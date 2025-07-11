from typing import List, Generator
from collections import namedtuple
import numpy as np
import os

from .base_reader import BaseReader
from ...core.frame import Frame
from ...config.settings import Settings

FrameIndex = namedtuple('FrameIndex', ['frame_id', 'num_nodes', 'lattice', 'byte_offset'])

class LAMMPSReader(BaseReader):
    """
    Reader for LAMMPS trajectory files.
    """
    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)

    def detect(self, filepath: str) -> bool:
        """
        Detects if the file is supported by this reader.

        Returns:
            bool: True if the file is supported, False otherwise.
        """
        if filepath.lower().endswith('.lammpstrj'):
            return True
        if filepath.lower().endswith('.lammps'):
            return True
        if filepath.lower().endswith('.data'):
            return True
        return False

    def scan(self) -> List[FrameIndex]:
        """
        Scans the LAMMPS trajectory file efficiently to index frames.

        This method reads the file sequentially to locate the start of each
        frame and parse its header. It uses a single buffered reader for
        efficient I/O and stores the byte offset for fast seeking later.
        """
        self.frame_indices = []
        self.num_frames = 0
        
        try:
            with open(self.filename, 'r') as f:
                while True:
                    frame_start_offset = f.tell()
                    
                    # Read the first line of the frame block
                    line = f.readline()
                    if not line:
                        break # End of file
                    
                    # Check if the line is the 'ITEM: TIMESTEP' header
                    if 'ITEM: TIMESTEP' not in line:
                        # This could indicate a malformed file or we are out of sync.
                        # For now, we'll assume it means the end of valid frames.
                        break
                        
                    try:
                        # --- Parse Frame Header ---
                        f.readline() # Timestep value
                        f.readline() # ITEM: NUMBER OF NODES
                        num_nodes = int(f.readline().strip())
                        
                        f.readline() # ITEM: BOX BOUNDS
                        
                        # Read lattice vectors
                        lattices = [f.readline().strip() for _ in range(3)]
                        lattice = []
                        for i in range(3):
                            line_parts = lattices[i].split()
                            low, high = float(line_parts[0]), float(line_parts[1])
                            lii = high - low
                            if i == 0:
                                lattice.append([lii, 0.0, 0.0])
                            elif i == 1:
                                lattice.append([0.0, lii, 0.0])
                            else:
                                lattice.append([0.0, 0.0, lii])
                        lattice = np.array(lattice)
                        
                        # Read node property header
                        self.columns = {col: i for i, col in enumerate(f.readline().strip().split()[2:])}
                        
                        # Skip the atomic data lines to get to the next frame
                        for _ in range(num_nodes):
                            f.readline()

                    except (ValueError, IndexError) as e:
                        raise IOError(
                            f"Failed to parse LAMMPS frame header at byte offset {frame_start_offset} in {self.filename}. "
                            f"Error: {e}"
                        )

                    # Store the indexed frame information
                    frame_index = FrameIndex(
                        frame_id=self.num_frames,
                        num_nodes=num_nodes,
                        lattice=lattice,
                        byte_offset=frame_start_offset
                    )
                    self.frame_indices.append(frame_index)
                    self.num_frames += 1

        except FileNotFoundError:
            raise
        except Exception as e:
            raise IOError(f"Error scanning trajectory file {self.filename}: {e}")

        if self.verbose:
            print(f"Scanned {self.num_frames} frames in {self.filename}")

        self.is_indexed = True
        return self.frame_indices

    def parse(self, frame_id: int) -> Generator[Frame, None, None]:
        """
        Parses the trajectory file, get node data and yields frames.

        Yields:
            Frame: A data structure representing a frame.
        """
        
        if not self.is_indexed:
            self.scan()

        frame_index = self.frame_indices[frame_id]
        
        with open(self.filename, 'r') as f:
            f.seek(frame_index.byte_offset)
        
            num_nodes = frame_index.num_nodes
            lattice = frame_index.lattice

            # Skip 9 header lines
            for _ in range(9):
                f.readline()
            
            c_type = self.columns['type']
            c_x = self.columns['x']
            c_y = self.columns['y']
            c_z = self.columns['z']
            
            symbols = []
            positions = []
            # Read node lines
            for _ in range(num_nodes):
                node_line = f.readline().strip()
                try:
                    parts = node_line.split()
                    symbol = parts[c_type]
                    x = float(parts[c_x])
                    y = float(parts[c_y])
                    z = float(parts[c_z])
                    symbols.append(symbol)
                    positions.append(np.array([x, y, z]))
                except ValueError:
                    raise ValueError("Node line must have 4 values: symbol, x, y, z")

            data = {
                'symbol': symbols,
                'position': positions
            }
            
            yield Frame(frame_id=frame_id, _data=data, lattice=lattice, nodes=[])
                
                

