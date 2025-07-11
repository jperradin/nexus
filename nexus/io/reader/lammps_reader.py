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
        Scans the trajectory file.
        Initializes Frame objects with the chunk locations of each frame.
        Parse the header to store the number of nodes, the lattice and other informations.

        Returns:
            List[FrameIndex]: A list of FrameIndex objects.
        """
        if not self.mmaped_file:
            with open(self.filename, 'rb') as f:
                self.mmaped_file = memoryview(f.read())
        self.frame_offsets = []
        self.frame_sizes = []
        self.num_frames = 0
        self.frame_indices = []

        try:
            file_size = os.path.getsize(self.filename)
            with open(self.filename, 'r') as f:
                offset = 0
                while offset < file_size:
                    # Store position at the start of the frame
                    frame_start = offset

                    # Skip 3 header lines
                    f.readline()
                    f.readline()
                    f.readline()
                    offset += 3

                    # Read number of nodes line
                    num_nodes_line = f.readline().strip()
                    if not num_nodes_line:
                        break # End of file
                    try:
                        num_nodes = int(num_nodes_line)
                    except ValueError:
                        raise ValueError("Number of nodes must be an integer")
                    offset += 1

                    # Skip 1 header line
                    f.readline()
                    offset += 1

                    # Read lattice line
                    lattices = [f.readline().strip() for _ in range(3)]
                    try:
                        lattice = []
                        for i in range(3):
                            line = lattices[i].split()
                            low = float(line[0])
                            high = float(line[1])
                            lii = high - low
                            if i == 0:
                                lattice.append([lii, 0.0, 0.0])
                            elif i == 1:
                                lattice.append([0.0, lii, 0.0])
                            elif i == 2:
                                lattice.append([0.0, 0.0, lii])
                        lattice = np.array(lattice)
                    except ValueError:
                        raise ValueError("Lattice must be a 3x3 matrix")
                    offset += 3

                    # Read node lines
                    self.columns = f.readline().strip().split()
                    self.columns.remove("ITEM:")
                    self.columns.remove("NODES")
                    columns = {col: i for i, col in enumerate(self.columns)}
                    self.columns = columns
                    offset += 1
                    for _ in range(num_nodes):
                        f.readline()
                        offset += 1
                    
                    # Store position at the end of the frame
                    self.frame_indices.append(FrameIndex(frame_id=self.num_frames, num_nodes=num_nodes, lattice=lattice, byte_offset=frame_start))
                    self.num_frames += 1
                    self.frame_sizes.append(offset - frame_start)
                    self.frame_offsets.append(frame_start)

        except Exception as e:
            raise Exception(f"Error scanning trajectory file: {str(e)}")

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
        
        frame_index = self.frame_indices[frame_id]
        with open(self.filename, 'r') as f:
            self.seek_to_line(frame_index.byte_offset) 
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
                
                

