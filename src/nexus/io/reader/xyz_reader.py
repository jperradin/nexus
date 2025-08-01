from typing import List, Generator
from collections import namedtuple
from colorama import Fore, Style
import numpy as np
import os

from .base_reader import BaseReader
from ...core.frame import Frame
from ...config.settings import Settings

FrameIndex = namedtuple(
    "FrameIndex", ["frame_id", "num_nodes", "lattice", "byte_offset"]
)


class XYZReader(BaseReader):
    """
    Reader for XYZ trajectory files.
    """

    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)

    def detect(self, filepath: str) -> bool:
        """
        Detects if the file is supported by this reader.

        Returns:
            bool: True if the file is supported, False otherwise.
        """
        return filepath.lower().endswith(".xyz")

    def scan(self) -> List[FrameIndex]:
        """
        Scans the trajectory file efficiently to index frames.

        This method reads the file sequentially to locate the start of each
        frame and parse its header. It uses a single buffered reader for
        efficient I/O and stores the byte offset of each frame for fast
        seeking later. The `mmaped_file` attribute is no longer needed for
        this operation.
        """
        self.frame_indices = []
        self.num_frames = 0

        try:
            # Use a single `with` statement for robust file handling.
            # Python's file object is already buffered and efficient.
            with open(self.filename, "r") as f:
                while True:
                    # Record the starting position of the potential frame.
                    frame_start_offset = f.tell()

                    num_nodes_line = f.readline()
                    if not num_nodes_line:
                        break  # End of file

                    header_line = f.readline()

                    try:
                        num_nodes = int(num_nodes_line.strip())

                        # Extract lattice information from the header
                        lattice_str = header_line.split('Lattice="')[1].split('"')[0]
                        parts = [float(p) for p in lattice_str.split()]
                        lattice = np.array(parts).reshape(3, 3)

                    except (ValueError, IndexError) as e:
                        # Provides more context if a header is malformed.
                        raise IOError(
                            f"Failed to parse frame header at byte offset {frame_start_offset} in {self.filename}. "
                            f"Ensure all frames have a number of nodes and a valid Lattice string. Error: {e}"
                        )

                    # Skip the atomic data to find the next frame's header
                    for _ in range(num_nodes):
                        f.readline()

                    # Store the indexed frame information.
                    frame_index = FrameIndex(
                        frame_id=self.num_frames,
                        num_nodes=num_nodes,
                        lattice=lattice,
                        byte_offset=frame_start_offset,
                    )
                    self.frame_indices.append(frame_index)
                    self.num_frames += 1

        except FileNotFoundError:
            raise
        except Exception as e:
            # Catch other potential I/O errors
            raise IOError(f"Error scanning trajectory file {self.filename}: {e}")

        if self.verbose:
            message = (
                Fore.LIGHTBLUE_EX
                + rf"""
   Scanned {self.num_frames} frames in {self.filename}
   Found {num_nodes} nodes
"""
                + Style.RESET_ALL
            )
            print(message)

        self.is_indexed = True
        # The `parse` method will use these FrameIndex objects to seek directly
        # to the correct position in the file to read a specific frame.
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

        with open(self.filename, "r") as f:
            f.seek(frame_index.byte_offset)

            num_nodes = frame_index.num_nodes
            lattice = frame_index.lattice

            # Skip 2 header lines
            f.readline()
            f.readline()

            symbols = []
            positions = []
            # Read node lines
            for _ in range(num_nodes):
                node_line = f.readline().strip()
                try:
                    parts = node_line.split()
                    symbol = parts[0]
                    x, y, z = map(float, parts[1:4])

                    symbols.append(symbol)
                    positions.append(np.array([x, y, z]))
                except ValueError:
                    raise ValueError("Node line must have 4 values: symbol, x, y, z")

            data = {"symbol": symbols, "position": positions}

            yield Frame(
                frame_id=frame_id,
                _data=data,
                lattice=lattice,
                nodes=[],
                _settings=self._settings,
            )
