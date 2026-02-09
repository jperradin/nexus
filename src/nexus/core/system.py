import numpy as np
from typing import List, Optional, Generator

from ..io.reader.base_reader import BaseReader
from .frame import Frame
from ..config.settings import Settings  # Import the Settings class


class System:
    """
    Manages trajectory data and provides frame-level access through a file reader.

    Wraps a reader with lazy frame iteration. On initialization it configures the
    reader's filename from settings and triggers a file scan to index frame offsets.
    Frames can then be accessed individually or iterated over as a generator.

    Attributes:
        reader (BaseReader): The file reader used to parse trajectory data.
        settings (Settings): Configuration parameters controlling frame range and file
            location.
        current_frame (Optional[Frame]): The most recently loaded frame, or None if no
            frame has been loaded yet.
    """

    def __init__(self, reader: BaseReader, settings: Settings):
        """
        Initialize the system with a reader and settings, then scan the trajectory file.

        Assigns the file location from settings to the reader and triggers the reader's
        ``scan()`` method to index frame byte offsets.

        Args:
            reader (BaseReader): The file reader instance to use for parsing frames.
            settings (Settings): Configuration settings containing file location and
                frame range.
        """
        self.reader: BaseReader = reader
        self.settings: Settings = settings
        self.current_frame: Optional[Frame] = None
        self._current_frame_index: Optional[int] = None # Index of the current frame
        self._num_frames: Optional[int] = None  # Cache for number of frames
        
        # Set the filename in the reader
        self.reader.filename = self.settings.file_location
        
        # Scan the file to initialize the reader
        if hasattr(self.reader, 'scan'):
            self.reader.scan()

    def load_frame(self, frame_index: int) -> bool:
        """
        Load a specific frame from the trajectory file into ``current_frame``.

        Validates the frame index against the range defined in settings, then delegates
        to the reader's ``parse()`` method.

        Args:
            frame_index (int): Zero-based index of the frame to load.

        Returns:
            bool: True if the frame was successfully loaded, False otherwise.

        Raises:
            ValueError: If ``frame_index`` is negative.
        """

        if frame_index < 0:
            raise ValueError("Frame index cannot be negative.")

        # Check the range from settings.
        start_frame, end_frame = self.settings.range_of_frames  # Unpack the tuple
        if not (start_frame <= frame_index <= (end_frame if end_frame != -1 else float('inf'))):
            print(f"Frame index {frame_index} is out of range specified in settings ({start_frame}-{end_frame}).")
            return False

        # Use the parse method to get the frame
        if hasattr(self.reader, 'parse'):
            try:
                # Get the frame using the parse method
                frame_generator = self.reader.parse(frame_index)
                frame = next(frame_generator)
                self.current_frame = frame
                self._current_frame_index = frame_index
                return True
            except (StopIteration, IndexError, ValueError) as e:
                print(f"Error loading frame {frame_index}: {str(e)}")
                return False
        
        print(f"Frame {frame_index} not found in trajectory.")
        return False


    def get_frame(self, frame_index: int) -> Optional[Frame]:
        """
        Retrieve a specific frame, loading it if necessary.

        Args:
            frame_index (int): Zero-based index of the frame to retrieve.

        Returns:
            Optional[Frame]: The loaded frame, or None if loading failed.
        """
        if self.load_frame(frame_index):
            return self.current_frame
        return None

    def get_num_frames(self) -> int:
        """
        Return the total number of frames in the trajectory.

        Uses the reader's ``num_frames`` attribute if available, otherwise counts frames
        by iterating through the trajectory. The result is cached for subsequent calls.

        Returns:
            int: Total number of frames, or 0 if an error occurs.
        """
        # First, check if we already calculated the number of frames
        if self._num_frames is not None:
            return self._num_frames
        
        # Use the reader's num_frames attribute if available
        if hasattr(self.reader, 'num_frames') and self.reader.num_frames > 0:
            self._num_frames = self.reader.num_frames
            return self._num_frames
            
        # Otherwise, count frames by iterating through them
        count = 0
        for _ in self.iter_frames():
            count += 1
        self._num_frames = count
        return count

    def iter_frames(self) -> Generator[Frame, None, None]:
        """
        Yield frames one at a time over the configured range.

        Generator-based iteration that avoids loading the entire trajectory into memory.
        Uses the reader's indexed frame offsets when available, falling back to sequential
        ``load_frame()`` calls otherwise. Respects the frame range defined in settings.

        Yields:
            Frame: The next frame in the trajectory.
        """
        start_frame, end_frame = self.settings.range_of_frames
        
        # If the reader has frame_indices, use them to iterate through frames
        if hasattr(self.reader, 'frame_indices') and self.reader.frame_indices:
            for frame_id in range(start_frame, min(end_frame + 1 if end_frame != -1 else float('inf'), len(self.reader.frame_indices))):
                try:
                    frame_generator = self.reader.parse(frame_id)
                    frame = next(frame_generator)
                    yield frame
                except (StopIteration, IndexError, ValueError) as e:
                    print(f"Error loading frame {frame_id}: {str(e)}")
                    continue
        else:
            # Fallback to loading frames one by one
            for frame_id in range(start_frame, end_frame + 1 if end_frame != -1 else float('inf')):
                if self.load_frame(frame_id):
                    yield self.current_frame  # type: ignore
                else:
                    break  # Stop if we can't load a frame


    def __iter__(self) -> 'System':
        """Reset the frame index and return self as an iterator."""
        self._current_frame_index = self.settings.range_of_frames[0]
        return self


    def __next__(self) -> Frame:
        """Return the next frame in the iteration sequence."""

        if self._current_frame_index is None:  # First call to next
             self._current_frame_index = self.settings.range_of_frames[0] # Initialize if needed.

        start_frame, end_frame = self.settings.range_of_frames

        if end_frame != -1 and self._current_frame_index > end_frame: #check if end frame is reach
            raise StopIteration

        if self._current_frame_index < self.get_num_frames():
            if self.load_frame(self._current_frame_index):
                self._current_frame_index += 1
                return self.current_frame  # type: ignore
            else: # If load frame return False
                raise StopIteration
        else: # if current frame index is greater than the number of frames.
            raise StopIteration
