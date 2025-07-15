import numpy as np
from typing import List, Optional, Generator

from ..io.reader.base_reader import BaseReader
from .frame import Frame
from ..config.settings import Settings  # Import the Settings class


class System:
    """
    Manages the atomic system, trajectory data, and interaction with file readers.

    Attributes:
        reader (BaseReader): The file reader used to load data.
        settings (Settings): The settings object containing configuration parameters.
        current_frame (Optional[Frame]): The currently loaded frame.  None if no frame is loaded.
    """

    def __init__(self, reader: BaseReader, settings: Settings):
        """
        Initializes the System object.

        Args:
            reader (BaseReader): The file reader instance to use.
            settings (Settings): The settings object.
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
        Loads a specific frame from the trajectory file.

        Args:
            frame_index (int): The index of the frame to load (0-based).

        Returns:
            bool: True if the frame was successfully loaded, False otherwise.
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
        """Retrieves a specific frame, loading it if necessary.

        Args:
            frame_index: The index of the frame to retrieve.

        Returns:
            The Frame object, or None if the frame could not be loaded.
        """
        if self.load_frame(frame_index):
            return self.current_frame
        return None

    def get_num_frames(self) -> int:
        """
        Gets the total number of frames in the trajectory.
        If the value is already calculated, it's directly returned.

        Returns:
           int: The total number of frames, 0 if an error occurs
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
        Iterates through the frames of the trajectory, yielding one Frame at a time.
        This is a generator, avoiding loading the entire trajectory into memory.
        It respects the range of frames defined in settings.

        Yields:
            Frame: The next Frame object in the trajectory.
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
        """
        Make the System object itself iterable.
        """
        # Reset the frame index for iteration.
        self._current_frame_index = self.settings.range_of_frames[0]
        return self


    def __next__(self) -> Frame:
        """
        Returns the next frame during iteration.
        """

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