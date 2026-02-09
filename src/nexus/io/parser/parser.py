from typing import List, Dict
import os

class Parser:
    """
    Discovers and indexes trajectory files in a directory for batch processing.

    Scans a directory for files matching a given format extension and reads an
    accompanying ``info.csv`` file to associate metadata (e.g., concentration,
    temperature) with each trajectory file.

    Attributes:
        file_location (str): Path to the directory or a file within it.
        format (str): File extension to match (e.g., ``"xyz"``).
        files (List[str]): Sorted list of discovered trajectory file paths.
        infos (Dict[str, List[float]]): Column-keyed metadata from ``info.csv``.
    """

    def __init__(self, file_location: str, format: str):
        """
        Initialize the parser, discover files, and load metadata.

        Args:
            file_location (str): Path to a directory or a file within one.
            format (str): File extension to filter by (without the leading dot).

        Raises:
            ValueError: If ``file_location`` does not exist.
        """
        self.file_location: str = file_location
        self.format: str = format
        self.files: List[str] = []
        self.infos: Dict[str, List[float]] = {}
        if os.path.exists(self.file_location):
            self.parse()
            self.parse_infos()
        else:
            raise ValueError(f"File location {self.file_location} does not exist")

    def parse(self) -> List[str]:
        """
        Discover trajectory files matching the configured format in the directory.

        If ``file_location`` is a file, its parent directory is scanned. If it is a
        directory, it is scanned directly. Results are stored in ``self.files``.

        Returns:
            List[str]: Sorted list of matching file paths.
        """
        if os.path.isfile(self.file_location):
            # take parent directory instead
            parent_directory = os.path.dirname(self.file_location)
            files = []
            for file in os.listdir(parent_directory):
                if file.endswith("." + self.format):
                    files.append(os.path.join(parent_directory, file))
            files.sort()
            self.files = files
        elif os.path.isdir(self.file_location):
            files = []
            for file in os.listdir(self.file_location):
                if file.endswith("." + self.format):
                    files.append(os.path.join(self.file_location, file))
            files.sort()
            self.files = files

    def parse_infos(self) -> None:
        """
        Load metadata from an ``info.csv`` file alongside the trajectory files.

        The CSV file must have a header row of column names and one data row per
        trajectory file. A ``project_name`` column is kept as strings; all other
        columns are parsed as floats.

        Raises:
            ValueError: If the file is missing, empty, or has a row count mismatch.
        """
        info_file = os.path.join(self.file_location, "info.csv")
        if os.path.exists(info_file):
            with open(info_file, "r") as f:
                lines = f.readlines()
                if len(lines) == 0:
                    raise ValueError("info.csv is empty")
                elif len(lines)-1 != len(self.files):
                    raise ValueError("info.csv does not have the same number of lines as trajectory files")
                
                keys = lines[0].strip().split(",")
                for k in keys:
                    self.infos[k] = []
                for line in lines[1:]:
                    values = line.strip().split(",")
                    for k, v in zip(keys, values):
                        if k == "project_name":
                            self.infos[k].append(v)
                        else:
                            self.infos[k].append(float(v))
        else:
            raise ValueError(f"info.csv not found in the directory : {self.file_location}")

    def get_files(self) -> List[str]:
        """
        Return the list of discovered trajectory files, scanning if needed.

        Returns:
            List[str]: Sorted list of trajectory file paths.
        """
        if not self.files:
            self.parse()
        return self.files

    def get_infos(self) -> Dict[str, List[float]]:
        """
        Return the metadata dictionary, loading from CSV if needed.

        Returns:
            Dict[str, List[float]]: Column-keyed metadata values.
        """
        if not self.infos:
            self.parse_infos()
        return self.infos
        
            
                    
                        
                
        