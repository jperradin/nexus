from typing import List, Dict
import os

class Parser:
    def __init__(self, file_location: str, format: str):
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
        # checktype of file_location
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
        # fetch a info.csv file containing the data associated with trajectory files
        # info.csv should contain any information that can be used to plot the results
        # ie probability, density, temperature, pressure, etc.
        # info.csv should have the same number of lines as trajectory files
        # info.csv should have the name of the columns at the first line separated by commas
        # info.csv should have the values of the columns separated by commas
        # info.csv should have a column named "project_name" with the name of the project for each trajectory file
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
        if not self.files:
            self.parse()
        return self.files

    def get_infos(self) -> Dict[str, List[float]]:
        if not self.infos:
            self.parse_infos()
        return self.infos
        
            
                    
                        
                
        