import os
import json
from dataclasses import asdict
from datetime import datetime

from ...config.settings import Settings
from ...utils.performance import Performance
from ...io.writer.base_writer import BaseWriter

class PerformanceWriter(BaseWriter):
    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self._settings: Settings = settings

    def write(self, performance: Performance) -> None:
        # Convert datetime objects to strings for JSON serialization
        def serialize_datetime(obj):
            if isinstance(obj, datetime):
                return obj.isoformat()
            return obj
        
        # Convert performance data to dict
        perf_dict = asdict(performance)
        
        # Convert datetime objects
        perf_dict["timestamp"] = serialize_datetime(perf_dict["timestamp"])
        for entry in perf_dict["history"]:
            entry["timestamp"] = serialize_datetime(entry["timestamp"])
        
        # Save to file
        perf_file = os.path.join(self._settings.export_directory, f"performance_{performance.name}.json")
        with open(perf_file, "w") as f:
            json.dump(perf_dict, f, indent=2, default=serialize_datetime)