from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
from datetime import datetime


@dataclass
class Performance:
    """
    A dataclass representing performance metrics.
    
    This class can be used to store, analyze, and compare performance data
    for different operations or components in your system.
    """
    
    # Required fields
    id: str
    name: str
    timestamp: datetime = field(default_factory=datetime.now)
    
    # Optional performance metrics
    execution_time_ms: Optional[float] = None
    memory_usage_mb: Optional[float] = None
    cpu_usage_percent: Optional[float] = None
    
    # Additional metrics can be stored in a dictionary
    metrics: Dict[str, Any] = field(default_factory=dict)
    
    # Historical data points
    history: List[Dict[str, Any]] = field(default_factory=list)
    
    def add_metric(self, name: str, value: Any) -> None:
        """Add a custom metric to the performance data."""
        self.metrics[name] = value
    
    def record_history(self) -> None:
        """
        Record the current state in history.
        Useful for tracking performance changes over time.
        """
        current_state = {
            "timestamp": datetime.now(),
            "execution_time_ms": self.execution_time_ms,
            "memory_usage_mb": self.memory_usage_mb,
            "cpu_usage_percent": self.cpu_usage_percent,
            "metrics": self.metrics.copy()
        }
        self.history.append(current_state)
    
    def get_average_execution_time(self) -> Optional[float]:
        """Calculate the average execution time from history."""
        times = [entry["execution_time_ms"] for entry in self.history 
                if entry["execution_time_ms"] is not None]
        if not times:
            return None
        return sum(times) / len(times)
    
    def __str__(self) -> str:
        """String representation of performance data."""
        return (f"Performance '{self.name}' (ID: {self.id})\n"
                f"Timestamp: {self.timestamp}\n"
                f"Execution time: {self.execution_time_ms} ms\n"
                f"Memory usage: {self.memory_usage_mb} MB\n"
                f"CPU usage: {self.cpu_usage_percent}%\n"
                f"Additional metrics: {self.metrics}")