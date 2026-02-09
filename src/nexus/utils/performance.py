from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
from datetime import datetime


@dataclass
class Performance:
    """
    Stores and tracks performance metrics for pipeline operations.

    Records execution time, memory usage, and CPU usage for a named operation.
    Supports custom metrics and historical snapshots for trend analysis.

    Attributes:
        id (str): Unique identifier for this performance record.
        name (str): Human-readable name of the tracked operation.
        timestamp (datetime): Creation time of this record.
        execution_time_ms (Optional[float]): Execution time in milliseconds.
        memory_usage_mb (Optional[float]): Peak memory usage in megabytes.
        cpu_usage_percent (Optional[float]): CPU utilization percentage.
        metrics (Dict[str, Any]): Custom key-value metrics.
        history (List[Dict[str, Any]]): Snapshots of past metric states.
    """

    id: str
    name: str
    timestamp: datetime = field(default_factory=datetime.now)

    execution_time_ms: Optional[float] = None
    memory_usage_mb: Optional[float] = None
    cpu_usage_percent: Optional[float] = None

    metrics: Dict[str, Any] = field(default_factory=dict)

    history: List[Dict[str, Any]] = field(default_factory=list)

    def add_metric(self, name: str, value: Any) -> None:
        """
        Add a custom metric to the performance data.

        Args:
            name (str): Name of the metric.
            value (Any): Value of the metric.
        """
        self.metrics[name] = value

    def record_history(self) -> None:
        """Record the current metric state as a snapshot in the history list."""
        current_state = {
            "timestamp": datetime.now(),
            "execution_time_ms": self.execution_time_ms,
            "memory_usage_mb": self.memory_usage_mb,
            "cpu_usage_percent": self.cpu_usage_percent,
            "metrics": self.metrics.copy()
        }
        self.history.append(current_state)
    
    def get_average_execution_time(self) -> Optional[float]:
        """
        Calculate the average execution time from recorded history.

        Returns:
            Optional[float]: Mean execution time in milliseconds, or None if no
                history entries contain execution time data.
        """
        times = [entry["execution_time_ms"] for entry in self.history 
                if entry["execution_time_ms"] is not None]
        if not times:
            return None
        return sum(times) / len(times)
    
    def __str__(self) -> str:
        """Return a human-readable summary of the performance record."""
        return (f"Performance '{self.name}' (ID: {self.id})\n"
                f"Timestamp: {self.timestamp}\n"
                f"Execution time: {self.execution_time_ms} ms\n"
                f"Memory usage: {self.memory_usage_mb} MB\n"
                f"CPU usage: {self.cpu_usage_percent}%\n"
                f"Additional metrics: {self.metrics}")