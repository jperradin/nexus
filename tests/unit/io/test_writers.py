"""Unit tests for the output writers and WriterFactory.

PerformanceWriter is the one that records wall-time / RSS memory / CPU per run;
its JSON output is what the benchmark memory reporting (roadmap Phase 4) reads.
"""

import json

from nexus.io.writer.clusters_writer import ClustersWriter
from nexus.io.writer.logs_writer import LogsWriter
from nexus.io.writer.performance_writer import PerformanceWriter
from nexus.io.writer.writer_factory import WriterFactory
from nexus.utils.performance import Performance


class TestWriterFactory:
    def test_get_performance_writer(self, default_settings):
        w = WriterFactory(default_settings).get_writer("PerformanceWriter")
        assert isinstance(w, PerformanceWriter)

    def test_get_logs_writer(self, default_settings):
        assert isinstance(WriterFactory(default_settings).get_writer("LogsWriter"), LogsWriter)

    def test_get_clusters_writer(self, default_settings):
        assert isinstance(WriterFactory(default_settings).get_writer("ClustersWriter"), ClustersWriter)

    def test_unknown_returns_none(self, default_settings):
        assert WriterFactory(default_settings).get_writer("Nope") is None


class TestPerformanceModel:
    def test_add_metric(self):
        p = Performance(id="1", name="run")
        p.add_metric("frames", 10)
        assert p.metrics["frames"] == 10

    def test_record_history_snapshots(self):
        p = Performance(id="1", name="run", execution_time_ms=50.0)
        p.record_history()
        assert len(p.history) == 1
        assert p.history[0]["execution_time_ms"] == 50.0

    def test_average_execution_time(self):
        p = Performance(id="1", name="run")
        p.execution_time_ms = 100.0
        p.record_history()
        p.execution_time_ms = 200.0
        p.record_history()
        assert p.get_average_execution_time() == 150.0


class TestPerformanceWriter:
    def test_writes_json_with_metrics(self, tmp_path, settings_factory, two_atoms_path):
        settings = settings_factory(two_atoms_path, export_directory=str(tmp_path))
        perf = Performance(
            id="1",
            name="bench",
            execution_time_ms=120.0,
            memory_usage_mb=42.5,
            cpu_usage_percent=75.0,
        )
        perf.record_history()

        PerformanceWriter(settings).write(perf)

        out = tmp_path / "performance_bench.json"
        assert out.exists()
        data = json.loads(out.read_text())
        assert data["name"] == "bench"
        assert data["memory_usage_mb"] == 42.5
        assert data["execution_time_ms"] == 120.0
        assert len(data["history"]) == 1
        # Timestamps must be serialized to ISO strings, not left as datetime.
        assert isinstance(data["timestamp"], str)
        assert isinstance(data["history"][0]["timestamp"], str)
