"""End-to-end integration tests for the full nexus.main pipeline.

Runs main(settings) on the small fixtures and checks that the analyzer output
files are written with the expected headline values. This exercises the whole
chain: reader -> system -> strategy -> union-find -> analyzers -> writers.

main() appends project_name to export_directory in place, so output lands in
``<tmp>/test`` (project_name defaults to "test" in the settings factory).
"""

import pytest

from nexus import main


def output_dir(tmp_path):
    return tmp_path / "test"


class TestPercolatingPipeline:
    def test_runs_and_writes_outputs(self, tmp_path, percolating_path, settings_factory):
        settings = settings_factory(percolating_path, export_directory=str(tmp_path))
        main(settings)

        out = output_dir(tmp_path)
        assert out.is_dir()
        # A representative selection of analyzer outputs must exist and be non-empty.
        for name in (
            "largest_cluster_size.dat",
            "percolation_probability.dat",
            "order_parameter.dat",
            "concentrations.dat",
        ):
            f = out / name
            assert f.exists() and f.stat().st_size > 0

    def test_largest_cluster_is_full_lattice(self, tmp_path, percolating_path, settings_factory):
        settings = settings_factory(percolating_path, export_directory=str(tmp_path))
        main(settings)
        content = (output_dir(tmp_path) / "largest_cluster_size.dat").read_text()
        assert "A-A" in content
        assert "27" in content  # 27-atom spanning cluster

    def test_percolation_detected(self, tmp_path, percolating_path, settings_factory):
        settings = settings_factory(percolating_path, export_directory=str(tmp_path))
        main(settings)
        # Data rows look like: A-A,<conc>,<Pi>,<std>,<err>
        rows = [
            line for line in (output_dir(tmp_path) / "percolation_probability.dat").read_text().splitlines()
            if line.startswith("A-A")
        ]
        assert rows
        pi = float(rows[0].split(",")[2])
        assert pi == 1.0


class TestFinitePipeline:
    def test_size_distribution_written(self, tmp_path, two_clusters_path, settings_factory):
        settings = settings_factory(two_clusters_path, export_directory=str(tmp_path))
        main(settings)
        dist = output_dir(tmp_path) / "cluster_size_distribution-A-A.dat"
        assert dist.exists()
        rows = [l for l in dist.read_text().splitlines() if l.startswith("A-A")]
        # Two clusters of size 3 -> a row for size 3 with count 2.
        size3 = [r for r in rows if r.split(",")[2] == "3"]
        assert size3
        assert int(float(size3[0].split(",")[3])) == 2


class TestErrorHandling:
    def test_missing_file_raises_runtime_error(self, tmp_path, settings_factory):
        settings = settings_factory("/nonexistent/missing.xyz", export_directory=str(tmp_path))
        with pytest.raises(RuntimeError):
            main(settings)
