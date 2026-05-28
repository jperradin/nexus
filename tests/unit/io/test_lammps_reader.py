"""Unit tests for nexus.io.reader.lammps_reader.LAMMPSReader.

detect(), scan(), and parse() are exercised against a small two-frame dump.
"""

import numpy as np
import pytest

from nexus.io.reader.lammps_reader import LAMMPSReader


class TestDetect:
    def test_accepts_lammps_extensions(self, default_settings):
        r = LAMMPSReader(default_settings)
        assert r.detect("dump.lammpstrj") is True
        assert r.detect("system.lammps") is True
        assert r.detect("config.data") is True

    def test_rejects_xyz(self, default_settings):
        assert LAMMPSReader(default_settings).detect("traj.xyz") is False


class TestScan:
    def test_two_frames(self, lammps_path, settings_factory):
        r = LAMMPSReader(settings_factory(lammps_path))
        idx = r.scan()
        assert len(idx) == 2
        assert r.num_frames == 2

    def test_frame_metadata(self, lammps_path, settings_factory):
        r = LAMMPSReader(settings_factory(lammps_path))
        fi = r.scan()[0]
        assert fi.num_nodes == 2
        # Box bounds 0..10 in each direction -> orthorhombic 10x10x10.
        np.testing.assert_allclose(np.diag(fi.lattice), [10.0, 10.0, 10.0])

    def test_columns_mapping(self, lammps_path, settings_factory):
        r = LAMMPSReader(settings_factory(lammps_path))
        r.scan()
        for col in ("type", "x", "y", "z"):
            assert col in r.columns

    def test_missing_file_raises(self, settings_factory):
        r = LAMMPSReader(settings_factory("/nonexistent/missing.lammpstrj"))
        with pytest.raises(FileNotFoundError):
            r.scan()


class TestParse:
    def test_parse_yields_frame_with_data(self, lammps_path, settings_factory):
        r = LAMMPSReader(settings_factory(lammps_path))
        r.scan()
        frame = next(r.parse(0))
        assert frame.frame_id == 0
        assert frame._data["symbol"] == ["A", "A"]
        np.testing.assert_allclose(
            np.array(frame._data["position"]), [[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]]
        )

    def test_parse_scans_if_not_indexed(self, lammps_path, settings_factory):
        r = LAMMPSReader(settings_factory(lammps_path))
        frame = next(r.parse(1))  # no explicit scan()
        assert r.is_indexed is True
        np.testing.assert_allclose(
            np.array(frame._data["position"])[0], [1.5, 1.0, 1.0]
        )

    def test_parse_lattice(self, lammps_path, settings_factory):
        r = LAMMPSReader(settings_factory(lammps_path))
        frame = next(r.parse(0))
        np.testing.assert_allclose(np.diag(frame.lattice), [10.0, 10.0, 10.0])
