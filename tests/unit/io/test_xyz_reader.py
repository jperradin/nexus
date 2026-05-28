"""Unit tests for nexus.io.reader.xyz_reader.XYZReader."""

import os

import numpy as np
import pytest

from nexus.io.reader.xyz_reader import XYZReader


def reader_for(path, settings_factory, **kwargs):
    return XYZReader(settings_factory(path, **kwargs))


class TestDetect:
    def test_accepts_xyz_extension(self, default_settings):
        r = XYZReader(default_settings)
        assert r.detect("trajectory.xyz") is True
        assert r.detect("UPPER.XYZ") is True

    def test_rejects_other_extensions(self, default_settings):
        r = XYZReader(default_settings)
        assert r.detect("data.lammpstrj") is False
        assert r.detect("notes.txt") is False


class TestScan:
    def test_single_frame_count(self, two_atoms_path, settings_factory):
        r = reader_for(two_atoms_path, settings_factory)
        idx = r.scan()
        assert len(idx) == 1
        assert r.num_frames == 1

    def test_two_frame_count(self, two_frames_path, settings_factory):
        r = reader_for(two_frames_path, settings_factory)
        idx = r.scan()
        assert len(idx) == 2

    def test_frame_index_fields(self, two_atoms_path, settings_factory):
        r = reader_for(two_atoms_path, settings_factory)
        fi = r.scan()[0]
        assert fi.frame_id == 0
        assert fi.num_nodes == 2
        np.testing.assert_allclose(np.diag(fi.lattice), [10.0, 10.0, 10.0])

    def test_sets_is_indexed(self, two_atoms_path, settings_factory):
        r = reader_for(two_atoms_path, settings_factory)
        r.scan()
        assert r.is_indexed is True

    def test_malformed_header_raises(self, malformed_path, settings_factory):
        r = reader_for(malformed_path, settings_factory)
        with pytest.raises(IOError):
            r.scan()

    def test_missing_file_raises(self, settings_factory):
        r = reader_for("/nonexistent/path/missing.xyz", settings_factory)
        with pytest.raises(FileNotFoundError):
            r.scan()


class TestParse:
    def test_yields_frame_with_correct_id(self, two_frames_path, settings_factory):
        r = reader_for(two_frames_path, settings_factory)
        r.scan()
        frame = next(r.parse(1))
        assert frame.frame_id == 1

    def test_parse_scans_if_not_indexed(self, two_atoms_path, settings_factory):
        r = reader_for(two_atoms_path, settings_factory)
        # No explicit scan(); parse must trigger it.
        frame = next(r.parse(0))
        assert r.is_indexed is True
        assert frame is not None

    def test_raw_data_symbols_and_positions(self, two_atoms_path, settings_factory):
        r = reader_for(two_atoms_path, settings_factory)
        frame = next(r.parse(0))
        assert frame._data["symbol"] == ["A", "A"]
        np.testing.assert_allclose(
            np.array(frame._data["position"]), [[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]]
        )

    def test_lattice_parsed_into_frame(self, two_clusters_path, settings_factory):
        r = reader_for(two_clusters_path, settings_factory)
        frame = next(r.parse(0))
        np.testing.assert_allclose(np.diag(frame.lattice), [20.0, 20.0, 20.0])

    def test_lattice_written_back_to_settings(self, two_atoms_path, settings_factory):
        settings = settings_factory(two_atoms_path)
        r = XYZReader(settings)
        next(r.parse(0))
        # Not using a custom lattice -> reader stores the parsed lattice in settings.
        np.testing.assert_allclose(np.diag(settings.lattice.lattice), [10.0, 10.0, 10.0])

    def test_second_frame_distinct_positions(self, two_frames_path, settings_factory):
        r = reader_for(two_frames_path, settings_factory)
        f0 = next(r.parse(0))
        f1 = next(r.parse(1))
        np.testing.assert_allclose(np.array(f0._data["position"])[0], [1.0, 1.0, 1.0])
        np.testing.assert_allclose(np.array(f1._data["position"])[0], [1.5, 1.0, 1.0])
