"""Unit tests for nexus.io.reader.reader_factory.ReaderFactory."""

import pytest

from nexus.io.reader.lammps_reader import LAMMPSReader
from nexus.io.reader.reader_factory import ReaderFactory
from nexus.io.reader.xyz_reader import XYZReader


class TestGetReader:
    def test_xyz_file_selects_xyz_reader(self, two_atoms_path, settings_factory):
        factory = ReaderFactory(settings_factory(two_atoms_path))
        assert isinstance(factory.get_reader(), XYZReader)

    def test_missing_file_raises(self, settings_factory):
        factory = ReaderFactory(settings_factory("/nonexistent/missing.xyz"))
        with pytest.raises(ValueError):
            factory.get_reader()

    def test_lammps_extension_selects_lammps_reader(self, tmp_path, settings_factory):
        path = tmp_path / "traj.lammpstrj"
        path.write_text("dummy\n")
        factory = ReaderFactory(settings_factory(str(path)))
        assert isinstance(factory.get_reader(), LAMMPSReader)

    def test_unsupported_extension_returns_none(self, tmp_path, settings_factory):
        path = tmp_path / "notes.txt"
        path.write_text("dummy\n")
        factory = ReaderFactory(settings_factory(str(path)))
        assert factory.get_reader() is None
