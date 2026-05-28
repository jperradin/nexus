"""Unit tests for nexus.core.system.System (lazy frame access over a reader)."""

import pytest

from nexus.core.frame import Frame
from nexus.core.system import System
from nexus.io.reader.xyz_reader import XYZReader


def make_system(path, settings_factory, **kwargs):
    settings = settings_factory(path, **kwargs)
    return System(XYZReader(settings), settings), settings


class TestFrameCount:
    def test_num_frames(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        assert system.get_num_frames() == 2

    def test_num_frames_cached(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        first = system.get_num_frames()
        # Mutate the reader; cached value must be returned unchanged.
        system.reader.num_frames = 99
        assert system.get_num_frames() == first


class TestLoadFrame:
    def test_valid_index_returns_true(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        assert system.load_frame(0) is True
        assert isinstance(system.current_frame, Frame)
        assert system.current_frame.frame_id == 0

    def test_negative_index_raises(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        with pytest.raises(ValueError):
            system.load_frame(-1)

    def test_out_of_settings_range_returns_false(self, two_frames_path, settings_factory):
        system, settings = make_system(two_frames_path, settings_factory)
        settings.range_of_frames = (0, 0)  # only frame 0 allowed
        assert system.load_frame(1) is False

    def test_beyond_file_returns_false(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        # Index in settings range but past the actual frames -> parse IndexError -> False.
        assert system.load_frame(5) is False


class TestGetFrame:
    def test_returns_frame(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        frame = system.get_frame(1)
        assert isinstance(frame, Frame)
        assert frame.frame_id == 1

    def test_missing_returns_none(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        assert system.get_frame(99) is None


class TestIteration:
    def test_iter_frames_yields_all(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        frames = list(system.iter_frames())
        assert len(frames) == 2
        assert [f.frame_id for f in frames] == [0, 1]

    def test_iterator_protocol(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        collected = [f.frame_id for f in system]
        assert collected == [0, 1]

    def test_iterator_stops(self, two_frames_path, settings_factory):
        system, _ = make_system(two_frames_path, settings_factory)
        it = iter(system)
        next(it)
        next(it)
        with pytest.raises(StopIteration):
            next(it)
