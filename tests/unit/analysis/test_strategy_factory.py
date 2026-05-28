"""Unit tests for nexus.analysis.strategy_factory.StrategyFactory selection logic."""

from nexus.analysis.strategy_factory import StrategyFactory
from nexus.analysis.strategies.bond_strategy import BondingStrategy
from nexus.analysis.strategies.coordination_strategy import CoordinationStrategy
from nexus.analysis.strategies.distance_strategy import DistanceStrategy
from nexus.analysis.strategies.shared_strategy import SharedStrategy


def factory(frame, settings):
    return StrategyFactory(frame, settings)


class TestStrategySelection:
    def test_distance_default(self, two_atom_frame, settings_factory, two_atoms_path):
        settings = settings_factory(two_atoms_path)
        strat = factory(two_atom_frame, settings).get_strategy(settings)
        assert isinstance(strat, DistanceStrategy)

    def test_bond_criterion(self, two_atom_frame, settings_factory, two_atoms_path):
        settings = settings_factory(two_atoms_path)
        settings.clustering.criterion = "bond"
        strat = factory(two_atom_frame, settings).get_strategy(settings)
        assert isinstance(strat, BondingStrategy)

    def test_coordination_number(self, two_atom_frame, settings_factory, two_atoms_path):
        settings = settings_factory(two_atoms_path)
        settings.clustering.with_coordination_number = True
        strat = factory(two_atom_frame, settings).get_strategy(settings)
        assert isinstance(strat, CoordinationStrategy)

    def test_shared_takes_precedence_over_coordination(self, two_atom_frame, settings_factory, two_atoms_path):
        settings = settings_factory(two_atoms_path)
        settings.clustering.with_coordination_number = True
        settings.clustering.with_number_of_shared = True
        strat = factory(two_atom_frame, settings).get_strategy(settings)
        assert isinstance(strat, SharedStrategy)


class TestRegistration:
    def test_all_four_registered(self, two_atom_frame, settings_factory, two_atoms_path):
        settings = settings_factory(two_atoms_path)
        f = factory(two_atom_frame, settings)
        for name in ("DistanceStrategy", "BondingStrategy", "SharedStrategy", "CoordinationStrategy"):
            assert name in f._strategies
