"""Analyzer validation on real amorphous-SiO2 data across a density sweep.

Runs the full pipeline on 5-frame samples of 384-atom v-SiO2 at three densities
(thesis dataset, B-part) using the production config (bond Si-O-Si,
coordination-O range [4,6], alternating), and checks the written analyzer
outputs.

The three densities span the tetrahedral -> octahedral coordination transition
of silica under compression, giving a physical oracle:

  * 2.212 g/cc (ambient): a fully connected corner-sharing SiO_4 network ->
    SiO_4-SiO_4 percolates, no 5/6-coordinated Si.
  * 3.006 g/cc (transition): SiO_4 still percolates; 5/6-coordination emerging
    but too sparse to percolate.
  * 4.109 g/cc (compressed): SiO_5 / SiO_6 networks dominate and percolate while
    sparse SiO_4 no longer does.

Two assertion layers:
  * Golden values: exact ensemble results for integer-derived analyzers
    (concentration, percolation, order parameter, largest/average size). These
    are deterministic and platform-independent (counts and their ratios).
  * Physics: percolation set per density and the monotonic rise of 6-coordinated
    Si with density. Gyration/correlation involve floating point and are only
    sanity-checked.
"""

import os
import tempfile

import pytest

from nexus import SettingsBuilder, main
import nexus.config.settings as c

CONNECTIVITIES = ["SiO_4-SiO_4", "SiO_4-SiO_5", "SiO_5-SiO_5", "SiO_5-SiO_6", "SiO_6-SiO_6"]

# Per-density expectations. ``percolating`` = connectivities with Pi == 1.0;
# ``conc`` = golden concentration values (subset checked exactly).
DENSITIES = {
    "2.212": {
        "file": "sio2_384_dens2.212.xyz",
        "percolating": {"SiO_4-SiO_4"},
        "conc": {"SiO_4-SiO_4": 1.0, "SiO_5-SiO_5": 0.0, "SiO_6-SiO_6": 0.0},
    },
    "3.006": {
        "file": "sio2_384_dens3.006.xyz",
        "percolating": {"SiO_4-SiO_4", "SiO_4-SiO_5"},
        "conc": {
            "SiO_4-SiO_4": 0.6671875,
            "SiO_4-SiO_5": 0.771875,
            "SiO_5-SiO_5": 0.2984375,
            "SiO_6-SiO_6": 0.015625,
        },
    },
    "4.109": {
        "file": "sio2_384_dens4.109.xyz",
        "percolating": {"SiO_5-SiO_5", "SiO_5-SiO_6", "SiO_6-SiO_6"},
        "conc": {"SiO_4-SiO_4": 0.034375, "SiO_5-SiO_6": 0.9203125, "SiO_6-SiO_6": 0.465625},
    },
}

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def run_pipeline(path):
    """Run the production SiO2 config on *path* and return the output directory."""
    tmp = tempfile.mkdtemp()
    general = c.GeneralSettings(
        project_name="sio2",
        export_directory=tmp,
        file_location=path,
        range_of_frames=(0, -1),
        apply_pbc=True,
        verbose=False,
    )
    clustering = c.ClusteringSettings(
        criterion="bond",
        node_types=["Si", "O"],
        node_masses=[28.0855, 15.9994],
        connectivity=["Si", "O", "Si"],
        cutoffs=[c.Cutoff(type1="Si", type2="O", distance=2.30)],
        with_coordination_number=True,
        coordination_mode="O",
        coordination_range=[4, 6],
        with_alternating=True,
    )
    settings = (
        SettingsBuilder()
        .with_general(general)
        .with_lattice(c.LatticeSettings(apply_custom_lattice=False))
        .with_clustering(clustering)
        .with_analysis(c.AnalysisSettings(with_all=True))
        .build()
    )
    main(settings)  # appends project_name to export_directory
    return os.path.join(tmp, "sio2")


def parse_dat(path):
    """Parse an analyzer .dat into {connectivity: [float, ...]} (columns after the label)."""
    rows = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        rows[parts[0]] = [float(x) for x in parts[1:]]
    return rows


@pytest.fixture(scope="module")
def all_outputs():
    """Run the pipeline once per density; return {density_name: output_dir}."""
    return {name: run_pipeline(os.path.join(DATA_DIR, info["file"])) for name, info in DENSITIES.items()}


def dat(all_outputs, name, filename):
    return parse_dat(os.path.join(all_outputs[name], filename))


# --------------------------------------------------------------------------- #
# Per-density golden + percolation checks
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", list(DENSITIES))
class TestPerDensity:
    def test_all_connectivities_present(self, all_outputs, name):
        rows = dat(all_outputs, name, "concentrations.dat")
        assert set(CONNECTIVITIES).issubset(rows)

    def test_concentration_golden(self, all_outputs, name):
        rows = dat(all_outputs, name, "concentrations.dat")
        for conn, val in DENSITIES[name]["conc"].items():
            assert rows[conn][0] == pytest.approx(val, rel=1e-9, abs=1e-12)

    def test_concentrations_valid_fractions(self, all_outputs, name):
        rows = dat(all_outputs, name, "concentrations.dat")
        for conn in CONNECTIVITIES:
            assert 0.0 <= rows[conn][0] <= 1.0

    def test_percolation_set(self, all_outputs, name):
        rows = dat(all_outputs, name, "percolation_probability.dat")
        expected = DENSITIES[name]["percolating"]
        for conn in CONNECTIVITIES:
            assert rows[conn][1] == (1.0 if conn in expected else 0.0)

    def test_order_parameter_matches_percolation(self, all_outputs, name):
        perc = dat(all_outputs, name, "percolation_probability.dat")
        order = dat(all_outputs, name, "order_parameter.dat")
        for conn in CONNECTIVITIES:
            if perc[conn][1] == 1.0:
                assert order[conn][1] > 0.0
            else:
                assert order[conn][1] == 0.0

    def test_average_size_nonnegative(self, all_outputs, name):
        # <S> excludes percolating clusters but finite clusters of the same
        # connectivity may remain, so it is only guaranteed non-negative.
        avg = dat(all_outputs, name, "average_cluster_size.dat")
        for conn in CONNECTIVITIES:
            assert avg[conn][1] >= 0.0

    def test_five_frames_averaged(self, all_outputs, name):
        text = open(os.path.join(all_outputs[name], "concentrations.dat")).read()
        assert "# Frames averaged: 5" in text


# --------------------------------------------------------------------------- #
# Cross-density physics: the coordination transition
# --------------------------------------------------------------------------- #
class TestCoordinationTransition:
    def test_six_coordination_rises_with_density(self, all_outputs):
        c6 = lambda n: dat(all_outputs, n, "concentrations.dat")["SiO_6-SiO_6"][0]
        assert c6("2.212") < c6("3.006") < c6("4.109")

    def test_four_coordination_falls_with_density(self, all_outputs):
        c4 = lambda n: dat(all_outputs, n, "concentrations.dat")["SiO_4-SiO_4"][0]
        assert c4("2.212") > c4("3.006") > c4("4.109")

    def test_percolation_pathway_shifts(self, all_outputs):
        # SiO_4 percolates at ambient density but not under high compression;
        # SiO_6 is the reverse.
        perc = lambda n, conn: dat(all_outputs, n, "percolation_probability.dat")[conn][1]
        assert perc("2.212", "SiO_4-SiO_4") == 1.0
        assert perc("4.109", "SiO_4-SiO_4") == 0.0
        assert perc("2.212", "SiO_6-SiO_6") == 0.0
        assert perc("4.109", "SiO_6-SiO_6") == 1.0


# --------------------------------------------------------------------------- #
# High-density detailed golden + float sanity
# --------------------------------------------------------------------------- #
class TestHighDensityGolden:
    def test_largest_cluster_sizes(self, all_outputs):
        rows = dat(all_outputs, "4.109", "largest_cluster_size.dat")
        assert rows["SiO_5-SiO_6"][1] == pytest.approx(117.8, rel=1e-6)
        assert rows["SiO_5-SiO_5"][1] == pytest.approx(58.8, rel=1e-6)
        assert rows["SiO_4-SiO_4"][1] == pytest.approx(2.4, rel=1e-6)

    def test_average_finite_cluster_golden(self, all_outputs):
        rows = dat(all_outputs, "4.109", "average_cluster_size.dat")
        assert rows["SiO_4-SiO_4"][1] == pytest.approx(2.2057142857142855, rel=1e-9)

    def test_percolating_connectivities_excluded_from_average(self, all_outputs):
        # At high density the percolating connectivities consist solely of the
        # giant spanning cluster, so <S> (finite clusters only) is exactly 0.
        rows = dat(all_outputs, "4.109", "average_cluster_size.dat")
        for conn in DENSITIES["4.109"]["percolating"]:
            assert rows[conn][1] == 0.0

    def test_correlation_length_nonnegative(self, all_outputs):
        rows = dat(all_outputs, "4.109", "correlation_length.dat")
        for conn in CONNECTIVITIES:
            assert rows[conn][1] >= 0.0

    def test_gyration_files_nonnegative(self, all_outputs):
        f = os.path.join(all_outputs["4.109"], "gyration_radius_distribution-SiO_4-SiO_4.dat")
        rows = [l for l in open(f) if l.strip() and not l.startswith("#")]
        assert rows
        assert all(float(r.split(",")[3]) >= 0.0 for r in rows)
