from __future__ import annotations

import csv
import io
import tempfile
import unittest
from contextlib import redirect_stdout
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import SIMULATION_exoplanet_orbit_temperature as simulation


ROOT = Path(__file__).resolve().parents[1]
CATALOGUE = ROOT / "exoplanet_orbital_quarter_data.csv"


class CatalogueTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.planets = simulation.load_planets(CATALOGUE, forced_tilt_deg=None)

    def test_bundled_catalogue_is_selectable_and_deduplicated(self) -> None:
        self.assertEqual(len(self.planets), 6284)
        identities = {(planet.name.casefold(), planet.host.casefold()) for planet in self.planets}
        self.assertEqual(len(identities), len(self.planets))

    def test_planet_can_be_selected_by_name_or_number(self) -> None:
        by_name = simulation.resolve_planet_selector(self.planets, "Kepler-1988 b")
        by_number = simulation.resolve_planet_selector(self.planets, by_name.menu_index)
        self.assertIs(by_name, by_number)

    def test_default_menu_shows_every_loaded_planet(self) -> None:
        sample = [
            replace(self.planets[0], menu_index=index, name=f"Visible planet {index}")
            for index in range(1, 82)
        ]
        output = io.StringIO()

        with redirect_stdout(output):
            simulation.print_planet_menu(sample)

        self.assertIn("Visible planet 81", output.getvalue())
        self.assertNotIn("more not shown", output.getvalue())

    def test_loader_prefers_kepler_consistent_archival_candidate(self) -> None:
        planet = simulation.resolve_planet_selector(self.planets, "TOI-1694 c")

        self.assertIsNotNone(planet.source_period_days)
        self.assertLess(
            abs(planet.source_period_days / planet.period_days - 1.0),
            0.001,
        )


class ProfileTests(unittest.TestCase):
    def setUp(self) -> None:
        self.earth = simulation.earth_sun_planet()

    def test_profile_covers_and_closes_one_orbit(self) -> None:
        profile = simulation.orbital_temperature_profile(self.earth, samples=361)

        self.assertAlmostEqual(
            self.earth.period_days,
            simulation.period_from_a_if_missing(1.0, 1.0, 1.0),
            places=12,
        )
        self.assertEqual(len(profile), 361)
        self.assertEqual(profile[0]["orbit_fraction"], 0.0)
        self.assertEqual(profile[-1]["orbit_fraction"], 1.0)
        self.assertAlmostEqual(profile[0]["distance_au"], profile[-1]["distance_au"], places=12)
        self.assertAlmostEqual(profile[0]["global_temp_k"], profile[-1]["global_temp_k"], places=10)
        self.assertGreater(profile[0]["global_temp_k"], profile[len(profile) // 2]["global_temp_k"])

    def test_profile_rejects_too_few_samples(self) -> None:
        with self.assertRaises(ValueError):
            simulation.orbital_temperature_profile(self.earth, samples=2)

    def test_processed_reference_row_is_reconstructed_consistently(self) -> None:
        planet = simulation.row_to_exoplanet(
            {
                "planet": "Reference b",
                "host_star": "Reference",
                "orbital_period_years": "1",
                "current_distance_au": "0.08",
                "orbit_fraction": "0",
                "eccentricity": "0.2",
                "host_mass_solar": "1",
                "exoplanet_global_temp_k": "900",
            },
            row_index=0,
            menu_index=1,
            forced_tilt_deg=None,
        )

        self.assertIsNotNone(planet)
        assert planet is not None
        profile = simulation.orbital_temperature_profile(planet, samples=5)
        self.assertAlmostEqual(planet.semi_major_axis_au, 0.1, places=12)
        self.assertAlmostEqual(profile[0]["distance_au"], 0.08, places=12)
        self.assertAlmostEqual(profile[0]["global_temp_k"], 900.0, places=9)
        self.assertAlmostEqual(planet.source_period_days, simulation.DAYS_PER_YEAR)
        self.assertNotAlmostEqual(planet.source_period_days, planet.period_days, places=1)

    def test_independent_axis_candidate_beats_period_only_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            csv_path = Path(temporary_directory) / "candidates.csv"
            with csv_path.open("w", encoding="utf-8", newline="") as file:
                writer = csv.DictWriter(
                    file,
                    fieldnames=[
                        "planet",
                        "host_star",
                        "period_days",
                        "semi_major_axis_au",
                        "host_mass_solar",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "planet": "Candidate b",
                        "host_star": "Candidate",
                        "period_days": "365.25",
                        "semi_major_axis_au": "",
                        "host_mass_solar": "1",
                    }
                )
                writer.writerow(
                    {
                        "planet": "Candidate b",
                        "host_star": "Candidate",
                        "period_days": "365.25",
                        "semi_major_axis_au": "1",
                        "host_mass_solar": "1",
                    }
                )

            planets = simulation.load_planets(csv_path, forced_tilt_deg=None)
            self.assertEqual(len(planets), 1)
            self.assertEqual(planets[0].row_index, 1)

    def test_individual_and_comparison_outputs_are_written(self) -> None:
        comparison = simulation.Exoplanet(
            row_index=0,
            menu_index=2,
            name="Test b",
            host="Test star",
            period_days=100.0,
            semi_major_axis_au=0.5,
            eccentricity=0.2,
            eccentricity_min=0.2,
            eccentricity_max=0.2,
            planet_mass_earth=2.0,
            host_mass_solar=1.0,
            host_luminosity_solar=0.8,
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = simulation.generate_profile_outputs(
                [self.earth, comparison],
                Path(temporary_directory),
                samples=25,
                comparison_planets=[self.earth, comparison],
            )

            self.assertEqual(len(paths), 5)
            self.assertEqual(
                paths[-1].name,
                "temperature_profile_comparison_0001-0002.png",
            )
            for path in paths:
                self.assertTrue(path.is_file(), path)
                self.assertGreater(path.stat().st_size, 0, path)


if __name__ == "__main__":
    unittest.main()
