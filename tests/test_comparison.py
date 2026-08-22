from __future__ import annotations

import unittest

import matplotlib

matplotlib.use("Agg")

import orbital_data_exoplanets_and_earth as comparison


class AllPlanetsComparisonTests(unittest.TestCase):
    def test_non_earth_like_planet_is_retained(self) -> None:
        earth_rows = {
            0: {
                "case": "midpoint_eccentricity",
                "orbital_period_years": "1",
                "current_distance_au": "1",
                "orbital_speed_m_s": "29780",
                "host_mass_solar": "1",
                "planet_mass_earth": "1",
                "global_temp_k": "255",
                "north_pole_temp_k": "235",
                "equator_temp_k": "275",
                "south_pole_temp_k": "235",
            }
        }
        exoplanet_rows = [
            {
                "planet": "Extreme b",
                "host_star": "Extreme",
                "orbit_quarter": "0",
                "orbit_fraction": "0",
                "orbital_period_years": "1000",
                "current_distance_au": "500",
                "orbital_speed_m_s": "100",
                "host_mass_solar": "50",
                "planet_mass_earth": "1000000",
                "eccentricity": "0.9",
                "exoplanet_global_temp_k": "5000",
                "exoplanet_north_pole_temp_k": "4900",
                "exoplanet_equator_temp_k": "5100",
                "exoplanet_south_pole_temp_k": "4900",
            }
        ]

        results = comparison.build_results(exoplanet_rows, earth_rows)

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]["planet"], "Extreme b")
        self.assertNotIn("passed_earth_like_filter", results[0])
        self.assertNotIn("similarity_score", results[0])

    def test_only_exact_orbital_quarters_are_compared(self) -> None:
        earth_row = {
            "case": "midpoint_eccentricity",
            "orbital_period_years": "1",
            "current_distance_au": "1",
            "orbital_speed_m_s": "29780",
            "host_mass_solar": "1",
            "planet_mass_earth": "1",
            "global_temp_k": "255",
            "north_pole_temp_k": "235",
            "equator_temp_k": "275",
            "south_pole_temp_k": "235",
        }
        exoplanet_rows = [
            {"planet": "Fractional b", "orbit_quarter": "1.9"},
            {"planet": "Negative b", "orbit_quarter": "-0.5"},
            {"planet": "Out of range b", "orbit_quarter": "4"},
            {"planet": "Quarter 1 b", "orbit_quarter": "1"},
        ]

        results = comparison.build_results(exoplanet_rows, {0: earth_row, 1: earth_row})

        self.assertEqual([row["planet"] for row in results], ["Quarter 1 b"])


if __name__ == "__main__":
    unittest.main()
