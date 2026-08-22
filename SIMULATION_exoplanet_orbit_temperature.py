"""Full-orbit exoplanet temperature profiles, comparisons, and Tk animation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

try:
    import tkinter as tk
except ImportError:  # Profile generation can still run headlessly without Tk.
    tk = None  # type: ignore[assignment]

import matplotlib

if tk is None:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt

# Physical constants
G = 6.67430e-11
AU_M = 1.495978707e11
M_SUN_KG = 1.98847e30
M_EARTH_KG = 5.9722e24
M_JUPITER_EARTH = 317.8
L_SUN_W = 3.828e26
SIGMA = 5.670374419e-8
DAYS_PER_YEAR = 365.25
SECONDS_PER_DAY = 86_400.0

# Display constants
WINDOW_WIDTH = 1320
WINDOW_HEIGHT = 820
SIDEBAR_WIDTH = 420
FPS_MS = 33
ORBIT_POINTS = 420
LATITUDE_BANDS = [-90, -60, -30, 0, 30, 60, 90]
PLANET_DISPLAY_SCALE = 1.5
STAR_DISPLAY_SCALE = 2.0
DEFAULT_PROFILE_SAMPLES = 361
DEFAULT_PROFILE_DIRECTORY = Path("orbital_temperature_profiles")

PROFILE_COLUMNS = [
    "planet",
    "host_star",
    "orbit_fraction",
    "time_days",
    "model_orbital_period_days",
    "source_orbital_period_days",
    "true_anomaly_deg",
    "distance_au",
    "eccentricity",
    "orbital_speed_m_s",
    "global_temp_k",
    "north_pole_temp_k",
    "equator_temp_k",
    "south_pole_temp_k",
]

# Data structures
@dataclass
class Exoplanet:
    row_index: int
    menu_index: int
    name: str
    host: str
    period_days: float
    semi_major_axis_au: float
    eccentricity: float
    eccentricity_min: float
    eccentricity_max: float
    planet_mass_earth: float
    host_mass_solar: float
    host_luminosity_solar: float
    albedo: float = 0.30
    axial_tilt_deg: float = 0.0
    source_period_days: float | None = None

    @property
    def host_mass_kg(self) -> float:
        return self.host_mass_solar * M_SUN_KG

    @property
    def planet_mass_kg(self) -> float:
        return self.planet_mass_earth * M_EARTH_KG

    @property
    def luminosity_w(self) -> float:
        return self.host_luminosity_solar * L_SUN_W

    @property
    def mu(self) -> float:
        return G * (self.host_mass_kg + self.planet_mass_kg)

    @property
    def has_varying_eccentricity(self) -> bool:
        return abs(self.eccentricity_max - self.eccentricity_min) > 1e-6


# ------------------------------------------------------------
# Built-in Earth/Sun fallback
# ------------------------------------------------------------

def earth_sun_planet() -> Exoplanet:
    """
    Built-in Earth/Sun example.

    Earth's eccentricity is not constant over geological time. It varies roughly
    between about 0.005 and 0.058 due to Milankovitch cycles. The animation uses
    those bounds but compresses the variation into a 30-second visual cycle.
    """
    return Exoplanet(
        row_index=-1,
        menu_index=1,
        name="Earth",
        host="Sun",
        period_days=period_from_a_if_missing(1.0, 1.0, 1.0),
        semi_major_axis_au=1.0,
        eccentricity=0.0167,
        eccentricity_min=0.005,
        eccentricity_max=0.058,
        planet_mass_earth=1.0,
        host_mass_solar=1.0,
        host_luminosity_solar=1.0,
        albedo=0.30,
        axial_tilt_deg=23.44,
        source_period_days=365.25,
    )

# CSV loading and cleaning
def parse_float(value: object, default: float = math.nan) -> float:
    if value is None:
        return default
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null", "--"}:
        return default
    try:
        return float(text)
    except ValueError:
        return default


def clamp_eccentricity(e: float) -> float:
    if not math.isfinite(e):
        return 0.0
    return min(max(e, 0.0), 0.95)

def first_valid_float(row: dict[str, str], names: list[str], default: float = math.nan) -> float:
    for name in names:
        value = parse_float(row.get(name), math.nan)
        if math.isfinite(value):
            return value
    return default

def host_luminosity_from_row(row: dict[str, str], albedo: float = 0.30) -> float:
    """
    Prefer st_lum_solar if your earlier script created it.
    Otherwise use NASA's st_lum, which is log10(L_star / L_sun).
    """
    lum_solar = parse_float(row.get("st_lum_solar"), math.nan)
    if math.isfinite(lum_solar) and lum_solar > 0:
        return lum_solar

    st_lum_log = parse_float(row.get("st_lum"), math.nan)
    if math.isfinite(st_lum_log):
        return 10.0 ** st_lum_log

    # Comparison/quarter-data files may contain a modelled temperature and
    # distance but no luminosity. Rearranging the same equilibrium-temperature
    # equation used by the simulator keeps those bundled files self-consistent.
    temperature_k = first_valid_float(
        row,
        ["exoplanet_global_temp_k", "global_temp_k"],
        math.nan,
    )
    distance_au = parse_float(row.get("current_distance_au"), math.nan)
    if (
        math.isfinite(temperature_k)
        and temperature_k > 0
        and math.isfinite(distance_au)
        and distance_au > 0
        and 0 <= albedo < 1
    ):
        distance_m = distance_au * AU_M
        luminosity_w = (
            temperature_k ** 4
            * 16.0
            * math.pi
            * SIGMA
            * distance_m ** 2
            / (1.0 - albedo)
        )
        return luminosity_w / L_SUN_W

    return 1.0

def eccentricity_from_row(row: dict[str, str]) -> float:
    """
    Prefer the column made by the statistical eccentricity program.
    Then fall back to the expected value, then the raw NASA value.
    """
    eccentricity = first_valid_float(
        row,
        ["eccentricity_for_analysis", "expected_eccentricity", "pl_orbeccen", "eccentricity"],
        default=0.0,
    )
    return clamp_eccentricity(eccentricity)

def eccentricity_limits_from_row(row: dict[str, str], base_e: float) -> tuple[float, float]:
    """
    Uses explicit e_min/e_max columns if available.

    If those do not exist, uses NASA-style uncertainty columns if they are
    present. In the NASA archive, lower uncertainties are often stored as
    negative values, so both signs are handled.

    If no limits exist, returns a static range: (base_e, base_e).
    """
    explicit_min = first_valid_float(
        row,
        ["eccentricity_min", "e_min", "pl_orbeccen_min"],
        default=math.nan,
    )
    explicit_max = first_valid_float(
        row,
        ["eccentricity_max", "e_max", "pl_orbeccen_max"],
        default=math.nan,
    )

    if math.isfinite(explicit_min) and math.isfinite(explicit_max):
        e_min = clamp_eccentricity(explicit_min)
        e_max = clamp_eccentricity(explicit_max)
        return min(e_min, e_max), max(e_min, e_max)

    err1 = parse_float(row.get("pl_orbeccenerr1"), math.nan)
    err2 = parse_float(row.get("pl_orbeccenerr2"), math.nan)

    candidates = [base_e]
    if math.isfinite(err1):
        candidates.append(base_e + abs(err1))
    if math.isfinite(err2):
        # err2 is commonly negative. abs(err2) also handles positive lower errors.
        candidates.append(base_e - abs(err2))

    candidates = [clamp_eccentricity(value) for value in candidates if math.isfinite(value)]

    if len(candidates) >= 2:
        return min(candidates), max(candidates)

    return base_e, base_e

def planet_mass_from_row(row: dict[str, str]) -> float:
    mass_earth = first_valid_float(row, ["planet_mass_earth", "pl_bmasse"], math.nan)
    if math.isfinite(mass_earth) and mass_earth > 0:
        return mass_earth

    mass_jupiter = parse_float(row.get("pl_bmassj"), math.nan)
    if math.isfinite(mass_jupiter) and mass_jupiter > 0:
        return mass_jupiter * M_JUPITER_EARTH

    return 1.0

def period_from_a_if_missing(
    semi_major_axis_au: float,
    host_mass_solar: float,
    planet_mass_earth: float = 0.0,
) -> float:
    """
    Two-body Kepler period for a known semi-major axis and component masses.
    """
    if semi_major_axis_au <= 0 or host_mass_solar <= 0 or planet_mass_earth < 0:
        return math.nan
    semi_major_axis_m = semi_major_axis_au * AU_M
    total_mass_kg = host_mass_solar * M_SUN_KG + planet_mass_earth * M_EARTH_KG
    return (
        2.0
        * math.pi
        * math.sqrt(semi_major_axis_m ** 3 / (G * total_mass_kg))
        / SECONDS_PER_DAY
    )

def a_from_period_if_missing(
    period_days: float,
    host_mass_solar: float,
    planet_mass_earth: float = 0.0,
) -> float:
    """
    Two-body semi-major axis for a known period and component masses.
    """
    if period_days <= 0 or host_mass_solar <= 0 or planet_mass_earth < 0:
        return math.nan
    period_seconds = period_days * SECONDS_PER_DAY
    total_mass_kg = host_mass_solar * M_SUN_KG + planet_mass_earth * M_EARTH_KG
    semi_major_axis_m = (
        G * total_mass_kg * period_seconds ** 2 / (4.0 * math.pi ** 2)
    ) ** (1.0 / 3.0)
    return semi_major_axis_m / AU_M


def a_from_reference_radius_if_available(
    row: dict[str, str],
    eccentricity: float,
) -> float:
    """Recover ``a`` from a processed row's radius and elapsed-time phase."""
    radius_au = parse_float(row.get("current_distance_au"), math.nan)
    orbit_fraction = parse_float(row.get("orbit_fraction"), math.nan)
    if (
        not math.isfinite(radius_au)
        or radius_au <= 0
        or not math.isfinite(orbit_fraction)
    ):
        return math.nan

    mean_anomaly = 2.0 * math.pi * (orbit_fraction % 1.0)
    eccentric_anomaly = solve_kepler(mean_anomaly, eccentricity)
    radius_factor = 1.0 - eccentricity * math.cos(eccentric_anomaly)
    if radius_factor <= 0:
        return math.nan
    return radius_au / radius_factor

def row_to_exoplanet(
    row: dict[str, str],
    row_index: int,
    menu_index: int,
    forced_tilt_deg: float | None,
) -> Exoplanet | None:
    name = (row.get("pl_name") or row.get("planet") or f"row_{row_index}").strip()
    host = (row.get("hostname") or row.get("host_star") or "Unknown host").strip()

    host_mass_solar = first_valid_float(row, ["st_mass", "host_mass_solar"], 1.0)
    if not math.isfinite(host_mass_solar) or host_mass_solar <= 0:
        host_mass_solar = 1.0

    planet_mass_earth = planet_mass_from_row(row)
    base_e = eccentricity_from_row(row)

    source_period_days = first_valid_float(row, ["pl_orbper", "period_days"], math.nan)
    if not math.isfinite(source_period_days) or source_period_days <= 0:
        period_years = parse_float(row.get("orbital_period_years"), math.nan)
        if math.isfinite(period_years) and period_years > 0:
            source_period_days = period_years * DAYS_PER_YEAR

    period_days = source_period_days

    semi_major_axis_au = first_valid_float(row, ["pl_orbsmax", "semi_major_axis_au"], math.nan)

    if not math.isfinite(semi_major_axis_au) or semi_major_axis_au <= 0:
        semi_major_axis_au = a_from_reference_radius_if_available(row, base_e)

    if not math.isfinite(period_days) or period_days <= 0:
        period_days = period_from_a_if_missing(
            semi_major_axis_au,
            host_mass_solar,
            planet_mass_earth,
        )

    if not math.isfinite(semi_major_axis_au) or semi_major_axis_au <= 0:
        semi_major_axis_au = a_from_period_if_missing(
            period_days,
            host_mass_solar,
            planet_mass_earth,
        )

    if not math.isfinite(period_days) or period_days <= 0:
        return None
    if not math.isfinite(semi_major_axis_au) or semi_major_axis_au <= 0:
        return None

    # The selected orbital scale and masses define the internally consistent
    # model period. Keep a conflicting catalogue period as source metadata.
    period_days = period_from_a_if_missing(
        semi_major_axis_au,
        host_mass_solar,
        planet_mass_earth,
    )

    axial_tilt_deg = first_valid_float(
        row,
        ["axial_tilt_deg", "obliquity_deg", "pl_obliq"],
        default=0.0,
    )
    if forced_tilt_deg is not None:
        axial_tilt_deg = forced_tilt_deg
    if not math.isfinite(axial_tilt_deg):
        axial_tilt_deg = 0.0

    e_min, e_max = eccentricity_limits_from_row(row, base_e)

    albedo = first_valid_float(row, ["albedo", "planet_albedo"], 0.30)
    if not math.isfinite(albedo) or not 0 <= albedo < 1:
        albedo = 0.30

    return Exoplanet(
        row_index=row_index,
        menu_index=menu_index,
        name=name,
        host=host,
        period_days=period_days,
        semi_major_axis_au=semi_major_axis_au,
        eccentricity=base_e,
        eccentricity_min=e_min,
        eccentricity_max=e_max,
        planet_mass_earth=planet_mass_earth,
        host_mass_solar=host_mass_solar,
        host_luminosity_solar=host_luminosity_from_row(row, albedo),
        albedo=albedo,
        axial_tilt_deg=axial_tilt_deg,
        source_period_days=(
            source_period_days
            if math.isfinite(source_period_days) and source_period_days > 0
            else None
        ),
    )


def source_period_consistency_score(planet: Exoplanet) -> float:
    """Smaller is better; compare source and dynamically consistent periods."""
    source_period = planet.source_period_days
    if source_period is None or not math.isfinite(source_period) or source_period <= 0:
        return math.inf
    return abs(math.log(source_period / planet.period_days))


def candidate_selection_rank(
    row: dict[str, str],
    planet: Exoplanet,
) -> tuple[int, float]:
    """Prefer independently checkable orbital scales, then period consistency."""
    explicit_axis = first_valid_float(
        row,
        ["pl_orbsmax", "semi_major_axis_au"],
        math.nan,
    )
    reference_radius = parse_float(row.get("current_distance_au"), math.nan)
    reference_phase = parse_float(row.get("orbit_fraction"), math.nan)
    has_independent_scale = (
        math.isfinite(explicit_axis)
        and explicit_axis > 0
    ) or (
        math.isfinite(reference_radius)
        and reference_radius > 0
        and math.isfinite(reference_phase)
    )
    has_source_period = (
        planet.source_period_days is not None
        and math.isfinite(planet.source_period_days)
        and planet.source_period_days > 0
    )

    if has_independent_scale and has_source_period:
        return 0, source_period_consistency_score(planet)
    if has_independent_scale:
        return 1, 0.0
    if has_source_period:
        return 2, 0.0
    return 3, 0.0


def load_planets(csv_path: Path, forced_tilt_deg: float | None) -> list[Exoplanet]:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file)
        planets: list[Exoplanet] = []
        planet_positions: dict[tuple[str, str], int] = {}
        best_ranks: list[tuple[int, float]] = []
        for row_index, row in enumerate(reader):
            planet = row_to_exoplanet(row, row_index, len(planets) + 1, forced_tilt_deg)
            if planet is None:
                continue

            key = (planet.name.casefold(), planet.host.casefold())
            rank = candidate_selection_rank(row, planet)
            position = planet_positions.get(key)
            if position is None:
                planet_positions[key] = len(planets)
                planets.append(planet)
                best_ranks.append(rank)
            elif (
                rank[0] < best_ranks[position][0]
                or (
                    rank[0] == best_ranks[position][0]
                    and rank[1] < best_ranks[position][1] - 1e-12
                )
            ):
                planet.menu_index = planets[position].menu_index
                planets[position] = planet
                best_ranks[position] = rank
    return planets

def print_planet_menu(planets: list[Exoplanet], limit: int | None = None) -> None:
    print("\nChoose a planet:")
    print("number | planet | host | model P/d | source P/d | a / AU | e")
    print("-" * 108)
    displayed_planets = planets if limit is None else planets[:limit]
    for planet in displayed_planets:
        if planet.has_varying_eccentricity:
            e_text = f"{planet.eccentricity_min:.3f}-{planet.eccentricity_max:.3f}"
        else:
            e_text = f"{planet.eccentricity:.3f}"
        source_period_text = (
            f"{planet.source_period_days:10.3f}"
            if planet.source_period_days is not None
            else f"{'-':>10s}"
        )
        print(
            f"{planet.menu_index:6d} | "
            f"{planet.name[:24]:24s} | "
            f"{planet.host[:20]:20s} | "
            f"{planet.period_days:9.3f} | "
            f"{source_period_text} | "
            f"{planet.semi_major_axis_au:7.4f} | "
            f"{e_text:>11s}"
        )
    if limit is not None and len(planets) > limit:
        print(f"... {len(planets) - limit} more not shown. Use --choice N to choose one directly.")
    print()

# Orbital mechanics
def solve_kepler(mean_anomaly: float, eccentricity: float, iterations: int = 12) -> float:
    E = mean_anomaly if eccentricity < 0.8 else math.pi
    for _ in range(iterations):
        f = E - eccentricity * math.sin(E) - mean_anomaly
        fp = 1.0 - eccentricity * math.cos(E)
        if abs(fp) < 1e-12:
            break
        E -= f / fp
    return E

def true_anomaly_from_time(time_days: float, planet: Exoplanet, eccentricity: float) -> float:
    M = 2.0 * math.pi * ((time_days / planet.period_days) % 1.0)
    E = solve_kepler(M, eccentricity)
    numerator = math.sqrt(1.0 + eccentricity) * math.sin(E / 2.0)
    denominator = math.sqrt(1.0 - eccentricity) * math.cos(E / 2.0)
    return 2.0 * math.atan2(numerator, denominator)

def relative_position_au(
    planet: Exoplanet,
    true_anomaly: float,
    eccentricity: float,
) -> tuple[float, float, float]:
    """
    Returns the planet position relative to a stationary host star at one focus.
    """
    a = planet.semi_major_axis_au
    r = a * (1.0 - eccentricity ** 2) / (1.0 + eccentricity * math.cos(true_anomaly))
    x = r * math.cos(true_anomaly)
    y = r * math.sin(true_anomaly)
    return x, y, r

def orbital_speed_m_s(planet: Exoplanet, radius_au: float) -> float:
    r_m = radius_au * AU_M
    a_m = planet.semi_major_axis_au * AU_M
    return math.sqrt(planet.mu * (2.0 / r_m - 1.0 / a_m))

def periapsis_au(planet: Exoplanet, eccentricity: float) -> float:
    return planet.semi_major_axis_au * (1.0 - eccentricity)

def apoapsis_au(planet: Exoplanet, eccentricity: float) -> float:
    return planet.semi_major_axis_au * (1.0 + eccentricity)

# Temperature model
def global_equilibrium_temperature_k(planet: Exoplanet, radius_au: float) -> float:
    """
    Black-body equilibrium temperature with full surface redistribution:
        T = [ L(1-A)/(16 pi sigma r^2) ]^(1/4)
    """
    r_m = radius_au * AU_M
    return ((planet.luminosity_w * (1.0 - planet.albedo)) / (16.0 * math.pi * SIGMA * r_m ** 2)) ** 0.25

def substellar_latitude_deg(planet: Exoplanet, time_days: float) -> float:
    """
    Simple seasonal tilt model. If axial tilt is unavailable, this is zero.
    """
    orbital_phase = 2.0 * math.pi * ((time_days / planet.period_days) % 1.0)
    return planet.axial_tilt_deg * math.sin(orbital_phase)

def local_temperature_k(global_temp_k: float, latitude_deg: float, substellar_lat_deg: float) -> float:
    """
    Educational regional-temperature approximation.

    This is not a full atmospheric model. It shifts the warmest band toward the
    substellar latitude and cools polar regions.
    """
    angular_distance = abs(latitude_deg - substellar_lat_deg)
    seasonal_heating = 16.0 * max(math.cos(math.radians(angular_distance)), 0.0)
    polar_cooling = 13.0 * (abs(latitude_deg) / 90.0) ** 1.15
    return global_temp_k + seasonal_heating - polar_cooling


def calculate_orbital_state(
    planet: Exoplanet,
    time_days: float,
    eccentricity: float | None = None,
) -> dict[str, float]:
    """Calculate orbital, speed, and temperature values at one simulated time."""
    active_eccentricity = planet.eccentricity if eccentricity is None else clamp_eccentricity(eccentricity)
    true_anomaly = true_anomaly_from_time(time_days, planet, active_eccentricity)
    x_au, y_au, radius_au = relative_position_au(planet, true_anomaly, active_eccentricity)
    speed_m_s = orbital_speed_m_s(planet, radius_au)
    global_temp_k = global_equilibrium_temperature_k(planet, radius_au)
    substellar_latitude = substellar_latitude_deg(planet, time_days)
    local_temperatures = {
        latitude: local_temperature_k(global_temp_k, latitude, substellar_latitude)
        for latitude in LATITUDE_BANDS
    }

    return {
        "nu": true_anomaly,
        "x": x_au,
        "y": y_au,
        "r": radius_au,
        "speed": speed_m_s,
        "global_k": global_temp_k,
        "sub_lat": substellar_latitude,
        "north_pole_k": local_temperatures[90],
        "equator_k": local_temperatures[0],
        "south_pole_k": local_temperatures[-90],
        "peri_au": periapsis_au(planet, active_eccentricity),
        "apo_au": apoapsis_au(planet, active_eccentricity),
        **{f"lat_{latitude}": temp for latitude, temp in local_temperatures.items()},
    }


def orbital_temperature_profile(
    planet: Exoplanet,
    samples: int = DEFAULT_PROFILE_SAMPLES,
    eccentricity: float | None = None,
) -> list[dict[str, object]]:
    """
    Sample a planet at evenly spaced times from periapsis through one full orbit.

    Both endpoints are included, so orbit fractions 0 and 1 describe the same
    physical position and make the plotted profile visibly close on itself.
    """
    if samples < 3:
        raise ValueError("A temperature profile needs at least three samples.")

    active_eccentricity = planet.eccentricity if eccentricity is None else clamp_eccentricity(eccentricity)
    profile: list[dict[str, object]] = []

    for sample_index in range(samples):
        orbit_fraction = sample_index / (samples - 1)
        time_days = orbit_fraction * planet.period_days
        state = calculate_orbital_state(planet, time_days, active_eccentricity)
        profile.append(
            {
                "planet": planet.name,
                "host_star": planet.host,
                "orbit_fraction": orbit_fraction,
                "time_days": time_days,
                "model_orbital_period_days": planet.period_days,
                "source_orbital_period_days": (
                    planet.source_period_days
                    if planet.source_period_days is not None
                    else ""
                ),
                "true_anomaly_deg": math.degrees(state["nu"]) % 360.0,
                "distance_au": state["r"],
                "eccentricity": active_eccentricity,
                "orbital_speed_m_s": state["speed"],
                "global_temp_k": state["global_k"],
                "north_pole_temp_k": state["north_pole_k"],
                "equator_temp_k": state["equator_k"],
                "south_pole_temp_k": state["south_pole_k"],
            }
        )

    return profile


def safe_filename(value: str) -> str:
    """Return a portable, human-readable filename component."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", value.strip()).strip("-._")
    return cleaned.lower() or "planet"


def profile_file_stem(planet: Exoplanet) -> str:
    return f"{planet.menu_index:04d}_{safe_filename(planet.name)}"


def write_temperature_profile(path: Path, profile: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=PROFILE_COLUMNS)
        writer.writeheader()
        writer.writerows(profile)


def make_temperature_profile_graph(
    planet: Exoplanet,
    profile: list[dict[str, object]],
    output_path: Path,
) -> None:
    """Save temperature and distance curves covering one complete orbit."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fractions = [float(row["orbit_fraction"]) for row in profile]

    figure, (temperature_axis, distance_axis) = plt.subplots(
        2,
        1,
        figsize=(10, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [2.0, 1.0]},
    )

    temperature_series = [
        ("Global", "global_temp_k", "#202020", 2.3),
        ("North pole", "north_pole_temp_k", "#2878b5", 1.4),
        ("Equator", "equator_temp_k", "#e07a1f", 1.4),
        ("South pole", "south_pole_temp_k", "#3ca370", 1.4),
    ]
    for label, column, colour, width in temperature_series:
        temperature_axis.plot(
            fractions,
            [float(row[column]) for row in profile],
            label=label,
            color=colour,
            linewidth=width,
        )

    if planet.has_varying_eccentricity:
        minimum_profile = orbital_temperature_profile(
            planet,
            len(profile),
            planet.eccentricity_min,
        )
        maximum_profile = orbital_temperature_profile(
            planet,
            len(profile),
            planet.eccentricity_max,
        )
        lower = [
            min(float(low["global_temp_k"]), float(high["global_temp_k"]))
            for low, high in zip(minimum_profile, maximum_profile)
        ]
        upper = [
            max(float(low["global_temp_k"]), float(high["global_temp_k"]))
            for low, high in zip(minimum_profile, maximum_profile)
        ]
        temperature_axis.fill_between(
            fractions,
            lower,
            upper,
            color="#777777",
            alpha=0.18,
            label="Global range from eccentricity limits",
        )

    temperature_axis.set_ylabel("Temperature / K")
    temperature_axis.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
    temperature_axis.legend(loc="best", ncols=2)

    distance_axis.plot(
        fractions,
        [float(row["distance_au"]) for row in profile],
        color="#7554a5",
        linewidth=2.0,
    )
    distance_axis.set_xlabel("Orbital phase (fraction of one period)")
    distance_axis.set_ylabel("Star distance / AU")
    distance_axis.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
    distance_axis.set_xlim(0.0, 1.0)

    eccentricity_text = f"e = {planet.eccentricity:.4f}"
    period_text = f"model P = {planet.period_days:.4g} days"
    if (
        planet.source_period_days is not None
        and abs(math.log(planet.source_period_days / planet.period_days)) > 1e-3
    ):
        period_text += f", source P = {planet.source_period_days:.4g} days"
    figure.suptitle(
        f"{planet.name} around {planet.host}: full-orbit temperature profile\n"
        f"{period_text}, a = {planet.semi_major_axis_au:.4g} AU, {eccentricity_text}"
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def make_temperature_comparison_graph(
    profiles: list[tuple[Exoplanet, list[dict[str, object]]]],
    output_path: Path,
) -> None:
    """Overlay nominal global-temperature profiles on a common orbital phase."""
    if len(profiles) < 2:
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(10, 6))
    for planet, profile in profiles:
        axis.plot(
            [float(row["orbit_fraction"]) for row in profile],
            [float(row["global_temp_k"]) for row in profile],
            linewidth=2.0,
            label=f"{planet.name} ({planet.host})",
        )

    axis.set_xlabel("Orbital phase (fraction of each planet's period)")
    axis.set_ylabel("Global equilibrium temperature / K")
    axis.set_title("Full-orbit exoplanet temperature comparison")
    axis.set_xlim(0.0, 1.0)
    axis.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
    axis.legend(loc="best")
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def comparison_graph_filename(planets: list[Exoplanet]) -> str:
    """Build a deterministic comparison filename without unbounded path length."""
    choices = "-".join(f"{planet.menu_index:04d}" for planet in planets)
    if len(choices) > 80:
        digest = hashlib.sha256(choices.encode("ascii")).hexdigest()[:12]
        choices = f"{len(planets)}-planets-{digest}"
    return f"temperature_profile_comparison_{choices}.png"


def generate_profile_outputs(
    planets: list[Exoplanet],
    output_directory: Path,
    samples: int,
    *,
    comparison_planets: list[Exoplanet] | None = None,
) -> list[Path]:
    """Write one CSV and graph per planet, plus an optional comparison graph."""
    if samples < 3:
        raise ValueError("--profile-samples must be at least 3.")

    output_directory.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    compared: list[Exoplanet] = []
    comparison_choices: set[int] = set()
    for planet in comparison_planets or []:
        if planet.menu_index not in comparison_choices:
            comparison_choices.add(planet.menu_index)
            compared.append(planet)

    planets_to_generate = list(planets)
    scheduled_choices = {planet.menu_index for planet in planets_to_generate}
    for planet in compared:
        if planet.menu_index not in scheduled_choices:
            scheduled_choices.add(planet.menu_index)
            planets_to_generate.append(planet)

    generated_choices: set[int] = set()
    comparison_profile_cache: dict[int, list[dict[str, object]]] = {}

    for planet in planets_to_generate:
        if planet.menu_index in generated_choices:
            continue
        generated_choices.add(planet.menu_index)
        profile = orbital_temperature_profile(planet, samples)
        stem = profile_file_stem(planet)
        csv_path = output_directory / f"{stem}_temperature_profile.csv"
        graph_path = output_directory / f"{stem}_temperature_profile.png"
        write_temperature_profile(csv_path, profile)
        make_temperature_profile_graph(planet, profile, graph_path)
        generated.extend([csv_path, graph_path])
        if planet.menu_index in comparison_choices:
            comparison_profile_cache[planet.menu_index] = profile

    comparison_profiles = [
        (planet, comparison_profile_cache[planet.menu_index])
        for planet in compared
    ]

    if len(comparison_profiles) >= 2:
        comparison_path = output_directory / comparison_graph_filename(compared)
        make_temperature_comparison_graph(comparison_profiles, comparison_path)
        generated.append(comparison_path)

    return generated

def temperature_to_color(temp_k: float, min_k: float, max_k: float) -> str:
    """
    Blue-white-orange-red gradient for planet bands.
    """
    temperature_span = max(max_k - min_k, 1e-9)
    anchors = [
        (min_k, (180, 220, 255)),
        (min_k + temperature_span * 0.38, (230, 245, 255)),
        (min_k + temperature_span * 0.55, (255, 220, 140)),
        (min_k + temperature_span * 0.75, (255, 145, 70)),
        (max_k, (215, 45, 45)),
    ]
    value = max(min_k, min(max_k, temp_k))
    for (ta, ca), (tb, cb) in zip(anchors, anchors[1:]):
        if ta <= value <= tb:
            factor = 0.0 if tb == ta else (value - ta) / (tb - ta)
            r = round(ca[0] + factor * (cb[0] - ca[0]))
            g = round(ca[1] + factor * (cb[1] - ca[1]))
            b = round(ca[2] + factor * (cb[2] - ca[2]))
            return f"#{r:02x}{g:02x}{b:02x}"
    return "#d72d2d"

# Tkinter animation engine
class ExoplanetOrbitApp:
    def __init__(
        self,
        planets: list[Exoplanet],
        selected: Exoplanet,
        orbit_seconds: float,
        eccentricity_cycle_seconds: float,
        duration_seconds: float,
    ) -> None:
        self.planets = planets
        self.selected = selected
        self.orbit_seconds = max(orbit_seconds, 0.1)
        self.eccentricity_cycle_seconds = max(eccentricity_cycle_seconds, 1.0)
        self.duration_seconds = max(duration_seconds, 0.0)
        self.start_time = time.perf_counter()
        self.paused = False
        self.pause_started_at = 0.0
        self.total_pause_time = 0.0
        self.finished = False

        self.root = tk.Tk()
        self.root.title("Exoplanet orbital-temperature simulation")
        self.root.geometry(f"{WINDOW_WIDTH}x{WINDOW_HEIGHT}")
        self.root.configure(bg="#08111f")

        self.canvas = tk.Canvas(
            self.root,
            width=WINDOW_WIDTH - SIDEBAR_WIDTH,
            height=WINDOW_HEIGHT,
            bg="#050a14",
            highlightthickness=0,
        )
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.sidebar = tk.Frame(self.root, width=SIDEBAR_WIDTH, bg="#0d1424")
        self.sidebar.pack(side=tk.RIGHT, fill=tk.Y)
        self.sidebar.pack_propagate(False)

        self.metrics: dict[str, tk.Label] = {}
        self._build_sidebar()

        self.max_extent_au = max(
            p.semi_major_axis_au * (1.0 + max(p.eccentricity, p.eccentricity_max))
            for p in self.planets
        )
        self.max_extent_au = max(self.max_extent_au, 0.01)

        self.temp_min, self.temp_max = self._estimate_temperature_range()

        self.root.bind("<space>", lambda event: self.toggle_pause())
        self.root.bind("<Escape>", lambda event: self.root.destroy())

    def _build_sidebar(self) -> None:
        title = tk.Label(
            self.sidebar,
            text="Exoplanet Orbit Model",
            fg="white",
            bg="#0d1424",
            font=("Arial", 18, "bold"),
        )
        title.pack(anchor="w", padx=20, pady=(20, 8))

        description = (
            f"{self.selected.host} is the host star. {self.selected.name} is the exoplanet being displayed.\n"
            "If eccentricity limits exist, the eccentricity oscillates between them."
        )
        tk.Label(
            self.sidebar,
            text=description,
            wraplength=360,
            justify="left",
            fg="#d9e8ff",
            bg="#0d1424",
        ).pack(anchor="w", padx=20, pady=(0, 12))

        for key in [
            "Selected planet",
            "Host star",
            "Model orbital period",
            "Source orbital period",
            "Semi-major axis",
            "Eccentricity",
            "Current periapsis",
            "Current apoapsis",
            "Current distance",
            "Orbital speed",
            "Global temp",
            "North pole temp",
            "Equator temp",
            "South pole temp",
        ]:
            self._add_metric_row(key)

        legend_text = (
            "Space: pause/resume\n"
            "Esc: quit\n"
        )
        tk.Label(
            self.sidebar,
            text=legend_text,
            wraplength=360,
            justify="left",
            fg="#9fb5d1",
            bg="#0d1424",
        ).pack(anchor="w", padx=20, pady=16)

    def _add_metric_row(self, key: str) -> None:
        frame = tk.Frame(self.sidebar, bg="#172138")
        frame.pack(fill=tk.X, padx=18, pady=3)
        tk.Label(
            frame,
            text=key,
            fg="white",
            bg="#172138",
            font=("Arial", 9, "bold"),
        ).pack(anchor="w", padx=10, pady=(5, 0))
        value = tk.Label(
            frame,
            text="-",
            fg="#9fd3ff",
            bg="#172138",
            font=("Consolas", 9),
        )
        value.pack(anchor="w", padx=10, pady=(1, 5))
        self.metrics[key] = value

    def _estimate_temperature_range(self) -> tuple[float, float]:
        temps = []
        for planet in self.planets:
            for e in [planet.eccentricity_min, planet.eccentricity, planet.eccentricity_max]:
                peri = periapsis_au(planet, e)
                apo = apoapsis_au(planet, e)
                if peri > 0:
                    temps.append(global_equilibrium_temperature_k(planet, peri))
                if apo > 0:
                    temps.append(global_equilibrium_temperature_k(planet, apo))
        if not temps:
            return 150.0, 400.0
        low = min(temps) - 30.0
        high = max(temps) + 30.0
        if abs(high - low) < 1.0:
            high = low + 1.0
        return low, high

    def toggle_pause(self) -> None:
        if self.finished:
            return
        if not self.paused:
            self.paused = True
            self.pause_started_at = time.perf_counter()
        else:
            self.paused = False
            self.total_pause_time += time.perf_counter() - self.pause_started_at
            self._tick()

    def elapsed_real_seconds(self) -> float:
        if self.paused:
            return self.pause_started_at - self.start_time - self.total_pause_time
        return time.perf_counter() - self.start_time - self.total_pause_time

    def current_eccentricity(self, planet: Exoplanet, elapsed_seconds: float) -> float:
        """
        Oscillates between e_min and e_max if a range exists.

        The oscillation is deliberately exaggerated in time for visual clarity:
        e_min -> e_max -> e_min over eccentricity_cycle_seconds.
        """
        if not planet.has_varying_eccentricity:
            return planet.eccentricity

        midpoint = 0.5 * (planet.eccentricity_min + planet.eccentricity_max)
        amplitude = 0.5 * (planet.eccentricity_max - planet.eccentricity_min)

        # Starts at the lower limit, reaches the upper limit halfway through.
        phase = 2.0 * math.pi * (elapsed_seconds / self.eccentricity_cycle_seconds) - math.pi / 2.0
        return clamp_eccentricity(midpoint + amplitude * math.sin(phase))

    def simulated_time(self, elapsed_seconds: float) -> tuple[float, float]:
        """
        Converts real animation time into simulated days.

        The selected planet completes exactly one orbit every orbit_seconds.
        Other planets in the same system move according to their real period
        relative to the selected planet's simulated days.
        """
        selected_orbit_count = elapsed_seconds / self.orbit_seconds
        sim_days = selected_orbit_count * self.selected.period_days
        return selected_orbit_count, sim_days

    def transform(self, x_au: float, y_au: float) -> tuple[float, float]:
        width = max(self.canvas.winfo_width(), 1)
        height = max(self.canvas.winfo_height(), 1)
        scale = 0.43 * min(width, height) / self.max_extent_au
        return width * 0.5 + x_au * scale, height * 0.5 - y_au * scale

    def orbit_points(self, planet: Exoplanet, eccentricity: float) -> Iterable[tuple[float, float]]:
        for i in range(ORBIT_POINTS + 1):
            nu = 2.0 * math.pi * i / ORBIT_POINTS
            x, y, _ = relative_position_au(planet, nu, eccentricity)
            yield x, y

    def planet_state(self, planet: Exoplanet, sim_days: float, eccentricity: float) -> dict[str, float]:
        return calculate_orbital_state(planet, sim_days, eccentricity)

    def draw_orbit(
        self,
        planet: Exoplanet,
        eccentricity: float,
        colour: str,
        width: int,
        dash: tuple[int, int] | None,
    ) -> None:
        points: list[float] = []
        for x, y in self.orbit_points(planet, eccentricity):
            px, py = self.transform(x, y)
            points.extend([px, py])
        if dash is None:
            self.canvas.create_line(*points, fill=colour, width=width, smooth=True)
        else:
            self.canvas.create_line(*points, fill=colour, dash=dash, width=width, smooth=True)

    def draw_planet_disc(self, planet: Exoplanet, state: dict[str, float], px: float, py: float, selected: bool) -> None:
        radius = int((15 if selected else 8) * PLANET_DISPLAY_SCALE)
        band_step = max(1, int(radius / 8))
        sub_lat = state["sub_lat"]
        global_k = state["global_k"]

        for pixel_y in range(int(py - radius), int(py + radius) + 1, band_step):
            dy = pixel_y - py
            if abs(dy) > radius:
                continue
            chord = math.sqrt(max(radius ** 2 - dy ** 2, 0.0))
            latitude = 90.0 * (dy / radius)
            temp = local_temperature_k(global_k, latitude, sub_lat)
            colour = temperature_to_color(temp, self.temp_min, self.temp_max)
            self.canvas.create_line(px - chord, pixel_y, px + chord, pixel_y, fill=colour, width=band_step)

        outline = "#ffffff" if selected else "#c8d8ed"
        self.canvas.create_oval(px - radius, py - radius, px + radius, py + radius, outline=outline, width=2)

        if selected:
            axis_angle = math.radians(90.0 - planet.axial_tilt_deg)
            dx = math.cos(axis_angle) * radius * 0.95
            dy = -math.sin(axis_angle) * radius * 0.95
            self.canvas.create_line(px - dx, py - dy, px + dx, py + dy, fill="#e8f3ff", width=2)
            self.canvas.create_text(px, py + radius + 13, text=planet.name, fill="#e8f3ff", font=("Arial", 9, "bold"))

    def draw_orbit_limit_labels(self, planet: Exoplanet) -> None:
        if not planet.has_varying_eccentricity:
            return

        min_apo_x, min_apo_y = self.transform(apoapsis_au(planet, planet.eccentricity_min), 0.0)
        max_apo_x, max_apo_y = self.transform(apoapsis_au(planet, planet.eccentricity_max), 0.0)
        self.canvas.create_text(
            min_apo_x + 8,
            min_apo_y - 12,
            text="min e orbit",
            fill="#b9d7ff",
            font=("Arial", 9),
            anchor="w",
        )
        self.canvas.create_text(
            max_apo_x + 8,
            max_apo_y + 12,
            text="max e orbit",
            fill="#ff9f9f",
            font=("Arial", 9),
            anchor="w",
        )

    def draw_frame(self) -> None:
        self.canvas.delete("all")
        elapsed = self.elapsed_real_seconds()
        selected_orbit_count, sim_days = self.simulated_time(elapsed)

        eccentricities = {
            planet.menu_index: self.current_eccentricity(planet, elapsed)
            for planet in self.planets
        }
        states = {
            planet.menu_index: self.planet_state(planet, sim_days, eccentricities[planet.menu_index])
            for planet in self.planets
        }

        # Draw a thin grey dashed semi-major-axis guide.
        selected_e = eccentricities[self.selected.menu_index]
        centre_x = -self.selected.semi_major_axis_au * selected_e
        peri_x = self.selected.semi_major_axis_au * (1.0 - selected_e)
        x1, y1 = self.transform(centre_x, 0.0)
        x2, y2 = self.transform(peri_x, 0.0)
        self.canvas.create_line(x1, y1, x2, y2, fill="#8f98a8", dash=(3, 6), width=1)

        self.draw_orbit(self.selected, selected_e, "#f4f8ff", 3, (8, 6))

        # Draw fixed host star.
        sx, sy = self.transform(0.0, 0.0)
        star_outer_radius = int(18 * STAR_DISPLAY_SCALE)
        star_inner_radius = int(7 * STAR_DISPLAY_SCALE)
        self.canvas.create_oval(
            sx - star_outer_radius,
            sy - star_outer_radius,
            sx + star_outer_radius,
            sy + star_outer_radius,
            fill="#ffca45",
            outline="",
        )
        self.canvas.create_oval(
            sx - star_inner_radius,
            sy - star_inner_radius,
            sx + star_inner_radius,
            sy + star_inner_radius,
            fill="#fff2a6",
            outline="",
        )
        self.canvas.create_text(
            sx,
            sy + star_outer_radius + 13,
            text=self.selected.host,
            fill="#ffe9a8",
            font=("Arial", 10, "bold"),
        )

        for planet in self.planets:
            state = states[planet.menu_index]
            selected = planet.menu_index == self.selected.menu_index
            px, py = self.transform(state["x"], state["y"])
            self.draw_planet_disc(planet, state, px, py, selected)

        selected_state = states[self.selected.menu_index]
        selected_e = eccentricities[self.selected.menu_index]
        self.update_sidebar(selected_state, selected_e)

        self.canvas.create_text(
            22,
            24,
            text="Exoplanet system orbit and temperature model",
            fill="#d6e7ff",
            font=("Arial", 13, "bold"),
            anchor="w",
        )
        self.canvas.create_text(
            22,
            45,
            text="Colours show the educational latitude-temperature approximation",
            fill="#9fb5d1",
            font=("Arial", 10),
            anchor="w",
        )

    def update_sidebar(
        self,
        state: dict[str, float],
        selected_e: float,
    ) -> None:
        self.metrics["Selected planet"].config(text=f"{self.selected.name}  [choice {self.selected.menu_index}]")
        self.metrics["Host star"].config(text=self.selected.host)
        self.metrics["Model orbital period"].config(
            text=f"{self.selected.period_days / DAYS_PER_YEAR:,.6f} years"
        )
        if self.selected.source_period_days is None:
            source_period_text = "not supplied"
        else:
            source_period_text = (
                f"{self.selected.source_period_days / DAYS_PER_YEAR:,.6f} years"
            )
        self.metrics["Source orbital period"].config(text=source_period_text)
        self.metrics["Semi-major axis"].config(text=f"{self.selected.semi_major_axis_au:.6f} AU")
        self.metrics["Eccentricity"].config(text=f"{selected_e:.6f}")
        self.metrics["Current periapsis"].config(text=f"{state['peri_au']:.6f} AU")
        self.metrics["Current apoapsis"].config(text=f"{state['apo_au']:.6f} AU")
        self.metrics["Current distance"].config(text=f"{state['r']:.6f} AU")
        self.metrics["Orbital speed"].config(text=f"{state['speed'] / 1000.0:.3f} km/s")
        self.metrics["Global temp"].config(text=f"{state['global_k']:.2f} K")
        self.metrics["North pole temp"].config(text=f"{state['north_pole_k']:.2f} K")
        self.metrics["Equator temp"].config(text=f"{state['equator_k']:.2f} K")
        self.metrics["South pole temp"].config(text=f"{state['south_pole_k']:.2f} K")

    def _tick(self) -> None:
        if self.paused:
            return

        elapsed = self.elapsed_real_seconds()
        if self.duration_seconds > 0 and elapsed >= self.duration_seconds:
            self.finished = True

        self.draw_frame()

        if self.finished:
            self.root.after(100, self.root.destroy)
        else:
            self.root.after(FPS_MS, self._tick)

    def run(self) -> None:
        self._tick()
        self.root.mainloop()

# Main program
def matching_planets(planets: list[Exoplanet], query: str) -> list[Exoplanet]:
    """Return planets whose name, host, or combined label contains a query."""
    normalised_query = query.strip().casefold()
    if not normalised_query:
        return list(planets)
    return [
        planet
        for planet in planets
        if normalised_query in planet.name.casefold()
        or normalised_query in planet.host.casefold()
        or normalised_query in f"{planet.name} @ {planet.host}".casefold()
    ]


def resolve_planet_selector(planets: list[Exoplanet], selector: str | int) -> Exoplanet:
    """Resolve a 1-based menu number, exact name, or unique name/host search."""
    text = str(selector).strip()
    if not text:
        raise ValueError("Planet selectors cannot be empty.")

    if text.isdecimal():
        choice = int(text)
        for planet in planets:
            if planet.menu_index == choice:
                return planet
        raise ValueError(f"Choice {choice} was not found. Choose a number from 1 to {len(planets)}.")

    folded = text.casefold()
    exact_label_matches = [
        planet
        for planet in planets
        if folded == f"{planet.name} @ {planet.host}".casefold()
    ]
    if len(exact_label_matches) == 1:
        return exact_label_matches[0]

    exact_name_matches = [planet for planet in planets if folded == planet.name.casefold()]
    if len(exact_name_matches) == 1:
        return exact_name_matches[0]

    matches = matching_planets(planets, text)
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValueError(f"No planet or host matched {text!r}.")

    examples = ", ".join(
        f"{planet.menu_index}: {planet.name} @ {planet.host}"
        for planet in matches[:5]
    )
    suffix = " ..." if len(matches) > 5 else ""
    raise ValueError(
        f"Selector {text!r} matched {len(matches)} planets. "
        f"Use a menu number or the exact 'planet @ host' label. Matches: {examples}{suffix}"
    )


def choose_selected_planet(
    planets: list[Exoplanet],
    requested_choice: int | None,
    requested_planet: str | None = None,
) -> Exoplanet:
    """
    Selects exactly one planet.

    Direct selectors are resolved immediately. Otherwise, the terminal accepts a
    global menu number, exact name, or unique planet/host search.
    """
    if requested_choice is not None:
        return resolve_planet_selector(planets, requested_choice)

    if requested_planet is not None:
        return resolve_planet_selector(planets, requested_planet)

    if len(planets) == 1 and planets[0].name == "Earth" and planets[0].host == "Sun":
        print("No exoplanet CSV was supplied, so the program will use Earth orbiting the Sun.")
        return planets[0]

    print_planet_menu(planets)
    while True:
        raw = input(
            "Enter a global choice number, an exact planet name, or search text: "
        ).strip()
        if not raw:
            print("Please enter a choice or search text.")
            continue
        try:
            return resolve_planet_selector(planets, raw)
        except ValueError as error:
            matches = matching_planets(planets, raw)
            print(error)
            if matches:
                print_planet_menu(matches)

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate full-orbit temperature profiles, compare planets, and optionally "
            "animate one selected exoplanet."
        )
    )
    parser.add_argument(
        "csv_file",
        type=Path,
        nargs="?",
        default=None,
        help="Optional CSV file containing exoplanet data. If omitted, Earth/Sun is used.",
    )
    parser.add_argument(
        "--choice",
        type=int,
        default=None,
        help="Planet menu number to choose directly. This is 1 to n, not the CSV row index.",
    )
    parser.add_argument(
        "--planet",
        default=None,
        help="Planet selector: exact name, 'planet @ host', or a unique name/host search.",
    )
    parser.add_argument(
        "--compare",
        action="append",
        default=[],
        metavar="SELECTOR",
        help=(
            "Add a planet to the comparison graph by menu number or name. "
            "Repeat this option to compare several planets."
        ),
    )
    parser.add_argument(
        "--list-planets",
        nargs="?",
        const="",
        default=None,
        metavar="FILTER",
        help="List selectable planets (optionally filtered by planet or host text) and exit.",
    )
    parser.add_argument(
        "--row",
        type=int,
        default=None,
        help="Deprecated alias for --choice. Use --choice instead.",
    )
    parser.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--orbit-seconds",
        type=float,
        default=60.0,
        help="Visual seconds per selected-planet revolution (default: 60; minimum: 0.1).",
    )
    parser.add_argument(
        "--eccentricity-cycle-seconds",
        type=float,
        default=30.0,
        help="Seconds for the visual min-max-min eccentricity cycle (default: 30; minimum: 1).",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=0.0,
        help="Animation duration in real seconds (default: 0, continuous).",
    )
    parser.add_argument("--tilt", type=float, default=None, help="Override axial tilt in degrees for all planets")
    parser.add_argument(
        "--profile-samples",
        type=int,
        default=DEFAULT_PROFILE_SAMPLES,
        help=(
            "Samples per complete orbit, including both endpoints "
            f"(default: {DEFAULT_PROFILE_SAMPLES}; minimum: 3)."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_PROFILE_DIRECTORY,
        help=f"Directory for profile CSVs and graphs (default: {DEFAULT_PROFILE_DIRECTORY}).",
    )
    parser.add_argument(
        "--all-profiles",
        action="store_true",
        help="Generate an individual profile CSV and graph for every loaded, deduplicated planet.",
    )
    parser.add_argument(
        "--no-gui",
        action="store_true",
        help="Generate profile outputs without opening the Tk animation.",
    )
    args = parser.parse_args()

    requested_choice = args.choice if args.choice is not None else args.row

    if args.choice is not None and args.row is not None:
        parser.error("Use either --choice or the deprecated --row alias, not both.")
    if requested_choice is not None and args.planet is not None:
        parser.error("Use either --choice/--row or --planet, not both.")
    if args.profile_samples < 3:
        parser.error("--profile-samples must be at least 3.")
    if not math.isfinite(args.orbit_seconds) or args.orbit_seconds < 0.1:
        parser.error("--orbit-seconds must be a finite number of at least 0.1.")
    if (
        not math.isfinite(args.eccentricity_cycle_seconds)
        or args.eccentricity_cycle_seconds < 1.0
    ):
        parser.error("--eccentricity-cycle-seconds must be a finite number of at least 1.")
    if not math.isfinite(args.duration) or args.duration < 0:
        parser.error("--duration must be a finite, non-negative number.")
    if args.tilt is not None and not math.isfinite(args.tilt):
        parser.error("--tilt must be finite.")

    if args.csv_file is None:
        earth = earth_sun_planet()
        if args.tilt is not None:
            earth.axial_tilt_deg = args.tilt
        planets = [earth]
    else:
        if not args.csv_file.exists():
            parser.error(f"Could not find CSV file: {args.csv_file}")
        else:
            planets = load_planets(args.csv_file, forced_tilt_deg=args.tilt)
            if not planets:
                parser.error(f"No usable exoplanets were loaded from {args.csv_file}.")

    if args.list_planets is not None:
        matches = matching_planets(planets, args.list_planets)
        if not matches:
            print(f"No planets matched {args.list_planets!r}.")
            return
        print_planet_menu(matches, limit=len(matches))
        print(f"Listed {len(matches)} of {len(planets)} loaded planets.")
        return

    if not args.no_gui and tk is None:
        parser.error("Tkinter is unavailable. Install Tk support or rerun with --no-gui.")

    try:
        if args.all_profiles and requested_choice is None and args.planet is None and not args.compare:
            selected = planets[0]
        else:
            selected = choose_selected_planet(planets, requested_choice, args.planet)
        comparison_planets = [selected]
        for selector in args.compare:
            comparison = resolve_planet_selector(planets, selector)
            if comparison.menu_index not in {planet.menu_index for planet in comparison_planets}:
                comparison_planets.append(comparison)
    except ValueError as error:
        parser.error(str(error))

    profile_planets = planets if args.all_profiles else comparison_planets
    print(f"Generating full-orbit outputs for {len(profile_planets)} planet(s)...")
    generated_paths = generate_profile_outputs(
        profile_planets,
        args.output_dir,
        args.profile_samples,
        comparison_planets=comparison_planets,
    )

    print(f"\nLoaded selectable planets: {len(planets)}")
    print(f"Primary planet: {selected.name} around {selected.host}")
    if len(comparison_planets) > 1:
        print("Compared profiles: " + ", ".join(planet.name for planet in comparison_planets))
    print(f"Generated {len(generated_paths)} profile file(s) in {args.output_dir}")
    for path in generated_paths:
        print(f"  {path}")

    if args.no_gui:
        return

    app = ExoplanetOrbitApp(
        planets=[selected],
        selected=selected,
        orbit_seconds=args.orbit_seconds,
        eccentricity_cycle_seconds=args.eccentricity_cycle_seconds,
        duration_seconds=args.duration,
    )
    app.run()


if __name__ == "__main__":
    main()
