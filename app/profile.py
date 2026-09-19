"""Instrument-profile simulation of the neutral isotope stick spectrum.

The neutral formula's isotope envelope (the *stick spectrum*) is broadened
into a continuous instrument profile: every stick at mass m is replaced by
a Gaussian kernel whose full width at half maximum (FWHM) follows the
resolving power R as ``FWHM = m / R``, i.e. ``sigma = FWHM / (2*sqrt(2*ln2))``.
Each kernel is accumulated onto a uniform mass grid within +/-4 sigma of
its centre; a kernel extending beyond the ``[mass_min, mass_max]`` interval
is truncated — its in-interval contribution is kept as-is and the
out-of-interval tail is dropped (no renormalisation). The accumulated
profile is finally normalised to a maximum intensity of 100.

When a range bound is omitted, the outermost +/-4-sigma boundary over all
sticks is used, so the default grid covers every kernel's support. The
grid is capped at a configurable number of points; requests needing more
points are rejected with a locatable error.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .enrichment import normalise_base_peak
from .envelope import Peak, isotope_envelope
from .parser import ParsedFormula

# Gaussian kernels are accumulated within +/- this many sigma of their centre.
KERNEL_WIDTH_SIGMAS = 4.0

# Gaussian sigma per unit FWHM: 1 / (2*sqrt(2*ln 2)).
_FWHM_TO_SIGMA = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))

# Float-jitter tolerance when locating grid indices on exact boundaries.
_EPS = 1e-9


class ProfileError(ValueError):
    """Validation error for profile-simulation parameters, carrying the
    name of the offending request field (``None`` when not field-specific)."""

    def __init__(self, message: str, field: str | None = None) -> None:
        super().__init__(message)
        self.message = message
        self.field = field


@dataclass
class ProfileResult:
    # Grid bounds actually used (defaults applied).
    mass_min: float
    mass_max: float
    # (mass, intensity) samples on the uniform grid, ascending mass.
    points: list[Peak]


def gaussian_sigma(mass: float, resolving_power: float) -> float:
    """Sigma of the Gaussian kernel for a stick at ``mass``: the resolving
    power sets ``FWHM = mass / resolving_power`` and a Gaussian's FWHM is
    ``2*sqrt(2*ln 2)`` sigma."""
    return mass / resolving_power * _FWHM_TO_SIGMA


def instrument_profile(
    sticks: list[Peak],
    resolving_power: float,
    mass_step: float,
    mass_min: float | None,
    mass_max: float | None,
    max_points: int,
) -> ProfileResult:
    """Broaden ``sticks`` into a Gaussian instrument profile on a uniform grid.

    Grid points sit at ``mass_min + k * mass_step``; an omitted bound
    defaults to the corresponding outermost +/-4-sigma boundary over all
    sticks. Raises :class:`ProfileError` for an empty stick list, an
    empty/inverted range, a grid exceeding ``max_points``, or a range no
    kernel contributes to (no valid sampling points).
    """
    if not sticks:
        raise ProfileError("Stick spectrum is empty; nothing to profile")

    kernels = [
        (mass, intensity, gaussian_sigma(mass, resolving_power))
        for mass, intensity in sticks
    ]
    lo = mass_min
    if lo is None:
        lo = min(m - KERNEL_WIDTH_SIGMAS * sigma for m, _, sigma in kernels)
    hi = mass_max
    if hi is None:
        hi = max(m + KERNEL_WIDTH_SIGMAS * sigma for m, _, sigma in kernels)
    if not lo < hi:
        # Only reachable with explicit bounds; the defaults always span at
        # least one kernel's 8-sigma support.
        raise ProfileError(
            f"mass range is empty: mass_min ({lo}) must be smaller than "
            f"mass_max ({hi})",
            "mass_min" if mass_min is not None and mass_max is None else "mass_max",
        )

    n_points = math.floor((hi - lo) / mass_step + _EPS) + 1
    if n_points > max_points:
        raise ProfileError(
            f"mass range {lo}..{hi} with mass_step {mass_step} needs "
            f"{n_points} grid points, exceeding the limit of {max_points}; "
            f"increase mass_step or narrow the range",
            "mass_step",
        )

    grid = [0.0] * n_points
    for mass, intensity, sigma in kernels:
        half = KERNEL_WIDTH_SIGMAS * sigma
        first = max(0, math.ceil((mass - half - lo) / mass_step - _EPS))
        last = min(n_points - 1, math.floor((mass + half - lo) / mass_step + _EPS))
        for k in range(first, last + 1):
            z = (lo + k * mass_step - mass) / sigma
            grid[k] += intensity * math.exp(-0.5 * z * z)

    top = max(grid)
    if top <= 0.0:
        # Only reachable with an explicit range that misses every kernel's
        # +/-4-sigma support; the defaults always cover all of them.
        if hi <= min(m for m, _, _ in kernels):
            field = "mass_max"
        elif lo >= max(m for m, _, _ in kernels):
            field = "mass_min"
        else:
            field = None
        raise ProfileError(
            f"no stick contributes within {lo}..{hi}: the range misses "
            f"every kernel's +/-{KERNEL_WIDTH_SIGMAS:g}-sigma support",
            field,
        )

    scale = 100.0 / top
    points = [(lo + k * mass_step, value * scale) for k, value in enumerate(grid)]
    return ProfileResult(lo, hi, points)


def simulate_profile(
    parsed: ParsedFormula,
    resolving_power: float,
    mass_step: float,
    mass_min: float | None,
    mass_max: float | None,
    merge_tolerance: float,
    prune: float,
    max_points: int,
) -> ProfileResult:
    """Instrument profile of the neutral isotope envelope of ``parsed``.

    The envelope is convolved with the existing machinery and normalised to
    base peak 100 (the conventional stick-spectrum scale), then broadened
    by :func:`instrument_profile`.
    """
    envelope = isotope_envelope(parsed, merge_tolerance, prune)
    sticks = normalise_base_peak(envelope)
    return instrument_profile(
        sticks, resolving_power, mass_step, mass_min, mass_max, max_points
    )
