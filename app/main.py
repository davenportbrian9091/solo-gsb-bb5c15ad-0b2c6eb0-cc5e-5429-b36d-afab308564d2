"""FastAPI application: mass-spectrometry isotope-envelope precompute API."""

from __future__ import annotations

import math
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field, field_validator

from . import __version__, db
from .adduct import AdductError, apply_adduct, parse_adduct
from .config import settings
from .enrichment import (
    LabelCombinationError,
    PeakListError,
    fit_dual_enrichment,
    fit_enrichment,
    normalise_base_peak,
    validate_observed_peaks,
)
from .envelope import (
    ELECTRON_MASS_DA,
    ion_envelope,
    isotope_envelope,
    monoisotopic_mass,
    top_peaks,
)
from .identify import IonCandidatePrediction, predict_ion_candidate, score_ion_candidate
from .isotopes import ISOTOPES
from .mixture import fit_ion_mixture, score_mixture
from .parser import FormulaError, LabelError, parse_formula, parse_label
from .profile import ProfileParameterError, simulate_profile


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Dependencies are installed at image build time; the SQLite schema is
    # created here so a fresh container is ready to serve immediately.
    db.init_db(settings.db_path)
    yield
    db.close_db()


app = FastAPI(
    title="MS Precompute API",
    version=__version__,
    description=(
        "Computes natural-abundance isotope envelopes (mass/intensity "
        "polynomial convolution) and monoisotopic masses for molecular "
        "formulas, with optional pure-isotope labels such as [13C]; also "
        "estimates single-isotope enrichment by fitting candidate "
        "abundances against observed peaks, jointly fits two isotopes of "
        "different elements on a 2-D abundance grid, fits pairwise "
        "0-100% mixtures of candidate adduct/charge ion envelopes, and "
        "simulates Gaussian instrument profiles of neutral stick spectra."
    ),
    lifespan=lifespan,
)


class EnvelopeRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    top_n: int = Field(
        settings.top_n_default,
        ge=1,
        le=settings.top_n_max,
        description="Number of most intense peaks to return",
    )


class PeakOut(BaseModel):
    mass: float
    intensity: float


class EnvelopeResponse(BaseModel):
    formula: str
    composition: dict[str, int]
    labels: dict[str, int]
    monoisotopic_mass: float
    base_peak_intensity: float
    peak_count: int
    peaks: list[PeakOut]
    cached: bool


class IonEnvelopeRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Neutral molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    adduct: str = Field(
        ...,
        description=(
            "Adduct made of consecutive '+/-[coefficient]formula' fragments; "
            "the coefficient multiplies the whole fragment, e.g. '+H', "
            "'+2H-H2O' or '+Na-H'"
        ),
        examples=["+H"],
    )
    charge: int = Field(
        ...,
        ge=-10,
        le=10,
        description="Signed ion charge; must be -10..-1 or 1..10 (zero is not allowed)",
        examples=[1],
    )
    top_n: int = Field(
        settings.top_n_default,
        ge=1,
        le=settings.top_n_max,
        description="Number of most intense peaks to return",
    )

    @field_validator("charge")
    @classmethod
    def charge_must_be_nonzero(cls, value: int) -> int:
        if value == 0:
            raise ValueError("charge must be non-zero: use -10..-1 or 1..10")
        return value


class MzPeakOut(BaseModel):
    mz: float
    intensity: float


class IonEnvelopeResponse(BaseModel):
    formula: str
    adduct: str
    charge: int
    composition: dict[str, int]
    labels: dict[str, int]
    monoisotopic_mz: float
    base_peak_intensity: float
    peak_count: int
    peaks: list[MzPeakOut]


class NeutralProfileRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Neutral molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    resolving_power: float = Field(
        ...,
        gt=0,
        description=(
            "Instrument resolving power R = m/FWHM; each stick at mass m is "
            "broadened with FWHM = m/R"
        ),
        examples=[100000],
    )
    mass_step: float = Field(
        ...,
        gt=0,
        description="Uniform output-grid spacing in Da (must be > 0)",
        examples=[0.001],
    )
    mass_min: float | None = Field(
        None,
        gt=0,
        description=(
            "Inclusive lower grid bound in Da; defaults to the smallest "
            "mass covered by any stick's +/-4 sigma kernel"
        ),
    )
    mass_max: float | None = Field(
        None,
        gt=0,
        description=(
            "Inclusive upper grid bound in Da; defaults to the largest mass "
            "covered by any stick's +/-4 sigma kernel"
        ),
    )

    @field_validator("resolving_power", "mass_step", "mass_min", "mass_max")
    @classmethod
    def values_must_be_finite(cls, value: float | None) -> float | None:
        if value is not None and not math.isfinite(value):
            raise ValueError("value must be a finite number")
        return value


class NeutralProfilePointOut(BaseModel):
    mass: float
    intensity: float


class NeutralProfileResponse(BaseModel):
    formula: str
    resolving_power: float
    mass_min: float
    mass_max: float
    mass_step: float
    point_count: int
    profile: list[NeutralProfilePointOut]


class ObservedPeak(BaseModel):
    mass: float = Field(..., description="Observed peak mass (Da)")
    intensity: float = Field(..., ge=0, description="Observed peak intensity")


class EnrichmentFitRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    label: str = Field(
        ...,
        description="Single isotope label to estimate, e.g. '13C'",
        examples=["13C"],
    )
    peaks: list[ObservedPeak] = Field(
        ...,
        max_length=200,
        description="Observed peaks, sorted by strictly increasing mass",
    )
    step_percent: float = Field(
        1.0,
        ge=0.5,
        le=10.0,
        description="Abundance grid step in percentage points (0.5-10)",
    )
    mass_tolerance: float = Field(
        ...,
        gt=0,
        description="Mass tolerance (Da) for pairing observed/predicted peaks",
    )


class PeakPairOut(BaseModel):
    observed_mass: float | None
    observed_intensity: float | None
    predicted_mass: float | None
    predicted_intensity: float | None


class EnrichmentFitResponse(BaseModel):
    formula: str
    label: str
    best_abundance_percent: float
    mse: float
    search_points: int
    pair_count: int
    pairs: list[PeakPairOut]


class DualLabelSpec(BaseModel):
    label: str = Field(
        ...,
        description="Isotope label to estimate, e.g. '13C'",
        examples=["13C"],
    )
    step_percent: float = Field(
        1.0,
        ge=0.5,
        le=10.0,
        description="Abundance grid step in percentage points (0.5-10)",
    )


class DualEnrichmentFitRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    labels: list[DualLabelSpec] = Field(
        ...,
        min_length=2,
        max_length=2,
        description="Exactly two isotope labels, which must target different elements",
        examples=[[{"label": "13C", "step_percent": 1.0}, {"label": "15N", "step_percent": 1.0}]],
    )
    peaks: list[ObservedPeak] = Field(
        ...,
        max_length=200,
        description="Observed peaks, sorted by strictly increasing mass",
    )
    mass_tolerance: float = Field(
        ...,
        gt=0,
        description="Mass tolerance (Da) for pairing observed/predicted peaks",
    )


class BestAbundanceOut(BaseModel):
    label: str
    best_abundance_percent: float


class DualEnrichmentFitResponse(BaseModel):
    formula: str
    best_abundances: list[BestAbundanceOut]
    mse: float
    search_points: int
    pair_count: int
    pairs: list[PeakPairOut]


class ObservedMzPeak(BaseModel):
    mz: float = Field(..., description="Observed mass-to-charge ratio")
    intensity: float = Field(..., ge=0, description="Observed peak intensity")


class IonCandidateSpec(BaseModel):
    adduct: str = Field(
        ...,
        description=(
            "Adduct made of consecutive '+/-[coefficient]formula' fragments, "
            "e.g. '+H', '+2H-H2O' or '+Na-H'"
        ),
        examples=["+H"],
    )
    charge: int = Field(
        ...,
        ge=-10,
        le=10,
        description="Signed ion charge; must be -10..-1 or 1..10 (zero is not allowed)",
        examples=[1],
    )

    @field_validator("charge")
    @classmethod
    def charge_must_be_nonzero(cls, value: int) -> int:
        if value == 0:
            raise ValueError("charge must be non-zero: use -10..-1 or 1..10")
        return value


class IonCandidatePairOut(BaseModel):
    observed_mz: float | None
    observed_intensity: float | None
    predicted_mz: float | None
    predicted_intensity: float | None


class IonCandidateOut(BaseModel):
    index: int = Field(..., description="0-based position of the candidate in the request")
    valid: bool
    adduct: str
    charge: int
    mse: float | None = Field(None, description="Normalised intensity MSE (null for invalid candidates)")
    monoisotopic_mz: float | None = None
    composition: dict[str, int] | None = None
    labels: dict[str, int] | None = None
    pair_count: int | None = None
    pairs: list[IonCandidatePairOut] | None = None
    error: str | None = Field(
        None, description="Machine-readable error code for invalid candidates"
    )
    message: str | None = Field(
        None, description="Human-readable reason why the candidate is invalid"
    )
    position: int | None = Field(
        None, description="0-based character position in the adduct string"
    )


class IonIdentifyBestOut(BaseModel):
    index: int = Field(..., description="0-based index of the best candidate in the request")
    adduct: str
    charge: int
    mse: float


class IonIdentifyResponse(BaseModel):
    formula: str
    mass_tolerance: float
    observed_peak_count: int
    candidates: list[IonCandidateOut]
    best_candidate: IonIdentifyBestOut | None = None


class IonIdentifyRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Neutral molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    peaks: list[ObservedMzPeak] = Field(
        ...,
        min_length=1,
        max_length=200,
        description="Observed m/z peaks, 1-200 entries, sorted by strictly increasing m/z",
    )
    mass_tolerance: float = Field(
        ...,
        gt=0,
        description="m/z tolerance for pairing observed/predicted peaks",
    )
    candidates: list[IonCandidateSpec] = Field(
        ...,
        min_length=1,
        max_length=20,
        description="1-20 candidate adduct/charge combinations to score",
    )


class IonMixtureFitRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Neutral molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    peaks: list[ObservedMzPeak] = Field(
        ...,
        min_length=1,
        max_length=200,
        description="Observed m/z peaks, 1-200 entries, sorted by strictly increasing m/z",
    )
    mass_tolerance: float = Field(
        ...,
        gt=0,
        description=(
            "m/z tolerance for fusing overlaid mixture peaks and for "
            "pairing observed/predicted peaks"
        ),
    )
    candidates: list[IonCandidateSpec] = Field(
        ...,
        min_length=2,
        max_length=10,
        description="2-10 candidate adduct/charge combinations to mix pairwise",
    )
    step_percent: float = Field(
        1.0,
        ge=0.5,
        le=10.0,
        description="Mixture-ratio grid step in percentage points (0.5-10)",
    )


class IonMixtureCandidateOut(BaseModel):
    index: int = Field(..., description="0-based position of the candidate in the request")
    valid: bool
    adduct: str
    charge: int
    monoisotopic_mz: float | None = None
    error: str | None = Field(
        None, description="Machine-readable error code for invalid candidates"
    )
    message: str | None = Field(
        None, description="Human-readable reason why the candidate is invalid"
    )
    position: int | None = Field(
        None, description="0-based character position in the adduct string"
    )


class IonMixtureCombinationOut(BaseModel):
    rank: int = Field(..., description="1-based rank by ascending MSE")
    first_index: int = Field(..., description="0-based index of the lower-index candidate")
    second_index: int = Field(..., description="0-based index of the higher-index candidate")
    first_adduct: str
    first_charge: int
    second_adduct: str
    second_charge: int
    first_fraction_percent: float = Field(
        ..., description="Mixture ratio of the first candidate, 0-100%"
    )
    second_fraction_percent: float = Field(
        ..., description="Mixture ratio of the second candidate, 0-100%"
    )
    mse: float


class IonMixtureBestOut(BaseModel):
    first_index: int
    second_index: int
    first_adduct: str
    first_charge: int
    second_adduct: str
    second_charge: int
    first_fraction_percent: float
    second_fraction_percent: float
    mse: float
    pair_count: int
    pairs: list[IonCandidatePairOut]


class IonMixtureFitResponse(BaseModel):
    formula: str
    mass_tolerance: float
    step_percent: float
    observed_peak_count: int
    valid_candidate_count: int
    combination_count: int
    candidates: list[IonMixtureCandidateOut]
    combinations: list[IonMixtureCombinationOut]
    best_combination: IonMixtureBestOut | None = None


@app.exception_handler(FormulaError)
async def formula_error_handler(request: Request, exc: FormulaError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": "invalid_formula",
                "message": exc.message,
                "position": exc.position,
            }
        },
    )


@app.exception_handler(AdductError)
async def adduct_error_handler(request: Request, exc: AdductError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": "invalid_adduct",
                "message": exc.message,
                "position": exc.position,
            }
        },
    )


@app.exception_handler(LabelError)
async def label_error_handler(request: Request, exc: LabelError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": "invalid_label",
                "message": exc.message,
                "position": exc.position,
            }
        },
    )


@app.exception_handler(PeakListError)
async def peak_list_error_handler(request: Request, exc: PeakListError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": "invalid_peaks",
                "message": exc.message,
                "index": exc.index,
            }
        },
    )


@app.exception_handler(LabelCombinationError)
async def label_combination_error_handler(
    request: Request, exc: LabelCombinationError
) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": "invalid_labels",
                "message": exc.message,
                "index": exc.index,
            }
        },
    )


@app.exception_handler(ProfileParameterError)
async def profile_parameter_error_handler(
    request: Request, exc: ProfileParameterError
) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": exc.error,
                "message": exc.message,
                "field": exc.field,
            }
        },
    )


@app.get("/health")
def health() -> dict:
    return {"status": "ok", "version": __version__}


@app.get("/api/v1/elements")
def elements() -> dict:
    """List supported element symbols and their natural isotopes."""
    return {
        symbol: [{"mass": m, "abundance": a} for m, a in isotopes]
        for symbol, isotopes in sorted(ISOTOPES.items())
    }


@app.post("/api/v1/envelope", response_model=EnvelopeResponse)
def envelope(req: EnvelopeRequest) -> dict:
    """Compute the isotope envelope and monoisotopic mass of a formula."""
    try:
        parsed = parse_formula(req.formula)
    except FormulaError as exc:
        db.log_request(req.formula, req.top_n, ok=False, message=exc.message)
        raise

    cached = db.get_cached(req.formula, req.top_n)
    if cached is not None:
        db.log_request(req.formula, req.top_n, ok=True, message="cache hit")
        return {**cached, "cached": True}

    env = isotope_envelope(parsed, settings.merge_tolerance, settings.prune_threshold)
    peaks = top_peaks(env, req.top_n)
    result = {
        "formula": req.formula,
        "composition": parsed.composition,
        "labels": {f"[{mass_no}{symbol}]": count for (symbol, mass_no), count in parsed.labels.items()},
        "monoisotopic_mass": round(monoisotopic_mass(parsed), 6),
        "base_peak_intensity": 100.0,
        "peak_count": len(peaks),
        "peaks": [
            {"mass": round(m, 6), "intensity": round(i, 6)} for m, i in peaks
        ],
        "cached": False,
    }
    db.put_cache(req.formula, req.top_n, {k: v for k, v in result.items() if k != "cached"})
    db.log_request(req.formula, req.top_n, ok=True, message="computed")
    return result


@app.post("/api/v1/ion-envelope", response_model=IonEnvelopeResponse)
def ion_envelope_endpoint(req: IonEnvelopeRequest) -> dict:
    """Compute the adducted ion isotope envelope as m/z peaks.

    The neutral formula is parsed, the adduct fragments are applied to the
    atom counts, and the envelope is convolved with the neutral peak merge
    tolerance scaled by |charge|; ionic masses are then converted to m/z via
    ``(mass - charge * m_e) / |charge|``.
    """
    try:
        parsed = parse_formula(req.formula)
        pieces = parse_adduct(req.adduct)
        adjusted = apply_adduct(parsed, pieces)
    except (FormulaError, AdductError) as exc:
        db.log_request(req.formula, req.top_n, ok=False, message=exc.message)
        raise

    # Larger charge states spread isotope peaks more tightly in m/z, so the
    # neutral-mass merge tolerance scales with |charge| before conversion.
    merge_tolerance = abs(req.charge) * settings.merge_tolerance
    env = isotope_envelope(adjusted, merge_tolerance, settings.prune_threshold)
    peaks = top_peaks(env, req.top_n)
    mz_peaks = ion_envelope(peaks, req.charge)
    mono_mz = (
        monoisotopic_mass(adjusted) - req.charge * ELECTRON_MASS_DA
    ) / abs(req.charge)

    result = {
        "formula": req.formula,
        "adduct": req.adduct,
        "charge": req.charge,
        "composition": adjusted.composition,
        "labels": {
            f"[{mass_no}{symbol}]": count
            for (symbol, mass_no), count in adjusted.labels.items()
        },
        "monoisotopic_mz": round(mono_mz, 6),
        "base_peak_intensity": 100.0,
        "peak_count": len(mz_peaks),
        "peaks": [
            {"mz": round(mz, 6), "intensity": round(i, 6)} for mz, i in mz_peaks
        ],
    }
    db.log_request(req.formula, req.top_n, ok=True, message="ion envelope computed")
    return result


@app.post("/api/v1/neutral-profile", response_model=NeutralProfileResponse)
def neutral_profile(req: NeutralProfileRequest) -> dict:
    """Simulate the instrument profile of a neutral-molecule stick envelope.

    The formula is parsed and convolved with the existing natural-abundance
    envelope machinery (no ``top_n`` truncation), then every stick at mass m
    is broadened into a Gaussian with FWHM = m/resolving_power and sampled on
    one uniform grid (spacing ``mass_step``) within +-4 sigma; kernels are
    clipped to the requested window but keep their in-window contribution.
    The summed profile is normalised to a peak intensity of 100.
    """
    try:
        parsed = parse_formula(req.formula)
    except FormulaError as exc:
        db.log_request(req.formula, None, ok=False, message=exc.message)
        raise

    env = isotope_envelope(parsed, settings.merge_tolerance, settings.prune_threshold)
    result = simulate_profile(
        env,
        req.resolving_power,
        req.mass_step,
        req.mass_min,
        req.mass_max,
    )
    db.log_request(req.formula, None, ok=True, message="neutral profile computed")
    return {
        "formula": req.formula,
        "resolving_power": req.resolving_power,
        "mass_min": result.mass_min,
        "mass_max": result.mass_max,
        "mass_step": result.mass_step,
        "point_count": result.point_count,
        "profile": [
            {"mass": mass, "intensity": round(intensity, 6)}
            for mass, intensity in result.profile
        ],
    }


def _build_pairs(
    observed: list[tuple[float, float]],
    fit_pairs: list[tuple[int | None, int | None]],
    predicted: list[tuple[float, float]],
) -> list[dict]:
    """Serialise the best-fit pairing, sorted by mass (observed first)."""
    obs = normalise_base_peak(observed)
    pairs = [
        {
            "observed_mass": round(obs[i][0], 6) if i is not None else None,
            "observed_intensity": round(obs[i][1], 6) if i is not None else None,
            "predicted_mass": round(predicted[j][0], 6) if j is not None else None,
            "predicted_intensity": round(predicted[j][1], 6) if j is not None else None,
        }
        for i, j in fit_pairs
    ]
    pairs.sort(
        key=lambda p: p["observed_mass"]
        if p["observed_mass"] is not None
        else p["predicted_mass"]
    )
    return pairs


def _build_mz_pairs(
    observed: list[tuple[float, float]],
    fit_pairs: list[tuple[int | None, int | None]],
    predicted: list[tuple[float, float]],
) -> list[dict]:
    """Serialise an ion-identification pairing on the m/z axis, sorted by
    m/z (observed m/z first for unmatched-predicted placeholders)."""
    pairs = [
        {
            "observed_mz": round(observed[i][0], 6) if i is not None else None,
            "observed_intensity": round(observed[i][1], 6) if i is not None else None,
            "predicted_mz": round(predicted[j][0], 6) if j is not None else None,
            "predicted_intensity": round(predicted[j][1], 6) if j is not None else None,
        }
        for i, j in fit_pairs
    ]
    pairs.sort(
        key=lambda p: p["observed_mz"]
        if p["observed_mz"] is not None
        else p["predicted_mz"]
    )
    return pairs


@app.post("/api/v1/enrichment-fit", response_model=EnrichmentFitResponse)
def enrichment_fit(req: EnrichmentFitRequest) -> dict:
    """Estimate the abundance of a single isotope label by fitting the
    computed envelope against observed peaks over a 0-100% grid."""
    parsed = parse_formula(req.formula)
    symbol, mass_number = parse_label(req.label)

    observed = [(p.mass, p.intensity) for p in req.peaks]
    validate_observed_peaks(observed)

    fit = fit_enrichment(
        parsed,
        symbol,
        mass_number,
        observed,
        step_percent=req.step_percent,
        mass_tolerance=req.mass_tolerance,
        merge_tolerance=settings.merge_tolerance,
        prune=settings.prune_threshold,
    )

    pairs = _build_pairs(observed, fit.pairs, fit.predicted)
    return {
        "formula": req.formula,
        "label": req.label,
        "best_abundance_percent": round(fit.best_abundance_percent, 6),
        "mse": fit.mse,
        "search_points": fit.search_points,
        "pair_count": len(pairs),
        "pairs": pairs,
    }


@app.post(
    "/api/v1/dual-enrichment-fit", response_model=DualEnrichmentFitResponse
)
def dual_enrichment_fit(req: DualEnrichmentFitRequest) -> dict:
    """Jointly estimate the abundances of two isotope labels (on different
    elements) by fitting the computed envelope against observed peaks over
    a 2-D 0-100% grid."""
    parsed = parse_formula(req.formula)
    parsed_labels = [parse_label(spec.label) for spec in req.labels]

    observed = [(p.mass, p.intensity) for p in req.peaks]
    validate_observed_peaks(observed)

    # Labels targeting the same element are rejected here with a locatable
    # index; parse_label above already covered malformed/unknown labels.
    if parsed_labels[0][0] == parsed_labels[1][0]:
        raise LabelCombinationError(
            f"The two labels must target different elements, but "
            f"'{req.labels[1].label}' targets element "
            f"'{parsed_labels[1][0]}' already named by "
            f"'{req.labels[0].label}'",
            1,
        )

    fit = fit_dual_enrichment(
        parsed,
        parsed_labels[0],
        parsed_labels[1],
        observed,
        step_percent=(req.labels[0].step_percent, req.labels[1].step_percent),
        mass_tolerance=req.mass_tolerance,
        merge_tolerance=settings.merge_tolerance,
        prune=settings.prune_threshold,
    )

    pairs = _build_pairs(observed, fit.pairs, fit.predicted)
    return {
        "formula": req.formula,
        "best_abundances": [
            {"label": req.labels[k].label, "best_abundance_percent": round(pct, 6)}
            for k, pct in enumerate(fit.best_abundances_percent)
        ],
        "mse": fit.mse,
        "search_points": fit.search_points,
        "pair_count": len(pairs),
        "pairs": pairs,
    }


@app.post("/api/v1/ion-identify", response_model=IonIdentifyResponse)
def ion_identify(req: IonIdentifyRequest) -> dict:
    """Score candidate adduct/charge combinations against observed m/z peaks.

    Each candidate reuses the adduct parsing and ionic envelope computation
    of ``/api/v1/ion-envelope``; predicted and observed (base-peak
    normalised) peaks are paired one-to-one, smallest within-tolerance mass
    difference first, and unmatched peaks count at zero intensity in the
    normalised MSE. Invalid candidates are reported in place with their
    request index and the same position-aware adduct error; valid candidates
    are returned sorted by MSE with input-order tie-break.
    """
    parsed = parse_formula(req.formula)

    observed = [(p.mz, p.intensity) for p in req.peaks]
    validate_observed_peaks(observed)
    obs = normalise_base_peak(observed)

    valid_results: list[dict] = []
    invalid_results: list[dict] = []
    for idx, spec in enumerate(req.candidates):
        try:
            prediction = predict_ion_candidate(
                parsed,
                spec.adduct,
                spec.charge,
                settings.merge_tolerance,
                settings.prune_threshold,
            )
        except AdductError as exc:
            # Bad candidates are part of the 200 response (with their index
            # and the existing position-aware error), not a request failure.
            invalid_results.append(
                {
                    "index": idx,
                    "valid": False,
                    "adduct": spec.adduct,
                    "charge": spec.charge,
                    "mse": None,
                    "error": "invalid_adduct",
                    "message": exc.message,
                    "position": exc.position,
                }
            )
            continue

        score = score_ion_candidate(idx, prediction, obs, req.mass_tolerance)
        pairs = _build_mz_pairs(obs, score.pairs, score.predicted)
        valid_results.append(
            {
                "index": idx,
                "valid": True,
                "adduct": spec.adduct,
                "charge": spec.charge,
                "mse": score.mse,
                "monoisotopic_mz": round(score.monoisotopic_mz, 6),
                "composition": score.composition,
                "labels": score.labels,
                "pair_count": len(pairs),
                "pairs": pairs,
                "error": None,
                "message": None,
                "position": None,
            }
        )

    # Stable ordering: ascending MSE; an MSE tie keeps the request order.
    # Invalid candidates follow the scored ones, also in request order.
    valid_results.sort(key=lambda c: (c["mse"], c["index"]))
    candidates = valid_results + invalid_results

    best = None
    if valid_results:
        winner = valid_results[0]
        best = {
            "index": winner["index"],
            "adduct": winner["adduct"],
            "charge": winner["charge"],
            "mse": winner["mse"],
        }

    db.log_request(
        req.formula,
        None,
        ok=bool(valid_results),
        message=(
            f"ion identify: {len(valid_results)} valid, "
            f"{len(invalid_results)} invalid candidates"
        ),
    )
    return {
        "formula": req.formula,
        "mass_tolerance": req.mass_tolerance,
        "observed_peak_count": len(observed),
        "candidates": candidates,
        "best_candidate": best,
    }


@app.post("/api/v1/ion-mixture-fit", response_model=IonMixtureFitResponse)
def ion_mixture_fit(req: IonMixtureFitRequest) -> dict:
    """Fit pairwise mixtures of candidate ion envelopes to observed m/z peaks.

    Each candidate reuses the adduct parsing and ionic envelope computation
    of ``/api/v1/ion-envelope`` (base-peak normalised). For every unordered
    pair of valid candidates — the lower request index acting as the first
    candidate — and every first-candidate ratio on a 0-100% grid
    (``step_percent``), the two envelopes are scaled by their ratio,
    overlaid (peaks within ``mass_tolerance`` fused) and re-normalised to
    base peak 100, then paired one-to-one with the observed peaks; the MSE
    over all pairs, unmatched peaks on both sides included, ranks the
    combinations. Ties resolve to the smaller candidate index pair, then
    the lower first-candidate ratio. Invalid candidates are reported in
    place with the same position-aware error as ``/api/v1/ion-identify``
    and their combinations are skipped; with fewer than two valid
    candidates the combination list is empty and the best combination null.
    """
    parsed = parse_formula(req.formula)

    observed = [(p.mz, p.intensity) for p in req.peaks]
    validate_observed_peaks(observed)
    obs = normalise_base_peak(observed)

    candidates: list[dict] = []
    valid: list[tuple[int, IonCandidatePrediction]] = []
    for idx, spec in enumerate(req.candidates):
        try:
            prediction = predict_ion_candidate(
                parsed,
                spec.adduct,
                spec.charge,
                settings.merge_tolerance,
                settings.prune_threshold,
            )
        except AdductError as exc:
            # Bad candidates are part of the 200 response (with their index
            # and the existing position-aware error); every combination
            # involving them is skipped.
            candidates.append(
                {
                    "index": idx,
                    "valid": False,
                    "adduct": spec.adduct,
                    "charge": spec.charge,
                    "monoisotopic_mz": None,
                    "error": "invalid_adduct",
                    "message": exc.message,
                    "position": exc.position,
                }
            )
            continue
        valid.append((idx, prediction))
        candidates.append(
            {
                "index": idx,
                "valid": True,
                "adduct": spec.adduct,
                "charge": spec.charge,
                "monoisotopic_mz": round(prediction.monoisotopic_mz, 6),
                "error": None,
                "message": None,
                "position": None,
            }
        )

    combinations: list[dict] = []
    best: dict | None = None
    if len(valid) >= 2:
        scored = fit_ion_mixture(valid, obs, req.step_percent, req.mass_tolerance)
        by_index = dict(valid)
        for rank, combo in enumerate(scored, start=1):
            first = by_index[combo.first_index]
            second = by_index[combo.second_index]
            combinations.append(
                {
                    "rank": rank,
                    "first_index": combo.first_index,
                    "second_index": combo.second_index,
                    "first_adduct": first.adduct,
                    "first_charge": first.charge,
                    "second_adduct": second.adduct,
                    "second_charge": second.charge,
                    "first_fraction_percent": round(combo.first_fraction_percent, 6),
                    "second_fraction_percent": round(
                        100.0 - combo.first_fraction_percent, 6
                    ),
                    "mse": combo.mse,
                }
            )

        winner = scored[0]
        first = by_index[winner.first_index]
        second = by_index[winner.second_index]
        # Recompute the winning mixture once to expose its pairing detail.
        fit_pairs, _, mixed = score_mixture(
            first.peaks,
            second.peaks,
            winner.first_fraction_percent / 100.0,
            obs,
            req.mass_tolerance,
        )
        pairs = _build_mz_pairs(obs, fit_pairs, mixed)
        best = {
            "first_index": winner.first_index,
            "second_index": winner.second_index,
            "first_adduct": first.adduct,
            "first_charge": first.charge,
            "second_adduct": second.adduct,
            "second_charge": second.charge,
            "first_fraction_percent": round(winner.first_fraction_percent, 6),
            "second_fraction_percent": round(
                100.0 - winner.first_fraction_percent, 6
            ),
            "mse": winner.mse,
            "pair_count": len(pairs),
            "pairs": pairs,
        }

    db.log_request(
        req.formula,
        None,
        ok=best is not None,
        message=(
            f"ion mixture fit: {len(valid)} valid candidates, "
            f"{len(combinations)} combinations"
        ),
    )
    return {
        "formula": req.formula,
        "mass_tolerance": req.mass_tolerance,
        "step_percent": req.step_percent,
        "observed_peak_count": len(observed),
        "valid_candidate_count": len(valid),
        "combination_count": len(combinations),
        "candidates": candidates,
        "combinations": combinations,
        "best_combination": best,
    }
