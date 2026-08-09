# The Wavenumber-1 Eccentric Shadow: A Structural Precursor for Tropical Cyclone Track Direction

**Author:** tygtDc (Deep Research)
**Affiliation:** Independent Research
**Date:** 2026-08-09
**Version:** 1.0.0

---

## Abstract

Tropical cyclone (TC) track forecasting relies on environmental steering, yet the *internal structural asymmetry* of the storm itself is rarely used as a predictive signal. We introduce the **Wavenumber-1 (WN1) Eccentric Shadow** operator: the phase of the azimuthal wavenumber-1 component of the radial wind in the inner ring (r = 1–5° from the storm center) at 500 hPa. We show that this phase points toward the future direction of TC motion with a mean angular error of **23.9°** over 76 pseudo-real-time samples from 6 typhoons (random baseline: 90°; skill ≈ 73%), evaluated at a 12-hour forward lead. With an uncertainty (UQ) filter on structural purity (ellipticity ellipt ≤ 0.4), the error drops to **21.1°** with **90–91%** of samples within 45°. The signal degrades slowly with lead time: 20.1° at 6 h to 26.6° at 48 h (≈ 0.16°/h), indicating a usable forecast window of at least 48 h, not merely 12 h. We compare three pressure levels: 850 hPa (52.0°), 200 hPa (31.8°), and 500 hPa (23.9°); 500 hPa — the classic steering level — is optimal, while 850 hPa is contaminated by convective coupling near the inner core. The operator responds to the three standard critiques of structural track signals (simultaneity, downwind sampling, vertical wind shear) with data, and a phase-jump diagnostic (Δφ) correlates with observed turning (ρ = 0.589, n = 69). The method is deployed as a fully automated, open-source GitHub Actions pipeline that measures WN1 phase every 6 h for active typhoons using real-time GFS analyses (NOMADS). The WN1 shadow is proposed as a *complementary* second-opinion signal for operational forecasting — particularly valuable as an early turning alarm — not as a standalone replacement for dynamical models.

---

## 1. Introduction

Tropical cyclone motion is governed to first order by the environmental steering flow (Chan & Gray, 1982). Operational track models assimilate four-dimensional observations into full-physics dynamical systems with ensemble members, achieving 24 h track errors of roughly 100–150 km. Yet two intrinsic weaknesses remain: (i) the prediction of **recurvature/turning**, where the subtropical high evolves and steering changes, and (ii) episodes of **rapid structural change** (rapid intensification, RI; eyewall replacement cycles, ERC) during which internal dynamics temporarily dominate over environmental advection.

The guiding idea of this work is that the *asymmetric structure of the TC itself* — not just the environment around it — encodes information about its imminent motion. A mature TC is not a symmetric vortex: its inner-core convection is organized into a wavenumber-1 (WN1) asymmetry whose orientation is physically linked to the vortex's translation and to the deep-layer steering. If this asymmetry can be measured cleanly, its azimuthal phase may serve as a **structural shadow** of the storm's future heading.

We define the WN1 Shadow operator as the phase of the azimuthal wavenumber-1 Fourier mode of the radial wind on a ring 1–5° from the center at 500 hPa. We then:

1. Validate its *forward* (not merely simultaneous) skill on 76 samples from 6 western North Pacific typhoons using ERA5 reanalyses;
2. Quantify the uncertainty-controlled operating envelope via a structural-purity filter;
3. Establish the lead-time degradation curve to define the usable forecast window;
4. Compare pressure levels to select the physically optimal one;
5. Test a phase-jump turning diagnostic;
6. Deploy the operator in a fully automated, real-time GitHub Actions pipeline fed by GFS analyses.

Throughout, we follow an honest-research protocol: all conclusions are labeled with evidence level, all data provenance is documented, and all numbers are recomputed from the published JSON artifacts before release.

---

## 2. Data and Methods

### 2.1 Data

| Dataset | Variables | Resolution | Source |
|---|---|---|---|
| ERA5 reanalysis (Hersbach et al., 2020) | u, v at 850 / 500 / 200 hPa | 0.25° | Copernicus Climate Data Store (CDS) |
| TC best-track (6 typhoons) | center position, wind | 6-hourly / 12-hourly | JTWC / HKO composite (`tracks_6typhoon.json`) |
| GFS operational analysis | u, v at 500 hPa | 0.25° | NOMADS (NCEP) |
| Active TC list (live mode) | names, positions | ~6-hourly advisories | cyclocane (JTWC advisories) |

**Sample construction.** Six typhoons with full life cycles were selected: **HATO (2017)**, **MANGKHUT (2018)**, **MERANTI (2016)**, **HAGIBIS (2019)**, **SAOLA (2023)**, and **GONI (2020)**. For each 12-hourly analysis time *t*, we computed the WN1 phase from ERA5 and compared it with the observed motion direction over the subsequent 12 h (forward validation). From 360 raw time steps, quality filters retained **76 typhoon-grade samples**: wind ≥ 40 kt (TC intensity, excluding background/weak fields) and WN1 amplitude ≥ 1.0 m/s (excluding phase noise). A pseudo-real-time protocol was enforced: only information available at time *t* was used to predict motion *t → t+12h*; no future data leaked into the predictor.

### 2.2 The WN1 Shadow Operator

At a given pressure level, on the ring of radius 1–5° around the TC center:

1. Compute the radial wind $v_r$ at each azimuth;
2. Take the azimuthal Fourier decomposition; extract the wavenumber-1 component;
3. Define the **phase φ** as the azimuth angle of the WN1 radial-wind dipole maximum (inward/outward asymmetry axis);
4. Define the **ellipticity** ellipt = |A₂/A₁| (wavenumber-2 amplitude / wavenumber-1 amplitude), a structural-purity metric: low ellipt indicates a clean, single-dipole asymmetry; high ellipt indicates a distorted or multi-lobed structure where the WN1 phase is unreliable.

The physical hypothesis: the WN1 radial-wind asymmetry at the steering level is the "shadow" cast by the deep steering flow onto the vortex structure; its phase therefore points in the direction the storm is being (and will be) advected.

### 2.3 Forward Validation

For each sample at time *t* with WN1 phase φ(t):

- **Observed motion direction** ψ(t → t+12h): azimuth of the displacement vector over the next 12 h;
- **Angular error** ε = circular distance(φ(t), ψ(t→t+12h));
- **Random baseline**: 2,000 Monte Carlo draws of uniform phase in [0°, 360°) → expected mean error ≈ 90°.

### 2.4 Uncertainty (UQ) Filter

Following the structural-purity rationale, we apply a threshold **ellipt ≤ 0.4** (a WN1-dominated, clean dipole). Samples above threshold are flagged untrustworthy and excluded from the "trusted" statistics. Additional diagnostics (phase jump Δφ, shape transition) were tested but did not add residual discrimination beyond ellipt ≤ 0.4 (see §3.5).

### 2.5 Reproducibility

All analysis scripts (`wn1_forward_validation.py`, `wn1_500_vs_200.py`, `wn1_uq_shape.py`, `wn1_lead_time_scan.py`, `wn1_track.py`) and result JSONs are published alongside this paper (GitHub + Zenodo). Raw ERA5/GFS fields are not packaged (size), but download scripts and provenance allow full re-pull.

---

## 3. Results

### 3.1 Forward skill (12 h lead, 76 samples)

| Metric | Value |
|---|---|
| Samples (6 typhoons) | 76 |
| Mean angular error (all) | **23.9°** (500 hPa) |
| Median angular error | ~18° |
| Fraction < 45° (all) | 82% |
| Random baseline | 90° |
| Skill vs. random | ≈ 73% |

Per-storm mean errors: HATO 9.9° (n=5), MANGKHUT 16.0° (n=18), MERANTI 15.4° (n=10), HAGIBIS 14.9° (n=15), GONI 29.1° (n=10), SAOLA 44.9° (n=18). Four of six storms average below 20°; SAOLA and GONI are the two weak cases (see §3.5).

### 3.2 The three critiques — answered with data

| Critique | Claim | Data response |
|---|---|---|
| **Simultaneity** (φ describes the present, not the future) | "You measured the phase at peak and the motion at peak — that's diagnosis, not forecast." | Forward validation uses φ(t) → motion(t→t+12h). Mean error 23.9° ≪ random 90°; moreover forward error is not larger than simultaneous error, and the lead-time curve (§3.4) degrades only slowly to 48 h — a pure diagnostic would collapse with lead. |
| **Downwind sampling** (all 6 storms sampled near peak) | "You only measured storms while they were moving steadily; turning cases would fail." | Turning segments (Δ > 20°/12h): mean 33.4°, median 19.7° (200 hPa); the sharpest turns do not collapse. The two genuinely weak storms (SAOLA, GONI) *include* the strongest turning/looping behavior — precisely the regime that defeats naive sampling. |
| **Vertical wind shear (VWS) confounding** | "The WN1 phase is just a VWS shadow." | High-VWS samples (≥10 m/s) are *more* accurate than low-VWS (24.9° vs 33.4°, 200 hPa); the residual of phase in the VWS coordinate frame is R = 0.963, inconsistent with a shear-alignment explanation. |

### 3.3 Pressure-level comparison (selecting the optimal level)

| Level | All mean error | ellipt ≤ 0.4 mean | Fraction < 45° | Notes |
|---|---|---|---|---|
| 850 hPa | 52.0° | 33.4° | 73% | Too close to inner core; convective heating asymmetry contaminates the phase |
| 200 hPa | 31.8° | 21.4° | 91% | Good; outflow-level signal |
| **500 hPa** | **23.9°** | **21.1°** | **90–91%** | **Optimal: classic steering level** |

The 850 hPa hypothesis (low-level steering should be more physical) is **rejected**: its phase differs from 200 hPa by a mean of 70.7°, only 13% of samples agree within 30°, and it performs worse at every aggregate metric. The two levels carry independent physical information rather than the same steering flow at different heights. 500 hPa wins or ties on all six storms and becomes the final operational level (200 hPa retained as fallback).

### 3.4 Lead-time degradation (forecast window)

Using the trusted subset (ellipt ≤ 0.4, n = 43):

| Lead | 6 h | 12 h | 24 h | 48 h |
|---|---|---|---|---|
| Mean error | 20.1° | 21.4° | 23.5° | 26.6° |
| Slope | — | — | — | ≈ 0.16°/h |

Degradation from 6 h to 48 h is only +6.5° — the WN1 phase tracks a *persistent* motion direction (deep steering is quasi-stationary over 2 days). The signal is therefore not a "12 h short-lived" quantity; it supports **24–48 h trend forecasting**. Turning segments show faster degradation (+7° over 48 h), as expected from the intrinsic unpredictability of the motion change itself, but remain below 45° at 24 h.

### 3.5 Failure modes and the UQ envelope

Two storms dominate the error budget, with distinct, detectable signatures:

- **GONI (2020)** — a rapid-intensification case. During RI, internal dynamics dominate over environmental advection; the WN1 phase is temporarily "captured" by the structural dipole of the intensifying core (10 samples, 29.1°). The failure is identifiable: high ellipticity / phase instability during RI onset.
- **SAOLA (2023)** — a looping storm. The WN1 phase lags behind rapid direction changes (up to 108°/12h turns); it cannot keep up with the evolving steering during loops (18 samples, 44.9°). Signature: phase jump Δφ and strong ellipticity excursions.

Both failure modes carry *diagnostic signatures* that the UQ filter (ellipt ≤ 0.4) largely captures: after filtering, trusted samples reach 21.1° / 90% < 45°. The filter therefore defines an **operating envelope**: when the storm structure is a clean WN1 dipole, the phase is a reliable directional proxy; when it is not, the operator flags itself as untrustworthy rather than silently failing.

### 3.6 Turning alarm (phase-jump diagnostic)

A complementary diagnostic uses the **phase jump Δφ** between consecutive 12-h analyses as an early-warning signal for turning. Over n = 69 consecutive pairs:

| Diagnostic | Value |
|---|---|
| Pearson ρ(Δφ, Δmotion) | **0.589** |
| Mean |Δmotion| when Δφ ≥ 15° | 17.7° / 12 h |
| Mean |Δmotion| when Δφ < 15° | 9.2° / 12 h |

A large phase jump is associated with a roughly 2× larger subsequent direction change — a physically interpretable, cheap early-warning trigger for recurvature, complementary to official "medium-confidence" turning forecasts.

### 3.7 Real-time deployment (GFS + GitHub Actions)

The operator is deployed as an open, fully automated pipeline (GitHub Actions, 6-hourly):

1. **cyclocane** → active TC list (JTWC advisories);
2. **NOMADS GFS** 0.25° 500 hPa analysis (nearest available cycle, with cycle-fallback logic that prefers the most recent published analysis rather than a stale previous-day cycle);
3. Compute WN1 phase + ellipticity per storm;
4. Append to `wn1_history.json` (idempotent, storm × cycle de-duplicated);
5. Auto-commit and push.

Live results (2026-08-09) include CHAN-HOM (WN1 219.7° SW, ellipt 0.086 — trusted) and PEILOU (67.8° ENE, ellipt 0.72 — flagged), demonstrating both the trusted and flagged regimes in real time. Latency matches operational reality: GFS analyses are available 3–6 h after nominal cycle time, comparable to official forecast cycles.

---

## 4. Discussion

### 4.1 Honest engineering assessment

The WN1 shadow is **not** proposed as a standalone track forecast system. Translating the directional error to distance: for a fast storm moving 480 km in 24 h, 23.9° ≈ 190 km cross-track, which is 30–60% worse than official 24 h errors (~100–150 km). For slow-moving storms the same angular error corresponds to a much smaller distance — the comparison is speed-dependent.

What the operator *does* provide is an independent, physically orthogonal signal:

| Aspect | Official dynamical models | WN1 shadow |
|---|---|---|
| Main error source | Environmental field prediction (subtropical high evolution) | Structural change (RI/ERC, looping) |
| Weakest regime | Recurvature ("medium confidence") | Structure-unstable periods (caught by UQ) |
| Strongest regime | Steady-track persistence | Turning alarm, structural forewarning |
| Data cost | Satellites + aircraft + supercomputers | Free public GFS + cyclocane |
| Automation | Institutional staff | Fully automated (GitHub Actions) |

Because the error sources are largely orthogonal, an optimal (inverse-variance) blend is estimated to improve official 24 h error from ~125 km to ~104 km (**~15–20% ideal; ~5–10% realistic** after accounting for shared GFS wind-field error sources). This is an *estimate, not a measured value*; blending weights must be regressed from live accumulated data (20–30 storm cases ≈ one season).

### 4.2 The genuine value: turning alarms

The highest-value operational use is the **shadow-deviation alarm**: when ellipt ≤ 0.4 (trusted structure) and the WN1 phase differs from the official forecast direction by > 30°, flag the disagreement as an independent second opinion. The value of such a signal is not a 5% average-error improvement; it is avoiding a single catastrophic misjudgment during recurvature — exactly the regime where official models express "medium confidence."

### 4.3 Limitations

- Sample size: 76 samples / 6 storms; sharp-turn (Δ > 40°) subset n = 3 — statistically weak;
- SAOLA/GONI weak cases suggest regime dependence that is signature-detectable but not yet mechanistically closed;
- Track data are sparse in places (e.g., HATO contributes only 5 qualified samples);
- 500 hPa level: GFS 0.25° cannot resolve the inner core (eye wall ~50 km), so the operator relies on the outer-ring structure — this is a feature (robustness) and a limitation (no core dynamics);
- The forward validation is pseudo-real-time on reanalyses; live GFS validation accumulates via the automated pipeline;
- Blend-with-official numbers are estimates pending live data.

---

## 5. Conclusion

The WN1 Eccentric Shadow — the phase of the wavenumber-1 radial-wind asymmetry on the inner ring at 500 hPa — is a structurally grounded, forward-validated, uncertainty-controlled directional precursor for tropical cyclone motion:

- **12 h forward mean error 23.9°** (76 samples, 6 typhoons) vs. 90° random baseline;
- **Trusted envelope (ellipt ≤ 0.4): 21.1°, 90% < 45°**;
- **Usable window ≥ 48 h** (26.6° at 48 h, slope ≈ 0.16°/h);
- **500 hPa optimal** (850 rejected, 200 fallback);
- **Phase-jump turning alarm** (ρ = 0.589);
- **Fully automated real-time deployment** on free public data (GitHub Actions + GFS + cyclocane).

We position it as a **complementary second-opinion signal** for operational forecasting — an independent, free, automated structural shadow with a well-defined trust envelope — with its highest value in **recurvature alarms**. The operator does not replace dynamical models; it adds an orthogonal measurement channel that official forecasts do not use.

---

## Data Availability

- **Analysis scripts** (all): `wn1_forward_validation.py`, `wn1_500_vs_200.py`, `wn1_850_vs_200.py`, `wn1_uq_shape.py`, `wn1_lead_time_scan.py`, `wn1_dphi_turning.py`, `wn1_track.py` (GitHub repo `MMDR10/typhoon-wn1-shadow`);
- **Result JSONs**: `wn1_forward_validation.json`, `wn1_500_vs_200.json`, `wn1_850_vs_200.json`, `wn1_uq_shape_v2.json`, `wn1_lead_time.json`, `wn1_dphi_turning.json`, `wn1_history.json`;
- **Figures**: `wn1_3levels_CDF.png`, `wn1_500_vs_200.png`, `wn1_850_vs_200.png`, `lead_time_curve.png`, `wn1_dphi_analysis.png`;
- **Raw ERA5 / GFS fields**: not packaged (size); reproducible via download scripts; provenance documented in §2.1.

## References

1. Chan, J. C. L., & Gray, W. M. (1982). Tropical cyclone movement and surrounding flow relationships. *Monthly Weather Review*, 110(10), 1354–1374.
2. Hersbach, H., et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999–2049.
3. NCEP. Global Forecast System (GFS) 0.25° operational analyses. NOMADS, NOAA/NCEP.
4. Joint Typhoon Warning Center (JTWC) best-track / advisories; HKO tropical cyclone database.
5. cyclocane — real-time TC advisory aggregator (https://www.cyclocane.com/).

---

*All conclusions labeled per evidence-level protocol: ✅ cross-validated (≥2 independent sources/analyses), ⚠️ single-source/limited-sample, ❌ unverified. See individual report artifacts in the repository for full evidence tables.*
