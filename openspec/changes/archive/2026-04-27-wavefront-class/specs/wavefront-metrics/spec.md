# wavefront-metrics Specification

## Purpose

Compute wavefront error metrics for optical quality assessment.

## ADDED Requirements

### Requirement: RMS Calculation

The system MUST compute the root-mean-square wavefront error.

#### Scenario: Compute RMS

- GIVEN `Wavefront` instance with phase `phi`
- WHEN user calls `rms = wf.computeRMS()`
- THEN `rms` SHALL be a scalar in the same units as `phi` (radians or meters)

#### Scenario: RMS from Fitted Coefficients

- GIVEN `Wavefront` instance
- WHEN user calls `rms = wf.computeRMS(coeffs)`
- THEN `rms` SHALL be computed from the fitted wavefront only (excluding residual)

### Requirement: Peak-to-Valley (PV) Calculation

The system MUST compute peak-to-valley wavefront error.

#### Scenario: Compute PV

- GIVEN `Wavefront` instance with phase `phi`
- WHEN user calls `pv = wf.computePV()`
- THEN `pv` SHALL equal `max(phi) - min(phi)`

### Requirement: Strehl Ratio Estimation

The system MUST estimate Strehl ratio from wavefront error using the Maréchal approximation.

#### Scenario: Compute Strehl

- GIVEN `Wavefront` instance with RMS phase error `sigma` in radians
- WHEN user calls `strehl = wf.computeStrehl()`
- THEN `strehl` SHALL be `exp(-sigma^2)` for small sigma
- AND `strehl` SHALL be between 0 and 1

#### Scenario: Strehl with Known WFE

- GIVEN `Wavefront` instance with phase RMS `sigma = 0.1` radians
- WHEN user calls `strehl = wf.computeStrehl()`
- THEN `strehl` SHALL be approximately `exp(-0.1^2) ≈ 0.990`

### Requirement: Wavefront Error Structure

The system MAY provide a structured output with all metrics.

#### Scenario: Get Metrics Struct

- GIVEN `Wavefront` instance
- WHEN user calls `metrics = wf.getMetrics()`
- THEN `metrics` SHALL be a struct with fields: `rms`, `pv`, `strehl`, `residualRMS`
