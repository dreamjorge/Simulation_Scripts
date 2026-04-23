# zernike-fitting Specification

## Purpose

Fit Zernike polynomials to wavefront data for modal decomposition and optical quality assessment.

## ADDED Requirements

### Requirement: Zernike Polynomial Computation

The system MUST compute Zernike polynomials using standard Noll indexing.

#### Scenario: Zernike Standard Indices

- GIVEN Zernike index `n = 4` (defocus)
- WHEN system computes `Z4(rho, theta)`
- THEN it SHALL use standard Noll definition: Z4 = sqrt(3)*(2*rho^2 - 1)

### Requirement: Fit Zernike Coefficients

The system MUST fit Zernike coefficients to a wavefront phase map using least-squares.

#### Scenario: Fit First 36 Terms

- GIVEN `Wavefront` instance with phase `phi`
- WHEN user calls `coeffs = wf.fitZernike(36)`
- THEN `coeffs` SHALL be a [36 x 1] column vector
- AND `coeffs(1)` SHALL be piston (mean offset)

#### Scenario: Fit Custom Number of Terms

- GIVEN `Wavefront` instance with phase `phi`
- WHEN user calls `coeffs = wf.fitZernike(nTerms)`
- THEN `coeffs` SHALL have exactly `nTerms` entries

### Requirement: Reconstruct Wavefront

The system MUST reconstruct the wavefront from fitted coefficients.

#### Scenario: Round-Trip Reconstruction

- GIVEN `Wavefront` instance
- WHEN user fits and reconstructs: `phi_fit = wf.fitZernike(36); phi_rec = wf.reconstructZernike(phi_fit)`
- THEN the residual error `norm(phi - phi_rec)` SHALL be small (RMS < 1e-10)

### Requirement: Standard Noll Set (36 Terms)

The system SHALL support at minimum these 36 Zernike terms:
1. Piston, 2. Tilt X, 3. Tilt Y, 4. Defocus, 5-6. Astigmatism, 7-8. Coma, 9-10. Spherical, and higher orders up to 36.

#### Scenario: Terms 1-36 Available

- GIVEN `fitZernike(36)` call
- WHEN terms are computed
- THEN all 36 standard Noll Zernike polynomials SHALL be correctly defined

### Requirement: Zernike Names

The system SHOULD provide names for each Zernike index for display purposes.

#### Scenario: Get Zernike Name

- GIVEN `Wavefront` instance
- WHEN user calls `name = wf.zernikeName(4)`
- THEN `name` SHALL be `'Defocus'`