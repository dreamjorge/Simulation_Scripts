# wavefront-extraction Specification

## Purpose

Extract wavefront (phase) information from complex optical field data.

## ADDED Requirements

### Requirement: Constructor with Field Data

The system MUST accept a complex field `[Ny x Nx]` array along with wavelength and optional grid metadata.

#### Scenario: Minimal Constructor

- GIVEN a complex field `E` of size [256, 256]
- WHEN user creates `wf = Wavefront(E, lambda)`
- THEN the system SHALL store `E` and `lambda`
- AND SHALL derive grid size from `size(E)`

#### Scenario: Constructor with Grid

- GIVEN complex field `E`, wavelength `lambda`, and GridUtils `grid`
- WHEN user creates `wf = Wavefront(E, lambda, grid)`
- THEN the system SHALL store grid metadata (dx, dy, Dx, Dy)

### Requirement: Getter Methods

The system MUST provide these getters:
- `getIntensity()` — returns `abs(E).^2`
- `getPhase()` — returns `angle(E)` wrapped to [-pi, pi]
- `getField()` — returns original complex field

#### Scenario: Get Intensity

- GIVEN `Wavefront` instance `wf`
- WHEN user calls `I = wf.getIntensity()`
- THEN `I` SHALL be [Ny x Nx] real array

#### Scenario: Get Phase

- GIVEN `Wavefront` instance `wf`
- WHEN user calls `phi = wf.getPhase()`
- THEN `phi` SHALL be [Ny x Nx] real array in [-pi, pi]

### Requirement: Grid Info Storage

The system SHOULD store these grid metadata properties:
- `dx`, `dy` — grid spacing in meters
- `Dx`, `Dy` — grid extent in meters
- `wavelength` — in meters

#### Scenario: Grid Metadata Preserved

- GIVEN `Wavefront(E, lambda, grid)` constructor
- WHEN `wf.dx` is accessed
- THEN it SHALL return `grid.dx` value