# wavefront-visualization Specification

## Purpose

Visualize wavefront phase maps and Zernike decomposition results.

## ADDED Requirements

### Requirement: Wavefront Phase Map Plot

The system MUST provide a method to plot the 2D wavefront phase map.

#### Scenario: Plot Wavefront Phase

- GIVEN `Wavefront` instance `wf`
- WHEN user calls `wf.plotWavefront()`
- THEN it SHALL create a figure with a 2D phase map
- AND colorbar SHALL show phase in radians or waves

### Requirement: Intensity Plot

The system MUST provide a method to plot the intensity distribution.

#### Scenario: Plot Intensity

- GIVEN `Wavefront` instance `wf`
- WHEN user calls `wf.plotIntensity()`
- THEN it SHALL create a figure with 2D intensity map
- AND axis labels SHALL be in meters

### Requirement: Zernike Coefficient Bar Chart

The system MUST provide a method to plot Zernike coefficients as a bar chart.

#### Scenario: Plot Zernike Coefficients

- GIVEN `Wavefront` instance with fitted coefficients
- WHEN user calls `wf.plotZernikeCoeffs(coeffs)`
- THEN it SHALL create a bar chart with coefficient values
- AND x-axis labels SHALL show Zernike names (Piston, Tilt X, ...)

### Requirement: Phase Slice Plot

The system MAY provide a 1D phase slice plot along a cross-section.

#### Scenario: Plot Phase Slice

- GIVEN `Wavefront` instance with phase `phi`
- WHEN user calls `wf.plotPhaseSlice('x', Ny/2)`
- THEN it SHALL plot `phi(Ny/2, :)` vs x-axis