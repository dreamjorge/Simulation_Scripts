# Paraxial Beam Simulation Project

This project demonstrates a **Component-Based Architecture** for simulating various types of paraxial beams, including **Gaussian** and **Hermite-Gaussian** beams. The architecture is designed for **reusability**, **extensibility**, and **testability**, allowing new beam types, services, and utilities to be easily added without modifying the core system.

## Project Structure

The project is organized into the following directories:
It seems there was an issue generating the download link again. However, you can copy the content below into a file named README.md in your project folder.
Paraxial Beam Simulation Project

This project demonstrates a Component-Based Architecture for simulating various types of paraxial beams, including Gaussian and Hermite-Gaussian beams. The architecture is designed for reusability, extensibility, and testability, allowing new beam types, services, and utilities to be easily added without modifying the core system.
Project Structure

The project is organized into the following directories:

````
ParaxialBeamSimulation/
│
├── entities/
│   ├── IBeam.m                   % Common interface for all beam types
│   ├── ParaxialBeam.m             % Gaussian Beam (base beam type)
│   └── HermiteGaussianBeam.m      % Hermite-Gaussian Beam (extends ParaxialBeam)
│
├── services/
│   └── BeamService.m              % Service for beam-related calculations
│
├── utils/
│   ├── RayTracing.m               % Utility for ray tracing and beam propagation
│   └── VisualizationUtils.m       % Plotting and visualization functions
│
├── app/
│   ├── main.m                     % Main simulation script for beam properties
│   └── main_ray_tracing.m         % Script for running ray tracing visualization
└── tests/
    └── test_beam_service.m        % Unit tests for the BeamService component
````

## Components

### 1. Beam Models (`entities/`)
- **`IBeam.m`**: Defines the common interface that all beam models must implement. This ensures that beam models like **Gaussian** and **Hermite-Gaussian** can be used interchangeably by the services and utilities.
- **`ParaxialBeam.m`**: A basic Gaussian beam model that calculates beam waist, radius of curvature, phase shift, and intensity profile.
- **`HermiteGaussianBeam.m`**: Extends `ParaxialBeam` to add support for higher-order Hermite-Gaussian modes.

### 2. BeamService (`services/`)
- **`BeamService.m`**: A service component that provides operations for calculating beam properties like waist size, phase shift, and intensity profile. This service works with any beam type implementing `IBeam`.

### 3. Ray Tracing and Visualization Utilities (`utils/`)
- **`RayTracing.m`**: A utility that traces the wavefront propagation of a beam and visualizes its phase front as it propagates through space.
- **`VisualizationUtils.m`**: Provides functions for visualizing beam intensity profiles and other beam properties.

### 4. Main Application Scripts (`app/`)
- **`main.m`**: Runs the beam simulation, calculating and plotting the properties of the selected beam.
- **`main_ray_tracing.m`**: Visualizes the wavefront propagation of the selected beam using ray tracing.

### 5. Unit Testing (`tests/`)
- **`test_beam_service.m`**: Unit tests to validate the `BeamService` component.

## Installation

To use this project:

1. **Clone the repository** to your local machine.
2. **Open the MATLAB environment**.
3. Ensure all paths to the project are added to MATLAB’s working directory.

## Usage

### Simulate Beam Properties (`main.m`)
1. Open `app/main.m` in MATLAB.
2. Set the beam type you want to simulate (Gaussian or Hermite-Gaussian).
3. Run the script to visualize the beam properties such as intensity and beam waist at a specific propagation distance.

### Perform Ray Tracing (`main_ray_tracing.m`)
1. Open `app/main_ray_tracing.m` in MATLAB.
2. Run the script to visualize the wavefront propagation of the beam over a range of distances.

## Example

### Example Code to Simulate a Gaussian Beam

```matlab
% Define beam parameters
wavelength = 633e-9;  % Wavelength in meters
waist = 1e-3;         % Waist size in meters
position = 0;         % Starting position (z = 0)

% Create a Gaussian beam
beam = ParaxialBeam(wavelength, waist, position);

% Pass the beam to the BeamService
beamService = BeamService(beam);

% Calculate properties at z = 1 meter
z = 1;
waist_at_z = beamService.calculateBeamWaist(z);
fprintf('Beam waist at z = %.2f m: %.2e m\\n', z, waist_at_z);

% Visualize intensity at z = 1 meter
grid_size = [-5, 5];
plotBeamIntensity(beam, z, grid_size);
```

### Extending the Project

To add a new beam type (e.g., Laguerre-Gaussian):

    Create a new class that implements the IBeam interface.
    Define the necessary beam properties such as beam waist, radius of curvature, phase shift, and intensity profile.
    Pass the new beam type into BeamService and reuse the ray tracing and visualization utilities.

## Testing

To run the unit tests:

    Open tests/test_beam_service.m.
    Run the test to verify that the beam service works as expected.

```matlab
test_beam_service;
```
### Future Improvements

    Add support for additional beam types like Laguerre-Gaussian or Bessel Beams.
    Introduce more advanced visualizations, including 3D rendering of wavefronts.
    Implement more complex beam interactions such as interference patterns and transformations.

License

This project is licensed under the MIT License.
