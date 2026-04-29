# Skill Registry for Simulation_Scripts

## Project Standards

### MATLAB/Octave Coding standards
- Use `classdef` for new components.
- Document methods with help headers.
- Use SI units for physical simulations.
- Ensure compatibility with Octave 11.1.0+ and MATLAB R2020b+.

### Testing Standards
- All new logic must have corresponding tests in `tests/`.
- Use `verifyEqual` or `assertEqual` with tolerances (`AbsTol`, `RelTol`) for floating point comparisons.

## Project Context
- **Language**: MATLAB/Octave (.m files)
- **Domain**: Optical physics simulation (Gaussian, Hermite-Gauss, Laguerre-Gauss beams)
- **Architecture**: Object-oriented MATLAB with classes

## User Skills (Global)
| Skill | Trigger Context |
|-------|-----------------|
| sdd-* | SDD workflow phases |
| judgment-day | "judgment day", "dual review" |
| branch-pr | Creating pull requests |
| issue-creation | Creating GitHub issues |
| receiving-code-review | Receiving code review feedback |
| requesting-code-review | Requesting code review |

## Project Conventions
- **Tests**: `tests/test_all.m` — run with `& "C:\Users\uidn7961\AppData\Local\Programs\GNU Octave\Octave-11.1.0\mingw64\bin\octave-cli.exe" --no-gui --eval "run('tests/test_all.m')"`
- **Octave Path**: `C:\Users\uidn7961\AppData\Local\Programs\GNU Octave\Octave-11.1.0\mingw64\bin\octave-cli.exe`
- **CI**: GitHub Actions (`.github/workflows/octave.yml`, `.github/workflows/matlab.yml`, `.github/workflows/release.yml`)
- **No linting/type-checking**: MATLAB/Octave doesn't have standard tooling
