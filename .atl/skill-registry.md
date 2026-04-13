# Skill Registry - Simulation_Scripts

Generated: 2026-04-07

## Project Context
- **Language**: MATLAB/Octave (.m files)
- **Domain**: Optical physics simulation (Gaussian, Hermite-Gauss, Laguerre-Gauss beams)
- **Architecture**: Object-oriented MATLAB with classes

## User Skills (Global)
Loaded from `~/.config/opencode/skills/`. Triggered by task context.

| Skill | Trigger Context |
|-------|-----------------|
| sdd-* | SDD workflow phases |
| judgment-day | "judgment day", "dual review" |
| branch-pr | Creating pull requests |
| issue-creation | Creating GitHub issues |
| receiving-code-review | Receiving code review feedback |
| requesting-code-review | Requesting code review |

## Project Conventions
- **Tests**: `tests/test_all.m` — run with `octave --no-gui --eval "run('tests/test_all.m')"`
- **CI**: GitHub Actions (`.github/workflows/octave.yml`)
- **No linting/type-checking**: MATLAB/Octave doesn't have standard tooling

## Notes
- This is a physics simulation project, not a typical software project
- MATLAB/Octave has no standard TDD framework
- Strict TDD Mode not applicable
