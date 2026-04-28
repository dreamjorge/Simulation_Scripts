# github-release-workflow Specification

## Purpose

Define a GitHub Actions workflow that builds Octave `.tar.gz` and MATLAB `.mltbx` packages on tagged releases and attaches them as release artifacts.

## Requirements

### Requirement: Workflow Trigger

The workflow SHALL trigger on Git tag pushes matching `refs/tags/v*`.

#### Scenario: Tag push triggers workflow

- GIVEN a commit is tagged `v2.0.0` and pushed
- WHEN `git push origin v2.0.0` is executed
- THEN GitHub Actions SHALL trigger the release workflow
- AND a Release SHALL be created on GitHub

### Requirement: Build Octave Package

The workflow SHALL build the `.tar.gz` Octave package from the tagged source.

#### Scenario: Octave package is built and attached

- GIVEN the workflow runs on tag `v2.0.0`
- THEN it SHALL produce `simulation_scripts-2.0.0.tar.gz`
- AND it SHALL be attached to the GitHub Release as an asset

### Requirement: Build MATLAB Toolbox

The workflow SHALL build the `.mltbx` MATLAB toolbox from the tagged source.

#### Scenario: MATLAB toolbox is built and attached

- GIVEN the workflow runs on tag `v2.0.0`
- THEN it SHALL produce `simulation_scripts-2.0.0.mltbx`
- AND it SHALL be attached to the GitHub Release as an asset

#### Scenario: MATLAB CI uses free open-source license

- GIVEN `matlab-actions/setup-matlab@v3` with `release: latest`
- WHEN the workflow runs on an open-source public repo
- THEN it SHALL use GitHub's hosted MATLAB license at no cost

### Requirement: Validate Before Release

The workflow SHALL run `portable_runner.m` before building packages.

#### Scenario: Tests pass before packaging

- GIVEN the tag is pushed
- WHEN the workflow begins
- THEN it SHALL first run `octave --no-gui --eval "setpaths; status = portable_runner(); if status ~= 0, error(...) end"`
- AND if tests fail, the workflow SHALL fail
- AND no release artifacts SHALL be built

### Requirement: Release Notes

The workflow SHALL generate release notes from commit messages since last tag.

#### Scenario: Release notes describe changes

- GIVEN commits since `v1.9.0`
- WHEN a release is created for `v2.0.0`
- THEN the release body SHALL list changes since `v1.9.0`
- AND it SHALL include a link to the full changelog

### Requirement: Artifact Naming

Release artifacts SHALL follow the naming convention `simulation_scripts-<VERSION>.{tar.gz,mltbx}`.

#### Scenario: Artifacts named correctly

- GIVEN version `2.0.0`
- THEN the Octave artifact SHALL be `simulation_scripts-2.0.0.tar.gz`
- AND the MATLAB artifact SHALL be `simulation_scripts-2.0.0.mltbx`