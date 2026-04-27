# Tasks: Octave CI Failure Propagation

## Phase 1: CI Command Fix

- [x] 1.1 Modify `.github/workflows/octave.yml` portable suite command to store `status = portable_runner()`.
- [x] 1.2 Add Octave-side `if status ~= 0, error(...) end` guard while keeping `tee portable-tests.log`.

## Phase 2: Static Verification

- [x] 2.1 Verify `.github/workflows/matlab.yml` remains unchanged and already checks status.
- [x] 2.2 Verify `tests/portable_runner.m` remains a function returning `totalFailed` without forced exit.

## Phase 3: Runtime Verification

- [x] 3.1 Check whether `octave` is executable in this session.
- [ ] 3.2 If available, run `octave --no-gui --eval "addpath('tests'); status = portable_runner(); if status ~= 0, error('portable_runner failed with %d failing tests', status); end"`.
- [x] 3.3 If unavailable, document local verification blocker and command for user `uib95096`.

> Runtime note: current session user is `automotive-wan\uidn7961`; `where.exe octave` did not find Octave in PATH. MATLAB is in PATH, but this session previously failed license checkout. User reports Octave/MATLAB are usable under `uib95096`; run task 3.2 from that user context.
