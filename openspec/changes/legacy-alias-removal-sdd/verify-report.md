# Verify Report: legacy-alias-removal-sdd (Phase 1)

Date: 2026-04-22
Branch: `feat/legacy-alias-removal-sdd`

## 1) Baseline execution evidence

### 1.1 Guardrail

Command:

```powershell
octave-cli --no-gui --quiet --eval "setpaths; run('tests/modern/test_LegacyAliasGuardrail.m'); exit"
```

Result:

- PASS: no deprecated `Hankele*` aliases in `examples/canonical` and `tests/modern`

### 1.2 Portable suite

Command:

```powershell
octave-cli --no-gui --quiet --eval "set(0,'defaultfigurevisible','off'); setpaths; run('tests/portable_runner.m'); exit"
```

Result:

- `Tests Pasados: 31`
- `Tests Fallados: 0`
- `ESTADO: ÉXITO`

### 1.3 Temporary alias-removal portable run (Gate B.2 evidence)

Method:

- Temporarily renamed:
  - `legacy/compat/HankeleHermite.m` -> `.removed`
  - `legacy/compat/HankeleLaguerre.m` -> `.removed`
- Ran `tests/portable_runner.m`
- Restored both files automatically after run

Result:

- Initial run: `Tests Pasados: 28`, `Tests Fallados: 3`
  - failures isolated to alias-only legacy tests
- After transitioning legacy alias tests to migration assertions:
  - rerun result: `Tests Pasados: 31`, `Tests Fallados: 0`
  - full portable suite green with aliases temporarily removed

## 2) Repository reference scan (`Hankele*`)

Command:

```text
grep pattern: \bHankeleHermite\b|\bHankeleLaguerre\b
include: *.{m,md,yml,yaml,txt}
```

### 2.1 Expected references

- `legacy/compat/HankeleHermite.m`, `legacy/compat/HankeleLaguerre.m` (current adapters)
- `tests/legacy_compat/*` (legacy compatibility coverage)
- `docs/migration/*` and `legacy/compat/README.md` (migration/deprecation docs)
- legacy/historical examples outside canonical (`examples/Main*.m`, thesis/test scripts)
- `tests/modern/test_LegacyAliasGuardrail.m` (forbidden token list for enforcement)

### 2.2 Unexpected references

- None found in `src/`
- None found in `examples/canonical/`
- None found in `tests/modern/` (except guardrail token constants by design)

## 3) Readiness gates status

### Gate A — Usage

- ✅ PASS (for protected modern surfaces)
  - no `Hankele*` in `src/`, `examples/canonical/`, `tests/modern` runtime usage
- ⚠️ PENDING (external dependency signal)
  - tracking checklist created: `docs/migration/USAGE_SIGNAL_CHECKLIST.md`
  - announcement template prepared: `docs/migration/ALIAS_REMOVAL_ANNOUNCEMENT_TEMPLATE.md`

### Gate B — Test

- ✅ PASS
  - Baseline stability with aliases present: `portable_runner` green (31/0)
  - Temporary alias-removal run executed with full portable suite green (31/0)
  - Legacy alias tests transitioned to migration assertions so they pass both
    pre-removal and post-removal modes
  - Post-removal mode is explicitly gated by environment flag:
    `LEGACY_ALIAS_REMOVAL_MODE=1`
    - default mode (flag unset): missing aliases => test failure
    - removal mode (flag set): missing aliases => expected pass

### Gate C — Docs

- ✅ PASS
  - README and migration docs do not recommend alias usage for new code
  - breaking-change release note template documented in
    `docs/migration/ALIAS_REMOVAL_RELEASE_PLAN.md`

### Gate D — Release

- ✅ PASS
  - named milestone and target tag defined in
    `docs/migration/ALIAS_REMOVAL_RELEASE_PLAN.md`
    (`legacy-alias-removal-r1` / `v2026.05-legacy-alias-removal`)
  - warning-cycle evidence linked across consecutive checkpoints:
    2026-04-15 and 2026-04-22

## 4) Decision

Current status: **PARTIALLY READY FOR ALIAS REMOVAL**.

Blocking items:

1. Complete Usage gate item for external dependency signal (user-impact check).
