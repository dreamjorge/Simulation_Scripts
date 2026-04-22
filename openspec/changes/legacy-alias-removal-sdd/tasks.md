# Tasks: legacy-alias-removal-sdd

## Phase 1: Readiness evidence

- [x] 1.1 Run guardrail + portable suite and capture baseline results
- [x] 1.2 Search repository for remaining `Hankele*` references and classify expected vs unexpected
- [x] 1.3 Write `verify-report.md` with Usage/Test/Docs/Release gate status

## Phase 2: Alias removal implementation

- [ ] 2.1 Remove `legacy/compat/HankeleHermite.m`
- [ ] 2.2 Remove `legacy/compat/HankeleLaguerre.m`
- [ ] 2.3 Update any path/bootstrap docs that still suggest alias usage

## Phase 3: Test strategy transition

- [ ] 3.1 Replace alias-constructor compatibility tests with migration assertion tests (or archive alias-only tests)
- [ ] 3.2 Ensure `tests/portable_runner.m` remains green after transition
- [ ] 3.3 Confirm `tests/modern/test_LegacyAliasGuardrail.m` still passes

## Phase 4: Documentation and release closeout

- [ ] 4.1 Update `docs/migration/LEGACY_MIGRATION_PLAN.md` with completed removal status
- [ ] 4.2 Add/Update release checkpoint with removal commit hash and rollback command
- [ ] 4.3 Final review of unresolved migration risks
