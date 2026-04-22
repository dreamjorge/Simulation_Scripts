# Legacy Compatibility Layer

## Status

As of release track `legacy-alias-removal-r1`, deprecated aliases were removed:

- `HankeleHermite`
- `HankeleLaguerre`

Use modern classes directly:

- `HankelHermite`
- `HankelLaguerre`

## Migration Notes

See:

- `docs/migration/LEGACY_MIGRATION_PLAN.md`
- `docs/migration/ALIAS_REMOVAL_RELEASE_PLAN.md`

If a historical script still references `Hankele*`, update constructor/static
calls to `Hankel*` equivalents before running.
