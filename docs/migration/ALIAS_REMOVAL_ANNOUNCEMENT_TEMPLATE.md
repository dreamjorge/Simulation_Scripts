# Announcement Template — Legacy Alias Removal

## Summary

In release `v2026.05-legacy-alias-removal`, deprecated aliases
`HankeleHermite` and `HankeleLaguerre` will be removed.

## Why

These aliases were kept temporarily for legacy compatibility. Modern API usage
has been stable, and migration guardrails are now in place.

## What to change

Replace:

```matlab
HankeleHermite(...)  -> HankelHermite(...)
HankeleLaguerre(...) -> HankelLaguerre(...)
```

Static delegation calls:

```matlab
HankeleHermite.getPropagateCartesianRays(...)  -> HankelHermite.getPropagateCartesianRays(...)
HankeleLaguerre.getPropagateCylindricalRays(...) -> HankelLaguerre.getPropagateCylindricalRays(...)
```

## Timeline

- Announcement date: ____
- Effective removal release: `v2026.05-legacy-alias-removal`

## Support

If you still depend on these aliases, open an issue with:

1. script path(s)
2. minimal reproducible snippet
3. expected migration timeline
