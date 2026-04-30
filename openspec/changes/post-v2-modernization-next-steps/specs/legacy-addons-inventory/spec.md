# Legacy Addons Inventory Specification

## Purpose

Ensure `ParaxialBeams/Addons/` is classified before any migration, deletion, or vendored-code decision.

## Requirements

### Requirement: Addon Classification

Every top-level addon file and addon subdirectory MUST be classified before cleanup decisions are made.

#### Scenario: Inventory is created

- GIVEN files exist under `ParaxialBeams/Addons/`
- WHEN the addons inventory is produced
- THEN each top-level `.m` file and subdirectory MUST appear in `docs/ADDONS_INVENTORY.md`
- AND each entry MUST have exactly one classification.

### Requirement: Supported Classifications

Addon entries MUST use one of: runtime-required, plotting-only, vendored-third-party, removable-candidate, or needs-investigation.

#### Scenario: Unknown usage

- GIVEN addon usage cannot be proven from tests, examples, or source references
- WHEN the inventory is written
- THEN the entry MUST be marked `needs-investigation`
- AND it MUST NOT be deleted in this change.

### Requirement: No Deletion Without Follow-up SDD

This change MUST NOT delete addon files solely because they appear unused.

#### Scenario: Removable candidate found

- GIVEN an addon is classified as `removable-candidate`
- WHEN implementation completes
- THEN the file MUST remain present
- AND a follow-up SDD change SHOULD be referenced for removal.
