# Release branch strategy: orphaned `main`, history-visible `develop`

**Status: agreed, not yet executed.** Decided 2026-09-05, on
`release-R1.0-prep`. Execution is deferred to R1.0 ship time --
today's state is still prep, not a release.

## Problem

Users of a public release should see only the code base, not the
internal development trail (feature-branch churn, WIP commits,
internal-only dataset/personal-name references, etc.). Developers need
the opposite: full commit history so normal git operations
(blame, bisect, merge, rebase) keep working.

The naive way to hide history -- squash-merge every feature branch
into the public branch -- breaks a specific case: if a developer
has an old feature branch still open (unmerged) and pulls the public
branch into it to pick up other changes, git sees the squashed commit
as unrelated to the original commits it was made from (different
SHAs, same diffs). Any later merge attempt either duplicates the diff
or conflicts against it. This gets worse the longer that branch stays
open. As of this date there are 11 other feature branches sitting on
`origin` (`cube-stat`, `diagnostics-heuristics`, `doc-update`,
`encapsulate-rmsynth`, `multi-band-tomography`, `optimise-io`,
`optimise-openmp`, `overlap-host-device-compute`, `refactor-hpc`,
`rmclean-integration`, `rmtool-logging`) -- concrete instances of
exactly this risk, not a hypothetical.

## Model

- **Feature branches** merge into `develop` normally (ordinary merge
  commits, never squashed). `develop` is the full, permanent,
  developer-facing history -- every commit that ever mattered lives
  here, forever.
- **`develop` -> `main`** is the only place history gets hidden, and
  it is one-way: `main` is rebuilt from scratch as a single **orphan**
  commit (`git checkout --orphan`, no parent) containing a snapshot of
  `develop`'s tree at release time, then force-pushed over whatever
  `main` currently is. Nothing is lost by this: every commit `main`
  used to contain remains reachable via `develop`/whatever
  `release-*-prep` branch fed it (verified for the R1.0 case
  below), just no longer reachable by walking `main`'s own ancestry.
- Developers never need to pull *from* `main` to continue unmerged
  work -- they pull from `develop`, which has full history -- so the
  squash-pull failure mode above never triggers.

## Verified state at decision time (2026-09-05)

- `main` (local and `origin/main`) at `b158dca`, confirmed a **strict
  ancestor** of `release-R1.0-prep` (`git merge-base --is-ancestor
  main release-R1.0-prep`, exit 0, confirmed). No commit unique to `main`.
- `develop` (local and `origin/develop`) at `2e9e745`, 0 commits ahead
  of `release-R1.0-prep`'s own history (clean subset -- confirms
  `develop` has only ever received fast-forwards from this branch, no
  independent divergence).
- So the first cut of the orphaned `main` will not lose anything: the
  full pre-R1.0 development history stays completely intact on
  `develop`.

## Runbook, for whenever R1.0 (or any future release) ships

1. Make sure `develop` is fully caught up first (the normal
   commit -> push -> merge-to-develop -> push cycle already used all
   session -- nothing special here, just don't skip it before step 2).
2. Build the orphan snapshot from `develop`'s then-current tip:
   ```
   git checkout --orphan main-new develop
   git add -A
   git commit -m "R1.0 release"
   ```
3. Force-push it over `main`:
   ```
   git branch -D main            # local
   git branch -m main-new main
   git push origin main --force
   ```
4. **Tag it immediately**, before any future release repeats step 2-3
   and force-pushes over this commit again -- the tag is what keeps
   this exact snapshot reachable/inspectable afterward:
   ```
   git tag v1.0
   git push origin v1.0
   ```
5. For the *next* release (v1.1, v2.0, ...): repeat steps 2-4 exactly,
   always re-orphaning off `develop`'s then-current tip, always
   tagging before it can be overwritten again.

## Explicitly out of scope for this decision

The 11 other feature branches listed above are not touched by this
plan and were not reviewed for staleness/mergeability as part of this
discussion.
