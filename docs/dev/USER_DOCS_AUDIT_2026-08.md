# User-Facing Documentation Audit (pre-R1.0), 2026-08

**Status: awaiting review.** This is a checklist of findings from a
full audit of the user-facing doc set (README.md, QUICKSTART.md,
BUILD.md, docs/user/*.md, cfg/CONFIG_README.md) against the current
code on `release-R1.0-prep`, done via three independent research
passes plus direct spot-verification of the highest-severity items.

**No doc content has been changed yet.** Check off, edit, or annotate
any item below directly in this file, then hand it back for the
changes to actually be applied. Where a fix is a simple, unambiguous
correction, this file proposes exact replacement text; where it's a
real judgment call, it says so explicitly and leaves the decision
open.

Every finding below cites the doc's current claim and the code reality
it was checked against, file:line on both sides where applicable. The
highest-severity items (marked ✅ spot-checked) were independently
re-verified directly, not just taken from the research pass.

---

## A. Structural moves: ARCHITECTURE.md / PARALLELISM.md / DESIGN doc

Decided already (confirmed with you): a user doesn't need internal
design rationale to run the tools — the practical tuning guidance
these get cited for already stands on its own in
EXAMPLES.md/APP_REFERENCE.md.

- [x] **A1.** `git mv docs/user/ARCHITECTURE.md docs/dev/ARCHITECTURE.md`
- [x] **A2.** `git mv docs/user/PARALLELISM.md docs/dev/PARALLELISM.md`
- [x] **A3.** `git mv docs/user/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md docs/dev/ARCHIVED/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md`
      (archived, not just relocated — it already self-describes as
      superseded: "This section predates, and is superseded by, the
      `3.0` IO-efficiency milestone").

- [x] **A4. DONE.** ARCHITECTURE.md's "Host RAM tiling" section had a
      paragraph (added earlier this session specifically because
      EXAMPLES.md §4 needed somewhere to point for it) explaining
      `mem_frac_ram` is a fraction of *total* system RAM, not whatever
      happens to be free, and why that matters on a busy shared node.
      **Migrated**: added a new "## Shared Parameters" section to
      APP_REFERENCE.md (right after Contents, before tool 1) with this
      note stated once, and cross-referenced it from all 5 tools' own
      `mem_frac_ram` rows (`(see [Shared Parameters](#shared-parameters)
      above)`) rather than repeating the full explanation 5 times.

      **Your question — "what did you mean to check in PARALLELISM.md,
      what is the issue?" — resolved, not an issue after all.** I meant
      the same class of check as the `mem_frac_ram` one above: does
      PARALLELISM.md contain some OTHER genuinely user-actionable
      nugget (I specifically had in mind its "give compute all your
      cores; read/write don't need their own" thread-budget guidance)
      that exists ONLY there, and would go missing/unreachable once the
      file moves to `docs/dev/`? I hadn't checked yet when I wrote that
      line, so I flagged it rather than guess either way. Checked
      directly afterward: `docs/user/EXAMPLES.md:351` already states
      this exact guidance fully standalone — *"Always give
      `OMP_NUM_THREADS` all your cores — it's the only genuinely
      CPU-bound thread count here..."* — so nothing needed migrating
      for that one. No other such gaps found. Nothing left open here.

- [x] **A5. Update every downstream reference to the 3 moved files:**
  - [x] README.md's "Documentation" table: remove the ARCHITECTURE.md/
        PARALLELISM.md/DESIGN_CPU_GPU_TIMELINE... rows, add one line
        after the table: *"Internal architecture/engineering
        documentation lives under `docs/dev/`."* (mirrors BUILD.md:15's
        existing model: *"...kept for the record in
        [docs/dev/ARCHIVED/](docs/dev/ARCHIVED/)..."*)
  - [x] QUICKSTART.md §7 "Architecture notes" (lines 294-306) —
        currently:
        ```
        - docs/user/ARCHITECTURE.md -- master architecture document...
        - docs/user/PARALLELISM.md -- memory/execution decomposition...
        - docs/user/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md -- ...
        ```
        **Proposed:** replace with a single line pointing to
        `docs/dev/`, retitle the section (it's no longer about
        "architecture notes" once the pointer shrinks to one line) —
        keep the `scripts/plot_tile_async_swimlane.py` line that's
        also in this section (a real, useful tool, unrelated to the
        3 moved docs). This section also contains 2 of the 4 stale
        `io_write_threads` occurrences (lines 298, 404) — see C8,
        fix together.
  - [x] `docs/user/EXAMPLES.md` §4 "Memory tuning" — currently ends
        with a "Full mechanics... are in PARALLELISM.md and
        ARCHITECTURE.md's 'When io_overlap=y can be detrimental'
        section" pointer. Per the standing rule already established
        this session (no `docs/dev/` links from user-facing docs),
        **drop this pointer** rather than reroute it to `docs/dev/`
        — once A4's migration lands, APP_REFERENCE.md is the complete
        user-facing reference.
  - [x] `docs/user/EXAMPLES.md` §5 "I/O parallelism" — same treatment,
        same reasoning.
  - [x] PARALLELISM.md's own line 22-24 (after moving to
        `docs/dev/PARALLELISM.md`) currently reads: *"see
        `docs/user/ARCHITECTURE.md` ('Parallel read/write' and 'Async
        read/write overlap')..."* — two fixes needed together: (a) path
        becomes same-directory `ARCHITECTURE.md` (both now in
        `docs/dev/`), (b) "Parallel read/write" isn't one heading —
        the real headings are two separate ones, `#### Parallel read —
        io_read_threads` and `#### Parallel write — nwriters`.
  - [x] Fix the stale plain-text `planning/` mentions (a directory
        renamed to `docs/dev/` a while back; a prior link-checker pass
        only caught markdown *links*, not plain-text mentions) —
        PARALLELISM.md:468-469, 497-498. Moot for the DESIGN doc since
        it's being archived, but fix if quick while the file is open.

---

## B. Config-doc de-duplication: `cfg/CONFIG_README.md`

Decided already: trim to CONFIG_README.md's own stated scope (parser
mechanics, strict validation rules); APP_REFERENCE.md becomes sole
source of truth for the actual key list. This directly resolves the
resiQ/etc. drift below (C1) at the root instead of patching two copies.

- [x] **B1.** Remove CONFIG_README.md's "Required Keys" section
      (currently lines ~48-70) and any other key-by-key list that
      duplicates APP_REFERENCE.md's rm_synthesis tables — keep only
      the parser-mechanics content (strict validation rules, comma-list
      schema, output-filename convention) that's genuinely about *how
      the format works*, not *what every key does*.
- [x] **B2.** CONFIG_README.md's "## Adding New Config Variables"
      section (contributor instructions: "you must also update
      `src/rm_synthesis_mod.f90`: 1. Add handling in
      `read_cfg_keyval`...") is dev-only content, not user content.
      **Open decision — your call:** relocate to
      `docs/dev/` (e.g. a short CONTRIBUTING-style note)?

---

## C. Factual/accuracy fixes

Each of these has a clear proposed fix; none depend on the structural
moves above. Items marked ✅ were independently re-verified directly
this pass (not just taken from the research agents' reports).

- [x] **C1. ✅ APP_REFERENCE.md:90** — `resiQ`/`slopeQ`/`resiU`/
      `slopeU` described as "only meaningful when `remove_qu_bias=y`"
      (implying optional otherwise). **Code reality**
      (`src/rm_synthesis_mod.f90:2782-2793`): checked unconditionally
      in the same required-key chain as `path`/`infileQ` — missing any
      one is a hard parse error (`status=-143..-146`) regardless of
      `remove_qu_bias`. **Proposed fix:** change the doc to state
      these 4 keys are always required, with the residual/slope
      *values* only being meaningful (i.e. non-zero) when
      `remove_qu_bias=y`.

- [x] **C2. ✅ HIGH — APP_REFERENCE.md's rmclean_cubes section never
      documents `write_clean_diagnostics`** (default `y`,
      `src/rmclean_cubes.f90:1067`, parsed at line 1345-1346) or the 7
      files it writes by default (`src/rmclean_cubes.f90:2953-2989`
      area): `<outfile>.NITER.MAP.FITS`, `.STOP_REASON.MAP.FITS`,
      `.RESID_PEAK.MAP.FITS`, `.RESID_RMS.MAP.FITS`,
      `.N_COMPONENTS.MAP.FITS`, `.COMP_RM_SPREAD.MAP.FITS`,
      `.TOTAL_POL_FLUX.MAP.FITS`. The binary's own `--help` text
      already describes it correctly (`src/rmclean_cubes.f90:1529`:
      *"write_clean_diagnostics (default y): write 7 additional
      per-pixel 2D diagnostic maps..."*) — APP_REFERENCE.md just never
      picked it up. **Every `rmclean_cubes` run today writes 7 files
      the reference doc doesn't mention.** Proposed fix: add the key
      to the Parameters table and the 7 files to the Output files list
      in APP_REFERENCE.md's rmclean_cubes section.

- [x] **C3.** `subim_parfile` is an accepted-but-unused/vestigial
      `rm_synthesis` cfg key (`src/rm_synthesis_mod.f90:2077-2085`,
      default `'subimage.par'`) — parsed, never read anywhere else in
      the codebase, undocumented anywhere. **Open decision — your
      call:** (a) document it as deprecated/no-op in APP_REFERENCE.md,
      (b) leave undocumented (accept a user could set it and see no
      effect, no error), or (c) separately — outside this docs pass —
      consider actually removing it from the parser. Not proposing a
      default here since it touches code, not just docs.
      Document it as deprecated/no-op

- [x] **C4. ✅ QUICKSTART.md:19** — "each lives under `bin/` and is
      also symlinked to `bin/rm_synthesis`". **Code reality**
      (`Makefile:161`): `@cp -f $@ $(EXECUTABLE)` — a plain file copy,
      no `ln -s` anywhere in the Makefile. **Proposed fix:** change
      "symlinked to" → "copied to" (and "the symlink" later in the
      same sentence → "the copy").

- [x] **C5. ✅ QUICKSTART.md:30** — "The GPU binary is built with
      `-ffast-math -DUSE_GPU`." **Code reality**: only true for the
      `GPU_FC=gfortran` path (`Makefile:18`,
      `GPU_GNUFLAGS := ... -ffast-math ... -DUSE_GPU`). The
      default/auto-selected-first path, `nvfortran`
      (`Makefile:17`, `GPU_NVFLAGS := -cpp -O3 -mp=gpu -gpu=cc80,managed
      -DUSE_GPU`), has **no** `-ffast-math`. **Proposed fix:** qualify
      the claim — e.g. "The GPU binary always defines `-DUSE_GPU`;
      the `gfortran`-offload path additionally adds `-ffast-math`
      (the default `nvfortran` path does not)."

- [x] **C6.** QUICKSTART.md's Requirements table (§5, line ~274) lists
      Starlink AST + FFTW3 as needed "only for
      `reproject_cubes`/`convolve_cubes`/`match_cubes`" — omits
      `rmclean_cubes`, which also links FFTW3 (`Makefile:359`).
      BUILD.md already states this correctly ("`rmclean_cubes` needs
      only FFTW3 + CFITSIO... no Starlink AST dependency"). **Proposed
      fix:** add `rmclean_cubes` to that row's tool list (FFTW3 only,
      not AST, if the table can express that distinction — otherwise
      add a footnote).

- [x] **C7.** Same Requirements table's matplotlib/ffmpeg row names
      only `plot_tile_async_swimlane.py`/`animate_fits_cube.py`.
      `plot_rmclean_advisory.py` (referenced by APP_REFERENCE.md
      itself) and `benchmark_omp.py` also import matplotlib.
      **Proposed fix:** broaden the row to name all 4, or reword to
      "various `scripts/*.py` diagnostic/plotting tools."

- [x] **C8. ✅ QUICKSTART.md uses the non-existent key name
      `io_write_threads` in 4 places** (lines 111, 242, 298, 404).
      Confirmed directly: `grep -rn "case ('io_write_threads')" src/`
      → zero hits anywhere in the codebase. The real, current key
      (confirmed at `src/rm_synthesis_mod.f90:2671`,
      `src/rmclean_cubes.f90:1354`, and used correctly everywhere else
      in the docs) is **`nwriters`**. `io_write_threads` only survives
      in historical *source-code comments* describing an old,
      superseded design. **Proposed fix:** replace all 4 occurrences
      with `nwriters`. (Lines 298/404 get fixed as part of A5's
      QUICKSTART §7 rewrite anyway; 111 and 242 need a standalone fix.)

- [x] **C9.** `cfg/rmsynth.cfg` (APP_REFERENCE.md calls it "a full
      annotated template with every key shown in context") is missing
      `lsq_ref_mode`/`lsq_ref_fixed_value` entirely — real keys with a
      whole explanatory paragraph in APP_REFERENCE.md
      (lines 111-112: `zero|mid|centroid|min|max|fixed` semantics and
      the `LSQREF` output header). **Proposed fix:** add both to
      `cfg/rmsynth.cfg`'s "Output format (optional)" section, with
      inline comments matching the file's existing style.

- [x] **C10. ✅ EXAMPLES.md §2a is broken as written.** Confirmed
      directly:
      - EXAMPLES.md:76-78 prose: *"Run against this repo's own
        two-band test fixture (`tests/data/TEST.Q.FITSCUBE` +
        `tests/data/TEST_BAND2.Q.FITSCUBE`, which share the same sky
        grid by construction)"*, then the command
        `bin/rm_synthesis_release_cpu_omp
        cfg/rmsynth-e2e-multiband-smalltest.cfg`, claiming it
        "completes cleanly, no preprocessing step involved."
      - `cfg/rmsynth-e2e-multiband-smalltest.cfg:11` actually reads:
        `infileQ = TEST.Q.FITSCUBE,TEST_BAND2_MISMATCH.Q.FITSCUBE` —
        the *mismatched* fixture, not the matched one the prose names.
      - This is the exact same fixture pair `tests/run_tests.sh`
        section 15 (line ~924) uses specifically to *provoke* the
        `ERROR: RA WCS mismatch for band 2` message shown two sections
        later in EXAMPLES.md §2b.
      - Running §2a's example as written would raise that WCS-mismatch
        error, not "complete cleanly."
      - **No static, checked-in cfg file anywhere in the repo currently
        uses the genuinely-matched `TEST_BAND2.Q.FITSCUBE` pair** —
        `tests/run_tests.sh` section 15 builds that cfg inline via a
        heredoc at test time (lines ~878-898), it's not a committed
        file. Checked this directly this pass — confirmed no existing
        fixture to just repoint to.
      - **Proposed fix (needs your go-ahead, this is new content, not
        a wording tweak):** create a new committed cfg,
        e.g. `cfg/rmsynth-e2e-multiband-matched-smalltest.cfg`, mirroring
        `tests/run_tests.sh`'s own inline `mb_match_cfg` heredoc
        (lines 878-898) but as a real file with repo-relative paths;
        update EXAMPLES.md §2a to reference it; then actually run it
        to confirm it completes cleanly before marking this done (not
        just asserted). 
        Alright, go ahead!

- [x] **C11.** RELEASE_NOTES_1.0.md never mentions the demo GIF /
      `scripts/animate_fits_cube.py`, despite the GIF being generated
      from the exact flagship WALLABY/EMU validation run
      RELEASE_NOTES_1.0.md already describes at length in its
      Validation section. **Proposed fix:** add a short mention/link,
      e.g. in the "Highlights" or "Validation" section, pointing at
      the README's `## Demo` section or the GIF directly.

- [x] **C12. RESOLVED, not a discrepancy.** RELEASE_NOTES_1.0.md's
      "166/166" regression count was unconfirmed by any captured log
      in the repo at audit time (the only saved log predated the test
      section that could produce 166 instead of 165). **Ran the full
      suite fresh this pass (2026-08-28): Total 166, Pass 166, Fail 0,
      Skip 0, Warn 0, RESULT: ALL PASSED.** The doc's claim is
      currently accurate. No action needed — leaving checked as
      verified rather than removing, so the verification is on record.

- [x] **C13.** `cfg/ARCHIVED/README.md`'s category breakdown doesn't
      match the real files on disk: claims Statistics=4/General
      Purpose=40 (lines 18, 21); actual directory contents are
      Statistics=5 (`rmpspec_3C286.cfg`, `rmpspec_3C303.cfg`,
      `rmpspec_casa.cfg`, `rmstat_3C303.cfg`, `rmstat_casa.cfg`) and
      General Purpose=39 — the total of 63 happens to still be
      correct, only the split between those two categories is off by
      one each way. **Proposed fix:** change "4 files" → "5 files" and
      "40 files" → "39 files" at those two lines.

- [x] **C14.** RELEASE_NOTES_1.0.md:148 renders one doc-list entry
      with its raw path as link text
      (`[../../QUICKSTART.md](../../QUICKSTART.md)`) while the 3
      sibling entries just below use clean display text
      (`[TUTORIAL.md](TUTORIAL.md)`, etc.). **Proposed fix:** change to
      `[QUICKSTART.md](../../QUICKSTART.md)` for consistency.

- [x] **C15. Low priority, not a live bug — flagging for awareness
      only.** ARCHITECTURE.md has two headings that would produce the
      same anchor (`## Purpose` and `### Purpose`), and
      APP_REFERENCE.md repeats `### Invocation`/`### Parameters`/
      `### Output files`/`### Example` once per tool (5×). GitHub
      auto-suffixes duplicates (`-1`, `-2`, ...) so nothing currently
      breaks, but a future bare `#invocation`-style link would land on
      the wrong tool's section. No action proposed — just noting for
      future editors. (Moot for ARCHITECTURE.md if A1 is accepted,
      since it becomes a dev-only doc.)

---

## D. Flagged for your decision — not yet resolved

- [x] **D1.** QUICKSTART.md §5 and BUILD.md's own Requirements
      sections give the same `apt-get`/`brew` install instructions for
      the same 3 dependency groups (gfortran, CFITSIO, Starlink
      AST+FFTW3) — genuine duplication, reworded/reformatted (table vs.
      subsections) but the same content maintained in two places.
      **Proposed default:** trim QUICKSTART.md's copy to a short table
      + "see BUILD.md for full install detail," keeping BUILD.md as
      the canonical detailed reference. **Not yet confirmed with
      you** — accept, modify, or reject.
      Accepted!

---

## Confirmed clean (no action needed)

For completeness — these were checked and found correct, so they're
not re-litigated above:

- All markdown-link anchors across all 12 user-facing docs resolve
  correctly, including tricky em-dash-containing headings (e.g.
  EXAMPLES.md §2b/2c/2d). Zero broken links found.
- README.md itself has zero `docs/dev/` leakage already (no directory
  tree, no project-history section, its Documentation table only
  lists the 10 user docs).
- All CLI/cfg parameter tables for `reproject_cubes`, `convolve_cubes`,
  and `match_cubes` in APP_REFERENCE.md match their source parsers
  exactly.
- RM-CLEAN's stopping-criteria docs (`abs_flux_floor`/`auto_nsigma`)
  are accurate, and the old retired keys
  (`threshold=`/`threshold_snr=`/etc.) are genuinely gone from the
  code, not just the docs.
- Every file path and CLI argument referenced in TUTORIAL.md is
  correct and exists.
- EXAMPLES.md's and QUICKSTART.md's own "Contents" lists exactly match
  their real headings — no numbering/title drift.
- No references anywhere in the tracked repo to the abandoned 3D
  PyVista rendering approach.
