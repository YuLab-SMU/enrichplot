# AGENTS.md — enrichplot

Guidance for AI agents (and humans) working in this repository. Last verified 2026-09-21.

## What this package is

enrichplot is the ggplot2-based visualization layer of the clusterProfiler ecosystem
(Bioconductor). It draws results produced by ORA/GSEA tools: `enrichResult`,
`gseaResult`, `compareClusterResult` (classes from DOSE), and `nseaResult`,
`mnseaResult` (classes from enrichit). It contains almost no statistics — its job is
`S4 object -> fortify() -> data.frame -> ggplot`.

Key runtime deps: ggplot2 (>= 3.5, tested against 4.0.x), ggnewscale (multi-panel
color scales), ggtangle (cnetplot edges), aplot (patchwork-style gseaplot2),
ggtree (treeplot), igraph (emapplot graphs), yulab.utils (error handling:
`yulab_abort`/`yulab_warn`, `parse_ratio`), enrichit (classes + gene/term accessors
like `geneInCategory()`).

## Data flow and the fortify contract

1. Plot methods call `fortify(object, showCategory=..., by=..., split=...)` to get a
   plain data.frame. All column-name assumptions downstream come from here.
2. `fortify_internal()` in R/method-fortify.R is the shared core. For
   `compareClusterResult` with `by=`:
   - `"geneRatio"`: character ratios `"k/n"` are parsed to numeric; Cluster labels are
     relabeled to `"<Cluster> \n (<n>)"` (cluster size) unless multiple ONTOLOGY values
     are present; `"count"`: `parse_ratio` only; `"rowPercentage"`: produces a
     `Percentage` column.
3. Term selection is centralized in R/data_utils.R (`update_n()`, `select_terms()`,
   `get_term_mapping()`, `resolve_term_rows()`). `showCategory` accepts a number OR a
   vector of term IDs/labels/Descriptions.

## Term identity vs. display labels (easy to break)

- Display labels come from `get_term_mapping()`: Description, with duplicates
  disambiguated as `"Description [ID]"`.
- In the current devel line, `pairwise_termsim()`/`get_similarity_matrix()` key the
  termsim matrix by **term labels** (Descriptions), not IDs, for ALL methods (JC and
  Wang/Resnik alike). Do not "helpfully" remap label-keyed rows back to IDs — that
  exact bug (NA edges -> igraph "edge data frame contains NAs") was fixed once.
- Internally, selection uses stable IDs; labels are attached at the last moment.

## Repository map

| Area | Files |
|---|---|
| Shared core | R/method-fortify.R, R/data_utils.R, R/AllGenerics.R, R/color_utils.R (`set_enrichplot_color`), R/plot_utils.R (`plotting.clusterProfile`), R/utils.R |
| Base plots | dotplot.R, barplot.R, heatplot.R, volcano/manhattan/volplot.R, densityplot.R |
| Network plots | cnetplot.R, emapplot.R, emapplot_utilities.R, ssplot.R, goplot.R, treeplot.R, pairwise_termsim.R |
| GSEA views | gseaplot.R, ridgeplot.R, upsetplot.R, wordcloud.R |
| Mechanism plots (nsea/mnsea) | phaseplot.R, rewireplot.R, consensusmap.R, mechanismflow.R, mnsea-helpers.R, nsea-mechanism-helpers.R |

## Interface rules for rebuilt features

Some historical features were implemented in the old `ggraph`-based code and may be
partially missing or intentionally reshaped after the move to `ggtangle` and the
current plotting architecture. When restoring a genuinely useful missing feature, do
**not** treat "feature parity" as "copy the old interface verbatim".

Use these rules instead:

1. **Simple beats exhaustive.** The default path should solve the common use case with
   as few arguments as possible. Extra control is only worth adding when it improves
   real user workflows without making the main path harder to follow.
2. **Follow the current API, not the historical one.** Reintroduced functionality must
   match the naming, argument style, and semantics of today's `enrichplot` functions.
   Do not revive old argument bundles, legacy flag combinations, or backend-shaped
   parameters just because they existed before.
3. **One concept, one parameter.** Prefer a small number of orthogonal arguments over
   multiple overlapping switches. Avoid designs where users must coordinate several
   booleans or memorize hidden precedence rules.
4. **User-facing arguments describe plot semantics, not implementation details.**
   Expose concepts such as term selection, grouping, labels, colors, and panel type.
   Do not leak `ggraph`/`ggtangle`/`igraph` internals into the public API unless there
   is a compelling and user-comprehensible reason.
5. **Prefer current shared vocabulary.** Reuse the modern parameter conventions already
   present in the package (`showCategory`, `group`, `group_legend`, normalized measure
   names such as `Count`/`GeneRatio`/`Percentage`, etc.) instead of inventing
   plot-specific synonyms.
6. **Behavioral parity matters more than signature parity.** If an old feature is worth
   bringing back, preserve the useful outcome in a cleaner form rather than matching
   every historical argument name or edge-case behavior.
7. **Make complexity pay rent.** If a restored control would only serve a narrow or
   confusing legacy workflow, prefer leaving it out or folding it into a simpler,
   more general option. Missing low-value complexity is acceptable; user-hostile API
   shape is not.
8. **Test the user story, not just the code path.** Regression tests for rebuilt
   features should assert the intended user-visible behavior under the new interface,
   not merely that the legacy internal route can still be reached.

## Testing

- `testthat`, edition 3. Run: `Rscript -e 'setwd("<repo>"); devtools::load_all("."); testthat::test_local(".")'`.
- Mock objects live in tests/testthat/helper-mock-results.R
  (`mock_enrich_result()`, `mock_gsea_result()` via test-plotting-regression.R,
  `mock_comparecluster_result()`, `mock_nsea_result()`, `mock_mnsea_result()`).
  Extend these instead of hitting the internet; real-data tests use
  `skip_if_not_installed(...)`.
- **Always evaluate plots in tests**: `ggplot()` construction succeeding proves
  nothing. `geom_bar()` + mapped y, missing columns in a layered aes, and discrete
  data on a continuous scale all only explode at draw time. Standard helper:

  ```r
  expect_ggplot <- function(p) {
      expect_s3_class(p, "ggplot")
      expect_error(ggplot2::ggplot_build(p), NA)
  }
  ```
- Regression gate: every fix must keep the full suite green before commit.

## Development environment gotchas

- Test devel code with `devtools::load_all()`. To compare against the Bioc RELEASE
  (1.32.x), a built library exists at `/tmp/rlib-ep132`; set
  `.libPaths(c("/tmp/rlib-ep132", .libPaths()))` **before** any `library()` call,
  and print `packageVersion("enrichplot")` + the lib path in verification scripts so
  you know which tree actually ran.
- `R CMD build` has **no `-o` flag**; run it from the target dir with the repo path
  as argument. `R CMD check` with `_R_CHECK_FORCE_SUGGESTS_=false` to skip missing
  Suggests. The only expected WARNINGs in a dev checkout are the two
  "vignettes but no inst/doc" ones.
- Shell working directory is not reliable across tool invocations — put
  `setwd("<abs path>")` INSIDE `Rscript -e` scripts. Beware stale trees from earlier
  sessions: `/tmp/enrichplot_orig/`, `/tmp/epcheck/enrichplot.Rcheck/` (a build
  snapshot). Before running tests, confirm which tree you are in.

## Maintenance commands (Makefile)

Routine package upkeep is scripted in the repo `Makefile` — use these targets
instead of re-deriving commands (they encode the maintainer's workflow):

- Docs and checks: `make rd` (roxygenise), `make all` (= `rd check clean`),
  `make check` (`devtools::check()`; needs Suggests installed — if Suggests are
  missing, fall back to raw `R CMD check` with `_R_CHECK_FORCE_SUGGESTS_=false`),
  `make check2` (`R CMD check` on the built tarball), `make check-dontrun`
  (checks with `\dontrun{}` examples executed), `make bioccheck`.
- Build/install: `make build` (`devtools::build()`), `make build2`
  (`R CMD build --no-build-vignettes` from the parent dir), `make install`.
- Release prep: `make for-release` (= `rd check-dontrun clean readme`; also renders
  README.Rmd). `BIOCVER := RELEASE_3_23` at the top pins the release branch — bump
  it each April/October release cycle.
- Git workflow (Bioconductor two-remote setup): `make update` syncs local devel
  from `upstream/devel` + `origin/devel`; `make release` checks out the pinned
  release branch; `make push` pushes devel to BOTH `upstream` (Bioconductor git,
  shared state) and `origin` (GitHub) — do not run without explicit instruction;
  `make biocinit` adds the upstream remote. `make rmrelease` deletes the local
  release branch (destructive).

## Hard-won bug patterns (all of these were real bugs here)

1. **Scalar boolean predicates.** `if (x && grep("/", x[1]))` crashes with
   "missing value where TRUE/FALSE needed" when `grep` returns `integer(0)` or `NA`.
   Use `isTRUE(grepl(...))`. Same for any `if()` over data-derived conditions: guard
   length-0 and NA (`isTRUE()`, `!is.null() && length() > 0`).
2. **`linewidth`, not `size`, on line geoms** (geom_line/segment/path/hline/vline).
   ggplot2 deprecates `size` for lines. Two traps:
   - An inherited global `size` aesthetic (point-size scaling) leaks into line layers
     and triggers the same deprecation. `aes(size = NULL)` does **not** suppress it on
     ggplot2 4.0.3 — the fix is an explicit full aes with `inherit.aes = FALSE`.
   - `ggnewscale::new_scale_*()` internally runs `ggplot_build()`, so layer warnings
     surface during `+` itself, not at final draw — don't misread the stack trace.
3. **`geom_col()`, not `geom_bar()`, whenever `y` is mapped.** `geom_bar()` is
   stat_count and errors at build time with a mapped y.
4. **Misspelled arguments vanish into `...`** (e.g. `stringAsFactors=` was silently
   ignored for months). When adding params, grep all call sites; when reviewing,
   check arg spellings against the callee signature.
5. **S4 objects + dplyr verbs cross package boundaries.** `filter()`/`group_by()` on a
   `compareClusterResult` only work because clusterProfiler registers
   `filter.compareClusterResult` etc. Fine for user-facing code (users have
   clusterProfiler loaded), but don't rely on it in enrichplot-internal helpers —
   operate on the fortified data.frame or `object@compareClusterResult`.
6. **Column contracts differ per input class.** ORA results have character
   `"k/n"` ratios and no FoldEnrichment; GSEA results have numeric NES and no
   GeneRatio; compareClusterResult adds Cluster; mnseaResult adds layer columns.
   Anything that indexes a column must either be class-specific or guard existence.
7. **`showCategory` can be a number or a character vector** of IDs/labels/Descriptions
   — route through `update_n()`/`resolve_term_rows()`, never `x[seq_len(n)]` directly.
8. **Normalize user-facing measure aliases before aesthetic mapping.**
   `count`/`Count`, `geneRatio`/`GeneRatio`, and `rowPercentage`/`Percentage`
   should be canonicalized once before they reach `.data[[...]]`; otherwise
   the error only appears at `ggplot_build()` time as a missing-column failure.
9. **`dotplot2()` cannot assume `FoldEnrichment` already exists.** Many
   compareCluster ORA results only carry `GeneRatio` and `BgRatio`; derive
   `FoldEnrichment = GeneRatio / BgRatio` when possible, otherwise fail with a
   direct error instead of a `[[<-` replacement-length crash.

## Conventions

- roxygen with markdown; imports declared per-function via `@importFrom` at each
  site (package avoids NAMESPACE-wide imports).
- NEWS.md entries go under the current devel version heading, one `+` bullet per
  user-visible change, dated `(YYYY-MM-DD, Weekday)`.
- When adding back a missing capability from the pre-`ggtangle` era, prefer a small,
  user-friendly interface that matches the current package style over exact historical
  argument compatibility.
- Commit style: `fix:`, `test:`, `docs:` prefixes seen in history. Version bumps in
  DESCRIPTION are made by the maintainer, not agents — ask before committing one.
- Bioc timeline: release 1.32.0 = Bioc 3.23; devel is 1.99.x heading to 2.0.0 with
  Bioc 3.24 (code freeze mid-October 2026). Behavior changes must keep the release
  line's documented interface working.
- R style: 4-space indent, 80-ish columns, `snake_case` internals, `NULL` defaults
  over `missing()`.

## Related workspaces

- Nature Methods Correspondence revision materials (response letter, review
  comments) live outside this repo at
  `/home/wang/data/manuscript-in-preparation/enrichplot/revision/`. Letter claims
  about fixes must be verifiable against this repo's tests.
- The tutorial book source is at
  `/home/wang/data/source/clusterProfiler_family/biomedical-knowledge-mining-book`
  and may be edited when a tutorial example is wrong.
