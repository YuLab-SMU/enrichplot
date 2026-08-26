# nsea / mnsea Visualization Implementation Checklist

## Goal

Turn the design spec in `nsea-mnsea-visualization-methods-spec.md` into an execution-ready checklist for `enrichplot`.

This checklist has two top-level goals:

1. **Complete the inherited `gseaResult`-style plotting workflow for `nseaResult` and `mnseaResult`.**
2. **Implement a first batch of new mechanism-oriented methods that are specific to network enrichment results.**

## Working Rule

- **First complete the old workflow.**
- **Then add new methods.**
- **Prefer helper-first implementation over plot-first patching.**
- **Each batch should end with tests, NEWS, and a clean commit boundary.**

## Batch 0: Baseline Confirmation

### Goal

Confirm what is already done and what remains missing, so later batches do not re-open closed work.

### Environment Prerequisites

- [x] Verify against source via `pkgload::load_all(".")` or a fresh install of the working tree (the installed Bioconductor `enrichplot` may predate all `mnsea` code; do not trust `library(enrichplot)` alone)
- [x] Install suggested dependencies used by the tested plots, in particular `ggupset` (required by `upsetplot()` tests)

### Checklist

- [x] Reconfirm current `mnseaResult` support in:
  - `dotplot()`
  - `heatplot()`
  - `cnetplot()`
  - `emapplot()`
  - `gseaplot()`
  - `ridgeplot()`
  - `upsetplot()`
  - `ssplot()`
- [x] Reconfirm existing helper contracts:
  - `fortify.mnseaResult()`
  - `fortify_mnsea_contribution()`
  - `fortify_mnsea_subnetwork()`
- [x] Reconfirm current gaps and their exact mechanism:
  - `gseaplot2()` / `gsearank()` / `hplot()` are plain functions calling plain `gsInfo()`; no S3 dispatch, no `layer` support
  - `treeplot()` already has a `gseaResult` S4 method (subclasses dispatch), but `pairwise_termsim` has no `mnseaResult` method
  - `barplot()` has no `gseaResult`-family S3 method, so `nsea` / `mnsea` fall through to `graphics::barplot.default`

### Verification

- [x] `pkgload::load_all(".")` + `testthat::test_file("tests/testthat/test-mnsea-helpers.R")` stays green (with `ggupset` installed)
- [x] No regression in existing `mnsea` helper outputs

## Batch 1: Complete the Running-Score Family

### Goal

Finish the plots that should naturally follow from `gseaplot.mnseaResult()`.

### Prerequisite Refactor

- [x] Turn `gsInfo()` into an S3 generic (`gsInfo <- function(object, geneSetID, ...) UseMethod("gsInfo")`), keep the current body as `gsInfo.gseaResult()`, and let `gsInfo.mnseaResult()` / any future `gsInfo.nseaResult()` dispatch through it
- [x] Keep `get_gsdata()` calling the generic `gsInfo()` so multi-`geneSetID` paths pick up the right method automatically

### Target Functions

- [x] `gseaplot2()` for `nseaResult`
- [x] `gseaplot2()` for `mnseaResult`
- [x] `gsearank()` for `nseaResult`
- [x] `gsearank()` for `mnseaResult`
- [x] `hplot()` for `nseaResult`
- [x] `hplot()` for `mnseaResult`

### Implementation Notes

- [x] Reuse `gsInfo.mnseaResult()` through the new generic (never call it by name inside plot functions)
- [x] Do not invent a second ranked-score pipeline
- [x] For `mnseaResult`, keep the same score-space rule:
  - `layer = NULL` -> `collapsed_scores`
  - `layer = "<single-layer>"` -> `layer_scores[[layer]]`
- [x] Reject ambiguous multi-layer running-score requests unless the plot truly supports them
- [x] Document that `mnseaResult@geneList` is already the collapsed-score ranked list, so the `layer = NULL` default is semantically correct even before the refactor; the refactor is what adds single-layer and stable-ID support

### Likely Files

- [x] `R/gseaplot.R`
- [x] `man/gseaplot.Rd`
- [x] `tests/testthat/test-mnsea-helpers.R`

### Tests

- [x] `gsInfo()` dispatch test: `gsInfo(mnseaObj, id)` and `gsInfo(mnseaObj, id, layer = "rna")`
- [x] `gseaplot2()` smoke test for `nseaResult`
- [x] `gseaplot2()` smoke test for `mnseaResult`
- [x] `gseaplot2()` stable `geneSetID` resolution test
- [x] `gseaplot2()` single-layer score-space test
- [x] `gseaplot2()` multi-`geneSetID` + `layer = NULL` test (define expected collapsed-only behavior)
- [x] `gsearank()` smoke tests for `nsea` and `mnsea`
- [x] `hplot()` smoke tests for `nsea` and `mnsea` (implemented with base `ggplot2`, no `ggHoriPlot` required)
- [x] boundary test for invalid `layer`

### Exit Criteria

- [x] Running-score family methods are usable end-to-end for both `nseaResult` and `mnseaResult`
- [x] No duplicated score-selection logic remains
- [x] `gsInfo()` has exactly one dispatch path, and all running-score plots consume it

## Batch 2: Complete the Similarity / Tree Family

### Goal

Finish the plots that depend on pathway similarity and term clustering.

### Prerequisite: publicize `mnsea` pairwise similarity

- [x] Add / stabilize a `pairwise_termsim` method or equivalent public helper for `mnseaResult` so `treeplot()`, `emapplot()`, and `ssplot()` share one layer-aware similarity definition
- [x] Do not create a second, incompatible similarity definition for `mnsea`

### Target Functions

- [x] `treeplot()` for `nseaResult` (confirm existing `gseaResult` S4 dispatch + add smoke tests; no new method required)
- [x] `treeplot()` for `mnseaResult` (add layer-aware similarity semantics)
- [x] `barplot()` support for `nseaResult`
- [x] `barplot()` support for `mnseaResult`

### Implementation Notes

- [x] `treeplot()` already dispatches to `signature(x = "gseaResult")` for both subclasses; the `nsea` work item is tests/documentation, not a new method
- [x] Reuse current `emapplot` / `ssplot` similarity logic only as a stopgap until the public `pairwise_termsim` path lands
- [x] Add a `barplot.gseaResult` (or `nsea` / `mnsea`) S3 method / alias because the current `barplot.enrichResult` never matches `gseaResult` subclasses; this is a method gap, not a “compatibility check”

### Likely Files

- [x] `R/paired-similarity.R` or `R/pairwise_termsim.R` (similarity publicization)
- [x] `R/treeplot.R`
- [x] `R/barplot.R`
- [x] `man/treeplot.Rd`
- [x] `man/barplot.Rd`
- [x] `tests/testthat/test-mnsea-helpers.R`

### Tests

- [x] `treeplot()` smoke test for `nsea`
- [x] `treeplot()` smoke test for `mnsea`
- [x] `treeplot()` no-precomputed-termsim path for `mnsea`
- [x] `pairwise_termsim()` returns term × term structure for `mnsea` (the shared similarity contract)
- [x] single-pathway / two-pathway boundary tests
- [x] `barplot()` smoke tests for `nsea` and `mnsea` (prove they no longer hit `graphics::barplot.default`)

### Exit Criteria

- [x] Existing `gseaResult`-style workflow is effectively complete for `nseaResult` / `mnseaResult`
- [x] Remaining gaps are no longer “missing old plots”, only “new method work”
- [x] `treeplot()` / `emapplot()` / `ssplot()` consume the same public `mnsea` similarity helper

## Batch 3: Stabilize Mechanism-Oriented Helper Layer

### Goal

Create the reusable summary helpers needed by new methods.

### Target Helpers

- [x] `compute_rewiring_score()`
- [x] `classify_mechanism_state()`
- [x] `summarize_nsea_mechanism()`
- [x] `extract_rewiring_features()`

### Minimal Contracts

#### `compute_rewiring_score()`

- [x] Accept `nseaResult` or `mnseaResult`; for cross-network / cross-condition comparisons, also accept a named list and return a context column
- [x] Return one numeric rewiring score per pathway
- [x] Use the agreed transparent first version: `rewiring_score = 1 - leading-edge Jaccard overlap`, range `[0, 1]` (`0` = conserved, `1` = fully rewired); edge overlap / centrality drift are optional later composite terms behind the same function
- [x] For a single `nseaResult` without a reference, return `NA_real_` (or error) rather than a fake 0

#### `classify_mechanism_state()`

- [x] Use enrichment shift + rewiring score
- [x] Use written defaults (and an optional `thresholds` argument) in the first version; do not silently auto-tune thresholds
- [x] Return one of:
  - `conserved`
  - `rewired`
  - `context_specific`
  - `contradictory`

#### `summarize_nsea_mechanism()`

- [x] Return a term-level summary table with at least:
  - `ID`
  - `Description`
  - `NES`
  - `p.adjust`
  - `leading_edge_size`
  - `leading_edge_overlap`
  - `rewiring_score`
  - `centrality_shift`
  - `mechanism_class`
- [x] Document explicitly which columns are `NA` for a single `nseaResult` (no reference context) and which are computed from `mnsea` layer vs collapsed comparisons
- [x] Never fill missing rewiring/centrality values with 0; use `NA` + explicit error/argument messaging

#### `extract_rewiring_features()`

- [x] Require `layer` / `reference_layer` (or a second reference object) instead of guessing a default reference
- [x] Return a feature-level comparison table with at least:
  - `Feature`
  - `score`
  - `abs_score`
  - `sign`
  - `status`
- [x] `status` takes one of `shared` / `gained` / `lost` / `shifted`; error when no reference is supplied instead of returning all-`shared`

### Likely Files

- [x] `R/nsea-mechanism-helpers.R` or equivalent new helper file
- [x] roxygen docs for helper contracts
- [x] `tests/testthat/test-nsea-mechanism-helpers.R`

### Tests

- [x] contract tests for each helper
- [x] deterministic classification tests
- [x] `compute_rewiring_score()` single-`nseaResult` no-reference behavior (NA / error)
- [x] `extract_rewiring_features()` missing-reference error test
- [x] empty / single-term / single-layer boundary tests
- [x] invalid-layer and invalid-pathway tests

### Exit Criteria

- [x] All new mechanism plots can consume the same helper layer
- [x] Helper outputs are explicit enough to test without plotting

## Batch 4: Implement `phaseplot()`

### Goal

Deliver the first new method with clear `nsea` / `mnsea` identity.

### Method Definition

- [x] X-axis = enrichment shift
- [x] Y-axis = rewiring score
- [x] Size = leading-edge size or overlap
- [x] Color = significance or mechanism class

### Work Items

- [x] Add generic to `R/AllGenerics.R`
- [x] Create `R/phaseplot.R`
- [x] Define methods for:
  - `nseaResult`
  - `mnseaResult`
- [x] Reuse `summarize_nsea_mechanism()`
- [x] Make default axis labels explicit and readable

### Likely Files

- [x] `R/AllGenerics.R`
- [x] `R/phaseplot.R`
- [x] `man/phaseplot.Rd`
- [x] `tests/testthat/test-phaseplot.R`

### Tests

- [x] smoke test for `nseaResult`
- [x] smoke test for `mnseaResult`
- [x] mechanism class mapping test
- [x] size/color semantic tests
- [x] empty result and one-term boundary tests

### Exit Criteria

- [x] `phaseplot()` provides information not already available from `dotplot()` / `emapplot()`
- [x] The default plot already separates conserved vs rewired patterns in a readable way

### Enhancements (2026-08-23)

- [x] `phaseplot()` accepts `reference` to compute real `delta_NES`
- [x] `phaseplot()` supports `x_axis = "delta_NES"` / `"NES"` and `size_var = "leading_edge_size"` / `"leading_edge_overlap"`
- [x] `phaseplot()` errors clearly when `delta_NES` is requested without a reference

## Batch 5: Implement `rewireplot()`

### Goal

Deliver the first pathway-specific mechanism explanation plot.

### Scope (decided)

- [x] First version supports `mnseaResult` only; `nseaResult` requires an explicit second result object / reference and is deferred
- [x] Require `reference_layer` (no silent default) for `mnseaResult`

### Method Definition

- [x] Nodes = leading-edge or pathway-driving features
- [x] Edges = pathway-specific subnetwork
- [x] Feature status displayed as:
  - `shared`
  - `gained`
  - `lost`
  - `shifted`
- [x] Optional display of coupling-mediated edges for `mnsea`

### Work Items

- [x] Add generic to `R/AllGenerics.R`
- [x] Create `R/rewireplot.R`
- [x] Reuse:
  - `fortify_mnsea_subnetwork()`
  - `extract_rewiring_features()`
- [x] Define stable pathway selection rules
- [x] Decide whether to facet by layer or color by layer for the first version
- [x] Emit an explicit error when called on a bare `nseaResult` without a reference, instead of producing an all-`shared` plot

### Likely Files

- [x] `R/AllGenerics.R`
- [x] `R/rewireplot.R`
- [x] `man/rewireplot.Rd`
- [x] `tests/testthat/test-rewireplot.R`

### Tests

- [x] smoke test for `mnseaResult`
- [x] missing-`reference_layer` error test
- [x] stable `pathway_id` resolution test
- [x] feature-status mapping test
- [x] no-edge / no-coupling boundary tests
- [x] bare `nseaResult` without reference error test

### Exit Criteria

- [x] `rewireplot()` answers “same pathway name, same mechanism or not?”
- [x] pathway-specific rewiring evidence is readable without extra manual preprocessing

## Batch 6: Implement `consensusmap()`

### Goal

Add a multi-network / multi-layer overview plot for mechanism agreement and disagreement.

### Input Contract (decide first)

- [x] Accept a named list of results (`list(networkA = resA, networkB = resB, ...)` or `list(conditionA = resA, ...)`) as the comparison input
- [x] For a single `mnseaResult`, allow layer-as-context fallback; for a single `nseaResult`, error and ask for a list

### Work Items

- [x] Add generic and method file
- [x] Build a term × context summary matrix
- [x] Show enrichment strength and topology consistency together
- [x] Attach one mechanism class per pathway

### Dependencies

- [x] `summarize_nsea_mechanism()`
- [x] `classify_mechanism_state()`

### Exit Criteria

- [x] Users can quickly identify conserved, rewired, and context-specific pathways
- [x] Single-object inputs fail loudly instead of producing a one-column plot

### Enhancements (2026-08-23)

- [x] `consensusmap()` exposes both enrichment strength (`fill`) and topology consistency (`size`)
- [x] `consensusmap()` accepts `fill_var` / `size_var` / `label` / `reference`
- [x] `consensusmap()` uses the first list element as default reference for `delta_NES`
- [x] `compute_rewiring_score()` supports cross-object `reference` for real rewiring scores

## Batch 7: Implement `mechanismflow()`

### Goal

Add an evolution-style view for pathway state transitions across conditions or layers.

### Input Contract

- [x] Transitions come from a named list of results or a `mnseaResult` layer sequence; there is no implied time order on a single object
- [x] First version focuses on `mnseaResult` (layer sequence); `nseaResult` needs an explicit named list

### Work Items

- [x] Add generic and method file
- [x] Define state transitions between contexts
- [x] Choose a first rendering strategy:
  - sankey
  - river
- [x] Keep the first version focused on `mnseaResult` if needed

### Exit Criteria

- [x] The plot makes pathway state transitions easier to read than side-by-side `NES` comparisons

### Enhancements (2026-08-23)

- [x] `mechanismflow()` accepts `reference` and `flow_var`
- [x] `mechanismflow()` uses flow magnitude (NES / delta NES / leading-edge size) for line width and point size
- [x] `mechanismflow()` uses a stable mechanism-state ordering on the y axis

## Shared Verification Rules

### For Every Batch

- [x] Add or update focused tests
- [x] Run targeted verification first
- [x] Update `NEWS.md`
- [x] Regenerate docs if roxygen changes
- [ ] Keep commit boundaries narrow and readable

### Test Priority

- [ ] contract tests before visual tests
- [ ] helper tests before plot tests
- [ ] edge-case tests before “pretty plot” checks

## Suggested Commit Boundaries

### Commit Group A: Complete Old Plots

- [x] `refactor: turn gsInfo into an S3 generic with layer-aware dispatch`
- [x] `feat: add nsea gseaplot2 family support`
- [x] `feat: add mnsea gseaplot2 family support`
- [x] `feat: publicize mnsea pairwise similarity for treeplot`
- [x] `feat: add nsea mnsea treeplot support`
- [x] `feat: add gseaResult-family barplot methods for nsea and mnsea`

### Commit Group B: Add New Helper Layer

- [x] `feat: add nsea mechanism helper summaries`

### Commit Group C: Add New Methods

- [x] `feat: add phaseplot for nsea and mnsea`
- [x] `feat: add rewireplot for mnsea`
- [x] `feat: add consensusmap for network mechanism overview`
- [x] `feat: add mechanismflow for pathway state transitions`

## Recommended Immediate Next Step

If implementation resumes now, the next coding batch should be:

0. **set up `pkgload` / `ggupset` verification baseline**
1. **refactor `gsInfo()` to dispatch, then `gseaplot2()`**
2. **`gsearank()` / `hplot()`**
3. **publicize `pairwise_termsim` `mnsea` similarity, then `treeplot()`**
4. **add `barplot` family methods**

Only after that should the work move to:

5. **`phaseplot()`**
6. **`rewireplot()`**

This keeps the development path clean:

- first complete compatibility
- then establish new method identity
