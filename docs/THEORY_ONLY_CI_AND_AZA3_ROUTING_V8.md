# anto v8 — Theory-only active CI; Cirsium empirical code lives in aza3

**2026-10-19.** This is a repository-boundary cleanup, not a new biological result.

## User-directed ownership is now enforced by CI

- **anto:** cross-angiosperm mathematical questions: movement order, history-dependent latent availability, sharp higher-order resource probability bounds and generic stochastic colonization nulls.
- **aza3:** all *Cirsium*-specific biological source intake, public transcript and annotation mapping, DFR/ANS sequence homology triage, ancestor/introgression comparisons, same-individual floret/involucre/head modules and ecological fitness study plans.

The source and experimental transfer was established in [aza3 PR #45](https://github.com/zuizui0223/aza3/pull/45). The operational source analysis is now being consolidated into [aza3](https://github.com/zuizui0223/aza3/blob/main/analysis/run_anto_cirsium_empirical_sources_v3.py), including accession-to-specimen mapping, source DFR/ANS candidate family triage and organ evidence. A later PR or pending branch may change the link's reachability during integration; the reference owner remains aza3.

## CI now executes exclusively generic models

The remaining active workflow .github/workflows/synthetic-validation.yml runs:

- three-state generic W+/W−/P lineage CTMC and conditional first-entry tests;
- generic seed colonization, pollen/introgression as independent mechanism states and negative controls;
- exact first/pairwise/higher-order module inventory accessibility bounds;
- module-loss-by-region migration timing counterfactuals;
- labelled synthetic history discrimination.

**Not executed in anto CI anymore:** public Cirsium SRA ingest, Ryukyu young-leaf RNA source QC, 2020 organ expression data parsing, DFR/ANS source-transcript sequence analysis, Cirsium molecular reference quality tests, source supplement downloading, or laboratory capability/fitness claims. Those operations belong to the aza3 research workflow.

Legacy anto docs, source programs and source datasets under historical paths remain visible and versioned solely to audit the route by which the general theory originated. They must not be silently used as a second active empirical authority. Do not delete existing provenance to fake a new scientific discovery.

## Consequence for universal inference

The 3-resource sharp bound of 0–25% given identical first/pairwise statistics is a probability-theory constraint, **not** a universal observation in flowering plants. Species-specific paths (structural pathway preservation, regulatory loss, reactivation, immigration, introgression and new mutation) still require independent molecular, genealogy and developmental measurements, and cross-clade held-out replication.

An original research programme can still use *Cirsium* as the initial mechanistic example through aza3, then compare validated independent angiosperm clades under anto. Until several independent clades have genuine molecular pathway-competence measurements and historical source orientation, no universal evolutionary accessibility law is empirically established.

## Reproducibility

    python -m pip install numpy scipy
    python -m unittest discover -s tests -p 'test_model.py' -v
    python -m unittest discover -s tests -p 'test_colonization.py' -v
    python -m unittest discover -s tests -p 'test_joint_accessibility_v7.py' -v
    python scripts/run_joint_accessibility_v7.py

The validated biological-source-to-gene pipeline's CLI is provided under aza3, not required for anto CI.
