# ANTO v44 — Published 2014 Phlox C++ selection, pollen dispersal and assortative mating sensitivities

2026-11-17. This is a source-compiled deterministic sensitivity audit around one original-author example-input START vector, NOT a newly fitted ecological model, selection parameter, MCMC posterior or original natural evolutionary experiment.

## Authentic original source

Hopkins et al (2014), Current Biology, DOI 10.1016/j.cub.2014.07.027. Original PhloxHopkins/PhloxSelection source commit 60ff508170f1be592bbc6fe7ed46ac6e5b459026. Original 20 Git blobs plus immutable v42/v43 source/results rechecked. Eight original 2014 forward-model C++ files unchanged. ANTO's v43-authored diagnostic driver modified only by strict predefined source insertion for 7 scenarios. The original C++ adult flower-frequency and family/offspring probability functions and calcTotLike were used. Original program files were NOT redistributed.

Original-author example input: 32 populations, 11,379 adults, 371 mother-family rows, 2,496 offspring. Main publication states 12,024 / 382 / 2,581; gaps 645 / 11 / 85 are UNRESOLVED. Initial author vector: 5 transect borders [.279,.362,.612,.639,.362], assortative mating a=0, pollen geometric kernel spread m=.72, original zone fitness DB/DR/LR [.88,.91,.90] and [3.7,6.3,5.5] relative to LB=1. These are **initial model settings**, not source author's fitted posterior.

## Original source C++ actual results

For every admitted scenario, Delta log-likelihood = original source likelihood of scenario minus original author's UNOPTIMIZED start likelihood. Positive means higher likelihood for this specific source example and given vector, NOT that the scenario is scientifically favored under fitted conditions. No p-values, confidence bounds, AIC or BIC are assigned.

| Source intervention | Adult Delta logL | Offspring Delta logL | Total Delta logL | Status |
| --- | ---: | ---: | ---: | --- |
| baseline (exact original start) | 0 | 0 | 0 | admitted; exactly v43 |
| all DB/DR/LR fitness neutral to 1 | -9,655.426 | -994.213 | **-10,649.639** | admitted; all five transects total negative |
| phenotype-level epistasis set 0 | -261.291 | **+4.474** | **-256.817** | admitted; two of five transects total positive |
| pollen geometric spread m .72 to .36 | -432.213 | **+8.674** | **-423.539** | admitted; two of five transects total positive |
| pollen geometric spread m .72 to 1.44 | -437.574 | -89.802 | **-527.376** | admitted; all five transects total negative |
| same relative fitness in both habitat zones | — | — | — | BLOCKED; source model final change exceeded original convergence tolerance |
| same-flower-phenotype assortment a=0 to .25 | — | — | — | BLOCKED; source model generated a non-normalized / invalid maternal-color probability vector |

The 5 admitted cases include the baseline; **2 of 7 precommitted scenarios were technically blocked**. Do NOT replace blocked with zero model effect or reject their biological mechanisms. Both blocked conditions remain in frozen derived result. The original source's source-coupled genotype-frequency convergence was checked explicitly, not inferred only from iterations below maxGens. All admissible variants had normalized phenotype and conditional-family vectors, original source total likelihood equal separately summed adult + conditional-family original MakeLikelihood components, and exactly identical output across 2 fresh runs.

**Most consequential biological model-side observation:** adult-colour and offspring-family likelihood components may react in opposite directions. At these one-factor source START perturbations, removing phenotype fitness epistasis or reducing geometric pollen dispersal slightly raised conditional offspring likelihood, while reducing adult source likelihood. It would be invalid to infer that offspring empirically favor reduced migration or a zero-epistasis natural genotype: no model optimization, uncertainty, population-family holdouts or independent field data were used. The adult observations (11,379) and children (2,496) have different denominators, so raw component log-likelihood sizes alone do not rank ecological contributions.

## Exact technical provenance, failures and repair

Unchanged original World / Deme / selection / geometric pollen kernel / mating / likelihood sources were compiled with a temporary narrow compatibility shim for the original missing nr3.h (VecDoub aliases and old unqualified std::cerr), and only 299 completely whitespace-only padding lines were removed from a temporary copy of the original fourth family file. Original raw author files and all nonblank record bytes were unchanged. No mcmc.cpp / main.cpp / missing Rseeds.txt was compiled or invented. No 75,000-step original MCMC run; this is original model forward simulation at seven hypothetical fixed parameter points.

Initial v44 original-source run returned a valid complete source scenario derived receipt, but the subsequent CI verification step referred to an incorrectly prefixed JSON key, causing a failed check **after** source calculation. The workflow key reference was corrected without changing author data, original C++, scenario values, scoring or biological results. Source derived JSON was frozen from all six complete, ordered source-result chunks; the final source CI must rebuild it identically. Neither the initial workflow failure nor a blocked model setting means zero natural selection.

## Files and next source-backed step

- Pre-readout protocol: docs/ANTO_V44_ORIGINAL_CPP_ABLATION_PROTOCOL_20261117.md
- Seven fixed original-C++ vectors: data/external/anto_v44_source_cpp_fixed_scenario_vectors_20261117.json
- Actual source forward execution: scripts/analyze_Phlox_original_CPP_source_ablations_v44.py
- 19 new adversarial numerical and source tests: tests/test_Phlox_original_CPP_source_ablation_v44.py
- Full source-derived frozen JSON: results/anto_v44_original_CPP_seven_scenarios_likelihood_20261117.json
- CI workflow: .github/workflows/anto-v44-original-cpp-scenarios.yml

This is the first same-source and same-input cross-mechanism fixed-point *model sensitivity* in this ANTO branch. It is NOT an independently fitted source replication of original 2014 spatial selection or proof of natural ancestral pigmentation re-acquisition. A legitimate next step will need constrained parameter optimization, source cohort reconciliation and separate out-of-sample checks, with blocked scenario numerical behavior diagnosed before attempting quantitative inference. Keep Draft PR stacked on v43 PR #53; no merge before substantive review.
