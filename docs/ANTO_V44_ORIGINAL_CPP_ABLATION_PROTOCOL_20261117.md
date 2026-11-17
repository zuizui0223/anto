# ANTO v44 — Source-authenticated Phlox C++ forward-model component ablations

**2026-11-17. Frozen BEFORE any v44 ablation results.** An exploratory deterministic perturbation of the ORIGINAL 2014 author model and EXAMPLE input, **not** model re-fitting, a selection coefficient, or a new experimental data set.

## Original frozen inputs and author model
Hopkins et al (2014), Current Biology DOI 10.1016/j.cub.2014.07.027; original author PhloxHopkins/PhloxSelection SHA **60ff508170f1be592bbc6fe7ed46ac6e5b459026**, 20 exact Git blob hashes source-checked in ANTO v42. Baseline v43 actual original forward result must exactly reproduce 5-transect total logL **−9867.530126894539**, adult logL **−8098.322108745898**, family logL **−1769.2080181486392**, plus per-transect logL to 1e−7. Original author README says **example input**: 32 adult populations / 11,379 adult colours / 371 maternal families / 2,496 offspring. Main paper 12,024/382/2,581 difference 645/11/85 unresolved. DO NOT insert fabricated records or imply full original 75,000-step MCMC.

Fixed 13-vector: borders .279, .362, .612, .639, .362; same-phenotype mating a=0; migration/dispersal KERNEL parameter m=.72; first habitat zone relative DB/DR/LR fitness .88/.91/.90 and second 3.7/6.3/5.5 (LB=1). Original author Params::addUntransformedPars has a at vector index 5, m at 6, fitness at indices 7–12. Original source uses geometric tenGBinMigRate(m, spacing .5 km), not proportion of migrants.

## Exactly seven source-vector scenarios
1. baseline: original vector unchanged
2. fitness_neutral: all 6 DB/DR/LR weights set to 1.0 (LB remains 1)
3. fitness_uniform_zone0: zone 2's 3 fitness values equal zone 1's original .88/.91/.90 (no between-habitat relative fitness contrast)
4. epistasis_zero: fix each zone DB/LR and LB=1; reset DR=DB+LR−1 (zone 0 DR=.78; zone 1 DR=8.2). This zeroes phenotype-level additive epistasis at START, not necessarily whole genetic mating-epistasis mechanisms
5. dispersal_half: m=.36 (narrower geometric pollen kernel; not absence of pollen flow)
6. dispersal_double: m=1.44 (broader geometric pollen kernel; not 144% migration)
7. assortative_25pct: a=.25 (same-phenotype mating conditional model; not directly measured bee behavior)

One parameter dimension changed per source contrast, selected BEFORE computing. No postreadout substitution, tuning or conditional branch picking. No source C++ edits: reuse v43 ANTO-authored driver with deterministic, strict temporary source-text modifications in v44 generated copy; original eight C++ source translation units unchanged. Same minimal old-nr3 alias/cerr shim and 299 all-whitespace family-padding line removal only in temporary input copy as documented v43. Source checkout and all v42/v43 frozen receipts replay-checked.

## Scoring and numerical gates
For all 7 scenarios × 5 transects, record original author C++ source calcTotLike(), original MakeLikelihood() summed separately over (A) adult 4-colour counts and (B) conditional maternal colour × offspring 4-colour counts. Record Delta logL vs baseline for each component and total; number of 5 transects with positive/negative Delta, model-predicted aggregate adult/offspring 4-colour probabilities, iterations, final getMaxChg(). Absolute tolerance 1e−6 on likelihood decomposition/probability normalization, 1e−7 on baseline per-transect v43 replay. Verify actual getMaxChg <= original 1e−8+1e−14: source original contSim also has a 5,000-step likelihood-abort and max-generation gates, so generations less than 100,000 alone do NOT prove equilibrium.

Each scenario runs twice, deterministic outputs must match exactly. Source original zero-probability nonzero-count floor (-1e11) is NOT interpreted as biological likelihood. Scenarios with source numerical/technical exit or nonconvergence remain in output explicitly as BLOCKED_NUMERICAL_OR_ORIGINAL_CODE with the exact error category; no replacement or silent omission. Baseline failure means entire analysis blocked. Report source model conditional one-factor sensitivity, NOT empirical selection significance, independent replication, AIC, natural population MCMC fit, measured migration rate, or ancestral colour reversion. Any increase in logL at this arbitrary original START is merely a better local hypothetical parameter value, NOT proof that a mechanism was absent or the paper wrong.

## Biological interpretation
Selection, pollen spread and mating assortment may have different effects on adult regional colour clines vs mother-conditioned offspring fractions. Because these are co-determined in source model, the scenarios are *conditional perturbations* and cannot isolate real-world causal selection. Experimental source effect sizes must not be inferred from differences in unoptimized START likelihoods. Existing v26–v43 biological claims and source outcomes must remain unchanged.

## Deliverables
Source-frozen design and seven vectors, user-authored v44 source-run code, adversarial tests, original author source CI with 7-scenario readout, immutable derived JSON and review-ready summary. Stack Draft PR on v43 PR #53; never merge without source and biological review.
