# ANTO v48 — original 2024 Phlox two-training-source-population matched refits

**Scope frozen on 2026-11-19, before v48 model readout.** Motivation was chosen **after** v46 and v47 showed negative deamii transfer. Therefore this is a targeted **POSTHOC** challenge of the small-training-population explanation, not a new prospective proof. Preserve v46 original 20-source-population primary *NOT_ESTABLISHED* decision and v47 no-refit exploratory reference verbatim.

## Exact original biological input and unit

Goulet-Scott et al. (2024) Nature Communications 10.1038/s41467-024-49075-6 original author **PhloxHopkins/PhloxFieldAdaptation** Git commit d792fe58cce36721394adfd97e9f412396578e86, original Dataset_fitness_leaftraits_distance_lat_long.csv, Git blob SHA1 4c780b0abe644ac69dae86a0450eb598c160626b, 175377 bytes. v46 verified the source: 963 clonal ramet observations, 284 genotypes, original 20 populations, three common gardens, original original taxa amoena (8 pops), pilosa (9), deamii (3); same genotype appears across gardens but belongs to one source population. New original independent source studies in v48 = ZERO; this is a re-fit of the **same 2024 dataset**. No flower colour, pollinator, genotype-paternity, multi-generation fitness, other garden, other radiation or independently sampled taxon observations.

## Controlled mechanism of the diagnostic

Earlier v47 matched the number of **evaluated** original populations (3 each), but did not match the number of **training** source populations within taxon: the original v46 deamii holdout learned from 2 remaining deamii populations whereas the amoena and pilosa holdouts had 7 and 8 same-taxon training populations. The v48 counterfactual brings each comparison to exactly **two training source populations and one withheld population of the SAME taxon**.

For each taxon (amoena n=8, pilosa n=9, deamii n=3 original populations):
- Exhaust all unordered original source-population triples: C(8,3)=56, C(9,3)=84, C(3,3)=1 -> **141 overlapping 3-pop panels**.
- For each of the three populations in a panel, withhold it entirely as test, fit only on the other TWO same-taxon original populations, with every original genotype-clone confined to its own train/test source population, pooled across the three original common gardens. **423 complete population-heldout fits** across all panels. Other original taxa are **excluded from the training and testing of that panel**.
- Each new model fitted per panel-fold has a TRAINING-ONLY intercept for each of the THREE original gardens. The baseline B0 therefore always compares two same-taxon training-population means for each original garden. There is no target population mean or target outcome in training, normalization or model selection.
- Center response and leaf-PC1/PC2 using each original garden's training-only mean, with training-only feature scaling by the pooled **population variance SD** across the training ramet rows (not across source populations); reject any zero SD.
- Fixed ridge penalty lambda=10, only for PC slopes, independent of endpoint and taxon. No post-outcome penalty change.
- Exact two candidate model structures: **G** adds two shared leaf PC1/2 coefficients pooled over original gardens, and **E** adds two coefficients specific to each of the THREE gardens (six slopes). When only one taxon is admitted to a fit, v46 T collapses to G and v46 TE collapses to E; do not invent distinct fourth/fifth effective models.
- Primary log(1+fruits_2019) predictive MSE, secondary log(1+flowers_2019). Original preregistered log(1+biomass_2019) remains technically nonadmissible in v46 because 13 original values are negative (minimum -233); do not clip/drop/shift biomass.
- Aggregate three one-pop test MSE losses within each triple **with equal original source-population weights** to form the panel-level fractional model gain 1 - sum(test MSE model)/sum(test MSE B0). All source-population ramet counts remain as observed and variable; matching source-population number does NOT guarantee matching genotype/ramet sample sizes or trait diversity.

## Frozen outputs and decisions

All 141 panels, their 423 original population-heldout fold MSEs for B0/G/E and both outcomes, exact panel gains, 56/84/1 panels by taxon, original 20 source-population identities, exact source hash, and per-panel training counts must be retained in the original-source-derived JSON.

Compare the ONE deamii panel gain with all 56 amoena and 84 pilosa panels (primary 140 **overlapping** reference panels), reporting counts at-or-below deamii and each taxon's min/median/max. Also report proportion of 140 reference panels with negative gain, and deamii individual 3 fold changes. Separate results by G and E and fruit versus flower, **no favorable subset selection**.

**Interpretation precommitted before v48 readout:** if G and E deamii remain worse than all/most same-training-population controls, the *source-population-count mismatch* alone is not a sufficient explanation for the observed prediction problem, conditional on this model family and original data. If many balanced comparison panels are as bad or worse, earlier v47's 0/140 contrast was sensitive to unequal training-population count and must be downgraded. Neither outcome provides a valid p-value or causally identifies taxon, pollinator, genetics, leaf trait mechanism or natural trait regain. The 140 panels overlap, 2024 source gardens all same three, and model refits use original samples. No new universal transfer confirmation, no pooling with PhloxSelection log-likelihood, ANTO pigment studies, SDEA or izu-core.

## Reproducibility

Commit protocol and source registry before new model readout; standard-library Python source parser reuses v46's source-identity and clone checks with a new two-train-pop fold solver; freeze derived-only original-source JSON and methods, adversarial tests for each 141/423 expected source/pop/clone/fit/garden/denominator boundary; run original source hash-verified recomputation under dedicated read-only GitHub Actions and general CI. Remain Draft on branch stacked on v47 PR #57, not on main without review.
