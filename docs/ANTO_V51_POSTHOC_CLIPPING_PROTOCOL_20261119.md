# ANTO v51 supplementary diagnostic — POST-OUTCOME prediction-clipping audit

2026-11-19. THIS DIAGNOSTIC WAS DESIGNED **AFTER** original v51 first-source results revealed that M0 and M1 have numerically IDENTICAL farthest-two-site loss at all six original source sites. It is retrospective and cannot be labeled preregistered v51 confirmation or selected to rescue a primary result.

Frozen input: original Castellanos et al. 2023 Ulex Zenodo record 7761289 five exact MD5s and the **unmodified v51 source-derived result** stored in ANTO, Git blob SHA1 8f1d5435d2d14c2b2c999671ff9801cffab01a71. Holdouts, models and source eligibility exactly as registered in v51: 202 matched individuals, six sampling SITES within one genetically continuous population, not six independent genetic populations.

For every source sampling site and fixed FULL5, NEAR2, FAR2 original training-site membership:
- Refit identical source-verified M0 and M1 ridge(5), with source train-only scaling and prediction raw-score before clipping.
- For every original heldout individual, record whether its **unclipped** prediction is BELOW 0, ABOVE 1, or inside [0,1], for each model. Aggregate only counts and within-site proportions; do NOT commit raw individual outcomes, site-level pollinator joins, genotype identities or individual prediction values.
- Also check the derived mean squared errors after clipping equal the original frozen v51 FULL5/NEAR2/FAR2 MSEs within 1e−10. If they differ, no interpretation; fail source/score.
- Interpret M0=M1 FAR2 MSE equality only after inspecting boundary clipping; a zero original trait increment is not evidence of a natural biological zero effect.
- Preserve v51 primary absolute-risk finding and FULL5 original SDEA v15 source parity; this diagnostic cannot change originally chosen training-site source sets, ridge, models, response transformation, or eligibility.

A strong saturation pattern indicates the *statistical prediction extrapolation + enforced [0,1] bounds* dominates that out-of-domain comparison. It is not evidence that the real floral–fitness relationship is flat. This remains one original Ulex 2023 source study, one continuous population and six within-region sites.
