# ANTO v59 — Season transfer fails while water-treatment transfer remains positive in one Mertensia experiment

2026-11-26. **Retrospective same-source diagnostic, no additional independent original study.** Original v58 frozen result and code are unchanged. v59 protocol and eleven specified model/split cases were committed before v59 source outcome calculations. The original v58 plant-level cohort remains 105 records. The four original phenology weeks are not four independent published studies.

## Source authentication and design correction

Gallagher & Campbell (2021), DOI 10.3389/fevo.2021.641693, original Dryad 10.7280/D16D7Z v4 and Zenodo 4731053. The exact original workbook SHA256 152e09ab8bd0f49c2c5052e284fbaeb3f7362fd0f82998ed1a27001bd0e7250d and the original metadata DOCX SHA256 b41dd78899e1b451731917e923e739bbf1d0cbbd734867b1faa7d04713a8d968 were checked against both authoritative archives. All 80 v58 group/model MSE values reproduced exactly (maximum difference **0**). The source workbook has 115 distinct raw plantID tokens, and all 105 eligible v58 plantIDs are distinct; there are no repeated IDs across original weeks, arrays or water codes. This confirms absence of identical source plant-ID duplication in the scored cohort. It does not establish clone/genetic independence.

The original article explicitly reports that week 4's twelve plants ceased flowering mid-week and the remaining plants were re-randomized into three new arrays. Nominal array labels may therefore share exposure context; no duplicated response leakage has been demonstrated. The conservative new split retains the original 15 week1–3 arrays but holds out ALL week4 plants together (16 folds). A stronger temporal test holds out each whole week (four folds). Whole-water-code holdouts are three different folds with no heldout-treatment dummy coefficient; plants from the same weeks/arrays may remain in training.

Original dictionary verifies mVWC is mean volumetric water content and mVisRate is mean plant visitation rate; the exact mseed field definition was not located in dictionary snippets. Retain v58's log1p(source mseed) label, not an inferred outcome unit. All 105 original measured mVWC values were admitted, with original source-code A/C/D means 17.4333/13.0180/10.4080 respectively; do not infer treatment label names without source dictionary confirmation.

## Fixed context and split comparison

B = original flower count + corolla width/length + source treatment dummies (except unseen water-code splits), with optional numeric linear week and measured VWC; V adds original log1p visit rate. Same cohort, alpha=5, train-only standardization. Negative gain = higher prediction MSE after adding visits. Common weighting shown below is 105-plant-weighted MSE, not equal-fold averages (fold grains differ).

| Context | Holdout | B MSE | V MSE | B−V | positive folds | positive week means |
|---|---|---:|---:|---:|---:|---:|
| no week | nominal 20 arrays | .172495448 | .158718617 | +.013776831 | 13/20 | 2/4 |
| linear week | nominal 20 arrays | .169618611 | .163050905 | +.006567706 | 14/20 | 2/4 |
| linear week + VWC | nominal 20 arrays | .161955786 | .153702337 | +.008253448 | 13/20 | 2/4 |
| no week | conservative 16 groups | .180593538 | .162565283 | +.018028255 | 10/16 | 2/4 |
| linear week | conservative 16 groups | .174697303 | .177474238 | -.002776935 | 9/16 | 1/4 |
| linear week + VWC | conservative 16 groups | .175955271 | .176765313 | -.000810042 | 8/16 | 1/4 |
| no week | entire four weeks | .176663618 | .187634419 | -.010970801 | 1/4 | 1/4 |
| **linear week (primary)** | **entire four weeks** | **.185193958** | **.282788069** | **-.097594111** | **0/4** | **0/4** |
| linear week + VWC | entire four weeks | .181348183 | .259180257 | -.077832073 | 0/4 | 0/4 |
| linear week, treatment blind | entire 3 water codes | .171122583 | .157245052 | +.013877532 | 3/3 | 3/4 |
| linear week + VWC, treatment blind | entire 3 water codes | .153862195 | .145358047 | +.008504148 | 3/3 | 2/4 |

All three full-week contexts have negative overall gain; this is not merely a numerical consequence of one chosen linear week model. The primary whole-week equal-week gain is -0.092923968 (0/4 positive weeks), 105-plant-weighted gain -0.097594111. In contrast, secondary water-code-heldout predictions improve in each of three codes, with or without measured VWC. The training weeks overlap the evaluation weeks for water-code splits, so this is a **difference in predictive transfer domains**, not a controlled causal effect of switching climate axes.

The conservative week4 regrouping makes the linear-week visitor increment negative at common plant weights, with and without VWC. It remains positive if week is omitted. Do not claim every conservative model fails; model-context dependence is retained.

## Why original and new week effects differ

Primary whole-week folds, all four retained:

| withheld week | n | B MSE | V MSE | 2 mean(e*d) | mean(d^2) | B−V | visitor-range outside count |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 28 | .218442472 | .531676362 | -.093455236 | .219778654 | -.313233890 | 16 |
| 2 | 27 | .115016271 | .146376910 | -.018333169 | .013027471 | -.031360640 | 0 |
| 3 | 27 | .102742598 | .104433749 | +.000554316 | .002245467 | -.001691151 | 0 |
| 4 | 23 | .323890734 | .349300928 | -.016654892 | .008755301 | -.025410194 | 0 |

Week 1 has extreme visitor support mismatch: 16/28 observed log visitor rates beyond training range, plus all 28 target numeric week indices out of range. Even though training slope is positive (0.908714 per log visit-rate unit), the visitor-model correction mean is +.372558 when mean B residual is -.275403, indicating aggregate overprediction worsens. Week 2 and week 4 also have negative error-correction alignment. Week3 is different: alignment is weak positive but smaller than the correction penalty. The same pattern in original v58 replay (week3 alignment +.000744212, squared correction .004784115, plant-weighted gain -.004039903) identifies a weak informative correction rather than evidence of a biological visit-effect reversal. This is a prediction-error algebra identity, not a decomposition of ecological mediation.

Measured VWC attenuates the full-week failure but does not remove it: 0/4 weeks still positive. No-week full-week sensitivity also negative overall, with only week4 positive. These observations motivate an actual guild/season-context hypothesis, not an inferred pollen mechanism.

## Fixed observational comparison and decision

200 source visit-label permutations within original week × water code, seed 20261126, preserve week distributions and outcomes; original B predictions unchanged. Original equal-week gain -.092923968; shuffled median -.008133113, q2.5 -.108995420, q97.5 +.027788707, 191/200 >= observed. Both primary direction and placebo separation gates fail. **WHOLE_WEEK_VISITOR_TRANSFER_NOT_ESTABLISHED**. This is a descriptive reference, not a causal p-value or a significance test for a universal ecological null.

## Integration and limitations

This is a counterexample to promoting the v58 within-season array result or the three-study v52/v56/v58 directional agreement into a universal environmental-transfer law. Original v58 was a valid comparison under its stated nominal split; its test domain does not withhold phenology. Source v59 identifies a difference between within-season spatial transfer, water-treatment transfer, and entirely unseen seasonal contexts. It does not prove that water has no effect, that a pollinator community shift is causal, or that visitation is biologically unimportant. Future tests should distinguish measured Bombus/non-Bombus composition, predictor support and pollen transfer, with a new independent study before broad generalization. All raw original plant IDs/measurement rows are excluded from the repository.

Original source CI 45758964051 succeeded with 10 unit tests and complete new first readout; exact-head saved-result replay must be separately checked. Draft stacked on #68, main unchanged.
