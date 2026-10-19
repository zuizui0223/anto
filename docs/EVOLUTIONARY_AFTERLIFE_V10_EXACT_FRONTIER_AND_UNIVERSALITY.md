# ANTO v10 — The evolutionary afterlife of a lost visible trait

**Date:** 2026-10-19 JST. **Status:** exact finite synthetic model and cross-repository *research design*. No new original plant phenotype, allele, Cirsium ancestry, anthocyanin, pollinator, or viable-achene data have been collected. No proposal to interpret the same cross-repository source model as independent replication.

## Scientific question

**When visible expression has disappeared, which distinctions determine whether a trait can reappear and be expressed repeatedly before the generative resources are lost?**

The paper-scale hypothesis is *not* that reversals are possible. Trait reversals and retained developmental pathways are well established in previous literature: Forni et al. (2026), *Biological Reviews* DOI 10.1002/brv.70168, distinguish pleiotropy, plasticity and reticulation; Esfeld et al. (2018), *Current Biology* DOI 10.1016/j.cub.2018.10.019 experimentally described Petunia AN2 functional resurrection; Berardi et al. (2021), *Plant Cell* DOI 10.1093/plcell/koab114 documented a distinct Petunia pigmentation route. Petunia constitutes **one plant genus**, not two independent clades.

The proposed contribution must instead **measure the quantitative difference between hidden functional capacity, expression, repeated expression, and effective ecological establishment**, and externally identify when an ancestral resource versus new genetic input explains a return. Here the last outcome—ecological establishment—is explicitly **not** modeled or validated.

## Four non-equivalent biological events

1. **Visible absence**: no expressed trait at the observed time.
2. **Resource availability**: all jointly necessary functional gene/development modules are carried by the *same lineage* (not simply detectable homologs or a collection of marginal allele frequencies).
3. **First expression**: a visibly absent lineage makes an expressed phenotype at least once.
4. **Three consecutive positive snapshots**: the lineage is expressed for at least three consecutive model generations at some time before horizon H. This is a narrow persistence *proxy*, **not** actual evolutionary establishment, fixation, function, fitness, viable seed production, or indefinite trait retention.

The state also records **whole-toolkit loss** and whether the **first three-generation run** followed initial retained capacity or **post-loss full-toolkit acquisition/reassembly**. Immigrants already displaying the trait are **not** part of this model, because geographic recolonisation is not intrinsic regain.

## Exact finite model and novel-to-this-program cross-object interface

All founders start visibly negative. Initially probability theta has a complete joint functional toolkit; 1-theta lacks the complete toolkit. At generation t, declared probabilities are:

- a: latent joint-complete toolkit activates expression;
- r: an expressed lineage remains expressed for another generation;
- d: a lineage loses its joint-complete toolkit, replacing either latent or expressed state with the no-toolkit state;
- g: a lineage currently without the whole toolkit acquires/reassembles a new complete functional inventory in one step. **g is an abstract acquisition operation**, not a measured recombination, introgression, mutation, dispersal or genotype-reconstruction rate.

Require 0<=a,r,d,g<=1, a+d<=1 and r+d<=1. Both first-expression and three-in-a-row are first-passage events tracked with **exact Fraction probabilities**. Two provenance labels native/acquired remain distinct. A reached streak does **not** freeze the current biological state: the lineage can later return to white or lose capacity. The event indicator remains true, while *present expression* can decline.

Run: python -m unittest discover -s tests -p 'test_afterlife_frontier_v10.py' -v

Replay: python scripts/run_afterlife_frontier_v10.py --output results/afterlife_frontier_v10.json

A source-pinned EGC simulator is **not** imported; this is an independently specified generic conditional finite chain. Existing EGC-v10 comparison remains an external mathematical consistency check, not an independent biological dataset.

### Conditional closed-world recovery frontier

For a homogeneous rate regime with g=0, first-step equations at the latent, first-expressed and twice-expressed states give, for one initially complete toolkit:

P(ever 3 consecutive positives before irreparable loss)
  = a r^2 / { a r^2 + d [1+a(1+r)] }.

A population mixture with initial fraction theta has the above probability multiplied by theta.

Proof sketch: let x denote the eventual three-streak probability from latent, y from 1-expression, z from 2-expressions. Then
x=(1-a-d)x+a y;
y=r z+(1-r-d)x;
z=r+(1-r-d)x.
Solving gives the exact expression. Degenerate cases a=0 or r=0 are assigned zero. For d=0 and positive a,r, eventual achievement is certain *conditional on initial capacity*, with no implication about true fitness.

At fixed positive degradation d and fixed theta, the eventual event probability is ordered by

Q(a,r) = a r^2 / [1+a(1+r)],

**not by activation rate a alone**. This is an elementary model-conditional competing-risk comparison, not a newly established universal mathematical theorem. It directly motivates separating expression and sustained expression in future biology.

If g>0 and a,r>0 remain strictly positive forever, eventual three-streak probability tends to one even when theta=0. That is a consequence of a permanently available acquisition source and infinitely many repeated opportunities—not ecological evidence that loss is always reversible. At finite H, this model keeps g-dependent delay, original-resource versus reacquired-resource histories and low success probabilities.

### Sharp uncertainty from first/pairwise module measurements

ANT0 v7 already proves that identical single-module probability 1/2 and pairwise joint 1/4 are consistent with theta=P(all three retained) anywhere in [0,1/4].

Because v10 groups all initial non-complete inventories into a single loss class and only theta changes, its finite-horizon event probability has the exact affine form

P(event by H) = theta A_H + (1-theta) B_H,

where A_H=P(event | initially complete) and B_H=P(event | initially incomplete, possibly acquiring later). Thus its **sharp interval** given the first/pairwise moments is the min/max over theta endpoints, for this declared model only. This is mathematical *partial identification*, **not** an empirical CI.

With no acquisition, B_H=0. With g>0, B_H>0 may hold: **observing a coloured return does not identify whether all resources were retained throughout the white phase**. Actual partial inventories can have different acquisition/recombination behaviour; then the affine-bound result need not apply.

## New numerical contrasts (synthetic, fixed illustrative parameters)

All four homogeneous scenarios: theta=1/8; horizon 12 generations. Fast/slow comparison deliberately fixes degradation d=2/25 and changes activation and repeat rates:

| Synthetic regime | a | r | d | g | First visible by gen 1 | Ever three consecutive by gen 12 | Asymptotic event if constant regime |
|---|---:|---:|---:|---:|---:|---:|---:|
| Fast flash, weak repeat | 0.60 | 0.15 | 0.08 | 0 | **7.50%** | **0.708%** | 1.135% |
| Slower activation, strong repeat | 0.20 | 0.75 | 0.08 | 0 | **2.50%** | **5.698%** | 6.378% |
| Same slow/strong plus acquisition | 0.20 | 0.75 | 0.08 | 0.05 | **2.50%** | **17.001%** | 100% only under inexhaustible acquisition / unlimited generations |
| No activation, no acquisition | 0 | 0.75 | 0.08 | 0 | 0 | 0 | 0 |

**Concrete rank reversal:** the fast-flash regime has **three times** the probability of first expression at generation 1, yet about **one-eighth** the 12-generation probability of ever observing a three-generation streak. This result is a **constructed conditional counterexample**, not an estimate of how often reversals occur in angiosperms.

In the acquisition arm the 12-generation three-streak event is exactly decomposed into original initial-retained toolkit routes (**5.698%**) and post-loss acquired/reassembled toolkit routes (**11.303%**), with no double counting of the first streak.

A separate 12-generation exposure-order example fixes the same **six generations of favorable** and **six unfavorable** rates, with equal total exposure but reversed order:
- favorable first, then unfavorable: **7.0878%** event;
- unfavorable first, then favorable: **0.5677%** event.

This is an illustrative **noncommuting ecological history** under the chosen first-passage target, not a universal claim that early colonization is advantageous. Time order can matter because the state-transition process is noncommutative; its sign depends on the rates and chosen endpoint.

### Negative and counterfactual controls

- theta=0, g=0 => no return in this closure, even if favorable a/r are specified.
- theta=0, g>0 => delayed re-expression after a capacity-restoration step is possible, but **not** immediate appearance of a pre-pigmented immigrant.
- a=0 or r=0 => no first three-streak event regardless of capacity and acquisition.
- More first expression does not imply more three-streak events.
- Identical visible absence but distinct theta can change future outcome.
- Different environment order with identical exposure totals need not agree.
- Present-day white colour never tells us the actual genetic mechanism; a homology alignment does not measure functional toolkit.

All twelve tests in tests/test_afterlife_frontier_v10.py passed in a local Python environment using a minimal exact functional equivalent of the previously committed v7 feasibility-bound function; the **authoritative repository integration verification is the dedicated GitHub CI workflow**.

## Cross-repository ownership and non-replication rule

| Programme | Genuine science owned there | Proper comparison to v10 |
|---|---|---|
| anto v7/v9/v10 | exact joint-resource bounds, history proxy sufficiency, conditional first-passage frontier | Theoretical result only; no Cirsium pipeline |
| egwe pinned EGC experiments | allele-seeded expression from global trait zero, joint-zero absorption, 4x256 source-pinned finite trajectories | Compare *qualitative* event ranking and absoption classes. Not independent empirical corroboration |
| izu-core Model 3 | reproductive assurance / flower-investment ordered changes and genetic support loss | Specific plant reproductive-fitness/diploid simulation; cannot infer ANTO parameters or claim confirmed flower recolouring |
| aza3 | all Ryukyu Cirsium field, head, same-individual molecular, genomics, introgression and viable-achene inference | **Only owner of Cirsium empirical experiments**; assay all required functional modules together |
| chun | phenotype-colour phylogenetic histories, lineage/trait-context memory | Visible history can guide hypotheses; cannot fill unmeasured floral competence |
| fcp / mois | source-photographic flower hue, spatial turnover and global composition | Identify candidate variation only. Neither original photo is a genotype or evolutionary event |
| Petunia original studies | experimentally documented real pigment-gene resurrection/substitution | one mechanistic positive-control clade, not a set of independent angiosperm replicates |

## Falsifiable cross-clade empirical programme (NOT executed)

The primary question is whether a **joint functional toolkit on the same biological unit** improves prediction of lineage re-expression or future phenotype under controlled assays beyond an observed colour/historical dispersal label.

**Prospective admissibility**: require defined true population/clade, source/chronology of visible loss, functional gene/regulator identity with genotype *and flower-scoped activity*, population ancestry including gene flow and cytotype, and genuinely time-resolved or experimental return outcome. Without these, classify as HISTORY_ONLY, PHENOTYPE_ONLY, or FUNCTIONAL_CAPACITY_ONLY and stop before claiming historical reversal.

**External candidate panel**: Petunia (one positive-control genus), Ryukyu Cirsium (hypothesis-generating; unfinished), Antirrhineae, Mimulus/Erythranthe and Camellia as conditional candidates. At least **two independently admitted additional angiosperm clades beyond Cirsium and Petunia** are required before any general cross-angiosperm conclusion. Do not cherry-pick only documented reversal-positive cases; include competent white non-returners and lost-capacity controls.

**Minimum model comparison on held-out entire clades/populations**:
- M0: present floral phenotype, region, environment, taxon;
- M1: M0 plus declared dated ancestry/colonization history;
- M2: M1 plus same-plant/functionally verified joint toolkit and tissue-scoped expression competence;
- M3: M2 plus migrant/introgressed alleles and path-specific origin support.

Measure held-out log-loss/Brier and calibration for genuinely observed **floral phenotype regeneration/retention events**, without inferring a theoretical waiting-time constant from a fossil-age proxy. Hypothesis: if ancestral history is mainly an incomplete proxy of present joint competence, its incremental forecast gain will shrink after competent modules are actually measured. If it does not, *do not* declare the world non-Markov: inspect missed regulatory state, germline environment, source origin, and assay error.

**Cirsium handoff to aza3**: the active v4 DFR/ANS vs ANR/FLS panel narrows 5 preliminary transcripts to 2 reference-robust family matches, *without establishing their enzyme function or flower expression in Ryukyu samples*. Prioritize same-individual floret anthocyanin/HPLC or LC-MS, validated candidate flower regulatory activation, nuclear ancestry, cytotype and viable achene outcome. A newly coloured lineage does not prove internal reactivation if an already-pigmented migrant or introgressed allele could have arrived.

## Publication decision

- **Now**: publish as a **testable model/gating programme** in anto; add the results only as an explicitly synthetic conditional extension. The loss-persistence/reversal principle and elementary Markov/competing risks are already in the literature and do not warrant a stand-alone Ecology Letters/NEE originality claim merely by formalization.
- **Later, if independent data exist**: paper title *The evolutionary afterlife of lost traits: when does recoverability outlast phenotype?*; main data would need contrasting phenotypically lost clades with functionally measured retained genetic capacity, dated genealogy and real re-expression / reproduction outcome across external systems.
- **If the cross-clade gate fails**: publish the causal Cirsium same-individual/molecular results under aza3 if they independently warrant it; keep v10 as mathematical methods/support instead of inventing an angiosperm-wide law.

### Reproduce

python -m unittest discover -s tests -p 'test_afterlife_frontier_v10.py' -v

python scripts/run_afterlife_frontier_v10.py --output /tmp/anto_afterlife_frontier_v10.json

Both work without any empirical Cirsium file or new third-party Python dependency.
