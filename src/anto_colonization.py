"""Synthetic plant-lineage colonization and floral re-expression simulator.

A discrete-generation, finite-population *toy* model. Plants produce seeds, seeds
move between regions, pollen can move without seeds, and recruit establishment
is limited by region capacity. Re-expression, immigration and introgression
are logged separately. This is NOT a fitted model or a population-genetic
model of real Cirsium; 'W-' refers to strongly reduced intrinsic competence,
not irreversible loss.
"""

from collections import Counter
from dataclasses import dataclass, field
from typing import Mapping
import math
import numpy as np

STATES = ("P", "W+", "W-")
ROUTES = ("intrinsic_reexpression", "introgressive_reacquisition",
          "pigmented_immigrant", "pigmented_inheritance")


@dataclass(frozen=True)
class Plant:
    lineage: str
    state: str
    introgressed: bool = False

    def __post_init__(self):
        if self.state not in STATES or not self.lineage:
            raise ValueError("Each plant needs a valid state and lineage ID")


@dataclass
class Parameters:
    regions: tuple[str, ...]
    capacity: Mapping[str, int]
    seed_migration: Mapping[tuple[str, str], float] = field(default_factory=dict)
    pollen_delivery: Mapping[tuple[str, str], float] = field(default_factory=dict)
    fecundity: float = 2.5
    loss: float = 0.0           # P -> W+ per offspring
    reexpression: float = 0.0   # W+ -> P, intrinsic, per offspring
    degradation: Mapping[str, float] = field(default_factory=dict) # W+ -> W-
    repair: float = 0.0         # W- -> W+
    introgression: float = 0.0  # white maternal offspring -> P when P pollen accessible
    catastrophe: Mapping[str, float] = field(default_factory=dict)
    recruitment_weight: Mapping[tuple[str, str], float] = field(default_factory=dict)

    def __post_init__(self):
        if not self.regions or len(set(self.regions)) != len(self.regions):
            raise ValueError("Unique nonempty regions required")
        if set(self.capacity) != set(self.regions):
            raise ValueError("Specify each region capacity")
        if any(type(v) is not int or v <= 0 for v in self.capacity.values()):
            raise ValueError("Capacities must be positive integers")
        for p in (self.loss, self.reexpression, self.repair, self.introgression):
            _probability(p)
        if self.reexpression + max(self.degradation.values(), default=0) > 1 + 1e-12:
            raise ValueError("W+ reexpression and degradation compete")
        if not math.isfinite(self.fecundity) or self.fecundity < 0:
            raise ValueError("Fecundity must be finite and nonnegative")
        for region in list(self.degradation) + list(self.catastrophe):
            if region not in self.regions:
                raise ValueError("Unknown region in a rate table")
        for v in list(self.degradation.values()) + list(self.catastrophe.values()):
            _probability(v)
        migration_out = Counter()
        for network in (self.seed_migration, self.pollen_delivery):
            for (src, dest), probability in network.items():
                if src not in self.regions or dest not in self.regions or src == dest:
                    raise ValueError("Links must connect distinct valid regions")
                _probability(probability)
                if network is self.seed_migration:
                    migration_out[src] += probability
        if any(v > 1 + 1e-12 for v in migration_out.values()):
            raise ValueError("Seed migration probabilities must sum to <= 1")
        for (region, state), w in self.recruitment_weight.items():
            if region not in self.regions or state not in STATES or w <= 0 or not math.isfinite(w):
                raise ValueError("Recruitment weights must be positive and valid")


def _probability(value):
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("Probability must lie in [0,1]")


def visible(population):
    """Observer sees P or white: hidden competence must NOT leak."""
    return {r: {"pigmented": sum(p.state == "P" for p in plants),
                "white": sum(p.state != "P" for p in plants)}
            for r, plants in population.items()}


def latent_truth(population):
    """Simulator oracle only; not an admissible colour-only observation."""
    return {r: {"competent_white": sum(p.state == "W+" for p in plants),
                "degraded_white": sum(p.state == "W-" for p in plants),
                "introgressed": sum(p.introgressed for p in plants)}
            for r, plants in population.items()}


def _pollen_access(region, population, params):
    """Probability of encountering coloured pollen; donors may be elsewhere.

    Delivery is an exposure probability conditional on >=1 P donor (not
    a measured pollen or Mendelian transmission model).
    """
    sources = [src for (src, dst), p in params.pollen_delivery.items()
               if dst == region and any(a.state == "P" for a in population[src])]
    nonarrival = 1.0
    for src in sources:
        nonarrival *= 1.0 - params.pollen_delivery[(src, region)]
    return 1.0 - nonarrival


def _destination(src, params, rng):
    links = [(dst, p) for (a, dst), p in params.seed_migration.items() if a == src]
    u = rng.random()
    for dst, p in links:
        if u < p:
            return dst
        u -= p
    return src


def _offspring(mother, source, destination, pollen_probability, params, rng):
    """Return (child, immediate causal route for P recruitment or None)."""
    if mother.state != "P" and rng.random() < params.introgression * pollen_probability:
        # Idealized introgression transfers the colour-restoring genetic effect;
        # source lineage identity of the maternal ancestor remains tracked.
        return Plant(mother.lineage, "P", True), "introgressive_reacquisition"
    if mother.state == "P":
        if rng.random() < params.loss:
            return Plant(mother.lineage, "W+", mother.introgressed), None
        route = "pigmented_immigrant" if source != destination else "pigmented_inheritance"
        return mother, route
    if mother.state == "W+":
        u = rng.random()
        if u < params.reexpression:
            return Plant(mother.lineage, "P", mother.introgressed), "intrinsic_reexpression"
        if u < params.reexpression + params.degradation.get(destination, 0.0):
            return Plant(mother.lineage, "W-", mother.introgressed), None
        return mother, None
    if mother.state == "W-" and rng.random() < params.repair:
        return Plant(mother.lineage, "W+", mother.introgressed), None
    return mother, None


def simulate(params, initial, generations, seed=1729):
    """Return histories and accepted-recruitment route counts, not attempted births.

    Nonoverlapping generations; reproductive adults may be lost to catastrophe
    before making offspring; imported seed can recolonize empty regions.
    Source provenance is maternally tracked, with hybrid flag on descendants.
    """
    if generations < 0 or not isinstance(generations, int):
        raise ValueError("Generations must be a nonnegative integer")
    if set(initial) != set(params.regions):
        raise ValueError("Initial populations must specify every region")
    pop = {r: list(initial[r]) for r in params.regions}
    if any(len(pop[r]) > params.capacity[r] for r in params.regions):
        raise ValueError("Initial carrying capacity exceeded")
    rng = np.random.default_rng(seed)
    histories = [{"generation": 0, "regions": visible(pop), "oracle_latent": latent_truth(pop)}]
    events = []
    for t in range(1, generations + 1):
        survivors = {r: ([] if rng.random() < params.catastrophe.get(r, 0.0) else pop[r])
                     for r in params.regions}
        pollen = {r: _pollen_access(r, survivors, params) for r in params.regions}
        candidates = {r: [] for r in params.regions}
        for src in params.regions:
            for mother in survivors[src]:
                for _ in range(rng.poisson(params.fecundity)):
                    dst = _destination(src, params, rng)
                    child, route = _offspring(mother, src, dst, pollen[dst], params, rng)
                    candidates[dst].append((child, route, src != dst))
        next_pop = {}
        counts = {r: Counter() for r in params.regions}
        for dst in params.regions:
            pool = candidates[dst]
            k = min(params.capacity[dst], len(pool))
            if k == 0:
                recruits = []
            elif k == len(pool):
                recruits = pool
            else:
                # Weighted lottery sampling WITHOUT replacement, with
                # independent exponential keys. Large weight -> higher odds.
                rank = []
                for j, (child, _, _) in enumerate(pool):
                    w = params.recruitment_weight.get((dst, child.state), 1.0)
                    rank.append((-math.log(max(rng.random(), 1e-300)) / w, j))
                chosen = sorted(rank)[:k]
                recruits = [pool[j] for _, j in chosen]
            next_pop[dst] = [x[0] for x in recruits]
            for child, route, migrant in recruits:
                if route is not None:
                    counts[dst][route] += 1
                if migrant:
                    counts[dst]["seed_immigrant"] += 1
                    if not survivors[dst]:
                        counts[dst]["seed_recolonist"] += 1
        pop = next_pop
        histories.append({"generation": t, "regions": visible(pop), "oracle_latent": latent_truth(pop)})
        events.append({"generation": t, "regions": {r: dict(counts[r]) for r in params.regions}})
    return {"histories": histories, "events": events,
            "final": {r: [{"lineage": p.lineage, "state": p.state,
                           "introgressed": p.introgressed} for p in pop[r]]
                      for r in params.regions}}


def first_pigmented(result, region):
    """First generation with any P in region, or None."""
    for h in result["histories"]:
        if h["regions"][region]["pigmented"]:
            return h["generation"]
    return None


def order_experiment_probability(delta_source, delta_island, delay, horizon,
                                 r_island, n_founders=1):
    """Closed-form counterfactual for dispersal timing of white competent seeds.

    Early colonists spend delay generations on island; late colonists spend
    delay generations at source, then both face horizon generations on island.
    Re-expression is disabled during delay; no extinction, migration after
    arrival, introgression or repair. Outputs probability >=1 of reacquisition.
    """
    for v in (delta_source, delta_island, r_island):
        _probability(v)
    if r_island + delta_island > 1:
        raise ValueError("Competing transitions must sum to <= 1")
    if any(type(x) is not int or x < 0 for x in (delay, horizon)):
        raise ValueError("Delay/horizon must be nonnegative integers")
    if type(n_founders) is not int or n_founders <= 0:
        raise ValueError("Positive integer founder count required")
    p_gain = 0.0
    if r_island + delta_island > 0:
        p_gain = r_island / (r_island + delta_island) * (1 - (1-r_island-delta_island)**horizon)
    def at_least_one(p):
        return 1 - (1-p)**n_founders
    return {"early": at_least_one((1-delta_island)**delay * p_gain),
            "late": at_least_one((1-delta_source)**delay * p_gain),
            "p_gain_if_competent": p_gain}


def first_pigmented_route(result, region):
    """Competing first-arrival routes; a same-generation tie is reported.

    Observing coloured flowers alone does not supply these route labels. They
    are simulator-oracle event logs, and cannot be used as empirical labels
    without ancestry/molecular evidence.
    """
    time = first_pigmented(result, region)
    if time is None:
        return {"generation": None, "route": "not_observed"}
    if time == 0:
        return {"generation": 0, "route": "initially_present"}
    observed = result["events"][time-1]["regions"][region]
    causes = sorted(k for k in ("intrinsic_reexpression",
                                "introgressive_reacquisition",
                                "pigmented_immigrant") if observed.get(k, 0) > 0)
    label = causes[0] if len(causes) == 1 else ("simultaneous_multiple" if causes else "undetermined")
    return {"generation": time, "route": label, "concurrent_causes": causes}
