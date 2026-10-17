# First executable baseline — synthetic, not empirical

Run `python -m pip install numpy scipy`, then `python -m unittest discover -s tests -v` and `python src/anto_model.py`.

The model follows a *single tagged lineage* moving between regions, with 3 latent expression/competence states in each region. It is **not** a model of whole-region occupancy, colonization success, multi-lineage phylogenetic diversification, pollen flow, or independent phenotype origins. Region-specific degradation and migration need not commute. Tests verify probability conservation, observable-state non-identifiability, and a conditional reacquisition window.

**Scientific controls:** All example rates are invented for simulation; no biological data are loaded. The visible-colour observation maps both W+ and W- to white. Any claim that historical routes can be distinguished requires new data and explicit observation-model tests. The next implementation should add a separate occupancy/branching layer, genetic introgression and a likelihood-based identifiability benchmark rather than calling this toy model a fitted macroevolutionary analysis.

**Sources and boundaries:** aza3 supplies the Ryukyu hypothesis and proposed same-individual molecular tests; izu-core supplies model-dependent evolutionary-order insights; egwe supplies state-sufficiency and movement-operator caution; chun supplies recurrent but non-universal temporal colour memory. None independently establishes the anto central hypothesis.
