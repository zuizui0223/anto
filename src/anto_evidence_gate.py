"""Fail-closed, tier-specific empirical data admission for anto.

The presence of raw leaf-RNA data or a dated lineage split does NOT imply
floral re-expression competence or historical gain polarity. This tool
classifies source-study level evidence and makes missing measurements explicit.
"""
from pathlib import Path
import json

REQUIRED_BOOLEAN = (
    'flower_colour_recorded','geographic_occurrence_recorded',
    'dated_lineage_tree_available','colour_transition_direction_established',
    'colonization_direction_established','young_leaf_rna_registered',
    'same_individual_genomic_flower_link','floral_pigment_chemistry',
    'flower_stage_rna','functional_reexpression_test',
    'pigment_locus_gene_flow_resolved',
    'replicated_independent_reexpression_events',
)
ALLOWED_RESOLUTION = ('taxon','individual_with_geographic_confounding','experimental','population','individual')
MECHANISTIC = ('colour_transition_direction_established', 'floral_pigment_chemistry',
              'flower_stage_rna', 'functional_reexpression_test')
GEOGRAPHIC = ('geographic_occurrence_recorded', 'dated_lineage_tree_available',
              'colonization_direction_established')
FULL = MECHANISTIC + GEOGRAPHIC + ('same_individual_genomic_flower_link',
      'pigment_locus_gene_flow_resolved','replicated_independent_reexpression_events')


def load_registry(path):
    document=json.loads(Path(path).read_text(encoding='utf-8'))
    if document.get('registry_version')!='anto_empirical_evidence_v3':
        raise ValueError('Unexpected registry version')
    records=document.get('source_records')
    if not isinstance(records,list) or not records:
        raise ValueError('Expected source records')
    ids=[]
    for r in records:
        if not isinstance(r,dict) or not isinstance(r.get('id'),str) or not r['id']:
            raise ValueError('Every record needs a study identifier')
        ids.append(r['id'])
        for flag in REQUIRED_BOOLEAN:
            if type(r.get(flag)) is not bool:
                raise ValueError(f"Missing or non-boolean evidence {flag} for {r['id']}")
        if r.get('colour_resolution') not in ALLOWED_RESOLUTION:
            raise ValueError('Unrecognized colour resolution')
        if not isinstance(r.get('primary_url'),str) or not r['primary_url'].startswith('https://'):
            raise ValueError('A public source URL is required')
        if not isinstance(r.get('limitations'),str) or len(r['limitations'])<20:
            raise ValueError('Limitation provenance required')
        if not isinstance(r.get('evidence_source_urls'),list) or any(
                not isinstance(u,str) or not u.startswith('https://')
                for u in r['evidence_source_urls']):
            raise ValueError('Evidence URLs must be HTTPS')
        if r['functional_reexpression_test'] and not (
                r['floral_pigment_chemistry'] and r['flower_stage_rna']):
            raise ValueError('Mechanistic positive control requires matched biochemical and floral evidence')
    if len(set(ids))!=len(ids):
        raise ValueError('Duplicate study IDs')
    return records


def evaluate(r):
    missing_mechanistic=[a for a in MECHANISTIC if not r[a]]
    missing_geographic=[a for a in GEOGRAPHIC if not r[a]]
    missing_full=[a for a in FULL if not r[a]]
    mechanism_positive_control = not missing_mechanistic
    biogeographic_history_available = not missing_geographic
    # A historic lineage date without trait-event polarization is *not* enough.
    ready_for_integrated_accessibility = not missing_full
    if ready_for_integrated_accessibility:
        tier='INTEGRATED_ACCESSIBILITY_ELIGIBLE'
    elif mechanism_positive_control:
        tier='MOLECULAR_MECHANISM_ONLY'
    elif r['dated_lineage_tree_available'] and r['flower_colour_recorded']:
        tier='VISIBLE_LINEAGE_HISTORY_ONLY'
    elif r['flower_colour_recorded']:
        tier='VISIBLE_COLOUR_ONLY'
    else:
        tier='UNSUITABLE'
    return {
        'system': r['system'], 'id':r['id'],
        'admission_tier':tier,
        'mechanistic_regain_positive_control': mechanism_positive_control,
        'geographic_history_identified': biogeographic_history_available,
        'integrated_accessibility_eligible':ready_for_integrated_accessibility,
        'missing_mechanism_evidence':missing_mechanistic,
        'missing_geographic_evidence':missing_geographic,
        'missing_integrated_evidence':missing_full,
        'strongest_allowed_claim':(
            'Joint trajectory analysis, conditional on current validation gates' if ready_for_integrated_accessibility else
            'Molecular regain within this studied system only' if mechanism_positive_control else
            'Lineage split, taxon colour and geography; no inferred historical regain' if r['dated_lineage_tree_available'] else
            'Visible phenotype variation only'
        ),
    }


def report(path):
    records=load_registry(path)
    return {'evidence_gate_version':'v3',
            'integrated_accessibility_systems':sum(evaluate(r)['integrated_accessibility_eligible'] for r in records),
            'systems':[evaluate(r) for r in records]}


if __name__=='__main__':
    import sys
    source=(sys.argv[1] if len(sys.argv)>1 else Path(__file__).resolve().parents[1] / 'data' / 'empirical_sources_v3.json')
    print(json.dumps(report(source),indent=2,ensure_ascii=False))
