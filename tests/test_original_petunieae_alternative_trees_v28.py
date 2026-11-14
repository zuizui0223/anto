"""Synthetic/adversarial tests for actual-source-constrained v28 topology validation."""
from __future__ import annotations

import copy
import hashlib
import math
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import analyze_original_petunieae_alternative_trees_v28 as m
from analyze_original_petunieae_phylogeny_v27 import (
    Node, TreeGateError, parse_newick
)


def binary_nwk(names):
    if len(names)==1:
        return names[0]+":0.1"
    mid=len(names)//2
    return "("+binary_nwk(names[:mid])+","+binary_nwk(names[mid:])+"):0.2"


def fake_source():
    """43 external species/117 individual IDs, with 31+6 shared synthetic data."""
    focal=["FC"+f"{i:02d}" for i in range(37)]
    out=list(sorted(m.OUTGROUP_CODES))
    species=focal+out
    source_lines=[]
    individual=[]
    for idx,code in enumerate(species):
        # 6 one-sample focal codes, all source-OUT sample counts fixed at 3.
        count=1 if 31<=idx<37 else 3
        for j in range(1,count+1):
            sample=code+str(j)
            individual.append(sample)
            source_lines.append(sample+"\t"+("OUT" if code in m.OUTGROUP_CODES else code))
    assert len(individual)==117
    assert len(species)==43
    mapping,outgroups=m.sample_map("\n".join(source_lines)+"\n")
    tree1=parse_newick(binary_nwk(species)+";")
    tree2=parse_newick(binary_nwk(individual)+";")
    measured=set(focal[:31]+out)
    # External 43 species has six focal additional not measured; 2023 contains
    # measured 37 plus 22 other four-character symbols.
    other=["U"+f"{i:03d}" for i in range(22)]
    original=parse_newick("(BROW:0.2,"+binary_nwk(sorted(measured)+other)+");")
    current={}
    gene_names=set(x for genes in m.MODELS.values() for x in genes)
    for idx,code in enumerate(sorted(measured)+other):
        current[code]={
            "genes":{g:math.log1p(idx/12+(j+1)/7) for j,g in enumerate(sorted(gene_names))},
            "responses":{
                "anthocyanidins_log1p_total":math.log1p(1+(idx%10)),
                "flavonols_log1p_total":math.log1p(1+(idx%14)),
            }
        }
    assert len(current)==59
    return tree1,tree2,mapping,outgroups,original,current


class V28PublishedAlternativesAdversarial(unittest.TestCase):
    def test_external_pin_noneditable_and_author_tree_shapes(self):
        self.assertEqual(len(m.LATER_PIN),40)
        self.assertEqual(len(m.EXTERNAL),3)
        self.assertEqual(m.EXPECTED_SHARED_FOCAL,31)
        self.assertEqual(m.EXPECTED_SHARED_ALL,37)
        self.assertEqual(m.MAX_CLADE_SIZES,(8,15))
        self.assertEqual(m.ROOT_PRIMARY,"HUTE")
        for rec in m.EXTERNAL.values():
            self.assertTrue(m.pinned_url(rec).startswith(
                "https://raw.githubusercontent.com/pedrohpezzi/"))
            self.assertEqual(len(rec["git_blob_sha1"]),40)

    def test_git_blob_SHA_hashes_actual_bytes_including_prefix(self):
        data=b"Original published tree data\n"
        expected=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\x00"+data).hexdigest()
        self.assertEqual(m.git_blob_sha(data),expected)
        self.assertNotEqual(m.git_blob_sha(data+b"x"),expected)
        with self.assertRaisesRegex(TreeGateError,"BLOB_BYTES_REQUIRED"):
            m.git_blob_sha("not bytes")

    def test_wrong_original_external_bytes_fail_closed(self):
        def malformed(url):
            return b"(source bad);"
        with self.assertRaisesRegex(TreeGateError,"EXTERNAL_PUBLISHED_BYTE_LENGTH_MISMATCH"):
            m.load_published(malformed)
        class Meta:
            pass
        original=m.EXTERNAL["astral_species_tree"]
        modified=b"x"*original["bytes"]
        def fake_exact_length(url):
            return modified
        with self.assertRaisesRegex(TreeGateError,"EXTERNAL_PUBLISHED_BLOB_SHA_MISMATCH"):
            m.load_published(fake_exact_length)

    def test_speciesmap_exact_original_IDs_and_out_codes(self):
        astral,iq,mapping,out,orig,current=fake_source()
        self.assertEqual(len(mapping),117)
        self.assertEqual(len(set(mapping.values())),43)
        self.assertEqual(out,m.OUTGROUP_CODES)
        self.assertEqual(m.published_tree_admission(astral,iq,mapping,out)[
            "original_2024_IQtree_individual_tips"],117)

    def test_duplicate_or_missing_sample_mapping_fails(self):
        astral,iq,mapping,out,orig,current=fake_source()
        d=copy.deepcopy(mapping)
        d.pop(next(iter(d)))
        with self.assertRaisesRegex(TreeGateError,"IQTREE_INDIVIDUAL_SAMPLE_MAP_MISMATCH"):
            m.published_tree_admission(astral,iq,d,out)
        with self.assertRaisesRegex(TreeGateError,"SOURCE_SAMPLE_MAP_DUPLICATE"):
            m.sample_map("FC001\tFC00\nFC001\tFC00\n")
        with self.assertRaisesRegex(TreeGateError,"SOURCE_SAMPLE_MAP_SPECIES_DISAGREEMENT"):
            m.sample_map("FC001\tWRONG\n")

    def test_2024_ASTAL_and_IQtree_not_same_observation_grain(self):
        astral,iq,mapping,out,orig,current=fake_source()
        self.assertEqual(len(astral.tips()),43)
        self.assertEqual(len(iq.tips()),117)
        for rank in ("first","middle","last"):
            tr,select=m.collapse_to_one_sample_per_species(iq,mapping,rank)
            self.assertEqual(len(tr.tips()),43)
            self.assertEqual(set(tr.tips()),set(astral.tips()))
            self.assertEqual(len(select),43)
        with self.assertRaisesRegex(TreeGateError,"UNDECLARED_INDIVIDUAL_SELECTION"):
            m.collapse_to_one_sample_per_species(iq,mapping,"best_predictive")

    def test_wrong_2023_species_overlap_blocks_analysis(self):
        astral,iq,mapping,out,orig,current=fake_source()
        trees,panel=m.panels_from_sources(current,orig,astral,iq,mapping,out)
        self.assertEqual(len(panel["primary_common_focal_codes_sorted"]),31)
        self.assertEqual(len(panel["secondary_common_codes_sorted"]),37)
        d=copy.deepcopy(current)
        del d["FC00"]
        with self.assertRaisesRegex(TreeGateError,"REGISTERED_31_FOCAL_37_SHARED_PANELS_DRIFT"):
            m.panels_from_sources(d,orig,astral,iq,mapping,out)

    def test_strict_source_outgroup_anchors_and_reroot_coverage(self):
        astral,iq,mapping,out,orig,current=fake_source()
        trees,panel=m.panels_from_sources(current,orig,astral,iq,mapping,out)
        for key,tr in trees.items():
            self.assertTrue(set(panel["secondary_common_codes_sorted"]).issubset(
                set(tr.tips()) & set(current)))
            if key != m.ORIGINAL_TREE:
                self.assertEqual(set(tr.tips()) & set(current),
                                 set(panel["secondary_common_codes_sorted"]))
            for anchor in (m.ROOT_PRIMARY,m.ROOT_SECONDARY):
                rer=m.reroot_outgroup_edge(tr,anchor)
                self.assertEqual(set(rer.tips()),set(tr.tips()))
                subset=m.prune_to_taxa(rer,set(panel["primary_common_focal_codes_sorted"]))
                self.assertEqual(len(subset.tips()),31)
                self.assertEqual(set(subset.tips()),
                                 set(panel["primary_common_focal_codes_sorted"]))
        with self.assertRaisesRegex(TreeGateError,"ROOT_OUTGROUP_NOT_PRESENT"):
            m.reroot_outgroup_edge(trees[m.ORIGINAL_TREE],"NOT_PRESENT")

    def test_all_fixed_RNA_modules_and_topological_blocks_for_primary(self):
        astral,iq,mapping,out,orig,current=fake_source()
        trees,panel=m.panels_from_sources(current,orig,astral,iq,mapping,out)
        rr=m.evaluate_tree_panel(current,trees[m.ASTRAL_TREE],
            panel["primary_common_focal_codes_sorted"],m.ROOT_PRIMARY,"synthetic")
        self.assertEqual(set(rr),{"max_clade_8","max_clade_15"})
        for v in rr.values():
            self.assertEqual(v["n_species_tested_once"],31)
            self.assertEqual(len(v["endpoints"]),2)
            self.assertEqual(sum(v["fold_sizes_sorted"]),31)
            for endpoint in v["endpoints"].values():
                self.assertEqual(set(endpoint["models"]),set(m.MODELS))
                self.assertTrue(all(z["test_species_exactly_once"]==31
                    for z in endpoint["models"].values()))

    def test_test_chemical_outcome_does_not_select_the_tree_or_folds(self):
        astral,iq,mapping,out,orig,current=fake_source()
        trees,panel=m.panels_from_sources(current,orig,astral,iq,mapping,out)
        sub=m.prune_to_taxa(m.reroot_outgroup_edge(trees[m.ASTRAL_TREE],
                             m.ROOT_PRIMARY),
                             set(panel["primary_common_focal_codes_sorted"]))
        first=m.exact_monophyletic_partitions(sub,8)
        d=copy.deepcopy(current)
        for species in d.values():
            species["responses"]["anthocyanidins_log1p_total"]+=1000
        second=m.exact_monophyletic_partitions(sub,8)
        self.assertEqual(first,second)
        self.assertEqual(sum(map(len,first)),31)

    def test_no_selective_posthoc_primary_win_gate(self):
        astral,iq,mapping,out,orig,current=fake_source()
        trees,panel=m.panels_from_sources(current,orig,astral,iq,mapping,out)
        v={name:m.evaluate_tree_panel(current,trees[name],
             panel["primary_common_focal_codes_sorted"],m.ROOT_PRIMARY,name)
             for name in (m.ORIGINAL_TREE,m.ASTRAL_TREE,m.IQ_FIRST_TREE)}
        normal=m.positive_primary_gate(v)
        self.assertEqual(normal["predeclared_total_direction_checks"],12)
        d=copy.deepcopy(v)
        d[m.ASTRAL_TREE]["max_clade_8"]["endpoints"][
            "flavonols_log1p_total"]["models"]["branch_three"][
                "RMSE_relative_to_same_fold_training_mean"]=1.0001
        failed=m.positive_primary_gate(d)
        self.assertEqual(failed["verdict"],
                         "NOT_UNIFORMLY_PORTABLE_ACROSS_PUBLISHED_TREE_DEFINITIONS")
        self.assertTrue(any("flavonols" in str(x)
                            for x in failed["failed_predeclared_gates"]))

if __name__=="__main__":
    unittest.main()
