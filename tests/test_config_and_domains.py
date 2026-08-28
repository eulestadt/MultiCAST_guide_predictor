import unittest

import numpy as np

from structure_ko.config import StructureConfig, load_config
from structure_ko.structure import domains_from_pae


class ConfigTests(unittest.TestCase):
    def test_full_yaml_loads(self):
        cfg = load_config("examples/config.full.yaml")
        self.assertIn("lacZ", cfg.genes)
        self.assertEqual(cfg.cast.insertion_offset_bp, 49)
        self.assertEqual(cfg.cast.oligo_left, "TACTACTGCAAAGTAGCTGATAAC")
        self.assertEqual(cfg.organism.preset, "ecoli_k12")

    def test_cli_genes_override(self):
        cfg = load_config(genes=["lacZ"], organism="ecoli_k12")
        self.assertEqual(cfg.genes, ["lacZ"])
        self.assertTrue(str(cfg.organism.genome).endswith(".fna") or cfg.organism.genome.endswith(".fna"))


class DomainSplitTests(unittest.TestCase):
    def test_two_blocks_split_on_high_inter_pae(self):
        n = 120
        pae = np.full((n, n), 20.0)
        # two tight domains 0-50 and 70-120
        pae[0:50, 0:50] = 2.0
        pae[70:120, 70:120] = 2.0
        plddt = np.full(n, 90.0)
        plddt[50:70] = 40.0
        cfg = StructureConfig(min_domain_aa=30, pae_split_delta=4.0)
        domains = domains_from_pae(pae, plddt, cfg)
        self.assertGreaterEqual(len(domains), 2)
        self.assertLess(domains[0].end, 75)


if __name__ == "__main__":
    unittest.main()
