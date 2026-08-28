import unittest

from structure_ko.cast import scan_guides
from structure_ko.config import CastConfig, DEFAULT_OLIGO_LEFT, DEFAULT_OLIGO_RIGHT
from structure_ko.gene import GeneRecord


def _record(cds: str) -> GeneRecord:
    return GeneRecord(
        query="fake",
        name="fake",
        locus_tag="b0000",
        seqid="NC_TEST",
        start=1,
        end=len(cds),
        strand="+",
        cds=cds,
        protein="A" * (len(cds) // 3),
        uniprot="P00000",
        protein_id=None,
        product="test",
        organism="test",
    )


class CastOffsetTests(unittest.TestCase):
    def test_coding_strand_insertion_is_offset_past_guide_3prime(self):
        # 80 N's, then ACC + CC PAM + 32xA guide, then more N's
        prefix = "N" * 80
        match = "ACC" + "CC" + ("A" * 32)
        cds = prefix + match + ("N" * 200)
        cfg = CastConfig(insertion_offset_bp=49)
        guides = [g for g in scan_guides(_record(cds), cfg) if g.strand == "coding" and g.guide_sequence == "A" * 32]
        self.assertTrue(guides)
        g = guides[0]
        self.assertEqual(g.pam_region, "ACCCC")
        self.assertEqual(g.guide_3prime_cds, 80 + 36)
        self.assertEqual(g.insertion_nt, 80 + 36 + 49)
        self.assertEqual(g.insertion_aa, (80 + 36 + 49) // 3)
        self.assertTrue(g.oligo.startswith(DEFAULT_OLIGO_LEFT))
        self.assertTrue(g.oligo.endswith(DEFAULT_OLIGO_RIGHT))
        self.assertEqual(len(g.oligo), len(DEFAULT_OLIGO_LEFT) + 32 + len(DEFAULT_OLIGO_RIGHT))

    def test_template_strand_insertion_moves_toward_n_terminus(self):
        from Bio.Seq import Seq

        # Place the same PAM/guide on the reverse complement so scan_guides
        # reports a template-strand hit.
        rc_match = "ACC" + "CC" + ("T" * 32)
        prefix_rc = "N" * 40
        suffix_rc = "N" * 120
        rc = prefix_rc + rc_match + suffix_rc
        cds = str(Seq(rc).reverse_complement())
        cfg = CastConfig(insertion_offset_bp=49)
        hits = [g for g in scan_guides(_record(cds), cfg) if g.strand == "template" and g.guide_sequence == "T" * 32]
        self.assertTrue(hits)
        g = hits[0]
        self.assertLess(g.insertion_nt, g.guide_3prime_cds)


if __name__ == "__main__":
    unittest.main()
