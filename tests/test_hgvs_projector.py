# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os

import unittest

import pytest

import vvhgvs.dataproviders.uta

import vvhgvs.location
import vvhgvs.parser
import vvhgvs.projector
from support import CACHE

mode_txt = os.environ.get("HGVS_CACHE_MODE", None)

class TestHgvsProjector(unittest.TestCase):
    @classmethod
    def setUp(cls):
        cls.hdp = vvhgvs.dataproviders.uta.connect(mode=mode_txt, cache=CACHE)
        cls.alt_ac = "NC_000001.10"
        cls.alt_aln_method = "splign"
        cls.hp = vvhgvs.parser.Parser()

    def tst_forward_and_backward(self, v1, v2):
        pj = vvhgvs.projector.Projector(self.hdp, self.alt_ac, v1.ac, v2.ac, self.alt_aln_method, self.alt_aln_method)
        self.assertEqual(pj.project_variant_forward(v1), v2)
        self.assertEqual(pj.project_variant_backward(v2), v1)

    @pytest.mark.quick
    def test_rs201430561(self):
        # rs201430561 http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=201430561
        hgvs_c = ["NM_001197320.1:c.281C>T", "NM_021960.4:c.740C>T", "NM_182763.2:c.688+403C>T"]
        var_c = [self.hp.parse_hgvs_variant(h) for h in hgvs_c]
        self.tst_forward_and_backward(var_c[0], var_c[1])
        self.tst_forward_and_backward(var_c[0], var_c[2])
        self.tst_forward_and_backward(var_c[1], var_c[2])

    @pytest.mark.quick
    def test_bad_acs(self):
        hgvs_c = ["NM_001197320.1:c.281C>T", "NM_021960.4:c.740C>T", "NM_182763.2:c.688+403C>T"]
        var_c = [self.hp.parse_hgvs_variant(h) for h in hgvs_c]
        pj = vvhgvs.projector.Projector(self.hdp, self.alt_ac, var_c[0].ac, var_c[1].ac, self.alt_aln_method,
                                        self.alt_aln_method)
        # intentionally call p_v_f with variant on *destination* transcript, and vice versa
        with self.assertRaises(RuntimeError):
            pj.project_variant_forward(var_c[1])
        with self.assertRaises(RuntimeError):
            pj.project_variant_backward(var_c[0])


if __name__ == "__main__":
    unittest.main()
