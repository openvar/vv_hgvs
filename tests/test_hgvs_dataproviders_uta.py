# -*- coding: utf-8 -*-
"""Tests uta postgresql client"""

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import re
import unittest

from vvhgvs.exceptions import HGVSDataNotAvailableError
import vvhgvs.dataproviders.uta
import vvhgvs.edit
import vvhgvs.location
import vvhgvs.posedit
import vvhgvs.variantmapper
import vvhgvs.sequencevariant
from support import CACHE

mode_txt = os.environ.get("HGVS_CACHE_MODE", None)
class UTA_Base(object):
    def test_get_acs_for_protein_seq(self):
        exp = ["NP_001005405.1", "MD5_8fc09b1d9a38a8c55176a0fa922df227"]
        s = """
        mgccgcsggc gsgcggcgsg sggcgsgcgg cgssccvpic cckpvcccvp acscsscgsc
        ggskggcgsc gsskggcgsc gcsqsncckp ccsssgcgsf ccqsscskpc ccqssccqss
        cckpcccqss ccqsscfkpc ccqssccvpv ccqcki
        """

        s = re.sub(r"\s+", "", s.upper())
        curr_accs = self.hdp.get_acs_for_protein_seq(s)
        for prot_ac in exp:
            assert prot_ac in curr_accs

        exp = ["NP_071928.2", "MD5_ffb0d4adbd5e0b5d71678228b3696984"]
        s = """
        masetektha llqtcstesl isslglgafc lvadrllqfs tiqqndwlra lsdnavhcvi
        gmwswavvtg ikkktdfgei ilagflasvi dvdhfflags mslkaaltlp rrpflhcstv
        ipvvvltlkf tmhlfklkds wcflpwmlfi swtshhirdg irhglwicpf gktsplpfwl
        yviitsslph icsfvmyltg trqmmsskhg vridv
        """

        s = re.sub(r"\s+", "", s.upper())
        curr_accs = self.hdp.get_acs_for_protein_seq(s)
        for prot_ac in exp:
            assert prot_ac in curr_accs

    def test_get_gene_info(self):
        gene_info = self.hdp.get_gene_info("VHL")
        self.assertEqual("VHL", gene_info["hgnc"])
        self.assertEqual("3p25.3", gene_info["maploc"])
        assert len(gene_info) in [6,7] #cope with gene id in newer databases
    def test_get_gene_info_by_id(self):
        gene_info_id = self.hdp.get_gene_info_by_id("HGNC:12687")
        gene_info_symbol = self.hdp.get_gene_info("VHL")
        self.assertEqual(gene_info_id,gene_info_symbol)
        self.assertEqual("VHL", gene_info_id["hgnc"])
        self.assertEqual("HGNC:12687", gene_info_id["hgnc_id"])
        self.assertEqual("3p25.3", gene_info_id["maploc"])
        self.assertEqual(7, len(gene_info_id))


    def test_get_tx_exons(self):
        tx_exons = self.hdp.get_tx_exons("NM_000551.3", "NC_000003.11", "splign")
        self.assertEqual(3, len(tx_exons))

    def test_get_tx_exons_invalid_tx_ac(self):
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.get_tx_exons("NM_999999.9", "NC_000003.11", "splign")

    def test_get_tx_exons_invalid_alt_ac(self):
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.get_tx_exons("NM_000551.3", "NC_000999.9", "splign")

    def test_get_tx_exons_invalid_alt_aln_method(self):
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.get_tx_exons("NM_000551.3", "NC_000999.9", "best")

    def test_get_tx_for_gene(self):
        tig_id = self.hdp.get_tx_for_gene("VHL")
        assert(len(tig_id)> 14)
    def test_get_tx_for_gene_id(self):
        tig_symbol = self.hdp.get_tx_for_gene("VHL")
        tig_id = self.hdp.get_tx_for_gene_id("HGNC:12687")
        for id_tx_set in tig_id:
            assert id_tx_set in tig_symbol
        assert(len(tig_id)> 14)

    def test_get_tx_for_gene_invalid_gene(self):
        tig = self.hdp.get_tx_for_gene("GENE")
        self.assertEqual(0, len(tig))

    def test_get_tx_info(self):
        #replace AC_000143.1 (HuRef-removed) with NC_000011.10 (GRCh38)
        tx_info = self.hdp.get_tx_info("NM_000051.3", "NC_000011.10", "splign")
        self.assertEqual(385, tx_info["cds_start_i"])
        self.assertEqual(9556, tx_info["cds_end_i"])
        self.assertEqual("NC_000011.10", tx_info["alt_ac"])

    def test_get_tx_info_invalid_tx_ac(self):
        #replace AC_000143.1 (HuRef-removed) with NC_000011.10 (GRCh38)
        # we want to test a non-existent transcript in an existing chr
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.get_tx_info("NM_999999.9", "NC_000011.10", "splign")

    def test_get_tx_mapping_options(self):
        tx_mapping_options = self.hdp.get_tx_mapping_options("NM_000551.3")
        self.assertIn(["NM_000551.3", "NC_000003.11", "splign"], tx_mapping_options)
        self.assertIn(["NM_000551.3", "NC_000003.11", "blat"], tx_mapping_options)

    def test_get_tx_mapping_options_invalid(self):
        tx_info_options = self.hdp.get_tx_mapping_options("NM_999999.9")
        self.assertEqual(tx_info_options, [])


class Test_hgvs_dataproviders_uta_UTA_default(unittest.TestCase, UTA_Base):
    @classmethod
    def setUpClass(cls):
        cls.hdp = vvhgvs.dataproviders.uta.connect(mode=mode_txt, cache=CACHE)


class Test_hgvs_dataproviders_uta_UTA_default_with_pooling(unittest.TestCase, UTA_Base):
    @classmethod
    def setUpClass(cls):
        cls.hdp = vvhgvs.dataproviders.uta.connect(
            pooling=True, mode=mode_txt, cache=CACHE)


class TestUTACache(Test_hgvs_dataproviders_uta_UTA_default):
    def _create_cdna_variant(self):
        start = vvhgvs.location.SimplePosition(118898437)
        end = vvhgvs.location.SimplePosition(118898437)
        iv = vvhgvs.location.Interval(start=start, end=end)
        edit = vvhgvs.edit.NARefAlt(ref="C", alt="T")
        posedit = vvhgvs.posedit.PosEdit(pos=iv, edit=edit)
        genomic_variant = vvhgvs.sequencevariant.SequenceVariant(
            ac="NC_000011.9",
            type="g",
            posedit=posedit,
        )
        variantmapper = vvhgvs.variantmapper.VariantMapper(self.hdp)
        return variantmapper.g_to_c(genomic_variant, "NM_001164277.1")

    def test_deterministic_cache_results(self):
        """
        Check that identical request to the UTA yields the same results.
        """
        var1 = self._create_cdna_variant()
        var2 = self._create_cdna_variant()
        self.assertEqual(str(var1), str(var2))


if __name__ == "__main__":
    unittest.main()

# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
