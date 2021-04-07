# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import pprint
import re
import sys
import os

import unittest

import pytest

from vvhgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSParseError, HGVSInvalidVariantError, HGVSInvalidVariantError
from vvhgvs.enums import Datum
import vvhgvs.assemblymapper
import vvhgvs.dataproviders.uta
import vvhgvs.normalizer
import vvhgvs.parser
import vvhgvs.sequencevariant
import vvhgvs.validator
import vvhgvs.variantmapper
from support import CACHE
mode_txt = os.environ.get("HGVS_CACHE_MODE", None)

@pytest.mark.issues
class Test_Issues(unittest.TestCase):
    def setUp(self):
        self.hdp = vvhgvs.dataproviders.uta.connect(mode=mode_txt, cache=CACHE)
        self.vm = vvhgvs.variantmapper.VariantMapper(self.hdp, replace_reference=False)
        self.vm_rr = vvhgvs.variantmapper.VariantMapper(self.hdp, replace_reference=True)
        self.hp = vvhgvs.parser.Parser()
        self.hn = vvhgvs.normalizer.Normalizer(self.hdp)
        self.hv = vvhgvs.validator.IntrinsicValidator()
        self.am37 = vvhgvs.assemblymapper.AssemblyMapper(
            self.hdp, replace_reference=True, assembly_name='GRCh37', alt_aln_method='splign')
        self.am38 = vvhgvs.assemblymapper.AssemblyMapper(
            self.hdp, replace_reference=True, assembly_name='GRCh38', alt_aln_method='splign')
        self.vn = vvhgvs.normalizer.Normalizer(self.hdp, shuffle_direction=3, cross_boundaries=True)

    def test_424_430_nochange_parse_and_format(self):
        h = "NM_012.3:c.1="
        v = self.hp.parse_hgvs_variant(h)
        self.assertEqual("NM_012.3:c.1=", str(v))
        self.assertEqual("", v.posedit.edit.ref)
        self.assertEqual("", v.posedit.edit.alt)

        h = "NM_012.3:c.1A="
        v = self.hp.parse_hgvs_variant(h)
        self.assertEqual("NM_012.3:c.1=", str(v))
        self.assertEqual("A", v.posedit.edit.ref)
        self.assertEqual("A", v.posedit.edit.alt)

    def test_459_exception_when_ac_nonexistent(self):
        bogus_ac = "NM_000000.99"
        with self.assertRaises(HGVSDataNotAvailableError):
            self.hdp.seqfetcher.fetch_seq(bogus_ac)

    def test_499_whole_gene_dup(self):
        # Verify that 1_*1dup works
        self.am37.c_to_p(self.hp.parse_hgvs_variant("NM_001637.3:c.1_*1dup"))

        # Now try -1_*1dup (essentially, this is issue #499)
        self.am37.c_to_p(self.hp.parse_hgvs_variant("NM_001637.3:c.-1_*1dup"))
