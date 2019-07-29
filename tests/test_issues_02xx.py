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


@pytest.mark.issues
class Test_Issues(unittest.TestCase):
    def setUp(self):
        self.hdp = vvhgvs.dataproviders.uta.connect(mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE)
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

    def test_260_raise_exception_when_mapping_bogus_variant(self):
        v = self.hp.parse_hgvs_variant("NM_000059.3:c.7790delAAG")
        with self.assertRaises(HGVSInvalidVariantError):
            self.am37.c_to_p(v)

    def test_285_partial_palindrome_inversion(self):
        # https://github.com/biocommons/hgvs/issues/285/
        # Inversion mapping is not palindrome aware
        v = self.hp.parse_hgvs_variant("NM_000088.3:c.589_600inv")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000088.3:c.590_599inv")

    def test_293_parser_attribute_assignment_error(self):
        # https://github.com/biocommons/hgvs/issues/293/
        var = self.hp.parse_hgvs_variant('NG_029146.1:g.6494delG')
        self.vn.normalize(var)    # per issue, should raise error
