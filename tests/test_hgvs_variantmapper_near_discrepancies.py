# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os

import pytest

from vvhgvs.exceptions import HGVSError, HGVSInvalidVariantError
import vvhgvs.dataproviders.uta
import vvhgvs.parser
import vvhgvs.assemblymapper
from support import CACHE

# 14 tests disabled due to revert of partly fixed gap issues, variant
# validator handles this for now so if you intend to fix both need to
# change at once
tests_fn = "tests/data/proj-near-disc.tsv"


def read_tests(fn):
    """read tests from tsv file, return iterator of dicts"""
    fh = open(fn)
    for line in fh:
        line = line.strip()
        if line == "":
            continue
        if line.startswith("#"):
            continue
        dt, lt, var, exp = line.split()
        yield {"disc_type": dt, "loc_type": lt, "variant": var, "expected": exp}

mode_txt = os.environ.get("HGVS_CACHE_MODE", None)
hp = vvhgvs.parser.Parser()
hdp = vvhgvs.dataproviders.uta.connect(mode=mode_txt, cache=CACHE)
# TODO: Use variantmapper instead of assemblymapper
am38 = vvhgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh38', normalize=False)

tests = list(read_tests(tests_fn))
params = [(t["variant"], t["expected"], t["loc_type"] + " " + t["disc_type"]) for t in tests]


@pytest.mark.pnd
@pytest.mark.parametrize("variant,expected,description", params)
def test_projection_near_discrepancies(variant, expected, description):
    var_n = hp.parse_hgvs_variant(variant)
    var_g = am38.n_to_g(var_n)
    assert expected == str(var_g), description


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
