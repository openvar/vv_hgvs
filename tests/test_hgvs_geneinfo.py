# -*- coding: utf-8 -*-
"""Tests vvta gene info code"""

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

hdp = vvhgvs.dataproviders.uta.connect(mode=None, cache=None)
class Test_UTA_Geneinfo(unittest.TestCase):
    def test_get_null_gene_info(self):
        test = "ThisShouldNotGetAGeneinfoResult"
        null_result = hdp.get_gene_info(test)
        self.assertEqual(None, null_result)
        null_result = hdp.get_gene_info_by_id(test)
        self.assertEqual(None, null_result)
        null_result = hdp.get_gene_info_by_alias(test)
        # we have no assurances that old symbols are unique
        # current ones and ids must be, so return array not
        # single result
        self.assertEqual(0, len(null_result))
    def test_get_gene_info(self):
        gene_info = hdp.get_gene_info("VHL")
        self.assertEqual("VHL", gene_info["hgnc"])
        self.assertEqual("HGNC:12687",gene_info["hgnc_id"])
        self.assertEqual("3p25.3", gene_info["maploc"])
        self.assertEqual(7, len(gene_info))
    def test_get_gene_info_by_id(self):
        gene_info_id = hdp.get_gene_info_by_id("HGNC:12687")
        gene_info_symbol = hdp.get_gene_info("VHL")
        self.assertEqual(gene_info_id,gene_info_symbol)
        self.assertEqual("VHL", gene_info_id["hgnc"])
        self.assertEqual("HGNC:12687", gene_info_id["hgnc_id"])
        self.assertEqual("3p25.3", gene_info_id["maploc"])
        self.assertEqual(7, len(gene_info_id))
    def test_get_gene_info_by_alias(self):
        # alias is actually historic symbol in new data
        gene_info = hdp.get_gene_info_by_alias("VHLP")
        gene_info = gene_info[0]
        self.assertEqual("VHLL", gene_info["hgnc"])
        self.assertEqual('HGNC:30666',gene_info["hgnc_id"])
        self.assertEqual("1q22", gene_info["maploc"])
        self.assertEqual(7, len(gene_info))
    def test_get_gene_info_multi_alias(self):
        # test get gene with 3 alias (or more) for end vs middle check
        #from hgnc gene file 2021/2
        #HGNC:34 ABCA4   ATP binding cassette subfamily A member 4       protein-coding gene     gene with protein product       Approved        1p22.1  01p22.1 "FFM|ARMD2|CORD3"       Stargardt disease       "STGD1|ABCR|RP19|STGD" ......
        gene_info_id = hdp.get_gene_info_by_id("HGNC:34")
        # get by id then loop to keep fwd compatiblity on changes
        print(gene_info_id['aliases'])
        for alias in re.split(',',gene_info_id['aliases']):
            result = hdp.get_gene_info_by_alias(alias)
            assert(result)
            assert(result[0] == gene_info_id)
    def test_get_gene_info_alias_like_other_curr(self):
        # test get gene by alias where alias=other gene's current
        # symbol. HTT was old symbol of HGNC:11050 now SLC6A4
        # HTT now symbol for, HGNC:4851 huntingtin
        gene_info_htt_old = hdp.get_gene_info_by_id("HGNC:11050")
        gene_info_htt_curr = hdp.get_gene_info_by_id("HGNC:4851")
        gene_info_htt_alias = hdp.get_gene_info_by_alias('HTT')
        gene_info_htt_symbol = hdp.get_gene_info('HTT')

        assert(gene_info_htt_alias[0] != gene_info_htt_symbol)
        self.assertEqual(gene_info_htt_old,gene_info_htt_alias[0])
        self.assertEqual(gene_info_htt_curr,gene_info_htt_symbol)
        self.assertEqual(gene_info_htt_curr['aliases'],'HD')

    def test_get_gene_info_multi_alias(self):
        # test fetch of mutiple results, where historic alias was used
        # for multiple genes, using AMY1 (amalyse alpha 1) set
        results = hdp.get_gene_info_by_alias('AMY1')
        for result in results:
            assert(result['hgnc'] in ['AMY1A','AMY1B','AMY1C'])
    def test_get_gene_info_alias_known_changes(self):
        # test a set of genes known to have changed alias
        symbol_change_data = [
            {'hgnc_id':'HGNC:4720','hgnc_symbol':'H1-6','old_symbol':'HIST1H1T'},
            {'hgnc_id':'HGNC:4775','hgnc_symbol':'H3C10','old_symbol':'HIST1H3H'},
            {'hgnc_id':'HGNC:4776','hgnc_symbol':'H3C2','old_symbol':'HIST1H3B'},
            {'hgnc_id':'HGNC:4768','hgnc_symbol':'H3C3','old_symbol':'HIST1H3C'},
            {'hgnc_id':'HGNC:20503','hgnc_symbol':'H3C14','old_symbol':'HIST2H3C'},
            {'hgnc_id':'HGNC:4764','hgnc_symbol':'H3-3A','old_symbol':'H3F3A'},
            {'hgnc_id':'HGNC:4765','hgnc_symbol':'H3-3B','old_symbol':'H3F3B'},
            {'hgnc_id':'HGNC:20510','hgnc_symbol':'H4-16','old_symbol':'HIST4H4'}]
        for data in symbol_change_data:
            alias_dat = hdp.get_gene_info_by_alias(data['old_symbol'])
            cur_dat = hdp.get_gene_info(data['hgnc_symbol'])
            self.assertEqual(alias_dat[0],cur_dat)
            self.assertEqual(alias_dat[0]['hgnc_id'],data['hgnc_id'])
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
