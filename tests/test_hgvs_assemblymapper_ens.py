# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os

import unittest

import pytest

from vvhgvs.exceptions import HGVSInvalidVariantError
import vvhgvs.dataproviders.uta
import vvhgvs.parser
import vvhgvs.variantmapper
import vvhgvs.assemblymapper
from support import CACHE

mode_txt = os.environ.get("HGVS_CACHE_MODE", None)
@pytest.mark.quick
class Test_VariantMapper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hdp = vvhgvs.dataproviders.uta.connect(mode=mode_txt, cache=CACHE)
        cls.am = vvhgvs.assemblymapper.AssemblyMapper(cls.hdp, alt_aln_method='genebuild')
        cls.am37 = vvhgvs.assemblymapper.AssemblyMapper(cls.hdp, assembly_name="GRCh37", alt_aln_method='genebuild')
        cls.hp = vvhgvs.parser.Parser()
    # for now do mixed test ENST00000617537(.4) in GRCh38 but not in 82 in any form
    def test_VariantMapper_quick(self):
        # From garcia.tsv:
        hgvs_g = 'NC_000007.14:g.36522056C>T'#"NC_000007.13:g.36561662C>T"
        hgvs_c = "ENST00000617537.5:c.1582G>A"#ENST00000380152.7
        hgvs_p = "ENSP00000483783.1:p.(Gly528Arg)"

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.am.g_to_c(var_g, "ENST00000617537.5")
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_c), hgvs_c)
        self.assertEqual(str(var_p), hgvs_p)

        with self.assertRaises(HGVSInvalidVariantError):
            self.am.c_to_p(self.hp.parse_hgvs_variant("ENSP00000369497.3:c.7790delAAG"))

        hgvs_c_37 = "ENST00000380152.3:c.7791A>G"
        hgvs_c_38 = "ENST00000380152.7:c.7791A>G"
        hgvs_p = "ENSP00000369497.3:p.(Lys2597=)"

        var_c_37 = self.hp.parse_hgvs_variant(hgvs_c_37)
        var_c_38 = self.hp.parse_hgvs_variant(hgvs_c_38)
        var_p_37 = self.am.c_to_p(var_c_37)
        var_p_38 = self.am.c_to_p(var_c_38)

        self.assertEqual(str(var_p_37), hgvs_p)
        self.assertEqual(str(var_p_38), hgvs_p)

        hgvs_c_37 = "ENST00000196061.4:c.1594_1596del"
        hgvs_c_38 = "ENST00000196061.5:c.1594_1596del"
        hgvs_p = "ENSP00000196061.4:p.(Glu532del)"

        var_c_37 = self.hp.parse_hgvs_variant(hgvs_c_37)
        var_c_38 = self.hp.parse_hgvs_variant(hgvs_c_38)
        var_p_37 = self.am.c_to_p(var_c_37)
        var_p_38 = self.am.c_to_p(var_c_38)

        self.assertEqual(str(var_p_37), hgvs_p)
        self.assertEqual(str(var_p_38), hgvs_p)

        hgvs_c_37 = "ENST00000304636.3:c.2490_2516del"
        hgvs_c_38 = "ENST00000304636.8:c.2490_2516del"
        hgvs_p = "ENSP00000304408.3:p.(Glu832_Gly840del)"

        var_c_37 = self.hp.parse_hgvs_variant(hgvs_c_37)
        var_c_38 = self.hp.parse_hgvs_variant(hgvs_c_38)
        var_p_37 = self.am.c_to_p(var_c_37)
        var_p_38 = self.am.c_to_p(var_c_38)

        self.assertEqual(str(var_p_37), hgvs_p)
        self.assertEqual(str(var_p_38), hgvs_p)

        hgvs_c_37 = "ENST00000511181.1:c.376="
        hgvs_c_38 = "ENST00000511181.5:c.376="
        hgvs_g_37 = "NC_000005.9:g.74011484="
        hgvs_g_38 = "NC_000005.10:g.74715659="


        var_c_37 = self.hp.parse_hgvs_variant(hgvs_c_37)
        var_c_38 = self.hp.parse_hgvs_variant(hgvs_c_38)
        var_p_37 = self.am37.c_to_g(var_c_37)
        var_p_38 = self.am.c_to_g(var_c_38)

        self.assertEqual(str(var_p_37), hgvs_g_37)
        self.assertEqual(str(var_p_38), hgvs_g_38)

        #hgvs_c_37 = "ENST00000601016.1:c.-8_-3inv6"
        hgvs_c_38 = "ENST00000601016.6:c.-8_-3inv6"
        #actual position on NW_003871103.3
        #hgvs_g_37 = "NC_000023.11:g.153640173_153640178inv"
        hgvs_g_38 = "NC_000023.11:g.154411836_154411841inv"

        var_c_38 = self.hp.parse_hgvs_variant(hgvs_c_38)
        var_g_38 = self.am.c_to_g(var_c_38)

        self.assertEqual(str(var_g_38), hgvs_g_38)

        # no Grch37 ver
        hgvs_c = "ENST00000622030.2:c.89_91inv3"
        hgvs_g = "NC_000002.12:g.31580810_31580812inv"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_g = self.am.c_to_g(var_c)

        self.assertEqual(str(var_g), hgvs_g)

        # not in 82, may add later, with intermediates 
        hgvs_c = "ENST00000617537.5:c.1582_1583inv"
        hgvs_p = "ENSP00000483783.1:p.(Gly528Pro)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c_37 = "ENST00000261866.7:c.-20_*20inv"
        hgvs_c_38 = "ENST00000261866.12:c.-20_*20inv"
        hgvs_p = "ENSP00000261866.7:p.?"

        # ncbi RNA is longer, EBI ver triggers out of bounds error
        #var_c_37 = self.hp.parse_hgvs_variant(hgvs_c_37)
        #var_c_38 = self.hp.parse_hgvs_variant(hgvs_c_38)
        #var_p_37 = self.am.c_to_p(var_c_37)
        #var_p_38 = self.am.c_to_p(var_c_38)

        #self.assertEqual(str(var_p_37), hgvs_p)
        #self.assertEqual(str(var_p_38), hgvs_p)

        '''
    def test_projection_at_alignment_discrepancy(self):
        hgvs_g = "NC_000019.10:g.50378564_50378565insTG"
        hgvs_n = "NM_007121.5:n.796_798delinsTG"
        # NM_007121.5 not in database
        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_n = self.am.g_to_n(var_g, 'NM_007121.5')

        self.assertEqual(str(var_n), hgvs_n)

        hgvs_g = "NC_000007.14:g.149779575delC"
        hgvs_n = "NM_198455.2:n.1115_1116insAG"

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_n = self.am.g_to_n(var_g, 'NM_198455.2')

        self.assertEqual(str(var_n), hgvs_n)

        # issue-353
        hgvs_g = "NC_000012.11:g.122064775C>T"
        hgvs_c = "NM_032790.3:c.127_128insTGCCAC"

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.am.g_to_c(var_g, 'NM_032790.3')

        self.assertEqual(str(var_c), hgvs_c)

        # issue-461
        hgvs_g = "NC_000002.11:g.73675227_73675228insCTC"
        hgvs_c = "NM_015120.4:c.1574_1576="

        var_g = self.hp.parse_hgvs_variant(hgvs_g)
        var_c = self.am37.g_to_c(var_g, 'NM_015120.4')

        self.assertEqual(str(var_c), hgvs_c)

        # issue-259
        hgvs_c = "NM_000116.4:c.-120_-119insT"
        hgvs_g = "NC_000023.10:g.153640061="

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_g = self.am37.c_to_g(var_c)

        self.assertEqual(str(var_g), hgvs_g)

        hgvs_c = "NM_000348.3:c.88del"
        hgvs_g = "NC_000002.11:g.31805882_31805883="

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_g = self.am37.c_to_g(var_c)

        self.assertEqual(str(var_g), hgvs_g)
'''# this needs recovering but is an issue pre uta->vvta switch
    def test_c_to_p_with_stop_gain(self):
        # issue-474
        #ENST00000361134.2 is GRCh 37 and 38 so use 
        # ENST00000538474, from same ccds to differ
        hgvs_c_37 = "ENST00000538474.1:c.1733_1735delinsTTT"
        hgvs_c_38 = "ENST00000538474.5:c.1733_1735delinsTTT"
        hgvs_p = "ENSP00000442397.1:p.(Pro578_Lys579delinsLeuTer)"

        var_c_37 = self.hp.parse_hgvs_variant(hgvs_c_37)
        var_c_38 = self.hp.parse_hgvs_variant(hgvs_c_38)
        var_p_37 = self.am.c_to_p(var_c_37)
        var_p_38 = self.am.c_to_p(var_c_38)

        self.assertEqual(str(var_p_37), hgvs_p)
        self.assertEqual(str(var_p_38), hgvs_p)

        # issue-492
        # test switched to ENST00000378505, but this is replaced by 
        # ENSP00000495537.1 later
        # neither new or old is in 99 so leave out GRCh38 in initial
        # test run (ens 82&99 only)
        hgvs_c = "ENST00000378505.2:c.2847_2848delAGinsCT"
        hgvs_p = "ENSP00000367766.2:p.(Glu949_Glu950delinsAspTer)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c = "ENST00000378505.2:c.2847_2848inv"
        hgvs_p = "ENSP00000367766.2:p.(Glu949_Glu950delinsAspTer)"

        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.am.c_to_p(var_c)

        self.assertEqual(str(var_p), hgvs_p)

        hgvs_c_37 = "ENST00000538474.1:c.1735A>T"
        hgvs_c_38 = "ENST00000538474.5:c.1735A>T"
        hgvs_p = "ENSP00000442397.1:p.(Lys579Ter)"

        var_c_37 = self.hp.parse_hgvs_variant(hgvs_c_37)
        var_c_38 = self.hp.parse_hgvs_variant(hgvs_c_38)
        var_p_37 = self.am.c_to_p(var_c_37)
        var_p_38 = self.am.c_to_p(var_c_38)

        self.assertEqual(str(var_p_37), hgvs_p)
        self.assertEqual(str(var_p_38), hgvs_p)

        hgvs_c_37 = "ENST00000538474.1:c.1795_*3delinsTAG"
        hgvs_c_38 = "ENST00000538474.5:c.1795_*3delinsTAG"
        hgvs_p = "ENSP00000442397.1:p.(Leu599Ter)"

        var_c_37 = self.hp.parse_hgvs_variant(hgvs_c_37)
        var_c_38 = self.hp.parse_hgvs_variant(hgvs_c_38)
        var_p_37 = self.am.c_to_p(var_c_37)
        var_p_38 = self.am.c_to_p(var_c_38)

        self.assertEqual(str(var_p_37), hgvs_p)
        self.assertEqual(str(var_p_38), hgvs_p)


class Test_RefReplacement(unittest.TestCase):
    test_cases = [
    # These casese attempt to test reference update in four dimensions:
    # currently only grch37 based it is also imposible to test this with ensembl data which is ref based
    # - variant type: n, c, g
    # - major mapping paths: c<->n, c<->g, n<->g
    # - variant class: sub, del, ins, delins, dup
    # - strand: +/-

    # ADRB2    │ NM_000024.5 │  239 │ 1481 │ NC_000005.9  │  1 │ 148206155,148208197 | 284=1X32=1X1724=
    # cseq = hdp.fetch_seq("NM_000024.5")
    # gseq = hdp.fetch_seq("NC_000005.9",148206155,148208197)
    # cseq[280:290] = "CAATAGAAGC"
    # gseq[280:290] = "CAATGGAAGC"
    #                      ^ @ n.285
    # These variants are in and around the first sub:
        {
            u"c": "NM_000024.5:c.42C>N",
            u"g": "NC_000005.9:g.148206436C>N",
            u"n": u"NM_000024.5:n.281C>N"
        },
        {
            u"c": "NM_000024.5:c.43A>N",
            u"g": "NC_000005.9:g.148206437A>N",
            u"n": u"NM_000024.5:n.282A>N"
        },
        {
            u"c": "NM_000024.5:c.44A>N",
            u"g": "NC_000005.9:g.148206438A>N",
            u"n": u"NM_000024.5:n.283A>N"
        },
        {
            u"c": "NM_000024.5:c.45T>N",
            u"g": "NC_000005.9:g.148206439T>N",
            u"n": u"NM_000024.5:n.284T>N"
        },
        {
            u"c": "NM_000024.5:c.46A>N",
            u"g": "NC_000005.9:g.148206440G>N",
            u"n": u"NM_000024.5:n.285A>N"
        },    # ref repl
        {
            u"c": "NM_000024.5:c.47G>N",
            u"g": "NC_000005.9:g.148206441G>N",
            u"n": u"NM_000024.5:n.286G>N"
        },
        {
            u"c": "NM_000024.5:c.48A>N",
            u"g": "NC_000005.9:g.148206442A>N",
            u"n": u"NM_000024.5:n.287A>N"
        },
        {
            u"c": "NM_000024.5:c.49A>N",
            u"g": "NC_000005.9:g.148206443A>N",
            u"n": u"NM_000024.5:n.288A>N"
        },
        {
            u"c": "NM_000024.5:c.50G>N",
            u"g": "NC_000005.9:g.148206444G>N",
            u"n": u"NM_000024.5:n.289G>N"
        },
        {
            u"c": "NM_000024.5:c.51C>N",
            u"g": "NC_000005.9:g.148206445C>N",
            u"n": u"NM_000024.5:n.290C>N"
        },

    # ins, del, delins, dup:
        {
            u"c": "NM_000024.5:c.46_47insNN",
            u"g": "NC_000005.9:g.148206440_148206441insNN",
            u"n": "NM_000024.5:n.285_286insNN"
        },
        {
            u"c": "NM_000024.5:c.45_47delTAG",
            u"g": "NC_000005.9:g.148206439_148206441delTGG",
            u"n": "NM_000024.5:n.284_286delTAG"
        },
        {
            u"c": "NM_000024.5:c.45_47delTAGinsNNNN",
            u"g": "NC_000005.9:g.148206439_148206441delTGGinsNNNN",
            u"n": "NM_000024.5:n.284_286delTAGinsNNNN"
        },
        {
            u"c": "NM_000024.5:c.45_47delTAGinsNNNN",
            u"g": "NC_000005.9:g.148206439_148206441delTGGinsNNNN",
            u"n": "NM_000024.5:n.284_286delTAGinsNNNN"
        },
        {
            u"c": "NM_000024.5:c.46dupA",
            u"g": "NC_000005.9:g.148206440dupG",
            u"n": "NM_000024.5:n.285dupA"
        },

    # IFNA16   │ NM_002173.2 │    6 │  576 │ NC_000009.11 │ -1 │  21216371, 21217310 | 691=2X246=
    # cseq = hdp.fetch_seq("NM_002173.2")
    # gseq = reverse_complement(hdp.fetch_seq("NC_000009.11",21216371,21217310))
    # cseq[685:695] = "AAATTTCAAA"
    # gseq[685:695] = "AAATTTTCAA"
    #                        ^^ @ n.692_693
    # These variants are in and around the 2X substitution
        {
            u"c": "NM_002173.2:c.*110A>N",
            u"g": "NC_000009.11:g.21216625T>N",
            u"n": u"NM_002173.2:n.686A>N"
        },
        {
            u"c": "NM_002173.2:c.*111A>N",
            u"g": "NC_000009.11:g.21216624T>N",
            u"n": u"NM_002173.2:n.687A>N"
        },
        {
            u"c": "NM_002173.2:c.*112A>N",
            u"g": "NC_000009.11:g.21216623T>N",
            u"n": u"NM_002173.2:n.688A>N"
        },
        {
            u"c": "NM_002173.2:c.*113T>N",
            u"g": "NC_000009.11:g.21216622A>N",
            u"n": u"NM_002173.2:n.689T>N"
        },
        {
            u"c": "NM_002173.2:c.*114T>N",
            u"g": "NC_000009.11:g.21216621A>N",
            u"n": u"NM_002173.2:n.690T>N"
        },
        {
            u"c": "NM_002173.2:c.*115T>N",
            u"g": "NC_000009.11:g.21216620A>N",
            u"n": u"NM_002173.2:n.691T>N"
        },
        {
            u"c": "NM_002173.2:c.*116C>N",
            u"g": "NC_000009.11:g.21216619A>N",
            u"n": u"NM_002173.2:n.692C>N"
        },    # ref repl
        {
            u"c": "NM_002173.2:c.*117A>N",
            u"g": "NC_000009.11:g.21216618G>N",
            u"n": u"NM_002173.2:n.693A>N"
        },    # ref repl
        {
            u"c": "NM_002173.2:c.*118A>N",
            u"g": "NC_000009.11:g.21216617T>N",
            u"n": u"NM_002173.2:n.694A>N"
        },
        {
            u"c": "NM_002173.2:c.*119A>N",
            u"g": "NC_000009.11:g.21216616T>N",
            u"n": u"NM_002173.2:n.695A>N"
        },

    # ins, del, delins, dup:
        {
            u"c": "NM_002173.2:c.*115_*117insNN",
            u"g": "NC_000009.11:g.21216618_21216620insNN",
            u"n": "NM_002173.2:n.691_693insNN"
        },
        {
            u"c": "NM_002173.2:c.*114_*117delTTCA",
            u"g": "NC_000009.11:g.21216618_21216621delGAAA",
            u"n": "NM_002173.2:n.690_693delTTCA"
        },
        {
            u"c": "NM_002173.2:c.*115_*117delTCAinsNN",
            u"g": "NC_000009.11:g.21216618_21216620delGAAinsNN",
            u"n": "NM_002173.2:n.691_693delTCAinsNN"
        },
        {
            u"c": "NM_002173.2:c.*115_*117delTCAinsNN",
            u"g": "NC_000009.11:g.21216618_21216620delGAAinsNN",
            u"n": "NM_002173.2:n.691_693delTCAinsNN"
        },
        {
            u"c": "NM_002173.2:c.*115_*117dupTCA",
            u"g": "NC_000009.11:g.21216618_21216620dupGAA",
            u"n": "NM_002173.2:n.691_693dupTCA"
        },

    # genomic variation within an indel discrepancy
    # NM_032790.3  c        108 >
    # NM_032790.3  n        301 > CGGGGAGCCCCCGGGGGCC------CCGCCACCGCCGCCGT
    #                             |||||||||||||||||||------||||||||||||||||
    # NC_000012.11 g  122064755 > CGGGGAGCCCCCGGGGGCCCCGCCACCGCCACCGCCGCCGT
    #                                              **********
        {
            u"c": "NM_032790.3:c.125_128delCCCC",
            u"g": "NC_000012.11:g.122064772_122064781delCCCCGCCACC",
            u"n": "NM_032790.3:n.318_321delCCCC"
        },
    ]

    @classmethod
    def setUpClass(cls):
        def _parse_rec(rec):
            rec["pv"] = {x: cls.hp.parse_hgvs_variant(rec[x]) for x in "cgn"}
            return rec

        cls.hdp = vvhgvs.dataproviders.uta.connect(mode=mode_txt, cache=CACHE)
        cls.am = vvhgvs.assemblymapper.AssemblyMapper(
            cls.hdp, replace_reference=True, assembly_name="GRCh37", alt_aln_method="genbuild")
        cls.hp = vvhgvs.parser.Parser()
        cls.tests = [_parse_rec(rec) for rec in cls.test_cases]

    def test_replace_reference_sequence(self):
        """AssemblyMapper: Replace invalid reference sequence"""

        for rec in self.tests:
            for x in "cgn":
                pv = rec["pv"][x]
                if pv.posedit.edit.ref:
                    # replace ref with junk
                    pv.posedit.edit.ref = "NNNNNN"
                self.am._replace_reference(pv)
                self.assertEqual(rec[x], pv.format(conf={'max_ref_length': None}))


@pytest.mark.quick
class Test_AssemblyMapper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        hdp = vvhgvs.dataproviders.uta.connect(mode=mode_txt, cache=CACHE)
        cls.hp = vvhgvs.parser.Parser()
        cls.am = vvhgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method="genebuild")
        cls.am38 = vvhgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh38", alt_aln_method="genebuild")

    def _test_mapping(self, hgvs_set,gen_ver='37'):
        """given list of variant strings, test all valid combinations of
        g<->n<->c<->p mappings

        """
        parsed_variants = [(hv, self.hp.parse_hgvs_variant(hv)) for hv in hgvs_set]
        hgvs = {v.type: hv for hv, v in parsed_variants}
        pvs = {v.type: v for hv, v in parsed_variants}
        assembbly_maper = self.am
        if gen_ver == '38':
            assembbly_maper = self.am38
        if "g" in pvs and "c" in pvs:
            self.assertEqual(hgvs["g"], str(assembbly_maper.c_to_g(pvs["c"])))
            self.assertEqual(hgvs["c"], str(assembbly_maper.g_to_c(pvs["g"], pvs["c"].ac)))
        if "g" in pvs and "n" in pvs:
            self.assertEqual(hgvs["g"], str(assembbly_maper.n_to_g(pvs["n"])))
            self.assertEqual(hgvs["n"], str(assembbly_maper.g_to_n(pvs["g"], pvs["n"].ac)))
        if "c" in pvs and "n" in pvs:
            self.assertEqual(hgvs["n"], str(assembbly_maper.c_to_n(pvs["c"])))
            self.assertEqual(hgvs["c"], str(assembbly_maper.n_to_c(pvs["n"])))
        if "c" in pvs and "p" in pvs:
            self.assertEqual(hgvs["p"], str(assembbly_maper.c_to_p(pvs["c"])))
            self.assertEqual(hgvs["p"], str(assembbly_maper.t_to_p(pvs["c"])))

    def test_SNV(self):
        """AssemblyMapper: smoketest with SNVs"""
        hgvs_set = [
            #"NC_000007.13:g.36561662C>T", "NM_001637.3:c.1582G>A", "NM_001637.3:n.1983G>A", "NP_001628.1:p.(Gly528Arg)"
            "NC_000007.14:g.36522056C>T", "ENST00000617537.5:c.1582G>A","ENSP00000483783.1:p.(Gly528Arg)"

        ]
        self._test_mapping(hgvs_set,'38')

    def test_intronic(self):
        """AssemblyMapper: smoketest with intronic SNVs"""
        hgvs_set = [
            #"NC_000010.10:g.89711873A>C", "NM_000314.4:c.493-2A>C", "NM_000314.4:n.1524-2A>C", "NP_000305.3:p.?"
            #"NC_000010.11:g.87952116A>C",'ENST00000371953.7:c.493-2A>C',"ENSP00000361021.3:p.?"
            "NC_000010.10:g.89711873A>C",'ENST00000371953.3:c.493-2A>C',"ENSP00000361021.3:p.?"
        ]

        self._test_mapping(hgvs_set)

    def test_t_to_p(self):
        #assert "non-coding" == str(self.am.t_to_p(self.hp.parse("NR_027676.1:n.3980del")))
        assert "ENSP00000369497.3:p.(Lys2597=)" == str(self.am.t_to_p(self.hp.parse("ENST00000380152.7:c.7791A>G")))



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
