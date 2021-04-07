#!/usr/bin/env python

import vvhgvs
import vvhgvs.dataproviders.uta
import vvhgvs.parser
import vvhgvs.sequencevariant
import vvhgvs.assemblymapper

hdp = vvhgvs.dataproviders.uta.connect()
hp = vvhgvs.parser.Parser()
evm = vvhgvs.assemblymapper.AssemblyMapper(hdp,
                                           replace_reference=True, assembly_name='GRCh37',
                                           alt_aln_method='splign')

#v = hp.parse_hgvs_variant("NM_000059.3:c.7790A>G")
#v = hp.parse_hgvs_variant("NM_000059.3:c.7790_7792delAAG")
v = hp.parse_hgvs_variant("NM_000059.3:c.7790delAAG")
evm.c_to_p(v)

