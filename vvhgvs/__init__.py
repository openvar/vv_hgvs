# -*- coding: utf-8 -*-
"""hgvs is a package to parse, format, and manipulate biological sequence
variants.  See https://github.com/biocommons/hgvs/ for details.

Example use:

>>> import vvhgvs.dataproviders.uta
>>> import vvhgvs.parser
>>> import vvhgvs.variantmapper

# start with these variants as strings
>>> hgvs_g, hgvs_c = "NC_000007.13:g.36561662C>T", "NM_001637.3:c.1582G>A"

# parse the genomic variant into a Python structure
>>> hp = vvhgvs.parser.Parser()
>>> var_g = hp.parse_hgvs_variant(hgvs_g)
>>> var_g
SequenceVariant(ac=NC_000007.13, type=g, posedit=36561662C>T)

# SequenceVariants are composed of structured objects, e.g.,
>>> var_g.posedit.pos.start
SimplePosition(base=36561662, uncertain=False)

# format by stringification
>>> str(var_g)
'NC_000007.13:g.36561662C>T'

# initialize the mapper for GRCh37 with splign-based alignments
>>> hdp = vvhgvs.dataproviders.uta.connect()
>>> am = vvhgvs.assemblymapper.AssemblyMapper(hdp,
...          assembly_name="GRCh37", alt_aln_method="splign",
...          replace_reference=True)

# identify transcripts that overlap this genomic variant
>>> transcripts = am.relevant_transcripts(var_g)
>>> sorted(transcripts)
['NM_001177506.1', 'NM_001177507.1', 'NM_001637.3']

# map genomic variant to one of these transcripts
>>> var_c = am.g_to_c(var_g, "NM_001637.3")
>>> var_c
SequenceVariant(ac=NM_001637.3, type=c, posedit=1582G>A)
>>> str(var_c)
'NM_001637.3:c.1582G>A'

# CDS coordinates use BaseOffsetPosition to support intronic offsets
>>> var_c.posedit.pos.start
BaseOffsetPosition(base=1582, offset=0, datum=Datum.CDS_START, uncertain=False)

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import re
import warnings

try:
    # Use stdlib importlib.metadata when available (Python 3.8+)
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    # For Python < 3.8 use the backport
    from importlib_metadata import version, PackageNotFoundError

from .config import global_config    # noqa (importing symbol)

logger = logging.getLogger(__name__)

_is_released_version = False
try:
    __version__ = version("vvhgvs")
    # TODO: match other valid release tags, like .post\d+
    if re.match(r"^\d+\.\d+\.\d+$", __version__) is not None:
        _is_released_version = True
except PackageNotFoundError:
    warnings.warn("can't get __version__ because vvhgvs package isn't installed", Warning)
    __version__ = None

# Enable DeprecationWarnings for the vvhgvs package
# N.B. The module name is provided as a regexp to the module *path*
warnings.filterwarnings('default', '', DeprecationWarning, '.*\\Wlib\\W.*\\Wvvhgvs\\W.*')

logger.info("vvhgvs " + str(__version__) + "; released: " + str(_is_released_version))

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
