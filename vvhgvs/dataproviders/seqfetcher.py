# -*- coding: utf-8 -*-
"""provides sequencing fetching from NCBI and Ensembl

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import os
import re
from threading import Lock

import bioutils.seqfetcher

from ..exceptions import HGVSDataNotAvailableError

_logger = logging.getLogger(__name__)


class SeqFetcher(object):
    """This class is intended primarily as a mixin for HGVS data providers
    that doen't otherwise have access to sequence data.  It uses the
    fetch_seq() function in this module to fetch sequences from
    several sources; see that function for details.

    >> sf = SeqFetcher()

    >> sf.fetch_seq('NP_056374.2',0,10)
    'MESRETLSSS'

    """
    
    def __init__(self, check_same_thread=False):
        # If HGVS_SEQREPO_DIR is defined, we use seqrepo for *all* sequences
        # Otherwise, we fall back to remote sequence fetching
        seqrepo_dir = os.environ.get("HGVS_SEQREPO_DIR")
        self.lock = Lock()

        if seqrepo_dir:
            from biocommons.seqrepo import SeqRepo
            sr = SeqRepo(seqrepo_dir, check_same_thread)

            def _fetch_seq_seqrepo(ac, start_i=None, end_i=None):
                return sr.fetch(ac, start_i, end_i)

            self.check_same_thread = check_same_thread
            self.fetcher = _fetch_seq_seqrepo
            self.source = "SeqRepo ({})".format(seqrepo_dir)
        else:
            quit("""
            V.V. usage can be quite heavy, variant validators "test_configuration.py" asserts that
            we should at least explicitly chose the location, therefore, for vvhgvs, disable silent
            public fallback, explicitly set a external seqrepo location if remote data is needed.
            """)
            self.fetcher = bioutils.seqfetcher.fetch_seq
            self.source = "bioutils.seqfetcher"
        _logger.info("Fetching sequences with " + self.source)

    def fetch_seq(self, ac, start_i=None, end_i=None):
        if self.check_same_thread is False:
            self.lock.acquire()
        try:
            fetched_seq = self.fetcher(ac, start_i, end_i)
        except Exception as ex:
            try:
                self.lock.release()
            except AttributeError:
                pass
            except RuntimeError:
                pass
            raise HGVSDataNotAvailableError("Failed to fetch {ac} from {self.source} ({ex})".format(
                ac=ac, ex=ex, self=self))
        else:
            try:
                self.lock.release()
            except AttributeError:
                pass
            except RuntimeError:
                pass
            return fetched_seq


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
