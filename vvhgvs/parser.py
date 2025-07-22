# -*- coding: utf-8 -*-
"""Provides parser for HGVS strings and HGVS-related conceptual
components, such as intronic-offset coordiates
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import copy
import re
import importlib.resources as resources  # ✅ Replacement for pkg_resources

import bioutils.sequences
import ometa.runtime
import parsley

from vvhgvs.exceptions import HGVSParseError

# The following imports are referenced by fully-qualified name in the
# vvhgvs grammar.
import vvhgvs.enums
import vvhgvs.edit
import vvhgvs.hgvsposition
import vvhgvs.location
import vvhgvs.posedit
import vvhgvs.sequencevariant


def parse_error_hash(self):
    """Define missing ParseError.__hash__()."""
    return hash((self.position, self.formatReason()))


parsley.ParseError.__hash__ = parse_error_hash


class Parser(object):
    """Provides comprehensive parsing of HGVS variant strings into structured Python representations."""

    def __init__(self, grammar_fn=None, expose_all_rules=False):
        if grammar_fn is None:
            # ✅ Replacing pkg_resources.resource_filename
            grammar_fn = resources.files(__package__ + "._data").joinpath("hgvs.pymeta")

        self._grammar_fn = grammar_fn
        with open(self._grammar_fn, "r") as f:
            grammar_text = f.read()

        self._grammar = parsley.makeGrammar(
            grammar_text, {
                "vvhgvs": vvhgvs,
                "bioutils": bioutils,
                "copy": copy
            })
        self._logger = logging.getLogger(__name__)
        self._expose_rule_functions(expose_all_rules)

    def parse(self, v):
        """parse HGVS variant `v`, returning a SequenceVariant"""
        return self.parse_hgvs_variant(v)

    def _expose_rule_functions(self, expose_all_rules=False):
        """add parse functions for public grammar rules"""

        def make_parse_rule_function(rule_name):
            def rule_fxn(s):
                try:
                    return self._grammar(s).__getattr__(rule_name)()
                except ometa.runtime.ParseError as exc:
                    raise HGVSParseError("{s}: char {exc.position}: {reason}".format(
                        s=s, exc=exc, reason=exc.formatReason()))
            rule_fxn.__doc__ = "parse string s using `%s' rule" % rule_name
            return rule_fxn

        exposed_rule_re = re.compile(r"hgvs_(variant|position)|(c|g|m|n|p|r)"
                                     r"_(edit|hgvs_position|interval|pos|posedit|variant)")
        exposed_rules = [m.replace("rule_", "") for m in dir(self._grammar._grammarClass) if m.startswith("rule_")]
        if not expose_all_rules:
            exposed_rules = [rule_name for rule_name in exposed_rules if exposed_rule_re.match(rule_name)]
        for rule_name in exposed_rules:
            att_name = "parse_" + rule_name
            rule_fxn = make_parse_rule_function(rule_name)
            self.__setattr__(att_name, rule_fxn)
        self._logger.debug("Exposed {n} rules ({rules})".format(n=len(exposed_rules), rules=", ".join(exposed_rules)))


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
