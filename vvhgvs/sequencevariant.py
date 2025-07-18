# -*- coding: utf-8 -*-
"""represents simple sequence-based variants

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import attr
import re

import vvhgvs.variantmapper
from vvhgvs.enums import ValidationLevel
from vvhgvs.utils.validation import validate_type_ac_pair


@attr.s(slots=True, repr=False,eq=False)
class SequenceVariant(object):
    """
    represents a basic HGVS variant.  The only requirement is that each
    component can be stringified; for example, passing pos as either a string
    or an vvhgvs.location.CDSInterval (for example) are both intended uses
    """
    ac = attr.ib()
    type = attr.ib()
    posedit = attr.ib()
    rel_ac = attr.ib(default='')

    def format(self, conf=None):
        """Formatting the stringification of sequence variants

        :param conf: a dict comprises formatting options. None is to use global settings.

        See :class:`vvhgvs.config`.
        """
        if self.posedit is not None:
            posedit = self.posedit.format(conf)
        else:
            posedit = "?"
        if self.ac is not None:
            return "{ac}:{type}.{posedit}".format(ac=self.ac, type=self.type, posedit=posedit)
        else:
            return "{type}.{posedit}".format(type=self.type, posedit=posedit)

    __str__ = format

    def __repr__(self):
        return "{0}({1})".format(self.__class__.__name__, ", ".join(
            (a.name + "=" + str(getattr(self, a.name))) for a in self.__attrs_attrs__))

    def __eq__(self,other):
        """Count transcripts with unset relative accessions as equal to any set acc

        This is done so that a simple t1->g->t2 for the same tx ac will asses t1==t2 as True
        """
        if not isinstance(other,SequenceVariant):
            return False
        if self.ac != other.ac or self.type != other.type or self.posedit != other.posedit:
            return False
        #if self.rel_ac and other.rel_ac and self.rel_ac != other.rel_ac:
        #    return False
        return True

    def __ne__(self,other):
        """Count transcripts with unset relative accessions as equal to any set acc

        This is done so that a simple t1->g->t2 for the same tx ac will asses t1==t2 as True
        """
        if not isinstance(other,SequenceVariant):
            return True
        if self.ac != other.ac or self.type != other.type or self.posedit != other.posedit:
            return True
        #if self.rel_ac and other.rel_ac and self.rel_ac != other.rel_ac:
        return False


    def fill_ref(self, hdp):
        hm = vvhgvs.variantmapper.VariantMapper(hdp)
        type = None
        if isinstance(self.posedit, vvhgvs.posedit.PosEdit) and isinstance(self.posedit.edit, vvhgvs.edit.Edit):
            type = self.posedit.edit.type
        if type in ["del", "delins", "identity", "dup", "repeat"] and self.posedit.edit.ref_s is None:
            hm._replace_reference(self)
        if type == "identity" and isinstance(self.posedit.edit, vvhgvs.edit.NARefAlt):
            self.posedit.edit.alt = self.posedit.edit.ref
        return self

    def validate(self):
        (res, msg) = (ValidationLevel.VALID, None)
        if self.ac and self.type:
            (res, msg) = validate_type_ac_pair(self.type, self.ac)
            if res == ValidationLevel.ERROR:
                return (res, msg)
        if self.posedit is None:
            return (res, msg)
        (pe_res, pe_msg) = self.posedit.validate()
        if pe_res == ValidationLevel.VALID:
            return (res, msg)
        return (pe_res, pe_msg)


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
