# -*- coding: utf-8 -*-
"""implements an vvhgvs data provider interface using UTA
(https://github.com/biocommons/uta)

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import contextlib
import inspect
import logging
import os
import re

import psycopg2
import psycopg2.extras
import psycopg2.pool

from bioutils.assemblies import make_ac_name_map
from bioutils.digests import seq_md5
from six.moves.urllib import parse as urlparse

import vvhgvs
from ..dataproviders.interface import Interface
from ..exceptions import HGVSError, HGVSDataNotAvailableError
from .seqfetcher import SeqFetcher
import six

_logger = logging.getLogger(__name__)


def _stage_from_version(version):
    """return "prd", "stg", or "dev" for the given version string.  A value is always returned"""
    if version:
        m = re.match(r"^(?P<xyz>\d+\.\d+\.\d+)(?P<extra>.*)", version)
        if m:
            return "stg" if m.group("extra") else "prd"
    return "dev"


def _get_uta_db_url():
    """returns UTA DB URL based on environment variables and code version

    * if UTA_DB_URL is set, use that
    * Otherwise, if _UTA_URL_KEY is set, use that as the name of a
      config file entry and use the corresponding URL
    * Otherwise, 

    """

    if "UTA_DB_URL" in os.environ:
        return os.environ["UTA_DB_URL"]

    if "_UTA_URL_KEY" in os.environ:
        url_key = os.environ["_UTA_URL_KEY"]
    else:
        quit("""
        V.V. usage can be quite heavy, variant validators "test_configuration.py" asserts that we 
        should at least explicitly chose the location, therefore, for vvhgvs, disable silent public 
        fallback, explicitly set a external url key if remote data is needed.
        """)
        sdlc = _stage_from_version(vvhgvs.__version__)
        url_key = "public_{sdlc}".format(sdlc=sdlc)
    return vvhgvs.global_config['uta'][url_key]


def connect(db_url=None, pooling=vvhgvs.global_config.uta.pooling, application_name=None, mode=None, cache=None):
    """Connect to a UTA database instance and return a UTA interface instance.

    :param db_url: URL for database connection
    :type db_url: string
    :param pooling: whether to use connection pooling (postgresql only)
    :type pooling: bool
    :param application_name: log application name in connection (useful for debugging; PostgreSQL only)
    :type application_name: str

    When called with an explicit db_url argument, that db_url is used for connecting.

    When called without an explicit argument, the function default is
    determined by the environment variable UTA_DB_URL if it exists, or
    vvhgvs.datainterface.uta.public_db_url otherwise.

    >>> hdp = connect()
    >>> hdp.schema_version()
    '1.1'

    The format of the db_url is driver://user:pass@host/database/schema (the same
    as that used by SQLAlchemy).  Examples:

    A remote public postgresql database:
        postgresql://anonymous:anonymous@uta.biocommons.org/uta/uta_20170707'

    A local postgresql database:
        postgresql://localhost/uta_dev/uta_20170707

    For postgresql db_urls, pooling=True causes connect to use a
    psycopg2.pool.ThreadedConnectionPool.
    """

    _logger.debug('connecting to ' + str(db_url) + '...')

    if db_url is None:
        db_url = _get_uta_db_url()

    url = _parse_url(db_url)
    if url.scheme == 'postgresql':
        conn = UTA_postgresql(url=url, pooling=pooling, application_name=application_name, mode=mode, cache=cache)
    else:
        # fell through connection scheme cases
        raise RuntimeError("{url.scheme} in {url} is not currently supported".format(url=url))
    _logger.info('connected to ' + str(db_url) + '...')
    return conn


class UTABase(Interface):
    required_version = "0.9"
    version_type = 'matview_version'
    # for the id quires we need at least uta 1.1 and vvta ver 0.6+ 
    # but no way to specify simply, without breaking whole set on lower ver num
    _queries = {
        "acs_for_protein_md5":"select ac from seq_anno where seq_id=%s",
        "gene_info":"select * from gene where hgnc=%s",
        "gene_info_by_id":"select * from gene where hgnc_id=%s",
        # to get alias (or prev symbol from alias column) use LIKE for
        # now.  'SIMILAR TO' tests slower, and like regex it has
        # safety issues "= ANY (string_to_array(aliases,','))" is the
        # same, if a little slower. If we split first and store as
        # array is 1/4 quicker, and does not have the problem of
        # having to repeat the input
        'gene_info_by_alias_symbol':'''
            SELECT *
            FROM gene
            WHERE
                aliases = %s
                OR aliases LIKE %s || ',%%'
                OR aliases LIKE '%%,' || %s || ',%%'
                OR aliases LIKE '%%,' || %s
            ''',
    # TODO: reconcile tx_exons query and build_tx_cigar
    # built_tx_cigar says it expects exons in transcript order, but this is genomic order.
        "tx_exons":"""
            select tx_ac, alt_ac,alt_aln_method,alt_strand,ord,tx_start_i,tx_end_i,alt_start_i,alt_end_i,cigar
            from tx_exon_aln_mv where tx_ac=%s and alt_ac=%s and alt_aln_method=%s 
            order by alt_start_i
            """,
    # This query should replace tx_exons in all new code, it is pre checked and should be a faster lookup vs
    # the list of results from tx_exons, at least for things like strand, the arrays are transcript order.
    # This contains extra details not used in the uta so that caching will prevent repeated lookups when 
    # this data is used by variant validator.
        "agg_exon_aln":"""
            SELECT 
            alt_strand,mapped_start,not_quite_cigar,mapped_end,
            cds_start_i, cds_end_i,
            transcript_exon_start_end,mapped_exon_start_end
            FROM full_tx_aln_w_nq_cigar_mv
            WHERE tx_ac=%s and alt_ac=%s and alt_aln_method=%s
            """,

        "tx_for_gene":"""
            select hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method
            from current_valid_mapped_transcript_per_gene_mv where hgnc=%s
            """,
        "tx_for_gene_id":"""
            select hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method
            from current_valid_mapped_transcript_per_gene_mv where hgnc_id=%s
            """,
        # note that this gives inclusive results for input coordinates
        # the </<= implicitly does 0 based (db) to 1 based (query) conversion
        "tx_for_region":"""
            select tx_ac,alt_ac,alt_strand,alt_aln_method,start_i,end_i
            from current_valid_mapped_transcript_spans_mv 
            where alt_ac=%s and alt_aln_method=%s and start_i < %s and %s <= end_i
            """,
        "tx_limits":"""
            SELECT ac, cds_start_i, cds_end_i, length, hgnc
            FROM transcript_lengths_mv
            WHERE ac=%s
            """,
        # compat query for old tx_identity_info will work with numeric indexing
        "tx_identity_info":"""
            SELECT ac as tx_ac, NULL AS alt_ac, NULL AS alt_aln_method, cds_start_i, cds_end_i, ARRAY[length] AS lengths, hgnc
            FROM transcript_lengths_mv
            WHERE ac=%s
            """,
        "tx_info":"""
            select hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method
            from all_mapped_transcript_mv
            where tx_ac=%s and alt_ac=%s and alt_aln_method=%s
            """,
        "tx_mapping_options": """
            select distinct tx_ac,alt_ac,alt_aln_method 
            from tx_exon_aln_mv where tx_ac=%s and cigar is not NULL
            """,
        "tx_seq":"select seq from seq S join seq_anno SA on S.seq_id=SA.seq_id where ac=%s",
        "tx_seq_anno":"select len,seq_id,descr from seq_anno join seq using(seq_id) where ac=%s",
        "tx_similar":"select * from tx_similarity_v where tx_ac1 = %s",
        "tx_to_pro":"select * from associated_accessions where tx_ac = %s order by pro_ac desc",
    }

    def __init__(self, url, mode=None, cache=None):
        self.url = url
        self.seqfetcher = SeqFetcher()
        if mode != 'run':
            self._connect()
        super(UTABase, self).__init__(mode, cache)

    def __str__(self):
        return ("{n} <data_version:{dv}; schema_version:{sv}; application_name={self.application_name};"
                " url={self.url}; sequences-from={sf}>").format(
                    n=type(self).__name__,
                    self=self,
                    dv=self.data_version(),
                    sv=self.schema_version(),
                    sf=os.environ.get("HGVS_SEQREPO_DIR", "seqfetcher"))

    def _fetchone(self, sql, *args):
        with self._get_cursor() as cur:
            cur.execute(sql, *args)
            return cur.fetchone()

    def _fetchall(self, sql, *args):
        with self._get_cursor() as cur:
            cur.execute(sql, *args)
            return cur.fetchall()

    ############################################################################
    # Queries

    def data_version(self):
        return self.url.schema

    def schema_version(self):
        return self._fetchone(f"select * from meta where key = '{self.version_type}'")['value']

    def get_seq(self, ac, start_i=None, end_i=None):
        return self.seqfetcher.fetch_seq(ac, start_i, end_i)

    def get_acs_for_protein_seq(self, seq):
        """
        returns a list of protein accessions for a given sequence.  The
        list is guaranteed to contain at least one element with the
        MD5-based accession (MD5_01234abc...def56789) at the end of the
        list.
        """
        md5 = seq_md5(seq)
        return [r['ac'] for r in self._fetchall(self._queries['acs_for_protein_md5'], [md5])] + ['MD5_' + md5]

    def get_gene_info(self, gene):
        """
        returns basic information about the gene.

        :param gene: HGNC gene name
        :type gene: str

        # database results
        hgnc    | ATM
        maploc  | 11q22-q23
        descr   | ataxia telangiectasia mutated
        summary | The protein encoded by this gene belongs to the PI3/PI4-kinase family. This...
        aliases | AT1,ATA,ATC,ATD,ATE,ATDC,TEL1,TELO1
        added   | 2014-02-04 21:39:32.57125

        """
        return self._fetchone(self._queries['gene_info'], [gene])
    #same as above but by id not symbol
    def get_gene_info_by_id(self, gene_id):
        return self._fetchone(self._queries['gene_info_by_id'], [gene_id])
    def get_gene_info_by_alias(self, gene_alias):
        return self._fetchall(
                self._queries['gene_info_by_alias_symbol'],
                [gene_alias,gene_alias,gene_alias,gene_alias]
                )

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        """
        return transcript exon info for supplied accession (tx_ac, alt_ac, alt_aln_method), or None if not found

        :param tx_ac: transcript accession with version (e.g., 'NM_000051.3')
        :type tx_ac: str

        :param alt_ac: specific genomic sequence (e.g., NC_000011.4)
        :type alt_ac: str

        :param alt_aln_method: sequence alignment method (e.g., splign, blat)
        :type alt_aln_method: str

        # tx_exons = db.get_tx_exons('NM_199425.2', 'NC_000020.10', 'splign')
        # len(tx_exons)
        3

        tx_exons have the following attributes::

            {
                'tes_exon_set_id' : 98390
                'aes_exon_set_id' : 298679
                'tx_ac'           : 'NM_199425.2'
                'alt_ac'          : 'NC_000020.10'
                'alt_strand'      : -1
                'alt_aln_method'  : 'splign'
                'ord'             : 2
                'tx_exon_id'      : 936834
                'alt_exon_id'     : 2999028
                'tx_start_i'      : 786
                'tx_end_i'        : 1196
                'alt_start_i'     : 25059178
                'alt_end_i'       : 25059588
                'cigar'           : '410='
            }

        For example:

        # tx_exons[0]['tx_ac']
        'NM_199425.2'

        """
        rows = self._fetchall(self._queries['tx_exons'], [tx_ac, alt_ac, alt_aln_method])
        if len(rows) == 0:
            raise HGVSDataNotAvailableError(
                "No tx_exons for (tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})".format(
                    tx_ac=tx_ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method))

        # TODO: Check that end == transcript sequence length (but length N/A in current hdp)
        ex0 = 0 if (rows[0]["alt_strand"] == 1) else -1
        if rows[ex0]["tx_start_i"] != 0:
            raise HGVSDataNotAvailableError("Alignment is incomplete; cannot use transcript for mapping"
                                            "(tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})".format(
                                                tx_ac=tx_ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method))
        return rows

    def get_agg_exon_aln(self, tx_ac, alt_ac, alt_aln_method):
        """
        return transcript alignment details for supplied (tx_ac, alt_ac, alt_aln_method), or None if not found
        pre-filtered for start =0 and contiguousness 
        return order = strand,cigar start offset, not quite cigar format alignment,mapped end pos,
        cds start,cds end,exon pos sets, exon mapping pos sets
        """
        return self._fetchone(self._queries['agg_exon_aln'], [tx_ac, alt_ac, alt_aln_method])


    def get_tx_for_gene(self, gene):
        """
        return transcript info records for supplied gene, in order of decreasing length

        :param gene: HGNC gene name
        :type gene: str
        """
        return self._fetchall(self._queries['tx_for_gene'], [gene])
    #same as above but by id not symbol
    def get_tx_for_gene_id(self, gene_id):
        return self._fetchall(self._queries['tx_for_gene_id'], [gene_id])

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        """
        return transcripts that overlap given region, inclusive

        :param str alt_ac: reference sequence (e.g., NC_000007.13)
        :param str alt_aln_method: alignment method (e.g., splign)
        :param int start_i: 5' bound of region (flipped on start_i > end_i)
        :param int end_i: 3' bound of region (flipped on start_i > end_i)
        """
        if start_i > end_i:
            start = end_i
            end_i = start_i
            start_i = start
        return self._fetchall(self._queries['tx_for_region'], [alt_ac, alt_aln_method, end_i, start_i])

    def get_tx_identity_info(self, tx_ac):
        """returns features associated with a single transcript.

        :param tx_ac: transcript accession with version (e.g., 'NM_199425.2')
        :type tx_ac: str

        # database output
        -[ RECORD 1 ]--+-------------
        tx_ac          | NM_199425.2
        alt_ac         | NM_199425.2
        alt_aln_method | transcript
        cds_start_i    | 283
        cds_end_i      | 1003
        lengths        | {707,79,410}
        hgnc           | VSX1

        """
        rows = self._fetchall(self._queries['tx_identity_info'], [tx_ac])
        if len(rows) == 0:
            raise HGVSDataNotAvailableError("No transcript definition for (tx_ac={tx_ac})".format(tx_ac=tx_ac))
        return rows[0]
    def get_tx_limits(self, tx_ac):
        """returns gene symbol and non alignment derived transcript features associated with a single transcript.
        same as get_tx_identity_info but,
        does not generate dummy values for alt_ac and method 
        does not wrap transcript length as lengths in array
        ac,cds_start_i,cds_end_i,length,hgnc
        """
        rows = self._fetchall(self._queries['tx_limits'], [tx_ac])
        if len(rows) == 0:
            raise HGVSDataNotAvailableError("No transcript definition for (tx_ac={tx_ac})".format(tx_ac=tx_ac))
        return rows[0]
       
    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        """return a single transcript info for supplied accession (tx_ac, alt_ac, alt_aln_method), or None if not found

        :param tx_ac: transcript accession with version (e.g., 'NM_000051.3')
        :type tx_ac: str

        :param alt_ac: specific genomic sequence (e.g., NC_000011.4)
        :type alt_ac: str

        :param alt_aln_method: sequence alignment method (e.g., splign, blat)
        :type alt_aln_method: str

        # database output
        -[ RECORD 1 ]--+------------
        hgnc           | ATM
        cds_start_i    | 385
        cds_end_i      | 9556
        tx_ac          | NM_000051.3
        alt_ac         | AC_000143.1
        alt_aln_method | splign

        """
        rows = self._fetchall(self._queries['tx_info'], [tx_ac, alt_ac, alt_aln_method])
        if len(rows) == 0:
            raise HGVSDataNotAvailableError(
                "No tx_info for (tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})".format(
                    tx_ac=tx_ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method))
        elif len(rows) == 1:
            return rows[0]
        else:
            raise HGVSError("Multiple ({n}) replies for tx_info(tx_ac="
                            "{tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})".format(
                                n=len(rows), tx_ac=tx_ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method))

    def get_tx_mapping_options(self, tx_ac):
        """Return all transcript alignment sets for a given transcript
        accession (tx_ac); returns empty list if transcript does not
        exist.  Use this method to discovery possible mapping options
        supported in the database

        :param tx_ac: transcript accession with version (e.g., 'NM_000051.3')
        :type tx_ac: str

        # database output
        -[ RECORD 1 ]--+------------
        hgnc           | ATM
        cds_start_i    | 385
        cds_end_i      | 9556
        tx_ac          | NM_000051.3
        alt_ac         | AC_000143.1
        alt_aln_method | splign
        -[ RECORD 2 ]--+------------
        hgnc           | ATM
        cds_start_i    | 385
        cds_end_i      | 9556
        tx_ac          | NM_000051.3
        alt_ac         | NC_000011.9
        alt_aln_method | blat

        """
        rows = self._fetchall(self._queries['tx_mapping_options'], [tx_ac])
        return rows

    def get_similar_transcripts(self, tx_ac):
        """Return a list of transcripts that are similar to the given
        transcript, with relevant similarity criteria.

        >> sim_tx = hdp.get_similar_transcripts('NM_001285829.1')
        >> dict(sim_tx[0])
        { 'cds_eq': False,
        'cds_es_fp_eq': False,
        'es_fp_eq': True,
        'tx_ac1': 'NM_001285829.1',
        'tx_ac2': 'ENST00000498907' }

        where:

        * cds_eq means that the CDS sequences are identical
        * es_fp_eq means that the full exon structures are identical
          (i.e., incl. UTR)
        * cds_es_fp_eq means that the cds-clipped portions of the exon
          structures are identical (i.e., ecluding. UTR)
        * Hint: "es" = "exon set", "fp" = "fingerprint", "eq" = "equal"

        "exon structure" refers to the start and end coordinates on a
        specified reference sequence. Thus, having the same exon
        structure means that the transcripts are defined on the same
        reference sequence and have the same exon spans on that
        sequence.

        """

        rows = self._fetchall(self._queries['tx_similar'], [tx_ac])
        return rows

    def get_pro_ac_for_tx_ac(self, tx_ac):
        """Return the (single) associated protein accession for a given transcript
        accession, or None if not found."""

        rows = self._fetchall(self._queries['tx_to_pro'], [tx_ac])
        try:
            return rows[0]['pro_ac']
        except IndexError:
            return None

    def get_tx_seq_anno(self,tx_ac):
        """Return the length, checksum type seq_id, and description for a
        transcript, in that order.
        """
        return self._fetchone(self._queries['tx_seq_anno'], [tx_ac])

    def get_assembly_map(self, assembly_name):
        """return a list of accessions for the specified assembly name (e.g., GRCh38.p5)

        """
        return make_ac_name_map(assembly_name)


class UTA_postgresql(UTABase):
    def __init__(self, url, pooling=vvhgvs.global_config.uta.pooling, application_name=None, mode=None, cache=None):
        if url.schema is None:
            raise Exception("No schema name provided in {url}".format(url=url))
        self.application_name = application_name
        self.pooling = pooling
        self._conn = None
        super(UTA_postgresql, self).__init__(url, mode, cache)

    def __del__(self):
        self.close()

    def close(self):
        if self.pooling:
            self._pool.closeall()
        else:
            if self._conn is not None:
                self._conn.close()

    def _connect(self):
        if self.application_name is None:
            st = inspect.stack()
            self.application_name = os.path.basename(st[-1][1])
        # fix for allowing '/' in host names (when actually UNIX socket files)
        host_or_socketfile = self.url.hostname.replace('%2F','/')
        conn_args = dict(
            host=host_or_socketfile,
            port=self.url.port,
            database=self.url.database,
            user=self.url.username,
            password=self.url.password,
            application_name=self.application_name + "/" + vvhgvs.__version__,
        )
        if self.pooling:
            _logger.info("Using UTA ThreadedConnectionPool")
            self._pool = psycopg2.pool.ThreadedConnectionPool(vvhgvs.global_config.uta.pool_min,
                                                              vvhgvs.global_config.uta.pool_max, **conn_args)
        else:
            self._conn = psycopg2.connect(**conn_args)
            self._conn.autocommit = True

        self._ensure_schema_exists()

    def _ensure_schema_exists(self):
        # N.B. On AWS RDS, information_schema.schemata always returns zero rows
        r = self._fetchone("select exists(SELECT 1 FROM pg_namespace WHERE nspname = %s)", [self.url.schema])
        if r[0]:
            return
        raise HGVSDataNotAvailableError("specified schema ({}) does not exist (url={})".format(
            self.url.schema, self.url))

    @contextlib.contextmanager
    def _get_cursor(self, n_retries=1):
        """Returns a context manager for obtained from a single or pooled
        connection, and sets the PostgreSQL search_path to the schema
        specified in the connection URL.

        Although *connections* are threadsafe, *cursors* are bound to
        connections and are *not* threadsafe. Do not share cursors
        across threads.

        Use this funciton like this::

            with hdp._get_cursor() as cur:
                # your code

        Do not call this function outside a contextmanager.

        """

        n_tries_rem = n_retries + 1
        while n_tries_rem > 0:
            try:

                conn = self._pool.getconn() if self.pooling else self._conn

                # autocommit=True obviates closing explicitly
                conn.autocommit = True

                cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
                cur.execute("set search_path = {self.url.schema};".format(self=self))

                yield cur

                # contextmanager executes these when context exits
                cur.close()
                if self.pooling:
                    self._pool.putconn(conn)

                break

            except psycopg2.OperationalError:

                _logger.warning("Lost connection to {url}; attempting reconnect".format(url=self.url))
                if self.pooling:
                    self._pool.closeall()
                self._connect()
                _logger.warning("Reconnected to {url}".format(url=self.url))

            n_tries_rem -= 1

        else:

            # N.B. Probably never reached
            raise HGVSError("Permanently lost connection to {url} ({n} retries)".format(url=self.url, n=n_retries))


class ParseResult(urlparse.ParseResult):
    """Subclass of url.ParseResult that adds database and schema methods,
    and provides stringification.

    """

    def __new__(cls, pr):
        return super(ParseResult, cls).__new__(cls, *pr)

    @property
    def database(self):
        path_elems = self.path.split("/")
        return path_elems[1] if len(path_elems) > 1 else None

    @property
    def schema(self):
        path_elems = self.path.split("/")
        return path_elems[2] if len(path_elems) > 2 else None

    def __str__(self):
        return self.geturl()


def _parse_url(db_url):
    """parse database connection urls into components

    UTA database connection URLs follow that of SQLAlchemy, except
    that a schema may be optionally specified after the database. The
    skeleton format is:

       driver://user:pass@host/database/schema

    >>> params = _parse_url("driver://user:pass@host:9876/database/schema")

    >>> params.scheme
    u'driver'

    >>> params.hostname
    u'host'

    >>> params.username
    u'user'

    >>> params.password
    u'pass'

    >>> params.database
    u'database'

    >>> params.schema
    u'schema'

    """
    return ParseResult(urlparse.urlparse(db_url))


if __name__ == "__main__":
    import doctest
    doctest.testmod()

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
