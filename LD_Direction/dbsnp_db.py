# -*- coding: utf-8 -*-
"""
Python SQLAlchemy Interface for dbSNP
"""
import os as _os
from subprocess import check_output as _check_output

from multiprocessing import Pool as _Pool
from multiprocessing import cpu_count as _cpu_count

from numpy import array_split as _array_split

from tqdm import tqdm

import sqlalchemy as sa

from sqlalchemy.ext.declarative import declarative_base as _base
from sqlalchemy.orm import sessionmaker as _sessionmaker


class DB(object):

    """A class for interacting with the database."""

    Base = _base()

    class Row(Base):

        """A simple SQLAlchemy dbSNP interface."""

        __tablename__ = 'dbSNP'

        id     = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
        name   = sa.Column(sa.String, index=True)
        chrom  = sa.Column(sa.String(length = 4), index=True)
        start  = sa.Column(sa.Integer, index=True)
        end    = sa.Column(sa.Integer, index=True)
        strand = sa.Column(sa.String(length=1), index=True)

        sa.Index('chrom_start_index', chrom, start)

        def __repr__(self):
            """Display simple information about the row."""
            return '{name}<{chrom}:{start}-{end}>'.format(
                name=self.name, chrom=self.chrom,
                start=self.start, end=self.end
            )

    def __init__(self, version='149', location='/godot/dbsnp'):
        """Connect to the database."""
        self.location = _os.path.abspath(location)
        self.version = version
        if not _os.path.isdir(self.location):
            raise ValueError('{} is not a directory'.format(location))
        self.file = '{}/dbsnp{}.db'.format(self.location, self.version)
        self.engine = sa.create_engine('sqlite:///{}'.format(self.file))
        self.Session = _sessionmaker(bind=self.engine)

    def get_session(self):
        """Return a session for this database."""
        return self.Session()

    def query(self, *args, **kwargs):
        """Create a pre-initialized query."""
        if not args:
            args = (self.Row,)
        session = self.get_session()
        return session.query(*args, **kwargs)

    def lookup_rsids(self, rsids: list([str])) -> list([Row]):
        """Return either one row or a list of rows by rsID.

        Parameters
        ----------
        rsids : list_of_str or str

        Returns
        -------
        list
            A list of DB.Row objects
        """
        if isinstance(rsids, str):
            rsids = [rsids]
        for i in rsids:
            assert isinstance(i, str)
            assert i.startswith('rs')
        return self.query().filter(self.Row.name.in_(rsids)).all()

    def lookup_location(self, chrom: str, start: int, end: int=None) -> Row:
        """Return a row by location.

        Only does one at a time, optimized for speed.

        Parameters
        ----------
        chrom : str
        start : int
        end : int, optional

        Returns
        -------
        DB.Row
        """
        if not isinstance(chrom, (str, int)):
            raise ValueError('chrom must be int or string')
        chrom = str(chrom)
        start = int(start)
        if end:
            end = int(end)
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        query =  self.query().filter(
            self.Row.chrom == chrom
        ).filter(
            self.Row.start == start
        )
        if end:
            query = query.filter(self.Row.end == end)
        query = query.with_hint(self.Row, 'USE INDEX chrom_start_index')
        return query.first()

    def lookup_locations(self, locs: dict(str=list)) -> list:
        """Return a row by location.

        One query for every chromosome.

        Parameters
        ----------
        loct : dict
            In format: {chrom: list of starts}

        Returns
        -------
        list of DB.Row
        """
        if not isinstance(list(locs.keys())[0], (str, int)):
            raise ValueError('chrom items must be int or string')
        c = {}
        for i in locs.keys():
            o = locs[i]
            i = str(i)
            if not i.startswith('chr'):
                i = 'chr' + i
            s = []
            for j in o:
                s.append(int(j))
            c[i] = s
        locs = c

        results = []
        for chrom in locs:
            query =  self.query().filter(
                self.Row.chrom == chrom
            ).filter(
                self.Row.start.in_(locs[chrom])
            )
            query = query.with_hint(self.Row, 'USE INDEX chrom_start_index')
            results += query.all()
        return results

    def initialize_db(self, f, commit_every=1000000):
        """Initialize the database, must exist already.

        Args:
            f (str): File to read.
            commit_every (int): How many rows to wait before commiting.
        """

        def write_to_db(records, insert):
            """Add extra columns in parallel and write to database.

            Args:
                records (list):  List of dicts of columns from initialize_db()
                insert: Insert object from initialize_db()
            """
            conn = self.engine.connect()
            chunks = _cpu_count()-1

            with _Pool(processes=chunks) as pool:
                results = []
                processing = []
                for chunk in _array_split(records, chunks):
                    processing.append(pool.apply_async(
                        _create_extra_columns,
                        (list(chunk),)
                    ))
                for in_process in processing:
                    results += in_process.get()

            if results:
                conn.execute(insert, results)

        a = input('Create a new database (will destroy existing one)? [y/N] ')
        if not a.upper().startswith('Y'):
            return

        # Create db tables, must be deleted first
        self.Base.metadata.create_all(self.engine)
        rows    = 0
        count   = commit_every
        dbsnp   = self.Row.__table__
        insert  = dbsnp.insert()
        db_len  = int(
            _check_output('wc -l {}'.format(f), shell=True)
            .decode().split(' ')[0]
        )
        records = []
        with open(f) as fin, tqdm(unit='rows', total=int(db_len)) as pbar:
            for line in fin:
                if line.startswith('track'):
                    continue
                chrom, start, end, name, _, strand = line.rstrip().split('\t')
                records.append(
                    {'name': name, 'chrom': chrom, 'start': start,
                     'end': end, 'strand': strand}
                )
                if count:
                    count -= 1
                else:
                    pbar.write('Writing {} records...'.format(commit_every))
                    write_to_db(records, insert)
                    pbar.write('Written')
                    count = commit_every-1
                    records = []
                rows += 1
                pbar.update()
        print('Writing final rows')
        write_to_db(records, insert)
        print('Done, {} rows written'.format(rows))

    @property
    def exists(self):
        """Check if the file exists."""
        return _os.path.isfile(self.file)


def _create_extra_columns(rows):
    """Add length and SNP bool to columns list.

    Args:
        rows (list): List of dicts of columns from initialize_db()

    Returns:
        list: Replacement dictionary list for initialize_db containing new cols
    """
    for row in rows:
        row['start']  = int(row['start'])
        row['end']    = int(row['end'])
        if row['end']-row['start'] > 1:
            rows.remove(row)
    return rows
