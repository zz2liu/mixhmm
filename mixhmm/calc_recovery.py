"""calc snp recoveries of state assignment.

create_SNP
insert_SNP
create_RegionCnv
insert_RegionCnv

create_SnpCnv
insert_SnpCnv_

create_SnpRecovery
recovs = recovery_from_SnpCnv
recovery_to_SnpRecovery(recovs)
"""
from __future__ import division
from py25 import partial, defaultdict, groupby, itemgetter, set_trace, warn, izip
import csv, sqlite3
from numpy import array, zeros, nan
from py_cna.cnv import Cnv
from py_cna.penncnv import parse_rawcnv as parse_penncnv
from py_cna.genoCN import states as genocna_states, copynumbers, imbalances

def create_SNP(cur, tname='SNP'):
    """create the SNP table
    """
    cur.execute('''create table %s
        (Chr, Loc integer,
        primary key (Chr, Loc))
        '''% tname)

def insert_SNP(cur, fname, cols=[1,2], tname='SNP'):
    """insert rows into SNP table.

    cols: colidxs of (Chr, Loc)
    """
    INSERT = '''insert into %s values (?, ?)
        ''' % tname
    rows = csv.reader(open(fname), delimiter='\t')
    rows.next() #ignore header
    for row in rows:
        row_ = [row[i] for i in cols]
        row_[1] = long(row_[1]) #convert the Loc
        try:
            cur.execute(INSERT, tuple(row_))
        except sqlite3.IntegrityError, e:
            warn('%s: %s' % (row, e))

def create_RegionCnv(cur, tname='RegionCnv'):
    """ create the table RegionCnv 
    """
    # debug: there must be a comma at the end
    create = '''create table %s
        (Chip, Chr, Start integer, End integer, State, Cn real, Im real,
        primary key (Chip, Chr, Start, End))
        ''' % tname
    print create
    cur.execute(create)

def each_row_from_penncnv(cnv_file, parse_penncnv=parse_penncnv):
    """yield each row of (chr, start, end, state, cn, im) from penncnv.
    """
    states = ['ignore',
            '0n',
            '1n',
            'FM',
            'FF',
            '3n',
            '4n'
            ]
    for row in parse_penncnv(cnv_file):
        chr, start, end, state, cn = map(row.__getitem__,
                ['chr', 'start', 'end', 'state', 'cn'])
        state = states[state]
        im = (0.5 if cn==2 else 0)
        yield chr, start, end, state, cn, im

def each_row_from_mixhmm(cnv_file):
    rows = csv.reader(cnv_file, delimiter='\t')
    rows.next()
    for row in rows:
        yield row

def each_row_from_genoCNA(cnv_file):
    rows = csv.reader(cnv_file, delimiter='\t')
    rows.next()
    for row in rows:
        chr, start, end, state = row[:4]
        si = int(state)-1
        yield (chr, long(start), long(end), genocna_states[si],
                copynumbers[si], imbalances[si])

def insert_RegionCnv(cur, chips, cnv_fnames,
        cols=[0, 1,2,-3, -2, -1], tname='RegionCnv',
        each_row=None):
    """
    cnv_fnames: each cnv tsv file with the following cols ...
    cols: the col idx of (chr, start, end, state, cn, im) in each file.
    each_row: a function(cnv_file) -> each row
    """
    INSERT = '''insert into %s
        values (?, ?, ?, ?, ?, ?, ?)
        ''' % tname
    cnv_files = map(open, cnv_fnames)
    for chip, cnv_file in zip(chips, cnv_files):
        print chip
        for row in each_row(cnv_file):
            chr, start, end, state, cn, im = [row[i] for i in cols]
            start, end = long(start), long(end)
            cn, im = float(cn), float(im)
            cur.execute(INSERT, (chip, chr, start, end, state, cn, im))

insert_RegionCnv_from_mixhmm =  partial(insert_RegionCnv,
    each_row=each_row_from_mixhmm)
insert_RegionCnv_from_penncnv = partial(insert_RegionCnv,
    each_row=each_row_from_penncnv)
insert_RegionCnv_from_genoCNA = partial(insert_RegionCnv,
    each_row=each_row_from_genoCNA)

def insert_RegionCnv_old(cur, chips, cnv_fnames,
        cols=[0, 1,2,-3, -2, -1], tname='RegionCnv'):
    """
    cnv_fnames: each cnv tsv file with the following cols ...
    cols: the col idx of (chr, start, end, state, cn, im) in each file.
    """
    INSERT = '''insert into %s
        values (?, ?, ?, ?, ?, ?, ?)
        ''' % tname
    cnv_files = map(open, cnv_fnames)
    for chip, cnv_file in zip(chips, cnv_files):
        print chip
        rows = csv.reader(cnv_file, delimiter='\t')
        rows.next()
        for row in rows:
            chr, start, end, state, cn, im = [row[i] for i in cols]
            start, end = long(start), long(end)
            cn, im = float(cn), float(im)
            cur.execute(INSERT, (chip, chr, start, end, state, cn, im))

def create_SnpCnv(cur, tname='SnpCnv'):
    create = '''create table %s
        (Chip, Chr, Loc integer, State, Cn real, Im real,
        primary key (Chip, Chr, Loc))
        ''' % tname
    cur.execute(create)

def get_chips(cur, chiplike='%', tname='RegionCnv', complementary=False):
    """return a dictintive 'chip' list"""
    if complementary:
        select = """select distinct chip
            from %s
            where chip not like %r
            """ % (tname, chiplike)
    else:
        select = """select distinct chip
            from %s
            where chip like %r
            """ % (tname, chiplike)
    print select
    cur.execute(select)
    return [e[0].encode() for e in cur]

def insert_SnpCnv_(cur, tSnp='SNP', tRegionCnv='RegionCnv', tname='SnpCnv', chiplike='%'):
    """using python to expand region cnv to snp cnvs

    chiplike: select a subset of chips (where chip like {chiplike})to expand.
    """
    insert = '''insert into %s
        values (?, ?, ?, ?, ?, ?)''' % tname

    def _select_snp_rows():
        select = """select Chr, Loc from %s
        order by Chr, Loc
        """ % tSnp
        cur.execute(select)
        return [(row[0].encode(), row[1]) for row in cur]

    def _select_region_vals(chip):
        select = """select Chr, Start, End, State, Cn, Im
            from %s
            where chip=%r
            """ % (tRegionCnv, chip)
        print select
        cur.execute(select)
        return [(e[0].encode(), e[1], e[2], [e[3].encode(), e[4], e[5]])
                for e in cur]

    chips = get_chips(cur, chiplike, tname=tRegionCnv)
    snp_rows = _select_snp_rows()
    cnv = Cnv(snp_rows)
    for chip in chips:
        region_vals = _select_region_vals(chip)
        snp_vals = cnv.snpValuesFromRegionValues(region_vals,
                default=['FM', 2, 0])
        for chr, block in groupby(snp_rows, itemgetter(0)):
            curr_vals = snp_vals[chr]
            for (chr, loc), (state, cn, im) in izip(block, curr_vals):
                cur.execute(insert, (chip, chr, loc, state, cn, im))

def view_Confusion(cur, exp_chip, field='State',
        tSnpCnv='SnpCnv'):

    view_exp = '''create temp table exp as
        select * from %s
        where Chip==%r and %s is not NULL and Chr not in ('6', '16')
        ''' % (tSnpCnv, exp_chip, field)

    view_comparison = '''create table vComparison_%(field)s as
        select %(tSnpCnv)s.Chip as Chip, 
            %(tSnpCnv)s.%(field)s as %(field)s,
            exp.%(field)s as Expect
        from %(tSnpCnv)s
        join exp using (Chr, Loc)
        where %(tSnpCnv)s.Chip != %(exp_chip)r
        ''' % locals()

    view_confusion = '''create view vConfusion_%(field)s as
        select Chip, %(field)s, Expect, count() as Frequency
        from vComparison_%(field)s
        group by Chip, %(field)s, Expect
        ''' % locals()

    cur.execute('drop table if exists exp')
    print view_exp
    cur.execute(view_exp)
    cur.execute('create index _idx on exp (Chr, Loc)')

    cur.execute('drop table if exists vComparison_%(field)s' % locals())
    print view_comparison
    cur.execute(view_comparison)
    cur.execute('''create index _idx_comp_%(field)s on vComparison_%(field)s
        (Chip, %(field)s, Expect)''' % locals())

    cur.execute('drop view if exists vConfusion_%(field)s' % locals())
    print view_confusion
    cur.execute(view_confusion)



#######
# old
def create_SnpRecovery_(cur, tname):
    create = '''create table %s
        (Chr, Chip, Value text, Match integer, Total integer, Recovery real,
        primary key (Chr, Chip, Value))
        ''' % tname
    cur.execute(create)

def insert_SnpRecovery_(cur, exp_chip, field='Cn',
        tSnpCnv='SnpCnv', tSnpRecovery='SnpRecovery', insert_or=''):
    cur.execute('''drop view if exists _exp''')
    cur.execute('''create view _exp as
        select * from %s
        where Chip==%r and %s is not NULL
        ''' % (tSnpCnv, exp_chip, field))
    for chip in get_chips(cur, exp_chip, complementary=True):
        cur.execute('''drop view if exists _obs''')
        cur.execute('''create view _obs as
            select * from SnpCnv
            where Chip==%r''' % (chip))
        insert = '''insert %(insert_or)s into %(tSnpRecovery)s
                (Chr, Chip, Value, Match, Total)
                select Chr, _obs.Chip, _exp.%(field)s,
                    sum(_exp.%(field)s==_obs.%(field)s),
                    count()
                from _exp
                    join _obs using (Chr, Loc)
                group by Chr, _exp.%(field)s
                ''' % locals()
        print insert
        cur.execute(insert)

def update_SnpRecovery(cur, tname='SnpRecovery'):
    update = '''update %s
        set Recovery= cast(Match as real)/Total''' % tname
    cur.execute(update)


def view_SnpRecoveryGlobal(cur, tname, vname=None, where=1):
    vname = vname or 'v%s_genomewide' % tname
    view = '''CREATE VIEW %s as
        select Chip, Value, sum(Match), sum(Total),
            cast(sum(Match) as real)/sum(Total)
        from %s
        where %s
        group by Chip, Value
        ''' % (vname, tname, where)
    print view
    cur.execute(view)






##
# old function
def recovery_from_db(cur, tExp, tObs, default=2):
    """
    tExp: with only one chip as exp
    tObs: with multiple chips as Obs
    """
    select = """select Loc, Im from %s"""
    # create the exp
    cur.execute(select % tExp)
    rows = list(cur)
    exp_ims = [val for loc, val in rows]
    exp_idxs = dict((loc, i) for i, (loc, val) in enumerate(rows))

    # get chips
    cur.execute("""select distinct chip from %s
        """ % tObs)
    chips = [row[0].encode() for row in cur]

    # for each chip
    result = {}
    for chip in chips:
        cur.execute("""select Loc, [Obs.Im] from %s where Chip=%r
            """ % (tObs, chip))
        obs_ims = [default] * len(exp_ims)
        for loc, obs_im in cur:
            obs_ims[exp_idxs[loc]] = obs_im
        result[chip] = snp_recovery(exp_ims, obs_ims)
    return result

def snp_recovery(exp, obs):
    """return the recovery rate for each state.

    - exp, obs: a list of states
    """
    refs = defaultdict(int)
    matches = defaultdict(int)
    for e, o in zip(exp, obs):
        refs[e] += 1
        if o == e: #match
            matches[e] += 1
    result = dict((e, (matches[e], refs[e]))
            for e in refs.iterkeys())
    return result

def recovery_from_SnpCnv(cur, exp_chip, field='Cn', tname='SnpCnv', default=2):
    """
    ref_chip: the chip used as reference
    field: the field used to calc recovery
    """
    select_exp = """select Loc, %s from %s where Chip==%r
        """ % (field, tname, exp_chip)
    # create the exp
    cur.execute(select_exp)
    rows = list(cur)
    exp_vals = [val for loc, val in rows]
    exp_idxs = dict((loc, i) for i, (loc, val) in enumerate(rows))

    # get other chips
    cur.execute("""select distinct chip from %s
        where Chip!=%r
        """ % (tname, exp_chip))
    chips = [row[0].encode() for row in cur]

    # for each chip
    result = {}
    for chip in chips:
        cur.execute("""select Loc, %s from %s where Chip=%r
            """ % (field, tname, chip))
        obs_vals = [default] * len(exp_vals)
        for loc, val in cur:
            obs_vals[exp_idxs[loc]] = val
        result[chip] = snp_recovery(exp_vals, obs_vals)
    return result

def recovery_from_SnpCnv_(cur, exp_chip, field='Cn', tname='SnpCnv', default=2):
    """
    ref_chip: the chip used as reference
    field: the field used to calc recovery
    """
    select_chip_chr = 'select Loc, %s from %s ' % (field, tname) +\
            'where Chip==%r and Chr=%r;'

    def get_other_chips():
        # get other chips
        cur.execute("""select distinct chip from %s
            where Chip!=%r
            """ % (tname, exp_chip))
        return [row[0].encode() for row in cur]

    def prepare_ref(chr):
        # create the exp
        cur.execute(select_chip_chr % (exp_chip, chr))
        rows = list(cur)
        exp_vals = [val for loc, val in rows]
        exp_idxs = dict((loc, i) for i, (loc, val) in enumerate(rows))
        return exp_vals, exp_idxs

    def chr_chip_recovery(chr, chip, ref_vals, loc_idxs):
        obs_vals = [default] * len(ref_vals)
        cur.execute(select_chip_chr % (chip, chr))
        for loc, val in cur:
            obs_vals[loc_idxs[loc]] = val
            return snp_recovery(ref_vals, obs_vals)

    chrs = map(str, range(1, 22+1))
    other_chips = get_other_chips()

    result = {} #{chr: {chip: {val: recov}}}
    for chr in chrs:
        result[chr] = {}
        ref_vals, loc_idxs = prepare_ref(chr)
        for chip in other_chips:
            result[chr][chip] = chr_chip_recovery(chr, chip, ref_vals, loc_idxs)
    return result


##
# performance is too poor to be useful
def insert_SnpCnv(cur, tname, vname, other_cond='rowid>-1'):
    insert = '''insert into %(vname)s
        select Chip, Chr, Loc, State, Cn, Im
        from SNP
            left join %(tname)s using (Chr)
        where Loc between Start and End
            and %(other_cond)s
        ''' % locals()
    print insert
    cur.execute(insert)

def create_SnpRecovery(cur, tname):
    create = '''create table %s
        (Chip text, Value text, Match integer, Total integer, Recovery real,
        primary key (Chip, Value))
        ''' % tname
    cur.execute(create)

def recovery_to_SnpRecovery_(cur, recovs, tname):
    """write the dct to a table file"""
    insert = """insert into %s
        values (?, ?, ?, ?, ?, ?)
        """ % tname
    for chr, chip_dct in recovs.iteritems():
        for chip, val_dct in chip_dct.iteritems():
            for val, (match, total) in val_dct.iteritems():
                cur.execute(insert, (chr, chip, val, match, total, match/total))

def recovery_to_SnpRecovery(cur, recovs, tname):
    """write the dct to a table file"""
    insert = """insert into %s
        values (?, ?, ?, ?, ?)
        """ % tname
    for chip, val_dct in recovs.iteritems():
        for val, (match, total) in val_dct.iteritems():
            cur.execute(insert, (chip, val, match, total, match/total))

