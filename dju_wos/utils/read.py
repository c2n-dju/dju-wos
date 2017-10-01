# -*- coding: utf-8 -*-

"""
/usr/bin/ipython3
import runpy
file_globals = runpy.run_path("read.py")

wos = sorted(file_globals['mismatchs'])[0]
pp = file_globals['pbUT'][wos]
p1 = file_globals['ciw'].publis[pp[0]]
p2 = file_globals['ciw'].publis[pp[1]]

def didi(kv):
    it = iter(kv)
    return dict(zip(it, it))

def mimi(p1, p2):
    d1 = didi(p1)
    d2 = didi(p2)
    dd1 = dict()
    dd2 = dict()
    keys = set(d1).intersection(d2)
    for k in keys:
        if d1[k] != d2[k]:
            dd1[k] = d1[k]
            dd2[k] = d2[k]
        d1.pop(k)
        d2.pop(k)
    return (dd1, dd2, d1, d2)

(dd1, dd2, d1, d2) = mimi(p1, p2)

c1 = dd1['CR'].split('\n')
c2 = dd2['CR'].split('\n')

dl = []
for (l1, l2) in zip(sorted(c1), sorted(c2)):
    if l1 != l2:
        dl.append((l1, l2))
"""

from functools import reduce

class CiwReader(object):

    """
    Comments TODO
    """
    def _read_fn(self, iline, line):
        i = line.find('FN ')
        if i < 0:
            self._error(iline, line, 'Usage Count (Last 180 Days)FN not found')
        self.header_residu = line[0:i]
        self.header_content = line[i+3:-1]
        self._read_line = self._read_vr
        # print('header_content = ' + self.header_content)

    def _read_vr(self, iline, line):
        if line[0:3] != 'VR ':
            self._error(iline, line, 'VR not found')
        self.version = line[3:-1]
        self._read_line = self._read_pt
        # print('version = ' + self.version)

    def _read_pt(self, iline, line):
        if line[0:3] == 'EF':
            self._seen_ef = True
            return
        if line == '\n':
            return
        if line[0:3] != 'PT ':
            self._error(iline, line, 'PT not found')
        self._publi = ['PT', line[3:]]
        self._read_line = self._read_next

    def _read_next(self, iline, line):
        if line[0] == ' ':
            if line[0:3] != '   ': # voir si la suite de clé à 3 caractères convient
                self._error(iline, line, 'bad next line')
            self._publi[-1] += line[3:]
        elif line == 'ER\n':
            if self._publi != []:
                self._publi[-1] = self._publi[-1][:-1]
            self.publis.append(self._publi)
            self._publi = []
            self._read_line = self._read_pt
        else:
            i = line.index(' ')
            key = line[:i]
            if line[i+1] == ' ':
                if key == 'FU':
                    pass
                elif key == 'RI':
                    pass
                else:
                    self._error(iline, line, "double white after key")
            self._publi[-1] = self._publi[-1][:-1]
            self._publi.append(key)
            self._publi.append(line[i+1:])

    def __init__(self, filename):
        self.publis = []
        self._filename = filename
        self._seen_ef = False
        self._read_line = self._read_fn
        f = open(self._filename, 'r')
        iline = 0
        while not self._seen_ef:
            line = f.readline()
            iline += 1
            self._read_line(iline, line)
        line = f.readline()
        if line != '':
            print('reste : ' + line)
        f.close()
    def _reconstruct_value(self, value):
        return '\n   '.join(value.split('\n')) + '\n'

    def _reconstruct_publi(self, publi):
        assert len(publi) % 2 == 0
        kvs = zip(publi[0::2], publi[1::2])
        r = reduce(lambda x, kv: x+kv[0]+' '+self._reconstruct_value(kv[1]), kvs, "")
        return r + 'ER\n\n'

    def _error(self, iline, line, message):
        raise Exception('Error: ' + message + ' in line ' + str(iline) + ": \"" + line + "\"")

    def reconstruct(self):
        ret = self.header_residu
        ret += 'FN ' + self.header_content + '\n'
        ret += 'VR ' + self.version + '\n'
        ret += reduce(lambda x, p: x+self._reconstruct_publi(p), self.publis, "")
        ret += 'EF'
        return ret

    def is_ok(self):
        f = open(self._filename, 'r')
        s1 = f.read()
        f.close()
        s2 = self.reconstruct()
        return s1 == s2

    def write_to(self, filename):
        f = open('/tmp/' + filename, 'w')
        f.write(self.reconstruct())
        f.close()


class Ciw(object):
    def __init__(self, ciws):
        self.header_content = ciws[0].header_content
        self.version = ciws[0].version
        self.publis = []
        for ciw in ciws:
            self.publis.extend(ciw.publis)


def keys(publi):
    ks = publi[0::2]
    ret = set(ks)
    assert len(ks) == len(ret) # sinon il y a des doublons C'est le cas avec RID
    return ret


def ok_doublons(ciw):
    ipub = -1
    for publi in ciw.publis:
        ipub += 1
        ks = publi[0::2]
        sks = set(ks)
        if len(ks) != len(sks):
            ksc = list(ks)
            ksc.sort()
            return (ipub, ksc)
    return None


def publi_all_keys(publi):
    return set(publi[0::2])


def publis_all_keys(ciw):
    return reduce(lambda s, publi: set.union(s, publi_all_keys(publi)), ciw.publis, set())

"""

PT Publication Type (conference, book, journal, book in series, or patent) : J S B -> existe toujours
UT Unique Article Identifier : WOS:000269196400018 -> existe toujours
DT Document Type : 'Article' 'Book Chapter' 'Correction' 'Editorial Material' 'Meeting Abstract' 'Proceeding Paper' 'Review'

TI Document Title
AB Abstract

J9 29-Character Source Abbreviation : 'APPL PHYS LETT'
JI ISO Source Abbreviation          : 'Appl. Phys. Lett.'
SO Publication Name                 : 'APPLIED PHYSICS LETTERS'
VL Volume
IS Issue
AR Article Number
BP Beginning Page
EP Ending Page
PG Page Count
PY Year Published
PD Publication Date : 'AUG 31' 'SEP-NOV'
DI Digital Object Identifier (DOI)

AU Authors          : 'Galopin, E'
AF Author Full Name : 'Galopin, Elisabeth'
C1 Author Address
RID ? Researcher ID ? Cf. http://www.researcherid.com Attention : champ répété

EM E-mail Address

CT Conference Title
CY Conference Date
HO Conference Host : 'Univ Zaragoza', 'Jaszowied Int Sch'

DE Author Keywords  : (498, 'silicon nanowires; silver nanoparticles; SERS; high sensitivity')
ID Keywords Plus®   : (498, 'SCATTERING; NANOPARTICLES; SERS; AG; MOLECULES')
SC Subject Category : (498, 'Science & Technology - Other Topics; Materials Science')
WC ? Domain ? :  (498, 'Nanoscience & Nanotechnology; Materials Science, Multidisciplinary')

NR Cited Reference Count
TC Times Cited

GA Document Delivery Number ? 726YT
Z9 ? number ?

FU Funding Agency and Grant Number
FX Funding Text (aussi acknowledgment)

BN ISBN
SN ISSN
PU Publisher : 'ROYAL SOC CHEMISTRY'
PA Publisher Address
PI Publisher City

RP Reprint Address

LA Language
BE ? Person ? (rare)
CL Conference Location : (très rare)
GP ? rare ? 'Optical Society of America' 'IEEE'
PN Part Number : rare [(78, 'Part 2'), (143, 'Part 2'), (375, 'Part 1-2')]
SE Book Series Title : 'Proceedings of SPIE'
SI Special Issue : 'SI'
SP Conference Sponsors : (413, 'Russian Fdn Basic Res, Russian Acad Sci')
SU Supplement ; rare (85, '2')

BS Book Series Subtitle
CR Cited References
ED Editors

VR Version Number
ER End of Record
EF End of File

MA Meeting Abstract Number (vide)
CA Group Authors (vide)
FN File Name (inexistant)

"""

def publis_by_keys(ciw):
    ret = dict()
    ipub = -1
    for publi in ciw.publis:
        ipub += 1
        for k, v in zip(publi[0::2], publi[1::2]):
            if not k in ret:
                ret[k] = []
            ret[k].append((ipub, v))
    return ret

def publis_by_PT(pbk):
    ret = dict()
    for i, v in pbk['PT']:
        if not v in ret:
            ret[v] = []
        ret[v].append(i)
    return ret

def publis_by_DT(pbk):
    ret = dict()
    for i, vs in pbk['DT']:
        for v in vs.split('; '):
            if not v in ret:
                ret[v] = []
            ret[v].append(i)
    return ret

def publis_by_OI(pbk):
    ret = dict()
    for i, vs in pbk['OI']:
        for v in vs.split(';'):
            v = v.strip(' \n')
            if len(v) == 0:
                continue
            lv = v.split('/')
            if len(lv) != 2:
                print('BAD ' + str(i) + ':' + v)
            vv = lv[1]
            if not vv in ret:
                ret[vv] = []
            ret[vv].append((lv[0].replace('\n', ' '), i))
    return ret

def publis_by_author(pbk, AUorAF):
    ret = dict()
    for i, v in pbk[AUorAF]:
        for au in v.split('\n'):
            if not au in ret:
                ret[au] = []
            ret[au].append(i)
    return ret

def publis_by_author_labo(pbk):
    ret1 = dict()
    ret2 = dict()
    for i, v in pbk['C1']:
        for al in v.split('\n'):
            if al[0] != '[':
                aus = ciw.publis[i][ciw.publis[i].index('AF') + 1]
                labo = al
                for a in aus.split('\n'):
                    if not a in ret1:
                        ret1[a] = []
                    if not labo in ret2:
                        ret2[labo] = []
                    ret1[a].append((i, labo))
                    ret2[labo].append((i, a))
            else:
                ic = al.index(']')
                if al[ic+1] != " ":
                    raise Exception('publi ' + str(i) + ', C1 = ' + al)
                labo = al[ic+2:]
                for a in al[1:ic].split('; '):
                    if not a in ret1:
                        ret1[a] = []
                    if not labo in ret2:
                        ret2[labo] = []
                    ret1[a].append((i, labo))
                    ret2[labo].append((i, a))
    return ret1, ret2


def publis_by_journal(pbk, JIorJ9):
    ret = dict()
    for i, v in pbk[JIorJ9]:
        if not v in ret:
            ret[v] = []
        ret[v].append(i)
    return ret

def publis_by_wos(pbk):
    ret = dict()
    for i, v in pbk['UT']:
        if not v in ret:
            ret[v] = []
        ret[v].append(i)
    return ret

def maxsizes(ciw, keys):
    ret = {}
    for publi in ciw.publis:
        for k, v in zip(publi[0::2], publi[1::2]):
            l = len(v)
            if (not k in ret) or (ret[k][0] < l):
                ret[k] = (l, v)
    return ret

"""
 http://images.webofknowledge.com/WOK46/help/WOS/h_fieldtags.html
 http://images.webofknowledge.com/WOKRS53B4/help/WOS/hs_wos_fieldtags.html
https://images.webofknowledge.com/WOKRS512B4/help/WOK/hs_alldb_fieldtags.html
 http://images.webofknowledge.com/WOKRS56B5/help/WOS/hs_wos_fieldtags.html
 http://images.webofknowledge.com/WOKRS58B4/help/WOK/hs_wos_fieldtags.html
"""

# Cf. (a) http://images.webofknowledge.com/WOKRS515B5/help/WOK/hs_wos_fieldtags.html
#     (b) http://images.webofknowledge.com/WOK46/help/WOS/h_fieldtags.html
#     (c) https://images.webofknowledge.com/images/help/WOS/hs_wos_fieldtags.html
#     (d) FP

wos_tags = """
a AB Abstract
a AE Patent Assignee
b AF Author Full Name
a AR Article Number
a AU Author(s) (in English)
a BA Book Author(s)
a BE Editor(s)
a BF Book Authors Full Name
a BN International Standard Book Number (ISBN)
a BP Beginning Page
a BS Book Series Subtitle
b C1 Author Address
a CA Group Author(s)
a CL Conference Location
b CR Cited References
a CT Conference Title
a CY Conference Date
a D2 Book Digital Object Identifier (DOI)
d DA Access Date
b DE Author Keywords
a DI Digital Object Identifier (DOI)
b DT Document Type
b ED Editors
a EF End of File
a EI Electronic International Standard Serial Number (eISSN)
b EM E-mail Address
a EP Ending Page
a ER End of Record
a FN File Name
a FT Foreign Title
b FU Funding Agency and Grant Number
b FX Funding Text
b GA Document Delivery Number
a GP Book Group Author(s)
b HO Conference Host
b ID Keywords Plus®
a IS Issue
b J9 29-Character Source Abbreviation
b JI ISO Source Abbreviation
b LA Language
a MA Meeting Abstract Number
b NR Cited Reference Count
d OA No or gold
a OI ORCID Identifier (Open Researcher and Contributor ID)
a PG Chapter Count (Book Citation Index)
b PA Publisher Address
a PD Publication Date
a PG Page Count
b PI Publisher City
c PM PubMed ID
a PN Part / Patent Number
a PT Publication Type (J=Journal; B=Book; S=Series; P=Patent)
b PU Publisher
a PY Publication Year
a RI ResearcherID Number
b RP Reprint Address
a SC Research Areas
a SE Book Series Title
a SI Special Issue
a SN International Standard Serial Number (ISSN)
a SO Full Source Title (includes title and subtitle)
a SP Conference Sponsor
a SU Supplement
a TC Times Cites in Web of Science Core Collection
a TI Title
c U1 Usage Count (Last 180 Days)
c U2 Usage Count (Since 2013)
a UT Accession Number / ISI Unique Article Identifier
a VL Volume
a VR Version Number
a WC Web of Science Category
a Z1 Title (in second language)
a Z2 Author(s) (in second language)
a Z3 Full Source Title (in second language) (includes title and subtitle)
a Z4 Abstract (in second language)
a Z8 Times Cited in Chinese Science Citation Database
a Z9 Total Times Cited (Web of Science, BIOSIS Citation Index, and Chinese Science Citation Database)
a ZB Times Cited in BIOSIS Citation Index
"""
wos_tags_dict = {}
for (k, v) in map(lambda x: (x[2:4], x[5:]), wos_tags.strip('\n').split('\n')):
    wos_tags_dict[k] = v

ciws = []

for f in ("marcoussis_1.ciw",
          "marcoussis_2013.ciw",
          "marcoussis_2.ciw",
          "marcoussis_3.ciw",
          "marcoussis_4.ciw",
          "marcoussis_5.ciw",
          "marcoussis_6.ciw",
          "marcoussis_7.ciw",
          "marcoussis_8.ciw",
          "marcoussis_2013.ciw",
          "2013-2014.ciw",
          "c2n_orsay.ciw",
          "c2n_special.ciw",
          "el-fond-00501-001000.ciw",
          "el-fond-00001-00500.ciw",
          "el-fond-00501-001092.ciw",
          "ief_orsay_00501-00622.ciw",
          "ief_orsay_00001-00500.ciw",
          "marcoussis-00001-00500.ciw",
          "marcoussis-01001-01500.ciw",
          "marcoussis-00501-01000.ciw",
          "marcoussis-01501-02000.ciw",
          "marcoussis-02501-03000.ciw",
          "marcoussis-02001-02500.ciw",
          "marcoussis-03001-03500.ciw",
          "marcoussis-03501-04000.ciw",
          "marcoussis-04001-04500.ciw",
          "marcoussis-05001-05123.ciw",
          "marcoussis-04501-05000.ciw",
          "oo-ief_orsay_00001-00025.ciw",
          "oo-umr_8622-00001-00002.ciw",
          "umr_8622-00001-00500.ciw",
          "umr8622-00001-00192.ciw",
          "umr_8622-00501-01000.ciw",
          "umr_8622-01001-01387.ciw",
         ):
    filename = '../data/' + f
    ciw = CiwReader(filename)
    print(str(ciw.is_ok()) + ' ' + filename)
    ciws.append(ciw)

ciw = Ciw(ciws)

pbk = publis_by_keys(ciw)
pbUT = publis_by_wos(pbk)
assert len(pbk['UT']) == len(ciw.publis)
publis = []
mismatchs = set()
for ut in sorted(pbUT.keys()):
    ps = pbUT[ut]
    pm = [ps[0]]
    publis.append(ciw.publis[ps[0]])
    for p in ps[1:]:
        mismatch = True
        for po in pm:
            if ciw.publis[p] == ciw.publis[po]:
                mismatch = False
        if mismatch:
            print(ut + " mismatch")
            publis.append(ciw.publis[p])
            pm.append(p)
            mismatchs.add(ut)
        else:
            print(ut + " match")

print(len(ciw.publis))
print(len(publis))
ciw.publis = publis

all_keys = list(publis_all_keys(ciw))
all_keys.sort()
print()
for k in all_keys:
    if k in wos_tags_dict:
        print(k + ' ' + wos_tags_dict[k])
    else:
        print(k + ' ' + '--- UNKNOWN ---')

pbk = publis_by_keys(ciw)
pbPT = publis_by_PT(pbk) # S J B
pbDT = publis_by_DT(pbk)

dts = sorted(pbDT.keys())
if dts != ['Article',
           'Biographical-Item',
           'Book Chapter',
           'Correction',
           'Correction, Addition',
           'Discussion',
           'Editorial Material',
           'Letter',
           'Meeting Abstract',
           'News Item',
           'Note',
           'Proceedings Paper',
           'Reprint',
           'Review',]:
    print()
    print('!! nouveaux DTs !!')
    print(sorted(pbDT.keys()))

print()
pbAU = publis_by_author(pbk, 'AU')
pbAF = publis_by_author(pbk, 'AF')
pbUT = publis_by_wos(pbk)
pbC1_au, pbC1_labo = publis_by_author_labo(pbk)
pbJ9 = publis_by_journal(pbk, 'J9')
pbJI = publis_by_journal(pbk, 'JI')
journals = list(pbJ9.keys())
journals.sort()
sizes = maxsizes(ciw, all_keys)
ks = list(sizes.keys())
ks.sort(key=lambda x: sizes[x][0])
for k in ks:
    print(k, sizes[k])
authors = list(pbAU.keys())
authors.sort(key=str.lower)
assert len(pbk['UT']) == len(ciw.publis)

print(len(ciw.publis))

poi = publis_by_OI(pbk)

names_from_oi = dict()
for oi in poi.keys():
    names = []
    for v in poi[oi]:
        name = v[0]
        if not name in names:
            names.append(name)
    names_from_oi[oi] = names

oi_doublons = []
for oi in names_from_oi.keys():
    names = names_from_oi[oi]
    if len(names) != 1:
        oi_doublons.append(oi)

for oi in oi_doublons:
    print(oi, names_from_oi[oi])

oi_from_name = dict()
for oi in poi.keys():
    for v in poi[oi]:
        name = v[0]
        if not name in oi_from_name:
            oi_from_name[name] = []
        if not oi in oi_from_name[name]:
            oi_from_name[name].append(oi)

oi_authors = list(oi_from_name.keys())
oi_authors.sort(key=str.lower)

for author in oi_authors:
    print(author, oi_from_name[author])

homonymes = dict()
for author in oi_authors:
    if len(oi_from_name[author]) != 1:
        homonymes[author] = oi_from_name[author]

print(homonymes)

no_doi = set(range(len(ciw.publis)))

for p in pbk['DI']:
    no_doi.remove(p[0])

no_doi = sorted(list(no_doi))

di = pbk['DI']
