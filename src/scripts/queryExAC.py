
import tabix
from collections import defaultdict

exacPops = ('', 'AFR', 'AMR', 'EAS',
            'FIN', 'NFE', 'OTH', 'SAS')

def getPosFilter(chromStr, pos, tb):
    records = tb.query(chromStr, pos-1, pos)
    retLs = []
    for record in records:
        retLs.append(record[6])
    return retLs

def getPosFilterWithRefAlt(chromStr, pos, tb):
    records = tb.query(chromStr, pos-1, pos)
    retLs = []
    for record in records:
        retLs.append( [record[i] for i in (3, 4, 6)] )
    return retLs

def getPosDataAll(chromStr, pos, tb):
    """Lowest frac is 1 / 122516 = 8.162e-06.
       Uses data from all popultaions.
       tb is a tabix object.
       Returns a list of ref, altLs, freqLs"""
#    print(chromStr, pos-1, pos)
    records = tb.query(chromStr, pos-1, pos)
    retLs = []
    for record in records:
        chrom, pos, rs, ref, altLsStr = record[0:5]
        afLs = [float(x) for x in record[-1].split('AF=')[1].split(';')[0].split(',')]
        retLs.append( (ref, altLsStr.split(','), afLs) )
    return retLs

def getPosDataEuro(chromStr, pos, tb):
    """Uses data from AC_NFE (non finn euro popultaion).
       tb is a tabix object.
       Returns a list of ref, altLs, freqLs"""
    records = tb.query(chromStr, pos-1, pos)
    retLs = []
    for record in records:
        chrom, pos, rs, ref, altLsStr = record[0:5]
        acLs = [float(x) for x in record[-1].split('AC_NFE=')[1].split(';')[0].split(',')]
        an = float( record[-1].split('AN_NFE=')[1].split(';')[0] )
        afLs = []
        for ac in acLs:
            if an:
                afLs.append(ac/an)
            else:
                afLs.append('NA')
        retLs.append( (ref, altLsStr.split(','), afLs) )
    return retLs

def getPosDataByPop(chromStr, pos, tb):
    if chromStr == 'MT':
        return {} #chromStr = 'M'
    records = tb.query(chromStr, pos-1, pos)
    retLs = defaultdict(list)
    for record in records:
        chrom, pos, rs, ref, altLsStr = record[0:5]
        fields = exacPops
        for field in fields:
            if not field:
                vcfAc, vcfAn = 'AC=','AN='
            else:
                vcfAc, vcfAn = 'AC_' + field + '=', 'AN_' + field + '='
                
            if vcfAc in record[-1]:
                acLs = [float(x) for x in record[-1].split(vcfAc)[1].split(';')[0].split(',')]
                an = float( record[-1].split(vcfAn)[1].split(';')[0] )
                afLs = []
                for ac in acLs:
                    if an:
                        afLs.append(ac/an)
                    else:
                        afLs.append('NA')
                if not field:
                    retLs['Total'].append( (ref, altLsStr.split(','), afLs) )
                else:
                    retLs[field].append( (ref, altLsStr.split(','), afLs) )
    return retLs
        # AC=1169,1;
        # AN=122970;
        
        # AC_AFR=6,0;
        # AN_AFR=10468;

        # AC_AMR=649,0;
        # AN_AMR=11552;

        # AC_EAS=366,0;
        # AN_EAS=8730;

        # AC_FIN=22,1;
        # AN_FIN=6738;

        # AC_NFE=73,0;
        # AN_NFE=67552;

        # AC_OTH=11,0;
        # AN_OTH=924;

        # AC_SAS=42,0;
        # AN_SAS=16552

        # AF=9.506e-03,8.132e-06;

def testEachPop():
    pos = 46615880
    chromStr = '22'
    exac_file = '/nas/nbl3/human_vars/ExAC/ExAC.r0.2.sites.vep.vcf.gz'
    tb = tabix.open(exac_file)
#    getPosDataByPop(chromStr, pos, tb)
    records = getPosDataByPop(chromStr, pos, tb)
    for pop in records:
        for r in records[pop]:
            print(pop, r)

def testDel():
    pos = 11827113
    chromStr = '1'
    exac_file = '/nas/nbl3/human_vars/ExAC/ExAC.r0.2.sites.vep.vcf.gz'
    tb = tabix.open(exac_file)
    records = getPosDataAll(chromStr, pos, tb)
    for record in records:
        print(record)

def testAll():
    pos = 46615880
    chromStr = '22'
    exac_file = '/nas/nbl3/human_vars/ExAC/ExAC.r0.2.sites.vep.vcf.gz'
    tb = tabix.open(exac_file)
    records = getPosDataAll(chromStr, pos, tb)
    for record in records:
        print(record)

def testEuro():
    pos = 46615880
    chromStr = '22'
    exac_file = '/nas/nbl3/human_vars/ExAC/ExAC.r0.2.sites.vep.vcf.gz'
    tb = tabix.open(exac_file)
    records = getPosDataEuro(chromStr, pos, tb)
    for record in records:
        print(record)

def testAFR():
    pos = 46615880
    chromStr = '22'
    exac_file = '/nas/nbl3/human_vars/ExAC/ExAC.r0.2.sites.vep.vcf.gz'
    tb = tabix.open(exac_file)
    r = getPosDataByPop(chromStr, pos, tb)
    for pop in r:
        print(pop, r[pop])

if __name__ == "__main__":
#    testAFR()
#    testEachPop()
    # testAll()
    # testEuro()
    testDel()
