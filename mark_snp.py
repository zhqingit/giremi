#!/usr/bin/python


def main(argv):

    snpF = ""
    varF = ""
    gF = ""

    try:
        opts, args = getopt.getopt(argv,"s:g:i:",["snpMark"])
    except getopt.GetoptError as err:
        sys.exit(2)

    for opt, arg in opts:
       if opt == '-h':
           print 'Usage: mark_snp.py -s -i'
           print ' -s, --snp    file    The SNP file (vcf format)'
           print ' -g, --gene   file    The ucsc table with gene annotation (vcf format)'
           print ' -i, --vcf    file    The variants file (vcf format)'

       elif opt == '-s':
           snpF = arg
       elif opt == '-i':
           varF = arg
       elif opt == '-g':
           gF = arg

    snp = {}
    with open(snpF,"r") as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith("#"):
                continue
            line = tmp.split("\t")
            chrom = line[0]
            pos = line[1]
            ref = line[3]
            alt = line[4]
            iid = chrom+"|"+pos+"|"+ref+"|"+alt
            snp[iid] = 1

    gstr = {}
    with open(gF,"r") as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith("#"):
                continue
            line = tmp.split("\t")
            strand = line[3]
            gene = line[12]
            gstr[gene] = strand

    ANN = re.compile(r".*ANN=(.*)")
    with open(varF,"r") as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith("#"):
                continue
            line = tmp.split("\t")
            chrom = line[0]
            pos = line[1]
            ref = line[3]
            alt = line[4]
            iid = chrom+"|"+pos+"|"+ref+"|"+alt

            info = line[7]

            ann = ""
            m = ANN.match(info)
            if m:
                ann = m.group(1)
            ll = ann.split(",")[0].split("|")
            gene = ll[3]

            infoP = chrom+"\t"+str(int(pos)-1)+"\t"+pos+"\t"+gene
            if iid in snp:
                infoP = infoP+"\t1"
            else:
                infoP = infoP+"\t0"
            
            if gene in gstr:
                print(infoP+"\t"+gstr[gene])
            
if __name__ == "__main__":
    main(sys.argv[1:])

