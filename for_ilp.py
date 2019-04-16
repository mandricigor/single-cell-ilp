#inputfile = "/home/igorm/workspace/recomb_2019/beagle.imputed.test1_true_igor.nick.vcf"
#inputfile = "/home/igorm/workspace/CliqueSNV/tmp/beagle.imputed.test1_true_igor.nick.vcf/test1_true_igor.nick.vcf"

#inputfile = "/home/igorm/workspace/recomb_2019/test1_true_igor.nick.vcf"
#outputfile = "/home/igorm/workspace/recomb_2019/beagle.imputed.test1_true_igor.nick.for_ilp.csv"


import sys


inputfile = sys.argv[1]
outputfile = sys.argv[2]


with open(inputfile) as f:
    a = f.readlines()
a = [x for x in a if not x.startswith("#")]
a = map(lambda x: x.strip().split()[9:], a)

#snps = []

#for i in range(len(a[0])):
#    snps.append([x[i] for x in a])

def lik_to_gen(x):
    x = x.split(":")[0]
    x = x.split("|")
    x = map(int, x)
    x = sum(x)
    if x == 2:
        x = 1
    return str(x)
    
    
output = []
with open(outputfile, "w") as f:
    for x in a:
        f.write("%s\n" % " ".join(map(lik_to_gen, x)))

