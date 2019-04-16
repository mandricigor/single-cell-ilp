
import sys
import random

# read the simulated reads vcf file
with open(sys.argv[1]) as f:
    vcf1 = f.readlines()


header = map(lambda x: x.strip(), vcf1[:3])
vcf1 = map(lambda x: x.strip().split(), vcf1[3:])




table = []

for line in vcf1:
    #line2 = line[:9]
    #line2[-1] = "GT"
    row = []
    for position in line[9:]:
        #print position, gen, "AAAAA"
        likelihoods = map(float, position.split(","))
        if likelihoods == [0.0, 0.0, 0.0]:
            row.append("3")
        else:
            maximum_likelihood = max(likelihoods)
            if likelihoods[0] == maximum_likelihood:
                row.append("0")
            elif likelihoods[1] == maximum_likelihood:
                row.append("1")
            elif likelihoods[2] == maximum_likelihood:
                row.append("2")
    #line2 += row
    #table.append(line2)
    table.append(row)


#for line in header:
#    print line

for line in table:
    print " ".join(line)




