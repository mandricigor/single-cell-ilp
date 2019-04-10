


import sys
import numpy as np


EPSILON = 0.001
DROPOUT = 0.1
CELL_COVERAGE_MEAN = 100
CELL_COVERAGE_STD = 20


input_csv = sys.argv[1]
output_vcf = sys.argv[2]
output_scite = sys.argv[3]



with open(input_csv) as f:
    a = f.readlines()
a = map(lambda x: x.strip().split(), a)
genotypes = map(lambda x: map(int, x), a)


n_cells = len(a)
n_muts = len(a[0])


likelihoods = []
scite_genotypes = []

cell_coverages = map(int, CELL_COVERAGE_STD * np.random.randn(n_cells) + CELL_COVERAGE_MEAN)


for cell in range(n_cells):
    cell_mutations = genotypes[cell]
    cell_coverage = cell_coverages[cell]
    cell_likelihoods = []
    cell_scite_genotypes = []
    for mutation in cell_mutations:
        mutation_coverage = np.random.poisson(cell_coverage)
        if mutation == 0:
            # generate dropout
            dr = 1 - DROPOUT ** 2
            x = np.random.random()
            if x > dr:
                # dropout
                cell_likelihoods.append("0.0,0.0,0.0")
                cell_scite_genotypes.append(3)
            else:
                # generate reads
                reads = np.random.binomial(1, EPSILON, mutation_coverage)
                num0s = len([x for x in reads if x == 0])
                num1s = mutation_coverage - num0s
                like1 = (np.log10(1 - EPSILON) * num0s) + (np.log10(EPSILON) * num1s)
                like2 = np.log10(0.5) * mutation_coverage
                like3 = (np.log10(1 - EPSILON) * num1s) + (np.log10(EPSILON) * num0s)
                cell_likelihoods.append(",".join(map(str, [like1, like2, like3])))
                if like1 == max([like1, like2, like3]):
                    cell_scite_genotypes.append(0)
                elif like2 == max([like1, like2, like3]):
                    cell_scite_genotypes.append(1)
                else:
                    cell_scite_genotypes.append(2)
        elif mutation == 1:
            dr1 = DROPOUT * (1 - DROPOUT)
            dr2 = DROPOUT * (1 - DROPOUT)
            dr3 = DROPOUT * DROPOUT
            dr4 = (1 - DROPOUT) * (1 - DROPOUT)
            x = np.random.random()
            if x < dr1:
                reads = np.random.binomial(1, EPSILON, mutation_coverage)
                num0s = len([x for x in reads if x == 0])
                num1s = mutation_coverage - num0s
                like1 = (np.log10(1 - EPSILON) * num0s) + (np.log10(EPSILON) * num1s)
                like2 = np.log10(0.5) * mutation_coverage
                like3 = (np.log10(1 - EPSILON) * num1s) + (np.log10(EPSILON) * num0s)
                cell_likelihoods.append(",".join(map(str, [like1, like2, like3])))
                if like1 == max([like1, like2, like3]):
                    cell_scite_genotypes.append(0)
                elif like2 == max([like1, like2, like3]):
                    cell_scite_genotypes.append(1)
                else:
                    cell_scite_genotypes.append(2)
            elif x < dr1 + dr2:
                reads = np.random.binomial(1, 1 - EPSILON, mutation_coverage)
                num0s = len([x for x in reads if x == 0])
                num1s = mutation_coverage - num0s
                like1 = (np.log10(1 - EPSILON) * num0s) + (np.log10(EPSILON) * num1s)
                like2 = np.log10(0.5) * mutation_coverage
                like3 = (np.log10(1 - EPSILON) * num1s) + (np.log10(EPSILON) * num0s)
                cell_likelihoods.append(",".join(map(str, [like1, like2, like3])))
                if like1 == max([like1, like2, like3]):
                    cell_scite_genotypes.append(0)
                elif like2 == max([like1, like2, like3]):
                    cell_scite_genotypes.append(1)
                else:
                    cell_scite_genotypes.append(2)
            elif x < dr1 + dr2 + dr3:
                cell_likelihoods.append("0.0,0.0,0.0")
                cell_scite_genotypes.append(3)
            else:
                reads = np.random.binomial(1, 0.5, mutation_coverage)
                num0s = len([x for x in reads if x == 0])
                num1s = mutation_coverage - num0s
                like1 = (np.log10(1 - EPSILON) * num0s) + (np.log10(EPSILON) * num1s)
                like2 = np.log10(0.5) * mutation_coverage
                like3 = (np.log10(1 - EPSILON) * num1s) + (np.log10(EPSILON) * num0s)
                cell_likelihoods.append(",".join(map(str, [like1, like2, like3])))
                if like1 == max([like1, like2, like3]):
                    cell_scite_genotypes.append(0)
                elif like2 == max([like1, like2, like3]):
                    cell_scite_genotypes.append(1)
                else:
                    cell_scite_genotypes.append(2)
    likelihoods.append(cell_likelihoods)
    scite_genotypes.append(cell_scite_genotypes)


likelihoods_transposed = []
for i in range(len(likelihoods[0])):
    likelihoods_transposed.append([x[i] for x in likelihoods])

scite_genotypes_transposed = []
for i in range(len(scite_genotypes[0])):
    scite_genotypes_transposed.append([x[i] for x in scite_genotypes])


with open(output_scite, "w") as f:
    for line in scite_genotypes_transposed:
        f.write("%s\n" % " ".join(map(str, line)))



with open(output_vcf, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("##fileDate=20180615\n")
    f.write("##source=PLINKv1.90\n")
    f.write("##contig=<ID=1,length=249225078>\n")
    f.write('##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
    f.write("%s\n" % "\t".join(map(lambda x: "ID%s" % x, range(1, len(likelihoods_transposed[0]) + 1))))
    i = 1
    for line in likelihoods_transposed:
        f.write("1\t%s\trs%s\tA\tC\t.\t.\tPR\tGL\t" % (12345 + 55 * i, i))
        f.write("%s\n" % "\t".join(line))
        i += 1




