



import sys
from itertools import *
from pprint import pprint
import cplex
from math import *
import networkx as nx
from cplex.exceptions import CplexSolverError


ALPHA = 0.1
BETA = 0.00001



matfile = sys.argv[1]
outfile = sys.argv[2]

with open(matfile) as f:
    matdata = f.readlines()

matdata = map(lambda x: x.strip().split(), matdata)

mutation_matrix = []
for line in matdata:
    mutation_matrix.append(map(int, line))





nrMutations = len(mutation_matrix) # number of rows - number of mutations
nrCells = len(mutation_matrix[0]) # number of cells

mutationPairs = list(combinations(range(nrMutations), 2))
cellTriples = list(combinations(range(nrCells), 3))



columns = []



for c2 in range(nrCells):
    c_column = [x[c2] for x in mutation_matrix]
    c_muts = []
    for order, elem in enumerate(c_column):
        if elem == 1:
	    c_muts.append(order)
    for i in range(len(c_muts)):
	for j in range(i + 1, len(c_muts)):
            m1_cells = []
	    m2_cells = []
	    for order1, elem1 in enumerate(mutation_matrix[c_muts[i]]):
		if elem1 == 1 and mutation_matrix[c_muts[j]][order1] == 0:
		    m1_cells.append(order1)
	    for order2, elem2 in enumerate(mutation_matrix[c_muts[j]]):
		if elem2 == 1 and mutation_matrix[c_muts[i]][order2] == 0:
		    m2_cells.append(order2)
	    m1 = c_muts[i]
	    m2 = c_muts[j]
	    for c1 in m1_cells:
		for c3 in m2_cells:
		    column = [0 for iii in range(nrCells * nrMutations)]
		    column[nrCells * m1 + c1] = 1
		    column[nrCells * m1 + c2] = 1
		    column[nrCells * m1 + c3] = 1
		    column[nrCells * m2 + c1] = 1
		    column[nrCells * m2 + c2] = 1
		    column[nrCells * m2 + c3] = 1
		    columns.append(column)







hitMatrix = []
for i in range(len(columns[0])):
    line = [x[i] for x in columns]
    hitMatrix.append(line)

hitNames = []
for i in range(nrMutations):
    for j in range(nrCells):
        hitNames.append((i, j))

hitZeroOne = []
for i in range(nrMutations):
    for j in range(nrCells):
        hitZeroOne.append(mutation_matrix[i][j])



#for name, line in zip(hitNames, hitMatrix):
#    print name, sum(line)




"""
emptyToRemove = []
for hitName, hitZo, hitLine in zip(hitNames, hitZeroOne, hitMatrix):
    if hitZo == 0 and sum(hitLine) > 0:
        # check which one to remove
        #print hitLine
        inDoubleUs = []
        for index, belongsToW in enumerate(hitLine):
            if belongsToW:
                inDoubleUs.append(index)
        #print inDoubleUs
        emptyCountDict = {}
        for index in inDoubleUs:
            w = [hitMatrix[i][index] for i in range(len(hitMatrix))]
            for hn, hzo, value in zip(hitNames, hitZeroOne, w):
                if hzo == 0 and value == 1 and hitName != hn:
                    if hn not in emptyCountDict:
                        emptyCountDict[hn] = 1
                    else:
                        emptyCountDict[hn] += 1
        print emptyCountDict, len(inDoubleUs)
        toRemove = False
        for x, y in emptyCountDict.items():
            if y >= len(inDoubleUs) and y > 1:
                toRemove = True
                break
        if toRemove:
            print "TOREMOVE", hitName, nrCells * hitName[0] + hitName[1]
            emptyToRemove.append(hitName)
        else:
            print "TOLEAVE", hitName


print "-------------------------------------------------------------"


nonEmptyToRemove = []
for hitName, hitZo, hitLine in zip(hitNames, hitZeroOne, hitMatrix):
    if hitZo == 1 and sum(hitLine) > 0:
        # check which one to remove
        #print hitLine
        inDoubleUs = []
        for index, belongsToW in enumerate(hitLine):
            if belongsToW:
                inDoubleUs.append(index)
        #print inDoubleUs
        nonEmptyCountDict = {}
        for index in inDoubleUs:
            w = [hitMatrix[i][index] for i in range(len(hitMatrix))]
            for hn, hzo, value in zip(hitNames, hitZeroOne, w):
                if hzo == 1 and value == 1 and hitName != hn:
                    if hn not in nonEmptyCountDict:
                        nonEmptyCountDict[hn] = 1
                    else:
                        nonEmptyCountDict[hn] += 1
        print nonEmptyCountDict, len(inDoubleUs)
        toRemove = False
        for x, y in nonEmptyCountDict.items():
            if y >= len(inDoubleUs) and y > 1:
                toRemove = True
                break
        if toRemove:
            print "TOREMOVE", hitName, nrCells * hitName[0] + hitName[1]
            nonEmptyToRemove.append(hitName)
        else:
            print "TOLEAVE", hitName







print "-------------------------------------------------------------"


nonEmptyToRemove2 = []
for hitName, hitZo, hitLine in zip(hitNames, hitZeroOne, hitMatrix):
    if hitZo == 1 and sum(hitLine) > 0:
        # check which one to remove
        #print hitLine
        inDoubleUs = []
        for index, belongsToW in enumerate(hitLine):
            if belongsToW:
                inDoubleUs.append(index)
        #print inDoubleUs
        nonEmptyCountDict2 = {}
        for index in inDoubleUs:
            w = [hitMatrix[i][index] for i in range(len(hitMatrix))]
            for hn, hzo, value in zip(hitNames, hitZeroOne, w):
                if hzo == 0 and value == 1 and hitName != hn:
                    if hn not in nonEmptyCountDict2:
                        nonEmptyCountDict2[hn] = 1
                    else:
                        nonEmptyCountDict2[hn] += 1
        print nonEmptyCountDict2, len(inDoubleUs)
        howManyEmptyCovered = 0
        for x, y in nonEmptyCountDict2.items():
            if y >= len(inDoubleUs) and y > 1:
                howManyEmptyCovered += 1
        if howManyEmptyCovered in [1, 2, 3]:
            print "TOREMOVE", hitName, nrCells * hitName[0] + hitName[1]
            nonEmptyToRemove2.append(hitName)
        else:
            print "TOLEAVE", hitName
"""




















cpx = cplex.Cplex()
cpx.set_results_stream("solution.txt")

cpx.parameters.mip.pool.absgap.set(0.0)
cpx.parameters.mip.pool.intensity.set(4)
cpx.parameters.mip.limits.populate.set(1)
cpx.parameters.timelimit.set(300)
#cpx.parameters.timelimit.set(360)

variables = []

# add variables
for name, line in zip(hitNames, hitMatrix):
    if 1 in set(line):
        name = "X#%s#%s" % name
        variables.append(name)
        cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])


# add constraints
for colNum, column in enumerate(columns):
    constraint = []
    for name, mutValue in zip(hitNames, column):
        if mutValue == 1:
            constraint.append(name)
    inds = map(lambda z: "X#%s#%s" % z, constraint)
    vals = [1 for z in constraint]
    cons = cplex.SparsePair(ind=inds, val=vals)
    cpx.linear_constraints.add( \
        lin_expr = [cons],\
        senses = ["G"],\
        rhs = [1],\
        names = ['cons-%s' % colNum]\
    )


"""
# remove empties to remove
for a, b in emptyToRemove:
    inds = ["X#%s#%s" % (a, b)]
    vals = [1]
    cons = cplex.SparsePair(ind=inds, val=vals)
    cpx.linear_constraints.add( \
        lin_expr = [cons],\
        senses = ["E"],\
        rhs = [0],\
        names = ['cons-%s' % colNum]\
    )


# remove empties to remove
for a, b in nonEmptyToRemove:
    inds = ["X#%s#%s" % (a, b)]
    vals = [1]
    cons = cplex.SparsePair(ind=inds, val=vals)
    cpx.linear_constraints.add( \
        lin_expr = [cons],\
        senses = ["E"],\
        rhs = [0],\
        names = ['cons-%s' % colNum]\
    )


# remove empties to remove
for a, b in nonEmptyToRemove2:
    inds = ["X#%s#%s" % (a, b)]
    vals = [1]
    cons = cplex.SparsePair(ind=inds, val=vals)
    cpx.linear_constraints.add( \
        lin_expr = [cons],\
        senses = ["E"],\
        rhs = [0],\
        names = ['cons-%s' % colNum]\
    )
"""





for var in variables:
    ind1, ind2 = tuple(map(int, var.split("#")[1:]))
    if mutation_matrix[ind1][ind2] == 0:
        coef = log((1 - BETA) / ALPHA)
    elif mutation_matrix[ind1][ind2] == 1:
        coef = log((1 - ALPHA) / BETA)
    cpx.objective.set_linear(var, coef)

cpx.objective.set_sense(cpx.objective.sense.minimize)

cpx.set_problem_type(cpx.problem_type.MILP)
cpx.write("program.txt", filetype="lp")
cpx.solve()

cpx.populate_solution_pool()


objVals = []
solutions = []


optimum = None

for ii in range(cpx.solution.pool.get_num()):

    objVals.append(cpx.solution.pool.get_objective_value(ii))

    solution = []
    for var in variables:
        if int(round(cpx.solution.pool.get_values(ii, var))) == 1:
            solution.append(var)
    solutions.append(solution)





minObjValue = min(objVals)
goodIndices = [x for x, y in enumerate(objVals) if y == minObjValue]




good_solutions = []
count = 0
for sol in solutions:
    if count in goodIndices:
        good_solutions.append(sol)
    count += 1
good_solutions2 = [good_solutions[0]]
for i in range(1, len(good_solutions)):
    new = True
    for x in good_solutions2:
        if set(x) == set(good_solutions[i]):
            new = False
            break
    if new:
        good_solutions2.append(good_solutions[i])



print good_solutions2


for sol in good_solutions2:
    nrZeros = 0
    nrOnes = 0
    for var in sol:
        ind1, ind2 = tuple(map(int, var.split("#")[1:]))
        if mutation_matrix[ind1][ind2] == 0:
            nrZeros += 1
        elif mutation_matrix[ind1][ind2] == 1:
            nrOnes += 1


change = []
for var in sol:
    ind = tuple(map(int, var.split("#")[1:]))
    change.append(ind)

# produce corrected matrix
corrected_matrix = []
for i in range(len(mutation_matrix)):
    line = []
    for j in range(len(mutation_matrix[0])):
        line.append(0)
    corrected_matrix.append(line)
for i in range(len(mutation_matrix)):
    for j in range(len(mutation_matrix[0])):
        if (i, j) in change:
            corrected_matrix[i][j] = 1 - mutation_matrix[i][j]
        else:
            corrected_matrix[i][j] = mutation_matrix[i][j]


# compute number zeros and ones
#for sol in good_solutions2:
#    nrZeros = 0
#    nrOnes = 0
#    for var in sol:
#        ind1, ind2 = tuple(map(int, var.split("#")[1:]))
#        if mutation_matrix[ind1][ind2] == 0:
#            nrZeros += 1
#        elif mutation_matrix[ind1][ind2] == 1:
#            nrOnes += 1


with open(outfile, "w") as f:
    for line in corrected_matrix:
        f.write("%s\n" % " ".join(map(str, line)))







