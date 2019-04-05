
# this is the class which solves orientation problem

import sys
import cplex
import math
import networkx as nx
from cplex.exceptions import CplexSolverError
from collections import Counter
from Bio.Seq import Seq
from Queue import Queue
from pprint import pprint
import matplotlib.pyplot as plt
import pygraphviz as PG
from math import log


ALPHA = 0.1
BETA = 0.001


def single_cell_phylogeny(matrix, nLoss=1):
    cpx = cplex.Cplex()
    cpx.set_results_stream("solution.txt")

    #cpx.parameters.mip.pool.absgap.set(0.0)
    #cpx.parameters.mip.pool.intensity.set(4)
    #cpx.parameters.mip.limits.populate.set(1000)
    #cpx.parameters.timelimit.set(360)

    #cpx.parameters.mip.pool.replace.set(1)
    #cpx.parameters.emphasis.mip.set(2)
    #cpx.parameters.mip.pool.capacity = 3

    # add root to the mutation matrix - the last row is all 0's
    #matrix.append([1] * len(matrix[0]))

    matrix = [[1] * len(matrix[0])] + matrix



    M = len(matrix) # number of rows - number of mutations
    K = len(matrix[0]) # number of cells


    # adding Xij variables (matrix M)
    for i in range(M):
        for j in range(M):
            name = "X#%s#%s" % (i, j)
            cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])

    # adding Ui
    for i in range(M):
        name = "U#%s" % i
        cpx.variables.add(lb=[0], ub=[1], types=["C"], names=[name])

    # adding Vj
    for j in range(M):
        name = "V#%s" % j
        cpx.variables.add(lb=[0], ub=[1], types=["C"], names=[name])

    for i in range(M):
        for j in range(M):
            if i < j:
                name = "A#%s#%s" % (i, j)
                cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])

                name = "B#%s#%s" % (i, j)
                cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])

                name = "C#%s#%s" % (i, j)
                cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])

                name = "D#%s#%s" % (i, j)
                cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])

    # adding Sik
    for i in range(M):
        for k in range(K):
            name = "S#%s#%s" % (i, k)
            cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])


    # adding Pik - single cell matrix we are looking for
    for i in range(M):
        for k in range(K):
            name = "P#%s#%s" % (i, k)
            cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])


    # adding Gikj (Column i of P is equal to product of matrix M and column i of S) ????????
    for i in range(M):
        for k in range(K):
            for j in range(M):
                name = "G#%s#%s#%s" % (i, k, j)
                cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])


    # adding Lij - mutation loss matrix (Lij = 1 iff mutation i is lost before j)
    for i in range(M):
        for j in range(M):
            name = "L#%s#%s" % (i, j)
            cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])

    # adding Rij - mutation-descendants matrix (R = L * M)
    for i in range(M):
        for j in range(M):
            name = "R#%s#%s" % (i, j)
            cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])

    # adding Cijh - products of L and R
    for i in range(M):
        for j in range(M):
            for h in range(M):
                name = "F#%s#%s#%s" % (i, j, h)
                cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])

    # adding Wij = Mij - Rij
    for i in range(M):
        for j in range(M):
            name = "W#%s#%s" % (i, j)
            cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])


    for i in range(M):
        for k in range(K):
            name = "H#%s#%s" % (i, k)
            cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])




    for i in range(M):
        # Ui <= Vi
        inds = ["U#%s" % i, "V#%s" % i]
        vals = [1, -1]
        cons = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [cons],\
            senses = ["L"],\
            rhs = [-0.01],\
            names = ['cons1-%s' % i]\
        )

        for j in range(M):
            if i < j:
                # A
                inds = ["A#%s#%s" % (i, j), "U#%s" % i, "U#%s" % j]
                vals = [-1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [-0.01],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )
                inds = ["U#%s" % i, "U#%s" % j, "A#%s#%s" % (i, j)]
                vals = [1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [1],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )

                # B
                inds = ["B#%s#%s" % (i, j), "U#%s" % j, "V#%s" % i]
                vals = [-1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [-0.01],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )
                inds = ["U#%s" % j, "V#%s" % i, "B#%s#%s" % (i, j)]
                vals = [1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [1],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )

                # C
                inds = ["C#%s#%s" % (i, j), "V#%s" % j, "U#%s" % i]
                vals = [-1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [-0.01],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )
                inds = ["V#%s" % j, "U#%s" % i, "C#%s#%s" % (i, j)]
                vals = [1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [1],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )

                # D
                inds = ["D#%s#%s" % (i, j), "V#%s" % i, "V#%s" % j]
                vals = [-1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [-0.01],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )
                inds = ["V#%s" % i, "V#%s" % j, "D#%s#%s" % (i, j)]
                vals = [1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [1],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )


                # A, B, C, D
                inds = ["A#%s#%s" % (i, j), "B#%s#%s" % (i, j), "D#%s#%s" % (i, j)]
                vals = [1, 1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [2],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )
                inds = ["A#%s#%s" % (i, j), "C#%s#%s" % (i, j), "D#%s#%s" % (i, j)]
                vals = [1, 1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["G"],\
                    rhs = [1],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )


                # Xij
                inds = ["X#%s#%s" % (i, j), "A#%s#%s" % (i, j),"D#%s#%s" % (i, j)]
                vals = [1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["G"],\
                    rhs = [0],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )

                inds = ["X#%s#%s" % (i, j), "A#%s#%s" % (i, j)]
                vals = [1, -1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [0],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )

                inds = ["X#%s#%s" % (i, j), "D#%s#%s" % (i, j)]
                vals = [1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [1],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )

                # Xji
                inds = ["X#%s#%s" % (j, i), "D#%s#%s" % (i, j), "A#%s#%s" % (i, j)]
                vals = [1, -1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["G"],\
                    rhs = [0],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )

                inds = ["X#%s#%s" % (j, i), "D#%s#%s" % (i, j)]
                vals = [1, -1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [0],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )

                inds = ["X#%s#%s" % (j, i), "A#%s#%s" % (i, j)]
                vals = [1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [1],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )

                inds = ["X#%s#%s" % (i, j), "X#%s#%s" % (j, i)]
                vals = [1, 1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [1],\
                    names = ['cons2-%s-%s' % (i, j)]\
                )




    # Xii must be set to 1
    for i in range(M):
        inds = ["X#%s#%s" % (i, i)]
        vals = [1]
        cons = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [cons],\
            senses = ["E"],\
            rhs = [1],\
            names = ['cons8.5-%s-%s' % (i, j)]\
        )


    # a cell must be attached to one of the mutation
    for k in range(K):
        inds = ["S#%s#%s" % (i, k) for i in range(M)]
        vals = [1 for i in range(M)]
        cons = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [cons],\
            senses = ["E"],\
            rhs = [1],\
            names = ['cons9-%s' % k]\
        )


    ################ Ensure number of leaves is at most number of cells
    for i in range(M):
        # Sum Xij + SUM Sij >= 2 for each i - > if it is a leaf - it must have a cell
        inds = ["S#%s#%s" % (i, k) for k in range(K)] + ["X#%s#%s" % (i, j) for j in range(M)]
        vals = [1 for k in range(K)] + [1 for j in range(M)]
        cons = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [cons],\
            senses = ["G"],\
            rhs = [2],\
            names = ['cons9-%s' % i]\
        )

    #########################################3



    #################### LOOOOSSSSSSSSSSS




    # main diagonal of L must be set to 0
    for i in range(M):
        inds = ["L#%s#%s" % (i, i)]
        vals = [1]
        cons = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [cons],\
            senses = ["E"],\
            rhs = [0],\
            names = ['cons10a-%s' % i]\
        )

    # no loss in the root
    for i in range(M):
        for j in range(M):
            if i == 0 or j == 0:
                inds = ["L#%s#%s" % (i, j)]
                vals = [1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["E"],\
                    rhs = [0],\
                    names = ['cons10b-%s-%s' % (i, j)]\
                )

    # all entries Lij must be not greater than Mij
    for i in range(M):
        for j in range(M):
            inds = ["L#%s#%s" % (i, j), "X#%s#%s" % (i, j)]
            vals = [1, -1]
            cons = cplex.SparsePair(ind=inds, val=vals)
            cpx.linear_constraints.add( \
                lin_expr = [cons],\
                senses = ["L"],\
                rhs = [0],\
                names = ['cons11-%s' % i]\
            )



    # in each column of R not more than two ones - allow only at most one loss for each cell
    for j in range(M):
        inds = ["R#%s#%s" % (i, j) for i in range(M)]
        vals = [1 for i in range(M)]
        cons = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [cons],\
            senses = ["L"],\
            rhs = [1],\
            names = ['cons11-%s' % k]\
        )



    # number of losses must be not greater than nLoss
    inds = ["L#%s#%s" % (i, j) for i in range(M) for j in range(M)]
    vals = [1 for i in range(M) for j in range(M)]
    cons = cplex.SparsePair(ind=inds, val=vals)
    cpx.linear_constraints.add( \
        lin_expr = [cons],\
        senses = ["E"],\
        rhs = [nLoss],\
        names = ['cons12']\
    )


    # add constraints for Cijh - products of Lih and Rhj ?????????????
    for i in range(M):
        for j in range(M):
            for h in range(M):
                inds = ["F#%s#%s#%s" % (i, j, h), "L#%s#%s" % (i, h)]
                vals = [1, -1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [0],\
                    names = ['cons13A-%s-%s-%s' % (i, j, k)]\
                )
                inds = ["F#%s#%s#%s" % (i, j, h), "X#%s#%s" % (h, j)]
                vals = [1, -1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [0],\
                    names = ['cons13B-%s-%s-%s' % (i, j, k)]\
                )
                inds = ["F#%s#%s#%s" % (i, j, h), "L#%s#%s" % (i, h), "X#%s#%s" % (h, j)]
                vals = [1, -1, -1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["G"],\
                    rhs = [-1],\
                    names = ['cons13C-%s-%s-%s' % (i, j, k)]\
                )


    # add constraint that Rij = sum Cijh
    for i in range(M):
        for j in range(M):
            inds = ["R#%s#%s" % (i, j)] +  ["F#%s#%s#%s" % (i, j, h) for h in range(M)]
            vals = [1] + [-1 for h in range(M)]
            cons = cplex.SparsePair(ind=inds, val=vals)
            cpx.linear_constraints.add( \
                lin_expr = [cons],\
                senses = ["E"],\
                rhs = [0],\
                names = ['cons14-%s' % i]\
            )


    # add equality that Wij = Mij - Rij
    for i in range(M):
        for j in range(M):
            inds = ["W#%s#%s" % (i, j), "X#%s#%s" % (i, j), "R#%s#%s" % (i, j)]
            vals = [1, -1, 1]
            cons = cplex.SparsePair(ind=inds, val=vals)
            cpx.linear_constraints.add( \
                lin_expr = [cons],\
                senses = ["E"],\
                rhs = [0],\
                names = ['cons15-%s-%s' % (i, j)]\
            )


    ##############################################



    # add constraints for Gikj
    for i in range(M):
        for k in range(K):
            for j in range(M):
                inds = ["G#%s#%s#%s" % (i, k, j), "W#%s#%s" % (i, j)]
                vals = [1, -1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [0],\
                    names = ['cons16A-%s-%s-%s' % (i, j, k)]\
                )
                inds = ["G#%s#%s#%s" % (i, k, j), "S#%s#%s" % (j, k)]
                vals = [1, -1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["L"],\
                    rhs = [0],\
                    names = ['cons16B-%s-%s-%s' % (i, j, k)]\
                )
                inds = ["G#%s#%s#%s" % (i, k, j), "W#%s#%s" % (i, j), "S#%s#%s" % (j, k)]
                vals = [1, -1, -1]
                cons = cplex.SparsePair(ind=inds, val=vals)
                cpx.linear_constraints.add( \
                    lin_expr = [cons],\
                    senses = ["G"],\
                    rhs = [-1],\
                    names = ['cons16C-%s-%s-%s' % (i, k, j)]\
                )

    # add constraint that Pik = sum Gikj
    for i in range(M):
        for k in range(K):
            inds = ["P#%s#%s" % (i, k)] +  ["G#%s#%s#%s" % (i, k, j) for j in range(M)]
            vals = [1] + [-1 for j in range(M)]
            cons = cplex.SparsePair(ind=inds, val=vals)
            cpx.linear_constraints.add( \
                lin_expr = [cons],\
                senses = ["E"],\
                rhs = [0],\
                names = ['cons17-%s-%s' % (i, k)]\
            )


    for i in range(M):
        inds = ["X#%s#%s" % (0, i)]
        vals = [1]
        cons = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [cons],\
            senses = ["E"],\
            rhs = [1],\
            names = ['cons18-%s' % i]\
        )


    for i in range(1, M):
        inds = ["X#%s#%s" % (i, 0)]
        vals = [1]
        cons = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [cons],\
            senses = ["E"],\
            rhs = [0],\
            names = ['cons18-%s' % i]\
        )

    """
    # xor objective
    for i in range(M):
        for k in range(K):
            inds = ["H#%s#%s" % (i, k), "P#%s#%s" % (i, k)]
            vals = [1, -1]
            cons = cplex.SparsePair(ind=inds, val=vals)
            cpx.linear_constraints.add( \
                lin_expr = [cons],\
                senses = ["L"],\
                rhs = [1 - matrix[i][k]],\
                names = ['cons17-%s-%s' % (i, k)]\
            )
            inds = ["H#%s#%s" % (i, k), "P#%s#%s" % (i, k)]
            vals = [1, -1]
            cons = cplex.SparsePair(ind=inds, val=vals)
            cpx.linear_constraints.add( \
                lin_expr = [cons],\
                senses = ["G"],\
                rhs = [matrix[i][k] - 1],\
                names = ['cons17-%s-%s' % (i, k)]\
            )
            inds = ["H#%s#%s" % (i, k), "P#%s#%s" % (i, k)]
            vals = [1, 1]
            cons = cplex.SparsePair(ind=inds, val=vals)
            cpx.linear_constraints.add( \
                lin_expr = [cons],\
                senses = ["G"],\
                rhs = [1 - matrix[i][k]],\
                names = ['cons17-%s-%s' % (i, k)]\
            )
            inds = ["H#%s#%s" % (i, k), "P#%s#%s" % (i, k)]
            vals = [1, 1]
            cons = cplex.SparsePair(ind=inds, val=vals)
            cpx.linear_constraints.add( \
                lin_expr = [cons],\
                senses = ["L"],\
                rhs = [1 + matrix[i][k]],\
                names = ['cons17-%s-%s' % (i, k)]\
            )

    # different objective
    for i in range(1, M):
        for k in range(K):
            cpx.objective.set_linear("H#%s#%s" % (i, k), 1)

    cpx.objective.set_sense(cpx.objective.sense.maximize)
    """


    # setting objective
    for i in range(1, M): # not counting the newly introduced root
        for k in range(K):
            if matrix[i][k] == 0:
                cpx.objective.set_linear("P#%s#%s" % (i, k), log((1 - BETA) / ALPHA))#     log((1 - ALPHA) / BETA))
            elif matrix[i][k] == 1:
                cpx.objective.set_linear("P#%s#%s" % (i, k), -log((1 - ALPHA) / BETA))#        -log((1 - BETA) / ALPHA))

    cpx.objective.set_sense(cpx.objective.sense.minimize)

    cpx.set_problem_type(cpx.problem_type.MILP)
    cpx.write("program.txt", filetype="lp")
    cpx.solve()
    cpx.populate_solution_pool()
    print cpx.solution.pool.get_num()

    print "OBJECTIVE FUNCTION VALUE:", float(cpx.solution.get_objective_value())# + sum(map(sum, matrix)))

    solutions = []
    mutationTreeMatrices = []
    lossMatrices = []
    attachmentsAll = []
    objVals = []

    
    for ii in range(cpx.solution.pool.get_num()):

        objVals.append(cpx.solution.pool.get_objective_value(ii))


        attachments = []
        for i in range(M):
            attRow = []
            for k in range(K):
                attRow.append(int(round(cpx.solution.pool.get_values(ii, "S#%s#%s" % (i, k)))))
            attachments.append(attRow)
        #print "ATTACHMENTS MATRIX"
        #pprint(attachments)



        losses = 0
        for i in range(M):
            for j in range(M):
                losses += int(round(cpx.solution.pool.get_values(ii, "L#%s#%s" % (i, j))))
        #print "ACTUAL NR LOSSES:", losses


        lossMatrix = []
        for i in range(M):
            lossRow = []
            for j in range(M):
                lossRow.append(int(round(cpx.solution.pool.get_values(ii, "L#%s#%s" % (i, j)))))
            lossMatrix.append(lossRow)
        #print "LOSS MATRIX"
        #pprint(lossMatrix)


        rMatrix = []
        for i in range(M):
            rRow = []
            for j in range(M):
                rRow.append(int(round(cpx.solution.pool.get_values(ii, "R#%s#%s" % (i, j)))))
            rMatrix.append(rRow)
        #print "R MATRIX"
        #pprint(rMatrix)



        mutationTreeMatrix = []
        for i in range(M):
            mutRow = []
            for j in range(M):
                mutRow.append(int(round(cpx.solution.pool.get_values(ii, "X#%s#%s" % (i, j)))))
            mutationTreeMatrix.append(mutRow)
        #print "MUTATION MATRIX"
        #pprint(mutationTreeMatrix)


        solution = []
        for i in range(1, M):
            solX = []
            for k in range(K):
                solX.append(int(round(cpx.solution.pool.get_values(ii, "P#%s#%s" % (i, k)))))
            solution.append(solX)
        print "Inferred matrix"
        pprint(solution)

        solutions.append(solution)
        mutationTreeMatrices.append(mutationTreeMatrix)
        lossMatrices.append(lossMatrix)
        attachmentsAll.append(attachments)


    print objVals

    minObjValue = min(objVals)
    goodIndices = [x for x, y in enumerate(objVals) if y == minObjValue]

    good_solutions = []
    count = 0
    for sol, mut, los, att in zip(solutions, mutationTreeMatrices, lossMatrices, attachmentsAll):
        if count in goodIndices:
            good_solutions.append((sol, mut, los, att))
        count += 1
    good_solutions2 = [good_solutions[0]]
    for i in range(1, len(good_solutions)):
        new = True
        for x in good_solutions2:
            if x == good_solutions[i]:
                new = False
                break
        if new:
            good_solutions2.append(good_solutions[i])
    print len(good_solutions2), "OPTIMAL SOLUTIONS"
    solutions = []
    mutationsTreeMatrices = []
    lossMatrices = []
    attachmentsAll = []
    for sol, mut, los, att in good_solutions2:
        solutions.append(sol)
        mutationTreeMatrices.append(mut)
        lossMatrices.append(los)
        attachmentsAll.append(att)
    return solutions, mutationTreeMatrices, lossMatrices, attachmentsAll


matfile = sys.argv[1]
nLoss = int(sys.argv[2])

with open(matfile) as f:
    matdata = f.readlines()

matdata = map(lambda x: x.strip().split(), matdata)

mutation_matrix = []
for line in matdata:
    mutation_matrix.append(map(int, line))


solutions, matrices, lossesAll, attachmentsAll = single_cell_phylogeny(mutation_matrix, nLoss)


index = 0
for solution, matrix, losses, attachments in zip(solutions, matrices, lossesAll, attachmentsAll):



    ancestry_dict = {}
    for i in range(len(matrix)):
        ancestry_dict[i] = []
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if i != j:
                if matrix[i][j] == 1:
                    ancestry_dict[i].append(j)

    print ancestry_dict, "ANCESTRY"


    parents = {}

    while True:
        leaves = [node for node in ancestry_dict if ancestry_dict[node] == []]
        if leaves == [0]:
            break
        for leaf in leaves:
            del ancestry_dict[leaf]
        # determine who is the parent
        possible_parents = []
        for leaf in leaves:
            possible_parents = [(node, len(ancestry_dict[node])) for node in ancestry_dict if leaf in ancestry_dict[node]]
            parent = min(possible_parents, key=lambda z: z[1])[0]
            parents[leaf] = parent
        for node in ancestry_dict:
            #lifs = [x for x in ancestry_dict[node] if x in leaves]
            #for l in lifs:
            #    parents[l] = node
            ancestry_dict[node] = [x for x in ancestry_dict[node] if x not in leaves]



    A = PG.AGraph(directed=True, strict=True)

    for u, v in parents.items():
        A.add_edge("M" + str(v), "M" + str(u))


    ######### uncomment
    for i in range(len(attachments[0])):
        for j in range(len(attachments)):
            if attachments[j][i] == 1:
                A.add_edge("M" + str(j), "S" + str(i + 1))



    # add losses
    loss_dict = {}
    for i in range(len(losses)):
        for j in range(len(losses[0])):
            if losses[i][j] == 1:
                if j not in loss_dict:
                    loss_dict[j] = []
                loss_dict[j].append(i)


    loss_freq_dict = {}
    for node1, nodes2 in loss_dict.items():
        for n in nodes2:
            if n not in loss_freq_dict:
                loss_freq_dict[n] = 0
            loss_freq_dict[n] += 1




    for node in loss_dict:
        parent =  A.in_neighbors("M" + str(node))[0]
        losses_to_insert = loss_dict[node]
        A.remove_edge(parent, "M" + str(node))
        for i in range(len(losses_to_insert)):
            if i == 0:
                if loss_freq_dict[losses_to_insert[i]] > 1:
                    suffix = "." + str(loss_freq_dict[losses_to_insert[i]])
                    loss_freq_dict[losses_to_insert[i]] -= 1
                else:
                    suffix = ""
                A.add_edge(parent, "L" + str(losses_to_insert[i]) + suffix)
                parent = "L" + str(losses_to_insert[i]) + suffix
            else:
                if loss_freq_dict[losses_to_insert[i]] > 1:
                    suffix = "." + str(loss_freq_dict[losses_to_insert[i]])
                    loss_freq_dict[losses_to_insert[i]] -= 1
                else:
                    suffix = ""
                A.add_edge(parent, "L" + str(losses_to_insert[i]) + suffix)
                parent = "L" + str(losses_to_insert[i]) + suffix
        A.add_edge("L" + str(losses_to_insert[-1]) + suffix, "M" + str(node))

    A.node_attr['style']='filled'
    for node in A.nodes():
        if node == "M0":
            node_attr = A.get_node(node)
            node_attr.attr['fillcolor']="yellow"
        else:
            if node.startswith("M"):
                node_attr = A.get_node(node)
                node_attr.attr['fillcolor']="#CCCCFF"
            elif node.startswith("L"):
                node_attr = A.get_node(node)
                node_attr.attr['fillcolor']="red"
            elif node.startswith("S"):
                node_attr = A.get_node(node)
                node_attr.attr['fillcolor']="green"


    where_to_save = sys.argv[1].split(".")[0]
    A.write(where_to_save + ".igor.solution_%s.dot" % index)

    A.layout(prog='dot')

    A.draw(where_to_save + ".igor.solution_%s.jpg" % index)


    ########################################## EVALUATE THE SOLUTION ###################

    if len(sys.argv) == 5: # here we have to supply the ground truth
        groundMfile = sys.argv[3]
        with open(groundMfile) as f:
            groundMfile = f.readlines()
        groundMfile = map(lambda x: x.strip().split(), groundMfile)
        groundM = []
        for line in groundMfile:
            groundM.append(map(int, line))

        # here we work with mutation_matrix and solution matrices
        # evaluate false positive and false negative
        trueMijOnes = 0
        mutMijOnes = 0
        solMijOnes = 0
        for mline, sline in zip(groundM, solution):
            for m, s in zip(mline, sline):
                if m == 1:
                    mutMijOnes += 1
                if s == 1:
                    solMijOnes += 1
                if m == 1 and s == 1:
                    trueMijOnes += 1
        sensitivity = round(trueMijOnes * 1.0 / mutMijOnes, 2)
        ppv = round(trueMijOnes * 1.0 / solMijOnes, 2)
        print "SENSITIVITY", sensitivity
        print "PPV", ppv



        alphas = 0
        betas = 0
        solution_mismatches = 0
        print(solution)
        print(groundM)
        print(mutation_matrix)
        for u, v in zip(solution, groundM):
            for uu, vv in zip(u, v):
                if uu != vv:
                    if uu == 1:
                        alphas += 1
                    else:
                        betas += 1
                    solution_mismatches += 1

        print "Solution %s: solution mismatches in ancestry matrix: %s, alphas: %s, betas: %s" % (index, solution_mismatches, alphas, betas)





        # evaluate mutations from mutation tree
        #pprint(matrix)
        with open(sys.argv[4]) as f:
            ancestFile = f.readlines()
        ancestFile = map(lambda x: x.strip().split(), ancestFile)
        cellAncestors = ancestFile[len(matrix) - 1:]
        ancestFile = ancestFile[:len(matrix) - 1]
        trueAncestDict = {}
        for i, line in enumerate(ancestFile):
            trueAncestDict[i + 1] = set(map(int, line))
        inferredAncestDict = {}
        for i in range(1, len(matrix[0])):
            anc = [matrix[j][i] for j in range(len(matrix))]
            #print anc
            ances = set()
            for u, v in enumerate(anc):
                if u != i: # do not count node as its own ancestor
                    if v == 1:
                        ances.add(u)
            if len(ances) > 1:
                ances.remove(0)
            inferredAncestDict[i] = ances

        nTruePar = 0
        nEstimPar = 0
        truePos = 0
        for i in range(1, len(matrix[0])):
            #print trueAncestDict[i], inferredAncestDict[i]
            nTruePar += len(trueAncestDict[i])
            nEstimPar += len(inferredAncestDict[i])
            truePos += len(trueAncestDict[i] & inferredAncestDict[i])
        print "SENSITIVITY ANCESTORS", round(truePos * 1.0 / nTruePar, 2)
        print "PPV ANCESTORS", round(truePos * 1.0 / nEstimPar, 2)
        #pprint(solution)


        trueCellAncestDict = {}
        for i, line in enumerate(cellAncestors):
            trueCellAncestDict[i] = set(map(int, line))

        inferredCellAncestDict = {}
        # for each cell check its ancestors
        for i in range(len(solution[0])):
            anc = [solution[j][i] for j in range(len(solution))]
            ances = set()
            for u, v in enumerate(anc):
                if v == 1:
                    ances.add(u + 1)
            if len(ances) == 0:
                ances.add(0)
            inferredCellAncestDict[i] = ances
        
        for i in range(len(solution[0])):
            #print trueCellAncestDict[i], inferredCellAncestDict[i]
            nTruePar += len(trueCellAncestDict[i])
            nEstimPar += len(inferredCellAncestDict[i])
            truePos += len(trueCellAncestDict[i] & inferredCellAncestDict[i])
        print "SENSITIVITY ALL", round(truePos * 1.0 / nTruePar, 2)
        print "PPV ANCESTORS", round(truePos * 1.0 / nEstimPar, 2)
        print "---------------------------------------------------------"
    index += 1


