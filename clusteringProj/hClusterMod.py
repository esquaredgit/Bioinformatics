# CSCI-B363 Final Project: Topic 2.1; Ethan Eldridge, Nick Frasco, Matt Kunin
# Data taken from DiRisi and colleagues: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28

from node import Node
import copy

# The actual stuff:

def hCluster(D, n):
    clusters = [Node([i], None, None) for i in range(1, n+1)]
    T = list()
    while len(clusters) > 1: 
        for r in D:
            print(r)
    #   find the two closest clusters Ci and Cj 
        minI = 0
        minJ = 1
        for i in range(0, len(D)):
            for j in range(0, len(D[0])):
                if i == j:
                    continue
                if D[i][j] < D[minI][minJ]: 
                    minI = i
                    minJ = j
        print("mins: " + str(minI) + " " + str(minJ))
    #   merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
    #   add a new node labeled by cluster Cnew to T
    #   connect node Cnew to Ci and Cj by directed edges
        iVal = copy.deepcopy(clusters[minI].getVal())
        jVal = copy.deepcopy(clusters[minJ].getVal())
        newVal = [i for i in iVal]
        newVal.extend([j for j in jVal])
        cNew = Node(newVal, clusters[minI], clusters[minJ])
        print(str(cNew))
        T.append(cNew)
    #   remove the rows and columns of D corresponding to Ci and Cj
        D = [D[i] for i in range(0, len(D)) if ((i != minI) and (i != minJ))]
        for i in range(0, len(D)):
            D[i] = [D[i][col] for col in range(0, len(D[i])) if ((col != minI) and (col != minJ))]
    #   remove Ci and Cj from Clusters
        clusters = [clusters[c] for c in range(0, len(clusters)) if ((c != minI) and (c != minJ))]
    #   add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
        newRow = []
        for c in clusters:
            num = 0
            for i in cNew.getVal():
                for j in c.getVal():
                    num += data[i-1][j-1]
            denom = len(cNew.getVal()) * len(c.getVal())
            dAvg = num/denom
            newRow.append(dAvg)
        newRow.append(0)
        print("newRow: " + str(newRow))
        # adding row
        D.append(newRow)
        for r in D: print()
        # adding col
        for i in range(0, len(D)-1):
            D[i].append(newRow[i])
    #   add Cnew to Clusters 
        clusters.append(cNew)
#   assign root in T as a node with no incoming edges
    T.append(clusters[0])
    return T

with open("answer.txt", "w") as file:
    for c in hCluster(data, N):
        line = ""
        for i in range(0, len(c.getVal())):
            if i == len(c.getVal())-1: 
                line += str(c.getVal()[i])
                continue
            line += str(c.getVal()[i]) + " "
        file.write(line + "\n")
