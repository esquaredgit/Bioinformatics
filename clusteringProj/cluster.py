# CSCI-B363 Final Project: Topic 2.1; Ethan Eldridge, Nick Frasco, Matt Kunin
# Data taken from DiRisi and colleagues: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28

from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plot
import numpy as np
import math
import sys

# Begin by cleaning up the expression data (removing any rows with missing values):

with open("GSE28_series_matrix.txt", "r") as file:
    exp = file.readlines()
    file.close()

col_names = exp[0]
col_names = col_names.replace("\n", "")
col_names = col_names.replace("\"", "")
col_names = col_names.split("\t")
del exp[0]

for i in range(0,len(exp)):
    exp[i] = exp[i].replace("\n", "")
    exp[i] = exp[i].split("\t")
    del exp[i][0]

# converting data to strings and removing rows with incomplete data
exp = [[float(string) for string in row] for row in exp if not ("" in row)]
"""
This is extra script for copying over the clean exp data:
exp = [[string for string in row] for row in exp if not ("" in row)]
with open("cleanExpMatrix.txt", "w") as file:
    first = "4 7\n"
    file.write(first)
    for row in exp:
        line = ""
        for point in row:
            line = line + str(point) + " "
        line = line + "\n"
        file.write(line)
"""

# Creating Distance Matrix using Euclidean distance (ETHAN)
distMat = np.zeros((len(exp), len(exp)))
counter = 0

zCounter = 0
continued = False
# with open("distanceMatrix.txt", "w") as file:
for row in range(0, len(distMat)):
    print("Processing distance matrix row: " + str(counter))
    # line = ""
    for col in range(0, len(distMat[0])):
        if (row == col): 
            distMat[row][col] = 0.0
            # line = line + "0.0 "
            zCounter += 1
            continued = True
            continue
        if (row > col):
            distMat[row][col] = distMat[col][row]
            continued = True
            continue
        if (row > col) or (row == col) and (continued == True): 
            print("WTFFF")
            sys.exit()
        continued = False
        zipped = list(zip(exp[row], exp[col]))
        under = sum([(t[1]-t[0])**2 for t in zipped])
        dist = math.sqrt(under)
        distMat[row][col] = dist
        # line = line + str(dist) + " "
    # file.write(line + "\n")
    counter += 1
print(distMat)

# Running hClustering (ETHAN)


# Creating Dendrogram visualization (ETHAN)
info = linkage(distMat, "average")
figure = plot.figure(figsize=(25, 10))
graph = dendrogram(info)
plot.show()
print(str(set(graph["color_list"])))

# Running Lloyd algorithm on Distance Matrix using k (maybe k = 4) from hierarchical clustering (NICK)
K = len(set(graph["color_list"]))

# Visualizing each cluster (line of best fit for each cluster AND all vectors in each cluster) (NICK)


# Isolating each cluster with activity AFTER the shift to use in the CSRE portion (create some visualization with these) (NICK)


# Gibbs sampling or something to find the CSRE for only those clusters (MATT)


# Create a simple visualization with each cluster and its consensus string (MATT)
