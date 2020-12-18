# CSCI-B363 Final Project: Topic 2.1; Ethan Eldridge, Nick Frasco, Matt Kunin
# Data taken from DiRisi and colleagues: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28

from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plot
import numpy as np
import math
import sys

from sklearn.cluster import KMeans 
from sklearn import metrics 
from scipy.spatial.distance import cdist 
import numpy as np 
import matplotlib.pyplot as plt 

# Begin by cleaning up the expression data (removing any rows with missing values):

with open("GSE28_series_matrix.txt", "r") as file:
    exp = file.readlines()
    file.close()


col_names = exp[0]
col_names = col_names.replace("\n", "")
col_names = col_names.replace("\"", "")
col_names = col_names.split("\t")
col_names = col_names[1:]
del exp[0]

print(col_names)

for i in range(0,len(exp)):
    exp[i] = exp[i].replace("\n", "")
    exp[i] = exp[i].split("\t")
    del exp[i][0]


# Lets find column averages to substitute for each empty val
colAvs = np.zeros(len(exp[0]))
for col in range(0, len(exp[0])):
    colSum = 0
    counter = 0
    for row in range(0, len(exp)):
        if exp[row][col] == "": continue
        counter += 1
        exp[row][col] = float(exp[row][col])
        colSum += exp[row][col]
    colAvs[col] = colSum/counter

# Substituting column average for blank spaces:
for row in range(0, len(exp)):
    for col in range(0, len(exp[0])):
        if exp[row][col] == "":
            exp[row][col] = colAvs[col]
print()
print("Averages: " + str(colAvs))
print()
print("First ten:" + str(exp[0:10]))
print()


col_names = np.asarray(col_names)
exp = np.asarray(exp)
print(exp[0:10])
print()

# switching columns to match time ordering in dataset:
exp[:, [1,0]] = exp[:, [0,1]]
col_names[1], col_names[0] = col_names[0], col_names[1]
colAvs[1], colAvs[0] = colAvs[0], colAvs[1]

exp[: , [2,3]] = exp[:, [3,2]]
col_names[2], col_names[3] = col_names[3], col_names[2]
colAvs[2], colAvs[3] = colAvs[3], colAvs[2]

exp[: , [1,2]] = exp[:, [2,1]]
col_names[1], col_names[2] = col_names[2], col_names[1]
colAvs[1], colAvs[2] = colAvs[2], colAvs[1]

exp[: , [3,4]] = exp[:, [4,3]]
col_names[3], col_names[4] = col_names[4], col_names[3]
colAvs[3], colAvs[4] = colAvs[4], colAvs[3]

exp[: , [2,3]] = exp[:, [3,2]]
col_names[2], col_names[3] = col_names[3], col_names[2]
colAvs[2], colAvs[3] = colAvs[3], colAvs[2]

exp[: , [4,5]] = exp[:, [5,4]]
col_names[4], col_names[5] = col_names[5], col_names[4]
colAvs[4], colAvs[5] = colAvs[5], colAvs[4]

exp[: , [3,4]] = exp[:, [4,3]]
col_names[3], col_names[4] = col_names[4], col_names[3]
colAvs[3], colAvs[4] = colAvs[4], colAvs[3]

print(col_names)
print(exp[0:10])


# converting data to strings and removing rows with incomplete data
# exp = [[float(string) for string in row] for row in exp if not ("" in row)]

# This is extra script for copying over the clean exp data:
"""
exp = [[string for string in row] for row in exp if not ("" in row)]
with open("cleanExpMatrixRm.txt", "w") as file:
    first = "4 7\n"
    file.write(first)
    for row in exp:
        line = ""
        for point in row:
            line = line + str(point) + " "
        line = line + "\n"
        file.write(line)


exp = [[string for string in row] for row in exp]
with open("cleanExpMatrixAv.txt", "w") as file:
    first = "5 7\n"
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
figure = plot.figure()
graph = dendrogram(info, truncate_mode="lastp")
plot.show()
print(str(set(graph["color_list"])))

# Running Lloyd algorithm on Distance Matrix using k (maybe k = 4) from hierarchical clustering (ETHAN)
distortions = [] 
inertias = [] 
mapping1 = {} 
mapping2 = {} 
K = range(1,10) 
  
for k in K: 
    #Building and fitting the model 
    kmeanModel = KMeans(n_clusters=k).fit(distMat) 
    kmeanModel.fit(distMat)     
      
    distortions.append(sum(np.min(cdist(distMat, kmeanModel.cluster_centers_, 
                      'euclidean'),axis=1)) / distMat.shape[0]) 
    inertias.append(kmeanModel.inertia_) 
  
    mapping1[k] = sum(np.min(cdist(distMat, kmeanModel.cluster_centers_, 
                 'euclidean'),axis=1)) / distMat.shape[0] 
    mapping2[k] = kmeanModel.inertia_ 

for key,val in mapping1.items(): 
    print(str(key)+' : '+str(val)) 

plt.plot(K, distortions, 'bx-') 
plt.xlabel('Values of K') 
plt.ylabel('Distortion') 
plt.title('The Elbow Method using Distortion') 
plt.show() 

for key,val in mapping2.items(): 
    print(str(key)+' : '+str(val)) 

plt.plot(K, inertias, 'bx-') 
plt.xlabel('Values of K') 
plt.ylabel('Inertia') 
plt.title('The Elbow Method using Inertia') 
plt.show() 

# Visualizing each cluster (line of best fit for each cluster AND all vectors in each cluster) (NICK)


# Isolating each cluster with activity AFTER the shift to use in the CSRE portion (create some visualization with these) (NICK)


# Gibbs sampling or something to find the CSRE for only those clusters (MATT)


# Create a simple visualization with each cluster and its consensus string (MATT)
