# Ethan Eldridge
import math
import random
import matplotlib.pyplot as plt
import sys

# Formatting 
with open("cleanExpMatrixAv.txt", "r") as file:
    data = file.readlines()

for i in range(0,len(data)):
    data[i] = data[i].replace("\n", "")
    data[i] = data[i].split()
    data[i] = [float(string) for string in data[i]]

K, M = int(data[0][0]), int(data[0][1])

del data[0]

for i in range(1, len(data)+1):
    data[i-1] = [i, data[i-1]]

print("k: " + str(K))
print("m: " + str(M))
for c in range(0, 10):
    print(data[c])
print()


count = 0

def lloyd(data, k, m):
    centers = data[0:k]
    clusters = [[c] for c in centers]
    while True:
        global count
        count += 1
        # assigning each data point to closest cluster:
        for point in data:
            distToCenters = []
            for c in centers:
                underRad = 0
                for i in range(0, m):
                    underRad += (point[1][i]-c[1][i])**2
                dist = math.sqrt(underRad)
                distToCenters.append(dist)
            assignedCenter = distToCenters.index(min(distToCenters))
            clusters[assignedCenter].append(point)
        # finding new cluster centers
        newCenters = []
        for c in clusters:
            grav = []
            for i in range(0, m):
                grav.append(sum([p[1][i] for p in c])/len(c))
            newCenters.append([0, grav])
        if [n[1] for n in newCenters] == [c[1] for c in centers]:
            break
        # assigning new centers to the clusters
        centers = newCenters
        clusters = [[c] for c in centers]
    return centers, clusters

output = lloyd(data, K, M)
rawCens = output[0]
rawClus = output[1]

# with open("clusterAnswers.txt", "w") as file:
#     file.write("k t 100")
#     for c in range(0, len(rawClus)):



"""
STILL NEED TO SORT OUT THE INDEXING STUFF BUT OTHER THAN THAT IT'S PRETTY MUCH GOOD
"""

cens = [c[1] for c in rawCens]
print(cens)
for c in rawClus:
    print(c[0][1])

for p in range(0, len(cens)):
    plt.plot([i for i in range(1, 8)], cens[p], label="Cluster " + str(p+1))
plt.yticks([-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3])
plt.legend(loc="lower left")
plt.show()

# print("Centers: " + str(cens))
# print()
# print("Clusters:")

clus = [[p[1] for p in c] for c in rawClus]

print("# of clusters: " + str(len(clus)))

# for p in range(0, 270):
#     plt.plot([i for i in range(1, 8)], clus[0][p])
# plt.show()

# for p in clus[1]:
#     plt.plot([i for i in range(1, 8)], p)
# plt.show()

# for p in clus[2]:
#     plt.plot([i for i in range(1, 8)], p)
# plt.show()


with open("clusterAnswers.txt", "w") as file:
    # Need indexes for clusters 2, 3, 5
    file.write("k " + str(len(clus[1])) + " 100\n")
    for i in range(0, len(rawClus[1])):
        file.write(str(rawClus[1][i][0])+"\n")

    file.write("k " + str(len(clus[2])) + " 100\n")
    for i in range(0, len(rawClus[2])):
        file.write(str(rawClus[2][i][0])+"\n")

    file.write("k " + str(len(clus[4])) + " 100\n")
    for i in range(0, len(rawClus[4])):
        file.write(str(rawClus[4][i][0])+"\n")

# with open("answer.txt", "w") as file:
#     for p in lloyd(data, K, M):
#         line = ""
#         for i in range(0, len(p)):
#             if i == len(p)-1: 
#                 line += str(round(p[i],3))
#                 continue
#             line += str(round(p[i],3)) + " "
#         file.write(line + "\n")