# -*- coding: utf-8 -*-
"""B363 Project.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1Rg0dlcRILtUmkNuqx3aZspnjdXHlR9-s

Lloyd Algorithm
"""

import numpy as np
from collections import defaultdict

open_file = open('/content/cleanExpMatrix.txt', 'r')
input_lines = open_file.read().splitlines()
k, m = [int(x) for x in input_lines[0].split()]
data = [[float(x) for x in line.split()] for line in input_lines[1:]]


def euclidean(x, y):
    point1 = np.array(x)
    point2 = np.array(y)
    dist = np.linalg.norm(point1 - point2)
    return dist


def center_p(point, centers):
    min_dist = float("Inf")
    for x in centers:
        current = euclidean(x, point)
        if current < min_dist:
            min_dist = current
            closest = x
    return closest


def c_m(cluster):
    m = len(cluster[0])
    center = [0] * m
    for point in cluster:
        for i in range(m):
            center[i] += point[i]
    center = [x / len(cluster) for x in center]
    return center


def lloyd_k_means(data, k):
    centers = data[:k]

    while True:
        cluster_rain = defaultdict(list)
        for point in data:
            center = center_p(point, centers)
            cluster_rain[tuple(center)].append(point)

        new_centers = [[]] * k
        for i in range(k):
            new_centers[i] = c_m(cluster_rain[tuple(centers[i])])

        if new_centers == centers:
            break
        centers = new_centers[:]
    return centers


centers = lloyd_k_means(data, k)
print("k =",k)
for center in centers:
    print(" ".join(map(str, center)))

"""For plotting line graphs"""

#print(centers)
centers = np.array(centers)
print(centers[:,0])


import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np


x = np.linspace(0, 2 * np.pi, 400)
plt.plot(x, np.sin(x))
plt.title("A Sine Curve")
plt.xlabel("x")
plt.ylabel("sin(x)");

"""For multiple plots"""

fig, axs = plt.subplots(2, 3)
x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)
axs[0, 0].plot(x, y)
axs[0, 0].set_title('Axis [0, 0]')
axs[0, 1].plot(x, y, 'tab:orange')
axs[0, 1].set_title('Axis [0, 1]')
axs[1, 0].plot(x, -y, 'tab:green')
axs[1, 0].set_title('Axis [1, 0]')
axs[1, 1].plot(x, -y, 'tab:red')
axs[1, 1].set_title('Axis [1, 1]')

for ax in axs.flat:
    ax.set(xlabel='x-label', ylabel='y-label')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
