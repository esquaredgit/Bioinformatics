import math 
lis1 = [10, 8, 10]
lis2 = [4, 8.5, 3]
zipped = list(zip(lis1, lis2))
under = sum([(t[1]-t[0])**2 for t in zipped])
dist = math.sqrt(under)

print(dist)