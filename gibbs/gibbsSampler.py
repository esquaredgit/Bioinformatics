# References/Citations:
# - https://machinelearningmastery.com/how-to-generate-random-numbers-in-python/
# - https://stackoverflow.com/questions/33359740/random-number-between-0-and-1-in-python
# - https://www.w3schools.com/python/python_file_open.asp
# - https://mail.python.org/pipermail/python-list/2001-February/105300.html

import random
from random import randint

def Random(t):
  return randint(0, t-1)

def GetRandomKMer(dnaSeq, profile, k):
  kmersInSequence = []
  for i in range(len(dnaSeq)):
    kmer = dnaSeq[i:i+k]
    if len(kmer) == k:
      kmersInSequence.append(kmer)
  
  kmerProbabilities = []

  for kmer in kmersInSequence:
    probabilityForKmer = 1
    i = 0
    for letter in kmer:
      row = -1
      if letter == 'A':
        row = 0
      elif letter == 'C':
        row = 1
      elif letter == 'G':
        row = 2
      elif letter == 'T':
        row = 3
      probabilityForKmer *= profile[row][i]
      i += 1
    
    kmerProbabilities.append(probabilityForKmer)

    # normalize for monte carlo method
    sum = 0
    for prob in kmerProbabilities:
      sum += prob

    for x in range(len(kmerProbabilities)):
      kmerProbabilities[x] /= sum

    for x in range(1, len(kmerProbabilities)):
      kmerProbabilities[x] += kmerProbabilities[x-1]

    rand = random.uniform(0, 1)
    
    z = 0
    for kmerProbability in kmerProbabilities:
      if rand < kmerProbability:
        return kmersInSequence[z]
      z += 1
    return None


def buildConsensusString(Motifs):
  consensusString = ""

  numOfCols = len(Motifs[0])

  for i in range(numOfCols):
    currentCol = []

    for motif in Motifs:
      if motif == None:
        continue
      currentCol.append(motif[i])

    letterFreqForCurCol = {}
    for letter in currentCol:
      if letter in letterFreqForCurCol:
        letterFreqForCurCol[letter] += 1
      else:
        letterFreqForCurCol[letter] = 0

    mostFreqLetterInCurrentCol = ""
    highestFreqInCurrentCol = 0

    for letter in letterFreqForCurCol:
      letterFreq = letterFreqForCurCol[letter]

      if letterFreq > highestFreqInCurrentCol:
        highestFreqInCurrentCol = letterFreq
        mostFreqLetterInCurrentCol = letter

    consensusString += mostFreqLetterInCurrentCol

  return consensusString

"""
def buildConsensusString(Motifs):
  consensusString = ""

  countMap = {}
  for i in range(len(Motifs[0])):
    countMap[i] = [0, 0, 0, 0]

  for motif in Motifs:
    for letterPos in range(len(Motifs[0])):
        if motif[letterPos] == 'A':
          countMap[letterPos][0] += 1
        if motif[letterPos] == 'T':
          countMap[letterPos][1] += 1
        if motif[letterPos] == 'C':
          countMap[letterPos][2] += 1
        if motif[letterPos] == 'G':
          countMap[letterPos][3] += 1

  for pos in countMap.values():
    letter = pos.index(max(pos))
    
    if letter == 0:
      consensusString += 'A'
    if letter == 1:
      consensusString += 'T'
    if letter == 2:
      consensusString += 'C'
    if letter == 3:
      consensusString += 'G'

  return consensusString
"""

"""
Motifs = []
Motifs.append('TCGGGGGTTTTT')
Motifs.append('CCGGTGACTTAC')
Motifs.append('ACGGGGATTTTC')
Motifs.append('TTGGGGACTTTT')
Motifs.append('AAGGGGACTTCC')
Motifs.append('TTGGGGACTTCC')
Motifs.append('TCGGGGATTCAT')
Motifs.append('TCGGGGATTCCT')
Motifs.append('TAGGGGAACTAC')
Motifs.append('TCGGGTATAACC')

print(buildConsensusString(Motifs))
"""

def Score(Motifs):
  # build consensus string
  # compute hamming distance between each Motif with the consensus string
  # add all hamming distances added together
  # hamming distances is number of mismatches
  consensusString = buildConsensusString(Motifs)

  score = 0
  for motif in Motifs:
    numberOfMismatches = 0
    i = 0
    for letter in motif:
      if letter != consensusString[i]:
        numberOfMismatches += 1
      i += 1
    score += numberOfMismatches
  return score

def getMotifsFrom(DNA, k):
  Motifs = []
  for sequence in DNA:
    kmersInSequence = []
    for i in range(len(sequence)):
      kmer = sequence[i:i+k]
      if len(kmer) == k:
        kmersInSequence.append(kmer)
    # pick random kmer from sequence array and add it to Motifs
    if len(kmersInSequence) > 0:
      Motifs.append(kmersInSequence[randint(0, len(kmersInSequence)-1)])
  return Motifs

def count(myVec, c):
  count = 0
  for i in range(len(myVec)):
    if myVec[i] == c:
      count += 1
  return count

def calculateCountMatrix(Motifs):
  numOfCols = len(Motifs[0])

  CountMatrix = []
  rowA = []
  rowC = []
  rowG = []
  rowT = []
  for i in range(numOfCols):
    rowA.append(0)
    rowC.append(0)
    rowG.append(0)
    rowT.append(0)
  CountMatrix.append(rowA)
  CountMatrix.append(rowC)
  CountMatrix.append(rowG)
  CountMatrix.append(rowT)

  for i in range(numOfCols):
    currentCol = []

    for motif in Motifs:
      currentCol.append(motif[i])

    numberOfA = count(currentCol, 'A')
    numberOfC = count(currentCol, 'C')
    numberOfG = count(currentCol, 'G')
    numberOfT = count(currentCol, 'T')

    CountMatrix[0][i] = numberOfA
    CountMatrix[1][i] = numberOfC
    CountMatrix[2][i] = numberOfG
    CountMatrix[3][i] = numberOfT

  return CountMatrix

def sumNums(myVec):
  count = 0
  for i in range(len(myVec)):
    count += myVec[i]
  return count

def calculateProfileMatrix(CountMatrix):
  numOfCols = len(CountMatrix[0])

  Profile = []

  rowA = []
  rowC = []
  rowG = []
  rowT = []
  for i in range(numOfCols):
    rowA.append(0.0)
    rowC.append(0.0)
    rowG.append(0.0)
    rowT.append(0.0)
  Profile.append(rowA)
  Profile.append(rowC)
  Profile.append(rowG)
  Profile.append(rowT)

  denominator = 0.0
  for i in range(len(CountMatrix)):
    denominator += CountMatrix[i][0]

  for i in range(len(CountMatrix)):
    for j in range(len(CountMatrix[0])):
      Profile[i][j] = CountMatrix[i][j] / denominator

  return Profile

"""
Motifs = []
Motifs.append('TAAC')
Motifs.append('GTCT')
Motifs.append('ACTA')
Motifs.append('AGGT')
CountMatrix = calculateCountMatrix(Motifs)
for n in range(len(CountMatrix)):
  for m in range(len(CountMatrix[0])):
    CountMatrix[n][m] += 1
Profile = calculateProfileMatrix(CountMatrix)
print(Profile)
"""

def inBounds(i, DNA):
  return i >= 0 and i <= len(DNA)-1

def GibbsSampler(DNA, k, t, N):
  Motifs = getMotifsFrom(DNA, k)
  BestMotifs = Motifs
  for j in range(N):
    i = Random(t) # 1 through t (between 0 and t-1 in python) random number generator

    while not inBounds(i, DNA):
      i = Random(t)

    NewMotifs = []
    x = 0
    for motif in Motifs:
      if x == i: # ignore DNA[i]
        x += 1
        continue
      NewMotifs.append(motif)
      x += 1

    CountMatrix = calculateCountMatrix(NewMotifs)

    for n in range(len(CountMatrix)):
      for m in range(len(CountMatrix[0])):
        CountMatrix[n][m] += 1

    Profile = calculateProfileMatrix(CountMatrix)

    Motifs[i] = GetRandomKMer(DNA[i], Profile, k)
    
    if Score(Motifs) < Score(BestMotifs):
      BestMotifs = Motifs
  return BestMotifs

def reverse(myStr):
    result = ""

    for c in myStr:
        result = c + result

    return result

def reverseComplement(myStr):
    result = ""
    
    myStr = reverse(myStr)

    for c in myStr:
        if c == 'A':
            result += 'T'
        elif c == 'T':
            result += 'A'
        elif c == 'C':
            result += 'G'
        elif c == 'G':
            result += 'C'

    return result

def parseSequence(seq, start, end):
    result = ""
    
    start -= 1
    end -= 1
    
    k = 500
    
    if start < end:
        firstPart = seq[0:start]
        if len(firstPart) <= k:
            result = seq[0:start]
        else:
            result = seq[start-k:start]
    elif start > end:
        myStr = reverse(seq[start+1:start+1+k])
        myStr = reverseComplement(myStr)
        result = myStr

    return result

parseSequence("ACTCGGCCGTACA", 9, 5) # reverse complement of CAT

def generateDNA(filename):
    result = []
    
    myMap = {} # index => chromosome, start_coord, end_coord
    
    e = open("data.txt", "r")

    currLine = e.readline()
    while currLine:
        data = currLine.split(' ')
        myKey = int(data[0])
        myData3 = data[3].replace("\n", "")
        myValue = [int(data[1]), int(data[2]), int(myData3)]
        myMap[myKey] = myValue
        
        currLine = e.readline()
    e.close()

    f = open(filename, "r")
  
    currLine = f.readline()
    while currLine:
        if int(currLine) in myMap.keys():
            chromosome = myMap[int(currLine)][0]
            start_coord = myMap[int(currLine)][1]
            end_coord = myMap[int(currLine)][2]

            myFile = str(chromosome) + ".txt"
            g = open('chromosomes/' + myFile, "r")

            chromosomeSequence = ""
            chromosomeCurrLine = g.readline()
            while chromosomeCurrLine:
                toAdd = chromosomeCurrLine
                toAdd = toAdd.replace("\n", "")
                chromosomeSequence += toAdd
                chromosomeCurrLine = g.readline()
                
            g.close()

            # extract the upstream or reverse complement of downstream and add to result
            parsedSequence = parseSequence(chromosomeSequence, start_coord, end_coord)
            result.append(parsedSequence)
            
        currLine = f.readline()
    f.close()

    return result

DNA1 = generateDNA("dna1.txt")
DNA2 = generateDNA("dna2.txt")
DNA3 = generateDNA("dna3.txt")

Gibbs1 = GibbsSampler(DNA1, 16, 1737, 100)
Gibbs2 = GibbsSampler(DNA2, 16, 462, 100)
Gibbs3 = GibbsSampler(DNA3, 16, 936, 100)

print(Gibbs1)
print(Gibbs2)
print(Gibbs3)

print(buildConsensusString(Gibbs1))
print(buildConsensusString(Gibbs2))
print(buildConsensusString(Gibbs3))
