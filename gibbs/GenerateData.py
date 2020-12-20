# Reference/Citation: https://www.w3schools.com/python/python_file_write.asp

def generateData():
    e = open("data.txt", "a")
    f = open("orf.txt", "r")
  
    currLine = f.readline()
    currLine = f.readline()
    myStrCount = 0
    while currLine:
        data = currLine.split('\t')

        if data[0] == "\n":
            break
        if data[1] != '': # ORF exists
            if data[4] == '' or data[6] == '' or data[7] == '':
                currLine = f.readline()
            else:
                myStr = data[0] + ' ' + data[4] + ' ' + data[6] + ' ' + data[7]
                myStrCount += 1

                if myStrCount >= 1:
                    e.write('\n')
                    e.write(myStr)
                else:
                    e.write(myStr)
            
                currLine = f.readline()
        else: # ORF doesn't exist
            currLine = f.readline()

    e.close()
    f.close()

    print("Data extracted")

generateData()
