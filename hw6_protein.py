"""
Protein Sequencing Project
Name:
Roll Number:
"""

import string

from numpy import char
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f = open(filename, "r") 
    lines = f.read().replace("\n","")
    return lines


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    list1=[]
    dna=dna[startIndex:].replace("T", "U")
    for i in range(0,len(dna),3):
        list1.append(dna[i:i+3])
    for i in list1:
        if i=="UAA"or i=="UGA" or i=="UAG":
          new=list1.index(i)
          return list1[:new+1]   
    return list1


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    dict1={}
    f = open(filename, "r")
    j = json.load(f)
    for k,v in j.items():
        for r in v:
            new=r.replace("T","U")
            dict1[new]=k
    return dict1


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    list1=[]
    for k in codons:
        for v in codonD:
            if k==v:
                list1.append(codonD[v])
                if list1[0]=="Met":
                 list1[0]="Start"
    return list1


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    new=readFile(dnaFilename)
    new1=makeCodonDictionary(codonFilename)
    list1=[]
    i=0
    unused_bias=0
    while i<len(new):
        list2=new[i:i+3]
        if list2=="ATG":
            r=dnaToRna(new,i)
            s=generateProtein(r,new1)
            list1.append(s)
            i+=3*len(s)
        else:
            i+=1
            unused_bias+=1
    print("total number of bases",len(new),"unused-base count",unused_bias,"total number of proteins synthesized",len(list1))
    return list1


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    list1=[]
    for i in proteinList1:
        if i in proteinList2:
            list1.append(i)
    return list1


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    list1=[]
    for i in proteinList:
        for j in i:
            list1.append(j)
    return list1


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    dict1={}
    for i in aaList:
        if i not in dict1:
            dict1[i]=aaList.count(i)   
    return dict1


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    new=combineProteins(proteinList1)
    new1=combineProteins(proteinList2)
    new2=aminoAcidDictionary(new)
    new3=aminoAcidDictionary(new1)
    list1=[]
    for i,j in new2.items():
        if i!="Start" and i!="Stop" and i not in list1:
            list1.append(i)
    for i,j in new3.items():
        if i!="Start" and i!="Stop" and i not in list1:
            list1.append(i)
    frequency=[]
    for i in list1:
        f1=0
        f2=0
        if i in new:
            f1=new2[i]/len(new) 
        if i in new1:
            f2=new3[i]/len(new1)
        diff=f2-f1
        if diff>cutoff or diff<-cutoff:
            frequency.append([i,f1,f2])
    
    return frequency


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    list1=[]
    for i in commonalities:
        i.remove("Start")
        i.remove("Stop")
        if len(i)>1:
            new='-'.join(i)
            list1.append([new])
        else:
            if i not in list1:
                list1.append(i)
    new=sorted(list1)
    for i in new:
        for j in i:
            print(j)
    for i in differences:
        print(i[0]+":"+str(round(i[1]*100,2))+"%"+ " in seq1 ,"+str(round(i[2]*100,2))+"%"+"in seq2")

    return 


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    #print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    #test.week1Tests()
    #print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    #runWeek1()
    #test.testReadFile()
    #test.testDnaToRna()
    #test.testMakeCodonDictionary()
    #test.testGenerateProtein()
    #test.testSynthesizeProteins()
    #test.testCommonProteins()
    #test.testCombineProteins()
    #test.testAminoAcidDictionary()
    #test.testFindAminoAcidDifferences()
    runWeek2()

    ## Uncomment these for Week 2 ##
    """
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    """

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
