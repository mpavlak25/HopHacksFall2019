import subprocess
import multiprocessing as mp
import os
import numpy as np
import cPickle
import saveMethods

#This uses the freely available (for academic work) LSAlign to perform rigid
#alignment (didn't have sufficient processing power on my laptop to perform
#flexible alignment), and then generating similarity scores. As per definitions
#of PCScore in https://zhanglab.ccmb.med.umich.edu/papers/2018_1.pdf and
#http://zhanglab.ccmb.med.umich.edu/LS-align/ReportParam.html the max is used
#for similarity.

toExecutable = "./LSsrc/LSalign"
toExamples = "./LS_3DData/"

if saveMethods.checkSave():
    moleculeSimDict = saveMethods.loadSave()
else:
    filenames = subprocess.check_output("ls "+toExamples,stderr = subprocess.STDOUT,shell=True).split("\n")[0:-1]
    with open("out.txt",'w') as fw:
        fw.write("")
    try:
        print("Aligning...")
        subprocess.check_output(toExecutable  +" "+ toExamples +filenames[0] + " "+toExamples +filenames[0] + "> out.txt",stderr = subprocess.STDOUT,shell=True)
    except:
        if os.path.getsize("./out.txt") == 0:
            print("An error occurred. Exiting...")
            exit()

    similarity = []
    with open("out.txt",'r') as fr:
        for line in fr:
            similarity.append(line.split())
    similarity = similarity[3:-5]
    keys = set(s[0] for s in similarity)

    print(similarity)
    moleculeSimDict = {}
    for row in similarity:
        if hash(row[0]) < hash(row[1]):
            moleculeSimDict[(row[0],row[1])] = (row[2],row[3])
        else:
            moleculeSimDict[(row[1],row[0])] = (row[2],row[3])

    moleculeSimDict = {k:max(map(float,v)) for k,v in moleculeSimDict.items()}
    cPickle.dump(moleculeSimDict,open("./Saves/savefile.pkl","wb"))
