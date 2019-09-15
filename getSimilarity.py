import subprocess
import multiprocessing as mp
import os
import numpy as np


toExecutable = "./LSsrc/LSalign"
toExamples = "./LS_Example/"

# print("ls "+toExamples)
filenames = subprocess.check_output("ls "+toExamples,stderr = subprocess.STDOUT,shell=True).split("\n")[0:-1]
with open("out.txt",'w') as fw:
    fw.write("")


print(filenames)
try:
    subprocess.check_output(toExecutable  +" "+ toExamples +filenames[0] + " "+toExamples +filenames[1] + "> out.txt",stderr = subprocess.STDOUT,shell=True)
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

kOrder = moleculeSimDict.keys()

np.ndarray()
