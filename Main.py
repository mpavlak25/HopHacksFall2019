import subprocess
import multiprocessing as mp
import numpy as np,math
from heapq import nlargest
from tqdm import tqdm
import csv
import saveMethods
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import ExtraTreesClassifier
import pandas as pd
import cPickle


#Taking in names mapping
with open("./InData/178 molecules with 3D structures.csv") as mapFile:
    r = csv.DictReader(mapFile)
    pertPubMap = {row["pert_id"]:row['pubchem_cid'] for row in r}

#checks for precalculated 3d similarity scores
if (saveMethods.loadSave()):
    moleculeSimDict = saveMethods.loadSave()
else:
    print("generate the mol2 data using getSimilarity.py... exiting...") #generate using getSimilarity.py
    exit(1)

#Prepared adverse drug reaction data
with open("./InData/FAERS 90-10.csv") as adrData:
    adrReader = csv.reader(adrData)
    categories = []
    pertToValDict = {}
    i = 0
    for lines in adrReader:
        if i == 0:
            labelNames = lines[1:len(lines)]
        else:
            pertToValDict[lines[0]] = map(float,lines[1:len(lines)])
        i += 1


pertOrder = list(pertToValDict.keys())
labels = [pertToValDict[k] for k in pertOrder]



def getKey(in1,in2):
    """Returns a consistent ordering of key of two elements, the one with the
    lesser hash in tuple at index 0.
    """
    i1 = pertPubMap[in1]
    i2 = pertPubMap[in2]
    if (hash(i1) < hash(i2)):
        return(i1,i2)
    else:
        return(i2,i1)


dataIn = []
for k1 in pertOrder:
    temp = []
    for k2 in pertOrder:
        temp.append(moleculeSimDict[getKey(k1,k2)])
    dataIn.append(temp)

#Opening prepared gene expression data.
with open("./InData/LINCS_Gene_Expression_signatures_178.csv") as GEData:
    GEReader = csv.reader(GEData)
    categories = []
    pertToGEValDict = {}
    i = 0
    for lines in GEReader:
        if i == 0:
             GELabelNames = lines[1:len(lines)]
        else:
                pertToGEValDict[lines[0]] = map(float,lines[1:len(lines)])
        i += 1


numD = len(dataIn)
dataInComparison = [[] for i in dataIn]
for i,k1 in enumerate(pertOrder):
    temp = []
    for k2 in pertOrder:

        dataIn[i].extend(pertToGEValDict[k2])
        dataInComparison[i].extend(pertToGEValDict[k2])
dataUnbiased = pd.DataFrame(dataInComparison) #this is before any 2d specific

print("\n\n")
print("3D data complete. Finishing 2D...")

#Getting the 2D fingerprints
with open("./InData/MACCS 178.csv") as MacData:
    MacReader = csv.reader(MacData)
    categoriesM = []
    pertToMacValDict = {}
    i = 0
    for lines in MacReader:
        if i == 0:
             MacLabelNames = lines[1:len(lines)]
        else:
                pertToMacValDict[lines[0]] = "".join(lines[1:len(lines)])
        i += 1

#adding 2d fingerprints for 2d data
for i,k1 in enumerate(pertOrder):
    dataInComparison[i].extend(pertToMacValDict[k1])


dataInput = pd.DataFrame(dataIn)
dataCompInput = pd.DataFrame(dataInComparison)
y = pd.DataFrame(labels)
ysplits = [y.iloc[:,i] for i in range(0,len(labelNames))]


ncores = mp.cpu_count()

# The following code is how we calculated the most strongly correlated adrs for
# use in classification. This will be saved for loading, it takes a while
# to run. If you choose to update save may need to update this.
#choosing not to comment it out solely to demonstrate indexing stays intact.
def corrFinder(i):
    split = ysplits[i]
    corrs = dataUnbiased.corrwith(split,axis=1).values.tolist()
    corrs = map(abs,corrs)
    return(i,sum(corrs))
c =corrFinder(18)

print("Corr data intact: " + str(c[1] > 39.12 and c[1] <39.13))
##This is the commented out code that the result is saved below.
# print("Generating the correlations. This will be resource intensive.")
# corrPool = mp.Pool(ncores)
# corrRes = list(tqdm(corrPool.imap(corrFinder,range(0,len(ysplits))),total = len(ysplits)))
# corrPool.close()
# corrPool.join()
#
# print(len(corrRes))
# features = nlargest(75,corrRes,key= lambda best: best[1])
# print(features)
# print("The above filters ")
# print([i[0] for i in features])
# [(18, 39.1215229740508), (951, 38.30917790758762), (777, 37.90343653707833),
# (600, 36.570675293419114), (515, 34.63535137905087), (516, 33.945624250690265),
# (386, 31.43589116630935), (517, 31.25738387630076), (673, 31.124840619481258),
# (529, 30.52846326122847), (160, 30.470445553889096), (633, 30.162501607702687),
# (341, 29.989534134400756), (331, 29.816319863871232), (433, 29.564464990711823),
# (1026, 29.388626314458794), (912, 29.382006534181212),(296, 28.921660418830207),
# (26, 28.555322534698277), (214, 28.151576678741915), (254, 27.3185413198628),
# (113, 27.07670456265395), (753, 26.8263687006867), (422, 26.457292058679137),
# (944, 26.29526845297295), (671, 26.28635919988787), (235, 26.079793456781868),
# (544, 25.93496808340084), (405, 25.926439574742517), (724, 25.920356724712402),
# (334, 25.801298706633734), (25, 25.604362712998377), (196, 25.302387011100922),
# (88, 25.190602870118248), (421, 24.91903353893394), (726, 24.506134848329793),
# (343, 24.502700716644185), (174, 24.37964570845944), (1073, 24.309651924735387),
# (396, 24.142122739156445), (617, 23.989008140998983), (244, 23.94408289678512),
# (22, 23.78729029255355), (969, 23.777963280049498), (966, 23.729877828804145),
# (835, 23.595816212077974), (225, 23.40807281101535), (996, 23.289185781203994),
# (656, 23.19752314938246), (217, 22.837560428090413), (63, 22.69336656935906),
# (953, 22.250149774028817), (920, 22.00002386171841), (800, 21.956601222252427),
# (891, 21.857212796563083), (1028, 21.84583331192912),(1033, 21.797758246581797),
# (178, 21.73792761023484), (659, 21.734258812996416), (35, 21.729679611456987),
# (650, 21.692004455195093), (342, 21.619715012752), (460, 21.570156780859584),
# (401, 21.48555771276663), (1067, 21.41910060105994), (186, 21.346528147144575),
# (950, 21.34538308032179), (324, 21.342499510101568), (452, 21.33888790680898),
# (629, 21.33553655755084), (84, 20.993097932718864), (576, 20.8906311353427),
# (270, 20.692883218789138), (555, 20.643156350450802),
# (917, 20.579089123929926)]

runIndices =[18, 951, 777, 600, 515, 516, 386, 517, 673, 529, 160, 633, 341,
331, 433, 1026, 912, 296, 26, 214, 254, 113, 753, 422, 944, 671, 235, 544, 405,
724, 334, 25, 196, 88, 421, 726, 343, 174, 1073, 396, 617, 244, 22, 969, 966,
835, 225, 996, 656, 217, 63, 953, 920, 800, 891, 1028, 1033, 178, 659, 35, 650,
342, 460, 401, 1067, 186, 950, 324, 452, 629, 84, 576, 270, 555, 917]

# If you'd like to see label names.
# for i in runIndices:
#     print(labelNames[i])
print("Retrieved categories.")


res = []

def train(i):
    """Takes in int i that represents index of label dataframe within ysplits.
    Meant for use with multiprocessing.
    """
    yCurr = ysplits[i]
    treesModel = ExtraTreesClassifier(n_estimators = 100)
    score = cross_val_score(treesModel, currdataInput, yCurr, cv=3, scoring='accuracy')
    return([i,treesModel,score])

print("Training 3D model...")
currdataInput = dataInput
procPool = mp.Pool(ncores)

res = list(tqdm(procPool.imap(train,runIndices),total = len(runIndices)))
procPool.close()
procPool.join()
print("Training 2D comparison...")
currdataInput = dataCompInput
procPool = mp.Pool(ncores)
resComp = list(tqdm(procPool.imap(train,runIndices),total = len(runIndices)))
procPool.close()
procPool.join()

tot = 0
totLast = 0.0
x3D = []
y3D = []

for i in res:
    print("Score 3D: " + str(i[0]) + " " + str(i[2]))
    x3D.append(i[0])
    x3D.append(i[0])
    x3D.append(i[0])
    y3D.append(i[2][0])
    y3D.append(i[2][1])
    y3D.append(i[2][2])
    temp = map(float,i[2])
    tot += temp[0]+temp[1] + temp[2]
    totLast += temp[2]
mean3D = tot/(3*len(res))
meanLastFold3D = totLast/len(res)
print("3D mean: " + str(mean3D))
print("3D mean of last fold: " + str(meanLastFold3D))
totcomp = 0
totLastComp = 0
x2D = []
y2D = []
for i in resComp:
    print("Score 2D: " + str(i[0]) + " " +str(i[2]))
    x2D.append(i[0])
    x2D.append(i[0])
    x2D.append(i[0])
    y2D.append(i[2][0])
    y2D.append(i[2][1])
    y2D.append(i[2][2])
    temp = map(float,i[2])
    totcomp += temp[0]+temp[1] + temp[2]
    totLastComp += temp[2]
mean2Dcomp = totcomp/(3*len(resComp))
mean2Dlast = totLastComp/len(resComp)
print("2D mean: " + str(mean2Dcomp))
print("2D mean of last fold: " + str(mean2Dlast))



cPickle.dump([x2D,y2D,x3D,y3D],open("./SavesTrained/toPlotSaves.pkl","wb"))
cPickle.dump([res,resComp],open("./SavesTrained/resandcompSavedTrees.pkl","wb"))
