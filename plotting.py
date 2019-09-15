import matplotlib.pyplot as plt
import pandas as pd
import cPickle
import numpy as np,math

#Fairly simple program for generating a quick plot

saved = cPickle.load(open("./SavesTrained/toPlotSaves.pkl","rb"))
x2D,y2D,x3D,y3D = saved[0],saved[1],saved[2],saved[3]

plt.scatter(x2D,y2D)
plt.scatter(x3D,y3D)
plt.plot([-500,1500],[1,1],color="r")
plt.xlim([-500,1500])
plt.ylim([0,1.2])

plt.title("ADR prediction accuracy; 2D vs 3D data")
plt.ylabel("Percent correct by tag from all 3 folds")
plt.xlabel("ADR tag index")

plt.show()
