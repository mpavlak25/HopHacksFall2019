import subprocess
import multiprocessing as mp
import os
import numpy as np
import cPickle


toExecutable = "./LSsrc/LSalign"
toExamples = "./LS_3DData/"

def checkSave():
    """Checks if a save file exists.
    """
    if os.path.isfile("./Saves/savefile.pkl") and os.path.getsize("./Saves/savefile.pkl") > 0:
        return raw_input("Use found save? (y/n): ") == "y"
    return False

def loadSave():
    """Loads a save file and returns the object.
    """
    with open("./Saves/savefile.pkl") as f:
        order = cPickle.load(f)
    return order
