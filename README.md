-------------------------------------------------------
Mitchell Pavlak
mpavlak1@jhu.edu

Brandon Liu
brandonmliu@gmail.com

Emmanuel Osikpa
eosikpa1@jhu.edu

Abbey Thrope
Abbey.thrope@gmail.com

# Prediction of ADRs using Extra Trees and 3D Molecule Data
##### *September 15th, 2019*

-------------------------------------------------------
__Part A:__
===========

The program is separated into four pieces. The primary components are Main.py,
getSimilarity.py, and saveMethods.py. plotting.py generates a matplotlib plot
to help show the results of the classifiers over all ADR tags.

Running steps vary based on whether or not a save file is available. If not,
run `python getSimilarity.py`. Keep all the folder structure the same as upon download.
Then run `python Main.py`. You can also run `python plotting.py`.

All code used in this project should be our own work or open-sourced code.
Some specific python libraries used include scikit-learn, subprocess,
multiprocessing, numpy, tqdm, pandas, cPickle.

Additionally, LSAlign was used to source PCScore and get rigid small molecule
alignment. LSAlign is open for any academic reuse or modification. For more
information about LSAlign, see: https://zhanglab.ccmb.med.umich.edu/papers/2018_1.pdf

See also the following:
https://academic.oup.com/bioinformatics/article/32/15/2338/1744048

__Part B:__
===========

The goal of this project was to leverage 3D molecular structure data as a
potential area for improvement of ADR prediction. We used LSAlign to align the
small molecules. Instead of taking more detailed measures of 3D shape, we
used a PCScore as developed in the top paper above (LSAlign). This metric allowed
us to quantify the similarity in 3D structure, bond activity etc between various
molecules. By training with a distance matrix of these similarities in addition
to gene expression data, we were able to see slight improvements in accuracy
over using the full 2d fingerprint of the molecule with GE data. Further, our
3D model required little over a third of the training time comparable 2D
models we created required.
