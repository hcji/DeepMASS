# DeepMASS
***
DeepMASS is a known-to-unknown metabolite identification workflow, which includes a deep-learning based model to predict structural similarity between the unknown metabolites and the known ones based on their MS/MS spectra and a rank method for picking out the possible candidate structures of the unknowns.

<div align="center">
<img src="https://github.com/hcji/DeepMASS/blob/master/support/figure.png" width=600 height=450 />
</div>

## Motivation
Metabolite identification is a great challenge restricting the metabolomics study. One of the main reasons is the limited number of high-quality MS/MS spectra of metabolites included in database. Using the transformational relationship and structural similarity between metabolites is a promising strategy to extend the number of metabolites can be identified by the existing database. Hence, we present the DeepMASS workflow, which is proved to be an effective way to predict the structural similarities between the unknown metabolite and the metabolites existing in the database.

## Depends
	Anaconda for python 3.6
	# RDKit
	conda install libboost=1.65.1
	conda install boost=1.65.1
	conda install boost-cpp=1.65.1
	conda install -c rdkit rdkit
	# Keras
	conda install keras
	# IsoSpec
	pip install IsoSpecPy == 1.0.7
	
## Installation
1. Download the model [here](https://www.researchgate.net/profile/Hongchao_Ji/publication/328822822_DeepMASS_Model_for_Deep_MSMS-Aided_Structural-similarity_Scoring_for_Unknown_Metabolites_Identification/data/5be4d6a5299bf1124fc41e39/model-40V.zip), then unzip it into **/model/40V** directory.
2. Clone the repo and put spectra files in **data/spectra** directory

## Note
### 2019/04/29
Under the requirement of the owner of experimental spectra (MetDNA dataset), the dataset has been removed. If you have already download the dataset, please keep it private, and the dataset can be only used to reproduce the results of DeepMASS paper. If you need the dataset for other use, please contact with the [owner](http://www.metabolomics-shanghai.org/software.php).   
### 2019/05/13
Since the experimental spectra has been removed, this package cannot be run directly. You can train your own model with your in-house database. Otherwise, you can wait some days. I m trying to include spectra from public databases.  

## Usage
Take **example.py** as an example for metabolite identification.   
If you want to train a model based on your in-house database, please put your spectra files in to **data/spectra** directory and run **test.py**.
