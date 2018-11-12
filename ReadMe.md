# DeepMASS
***
DeepMASS is a known-to-unknown metabolite identification workflow, which includes a deep-learning based model to predict structural similarity between the unknown metabolites and the known ones based on their MS/MS spectra and a rank method for picking out the possible candidate structures of the unknowns.

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
2. Download the repo and unzip data/spectra.zip

## Usage
Take **example.py** as an example for metabolite identification.   
If you want to train a model based on your in-house database, please put your spectra files in to **data/spectra** directory and run **test.py**.