# DeepMASS
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
	
## Installation
Download the repo and unzip data/spectra.zip

## Usage
Take run_one_example.py as an example for metabolite identification. If you want to train a model based on your in-house database, please put your spectra files in to data/spectra directory and run test.py

## Note
Due to the restrictions of file size, the repo in GitHub do not included the trained model. You can download them in [Gitee mirror](https://gitee.com/hcji/KMet)
