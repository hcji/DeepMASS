# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 12:42:39 2018

@author: jihon
"""

# transform structuredb to metfrag db type
import pandas as pd 
from rdkit import Chem

# transform structuredb to metfrag db type
metdb = pd.DataFrame(columns=['Identifier', 'MonoisotopicMass', 'MolecularFormula', 'InChIKey1', 'SMILES', 'Name', 'InChI'])
StructureDB = pd.read_table('data/Database/MsfinderStructureDB-VS12.esd')
for i in StructureDB.index:
  Identifier = 'MSF' + str(i)
  MonoisotopicMass = StructureDB['Exact mass'][i]
  MolecularFormula = StructureDB['Formula'][i]
  InChIKey1 = StructureDB['InChIkey'][i]
  SMILES = StructureDB['SMILES'][i]
  Name = StructureDB['ChEBI'][i] # use for judge whether the identification is right
  try:
      Mol = Chem.MolFromSmiles(SMILES)
      SMILES = Chem.MolToSmiles(Mol) # standardize smiles
      InChi = Chem.MolToInchi(Mol) # write inchi
  except:
      continue
  metdb.loc[i] = [Identifier, MonoisotopicMass, MolecularFormula, InChIKey1, SMILES, Name, InChi]
  
metdb.to_csv('support/metdb.csv', index=False)