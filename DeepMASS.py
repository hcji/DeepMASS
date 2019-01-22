# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 14:49:47 2018

@author: jihon
"""

import re
import gc
import os
import sys
import bisect
import numpy as np
# import numpy.fft as nfft
import pandas as pd
from scipy import signal
from scipy.sparse import csr_matrix, vstack, hstack
from IsoSpecPy import IsoSpecPy
from itertools import chain
import json
import requests
import pubchempy as pc
from bs4 import BeautifulSoup
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs
from keras.models import Model, model_from_json
from keras.layers import Dense, Input, Dropout, concatenate
from keras.optimizers import SGD
from sklearn.preprocessing import MinMaxScaler
from sklearn.externals import joblib
import matplotlib.pyplot as plt

from MassPlot import visualize_molecule

formulaDB = pd.read_table('data/Database/MsfinderFormulaDB-VS10.efd')
structureDB = pd.read_table('data/Database/MsfinderStructureDB-VS12.esd')

def progressBar(count,total,size):
    percent = float(count)/float(total)*100
    sys.stdout.write("\r" + str(int(count)).rjust(3,'0')+"/"+str(int(total)).rjust(3,'0') + ' [' + '='*int(percent/10)*size + ' '*(10-int(percent/10))*size + ']')

def get_score(pred, real):
    '''
    Task: 
        Evaluate the matching degree of DeepMASS score and FP score with dot product
        when calculating FP score, some molecule may raising error, *compare_structure*
        function will return -1 in such case. Here they are removed.
    Parameters:
        pred: array, DeepMASS score
        real: array, Fingerprint score
    '''
    keep = np.where(real >= 0)[0] 
    if len(keep)<1:
        return 0
    else:
        score = sum(np.array(pred[keep]) * np.array(real[keep])) / len(keep)
    return round(score, 2)


def read_ms(csv, precursor = None, norm = True):
    '''
    Task: 
        Read spectrum from csv file
    Parameters:
        csv: str, file path of spectrum
        precursor: folat, m/z of precursor
        norm: logic, whether to normalizate the spectrum or not
    '''
    spec = pd.read_csv(csv)
    spec = spec.iloc[:,range(2)]
    spec.columns = ['mz', 'intensity']
    if norm:
        spec['intensity'] = spec['intensity'] / max(spec['intensity'])
    if precursor is not None:
        keep = spec['mz'] < precursor + 1
        spec = spec.loc[keep]
    return spec
	

def isotope_pattern(formula, thres=0.99):
    '''
    Task: 
        Generate theoretical isotope distribution
    Parameters:
        formula: str, chemical formula
        thres: folat, abundance threshold
    '''
    if type(formula) is not str:
        raise ValueError('input formula must be a character')
    isotope = IsoSpecPy.IsoSpec.IsoFromFormula(formula, thres)
    isotope = isotope.getConfsNumpy()
    output = pd.DataFrame({'mass': isotope[0], 'intensity': np.exp(isotope[1])})
    return output
    

def compare_isotope(measured, expected, tolerance=0.001):
    '''
    Task: 
        Compare theoretical isotope distribution and measured isotope distribution
    Parameters:
        measured: DataFrame, measured isotope distribution
        expected: DataFrame, theoretical isotope distribution
        tolerance: float, m/z tolerance
    '''
    if (type(measured) is not pd.DataFrame) or (type(expected) is not pd.DataFrame):
        raise ValueError('input data must be pandas.DataFrame')
    measured['intensity'] = measured['intensity']/sum(measured['intensity'])
    expected['intensity'] = expected['intensity']/sum(expected['intensity'])
    expected_m0 = expected['intensity'][0]
    measured_m0 = measured['intensity'][0]
    expected_m1 = sum(expected['intensity'][(expected['mass'] > expected['mass'][0] - tolerance + 0.997) & (expected['mass'] < expected['mass'][0] + tolerance + 1.006)])
    measured_m1 = sum(measured['intensity'][(measured['mass'] > measured['mass'][0] - tolerance + 0.997) & (measured['mass'] < measured['mass'][0] + tolerance + 1.006)])
    expected_m2 = sum(expected['intensity'][(expected['mass'] > expected['mass'][0] - tolerance + 1.994) & (expected['mass'] < expected['mass'][0] + tolerance + 2.013)])
    measured_m2 = sum(measured['intensity'][(measured['mass'] > measured['mass'][0] - tolerance + 1.994) & (measured['mass'] < measured['mass'][0] + tolerance + 2.013)])
    expected_m3 = sum(expected['intensity'][(expected['mass'] > expected['mass'][0] - tolerance + 2.991) & (expected['mass'] < expected['mass'][0] + tolerance + 3.018)])
    measured_m3 = sum(measured['intensity'][(measured['mass'] > measured['mass'][0] - tolerance + 2.991) & (measured['mass'] < measured['mass'][0] + tolerance + 3.018)])
    score = (1 - abs(expected_m0 - measured_m0)) * (1 - abs(expected_m1 - measured_m1)) * (1 - abs(expected_m2 - measured_m2)) * (1 - abs(expected_m3 - measured_m3))
    return score


def ms2vec(ms, precision=0.01, maxmz=2000):
    '''
    Task: 
        Convert mass spectrum to (sparse) vector
    Parameters:
        ms: DataFrame, ms spectrum
        precision: float, related to binning size 
        maxmz: int, maximum m/z, related to length of vector
    '''
    # default, convert to sparse vector
    bits = int(1/precision)
    idx = round(ms['mz']*bits)
    val = ms['intensity']
    keep = idx < maxmz*bits
    idx = idx[keep]
    val = val[keep]
    vec = csr_matrix((val, (np.zeros(len(idx)), idx)), shape=(1, bits*maxmz))
    
    # old codes, convert to dense vector
    '''
    vec = np.zeros(10*maxmz)
    for i in range(len(idx)):
        if int(idx[i]) < len(vec):
            vec[int(idx[i])] = val[i]
    '''
    return vec


def formula2vec(formula, elements=['C', 'H', 'O', 'N', 'P', 'S']):
    '''
    Task: 
        Convert formula vector
    Parameters:
        formula: str, chemical formula
        elements: str list, elements
    '''
    formula_p = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    vec = np.zeros(len(elements))
    for i in range(len(formula_p)):
        ele = formula_p[i][0]
        num = formula_p[i][1]
        if num == '':
            num = 1
        else:
            num = int(num)
        if ele in elements:
            vec[elements.index(ele)] += num
    return np.array(vec) 


def fft_ccor(ts1, ts2, norm=True, output_type='sparse'):
    '''
    Task: 
        Calculate fft cross correlation for (sparse) vector
    Parameters:
        ts1: csr_matrix, signal 1
        ts2: csr_matrix, signal 2
        norm: logic, norm the output
        output_type: str, output sparse vector or dense vector
    '''
    
    # convert to dense vector
    if type(ts1) is csr_matrix:
        ts1 = ts1.A.squeeze()
    if type(ts2) is csr_matrix:
        ts2 = ts2.A.squeeze()
    
    # old code, use fft
    '''
    fft_ts1 = nfft.fft(ts1)
    fft_ts2 = nfft.fft(ts2)
    output = ((1 / (1. * len(ts1))) * nfft.ifft(fft_ts1 * np.conjugate(fft_ts2)).real)
    '''
    
    # use scipy.sinal.correlate
    output = signal.correlate(ts1, ts2, mode='same', method='fft')
    maxccor = max(output)
    if norm:
        output = output / maxccor
    if output_type == 'sparse':
        output[output < maxccor/10**6] = 0
        output = csr_matrix(output)
    
    return output

    
def sparse_pearsonr(A, B):
    '''
    Task: 
        Calculate pearson correlation for (sparse) vector
        **This function is not used any more in the current project**
    Parameters:
        A: csr_matrix, signal 1
        B: csr_matrix, signal 2
    '''
    # calculate base on definition.
    # since A & B are sparse, take mean(A) and mean(B) as 0
    a = A.multiply(B).sum()
    b = np.sqrt(A.multiply(A).sum()) * np.sqrt(B.multiply(B).sum()) + 10 ** -99
    res = a/b
    
    # convert to dense vector and use scipy.stat.pearsonr
    '''
    from scipy.stat import pearsonr
    AA = A.toarray()[0]
    BB = B.toarray()[0]
    res1, pval = pearsonr(AA, BB)
    abs(res-res1)
    '''
    return res


def sparse_cross_correlate(A, B):
    '''
    Task: 
        Calculate cross pearson correlation for (sparse) vector
        **This function is not used any more in the current project**
    Parameters:
        A: csr_matrix, signal 1
        B: csr_matrix, signal 2
    '''
    colA = A.tocoo().col
    colB = B.tocoo().col
    shifts = [colA-xx for xx in colB]
    shifts = np.unique(np.array(list(chain(*shifts))))
    shifts = shifts[abs(shifts) < A.tocoo().shape[1] / 2]
    ccors = []
    for s in shifts:
        if s < 0:
            ss = csr_matrix([], shape=(1, abs(s)))
            sB = hstack((B, ss))
            sB = sB.tocsr()[:, range(abs(s),sB.shape[1])]
        elif s == 0:
            sB = B
        elif s > 0:
            ss = csr_matrix([], shape=(1, abs(s)))
            sB = hstack((ss, B))
            sB = sB.tocsr()[:, range(A.tocoo().shape[1])]
        ccors.append(sparse_pearsonr(A, sB))
    idx = shifts + int(A.tocoo().shape[1]/2)
    res = csr_matrix((ccors, (np.zeros(len(idx)), idx)), shape=A.tocoo().shape)
    return res
  
    
def compare_structure(smiles1, smiles2, fp_type='Morgan', sim_type='Dice'):
    '''
    Task: 
        Compare structual similarity of two compound based on fingerprints.
    Parameters:
        smiles1: str, smiles of the compound 1
        smiles2: str, smiles of the compound 2
        fp_type: str, type of fingerprints
        sim_type: str, method for calculating similarity
    '''
    if fp_type == 'Morgan':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2, useFeatures=False)
    elif fp_type == 'MorganWithFeature':
        getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2, useFeatures=True)
    elif fp_type == 'MACCS':
        getfp = lambda smi: Chem.MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smi))
    elif fp_type == 'Topological':
        getfp = lambda smi: FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smi))
    elif fp_type == 'AtomPairs':
        getfp = lambda smi: Pairs.GetAtomPairFingerprint(Chem.MolFromSmiles(smi))
    
    try:
        fp1 = getfp(smiles1)
        fp2 = getfp(smiles2)
        if sim_type == 'Dice':
            sim_fp = DataStructs.DiceSimilarity(fp1, fp2)
        elif sim_type == 'Tanimoto':
            sim_fp = DataStructs.TanimotoSimilarity(fp1, fp2)
        elif sim_type == 'Cosine':
            sim_fp = DataStructs.CosineSimilarity(fp1, fp2)
        elif sim_type == 'Sokal':
            sim_fp = DataStructs.SokalSimilarity(fp1, fp2)
        elif sim_type == 'Russel':
            sim_fp = DataStructs.RusselSimilarity(fp1, fp2)
            
    except Exception as e:
        sim_fp = -1
    return sim_fp


def search_formula(mass, ppm):
    '''
    Task: 
        Search formula from formula database.
    Parameters:
        mass: float, exact mass of compound
        ppm: float, ppm
    '''
    mmin = mass - mass*ppm/10**6
    mmax = mass + mass*ppm/10**6
    lf = bisect.bisect_left(formulaDB['Exact mass'], mmin)
    rg = bisect.bisect_right(formulaDB['Exact mass'], mmax)
    formulas = list(formulaDB['Formula'][lf:rg])
    return formulas


def search_structure(formula):
    '''
    Task: 
        Search chemical structure from structural database.
    Parameters:
        formula: str, chemical formula
    '''
    structures = structureDB[structureDB['Formula']==formula]
    return structures


def search_pubchem(formula, timeout=999):
    '''
    Task: 
        Search chemical structure from pubchem.
    Parameters:
        formula: str, chemical formula
    '''
    # get pubchem cid based on formula
    cids = pc.get_cids(formula, 'formula', list_return='flat')
    idstring = ''
    smiles = []
    inchikey = []
    all_cids = []
    # search pubchem via formula with pug
    for i, cid in enumerate(cids):
        idstring += ',' + str(cid)
        if ((i%100==99) or (i==len(cids)-1)):
            url_i = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + idstring[1:(len(idstring))] + "/property/InChIKey,CanonicalSMILES/JSON"
            res_i = requests.get(url_i, timeout=timeout)
            soup_i = BeautifulSoup(res_i.content, "html.parser")
            str_i = str(soup_i)
            properties_i = json.loads(str_i)['PropertyTable']['Properties']
            idstring = ''
            for properties_ij in properties_i:
                smiles_ij = properties_ij['CanonicalSMILES']
                if smiles_ij not in smiles:
                    smiles.append(smiles_ij)
                    inchikey.append(properties_ij['InChIKey'])
                    all_cids.append(str(properties_ij['CID']))
                else:
                    wh = np.where(np.array(smiles)==smiles_ij)[0][0]
                    all_cids[wh] = all_cids[wh] + ', ' + str(properties_ij['CID'])
    
    result = pd.DataFrame({'InChIKey': inchikey, 'SMILES': smiles, 'PubChem': all_cids})
    return result
     

def get_data(pairs, spec_dir, includeY=True):
    '''
    Task: 
        get predictor (X) and label (y) for train/predict. 
    Parameters:
        pairs: DataFrame, DataFrame from csv file, take data/kegg_pairs_refine.csv as example.
        spec_dir: str, path of spectra files.
        includeY: logic, whether to generate the label. (True for train dataset, False for predict data)
    '''
    # data_X1 = []
    data_X1 = None
    data_X2 = []
    data_Y = []
    for i in pairs.index:
        progressBar(int(i+1), len(pairs.index), 1)
        try:
            spec1 = read_ms(spec_dir + '/' + pairs['kegg1'][i] + '.csv', pairs['mass1'][i]+1.0078)
            spec2 = read_ms(spec_dir + '/' + pairs['kegg2'][i] + '.csv', pairs['mass2'][i]+1.0078)
        except Exception as e:
            continue
        spec_vec1 = ms2vec(spec1)
        spec_vec2 = ms2vec(spec2)
        # spec_ccor = sparse_cross_correlate(spec_vec1, spec_vec2)
        spec_ccor = fft_ccor(spec_vec1, spec_vec2)
        
        formula1 = pairs['formula1'][i]
        formula2 = pairs['formula2'][i]
        mass1 = pairs['mass1'][i]
        mass2 = pairs['mass2'][i]
        mass_diff = np.array([mass1, mass2, abs(mass1-mass2)])
        formula_vec1 = formula2vec(formula1)
        formula_vec2 = formula2vec(formula2)
        formula_diff = list(abs(formula_vec1 - formula_vec2))
        
        if includeY:
            smi1 = pairs['smiles1'][i]
            smi2 = pairs['smiles2'][i]
            sim = compare_structure(smi1, smi2)
            if sim < 0:
                continue
            data_Y.append(sim)
        
        # data_X1.append(np.concatenate([spec_vec1, spec_vec2, spec_ccor]))
        data_X1 = vstack((data_X1, hstack((spec_vec1, spec_vec2, spec_ccor))))
        data_X2.append(np.concatenate([formula_vec1, formula_vec2, formula_diff, mass_diff]))
        
    # return np.array(data_X1), np.array(data_X2), np.array(data_Y)
    return csr_matrix(data_X1), np.array(data_X2), np.array(data_Y)


def delete_rows_csr(mat, indices):
    """
    The code is from stackoverflow:
        https://stackoverflow.com/questions/13077527/is-there-a-numpy-delete-equivalent-for-sparse-matrices
    Remove the rows denoted by ``indices`` form the CSR sparse matrix ``mat``.
    """
    if not isinstance(mat, csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    indices = list(indices)
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[indices] = False
    return mat[mask]


def new_data(mass, formula, spectrum_file, energy='40V'):
    '''
    Task: 
        get predictor (X) for new unknown spectra. 
    Parameters:
        mass: float, exact mass.
        formula: chemical formula.
        spec_dir: str, path of spectrum file.
        energy: only 40V is support until now.
    '''
    spec1 = read_ms(spectrum_file, mass+1.0078)
    spec_vec1 = ms2vec(spec1)
    mass1 = mass
    formula_vec1 = formula2vec(formula)
    
    known = pd.read_csv('data/known_compounds.csv')
    spec_dir = 'data/spectra/measured_spectra/' + energy
    # data_X1 = []
    data_X1 = None
    data_X2 = []
    for i in known.index:
        spec2 = read_ms(spec_dir + '/' + known['kegg'][i] + '.csv')
        spec_vec2 = ms2vec(spec2)
        # spec_ccor = sparse_cross_correlate(spec_vec1, spec_vec2)
        spec_ccor = fft_ccor(spec_vec1, spec_vec2)
        
        formula2 = known['formula'][i]
        mass2 = known['mass'][i]
        mass_diff = np.array([mass1, mass2, abs(mass1-mass2)])
        formula_vec2 = formula2vec(formula2)
        formula_diff = list(abs(formula_vec1 - formula_vec2))
        
        # data_X1.append(np.concatenate([spec_vec1, spec_vec2, spec_ccor]))
        data_X1 = vstack((data_X1, hstack((spec_vec1, spec_vec2, spec_ccor))))
        data_X2.append(np.concatenate([formula_vec1, formula_vec2, formula_diff, mass_diff]))
    
    return csr_matrix(data_X1), np.array(data_X2), known
    

def build_model(energy='40V', Test=False, Save=True):
    '''
    Task: 
        train DeepMASS model with simulated mass spectra. 
    Parameters:
        energy: only 40V is support until now.
        Test: logic, whether to split dataset for test.
        Save: logic, whether to save the model or not.
    '''
    kegg_pairs = pd.read_csv('data/kegg_pairs_refine.csv')
    random_pairs = pd.read_csv('data/random_pairs_simulated.csv')
    spec_dir = 'data/spectra/simulated_spectra/' + energy
    
    print('generate dataset')
    pairs = pd.DataFrame.append(kegg_pairs, random_pairs, ignore_index=True)
    data_X1, data_X2, data_Y = get_data(pairs, spec_dir, includeY=True)
    gc.collect()
    
    if Test:
        print('split dataset')
        test = np.random.choice(range(len(data_Y)), int(0.2*len(data_Y)))
        # train_X1 = np.delete(data_X1, test, axis=0)
        train_X1 = delete_rows_csr(data_X1, test)
        train_X2 = np.delete(data_X2, test, axis=0)
        train_Y = np.delete(data_Y, test, axis=0)
        test_X1 = data_X1[test,:]
        test_X2 = data_X2[test,:]
        test_Y = data_Y[test]
    else:
        train_X1 = data_X1
        train_X2 = data_X2
        train_Y = data_Y
    gc.collect()
    
    # scaler
    scaler = MinMaxScaler()
    scaler.fit(train_X2)
    train_X2 = scaler.transform(train_X2)
    if Test:
        test_X2 = scaler.transform(test_X2)
    gc.collect()

    # build model
    spec_input = Input(shape=(train_X1.shape[1],), sparse=True)
    # spec_input = Input(shape=(train_X1.shape[1],))
    formula_input = Input(shape=(train_X2.shape[1],))
    spec_output = Dense(128, activation='relu')(spec_input)
    formula_output = Dense(16, activation='relu')(formula_input)
    hidden = concatenate([formula_output, spec_output])
    hidden = Dense(32, activation='relu')(hidden)
    hidden = Dropout(0.3)(hidden)
    hidden = Dense(32, activation='relu')(hidden)
    hidden = Dropout(0.2)(hidden)
    predictions = Dense(1, activation='linear')(hidden)

    model = Model(inputs=[spec_input, formula_input], outputs=predictions)
    model.compile(optimizer=SGD(lr=0.05),loss='mse',metrics=['mae'])
    model.fit([train_X1, train_X2], train_Y, epochs=20)
    
    if Test:
        print(model.evaluate([test_X1, test_X2], test_Y, verbose=0))
        preds = model.predict([test_X1, test_X2])
        plt.figure(figsize=(8, 8))
        plt.plot(test_Y, preds, 'bo')
        plt.plot([0,1], [0,1], color ='red')
        plt.xlabel('similarity')
        plt.ylabel('prediction')
        plt.show()  
    
    if Save:
        joblib.dump(scaler, 'model/' + energy + '/scaler_X2.save') 
        model_json = model.to_json()  
        with open("model/" + energy + "/model_dnn_raw.json", "w") as json_file:  
            json_file.write(model_json)
        model.save_weights('model/' + energy + '/model_dnn_raw.h5')
    
    print('done')


def fine_tune(energy='40V', Test=False, Save=True):
    '''
    Task: 
        fine tune DeepMASS model with measured mass spectra. 
    Parameters:
        energy: only 40V is support until now.
        Test: logic, whether to split dataset for test.
        Save: logic, whether to save the model or not.
    '''
    kegg_pairs = pd.read_csv('data/kegg_pairs_refine.csv')
    random_pairs = pd.read_csv('data/random_pairs_measured.csv')
    spec_dir = 'data/spectra/measured_spectra/' + energy
    
    print('generate dataset')
    pairs = pd.DataFrame.append(kegg_pairs, random_pairs, ignore_index=True)
    data_X1, data_X2, data_Y = get_data(pairs, spec_dir, includeY=True)
    gc.collect()
    
    print('load model')
    json_file = open('model/' + energy + '/model_dnn_raw.json', 'r') 
    loaded_model_json = json_file.read() 
    json_file.close()  
    model = model_from_json(loaded_model_json)
    model.load_weights("model/" + energy + "/model_dnn_raw.h5")
    model.compile(optimizer=SGD(lr=0.02),loss='mse',metrics=['mae'])
    
    scaler = joblib.load('model/' + energy + '/scaler_X2.save')
    data_X2 = scaler.transform(data_X2)
    gc.collect()
    
    if Test:
        test = np.random.choice(range(len(data_Y)), int(0.2*len(data_Y)))
        # train_X1 = np.delete(data_X1, test, axis=0)
        train_X1 = delete_rows_csr(data_X1, test)
        train_X2 = np.delete(data_X2, test, axis=0)
        train_Y = np.delete(data_Y, test, axis=0)
        test_X1 = data_X1[test,:]
        test_X2 = data_X2[test,:]
        test_Y = data_Y[test]
        gc.collect()

        model.fit([train_X1, train_X2], train_Y, epochs=20)    
        print(model.evaluate([test_X1, test_X2], test_Y, verbose=0))
        preds = model.predict([test_X1, test_X2])
        plt.figure(figsize=(8, 8))
        plt.plot(test_Y, preds, 'bo')
        plt.plot([0,1], [0,1], color ='red')
        plt.xlabel('similarity')
        plt.ylabel('prediction')
        plt.show()
    else:
        model.fit([data_X1, data_X2], data_Y, epochs=20)        
  
    # save model
    if Save:
        model_json = model.to_json()  
        with open('model/' + energy + '/model_dnn_tuned.json', "w") as json_file:  
            json_file.write(model_json)
        model.save_weights('model/' + energy + '/model_dnn_tuned.h5')
    
    print('done')


def leave_one_out_test(energy='40V', database='structureDB'):
    '''
    Task: 
        leave one out test
    Parameters:
        energy: only 40V is support until now.
        database: search in which database? structureDB or pubchem.
    '''
    spec_dir = 'data/spectra/measured_spectra/' + energy
    knowns = pd.read_csv('data/known_compounds.csv')
    all_spectra = os.listdir(spec_dir)
    
    # load model
    json_file = open('model/' + energy + '/model_dnn_tuned.json', 'r') 
    loaded_model_json = json_file.read() 
    json_file.close()  
    model = model_from_json(loaded_model_json)
    model.load_weights("model/" + energy + "/model_dnn_tuned.h5")
    scaler = joblib.load('model/' + energy + '/scaler_X2.save')
    
    ranks = []
    totals = []
    for i in knowns.index:
        if knowns['kegg'][i] + '.csv' not in all_spectra:
            continue
        if database == 'pubchem':
            if np.isnan(knowns['pubchem'][i]):
                ranks.append(np.nan)
                totals.append(np.nan)
                continue
            else:
                true = str(int(knowns['pubchem'][i]))
            try:
                candidate = search_pubchem(knowns['formula'][i])
                candidate_cids = list(candidate['PubChem'])
                candidate_cids = [xx.split(', ') for xx in candidate_cids]
                candidate_smiles = list(candidate['SMILES'])
            except Exception as e:
                ranks.append(-1)
                totals.append(-1)
                continue
        elif database == 'structureDB':
            if type(knowns['chebi'][i]) is not str:
                ranks.append(np.nan)
                totals.append(np.nan)
                continue
            else:
                true = str(knowns['chebi'][i])
            candidate = search_structure(knowns['formula'][i])
            candidate_cids = list(candidate['ChEBI'])
            candidate_cids = [str(xx).split(';') for xx in candidate_cids]
            candidate_smiles = list(candidate['SMILES'])
            if len(candidate_smiles) == 0:
                ranks.append(-1)
                totals.append(-1)
                continue                
        
        a = knowns.iloc[i, range(4)]
        a = pd.concat([a] * (knowns.shape[0]-1), axis=1) 
        a = a.transpose()
        b = knowns.iloc[knowns.index != i, range(4)]
        pred_pairs = pd.concat([a.reset_index(drop=True), b.reset_index(drop=True)], axis=1)
        pred_pairs.columns = ['kegg1', 'smiles1', 'mass1', 'formula1', 'kegg2', 'smiles2', 'mass2', 'formula2']
        data_X1, data_X2, data_Y = get_data(pred_pairs, spec_dir, includeY=False)
        data_X2 = scaler.transform(data_X2)
        
        preds = model.predict([data_X1, data_X2])[:,0]
        neighbors = np.where(preds > 0.6)[0]
        if len(neighbors)<1:
            ranks.append(np.nan)
            totals.append(np.nan)
            continue
        
        smi_neighbors = pred_pairs['smiles2'][neighbors]
        scores = []
        for smi in candidate_smiles:
            this_sims = np.array([compare_structure(smi, x) for x in smi_neighbors])
            scores.append(get_score(preds[neighbors], this_sims))
            
        try:
            wh_true = np.where([true in xx for xx in candidate_cids])[0][0]
        except Exception as e:
            ranks.append(np.nan)
            totals.append(np.nan)
            continue
        
        '''
        true_sims = np.array([compare_structure(candidate_smiles[wh_ture], x) for x in smi_neighbors])
        plt.plot(preds[neighbors], true_sims, 'bo')
        plt.plot([0,1], [0,1], color ='red')
        '''
        
        rank = len(np.where(scores > scores[wh_true])[0])+1
        ranks.append(rank)
        totals.append(len(scores))
        print(str(knowns['kegg'][i]) + ': ' + str(rank) + '/' + str(len(scores)))
    
    result = pd.concat([knowns.reset_index(drop=True), pd.Series(ranks), pd.Series(totals)], axis=1)
    result.columns = ['kegg','smiles','mass','formula','pubchem', 'chebi','rank', 'totel']
    result.to_csv('result/search_result_' + database + '_'  + energy + '.csv', index=False)
    return result


def run_one_example(mass, formula, ms_file, energy='40V', thres = 0.5, database='structureDB'):
    '''
    Task: 
        run DeepMASS workflow with a new spectrum.
        **this is the only function for end user**
    Parameters:
        mass: float, exact mass.
        formula: str, chemical formula.
        ms_file: str, path of mass spectrum.
        energy: only 40V is support until now.
        thres: folat, 0-1, threshold of structural similarity of reference compound.
        database: search in which database? structureDB or pubchem.
    '''
    # load model
    json_file = open('model/' + energy + '/model_dnn_tuned.json', 'r') 
    loaded_model_json = json_file.read() 
    json_file.close()  
    model = model_from_json(loaded_model_json)
    model.load_weights("model/" + energy + "/model_dnn_tuned.h5")
    scaler = joblib.load('model/' + energy + '/scaler_X2.save')
    
    if database == 'structureDB':
        candidate = search_structure(formula)
    else:
        candidate = search_pubchem(formula)
    candidate_smiles = candidate['SMILES']
    data_X1, data_X2, known = new_data(mass, formula, ms_file, energy) 
    data_X2 = scaler.transform(data_X2)
    
    preds = model.predict([data_X1, data_X2])[:,0]
    neighbors = np.where(preds > thres)[0]
    if len(neighbors) < 1:
        print ('no neighbors are found')
    else:
        meta_neighbors = known.iloc[neighbors,]
        fig_neighbors = visualize_molecule(list(meta_neighbors['smiles']))
        meta_neighbors = pd.concat([pd.Series(preds[neighbors]), meta_neighbors.reset_index(drop=True), pd.Series(fig_neighbors)], axis=1)
        meta_neighbors = meta_neighbors.rename(index=str,columns={0:'similarity', 1:'Image'})
        smi_neighbors = meta_neighbors['smiles']
        
        scores = []
        for smi in candidate_smiles:
            this_sims = np.array([compare_structure(smi, x) for x in  smi_neighbors])
            scores.append(get_score(preds[neighbors], this_sims))
        
        fig_candidate = visualize_molecule(list(candidate['SMILES']))
        candidate = pd.concat([pd.Series(scores), candidate.reset_index(drop=True), pd.Series(fig_candidate)], axis=1)
        candidate = candidate.rename(index=str,columns={0:'Score', 1:'image'})
        order = np.argsort(-np.array(scores))
    
        return candidate.iloc[order,], meta_neighbors
    
    

if __name__ == '__main__':
    sim = compare_structure('CCC', 'CC')
    ms1 = read_ms('data/spectra/simulated_spectra/40V/C04640.csv')
    ms2 = read_ms('data/spectra/simulated_spectra/40V/C03373.csv')
    ms_vec1 = ms2vec(ms1)
    ms_vec2 = ms2vec(ms2)
    shifts = np.array(range(int(-ms_vec1.tocoo().shape[1]/2), int(ms_vec1.tocoo().shape[1]/2)))
    ccor = fft_ccor(ms_vec1.A.squeeze(), ms_vec2.A.squeeze()).A.squeeze()
    plt.plot(shifts, ccor)
    
    '''
    mass = 146.19
    formula = 'C6H14N2O2'
    spectrum_file = 'experiment/spectra/Lysine.csv'
    result, neighbors = run_one_example(mass, formula, spectrum_file, energy='40V', thres = 0.5, database='structureDB')
    '''