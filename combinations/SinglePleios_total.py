#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from itertools import combinations
#import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
import pickle
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import copy
from scipy.stats import sem, t
from scipy import mean
import re
import os

"""
CREATING DATABASE WITHDIFF PLEIOTROPIES
"""

#VARIABLES
PATH = "/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/Cancer/"
OUTPATH = sys.argv[2]
domains = {
'skin': list("skin disease"),
'skeletal': list("skeletal system disease"),
'respiratory': ["respiratory system disease", "Abnormality of the respiratory system"],
'reproductive': list("reproductive system disease"),
'cancer': ["neoplasm", "Meningioma"],
'immune': ["immune system disease"],
'infectious': ["infectious disease"],
'kidney': ["kidney disease", "rapid kidney function decline"],
'cardio': ["cardiovascular disease", "Abnormality of the cardiovascular system"],
'metabolic': ["digestive system disease", "metabolic disease", "Abnormality of metabolism/homeostasis"],
'headneck': ["Abnormality of head or neck", "head and neck disorder"],
'nervous': ["Abnormality of the nervous system", "nervous system disease"],
'eye': list("eye disease"),
'endocrine': list("endocrine system disease")}

catalog = pd.read_csv(PATH + "data/GWAS_Age_merged.csv", sep='\t', low_memory=False)#panda creation
LD = pd.read_csv("~/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/CEU_mod.csv", sep='\t')
SNPs_multiple_diseases = {}
out_file1 = OUTPATH + "SinglePleiotropiesOut.tsv"
out_file2 = OUTPATH + "SinglePleiosLoc.tsv"
age_threeshold = int(sys.argv[1])

## 1st strategy -> Single SNPs influencing more than one disease

def create_diseases_dictionary_for_each_SNP(catalog_fun, SNPs_dict):
    for i, row in catalog_fun.iterrows():
        if row['SNPS'] not in SNPs_dict.keys():
            SNPs_dict[row['SNPS']] = [row['DISEASE/TRAIT']]
        else:
            if row['DISEASE/TRAIT'] not in SNPs_dict[row['SNPS']]:
                SNPs_dict[row['SNPS']].append(row['DISEASE/TRAIT'])
    return SNPs_dict


def filter_pleiotropic_single_SNPs(SNPs_dict):
    pleiotropic_singles = {k:v for k,v in SNPs_dict.items() if len(SNPs_dict[k])>=2}
    return pleiotropic_singles

def get_subset_disease_Pleiotropies(catalog_fun, pleio_dict, domains):
    disease = []
    pleios = pd.DataFrame(columns=['SNP', 'Disease', 'Group',
    'PleioType', 'skin', 'skeletal', 'respiratory', 'reproductive', 'cancer',
    'immune', 'infectious', 'kidney', 'cardio', 'metabolic',
    'headneck', 'nervous', 'eye', 'endocrine'])
    for entry in pleio_dict:
        subset = catalog_fun.loc[catalog_fun['SNPS'] == entry]
        subset["PleioType"] = ""
        snps = subset["SNPS"].tolist()
        diseases = subset['Disease'].tolist()
        ages = subset["Onset"].tolist()
        groups = subset['Group'].tolist()
        current_dict = {}
        for key in domains:
            if any(value in groups for value in domains[key]):
                current_dict[key] = 1
            else:
                current_dict[key] = 0
        if subset.Risk_allele.nunique() == 1 :
            type = "Agonistic"
        else:
            type = "Antagonistic"
        pleiotype = define_early_late(age_threeshold, ages, type)
        disease.append(subset)
        pleios = pleios.append({'SNP': snps[0], 'Disease': set(diseases),
        'Group': set(groups), 'PleioType': pleiotype, 'skin': current_dict['skin'],
        'skeletal': current_dict['skeletal'], 'respiratory': current_dict['respiratory'],
        'reproductive': current_dict['reproductive'], 'cancer': current_dict['cancer'],
        'immune': current_dict['immune'], 'infectious': current_dict['infectious'],
        'kidney': current_dict['kidney'], 'cardio': current_dict['cardio'],
        'metabolic': current_dict['metabolic'], 'headneck': current_dict['headneck'],
        'nervous': current_dict['nervous'], 'eye': current_dict['eye'],
        'endocrine': current_dict['endocrine']},
        ignore_index=True)
    disease = pd.concat(disease)
    print(pleios)
    return disease, pleios

def define_early_late(threeshold, ages, type):
    if all(x < threeshold for x in ages):
        outcome = "Early-Early " + type
    elif all(x > threeshold for x in ages):
        outcome = "Late-Late " + type
    else:
        outcome = "Early-Late " + type
    return(outcome)


"""
def get_subset_Cancer_Pleiotropies(table, array):
    for i, row in table.iterrows():
        if row['Group'] == "neoplasm" or row['Group'] == "Meningioma":
            array.append(row['SNP'])
    table_subset1 = table[table.SNP.isin(array)]
    table_subset1 = table_subset1.sort_values(by=['SNP'])
    return table_subset1

def get_subset_immunitary_Pleiotropies(table, array):
    for i, row in table.iterrows():
        if row['Group'] == "immune system disease":
            array.append(row['SNP'])
    table_subset2 = table[table.SNP.isin(array)]
    table_subset2 = table_subset2.sort_values(by=['SNP'])
    return table_subset2

def agonistic_or_antagonistic(table):
    Agonistic = 0
    Antagonistic = 0
    for snp in table.SNP.unique().tolist():
        df = table.loc[table['SNP'] == snp]
        if len(df.RiskAll.unique().tolist()) == 1:
            Agonistic += 1
        else:
            Antagonistic +=1
    #while table.loc[df['column_name'] == some_value]
    return Agonistic, Antagonistic
"""


#def sadasfsafsaf():

#MAIN

SNPs_multiple_diseases = create_diseases_dictionary_for_each_SNP(catalog, SNPs_multiple_diseases)
#print(SNPs_multiple_diseases)
pleiotropic_singles = filter_pleiotropic_single_SNPs(SNPs_multiple_diseases)
#for key in pleiotropic_singles.keys():
    #print(key, pleiotropic_singles[key])
pleio_dataset, pleios = get_subset_disease_Pleiotropies(catalog, pleiotropic_singles, domains)
pleio_dataset.to_csv(out_file1, sep='\t')
pleios.to_csv(out_file2, sep='\t')

"""
SinglePleio_subset1 = get_subset_infectious_Pleiotropies(SinglePleio, SingleInfectSNPs)
SinglePleio_subset2 = get_subset_immunitary_Pleiotropies(SinglePleio, SingleInfectSNPs)
SinglePleio_subset2.to_csv(out_file, sep='\t')

Agonistic, Antagonistic = agonistic_or_antagonistic(SinglePleio_subset2)
print(Agonistic, Antagonistic)
"""
