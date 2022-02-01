# -*- coding: utf-8 -*-
################################################################################
'''Module "mc.py", January 2022, is part of
BETR-QWASI-LU by Recep Kaya Goktas
based on BETR-Global 4.0 by Matt MacLeod

This module implements functions to conduct Monte Carlo simulations'''
################################################################################
from numpy import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def LN_value(m,k): # RKG, 20.01.2022
    """produces a random value from lognormal distribution, 
    m : the nominal value
    k : the dispersion factor"""
    
    CV = sqrt((np.exp((np.log(k)/1.96)**2) - 1))
    s = CV*m

    sigma = sqrt(np.log(1 + s**2/m**2))
    mu = log(m) - 0.5*sigma**2
    
    s = random.lognormal(mu, sigma)
    
    return s

def mu(m,k): # RKG, 20.01.2022
    """calculates parameter mu of a lognormal distribution
    m : the nominal value
    k : the dispersion factor"""
    
    CV = sqrt((exp((log(k)/1.96)**2) - 1))
    s = CV*m

    sigma = sqrt(log(1 + s**2/m**2))
    mu = log(m) - 0.5*sigma**2
    
    return mu

def sigma(m,k): # RKG, 20.01.2022
    """calculates parameter sigma of a lognormal distribution
    m : the nominal value
    k : the dispersion factor"""
    
    CV = sqrt((exp((log(k)/1.96)**2) - 1))
    s = CV*m

    sigma = sqrt(log(1 + s**2/m**2))
        
    return sigma

def sr_display(df):
    # calculates the spearman correlation coefficient
    # and displays it as a heatmap
    r = df.corr(method="spearman")
    plt.figure(figsize=(10,6))
    heatmap = sns.heatmap(df.corr(), vmin=-1, vmax=1, annot=True)
    plt.title("Spearman Rank Correlation")
    return r

def sr2_display(df,rows,columns,fn):
    # calculates the squared spearman correlation coefficient
    # and displays it as a heatmap
    r = df.corr(method="spearman")
    r = r**2
    plt.figure(figsize=(len(columns),len(rows)))
    heatmap = sns.heatmap(r.loc[rows,columns], vmin=0, vmax=1, annot=True, cmap="YlOrRd")
    #plt.title("Squared Spearman Rank Correlation Coefficients")
    plt.title("Spearman Sıralama Korelasyon Katsayılarının Karesi")
    plt.savefig(fn,bbox_inches="tight",dpi=300)
    return r