################################################################################
'''Module "tempcorrect", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module calculates temperature corrected fields of Koa,Kaw,Kow and
degradation half-lifes for the selected chemical in all compartments at
all timesteps'''
################################################################################

from numpy import *
import inspect
from globalz import *
import pdb

def tempcorrect(par, compdict, chem):
    '''Creates temperature corrected partitioning properties and half-lives
    for the selectedd chemical in all regions in all compartments at all
    timesteps.
    <compdict> is the dicitionary containing the compartments used.
    <par> is the cells x timesteps record-array of all cell-parameters.
    <chem> is a dictionary containing (standard) chemical properties of
    one substance.

    Returns a dictionary with compartment IDs as keys,
    and [cells] x [timesteps] record-arrays with fields
    (k_reac, Kaw, Koa, Kow), containing reaction rate-constant k_reac (1/h)
    and phase partioning coefficients Kaw, Koa, Kow'''

    ## ATT : half-lifes in the input should be changed to rate-constants
    ##       then perfect persistence would be much easier to implement

    chem_spatial_dict={}
    for k in list(compdict.keys()):
        #print("in tempcorrect: k = ", k) # RKG
        ## initialize record-array
        pa=zeros(par.shape, dtype=dtype([('k_reac','f8'), ('Kaw','f8'),
                                              ('Koa','f8'),('Kow','f8')]))
        tempfield=par[compdict[k]['temp_variable']]
        #print("in tempcorrect: tempfield = ", tempfield) # RKG, 03.12.2021
        t0=chem['T0']
        #print("in tempcorrect: t0 = ", t0) # RKG, 03.12.2021
        ## reaction rate (1/h)
        kreac0=log(2)/chem[compdict[k]['halflife_variable']]
        #print("in tempcorrect: kreac0 = ", kreac0) # RKG, 03.12.2021
        ea=chem[compdict[k]['EA_variable']] # J/mol  
        #print("in tempcorrect: ea = ", type(ea), ea) # RKG, 03.12.2021
        #print("in tempcorrect: R = ", R) # RKG, 03.12.2021
        pa['k_reac']= kreac0*exp(ea/R*(1/t0-1/tempfield))
        ## partitioning coefficients
        dUoa=chem['DUoa'] # J/mol
        dUow=chem['DUow'] # J/mol
        dUaw=dUow-dUoa # J/mol
        Kaw0=10**chem['logKaw']
        Koa0=10**chem['logKoa']
        Kow0=10**chem['logKow']
        pa['Kaw']=Kaw0*exp(dUaw/R*(1/t0-1/tempfield)) 
        pa['Koa']=Koa0*exp(dUoa/R*(1/t0-1/tempfield))
        pa['Kow']=Kow0*exp(dUow/R*(1/t0-1/tempfield))
        chem_spatial_dict[k]=pa

    ########################
    ### OH-radical correction is commented out since BETR-QWASI does not have an
    ### air compartment - RKG, 01.12.2021
    ## weight degradation in air according to OH-radical concentration
    ## ATT : This special treatment of air1 and air2 should be generalized
    ## somehow.
    #if 1 in list(chem_spatial_dict.keys()):
    #    chem_spatial_dict[1]['k_reac']*=par['OHair1'] # upper air
    #if 2 in list(chem_spatial_dict.keys()):
    #    chem_spatial_dict[2]['k_reac']*=par['OHair2'] # lower air
    #if 8 in list(chem_spatial_dict.keys()):
    #    chem_spatial_dict[8]['k_reac']*=par['OHair1'] # upper air coarse aerosol
    #if 9 in list(chem_spatial_dict.keys()):
    #    chem_spatial_dict[9]['k_reac']*=par['OHair2'] # lower air fine aerosol 
    #if 10 in list(chem_spatial_dict.keys()):
    #    chem_spatial_dict[10]['k_reac']*=par['OHair2'] # lower air coarse aerosol
    ####################### end of the BETR-QWASI modification - RKG, 01.12.2021
    return(chem_spatial_dict)
