################################################################################
'''Module "mkflowsD", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module calculates D-values for inter-cell transport processes of a
particular model parametrization.'''
################################################################################
from numpy import *
import inspect
from globalz import *
import sys
import copy

def mkflowD(model):
    fdict=copy.deepcopy(model.flowdict)
    print("in mkflowD: fdict = ", fdict) # RKG, 03.12.2001
    for f in fdict.keys():
        print("in mkflowD: f = ", f) # RKG, 07.12.2001
        print("in mkflowD: fdict[f] = ", type(fdict[f]), shape(fdict[f]), fdict[f]) # RKG, 07.12.2001
        try: ## try-except added by RKG, 07.12.2021
            fromcells=fdict[f][:,0].astype(int)
            zvals=model.zdict[f[0]]['bulk'][fromcells-1,:]
            fdict[f][:,2:]=fdict[f][:,2:]*zvals
        except:
            pass # if there is only 1 CELL, RKG, 07.12.2021
        
    return(fdict)
