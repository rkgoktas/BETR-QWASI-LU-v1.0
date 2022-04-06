## This is a run for BETR-Reserach with annually changing env.-parameters ####
import os
import time
import BETRS
import importlib
importlib.reload(BETRS)
from BETRS import *
import pdb

t_s = time.time() # Start time

"""
SET OPTIONS FOR THE FAST SOLVER AND THE FLUX_INTEGRATION AND
PRIMARY/SECONDARY_EMISSION MODE
"""

use_odespy = False       #SWITCH ON/OFF FAST SOLVER
track_fluxes = False    #SWITCH ON/OFF FLUX INTEGRATION
track_se = False         #SWICH ON/OFF TRACKING OF SECONDARY EMISSIONS
use_correction = False   #SWICH ON/OFF flow with correction

## options
"""
Change to current Directory to ensure that relative paths are set correctly

Tends to cause problem under Windows otherwise

"""

abspath=os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

runID = ['BaP_dyn_SedCap_100E10W'] # output names 
years = [list(range(1,16))]*len(runID)    # range of modeling run (years)

#emisdir = ['Emission_BDE209_10_comparts']  # emission inventory ('Emissions/annual/')
emisfile = ['Emission_BaP_const.txt']*len(runID) # emission inventory ('Emissions/')
seasparfile = ['seasonal_parameters_QWASI_2BasinIzmit_const.txt']*len(runID)#  seasonally varying parameters ('Environment/)
constparfile = ['const_parameters_QWASI_2BasinIzmit.txt']*len(runID)  # seasonally constant parameters ('Environment/')
flowdirectory = ['QWASI_2BasinIzmit_Layered_ConstFlow']*len(runID)  # flows in the atmosphere, ocean and fresh water ('Flows/)

chemdata = ['chemicals_PAH.txt']*len(runID)  # chemical properties ('Chemicals/')
chemnr = [13]*len(runID)  # selection of chemical from chemical properties files
compfile = ['compartments_QWASI.txt']*len(runID)   # compartments used in the model ('Environment/') 

procfile = ['processes_QWASI.txt']*len(runID)   # processes used in the model  ('Processes/')
contfile = ['control_QWASI.txt']*len(runID)      # some options ('Control/')
solvfile = ['solvparams_default.txt']*len(runID)    # options for ODE solver  ('Solver/')
mkendfile = False


#for v in [years, emisdir, seasparfile, chemdata, chemnr, constparfile, compfile, flowdirectory, procfile, contfile, solvfile]:
for v in [years, emisfile, seasparfile, chemdata, chemnr, constparfile, compfile, flowdirectory, procfile, contfile, solvfile]:
    if len(v) != len(runID):
        sys.exit('Warning: one of your input lists is not of same length as number of runs specified')        
    

## now run the model

for i in range(0, len(runID)):
    print(('\n\nStarting run ' + runID[i]))
    ## model first year and write temporary result to text file
    print(('\n\nBETR run ' + runID[i] + ' for year ' + str(years[i][0])))
    m=Model(chemical = chemnr[i],
            run = runID[i],
            chemdb = chemdata[i],
            #seasonalparfile = os.path.join('annual', seasdir[i], str(years[i][0]) + '.txt'),
            seasonalparfile = seasparfile[i],    
            constantparfile = constparfile[i], 
            compartmentfile = compfile[i], 
            #flowdir = os.path.join('annual', flowdirectory[i], str(years[i][0])),
            flowdir = os.path.join(flowdirectory[i]),
            processfile = procfile[i], 
            controlfile=contfile[i],
            use_correction = use_correction, 
            track_flows = (track_fluxes or track_se),
            )
    #m.update_emissions(os.path.join('annual', emisdir[i], str(years[i][0]) + '.txt'))
    m.update_emissions(emisfile[i])    
    m.update_solver(solvparamfile = solvfile[i], initfile = 'initial_BaP_from_ss_constEm.txt')
    m.solve_dyn(1,use_odespy=use_odespy)

    if mkendfile == True:
        m.output_end_txt(filename = 'endstate.txt' + str(years[i][0]) + '.txt')
    result = m.dyn_res                                                          # save model result from year 1
    flux_rkg = m.flux_res # experimenting to get the whole fluxes written, RKG, 04.07.2014
    m.output_dyn_tmp(filename = runID[i] + '_tmp.txt')                             # save end state of year 1 as text file
    

    ## model year 2 to n
    for y in years[i][1:]:
        del(m)
        print(('\n\nBETR run ' + runID[i] + ' for year ' + str(y)))
        m=Model(chemical = chemnr[i], 
            run = runID[i],
            chemdb = chemdata[i],
            #seasonalparfile = os.path.join('annual', seasdir[i], str(y) + '.txt'),
            seasonalparfile = seasparfile[i],
            constantparfile = constparfile[i], 
            compartmentfile = compfile[i], 
            #flowdir = os.path.join('annual', flowdirectory[i], str(y)),
            flowdir = os.path.join(flowdirectory[i]),
            processfile = procfile[i], 
            controlfile=contfile[i],
            use_correction = use_correction,
            track_flows = (track_fluxes or track_se),
            )       
        #m.update_emissions(os.path.join('annual', emisdir[i], str(y) + '.txt'))
        m.update_emissions(emisfile[i]) 
        if y == 6:
            m.update_solver(solvparamfile = solvfile[i], initfile = 'BaP_endstate5_cap_100E10W.txt') # update solver with end state of capping scenario
        else:
            m.update_solver(solvparamfile = solvfile[i], initfile = runID[i] + '_tmp.txt') # update solver with end state of year - 1 from text file
        m.solve_dyn(1,use_odespy=use_odespy)
        if mkendfile == True:
            m.output_end_txt(filename = 'endstate.' + str(y) + '.txt')
        if y == years[i][len(years[i])-1]: # If this is the last year, save the end state.
            m.output_end_txt(filename = 'endstate.' + str(y) + '.txt')
        result = hstack((result, m.dyn_res[:, 1:13]))                           # save model result from year y
        flux_rkg = hstack((flux_rkg, m.flux_res)) # experimenting to get the whole fluxes written, RKG, 04.07.2014
        m.output_dyn_tmp(filename = runID[i] + '_tmp.txt')                         # save end state of year y as text file


    ## nc-file output
    m.dyn_res = result
    #m.output_dyn(filename = 'OUT', units = ['kg', 'kg_per_m3', 'Pa', 'mol', 'mol_per_m3'], netcdf = True, cpk = False)
    m.output_dyn(filename = 'OUT', units = ['mol', 'Pa', 'mol_per_m3', 'ng_per_m3'], netcdf = False, cpk = False)
    os.remove(os.path.join('Solver', runID[i] + '_tmp.txt'))                    # remove temporary file
    if track_fluxes:
        m.flux_res = flux_rkg # experimenting to get the whole fluxes written, RKG, 04.07.2014
        m.output_fluxes(netcdf=True)
    if track_se:
        m.output_se(cpk=False, netcdf=True,scenario="air")

    t_e = time.time() # End time

    print('Simulation completed!')
    print(runID)
    print('TOTAL SIMULATION TIME: %f minutes = %f hours.' % ( (t_e - t_s) / 60.0 , (t_e - t_s) / 3600.0 ))
    
    # Plotting the results
    m.output_dyn_QWASI(filename = 'OUT', units = ['mol', 'Pa', 'mol_per_m3', 'ng_per_m3', 'mg_per_m3', 'mg_per_kg', 'sed_ng_per_g', 'sed_mg_per_kg'] )  
    
