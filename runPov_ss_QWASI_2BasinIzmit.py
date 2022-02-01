## This is a run for BETR-Reserach with annually changing env.-parameters ####
import os
import time
import BETRS
import importlib
importlib.reload(BETRS)
from BETRS import *
import pdb
import csv

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

# Open the csv file to write the results for all runs
f_ALL = open('Output/' + 'ssResults.csv','w',newline='')
header = ['Chem_Name',\
          'Central_Dissolved_Upper', 'Central_Dissolved_Lower',\
          'Eastern_Dissolved_Upper', 'Eastern_Dissolved_Lower',\
          'Central_Particle_Upper', 'Central_Particle_Lower',\
          'Eastern_Particle_Upper', 'Eastern_Particle_Lower',\
          'Central_Sediment', 'Eastern_Sediment']
writer = csv.writer(f_ALL)
writer.writerow(header)   
              

runID = ['NaP_ss', 'AcNP_ss', 'AcN_ss', 'Flo_ss', 'PhA_ss', 'Ant_ss', 'Flrn_ss',\
         'Pir_ss', 'BaA_ss', 'Kri_ss', 'BbF_ss', 'BkF_ss','BaP_ss', 'IP_ss',\
         'DahA_ss', 'BghiP_ss',\
         'PCB28_ss', 'PCB153_ss', 'BDE-47', 'BDE-99', 'BDE-153', '12347-PCDD',\
         'Dibenzo-p-furan'] # output names 
years = [list(range(1,5))]*len(runID)    # range of modeling run (years)

#emisdir = ['Emission_BDE209_10_comparts']  # emission inventory ('Emissions/annual/')
emisfile = ['Emission_QWASI_100.txt']*len(runID) # emission inventory ('Emissions/')
seasparfile = ['seasonal_parameters_QWASI_2BasinIzmit.txt']*len(runID)#  seasonally varying parameters ('Environment/)
constparfile = ['const_parameters_QWASI_2BasinIzmit.txt']*len(runID)  # seasonally constant parameters ('Environment/')
flowdirectory = ['QWASI_2BasinIzmit_Layered']*len(runID)  # flows in the atmosphere, ocean and fresh water ('Flows/)

chemdata = ['chemicals.txt']*len(runID)  # chemical properties ('Chemicals/')
chemnr = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]  # selection of chemical from chemical properties files
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
    m.solve_ss()
    if mkendfile == True:
        m.output_end_txt(filename = 'endstate.txt' + str(years[i][0]) + '.txt')
    result = m.ss_res
    m.ss_res = result
    m.output_ss_QWASI(filename = 'ss', units = ['mol_per_m3', 'Pa','ng_per_m3', 'diss_ng_per_L', 'part_ng_per_mg', 'sed_ng_per_g'], cpk = True)
    m.output_ss_txt(filename = 'ss.txt')

    if track_fluxes:
        m.flux_res = flux_rkg # experimenting to get the whole fluxes written, RKG, 04.07.2014
        m.output_fluxes(netcdf=True)
    if track_se:
        m.output_se(cpk=False, netcdf=True,scenario="air")

    t_e = time.time() # End time

    print('Simulation completed!')
    print(runID)
    print('TOTAL SIMULATION TIME: %f minutes = %f hours.' % ( (t_e - t_s) / 60.0 , (t_e - t_s) / 3600.0 ))
    
    ## Write results to csv file 
        
    f_output = open('Output/'+runID[i]+'/ss_out.cpk','rb')
    result = pickle.load(f_output)
    print('\nWriting the result for the run ' + str(runID[i]) + ' to ' + f_ALL.name)
    # Extract and write the results to f_ALL
    Central_Dissolved_Upper = result[1]['diss_ng_per_L'][0]
    Central_Dissolved_Lower = result[2]['diss_ng_per_L'][0]
    Eastern_Dissolved_Upper = result[1]['diss_ng_per_L'][1]
    Eastern_Dissolved_Lower = result[2]['diss_ng_per_L'][1]
    Central_Particle_Upper = result[1]['part_ng_per_mg'][0]
    Central_Particle_Lower = result[2]['part_ng_per_mg'][0]
    Eastern_Particle_Upper = result[1]['part_ng_per_mg'][1]
    Eastern_Particle_Lower = result[2]['part_ng_per_mg'][1]
    Central_Sediment = result[3]['sed_ng_per_g'][0]
    Eastern_Sediment = result[3]['sed_ng_per_g'][1]
    data = [m.chemdict['Name'],\
            Central_Dissolved_Upper, Central_Dissolved_Lower,\
            Eastern_Dissolved_Upper, Eastern_Dissolved_Lower,\
            Central_Particle_Upper, Central_Particle_Lower,\
            Eastern_Particle_Upper, Eastern_Particle_Lower,\
            Central_Sediment, Eastern_Sediment]
    writer.writerow(data)
    
    ## Write selected D-values to a csv file
    proc_name = 'betr_diff_loss'
    comp_from = 1
    comp_to = 1
    Dvals = m.Dproc[comp_from, comp_to, proc_name]
    f_Dvals = open('Output/' + runID[i] + '/' + runID[i] + '_' + proc_name + '.csv','w', newline='')
    print('\nWriting ' + proc_name + ' D values for the run ' + str(runID[i]) + ' to ' + f_Dvals.name)
    Dwriter = csv.writer(f_Dvals)
    header = ['month', 'Central Basin', 'Eastern Basin']
    Dwriter.writerow(header)
    for j in range(1,m.nots+1):
        data = [j, Dvals[0][j-1], Dvals[1][j-1]]
        Dwriter.writerow(data)
    f_Dvals.close()
    
    ## Write selected Z-values to a csv file
    comp_no = 1
    z_name = 'air'
    Zvals = m.zdict[comp_no][z_name]
    f_Zvals = open('Output/' + runID[i] + '/' + "Z_" + str(comp_no) + '_' + z_name + '.csv','w', newline='')
    print('\nWriting ' + z_name + ' Z values for the comp ' + str(comp_no) + ' for the run ' + str(runID[i]) + ' to ' + f_Zvals.name)
    Zwriter = csv.writer(f_Zvals)
    header = ['month', 'Central Basin', 'Eastern Basin']
    Zwriter.writerow(header)
    for j in range(1,m.nots+1):
        data = [j, Zvals[0][j-1], Zvals[1][j-1]]
        Zwriter.writerow(data)
    f_Zvals.close()

    # Persistence Calculations
    ssem = sum(m.get_ssEmission())
    TotMass = sum(m.ss_res)
    Pov = TotMass/ssem
    f_Pov = open('Output/' + runID[i] + '/ss_Pov.txt','w')
    f_Pov.write("Chem TotalEmissions_mol_per_h Total_Mass_mol  Pov_hours Pov_days\n")
    f_Pov.write("{} {:g} {:g} {:g} {:g}".format(m.chemdict['Name'],ssem,TotMass,Pov,Pov/24))
    print("Chem TotalEmissions_mol_per_h Total_Mass_mol  Pov_hours Pov_days")
    print("{}   {:g}                     {:g}            {:g}      {:g}".format(m.chemdict['Name'],ssem,TotMass,Pov,Pov/24))
    f_Pov.close()    

f_ALL.close()
            
    

