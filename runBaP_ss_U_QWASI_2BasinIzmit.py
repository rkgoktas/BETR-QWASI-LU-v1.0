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
monte_carlo = True # SWITCH ON/OFF MONTE CARLO SIMULATIONS, RKG, 20.01.2022
mc_iter = 10000    # assign the number of Monte Carlo iterations, RKG, 20.01.2022
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
          
runID = ['BaP_ss_U'] # output names 
years = [list(range(1,5))]*len(runID)    # range of modeling run (years)

#emisdir = ['Emission_BDE209_10_comparts']  # emission inventory ('Emissions/annual/')
emisfile = ['emission_BaP.txt']*len(runID) # emission inventory ('Emissions/')
seasparfile = ['seasonal_parameters_QWASI_2BasinIzmit.txt']*len(runID)#  seasonally varying parameters ('Environment/)
constparfile = ['const_parameters_QWASI_2BasinIzmit.txt']*len(runID)  # seasonally constant parameters ('Environment/')
u_envpar = {'fp1':3, 'fp2':3, 'partsett':3, 'seddep':3, 'focs7':1.5, 'fs7':3,\
            'rhos7':1.5, 'sedburial':3, 'sedresup':3, 'h7':2, 'tair2':1.1,\
            'Gup':3, 'Glow':3} # define the uncertain environmental parameters and their k-values
flowdirectory = ['QWASI_2BasinIzmit_Layered']*len(runID)  # flows in the atmosphere, ocean and fresh water ('Flows/)

chemdata = ['chemicals_PAH.txt']*len(runID)  # chemical properties ('Chemicals/')
chemnr = [13]*len(runID)  # selection of chemical from chemical properties files
u_chempar = {'logKaw':1.5, 'logKow':1.1, 'halflife_sediment':3} # define the uncertain chemical parameters and their k-values
compfile = ['compartments_QWASI.txt']*len(runID)   # compartments used in the model ('Environment/') 

procfile = ['processes_QWASI.txt']*len(runID)   # processes used in the model  ('Processes/')
contfile = ['control_QWASI.txt']*len(runID)      # some options ('Control/')
solvfile = ['solvparams_default.txt']*len(runID)    # options for ODE solver  ('Solver/')
mkendfile = False

#for v in [years, emisdir, seasparfile, chemdata, chemnr, constparfile, compfile, flowdirectory, procfile, contfile, solvfile]:
for v in [years, emisfile, seasparfile, chemdata, chemnr, constparfile, compfile, flowdirectory, procfile, contfile, solvfile]:
    if len(v) != len(runID):
        sys.exit('Warning: one of your input lists is not of same length as number of runs specified')        

input_params = {} # dictionary that will store all the input values of uncertain parameters, RKG, 22.01.2022
for p in u_envpar.keys(): # environmental parameters
    input_params[p] = []  # assign the keys with empty lists as values, RKG, 22.01.2022
for p in u_chempar.keys(): # chemical parameters
    input_params[p] = []    

output_params = {} # dictionary that will store all the output values from Monte Carlo simulations, RKG, 22.01.2022
output_comps = [1,2,3] # compartment numbers whose MC outputs will be stored, RKG, 22.01.2022
for c in output_comps:
    output_params[c]= [] # assign the keys with empty lists as values, RKG, 22.01.2022


 
## now run the model

for i in range(0, len(runID)):
    
    print(('\n\nStarting run ' + runID[i]))
    if monte_carlo == False:
        mc_iter = 1
    for iter in range(mc_iter): # Monte Carlo simulations with mc_iter iterations
        print('\iteration = ', iter+1)
        ## model first year and write temporary result to text file
        #print(('\n\nBETR run ' + runID[i] + ' for year ' + str(years[i][0])))
        m=Model(chemical = chemnr[i],
                run = runID[i],
                chemdb = chemdata[i],
                #seasonalparfile = os.path.join('annual', seasdir[i], str(years[i][0]) + '.txt'),
                seasonalparfile = seasparfile[i],    
                constantparfile = constparfile[i],
                u_envpar = u_envpar,
                u_chempar = u_chempar,
                compartmentfile = compfile[i], 
                #flowdir = os.path.join('annual', flowdirectory[i], str(years[i][0])),
                flowdir = os.path.join(flowdirectory[i]),
                processfile = procfile[i], 
                controlfile=contfile[i],
                monte_carlo = monte_carlo,
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
        m.output_ss_QWASI(filename = 'ss', units = ['mol_per_m3', 'Pa','ng_per_m3', 'mg_per_m3', 'mg_per_kg','diss_ng_per_L', 'part_ng_per_mg', 'sed_ng_per_g', 'sed_mg_per_kg'], cpk = True)
        #m.output_ss_txt(filename = 'ss.txt')
            
        if monte_carlo: # store input params for this Monte Carlo iteration
            for param_name in u_envpar.keys():
                input_params[param_name].append(mean(m.par[param_name],axis=1))
            for param_name in u_chempar.keys():
                input_params[param_name].append(m.chemdict[param_name])
                
        f_output = open('Output/'+runID[i]+'/ss_out.cpk','rb')
        result = pickle.load(f_output)
        
        output_params[1].append(result[1]['mg_per_m3'])
        output_params[2].append(result[2]['mg_per_m3'])
        output_params[3].append(result[3]['sed_mg_per_kg'])
    
    # Summarizing Monte Carlo simulation results
    if monte_carlo:
        # open the file to write MC results
        f_mc = open('Output/' + runID[i] + '/mc.txt','w')
        # Input parameters
        print("INPUT PARAMETERS")
        f_mc.write("INPUT PARAMETERS\n")
        for param_name in u_envpar.keys():
            print(param_name, "k =", u_envpar[param_name])
            f_mc.write(param_name + " k = " + str(u_envpar[param_name]) + "\n")
            print("95% Confidence Intervals:")
            f_mc.write("95% Confidence Intervals:\n")
            p025 = percentile(input_params[param_name],2.5,axis=0)
            p975 = percentile(input_params[param_name],97.5,axis=0)
            print("Central basin: {:.2e} - {:.2e}".format(p025[0],p975[0]))
            print("Eastern basin: {:.2e} - {:.2e}".format(p025[1],p975[1]))
            f_mc.write("Central basin: {:.2e} - {:.2e}\n".format(p025[0],p975[0]))
            f_mc.write("Eastern basin: {:.2e} - {:.2e}\n".format(p025[1],p975[1]))
        for param_name in u_chempar.keys():
            print(param_name, "k =", u_chempar[param_name])
            f_mc.write(param_name + " k = " + str(u_chempar[param_name]) + "\n")
            print("95% Confidence Intervals:")
            f_mc.write("95% Confidence Intervals:\n")
            p025 = percentile(input_params[param_name],2.5,axis=0)
            p975 = percentile(input_params[param_name],97.5,axis=0)
            print("Central basin: {:.2e} - {:.2e}".format(p025,p975))
            print("Eastern basin: {:.2e} - {:.2e}".format(p025,p975))
            f_mc.write("Central basin: {:.2e} - {:.2e}\n".format(p025,p975))
            f_mc.write("Eastern basin: {:.2e} - {:.2e}\n".format(p025,p975))    
       
        print("-----------------------------")
        f_mc.write("---------------------------------\n")
      
        print("OUTPUT PARAMETERS")
        f_mc.write("OUTPUT PARAMETERS\n")
        

        #print(" 2.5th percentile = {:.2e}\n 97.5th percentile = {:.2e}".format(p025[0], p975[0]))
        print("95% Confidence Intervals:")
        f_mc.write("95% Confidence Intervals:\n")
        
        p025 = percentile(output_params[1],2.5,axis=0)
        p975 = percentile(output_params[1],97.5,axis=0)
        print("Upper Ocean (mg/m^3):")
        print("Central basin: {:.2e} - {:.2e}\n".format(p025[0],p975[0]))
        print("Eastern basin: {:.2e} - {:.2e}".format(p025[1],p975[1]))
        f_mc.write("Upper Ocean (mg/m^3):\n")
        f_mc.write("Central basin: {:.2e} - {:.2e}".format(p025[0],p975[0]))
        f_mc.write("Eastern basin: {:.2e} - {:.2e}\n".format(p025[1],p975[1]))
        
        p025 = percentile(output_params[2],2.5,axis=0)
        p975 = percentile(output_params[2],97.5,axis=0)
        print("Lower Ocean (mg/m^3):")
        print("Central basin: {:.2e} - {:.2e}".format(p025[0],p975[0]))
        print("Eastern basin: {:.2e} - {:.2e}".format(p025[1],p975[1]))
        f_mc.write("Lower Ocean (mg/m^3):\n")
        f_mc.write("Central basin: {:.2e} - {:.2e}\n".format(p025[0],p975[0]))
        f_mc.write("Eastern basin: {:.2e} - {:.2e}\n".format(p025[1],p975[1]))
             
        p025 = percentile(output_params[3],2.5,axis=0)
        p975 = percentile(output_params[3],97.5,axis=0)
        print("Sediment (mg/kg):")
        print("Central basin: {:.2e} - {:.2e}".format(p025[0],p975[0]))
        print("Eastern basin: {:.2e} - {:.2e}".format(p025[1],p975[1]))
        f_mc.write("Sediment (mg/kg):\n")
        f_mc.write("Central basin: {:.2e} - {:.2e}\n".format(p025[0],p975[0]))
        f_mc.write("Eastern basin: {:.2e} - {:.2e}\n".format(p025[1],p975[1]))        
        
        f_mc.close()
        
        # Sensitivity Analysis (Spearman Rank Correlation Matrix & Heatmap)
        import pandas as pd
        
        inp_arr = [] # list to store input params as arrays
        # convert input_params into arrays
        for k in input_params.keys():
            inp_arr.append(array(input_params[k]))
     
        # create a dataframe for input parameters
        inp_dfs = []
        for cell in range(m.nocells):
            inp_s = [] # list to store input params as series
            for j,k in zip(range(len(inp_arr)), input_params.keys()):
                if type(inp_arr[j][0]) == numpy.ndarray:                  
                    inp_s.append(pd.Series(inp_arr[j][:,cell], name = k))
                else:
                    inp_s.append(pd.Series(inp_arr[j][:], name = k))
            inp_dfs.append(pd.concat(inp_s, axis=1))    
        
        out_arr = [] # list to store output params as arrays
        # convert output params into arrays
        for k in output_params.keys():
            out_arr.append(array(output_params[k]))
        
        # create a dataframe for output parameters
        out_dfs = []
        for cell in range(m.nocells):
            out_s = [] # list to store output params as series
            for j,k in zip(range(len(out_arr)), output_params.keys()):
                out_s.append(pd.Series(out_arr[j][:,cell], name = k))
            out_dfs.append(pd.concat(out_s, axis=1))
            
        # combine input and output dataframes for each cell
        all_dfs = []
        for cell in range(m.nocells):
            all_dfs.append(pd.concat([inp_dfs[cell], out_dfs[cell]], axis=1))
            
        for cell in range(m.nocells):
            all_dfs[cell].to_excel('Output/' + runID[i] + '/' + 'mc_cell_' + str(cell) + '.xlsx', index=False)
        
        # Calculate Spearman correlation coefficients
        r = []
        for cell in range(m.nocells):
            r.append(all_dfs[cell].corr(method="spearman"))
        
        r1 = sr_display(all_dfs[0])
        r2 = sr_display(all_dfs[1])
        
        r1 = sr2_display(all_dfs[0], list(output_params.keys()), list(input_params.keys()), 'Output/' + runID[i] + '/' + "sr2_merkez.png")
        r2 = sr2_display(all_dfs[1], list(output_params.keys()), list(input_params.keys()), 'Output/' + runID[i] + '/' + "sr2_doğu.png")
 
 
                
            
            
        

    t_e = time.time() # End time

    print('\nSimulation completed!')
    print(runID)
    print('TOTAL SIMULATION TIME: %f minutes = %f hours.' % ( (t_e - t_s) / 60.0 , (t_e - t_s) / 3600.0 ))
    

           
    

