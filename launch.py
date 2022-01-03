#copyright 2021, Auton Systems LLC

import argparse

parser = argparse.ArgumentParser(description='run simulations')
parser.add_argument('-N','--N_Samples', type=int, default=100, help='Number of samples to run')
parser.add_argument('-n','--N_Workers', type=int, default=2, help='Number of workers to use')
parser.add_argument('-t','--TimeLimit', type=int, default=1800, help='Timeout (s) for each sample')
parser.add_argument('-p','--ClassPath', type=str, help='Path to MutAntiGen jar files')
args = parser.parse_args()

import numpy as np
from joblib import Parallel, delayed
import os
#import lhsmdu

launch_cmd = '''timeout '''+str(args.TimeLimit)+'''s java -Xmx6G -XX:ParallelGCThreads=8 -cp "'''+args.ClassPath+'''" Mutantigen'''

DefaultParams = {"burnin": 0 ,
"endDay": 1095,#(1095=3yrs)    #3650 ,
"printStep": 7 ,
"tipSamplingStartDay": 2000 ,
"tipSamplingEndDay": 3450 ,
"tipSamplingRate": 0.00125 ,
"tipSamplesPerDeme": 600 ,
"tipSamplingProportional": True ,
"treeProportion": 0.01 ,
"diversitySamplingCount": 1000 ,
"netauWindow": 100 ,
"repeatSim": False ,
"immunityReconstruction": False ,
"memoryProfiling": False ,
"yearsFromMK": 1.0 ,
"pcaSamples": False ,
"detailedOutput": False ,
"restartFromCheckpoint": False ,
"hostImmuneHistorySampleCount": 10000 ,
"fitSampleCount": 100 ,
"printFitSamplesStep": 100000000 ,
"fitnessLogValue": False ,
"demeCount": 1 ,
"demeNames": ["tropics"] ,
"initialNs": [40000000] ,
"birthRate": 0.000091 ,
"deathRate": 0.000091 ,
"swapDemography": True ,
"initialI": 7406 ,
"initialDeme": 1 ,
"initialPrR": 0.5088 ,
"beta": 0.5627 ,
"nu": 0.25 ,
"betweenDemePro": 0.0000 ,
"externalMigration": 200.0 ,
"transcendental": False ,
"immunityLoss": 0.0 ,
"initialPrT": 0.0 ,
"backgroundImmunity": False ,
"backgroundDistance": 0.2 ,
"demeBaselines": [1.] ,
"demeAmplitudes": [0.0] ,
"demeOffsets": [0.] ,
"phenotypeSpace": "mutLoad" ,
"lambda": 0.10 ,
"mutCost": 0.008 ,
"probLethal": 0.0 ,
"epsilon": 0.16 ,
"epsilonSlope": 0.0 ,
"lambdaAntigenic": 0.00075 ,
"meanAntigenicSize": 0.012 ,
"antigenicGammaShape": 2.0 ,
"thresholdAntigenicSize": 0.012 ,
"thresholdMutationClade": 25 ,
"CladeSamples": 4000,
"antigenicEvoStartDay": 0 ,
"cleanUpDistance": 0.2 ,
"demoNoiseScaler": 0.0 ,
"muPhenotype": 0.0 ,
"smithConversion": 0.1 ,
"homologousImmunity": 0.95 ,
"initialTraitA": -6. ,
"meanStep": 0.3 ,
"sdStep": 0.3 ,
"mut2D": False }

def WriteParamsFile(edits={}):
    f = open("parameters_load.yml","w")
    for k,v in DefaultParams.items():
        if k in edits:
            v = edits[k]
        if isinstance(v,str):
            vstr = '"'+v+'"'
        elif isinstance(v,bool):
            vstr = str(v).lower()
        else:
            vstr = str(v)
        f.write(k+": "+vstr+'\n')
    f.close()

def Launch(dat):
    try:
      id_, wd, params = dat
      os.chdir(wd)
      os.system("mkdir job"+str(id_))
      os.chdir(wd+"/job"+str(id_))
      WriteParamsFile(params)
      os.system(launch_cmd)
    except Exception as e:
      print(str(e))

sample_pars = [
("initialNs",6,6) , # constant
("initialPrR",0.1,0.1) # constant
]

#sample_pars = [
#("initialNs",4,8) ,
#("initialPrR",0.,0.5) ,
#("beta",0.3,0.6) ,
#("nu",0.15,0.25) ,
#("lambda",0.05,0.25) ,
#("mutCost",0.006,0.01) ,
#("epsilon",0.1,0.2) ,
#("lambdaAntigenic",0.00075,0.001)
#]

p0 = np.array([tup[1] for tup in sample_pars])
p1 = np.array([tup[2] for tup in sample_pars])

print("Sampling parameters...")
#l = lhsmdu.sample(len(sample_pars),args.N_Samples)
#l = np.array(l)
l = np.random.random((len(sample_pars),args.N_Samples))
d = p1-p0
#for i in range(len(d)):
#    assert d[i]>0

l = d[:,None]*l
l = l+p0[:,None]

params = [dict( [(sample_pars[i][0],l[i,k]) for i in range(len(sample_pars))] ) for k in range(l.shape[1])]
for p in params: 
    p["initialI"] = int( np.sqrt( 10**p["initialNs"] ) )
    p["initialNs"] = [int( 10**p["initialNs"] )]


print("Prepping jobs...")
wd = os.getcwd()
jobs = [(i,wd,params[i]) for i in range(len(params))]

# run jobs
Parallel(n_jobs=args.N_Workers)(delayed(Launch)(tup) for tup in jobs)

print("Done.")

