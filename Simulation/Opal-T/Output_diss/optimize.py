import os
import uuid
import numpy as np
import h5py
import opt_func as optH
from shutil import copyfile, rmtree
from scipy.stats import powerlaw, skew
from opal.opal import load_dataset


abs_path = os.getcwd()
ref_file = 'input.ref.in'
input_list = ['note.txt'] # This includes the initial particle distribution as well  
# note is an unused editable file to log the case
input_list.append(ref_file)

zstop = 12.00  # nominal location where the simulation ends  

Np = 30000  # number of particles

# Ranges of the parameters for adjusting the quadrupole magnet
# Myranges: TDC #2 kick (ScalTDC2), Q2 strength (T/m), Q2 position (m)

myRanges = [(40,50), 
            (-20,20),
            (0,550),
            (0,440),
            (0.8,1.2),
            (10,12),
            (-20,20),
            (10,12),
            (-20,20),
            (2e-3,15e-3)
           ]

NDIM = len(myRanges)
BOUND_UP = 1
BOUND_LOW = 0


def writeOPAL(X):
    '''
    Prepare input file for OPAL
    '''
    X = optH.parametersRestore(X, myRanges)

    parameters = {
                 '__gun_inj_amp' : X[0],
                 '__gun_inj_phi' : X[1],
                 '__Ifocu' : X[2],
                 '__Imain' : X[3],
                 '__Ibuck_scale' : X[4],
                 '__L1_Amp' : X[5],
                 '__L1_phi' : X[6],
                 '__L2_Amp' : X[7],
                 '__L2_phi' : X[8],
                 '__radiusLaser' : X[9]
                 }



    # Search and replace. Ref file should be prepared.
    f = open(ref_file,'r')
    deck = f.read()
    f.close()
    for i in parameters:
        deck = deck.replace(i, str(parameters[i]))
    f = open('input.in','w')
    f.write(deck)
    f.close()
#    os.system("python /lcrc/project/Bright-Beams/seongyeol/FY23/RFBT/2023_June_Study/1_140mT_Beta1.76m_Drivelinac_MOGA/formatter.py")


def removeDirs(hof):
    '''Remove all directories under abs_path'''
    os.chdir(abs_path)
    exclude_list = ['input','logs','__pycache__','archive','minW1','minEm']
    hof_list = [g.indID for g in hof]
    exclude_list = exclude_list + hof_list
    dir_list = [f for f in os.listdir('.') if os.path.isdir(f) and f not in exclude_list]
    for ind in dir_list:
        rmtree(ind)
                

####################################################################
###                      objective function                      ###
####################################################################

# cma always minimize the objective function
# use fScale to control maximize / minimize
# negative weight means that it is minimization problem
badFitness = (1000, 1000, 1000, 1000) #, 1000, 1000)
myWeights = (-1.0, -1.0 , -1.0, -1.0) # 
# want to minimize epsilon and minimize s_rms
# is -1 to maximize and 1 to minimize?
fScale = (1, 1, 1, 1)


def softConstraint (x,constraintValue,weight,softness):
   '''
   sets an upper boudary soft contraint
   '''
   return (x+weight*(1+np.tanh((x-constraintValue)/softness)))

def narrowContraint (x,constraintValue, weight, softness):
   '''
   sets a target contraint (use a band-pass filter)
   '''
   n=4
   return (weight*(1-np.exp(-((x-constraintValue)/softness)**n)))


def fitnessFunc(X):
    os.chdir(abs_path)
    dirName = X.indID
    os.mkdir(dirName)
    os.chdir(dirName)

    for f in input_list:
        copyfile(F'{abs_path}/{f}', F'{abs_path}/{dirName}/{f}')
    
    writeOPAL(np.asarray(X))
#    os.system(f"/opt/metis/el8/contrib/opal/OPAL-2022.1/bin/opal input.in > log")
    os.system(f"module load opal/opal-2022.1; opal input.in > opal.log.su")
#     os.system(f"module load opal/opal-2022.1; opal input.in > opal.log.su")
    
    # OPAL results
    # File name
    opalstatout = "input.stat"
    
    if not os.path.isfile(opalstatout):
        return [tuple(np.divide(badFitness, fScale)), (np.nan,)]

    # contraint / run completion
    ds = load_dataset('./', fname=opalstatout)
    s  = ds.getData('s')
#    if s[-1]<i0.99*zstop:     # simulation did not finish 
#        return [tuple(np.divide(badFitness, fScale)), (np.nan,)]
#
# record beam parameters:
#
    Nalive  = ds.getData('numParticles')[-1] 
    kinetic = ds.getData('energy')[-1]
    emit_x  = ds.getData('emit_x')[-1]
    emit_y  = ds.getData('emit_y')[-1]
    emit_s  = ds.getData('emit_s')[-1]
    rms_x   = ds.getData('rms_x')[-1]    
    rms_y   = ds.getData('rms_y')[-1]
    rms_s   = ds.getData('rms_s')[-1] 
    rms_ps  = ds.getData('rms_ps')[-1]

    # save other parameters 
    otherParams=(rms_x, rms_y, rms_s, rms_ps, emit_x, emit_y, emit_s, kinetic, Nalive)

    
    if Nalive < Np*0.8: # less tha 80% of intial num of macroparticle
        return [tuple(np.divide(badFitness, fScale)), otherParams]

#    betax = np.divide(rms_x,np.divide(emit_x,rms_ps))
#    if betax > 75:
#	return [tuble(np.divide(badFitness, fScale)), otherParams]

    #contraints
    Nfit = Np-Nalive
    Kfit = narrowContraint(kinetic, 40, 100, 10)

    #build fitness function:
    fitness = (emit_x, rms_x, rms_s, Nfit) # , Nfit, Kfit) #emit_x, emit_s, Nfit, Kfit)

    fitness = tuple(np.divide(fitness, fScale))


    # combine fitness and other parameters
    allFitness = [fitness, otherParams]


    return allFitness



def main():
    start_point_rd = np.random.uniform(0,1,NDIM)
    # Initial random start
    start_point = [start_point_rd[0]*myRanges[0][1], 
                   start_point_rd[1]*myRanges[1][1]                   ]

    dirName = uuid.uuid4().hex[:16].upper()
    dirName = "copies"
    os.chdir(abs_path)
    os.mkdir(dirName)
    os.chdir(dirName)
    # Prepare input and run OPAL
    for f in input_list:
        copyfile(F'{abs_path}/{f}', F'{abs_path}/{dirName}/{f}')
    writeOPAL(start_point)

if __name__ == "__main__":
    main()

