import os
import shutil
import numpy as np
import scipy.constants as sc

curdir = os.getcwd()
print(curdir)
cms = sc.c


#zstart = 3e-3
#zend = 8e-3
nsteps = 7


tstart = 0.0030*1e-9
tend = 0.0070*1e-9

ztot = tend*cms*1e3

print('tstart: ',tstart)
print('tend: ',tend)

nums = np.arange(nsteps)
tss = np.linspace(tstart,tend,nsteps)


#dirnames = [str(nsteps)+'steps/tstep_'+str(i) for i in nums]

print(nums,len(tss))

#for i in nums:
#    os.mkdir(dirnames[i])

for i in nums:
    print(i)
    with open(r'inputts.in','r') as file:
        data = file.read()
        data = data.replace('tdiff',str(tss[i]))
    with open(r'input'+str(i)+'.in','w') as file:
        file.write(data)
    with open(r'submit_metists.pbs','r') as file:
        data = file.read()
        data = data.replace('theory','theory'+str(i))
        data = data.replace('input.in','input'+str(i)+'.in')
    with open(r'submit_metis'+str(i)+'.pbs','w') as file:
        file.write(data)
    os.system(f'qsub submit_metis'+str(i)+'.pbs')























