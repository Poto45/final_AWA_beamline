import numpy as np
import os

curdir = os.getcwd()

nsteps = 7



nums = np.arange(nsteps)


for i in nums:
    print(i)
    os.chdir('theory'+str(i))
    os.rename('input'+str(i)+'.stat','input.stat')
    os.rename('input'+str(i)+'.h5','input.h5')
    os.rename('input'+str(i)+'_Monitors.stat','input_Monitors.stat')
    os.rename('input'+str(i)+'.in','input.in')
    os.chdir('data')
    os.rename('input'+str(i)+'_DesignPath.dat','input_DesignPath.dat')
    os.rename('input'+str(i)+'_ElementPositions.py','input_ElementPositions.py')
    os.rename('input'+str(i)+'_ElementPositions.txt','input_ElementPositions.txt')
    os.rename('input'+str(i)+'_ElementPositions.sdds','input_ElementPositions.sdds')
    os.rename('input'+str(i)+'_3D.opal','input_3D.opal')
    os.rename('input'+str(i)+'_DISTFT.dat','input_DISTFT.dat')
  

    os.chdir(curdir)









