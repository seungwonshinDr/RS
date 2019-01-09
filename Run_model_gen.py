import numpy as np
import PointMutGen as Rseq
import datetime
import time

starttime = time.time()
print 'Start calcultion at', datetime.datetime.now()

filename = 'Models'
Randombase = 8 #Total number of bases
mutantNumber = 5 #Number of Mutants


Rseqtime = time.time()
print 'Calculate Random clusters including Target and Mutants', datetime.datetime.now()
for i in range(1):
    title = Rseq.Rseq(Randombase, mutationpoint, mutantNumber)
    R = title.R()
    Mut = title.Mut()

    np.savetxt(filename+'_Target',R,delimiter = ',')
    np.savetxt(filename+'_Mutants',Mut,delimiter = ',')

print R
print Mut

Rseqdone = time.time()
print 'Calculation is done. It took', Rseqdone-Rseqtime, 'seconds'

