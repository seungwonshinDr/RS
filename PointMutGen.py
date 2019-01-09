import numpy as np
import random

class Rseq:
    def __init__(self, Randombase, mutationpoint, mutantNumber):
        self.Randombase = Randombase
        self.mutationpoint = mutationpoint
        self.mutantNumber = mutantNumber


    def R(self):
        self.R = np.zeros((1, self.Randombase))
        for i in range(self.Randombase):
            self.R[0,i] = random.randrange(1,5)
        #print self.R
        return self.R


    def Mut(self):
        self.Mut = np.zeros((self.mutantNumber, self.Randombase))
        Mutpoint = np.zeros((1,self.mutationpoint))
        Mut_single = np.zeros((1,self.Randombase))


        for k in range(self.mutantNumber):
            list = np.array(range(self.Randombase))
            Mutpoint = np.array(random.sample(list,self.mutationpoint))

            for m in range(self.Randombase):
                Mut_single[0,m] = self.R[0,m]
            for j in range(self.mutationpoint):
                if Mut_single[0, Mutpoint[j]] == 1:
                    Mut_single[0, Mutpoint[j]] = np.random.choice([2, 3, 4])
                elif Mut_single[0, Mutpoint[j]] == 2:
                    Mut_single[0, Mutpoint[j]] = np.random.choice([1, 3, 4])
                elif Mut_single[0, Mutpoint[j]] == 3:
                    Mut_single[0, Mutpoint[j]] = np.random.choice([1, 2, 4])
                elif Mut_single[0, Mutpoint[j]] == 4:
                    Mut_single[0, Mutpoint[j]] = np.random.choice([1, 2, 3])
                self.Mut[k,:] = Mut_single[0,:]

        return self.Mut
