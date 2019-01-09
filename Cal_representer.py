import numpy as np
import sympy
from sympy import Symbol

class gen:
    def __init__(self, strand, conc, Temp):
        self.strand = strand
        self.conc = conc
        self.Temp = Temp
        self.fitness = 0

    def rep(self):
        rep_strand = np.zeros((1,np.shape(self.strand)[1]))
        for i in range(np.shape(rep_strand)[1]):
            if self.strand[0,i] == 1:
                rep_strand[0,i] = 2
            elif self.strand[0,i] == 2:
                rep_strand[0,i] = 1
            elif self.strand[0,i] == 3:
                rep_strand[0,i] = 4
            else:
                rep_strand[0,i] = 3

        complex = np.concatenate((rep_strand,self.strand), axis = 0)
        #print complex

        hyb_code = np.zeros((1,np.shape(complex)[1]-1))
        for i in range(np.shape(complex)[1]-1):
            hyb_code[0,i] = complex[0,i] * 1000 + complex[0,i+1] * 100 +complex[1,i] * 10 + complex[1,i+1]

        #print hyb_code

        for i in range(np.shape(hyb_code)[1]):
            if hyb_code[0,i] == 1122:
                self.fitness += -7.9 - self.Temp * (-22.2)
            elif hyb_code[0,i] == 1221:
                self.fitness += -7.2 - self.Temp * (-20.4)
            elif hyb_code[0,i] == 1324:
                self.fitness += -7.8 - self.Temp * (-21)
            elif hyb_code[0,i] == 1423:
                self.fitness += -8.4 - self.Temp * (-22.4)
            elif hyb_code[0,i] == 2112:
                self.fitness += -7.2 - self.Temp * (-21.3)
            elif hyb_code[0,i] == 2211:
                self.fitness += -7.9 - self.Temp * (-22.2)
            elif hyb_code[0,i] == 2314:
                self.fitness += -8.5 - self.Temp * (-22.7)
            elif hyb_code[0,i] == 2413:
                self.fitness += -8.2 - self.Temp * (-22.2)
            elif hyb_code[0,i] == 3142:
                self.fitness += -8.2 - self.Temp * (-22.2)
            elif hyb_code[0,i] == 3241:
                self.fitness += -8.4 - self.Temp * (-22.4)
            elif hyb_code[0,i] == 3344:
                self.fitness += -8 - self.Temp * (-19.9)
            elif hyb_code[0,i] == 3443:
                self.fitness += -9.8 - self.Temp * (-24.4)
            elif hyb_code[0,i] == 4132:
                self.fitness += -8.5 - self.Temp * (-22.7)
            elif hyb_code[0,i] == 4231:
                self.fitness += -7.8 - self.Temp * (-21)
            elif hyb_code[0,i] == 4334:
                self.fitness += -10.6 - self.Temp * (-27.2)
            elif hyb_code[0,i] == 4433:
                self.fitness += -8 - self.Temp * (-19.9)
            else:
                self.fitness += 0
        #print self.fitness
        return rep_strand

    def concentration(self):
        dG = self.fitness
        K = sympy.exp(dG * 4184 / (-8.314) / 300)

        x = Symbol('x')
        conc_cand = np.array(sympy.solve((K*(x-self.conc)*(x-self.conc)-self.conc*(2*x-self.conc))))
        #print conc_cand
        conc = conc_cand.max()
        #print conc
        return conc
 
