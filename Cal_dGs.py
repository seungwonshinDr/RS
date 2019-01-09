import numpy as np
from itertools import product
class dGs:
    def __init__(self, strands, Temp):
        self.strands = strands
        self.Temp = Temp
        self.commons = np.array(list(product([1, 2, 3, 4], repeat= np.shape(self.strands)[1])))

    def popcode(self):
        return self.commons

    def dGs(self):
        hyb_code = np.zeros((np.shape(self.strands)[0], np.shape(self.commons)[0], np.shape(self.strands)[1] - 1))
        fitness = np.zeros((np.shape(self.commons)[0], np.shape(self.strands)[0]))
        #np.savetxt('popcodes_'+str(np.shape(self.strands)[1]), self.commons, delimiter= ',')

        for k in range(np.shape(self.strands)[0]):
            for i in range(np.shape(self.commons)[0]):
                for j in range(np.shape(self.commons)[1] - 1):
                    hyb_code[k, i, j] = self.strands[k, j] * 1000 + self.strands[k, j + 1] * 100 + self.commons[i, j] * 10 + self.commons[i, j + 1]

        for k in range(np.shape(self.strands)[0]):
            for i in range(np.shape(self.commons)[0]):
                for j in range(np.shape(self.commons)[1] - 1):
                    if hyb_code[k, i, j] == 1122:
                        fitness[i, k] += -7.9 - self.Temp * (-22.2)
                    elif hyb_code[k, i, j] == 1221:
                        fitness[i, k] += -7.2 - self.Temp * (-20.4)
                    elif hyb_code[k, i, j] == 1324:
                        fitness[i, k] += -7.8 - self.Temp * (-21)
                    elif hyb_code[k, i, j] == 1423:
                        fitness[i, k] += -8.4 - self.Temp * (-22.4)
                    elif hyb_code[k, i, j] == 2112:
                        fitness[i, k] += -7.2 - self.Temp * (-21.3)
                    elif hyb_code[k, i, j] == 2211:
                        fitness[i, k] += -7.9 - self.Temp * (-22.2)
                    elif hyb_code[k, i, j] == 2314:
                        fitness[i, k] += -8.5 - self.Temp * (-22.7)
                    elif hyb_code[k, i, j] == 2413:
                        fitness[i, k] += -8.2 - self.Temp * (-22.2)
                    elif hyb_code[k, i, j] == 3142:
                        fitness[i, k] += -8.2 - self.Temp * (-22.2)

                    elif hyb_code[k, i, j] == 3241:
                        fitness[i, k] += -8.4 - self.Temp * (-22.4)
                    elif hyb_code[k, i, j] == 3344:
                        fitness[i, k] += -8 - self.Temp * (-19.9)
                    elif hyb_code[k, i, j] == 3443:
                        fitness[i, k] += -9.8 - self.Temp * (-24.4)
                    elif hyb_code[k, i, j] == 4132:
                        fitness[i, k] += -8.5 - self.Temp * (-22.7)
                    elif hyb_code[k, i, j] == 4231:
                        fitness[i, k] += -7.8 - self.Temp * (-21)
                    elif hyb_code[k, i, j] == 4334:
                        fitness[i, k] += -10.6 - self.Temp * (-27.2)
                    elif hyb_code[k, i, j] == 4433:
                        fitness[i, k] += -8 - self.Temp * (-19.9)
                    else:
                        fitness[i, k] += 0
        return fitness
