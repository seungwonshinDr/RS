import numpy as np
import sympy
from sympy import Symbol, solve, re

class ratio:
    def __init__(self, dGs, concentration):
        self.dGs = dGs
        self.conc = concentration
        self.product = np.zeros((np.shape(self.dGs)[0],4))

    def cal_ratio(self):
        for i in range(np.shape(self.dGs)[0]):
            print i+1, np.shape(self.dGs)[0]
            dGa = self.dGs[i, 1]
            dGb = self.dGs[i, 2]
            a = self.conc[0,0]
            b = self.conc[0,1]
            c = a + b
            print dGa, dGb
            #print con_a, con_b, con_anti

            Ka = sympy.exp(dGa * 4184 / (-8.314) / 300)
            Kb = sympy.exp(dGb * 4184 / (-8.314) / 300)

            x = Symbol('x')
            y = Symbol('y')
            ans = np.array(solve([(x*(a+b+c-x-y)-Ka*(a-x)*(c-x-y)), (y*(a+b+c-x-y)-Kb*(b-y)*(c-x-y))], x, y))

            if np.shape(np.shape(ans))[0] > 1:
                re_ans = np.zeros((np.shape(ans)[0], np.shape(ans)[1]))
                for j in range(np.shape(ans)[0]):
                    for k in range(np.shape(ans)[1]):
                        re_ans[j, k] = re(ans[j, k])

                #print np.where(re_ans<0)

                re_ans = np.delete(re_ans, np.where(re_ans<0)[0],0)

                #print re_ans
                #print self.conc
                array = np.array(np.where(self.conc == self.conc.min())).max()
                re_point = np.array(np.where(re_ans[:, array] < self.conc.min()))
                #re_point = np.array(np.where(re_ans < self.conc.min()))
                sel_ans = np.array(re_ans[re_point[0,:].max(), :])
                #print re_ans
                #print sel_ans
                #print re_point[:,0].max()

            elif np.shape(ans)[0] == 2:
                re_ans = np.zeros((1, 2))
                for j in range(1):
                    re_ans[j] = re(ans[j])
                sel_ans = re_ans
            else:
                sel_ans = np.array([0, 0])

#            if sel_ans.min() < 0:
#                sel_ans = np.array([0, 0])

            #print re_ans
            print sel_ans

            self.product[i, 0] = self.dGs[i,0]
            self.product[i, 1] = sel_ans[0]
            self.product[i, 2] = sel_ans[1]
            self.product[i, 3] = sel_ans[0] + sel_ans[1]
            print self.product
        return self.product

