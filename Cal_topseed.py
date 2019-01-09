import numpy as np

class topseed:
    def __init__(self, dGs, NumRanker):
        self.Ranked = dGs
        self.NumRanker = NumRanker

    def rank(self):
        for i in range(np.shape(self.Ranked)[1]):
            self.Ranked = np.array(sorted(self.Ranked, key=lambda a_entry: a_entry[i]))
            self.Ranked[:,i] = 0
            for k in range(self.NumRanker):
                self.Ranked[k,i] = 1
        #print self.Ranked
        return self.Ranked

    def closeness(self):
        self.Closeness = np.zeros((np.shape(self.Ranked)[1],np.shape(self.Ranked)[1]))

        for i in range(np.shape(self.Ranked)[1]):
            for k in range(np.shape(self.Ranked)[1]):
                self.Closeness[i,k] = sum(self.Ranked[:,i]*self.Ranked[:,k])
        #print self.Closeness
        return self.Closeness

    def topseed(self):
        matrix = self.Closeness
        for i in range(np.shape(matrix)[0]):
            matrix[i,i] = 0
        print matrix
        Topseed_value = matrix.max()
        #print Topseed_value
        Topseed_cand = np.array(np.where(matrix == Topseed_value))

        Topseed = Topseed_cand[:,0]
        if Topseed[0] == Topseed[1]:
            Topseed[0] = -1
        print Topseed

        return Topseed
