import numpy as np
import sympy
from sympy import Symbol
from itertools import product
import datetime
import time
import scipy
from scipy import stats

import Cal_dGs as dGs
import Cal_topseed as topseed
import Cal_ratio as ratio
import Gen_representer as gen
#import Pearson as pear


filename = 'Analogs'
cluster = np.loadtxt(filename)

Temp = 0.3

NumRanker_topseed = 1000
NumRanker_ratio = 100
concentrations = np.ones((1,np.shape(cluster)[0]))

#topseed_list = np.zeros((1,2))

print cluster

for h in range(np.shape(cluster)[0]-1):
    cluster_init = dGs.dGs(cluster, Temp)
    cluster_dGs = cluster_init.dGs()
    pearson_list = np.zeros((np.shape(cluster_dGs)[1],np.shape(cluster_dGs)[1]))

    for i in range(np.shape(cluster_dGs)[1]):
        for j in range(np.shape(cluster_dGs)[1]):
            pearson_list[i,j] = np.array(scipy.stats.pearsonr(cluster_dGs[:, i], cluster_dGs[:, j]))[0]
    np.savetxt(filename+'_Pearson_K'+str(h), pearson_list, delimiter= ',')
    now = datetime.datetime.now()
    print now

    np.savetxt(filename+'_dGs_K'+str(h), cluster_dGs, delimiter=',')
    #np.savetxt(filename + '_cluster_K' + str(h), cluster, delimiter=',')

    topseed_init = topseed.topseed(cluster_dGs, NumRanker_topseed)
    topseed_rank = topseed_init.rank()
    topseed_close = topseed_init.closeness()
    topseed_topseed = topseed_init.topseed()
    if topseed_topseed.min() < 0:
        print 'More NumRanker needs!!'
        break

    #print np.shape(topseed_list), np.shape(topseed_topseed)
    #topseed_list = np.concatenate((topseed_list, topseed_topseed), axis = 0)


    np.savetxt(filename+'_rankers_K'+str(h), topseed_rank, delimiter=',')
    np.savetxt(filename + '_closeness_K' + str(h), topseed_close, delimiter=',')


    conc = np.zeros((1,2))
    conc[0,0] = concentrations[0,topseed_topseed[0]]
    conc[0,1] = concentrations[0,topseed_topseed[1]]



    strand_dGs = np.zeros((np.shape(cluster_dGs)[0],3))
    for i in range(np.shape(cluster_dGs)[0]):
        strand_dGs[i,0] = i
        strand_dGs[i,1] = cluster_dGs[i, topseed_topseed[0]]
        strand_dGs[i,2] = cluster_dGs[i, topseed_topseed[1]]
        #print cluster_dGs[i, topseed_topseed[0]], cluster_dGs[i, topseed_topseed[1]]
    #np.savetxt('test', strand_dGs, delimiter=',')

    selected_strands = np.zeros((NumRanker_ratio,3))
    strand_dGs_order = np.array(sorted(strand_dGs, key= lambda a_entry: a_entry[2]))

    for i in range(NumRanker_ratio):
        selected_strands[i, 0] = strand_dGs_order[i, 0]
        selected_strands[i, 1] = strand_dGs_order[i, 1]
        selected_strands[i, 2] = strand_dGs_order[i, 2]


    ratio_init = ratio.ratio(selected_strands, conc)
    ratio_cal = ratio_init.cal_ratio()
    ratio_cal_order = np.array(sorted(ratio_cal, key= lambda a_entry: a_entry[3], reverse= True))
    #print ratio_cal_order
    common_seq_code = ratio_cal_order[0,0]

    popcodes = cluster_init.popcode()

    sub_strand_1 = cluster[topseed_topseed[0],:]
    sub_strand_2 = cluster[topseed_topseed[1],:]
    print 'sub_strand1', sub_strand_1, conc[0,0]
    print 'sub_strand2',sub_strand_2, conc[0,1]

    common_seq = np.zeros((1,np.shape(popcodes)[1]))
    for i in range(np.shape(popcodes)[1]):
        common_seq[0,i] = popcodes[common_seq_code,i]


    rep_init = gen.gen(common_seq, ratio_cal_order[0,3], Temp)
    rep_gen = rep_init.rep()
    rep_conc = np.zeros((1,1))
    rep_conc[0,0] = rep_init.concentration()
    print 'representor:', rep_gen, rep_conc[0,0]

    print 'topseed:', topseed_topseed
    concentrations = np.delete(concentrations, topseed_topseed.min(),1)
    concentrations = np.delete(concentrations, topseed_topseed.max()-1, 1)
    concentrations = np.concatenate((concentrations, rep_conc), axis=1)
    print 'concentrations:', concentrations

    #print cluster
    cluster = np.delete(cluster, topseed_topseed.min(),0)
    cluster = np.delete(cluster, topseed_topseed.max()-1, 0)
    cluster = np.concatenate((cluster, rep_gen), axis = 0)
    #print cluster

    np.savetxt(filename + '_cluster_K' + str(h), cluster, delimiter=',')
    np.savetxt(filename + '_representer_K' + str(h), rep_gen, delimiter=',')
    np.savetxt(filename + '_conc_K' + str(h), concentrations, delimiter=',')

