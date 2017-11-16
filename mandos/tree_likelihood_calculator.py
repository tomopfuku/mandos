import math
from model_calc import *
import numpy as np
from numpy import *
from scipy.linalg import *
from scipy.optimize import fmin_powell
from scipy.optimize import fmin_bfgs
import random
import sys

#these functions contributed by Stephen Smith 

#move to datatype
position = {"A":[0],"C":[1],"G":[2],"T":[3],"-":[0,1,2,3],"-":[0,1,2,3],"N":[0,1,2,3],"Y":[1,3],"R":[0,2],"W":[0,3],"M":[0,1],"B":[2,1,3],"V":[2,0,1],"S":[2,1],"K":[2,3],"H":[0,3,1]}


LARGE = 10000000

def calc_mk_like(sitels,tree,seqs,optimize=False):
    site_liks = []
    for i in range(len(sitels)):
        cur_sites = sitels[i]
        cur_st = cur_sites[0]-1
        cur_end = cur_sites[1]
        if i+1 == 1:
            state_space = 2
        else:
            state_space = i+1
        match_tips_and_seqs(tree,seqs)
        
        # set rmatrix and base freqs with Mk rate = 1.0
        rmatrix,basefreq = set_multi_jc_rmatrix_basefreq(state_space,1.) 
        if optimize == False:
            curcost = calc_mult_tree_likelihood(tree,rmatrix,basefreq,range(cur_st,cur_end),False)
            site_liks.append(curcost)
        else:
            #set arbitrary starting brlens
            for n in tree.iternodes():
                n.length = 0.01
            
            #optimize brlens
            res = optimize_brlen(tree,range(cur_st,cur_end),rmatrix,basefreq)
            site_liks.append(res[1])
    return -sum(site_liks)

def match_tips_and_seqs(tree,seqs):
    lvs = [i for i in tree.iternodes() if i.label != ""]
    for i in lvs:
        test = False
        for j in seqs:
            if j.name == i.label:
                i.data['seq'] = j.seq#.upper()
                test = True
                break
        if test == False:
            print "can't find "+i.label+" in seqs"
            return False

def calc_mult_tree_likelihood(tree,rmatrix,basefreq,sites,verbose=False):
    numstates = len(basefreq)
    multposition = {}
    for i in range(numstates):
        multposition[str(i)] = [i]
    multposition["-"] = range(0,len(basefreq))
    q = calc_mult_q_matrix(rmatrix,basefreq)
    ps = {}#key is length, value is pmatrix
    loglike = 0
    for s in sites:
        tempretlike = 0
        for c in tree.iternodes(order = 1):
            #if 'seq' in c.data:
                #print c.label,c.data['seq'][s]
            if len(c.children) > 0: #this is an internal node
                c.data['probs'] = ones(len(basefreq))
                if 'seq' in c.data:
                    c.data['probs'] = zeros(len(basefreq))
                    #print multposition
                    positions = multposition[c.data['seq'][s]]
                    for j in positions:
                        c.data['probs'][j] = 1
                for j in range(len(basefreq)):
                    for i in c.children:
                        if i.length not in ps:
                            p = calc_p_matrix(q,i.length)
                            #print p
                            ps[i.length] = p
                            #print i.length,p
                        p = ps[i.length]
                        templike = 0
                        for k in range(len(basefreq)):
                            templike += (p[j][k] * i.data['probs'][k])
                        c.data['probs'][j] *= templike
                #print c.get_newick_repr(),c.data['probs']
            else: #this is a tip
                c.data['probs'] = zeros(len(basefreq))
                positions = multposition[c.data['seq'][s]]
                for j in positions:
                    c.data['probs'][j] = 1
        #the tempretlike will be the sum(root prob * basefreq) for all bases
        tempretlike = 0
        for i in range(len(basefreq)):
            tempretlike  += tree.data['probs'][i]*basefreq[i]
        if verbose:
            print "site",s,"log(L):",log(tempretlike),"like:",tempretlike
        loglike -= log(tempretlike)
    return loglike



def set_multi_jc_rmatrix_basefreq(numstates,rate):
    #rmatrix = ones((numstates,numstates))
    rmatrix = full((numstates,numstates),rate)
    fill_diagonal(rmatrix,0)
    basecomp = ones(numstates)*(1/float(numstates))  
    sc = sum(basecomp)
    basecomp /= sc
    return rmatrix,basecomp



"""
params contain the rmatrix and basefreq params in the form
  A C G T
A   1 2 3
C     4 5
G 
T 
the three base composition rates will be parameters
6, 7, 8. The last base composition (T) will be 1- sum(6,7,8)
and the total will be scaled to 1
def calc_params_rateparams(params,tree,sites,alpha=None,cats=None):
    for i in params: 
        if i < 0:
            return LARGE
    if sum(params[5:8]) > 1:
        return LARGE
    rmatrix,basefreq = set_gtr_rmatrix_basefreq(params)
    like = -1
    if alpha != None:
        like = calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats)
    else:
        like = calc_nuc_tree_likelihood(tree,rmatrix,basefreq,sites)
    #print like
    if like < 0 or isnan(like):
        return LARGE
    return like
"""
def calc_jcparams_rateparams(rate,tree,sites,states):
    if rate[0] < 0:
        return LARGE
    rmatrix,basefreq = set_multi_jc_rmatrix_basefreq(states,rate)
    rmatrix *= rate[0]
    like = calc_mult_tree_likelihood(tree,rmatrix,basefreq,sites)
    if like<0 or isnan(like):
        return LARGE
    return like

def optimize_jcrateparams(a,tree,sites,states):
    assert len(a) == 1
    res = fmin_bfgs(calc_jcparams_rateparams,a,args = (tree,sites,states),disp=False,full_output=True)
    return res

def calc_params_treebl(params,tree,sites,rmatrix,basefreq,alpha=None,cats=None):
    count = 0
    for i in tree.iternodes():
        if i != tree:
            if params[count] < 0:
                return LARGE
            i.length = params[count]
            count += 1
        #else:
            #print i.length
    like = -1
    if alpha != None:
        like = calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats)
    else:
        like = calc_mult_tree_likelihood(tree,rmatrix,basefreq,sites)
    #print tree.get_newick_repr(True)
    #print like
    if like < 0 or isnan(like):
        return LARGE
    return like

def optimize_brlen(tree,sites,rmatrix,basefreq,alpha=None,cats=None):
    blstart = []
    for i in tree.iternodes():
        if i != tree:
            blstart.append(i.length)
            blstart.append(random.uniform(0,0.3))

    res = fmin_bfgs(calc_params_treebl,blstart,args=(tree,sites,rmatrix,basefreq,alpha,cats),full_output=True,disp=False)
    #print tree.get_newick_repr(True)
    return res

def optimize_rateparams(a,tree,sites,alpha=None,cats=None):
    assert len(a) == 8
    res =  fmin_powell(calc_params_rateparams,a,args=(tree,sites,alpha,cats),full_output=True)
    return res

