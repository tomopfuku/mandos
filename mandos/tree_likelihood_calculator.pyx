import math
from model_calc import *
import numpy as np
cimport numpy as np
from numpy import *
from scipy.linalg import *
from scipy.optimize import fmin_powell
from scipy.optimize import fmin_bfgs
import random
import sys

#these functions are based on unpublished code written by Stephen Smith 

#move to datatype
position = {"A":[0],"C":[1],"G":[2],"T":[3],"-":[0,1,2,3],"-":[0,1,2,3],"N":[0,1,2,3],"Y":[1,3],"R":[0,2],"W":[0,3],"M":[0,1],"B":[2,1,3],"V":[2,0,1],"S":[2,1],"K":[2,3],"H":[0,3,1]}

DTYPE = np.double
ctypedef np.double_t DTYPE_t

cdef double LARGE = 10000000000.0

def calc_mk_like(sitels,tree,seqs,optimize=False,random_start = True,numstates=None):
    site_liks = []
    numparams = 0
    for i in sitels:
        cur_sites = sitels[i]
        cur_st = cur_sites[0]-1
        cur_end = cur_sites[1]
        if numstates == None:
            state_space = i
        else:
            state_space = numstates
        match_tips_and_seqs(tree,seqs)
        
        # set rmatrix and base freqs with Mk rate = 1.0
        rmatrix,basefreq = set_multi_jc_rmatrix_basefreq(state_space,1.) 
        if optimize == False:
            curcost = calc_mult_tree_likelihood(tree,rmatrix,basefreq, np.array(range(cur_st,cur_end)),False)
            site_liks.append(curcost)
        else:
            #set arbitrary starting brlens
            
            #optimize brlens
            res = optimize_brlen(tree,range(cur_st,cur_end),rmatrix,basefreq,random_start)
            site_liks.append(res[0][1])
            numparams += res[1]
    return [-sum(site_liks),numparams]

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

cdef calc_mult_tree_likelihood(object tree, np.ndarray[np.float64_t, ndim = 2] rmatrix, np.ndarray[np.float64_t, ndim = 1] basefreq, np.ndarray[np.int_t, ndim = 1] sites, bint verbose=False):
    cdef:
        unsigned int numstates,s,j,i,k
        dict multposition,ps
        double loglike,tempretlike,
        np.ndarray[np.int_t, ndim = 1] positions
        np.ndarray[np.float64_t, ndim = 2] q,p
        object c, child
    numstates = len(basefreq)
    multposition = {}
    
    for i in range(numstates):
        multposition[str(i)] = [i]
    multposition["-"] = range(0,len(basefreq))
    q = calc_mult_q_matrix(rmatrix,basefreq)
    ps = {}#key is length, value is pmatrix
    loglike = 0.0
    for s in sites:
        tempretlike = 0.0
        for c in tree.iternodes(order = 1):
            #if 'seq' in c.data:
                #print c.label,c.data['seq'][s]
            if len(c.children) > 0: #this is an internal node
                c.data['probs'] = ones(len(basefreq))
                if 'seq' in c.data:
                    c.data['probs'] = zeros(len(basefreq))
                    #print multposition
                    positions = np.array(multposition[c.data['seq'][s]])
                    #print positions
                    for j in positions:
                        c.data['probs'][j] = 1
                for j in range(len(basefreq)):
                    for child in c.children:
                        if child.length not in ps:
                            p = calc_p_matrix(q,child.length)
                            #print p
                            ps[child.length] = p
                            #print i.length,p
                        p = ps[child.length]
                        templike = 0
                        for k in range(len(basefreq)):
                            templike += (p[j][k] * child.data['probs'][k])
                        c.data['probs'][j] *= templike
                #print c.get_newick_repr(),c.data['probs']
            else: #this is a tip
                c.data['probs'] = zeros(len(basefreq))
                positions = np.array(multposition[c.data['seq'][s]])
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
    like = -1
    if alpha != None:
        like = calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats) # doesn't work right yet
    else:
        like = calc_mult_tree_likelihood(tree,rmatrix,basefreq,np.array(sites))
    if like < 0 or isnan(like):
        return LARGE
    return like

def optimize_brlen(tree,sites,rmatrix,basefreq,random_start=False,alpha=None,cats=None):
    blstart = []
    for i in tree.iternodes():
        if i != tree:
            if random_start == True:
                i.length = random.uniform(0,0.3)
            blstart.append(i.length)               
    print "optimizing "+str(len(blstart)) +" parameters under Mk"
    res = fmin_bfgs(calc_params_treebl,blstart,args=(tree,sites,rmatrix,basefreq,alpha,cats),full_output=True,disp=False)
    #print tree.get_newick_repr(True)
    return [res,len(blstart)]

def calc_params_treebl_uni(param,br_num,tree,sites,rmatrix,basefreq,alpha=None,cats=None):
    if param[0] < 0:
        return LARGE
    count = 0
    for i in tree.iternodes():
        if i != tree:
            if count == br_num:
                i.length = param[0]
            count += 1
    like = -1
    if alpha != None:
        like = calc_nuc_tree_likelihood_gamma(tree,rmatrix,basefreq,sites,alpha,cats)
    else:
        like = calc_mult_tree_likelihood(tree,rmatrix,basefreq,sites)
    if like < 0 or isnan(like):
        return LARGE
    return like

def optimize_brlen_uni(tree,sites,rmatrix,basefreq,alpha=None,cats=None,random_start=False):
    count = 0
    for i in tree.iternodes():
        if i == tree:
            continue
        if random_start:
            i.length = random.uniform(0.0001,0.3)
        else:
            blstart = [i.length]
        res = fmin_bfgs(calc_params_treebl_uni,blstart,args=(count,tree,sites,rmatrix,basefreq,alpha,cats),full_output=True,disp=False)
        i.length = res[0][0]
        print i.length
        count += 1
    return res

def optimize_rateparams(a,tree,sites,alpha=None,cats=None):
    assert len(a) == 8
    res =  fmin_powell(calc_params_rateparams,a,args=(tree,sites,alpha,cats),full_output=True)
    return res

