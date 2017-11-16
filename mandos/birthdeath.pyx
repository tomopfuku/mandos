import tree_utils2
import tree_reader2
import math
from scipy import optimize
import numpy as np
cimport numpy as np
import random
import sys
#from libcpp cimport bool

DTYPE = np.double
ctypedef np.double_t DTYPE_t
cdef double LARGE = 10000000000.0

cpdef double single_birth_death_loglike(object tree, double rho):
    cdef:
        #double ,treelik = 0
        list nodels = [node for node in tree.iternodes()]
        unsigned int i, N = len(nodels)
    for i in range(N):
        if nodels[i] == tree:
            continue
        


cpdef double hr97_loglike(object tree, double lam):
    cdef:
        double tf,tl,f,l,a,b,c,top,bot,brlik,treelik =0
        list nodels = [node for node in tree.iternodes()]
        unsigned int i, nfos, N = len(nodels)
    for i in range(N):
        if nodels[i] == tree:
            continue
        tf = nodels[i].parent.height *10000
        tl = nodels[i].height * 10000
        #print tf,tl,i.length,i.length+i.height
        #cdef int nfos 
        nfos = nodels[i].num_occurrences
        if nfos > 1:
            f = nodels[i].upper * 10000
            l = nodels[i].lower*10000
            a = math.log(math.pow(f-l,(nfos-2)))
            b = math.log(math.pow(lam,nfos))
            c = -lam*(tf-tl)
            top = a+b+c
            bot = math.log(math.factorial(nfos-2))
            brlik = top-bot
        elif nfos == 1: #or nfos == 0:
            brlik = math.log(lam)+(-lam*(tf-tl))
            #print math.log(lam*(math.exp(-lam*(tf-tl))))
        elif nfos == 0:
            brlik =  -lam*(tf-tl)
        treelik += brlik
    return treelik

cdef bint assign_node_heights(np.ndarray[DTYPE_t] h, object tree, bint fixed_root=False):
    cdef:
        bint bad = False
        list treels = [node for node in tree.iternodes()]
        unsigned int i,N = len(treels),count = 0
    for i in range(N):
        treels[i].height = h[count]
        count += 1
        if treels[i].parent == None: #set height+length for root
            #i.length = max(i.occurrences)-i.height #set length if root is a "tip"
            continue
        else:
            treels[i].length = treels[i].parent.height-treels[i].height
        #print i.height
        if treels[i].length < 0:
            bad = True
            return bad
        if treels[i].num_occurrences != 0:
            if treels[i].height > treels[i].lower or treels[i].height + treels[i].length < treels[i].upper:
                bad = True
                return bad
    return bad

#calculate likelhood at a single branch. use for wrapping in list comp/loop
def strat_branch_loglike(node,lam):
    tf = node.parent.height *10000
    tl = node.height * 10000
    #print tf,tl,i.length,i.length+i.height
    nfos = node.num_occurrences
    if nfos > 1:
        f = max(node.occurrences) * 10000
        l = min(node.occurrences)*10000
        a = math.log(math.pow(f-l,(nfos-2)))
        b = math.log(math.pow(lam,nfos))
        c = -lam*(tf-tl)
        top = a+b+c
        bot = math.log(math.factorial(nfos-2))
        brlik = top-bot
    elif nfos == 1:
        brlik = math.log(lam)+(-lam*(tf-tl))
    elif nfos == 0:
        brlik =  -lam*(tf-tl)
    return brlik

cpdef double calc_like_strat(np.ndarray[DTYPE_t] p, object tree, dict strat): #p[0] should be lambda, p[1:] is node heights
    cdef:
        unsigned int i,N = len(p)
        double val, LARGE = 1000000.0
    for i in range(N):
        if p[i] < 0:
            return LARGE
    #tree_utils.assign_node_nums(tree)
    cdef bint bad = assign_node_heights(p[1:],tree)

    if bad:
        return LARGE
    try:
        val = -hr97_loglike(tree,p[0])
    except:
        return LARGE
    return val

def optim_lambda_heights(object tree,dict strat,verbose = False):
    #tree_utils.assign_node_nums(tree)
    #print [i.height for i in tree.iternodes()]
    #cdef list start
    cdef list start_val = [0.1]
    cdef list nodes
    nodes = [i.height for i in tree.iternodes(order = 0)]
    x = start_val + nodes
    cdef np.ndarray[DTYPE_t] start
    cdef unsigned int N = len(x)
    start = np.array(x)
    opt = optimize.fmin_powell(calc_like_strat,start,args=(tree,strat),full_output = True, disp = verbose)
    #tree_utils.assign_node_heights(opt[0][1:],tree)
    return [tree,opt]

def calc_like_lambda(p,tree,strat): #p[0] should be lambda, p[1:] is node heights
    for i in p:
        if i < 0:
            return LARGE
    #if p[0] > 1:
    #    return LARGE
    try:
        val = -hr97_loglike(tree,p[0])
    except:
        return LARGE
    return val

def optim_lambda(tree,strat,verbose=False):
    start = [1.]
    opt = optimize.fmin_bfgs(calc_like_lambda,start,args=(tree,strat),full_output = True, disp = verbose)
    return [tree,opt]









