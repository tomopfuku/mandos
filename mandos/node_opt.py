import math
from scipy import optimize
import tree_reader2
import tree_utils2
import calc_bm_likelihood
import random
import sys
LARGE = 100000000000


def bm_sigsq_optim(tree,traits,rate=1):
    tree_utils2.assign_sigsq(tree,[rate])
    start = [rate]
    print start
    opt = optimize.fmin_bfgs(calc_like_sigsq, start, args = (tree,traits),full_output=False,disp=True)
    return [tree.get_newick_repr(True),opt]

def calc_like_sigsq(sigsq,tree,traits):
    if sigsq < 0:
        return LARGE
    tree_utils2.assign_sigsq(tree,sigsq)
    try:
        val = -calc_bm_likelihood.bm_prune(tree,traits)
    except:
        return LARGE
    #print ht,val
    return val

def bm_brlen_optim(tree,traits,method="bfgs"):
    tree_utils2.assign_branch_nums(tree)
    #tree_utils2.assign_sigsq(tree,[1.0])
    start = [i.length for i in tree.iternodes() if i != tree]
    #print len(start);sys.exit(0)
    if method.lower() == "bfgs" or method.lower() == "unconstrained":
        opt = optimize.fmin_bfgs(calc_like_brlens, start, args = (tree,traits),full_output=False,disp=True)
    elif method.lower() == "l-bfgs-b":
        bounds = [(0.0,100000.0)]*len(start) 
        opt = optimize.fmin_l_bfgs_b(calc_like_brlens, start, approx_grad = True,bounds = bounds, args = (tree,traits))
    return [tree.get_newick_repr(True),opt]

def fix_bad_brlens(node):
    for i in node.iternodes():
        if i.length < 0.001:
            i.length = 0.1

def calc_like_brlens(l,tree,traits):
    for i in l:
        if i < 0:
            return LARGE
    bad = tree_utils2.assign_brlens(l,tree)
    if bad:
        return LARGE
    try:
        ll = -calc_bm_likelihood.bm_prune(tree,traits)
    except:
        return LARGE
    return ll

def calc_like_single_brlen(bl,branch,tree,traits):
    if bl[0] < 0:
        return LARGE
    branch.length = bl
    try:
        ll = -calc_bm_likelihood.bm_prune(tree,traits)
    except:
        return LARGE
    return ll

def bm_single_brlen_optim(branch,tree,traits,verbose = False):
    start = [branch.length]
    if verbose == False:
        opt = optimize.fmin_bfgs(calc_like_single_brlen, start, args = (branch,tree,traits),full_output=False,disp=False)
    else:
        opt = optimize.fmin_bfgs(calc_like_single_brlen, start, args = (branch,tree,traits),full_output=False,disp=True)
    return [tree.get_newick_repr(True),opt]

def update_all_brlens(tree,traits):
    for i in tree.iternodes():
        if i == tree:
            continue
        bm_single_brlen_optim(i,tree,traits)
        

