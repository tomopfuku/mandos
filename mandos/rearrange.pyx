import tree_reader2,tree_utils2
import stratoML,calc_bm_likelihood
import sys
import math
from copy import deepcopy
from random import sample
import node_opt
import node2
#import cPickle as pickle

def get_random_node(nodels):
    samp = sample(nodels,1)[0]
    if len(samp.children[0].children) != 2 or len(samp.children[1].children) != 2:
        get_random_node(nodels)
        #samp = sample(nodels,1)[0]
    else:
        return samp

cpdef void random_nni(object tree):
    cdef:
        object node,i,samp
        str old_tree = tree.get_newick_repr(True)+";"
        list nni_set,nodels
    nodels = [node for node in tree.iternodes() if node.istip == False and len(node.children[0].children) == 2 and len(node.children[1].children) == 2]
    samp = sample(nodels,1)[0]
    s0 = sample(samp.children[0].children,1)[0]
    s1 = sample(samp.children[1].children,1)[0]
    p0 = s0.parent
    p1 = s1.parent
    p0.remove_child(s0)
    p1.remove_child(s1)
    p1.add_child(s0)
    p0.add_child(s1)
    if p0.length < 0 or p1.length < 0:
        print samp.get_newick_repr(True)
    #return tree

cpdef void random_nni_optim(object tree, dict traits):
    cdef:
        object node,i,samp
        str old_tree = tree.get_newick_repr(True)+";"
        list nni_set,nodels
    nodels = [node for node in tree.iternodes() if node.istip == False and len(node.children[0].children) == 2 and len(node.children[1].children) == 2]
    samp = sample(nodels,1)[0]
    s0 = sample(samp.children[0].children,1)[0]
    s1 = sample(samp.children[1].children,1)[0]
    p0 = s0.parent
    p1 = s1.parent
    p0.remove_child(s0)
    p1.remove_child(s1)
    p1.add_child(s0)
    p0.add_child(s1) 
    """
    for i in samp.iternodes(): # fix negative branch lengths
        if i.length < 0.001:
            i.length = 0.1
    """
    node_opt.bm_single_brlen_optim(s0,p1,traits)
    node_opt.bm_single_brlen_optim(s1,p0,traits)
    #return tree

cpdef void fix_singleton(object node): #remove singleton node resulting from some rearrangements
    if len(node.children) != 1:
        raise ValueError("Is the tree fully bifurcating?")
    
    c = node.children[0]
    p = node.parent
    c.length += node.length
    node.remove_child(c)
    c.parent = p
    p.children.append(c)
    p.length += node.length
    #node.label = "COFFEE"
    #p.label = "COFFEETIP"
    #p.children.remove(node)

cpdef void random_spr(object tree): #TODO: fix issue with children of root
    cdef:
        object regraft,node,i,prune,new,rpar
        double l
        str old_tree = tree.get_newick_repr(True)
        list nodels   
    nodels = [node for node in tree.iternodes() if node.istip == False and node != tree and node.parent != tree and node.parent.parent != None] #TODO see above and below 
    prune = sample(nodels,1)[0]
    print "prune", prune
    par = prune.parent
    gpar = par.parent
    print "gpar",gpar  #TODO why is this coming out as None sometimes?! 
    if gpar == None:
        par.label = "NEW"
        print par.get_newick_repr(True)
        print tree.get_newick_repr(True)
    psib = [node for node in par.children if node != prune][0]
    prune_clade = [node for node in prune.iternodes()]
    nodels = [node for node in tree.iternodes() if node != prune and node.parent != None and node != node.parent and node.parent != gpar and node != psib and node != gpar and node != tree and node.parent != tree and node not in prune_clade]
    regraft = sample(nodels,1)[0]
    print regraft.parent
    if regraft == None:
        regraft.label = "NEW"
        print regraft.get_newick_repr(True)
        print tree.get_newick_repr(True)
    print "rpar",regraft
    psib.length += par.length
    par.remove_child(psib)
    print "gpar",gpar
    
    gpar.add_child(psib)
    gpar.remove_child(par)
    l = regraft.length/2
    regraft.length -= l
    #regraft.parent.label = "REGRAFT"
    rpar = regraft.parent
    if rpar == None or regraft == None:
        regraft.label = "REGRAFT"
        prune.label = "PRUNE"
        par.label = "PAR"
        print regraft.get_newick_repr(True)
        print tree.get_newick_repr(True)

    print "rpar",rpar,"regraft",regraft
    print "prune",prune,"par",par
    rpar.add_child(par)
    par.length = l
    par.add_child(regraft)
    rpar.remove_child(regraft)
    

cpdef void spr(object prune, object regraft):
    cdef:
        object new = node2.Node(), c
    par = prune.parent
    par.remove_child(prune)
    fix_singleton(par)
    c = sample(regraft.children,1)[0]
    c.length -= c.length/2
    new.length = c.length
    regraft.add_child(new)
    regraft.remove_child(c)
    new.add_child(c)
    new.add_child(prune)
   
cpdef void spr_all_reattachments(object tree):
    cdef:
        list nodels
        object regraft,node,i,prune
    nodels = [node for node in tree.iternodes() if node.istip == False and node != tree]
    prune = sample(nodels,1)[0] #pick random node to prune
    start = tree.get_newick_repr(True)+";"
    for node in tree.iternodes():
        spr(prune,node) #reattach subtree at each node

    






