from tree_utils2 import * 
import tree_reader2
from scipy import optimize
import math
from numpy import random
import sys
import tree_utils2
import numpy as np


cdef double LARGE = 1000000000.0

cpdef double bm_prune(object tree, unsigned int mat_len, double sigsq = 1.0):
    cdef:
        double like,charst,tlike
        unsigned int site
        object node
        list prune = []
        str tip
    like = 0.0
    for site in range(mat_len):
        #if site > 4:
        #    continue
        prune = []
        #old_tree = tree.get_newick_repr(True)+";"
        #for j in tree.iternodes():
        #    j.old_length = j.length        
        for node in tree.iternodes():
            if node.istip:
                charst = node.cont_traits[site]
                #if charst == LARGE:#"?" or val == "-":
                #    prune.append(node.label)
                #    continue
                #node.old_length = node.length
                node.contrast_length = node.length
                node.charst = charst
        #for tip in prune:
        #    tree_utils2.prune_tip(tree,tip)
            #print tip,len(list(tree.leaves()))
        #for j in tree.iternodes():
        #    j.length = j.old_length

        tlike = 0.
        tlike = c_bm_prune(tree,sigsq)
        like += tlike
        #for j in tree.iternodes():
            #j.length = j.old_length
        #print "site:",site,tlike
        #tree = tree_reader2.read_tree_string(old_tree)
    #print like
    return like


cpdef double c_bm_prune(object tree, double sigsq = 1.0):
    cdef:
        double contrast,cur_var,curlike,temp_charst,temp_brlen,l_brlen,r_brlen,l_charst,r_charst
        double node_likes = 0.0
        #list child_charst 
        unsigned int i
        object j
    for j in tree.iternodes(order=1):
        if j.istip: 
            continue
        l_charst = j.children[0].charst
        l_brlen = j.children[0].contrast_length
        r_charst = j.children[1].charst
        r_brlen = j.children[1].contrast_length
        contrast = l_charst-r_charst 
        cur_var = l_brlen+r_brlen
        curlike =((-0.5)* ((math.log(2*math.pi*sigsq))+(math.log(cur_var))+(math.pow(contrast,2)/(sigsq*cur_var))))
        node_likes += curlike
        if j != tree:
            temp_charst = ((r_brlen*l_charst)+(l_brlen*r_charst))/(cur_var)
            temp_brlen = j.contrast_length+((l_brlen*r_brlen)/(l_brlen+r_brlen))
            j.charst = temp_charst
            j.contrast_length = temp_brlen
    return node_likes

"""
cpdef void opt_brlen_all(object tree):
    cdef:
        object node
    for node in tree.iternodes():
        return None
"""
"""
This will be used to prune to the tritomy root of an unrooted tree. 
This calculates the PICs and averaged branch lengths of any internal child branches so that we can estimate ML branch lengths analytically using a three-taxon tree. 
"""

cpdef void prune_to_root(object tree, unsigned int ntraits):
    cdef:
        double contrast,cur_var,temp_charst,temp_brlen,l_brlen,r_brlen,l_charst,r_charst
        #list child_charst 
        unsigned int i
        object j,node
        #double trait

    #for i in range(ntraits): 
    for j in tree.iternodes(order=1):
        if j.istip or j == tree: 
            j.contrast_length = j.length
            continue
        for i in range(ntraits):
            l_charst = j.children[0].cont_traits[i]
            l_brlen = j.children[0].contrast_length
            r_charst = j.children[1].cont_traits[i]
            r_brlen = j.children[1].contrast_length
            
            #calculate PIC for each trait and store in a vector for use at the root
            temp_charst = (((1/r_brlen)*r_charst)+((1/l_brlen)*l_charst))/((1/r_brlen)+(1/l_brlen))
            j.cont_traits[i] = temp_charst
        temp_brlen = 1/((1/r_brlen)+(1/l_brlen))
        temp_brlen += j.length
        j.contrast_length = temp_brlen


"""
this will do the BM prune to any node on an unrooted tree, leaving a tritomy tree
w/ BM adjusted branch lengths subtending
TODO: need to add traits to prune_to_urnode() and prune_to_rnode() 
"""
cpdef void prune_to_urnode(object newroot, object tree, unsigned int ntraits):
    cdef:
        object node,newroot_subtree,c,nrsubtree,c1,c2,sib
        double root_prune_length = 0.0
        unsigned int i
        double tr,bot
        double[:] temp_traits = np.zeros(ntraits,dtype=np.float)

    check_unrooted(len(tree.children))
    
    """
    go through each of the 3 children of the root and mark the subtree containing newroot
    also prune any branches subtending from newroot for later use
    """
    for c in tree.children:
        for node in c.iternodes():
            if node == newroot:
                nrsubtree = c
                c.marked = True
                c1 = node.children[0]
                c2 = node.children[1]
                prune_to_rnode(c1,ntraits)
                prune_to_rnode(c2,ntraits)
                break
    
    #now prune each of the other two subtrees
    for c in tree.children:
        if c.marked == True:
            nrsubtree = c
        elif c.marked == False:
            prune_to_rnode(c,ntraits)
            root_prune_length += 1/c.contrast_length
            for i in range(ntraits):
                temp_charst = ((1/c.contrast_length)*c.cont_traits[i])
                temp_traits[i] += temp_charst

    for i in range(ntraits):
        nrsubtree.cont_traits[i] = temp_traits[i]/root_prune_length #switch directons
    root_prune_length = 1/root_prune_length
    
    nrsubtree.contrast_length = nrsubtree.length
    nrsubtree.contrast_length += root_prune_length

    #now need to prune from tree root to newroot
    if nrsubtree != newroot:
        preorder_prune(nrsubtree,newroot,ntraits)


cpdef void preorder_prune(object nrsubtree, object newroot, unsigned int ntraits):
    cdef:
        object node,c
        double bot,v1,v2,x1,x2,temp_charst
        unsigned int i

    for node in nrsubtree.iternodes(order=0): #we can do this in a preorder traversal
        if node == nrsubtree:
            continue
        for c in node.parent.children:
            if c != node:
                if c.contrast_length == 0.:
                    c.contrast_length = c.length
                sib = c
        v1 = node.parent.contrast_length
        v2 = sib.contrast_length
        node.contrast_length = node.length
        bot = ((1/v1)+(1/v2))
        node.contrast_length += 1/bot
        for i in range(ntraits):
            x1 = node.parent.cont_traits[i]
            x2 = sib.cont_traits[i]
            temp_charst = (((1/v1)*x1)+((1/v2)*x2))/bot
            node.cont_traits[i] = temp_charst
        #print [tr for tr in node.cont_traits]
        #print node.contrast_length
        if node == newroot:
            break

cpdef void prune_to_rnode(object tree, unsigned int ntraits):
    cdef:
        object node
        unsigned int i
        double bot,temp_charst

    for node in tree.iternodes(order=1):
        node.contrast_length = node.length
        if node.istip:
            continue
        c1 = node.children[0]
        c2 = node.children[1]
        bot = ((1/c1.contrast_length)+(1/c2.contrast_length))
        node.contrast_length += 1/bot
        for i in range(ntraits):
            temp_charst = (((1/c1.contrast_length)*c1.cont_traits[i])+((1/c2.contrast_length)*c2.cont_traits[i]))/bot
            node.cont_traits[i] = temp_charst


def check_unrooted(nchildren):
    if nchildren != 3:
        raise AssertionError("Branch lengths need to be estimated on unrooted trees!")

cpdef object pick_new_root(object tree):
    cdef:
        object child

    for child in tree.children:
        if child.istip:
            continue
        else:
            return child


cpdef void iterate_lengths(object tree, unsigned int ntraits, unsigned int iterations=2):
    cdef:
        object node,first_root,curroot,newroot
        unsigned int i,max_nodes,count
        list nodes = []

    nodes = [node for node in tree.iternodes() if not node.istip]
    first_root = tree
    curroot = tree
    count = 0
    max_nodes = 10000
    max_nodes = tree.nnodes()
    for i in range(iterations):
        for count in range(len(nodes)):
            newroot = nodes[count]
            tree = newroot.reroot(tree)
            for node in tree.children:
                prune_to_rnode(node,ntraits)
            tritomy_ML(tree,ntraits)
    tree=first_root.reroot(tree)

cpdef void tritomy_ML(object tree, unsigned int ntraits) except *: 
    cdef:
        double x1,x2,x3,v1,v2,v3,temp_v1,temp_v2,temp_v3
        unsigned int i
   
    check_unrooted(len(tree.children))
    v1 = tree.children[0].length
    v2 = tree.children[1].length
    v3 = tree.children[2].length 
    temp_v1 = 0.0
    temp_v2 = 0.0
    temp_v3 = 0.0
    #ntraits = 1
    for i in range(ntraits): 
        x1 = tree.children[0].cont_traits[i]
        x2 = tree.children[1].cont_traits[i]
        x3 = tree.children[2].cont_traits[i]
        temp_v1 += ((x1-x2)*(x1-x3))
        temp_v2 += ((x2-x1)*(x2-x3))
        temp_v3 += ((x3-x1)*(x3-x2))
    if temp_v1 < 0.0:
        temp_v1 = 0.000001
        temp_v2 = 0.0
        temp_v3 = 0.0
        for i in range(ntraits):
            x1 = tree.children[0].cont_traits[i]
            x2 = tree.children[1].cont_traits[i]
            x3 = tree.children[2].cont_traits[i]
            temp_v2 += math.pow((x1-x2),2)
            temp_v3 += math.pow((x1-x3),2)
    elif temp_v2 < 0.0:
        temp_v1 = 0.0
        temp_v2 = 0.000001
        temp_v3 = 0.0
        for i in range(ntraits):
            x1 = tree.children[0].cont_traits[i]
            x2 = tree.children[1].cont_traits[i]
            x3 = tree.children[2].cont_traits[i]
            temp_v1 += math.pow((x2-x1),2)
            temp_v3 += math.pow((x2-x3),2)
    elif temp_v3 < 0.0:
        temp_v1 = 0.0
        temp_v2 = 0.0
        temp_v3 = 0.00001
        for i in range(ntraits):
            x1 = tree.children[0].cont_traits[i]
            x2 = tree.children[1].cont_traits[i]
            x3 = tree.children[2].cont_traits[i]
            temp_v1 += math.pow((x3-x1),2)
            temp_v2 += math.pow((x2-x3),2)
    v1 = temp_v1/ntraits
    v2 = temp_v2/ntraits
    v3 = temp_v3/ntraits
    tree.children[0].length = v1
    tree.children[1].length = v2
    tree.children[2].length = v3
    #print tree.get_newick_repr(True)
    #print v1,v2,v3
    v1 = v1-(tree.children[0].contrast_length-tree.children[0].length)
    v2 = v2-(tree.children[1].contrast_length-tree.children[1].length)
    v3 = v3-(tree.children[2].contrast_length-tree.children[2].length)
    if v1 <= 0:
        v1 = 0.0001
    elif v2 <= 0:
        v2 = 0.0001
    if v3 <= 0:
        v3 = 0.0001



cpdef void node_ML(object tree, unsigned int ntraits):
    cdef:
        double x1,x2,x3,v1,v2,v3,temp_v1,temp_v2,temp_v3
        unsigned int i
   
    v1 = tree.contrast_length
    v2 = tree.children[0].contrast_length
    v3 = tree.children[1].contrast_length 
    temp_v1 = 0.0
    temp_v2 = 0.0
    temp_v3 = 0.0
    #ntraits = 1
    for i in range(ntraits): 
        x1 = tree.cont_traits[i]
        x2 = tree.children[0].cont_traits[i]
        x3 = tree.children[1].cont_traits[i]
        temp_v1 += ((x1-x2)*(x1-x3))
        temp_v2 += ((x2-x1)*(x2-x3))
        temp_v3 += ((x3-x1)*(x3-x2))
    
    if temp_v1 < 0.0:
        temp_v1 = 0.0001
        temp_v2 = 0.0
        temp_v3 = 0.0
        for i in range(ntraits):
            x1 = tree.cont_traits[i]
            x2 = tree.children[0].cont_traits[i]
            x3 = tree.children[1].cont_traits[i]
            temp_v2 += math.pow((x1-x2),2)
            temp_v3 += math.pow((x1-x3),2)
    elif temp_v2 < 0.0:
        temp_v1 = 0.0
        temp_v2 = 0.0001
        temp_v3 = 0.0
        for i in range(ntraits):
            x1 = tree.cont_traits[i]
            x2 = tree.children[0].cont_traits[i]
            x3 = tree.children[1].cont_traits[i]
            temp_v1 += math.pow((x2-x1),2)
            temp_v3 += math.pow((x2-x3),2)
    elif temp_v3 < 0.0:
        temp_v1 = 0.0
        temp_v2 = 0.0
        temp_v3 = 0.0001
        for i in range(ntraits):
            x1 = tree.cont_traits[i]
            x2 = tree.children[0].cont_traits[i]
            x3 = tree.children[1].cont_traits[i]
            temp_v1 += math.pow((x3-x1),2)
            temp_v2 += math.pow((x2-x3),2)
    v1 = temp_v1/ntraits
    v2 = temp_v2/ntraits
    v3 = temp_v3/ntraits

    v1 = v1-(tree.contrast_length-tree.length)
    v2 = v2-(tree.children[0].contrast_length-tree.children[0].length)
    v3 = v3-(tree.children[0].contrast_length-tree.children[0].length)
    if v1 <= 0:
        v1 = 0.0001
    elif v2 <= 0:
        v2 = 0.0001
    if v3 <= 0:
        v3 = 0.0001
    tree.length = v1
    tree.children[0].length = v2
    tree.children[1].length = v3
    #print tree.parent.get_newick_repr(True)
    #print v1,v2,v3

cpdef double sigsqML(object tree): #tree must already have characters mapped to tips using match_traits_tips()
    cdef:
        unsigned int n,p
        list vals,x,t
        double ui,Vi,add,div,sig2
    n = len(tree.leaves())
    vals = [None]*(n-2)
    p = 0
    for i in tree.iternodes(order=1):
        if i.istip == False and i != tree:
            x = [j.charst for j in i.children]
            t = [j.length for j in i.children]
            ui = abs(x[0]-x[1])
            Vi = sum(t)
            vals[p] = (ui,Vi)
            add = (t[0]*t[1])/(t[0]+t[1])
            i.length = i.length + add
            p += 1
        if i == tree:
            t = [j.length for j in i.children]
            Vi = sum(t)
            V0 = (t[0]*t[1])/(t[0]+t[1])
    div = sum([math.pow(i[0],2)/i[1] for i in vals])+(0/V0)
    sig2 = (1./n) * div
    return sig2
