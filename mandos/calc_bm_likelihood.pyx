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
        if site > 4:
            continue
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
                node.old_length = node.length
                node.charst = charst
        #for tip in prune:
        #    tree_utils2.prune_tip(tree,tip)
            #print tip,len(list(tree.leaves()))
        #for j in tree.iternodes():
        #    j.length = j.old_length

        tlike = 0
        tlike = c_bm_prune(tree,sigsq)
        like += tlike
        for j in tree.iternodes():
            j.length = j.old_length
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
        """
        l_charst = j.children[0].charst
        l_brlen = j.children[0].length
        r_charst = j.children[1].charst
        r_brlen = j.children[1].length
        contrast = l_charst-r_charst 
        cur_var = l_brlen+r_brlen
        """
        #curlike =((-0.5)* ((math.log(2*math.pi*sigsq))+(math.log(cur_var))+(math.pow(contrast,2)/(sigsq*cur_var))))
        l = 0.
        for child in j.children:
            l += ((-0.5)*((math.log(2*math.pi*sigsq))+(math.log(child.length))+(math.pow(child.charst,2)/(sigsq*child.length))))
        print l
        node_likes += curlike
        if j != tree:
            temp_charst = ((r_brlen*l_charst)+(l_brlen*r_charst))/(cur_var)
            temp_brlen = j.length+((l_brlen*r_brlen)/(l_brlen+r_brlen))
            j.charst = temp_charst
            j.length = temp_brlen
    return node_likes

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
