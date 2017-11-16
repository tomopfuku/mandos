from tree_utils2 import * 
import tree_reader2
from scipy import optimize
import math
from numpy import random
import sys

cdef double LARGE = 1000000000.0

cpdef double bm_prune(object tree, dict traits):
    cdef:
        double like
    like = c_bm_prune(tree,traits)
    return like

cdef double c_bm_prune(object tree, dict traits, double sigsq = 1.0):
    cdef:
        double contrast,cur_var,curlike,temp_charst,temp_brlen,trait_likes = 0.0
        double node_likes = 0.0
        list child_charst 
        unsigned int i,mat_len
        object j
    mat_len = len(traits.values()[0])
    for i in range(mat_len):
        #match_traits_tips(tree,traits,i)
        for j in tree.iternodes():
            if j.istip:
                j.charst = traits[j.label][i]
            j.old_length = j.length
        for j in tree.iternodes(order=1):
            if j.istip == False and j != tree:
                child_charst = [k.charst for k in j.children]
                brlens = [k.length for k in j.children]
                contrast = child_charst[0]-child_charst[1]
                cur_var = brlens[0]+brlens[1]
                curlike =((-0.5)* ((math.log(2*math.pi*sigsq))+(math.log(cur_var))+(math.pow(contrast,2)/(sigsq*cur_var))))
                node_likes += curlike
                #temp_charst = (((1/brlens[0])*child_charst[0])+((1/brlens[1])*child_charst[1]))/((1/brlens[0])+(1/brlens[1]))
                temp_charst = ((brlens[1]*child_charst[0])+(brlens[0]*child_charst[1]))/(sum(brlens))
                temp_brlen = j.length+((brlens[0]*brlens[1])/(brlens[0]+brlens[1]))
                j.charst = temp_charst
                j.length = temp_brlen
            elif j == tree:
                child_charst = [k.charst for k in j.children]
                brlens = [k.length for k in j.children]
                contrast = child_charst[0]-child_charst[1]
                cur_var = brlens[0]+brlens[1]
                curlike =((-0.5)* ((math.log(2*math.pi*sigsq))+(math.log(cur_var))+(math.pow(contrast,2)/(sigsq*cur_var))))
                node_likes += curlike
        for j in tree.iternodes():
            j.length = j.old_length
        trait_likes += node_likes
    #print -sum(trait_likes)
    return trait_likes

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
