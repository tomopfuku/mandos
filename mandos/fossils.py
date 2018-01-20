import calc_bm_likelihood
import node2
import numpy as np
import sys

def place_fossils(tree, traits, fossil_names,ntraits):
    for fossil in fossil_names:
        new = node2.Node()
        new.label = fossil
        new.length = 0.01
        new.cont_traits = traits[new.label]
        new.istip = True
        newpar = node2.Node()
        newpar.length = 0.01
        newpar.add_child(new)
        newpar.cont_traits = np.zeros(ntraits,dtype=np.double)
        nodels = [node for node in tree.iternodes() if node != tree]
        curbest = -10000000.
        for node in nodels:
            curlike = 0.
            node.parent.add_child(newpar)
            node.parent.remove_child(node)
            newpar.add_child(node)
            newpar.length = node.length/2
            node.length = newpar.length

            #estimate brlens for placement and calculate the loglike
            calc_bm_likelihood.iterate_lengths(tree,ntraits,2)
            #for n in tree.children:
            rr = tree.root_tree()
            tree = rr[0]
            curlike = calc_bm_likelihood.bm_prune(tree,ntraits)
            print curlike,tree.get_newick_repr(True)
            tree = tree.unroot_tree(rr[1])
            if curlike > curbest:
                curbest = curlike
                curtree = tree.get_newick_repr(True)

            #now restore the original tree and move to the next node
            newpar.parent.add_child(node)
            newpar.children.remove(node)
            newpar.parent.remove_child(newpar)
            node.length = node.length*2
        print curtree

def unroot_tree(tree,oldroot):
    s = [i for i in tree.children if i!=oldroot][0]
    tree.remove_child(oldroot)
    tree.remove_child(s)
    oldroot.add_child(s)
    s.length+=oldroot.length
    oldroot.length = 0.
    return oldroot

def root_tree(tree):
    newroot = node2.Node()
    c1 = tree.children[0]
    tree.remove_child(c1)
    newroot.add_child(c1)
    newroot.add_child(tree)
    c1.length = c1.length/2
    tree.length = c1.length
    return (newroot,tree)
