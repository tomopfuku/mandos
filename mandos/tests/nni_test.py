from mandos import *
from copy import deepcopy

#nwk = open("/home/carolioness/documents/coding_stuff/mandos/examples/cetaceans/cetacean.part.phy.tre","r").readline().strip()
nwk = open("/home/carolioness/documents/current_projects/hominoid_stratoclad/SC_ml/hominin.MORPH.tre","r").readline().strip()
tree = tree_reader.read_tree_string(nwk)

#print tree_utils.nni_set(tree,tree)
newtree = deepcopy(tree)
"""
for node in newtree.iternodes():
    if node.istip == False and len(node.descendants()) > 2:
        #print len(tree_utils.nni_set(node,tree))#.nni_set()) #should be either 3 or 4
        print tree_utils.nni_set(node,newtree)
        print "\n\n"
"""
likefl = "/home/carolioness/documents/current_projects/hominoid_stratoclad/branch_support/TREELIKES.1764.rr.likelihoods.nni.trees"
stratfl = "/home/carolioness/documents/current_projects/hominoid_stratoclad/SC_ml/hom_occurrences.csv"

print branch_support.aBayes_precalc_output(tree,likefl,stratfl)
