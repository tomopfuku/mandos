import mandos
import sys

if len(sys.argv) != 3:
    print "usage: "+sys.argv[0]+ " <tree> <traits>"
    sys.exit(0)

tree = mandos.tree_utils2.read_tree(sys.argv[1])
traits = mandos.tree_utils2.read_continuous(sys.argv[2],tree) #this will return (ntax,ntraits)

print tree.get_newick_repr(True)

print mandos.calc_bm_likelihood.bm_prune(tree,traits[1])
mandos.calc_bm_likelihood.tritomy_ML(tree,traits[1])
print mandos.calc_bm_likelihood.bm_prune(tree,traits[1])
print tree.get_newick_repr(True)

