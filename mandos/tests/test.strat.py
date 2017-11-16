import sys
import tree_utils2
import stratoML

if len(sys.argv) != 3:
    print "usage: "+sys.argv[0]+ " <newick file> <range data file>"
    sys.exit(0)

tree = tree_utils2.read_tree(sys.argv[1])
ranges = tree_utils2.read_strat(sys.argv[2])
tree_utils2.match_strat(tree,ranges)
ranges = None

tree_utils2.init_heights_strat(tree)
#print tree.get_newick_repr(True)
#stratoML.collapse_non_overlapping_ranges(tree)

#print -stratoML.hr97_loglike(tree,lam=1.0)
#print tree.get_newick_repr(True)

opt = stratoML.optim_lambda_heights(tree,ranges)

#print opt
t1= tree.get_newick_repr(True)
t2= opt[0].get_newick_repr(True)
print opt[1][1]
print t1
