import mandos
import sys

if len(sys.argv) != 3:
    print "usage: "+sys.argv[0]+" <tree> <traits>"
    sys.exit(0)

tree = mandos.tree_utils2.read_tree(sys.argv[1])
traits = mandos.tree_utils2.read_continuous(sys.argv[2],tree) #this will return (ntax,ntraits)

mandos.calc_bm_likelihood.tritomy_ML(tree,traits[1])



