import mandos
import sys

if len(sys.argv) != 3:
    print "usage: "+sys.argv[0]+ " <tree> <traits>"
    sys.exit(0)

tree = mandos.tree_utils2.read_tree(sys.argv[1])
traitdic = mandos.tree_utils2.read_continuous(sys.argv[2],tree) 
traits = mandos.tree_utils2.map_continuous(tree,traitdic) #this will return (ntax,ntraits)


root = tree
for node in tree.iternodes():
    if node.label == "N":
        pnode = node
    elif node.label == "f":
        p2node = node.parent
    elif node.label == "a":
        rnode = node.parent

#tree = mandos.fossils.root_tree(tree)[0]
#print mandos.calc_bm_likelihood.bm_prune(tree,traits[1])
#mandos.calc_bm_likelihood.iterate_lengths(tree,traits[1],2)
#print mandos.calc_bm_likelihood.bm_prune(tree,traits[1])

mandos.fossils.place_fossils(tree,traitdic,["g"],traits[1])

#print tree.get_newick_repr(True)

