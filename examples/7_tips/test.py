import mandos
import sys

if len(sys.argv) != 3:
    print "usage: "+sys.argv[0]+ " <tree> <traits>"
    sys.exit(0)

tree = mandos.tree_utils2.read_tree(sys.argv[1])
traits = mandos.tree_utils2.read_continuous(sys.argv[2],tree) #this will return (ntax,ntraits)

root = tree
for node in tree.iternodes():
    if node.label == "N":
        pnode = node
    elif node.label == "f":
        p2node = node.parent
    elif node.label == "a":
        rnode = node.parent

print tree.get_newick_repr(True)

mandos.calc_bm_likelihood.iterate_lengths(tree,traits[1],1)
"""
for i in range(5):
    for node in tree.children:
        mandos.calc_bm_likelihood.prune_to_rnode(node,traits[1])
    mandos.calc_bm_likelihood.tritomy_ML(tree,traits[1])

    tree = pnode.parent.reroot(tree)
    for node in tree.children:
        mandos.calc_bm_likelihood.prune_to_rnode(node,traits[1])
    mandos.calc_bm_likelihood.tritomy_ML(tree,traits[1])

    tree = pnode.reroot(tree)
    for node in tree.children:
        mandos.calc_bm_likelihood.prune_to_rnode(node,traits[1])
    mandos.calc_bm_likelihood.tritomy_ML(tree,traits[1])

    tree = p2node.reroot(tree)
    for node in tree.children:
        mandos.calc_bm_likelihood.prune_to_rnode(node,traits[1])
    mandos.calc_bm_likelihood.tritomy_ML(tree,traits[1])

    tree = rnode.reroot(tree)
    for node in tree.children:
        mandos.calc_bm_likelihood.prune_to_rnode(node,traits[1])
    mandos.calc_bm_likelihood.tritomy_ML(tree,traits[1])
    tree = root.reroot(tree)
"""

"""
for node in tree.iternodes(order=1):
    if node == tree:
        for node in tree.children:
            mandos.calc_bm_likelihood.prune_to_rnode(node,traits[1])
        mandos.calc_bm_likelihood.tritomy_ML(tree,traits[1])
        continue
    elif node.istip:
        continue
    mandos.calc_bm_likelihood.prune_to_urnode(node,tree,traits[1])
    mandos.calc_bm_likelihood.node_ML(node,traits[1])
"""

print tree.get_newick_repr(True)

