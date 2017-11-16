import sys
import tree_utils2
import tree_reader2
import calc_bm_likelihood,node_opt
import rearrange
import node2

if len(sys.argv) != 3:
    print "usage: "+sys.argv[0]+ " <newick file> <continuous traits>"
    sys.exit(0)


def nni_test(tree,traits):
    node_opt.bm_brlen_optim(tree,traits,"bfgs")
    #node_opt.update_all_brlens(tree,traits)
    #tree.get_newick_repr(True)

    for i in tree.iternodes(True):
        if i.length < 0:
            print i,i.get_newick_repr(True)
    
    cur_best = calc_bm_likelihood.bm_prune(tree,traits)

    for i in range(0,10):
        last_tree = tree.get_newick_repr(True)+";"
        rearrange.random_nni(tree)
        #rearrange.random_nni_optim(tree,traits)
        for node in tree.iternodes():
            if node.length <0:
                print node.length,node.get_newick_repr(True)
        cur_like = calc_bm_likelihood.bm_prune(tree,traits)
        if cur_like < cur_best:
            tree = tree_reader2.read_tree_string(last_tree)
        elif cur_like > cur_best:
            cur_best = cur_like
    return cur_best
   
def spr_test(tree,traits):
    #node_opt.bm_brlen_optim(tree,traits,"bfgs")
    node_opt.update_all_brlens(tree,traits)
    for i in tree.iternodes():
        if i.length < 0.001:
            i.length = 0.1
    print tree.get_newick_repr(True)
    cur_best = calc_bm_likelihood.bm_prune(tree,traits)
    #rearrange.random_spr(tree)
    #cur_best = calc_bm_likelihood.bm_prune(tree,traits)
    #rearrange.random_spr(tree)
    for i in range(0,20):
        print i
        print tree.nnodes("tips")
        last_tree = tree.get_newick_repr(True)+";"
        rearrange.random_spr(tree)
        #print tree.get_newick_repr(True)
        cur_like = calc_bm_likelihood.bm_prune(tree,traits)
        print cur_like
        if cur_like < cur_best:
            print "reverting to previous tree"
            tree = tree_reader2.read_tree_string(last_tree)
        elif cur_like > cur_best:
            cur_best = cur_like
    return cur_best

if __name__ == "__main__":
    tree = tree_utils2.read_tree(sys.argv[1])
    traits = tree_utils2.read_traits(sys.argv[2])
    print tree.get_newick_repr(True)

    #cur_best = nni_test(tree,traits)
    best = spr_test(tree,traits)


    #print node_opt.bm_brlen_optim(tree,traits,"bfgs")
    #print tree.get_newick_repr(True)
    print best

         


#print rearrange.random_spr(tree).get_newick_repr()

"""
t1= tree.get_newick_repr(True)
t2= opt[0].get_newick_repr(True)
print opt[1][1]
print t1
"""
