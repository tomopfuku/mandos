import sys,os
from mandos import *

"""
if len(sys.argv) != 7:
    print "usage: " +sys.argv[0]+ "<file containing trees> <stratigraphic range file> <morphological alignment> <morphological partitions file> <right and left descendents of root node>"
    sys.exit(0)
"""
def search_trees(trees,strt,aln,parts,root_tax,likefl):
    """
    trees = sys.argv[1]
    strt = sys.argv[2]
    aln = sys.argv[3]
    parts = sys.argv[4]
    root_tax = sys.argv[5]
    likefl = sys.argv[6]
    """
    aln_liks = []
    treels = []

    for i in open(likefl,"r"):
        spls = i.strip().split("\t")
        nwk = spls[1].strip()
        lik = float(spls[0])
        aln_liks.append(lik)
        treels.append(nwk)

    best = None
    best_strat = None
    bls = []
    interval = 100
    prog = 0
    for i in range(len(aln_liks)):#aln_liks.keys():
        tree = tree_utils.read_tree(treels[i],nwk=True)
        ranges = tree_utils.read_strat(strt)
        tree_utils.match_strat(tree,ranges)
        tree_utils.init_heights_strat(tree)
        #print tree.get_newick_repr(True)
        #sys.exit(0)
        ranges = None
        strat_opt = stratoML.optim_lambda_heights(tree,ranges)
        strat_heights_tree = strat_opt[0].get_newick_repr(True)
        strat_like = -float(strat_opt[1][1])
        combined_like = strat_like + aln_liks[i]
        #print combined_like
        if best == None:
            best = (combined_like,strat_heights_tree)
            bls.append((combined_like,strat_heights_tree))
            print best
        else:
            if best[0] < combined_like:
                best = (combined_like,strat_heights_tree)
                bls.append((combined_like,strat_heights_tree))
                print bls
                print best
        if best_strat == None:
            best_strat = (strat_like,strat_heights_tree)
        else:
            if best_strat[0] < strat_like:
                best_strat = (strat_like,strat_heights_tree)
                print best_strat
        if len(aln_liks) % interval == 0:
            print str(100* prog/len(aln_liks))+"% done"
        prog+=1


    treefl = open(aln+".best_combined.tre","w")
    treefl.write(best[1])
    treefl.close()
    treefl = open(aln+".best_STRAT.tre","w")
    treefl.write(best_strat[1])
    treefl.close()
    treefl = open(aln+".searched.trees","w")
    for i in bls:
        treefl.write(str(i[0])+"\t"+i[1]+"\n")
    treefl.close()

    print "\n"
    print "#################################"
    print "----- Combined data ML tree -----"
    print "#################################"
    print best
    print "\n"
    print "#################################"
    print "----- Stratigraphic ML tree -----"
    print "#################################"
    print best_strat 

if __name__ == "__main__":
    main()
