import sys,os
import tree_utils
import stratoML
import phylo3
import tree_reader

if len(sys.argv) != 6:
    print "usage: " +sys.argv[0]+ "<file containing trees> <stratigraphic range file> <morphological alignment> <morphological partitions file> <right and left descendents of root node>"
    sys.exit(0)

def main():
    trees = sys.argv[1]
    strt = sys.argv[2]
    aln = sys.argv[3]
    parts = sys.argv[4]
    root_tax = sys.argv[5]

    #calculate likelihoods
    cmd = "iqtree-omp -s "+aln+" -st MORPH -nt 2 -spp "+parts+" -z " +trees
    os.system(cmd)
    likefl = parts+".trees"
    aln_liks = []

    #rrtreefl = open(trees+".rr","w")
    for i in open(likefl,"r"):
        spls = i.strip().split("]")
        nwk = spls[1]
        #rrtreefl.write(nwk+"\n")
        #tre = tree_reader.read_tree_string(nwk)
        #root_nodes = [i for i in tre.leaves() if i.label in root_tax]
        likespls = spls[0].strip().split("lh=")
        lik = float(likespls[1])
        aln_liks.append(lik)
    #likefl.close()
    #rrtreefl.close()
    cmd = "pxrr -t "+trees+" -g "+root_tax+ " > "+trees+".rr"
    os.system(cmd)

    rr = open(trees+".rr","r")
    treels = []
    for i in rr:
        treels.append(i.strip())
    rr.close()
    cmd = "rm "+trees+".rr"
    os.system(cmd)

    best = None
    best_strat = None
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
        if best == None:
            best = (combined_like,strat_heights_tree)
        else:
            if best[0] < combined_like:
                best = (combined_like,strat_heights_tree)
                print best
        if best_strat == None:
            best_strat = (strat_like,strat_heights_tree)
        else:
            if best_strat[0] < strat_like:
                best_strat = (strat_like,strat_heights_tree)
                print best_strat

    treefl = open(aln+"best_combined.tre","w")
    treefl.write(best[1])
    treefl.close()
    treefl = open(aln+"best_STRAT.tre","w")
    treefl.write(best_strat[1])
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

