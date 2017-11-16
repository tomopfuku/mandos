import sys,os
from mandos import *
from mandos.tree_likelihood_calculator import *
from copy import deepcopy

"""
this function is for adding a root to an SA tree where a tip has been made the root
"""

def tree_AIC(tree,log_likelihood,model_param_count):
    k = tree.nnodes("all") + model_param_count
    aic = (2*k) - (2*(log_likelihood))
    return aic

def add_root_to_SA_tree(tree):
    root = phylo3.Node()
    root.children.append(tree)
    root.label = ""
    root.height = tree.height + tree.length
    root.length = 0.0
    root.isroot = True
    tree.isroot = False
    tree.parent = root
    return root

def MK_likelihood(sitels,tree,seqs):
    site_liks = []
    for i in range(len(sitels)):
        cur_sites = sitels[i]
        cur_st = cur_sites[0]-1
        cur_end = cur_sites[1]
        #print cur_st,cur_end
        if i+1 == 1:
            state_space = 2
        else:
            state_space = i+1
        #print cur_st,cur_end
        #sites = seqs[1].seq[cur_st:cur_end]#len(seqs[0].seq)
        match_tips_and_seqs(tree,seqs)
        rmatrix,basefreq = set_multi_jc_rmatrix_basefreq(state_space,1.)
        #for n in tree.iternodes():
        #    n.length = 0.0
        curcost = calc_mult_tree_likelihood(tree,rmatrix,basefreq,range(cur_st,cur_end),False)
        res = optimize_brlen(tree,range(cur_st,cur_end),rmatrix,basefreq)
        for n in tree.iternodes():
            n.length = 0.01
        curcost = calc_mult_tree_likelihood(tree,rmatrix,basefreq,range(cur_st,cur_end),False)
        start = [0.1]
        res1 = optimize_jcrateparams(start,tree,range(cur_st,cur_end),state_space)
        rmatrix,basefreq = set_multi_jc_rmatrix_basefreq(state_space,res1[0])
        res = optimize_brlen(tree,range(cur_st,cur_end),rmatrix,basefreq)
        site_liks.append(res[1])
    return -sum(site_liks)

def main():
    if len(sys.argv) != 6:
        print "python "+sys.argv[0]+" <tree> <alignment> <taxa for testing> <partitions file> <stratigraphic ranges>"
        sys.exit(0)
    infile = open(sys.argv[1],"r")
    tree = tree_reader.read_tree_string(infile.readline())
    infile.close()
    seqs = tree_utils.read_phylip_file(sys.argv[2])
    taxa = sys.argv[3].split(",")
    parts = sys.argv[4]
    strat = tree_utils.read_strat(sys.argv[5])

    tree_utils.match_strat(tree,strat)
    tree_utils.init_heights(tree)
    tree = tree_utils.prune_SA(tree,taxa)
    like =  stratoML.hr97_loglike(tree,2.0)
    print like
    print tree_AIC(tree,like,1)
    prune_seqs = []

    for i in seqs:
        if i.name in taxa:
            #print i.get_fasta()
            prune_seqs.append(i)
    #prune_seqs[0].seq = "00000"
    #prune_seqs[1].seq = "11100"
    sitels=tree_utils.read_partition_file(parts)
    #sitels = [(1,2),(3,5)]
    #print tree.get_newick_repr(True)
    #resst = stratoML.optim_lambda(tree,strat)
    """
    print resst[1]
    print tree.nnodes("all")
    print tree.get_newick_repr(True)
    print [(i.height,i.length,i.label) for i in tree.iternodes()]
    #sys.exit(0)
    """
    #strat_lik = -resst[1][1]
    #print tree_AIC(tree,strat_lik,2)
    newtree = deepcopy(tree)
    """ 
    morpholike = MK_likelihood(sitels,newtree,prune_seqs)
    print "MK Log-likelihood: "+ str(morpholike)
    print "Strato Log-likelihood: "+str(strat_lik)
    tot = morpholike +(strat_lik)
    print "Log-likelihood: "+str(tot)
    aic = tree_AIC(tree,tot,2)#(2*6)-(2*tot)
    print "AIC: "+str(aic)
    aicc =  aic+((2*6.)*(7))/(391.-6-1)
    print "AICc: "+str(aicc)
    """
    print tree.get_newick_repr(True)
    tree_utils.collapse_non_overlapping_ranges(tree)
    tree = add_root_to_SA_tree(tree)

    #resst = stratoML.optim_lambda(tree,strat)
    print "\n\n"
    like =  stratoML.hr97_loglike(tree,2.0)
    print like
    print tree_AIC(tree,like,1)
    print tree.get_newick_repr(True)
    sys.exit(0)
    """
    print resst
    print tree.nnodes("all")
    print tree.get_newick_repr(True)
    print tree_AIC(tree,strat_lik,2)
    sys.exit(0)
    """
    strat_lik = -resst[1][1]
    newtree = deepcopy(tree)
    morpholike = MK_likelihood(sitels,newtree,prune_seqs)

    #print tree.get_newick_repr(True)
    print "MK Log-likelihood: "+ str(morpholike)
    print "Strato Log-likelihood: "+str(strat_lik)
    tot = morpholike+strat_lik
    print "Log-likelihood: "+str(tot)
    aic = tree_AIC(tree,tot,1)#2*5)-(2*tot)
    print "AIC: "+str(aic)
    aicc = aic+((2*5.)*(6))/(391.-5-1)
    print "AICc: "+str(aicc)

if __name__ == "__main__":
    main()
