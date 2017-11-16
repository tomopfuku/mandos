import tree_reader2,tree_utils2
import sys
from numpy.random import exponential

def sim_occurrences(tree,r):
    occurrences = {}
    for i in tree.iternodes():
        if i.istip:
            cur_length = i.length
            occurrences[i.label] = []
            i.old_length = i.length
            while (cur_length > 0):
                cur_length -= exponential(r)
                cur_time = i.height + cur_length
                if cur_length > 0:
                    occurrences[i.label].append(cur_time)
            if len(occurrences[i.label]) == 0:
                occurrences[i.label].append(i.height)
            elif i.height == 0.:
                occurrences[i.label].append(i.height)
            i.length = i.old_length
    return occurrences

if __name__ == "__main__": #run a test
    fl = open("/home/carolioness/documents/coding_stuff/mandos/examples/hominins/hominin.best_combined.tre","r")
    nwk = fl.readline()
    nwk = nwk.strip()
    tree = tree_reader2.read_tree_string(nwk)
    strat = tree_utils2.read_strat("/home/carolioness/documents/coding_stuff/mandos/examples/hominins/hom_occurrences.csv")
    tree_utils2.match_strat(tree,strat)
    tree_utils2.init_heights(tree)
    reps = 100
    for i in range(reps):
        fos = sim_occurrences(tree,0.2)
        outfl = open("simulated_ranges/"+str(i)+".fos","w")
        for i in fos.keys():
            cur_upper = max(fos[i])
            cur_lower = min(fos[i])
            cur_num_fos = len(fos[i])
            outfl.write(i+"\t"+str(cur_upper)+"\t"+str(cur_lower)+"\t"+str(cur_num_fos)+"\n")


