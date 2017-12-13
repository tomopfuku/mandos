import tree_reader2
import tree_likelihood_calculator as tlc
from node2 import Node
from sequence import Sequence
from copy import deepcopy
import sys
#import rearrange
import numpy as np

def aic(log_likelihood,model_param_count):
    k =  model_param_count
    aic = (2*k) - (2*(log_likelihood))
    return aic

def find_node_by_label(tree,label):
    for i in tree.iternodes():
        if i.label == label:
            return i

def get_clades(tree):
    clades = []
    for node in tree.iternodes():
        clades.append([c.label for c in node.leaves()])
    clades = clades
    return clades

def get_all_tree_NNIs(tree):
    for node in tree.iternodes():
        if node.istip or len(node.descendants()) < 3 or node == tree: 
            continue
        nnis = nni_set(node,tree)
        for i in nnis:
            print i+";"

def swap_two_tips(tip_labels,tree):
    #mrca = phylo3.getMRCA(tip_labels,tree)
    n1 = find_node_by_label(tree,tip_labels[0])
    n2 = find_node_by_label(tree,tip_labels[1])
    p1 = n1.parent
    p2 = n2.parent

    p1.remove_child(n1)
    p1.add_child(n2)
    p2.remove_child(n2)
    p2.add_child(n1)

def compare_trees(tree1,tree2):
    ls1 = []
    for node in tree1.iternodes():
        ls1.append((node,node.parent))
    ls2 = []
    for node in tree2.iternodes():
        ls2.append((node,node.parent))
    if ls1 == ls2:
        return True
    else:
        return False


# this returns all possible NNIs for a single bifurcating node with bifurcating children
def nni_set(innode,inroot):
    if len(innode.children) != 2 or len(innode.descendants()) < 3:
        print "this only works on bifurcating nodes that parent multiple subtrees (ie. does not lead to only terminal edges)"
        return None
    node = deepcopy(innode)
    root = deepcopy(inroot)
    subtrees = []

    done = False
    for each_node in root.iternodes():
        if [tip.label for tip in each_node.leaves()] != [tip.label for tip in node.leaves()]:
            continue
        for child in each_node.children:
            if child.istip == False:
                assert len(child.children) == 2
                for sub in child.children:
                 subtrees.append(sub)
        done = True
        single_node = each_node
        break

    subtrees += [i for i in single_node.children if i.istip] #add terminal subtree child --> 'c' in (a,b),c)) 
    assert len(subtrees) == 3 or len(subtrees) == 4

    nni_trees = []
    for c1 in subtrees:
        for c2 in subtrees:
            p1 = c1.parent
            p2 = c2.parent
            if c1 == c2 or p1 == p2: #can't swap subtrees with same parent 
                continue
            p1.remove_child(c1)
            p1.add_child(c2)
            p2.remove_child(c2)
            p2.add_child(c1)
            c1.parent = p2  #swap subtrees
            c2.parent = p1
            nni_trees.append(root.get_newick_repr())
    nni_trees = list(set(nni_trees)) #remove duplicates
    #print len(nni_trees)
    return nni_trees

def read_phylip_file(infilename):
    infile = open(infilename,"r")
    infile.readline()
    seqls = []
    for i in infile:
        ls = i.strip().split()
        tax = ls[0]
        seq = ls[1]
        tseq = Sequence()
        tseq.name = tax
        tseq.seq = seq
        seqls.append(tseq)
    return seqls


def prune_SA(tree,taxa): #taxa should be list of all taxa present in the subtree to be tested
    for i in tree.iternodes():
        if i.istip:
            continue
        if set([j.label for j in i.leaves()]) == set(taxa):
            for k in i.parent.children:
                if k != i:
                    i.parent.remove_child(k)
                    i.parent = None
            i.length = 0.0
            #add_root_to_SA_tree(i)
            return i#.parent

def make_ancestor(tree,tax_label):
    taxa = tax_label.split(",")
    for t in taxa:
        for node in tree.iternodes():
            if node.label == t:
                par = node.parent
                par.label = node.label
                par.upper = node.upper
                par.lower = node.lower
                if par.parent == None: # need to make a new root if the new SA is outgroup 
                    newroot = Node()
                    newroot.add_child(par)
                    tree = newroot
                    par.length += node.length
                    par.height = node.height
                    par.num_occurrences = node.num_occurrences
                    par.remove_child(node)
                    par.parent.height = par.height + par.length
                    par.istip = True
                else:
                    par.length += node.length
                    par.height = node.height
                    par.num_occurrences = node.num_occurrences
                    old_length = node.length
                    par.remove_child(node)
                    par.istip = True
                    par.children[0].length -= old_length 
    #print newroot.get_newick_repr(True) 
    return tree         


def read_partition_file(fl): #reads RAxML style partition file
    parts = open(fl,"r")
    sites = {}
    for i in parts:
        ls = i.strip().split("=")
        lss = ls[1].strip().split("-")
        numstates = int(ls[0].strip().split(",")[1].strip())
        st = int(lss[0])
        end = int(lss[1])
        sites[numstates] = (st,end)
        #sites.append((st,end))
    parts.close()
    return sites

def init_heights_strat(tree,fixed_root=False):
    for i in tree.iternodes(order=1):
        if i == tree and fixed_root == True:
            continue
        elif i.istip and i.num_occurrences != 0: #Set heights for tips
            if i.lower != 0.0:
                i.height = i.lower - 0.01
            else:
                i.height = i.lower
        elif i.istip == False and i.label == "":
            o = []
            start = False
            for j in i.children:
                if j.num_occurrences != 0:
                    if start == False:
                        o = [j.upper,j.lower]
                        #print o
                        start = True
                    else:
                        o = o+[j.upper,j.lower]
            if len(o) > 0:
                i.height = max(o+[j.height for j in i.children]) + 0.01
            else:
                i.height = max([j.height for j in i.children])+0.01
    for i in tree.iternodes():
        if i == tree:
            continue
        i.length = i.parent.height-i.height

def assign_node_nums(tree,tips = True,fixed_root = False):
    num = 0
    for i in tree.iternodes(order=0):
        if i.istip and tips == False:
            continue
        elif i == tree and fixed_root == True:
            continue
        else:
            i.number = num
            num += 1


def assign_branch_nums(tree):
    num = 0
    for i in tree.iternodes(order=0):
        if i == tree:
            continue
        i.number = num
        num +=1

def assign_brlens(l,tree):
    count = 0
    for i in tree.iternodes(order = 0):
        if i.parent == None:
            continue
        if l[i.number] < 0.00001 :
            return True
        i.length = l[i.number]
    return False

def assign_node_heights(h,tree,fixed_root=False):
    count = 0
    for i in tree.iternodes(order=0):
        i.height = h[count]
        count += 1
        if i.parent == None: #set height+length for root
            continue
        else:
            i.length = i.parent.height-i.height
        if i.length < 0:
            print "h"
            return True
        if i.num_occurrences != 0:
            if i.height > i.occurrences[1] or i.height + i.length < max(i.occurrences):
                print i.label,i.height,i.occurrences,i.parent.height
                return True
    return False

def tip_dates(tree,dates,root_height):
    d = {}
    for i in open(dates,"r"):
        spls = i.strip().split()
        d[spls[0]] = float(spls[1])
    for i in tree.iternodes():
        if i.istip == True:
            i.height = d[i.label]
        elif i.parent==None:
            i.height = root_height

def match_traits_tips(tree,traits,number):
    for i in tree.leaves():
        i.charst = traits[i.label][number]

def read_strat(stratfl):
    ran = {}
    for i in open(stratfl,"r"):
        spls = i.strip().split("\t")
        t = spls[0]
        s = [float(i) for i in spls[1:-1] if i != "NA"]
        if s[1] > s[0]:
            print "ERROR: Sampling range data is incorrectly formatted.\n"
            sys.exit(0)
        n = float(spls[-1])
        if s == []:
            ran[t] = [-1.,-1.]
        else:
            ran[t]=s
            ran[t].append(n)
    return ran

def match_strat(tree,strat):
    for i in tree.iternodes(order=1):
        if i.label != "":
            i.upper = strat[i.label][0]
            i.lower = strat[i.label][1]
            i.num_occurrences = strat[i.label][2]


def read_tree(in_tree,nwk=False):
    if nwk == False:
        nwk = open(in_tree,"r").readlines()[0].strip()
    else:
        nwk = in_tree
    tree = tree_reader2.read_tree_string(nwk)
    for i in tree.iternodes():
        i.old_length = i.length
    return tree

def read_continuous(traitfl,tree): #should be phylip formatted
    traits = {}
    op = open(traitfl,"r")
    firstline = op.readline()
    fspls = firstline.strip().split()
    try:
        ntax = int(fspls[0])
        ntraits = int(fspls[1])
    except:
        print "ERROR: trouble reading the header. continuous traits should be phylip formatted and tab separated"
        print "can't read:"
        print firstline
        sys.exit(0)
    for i in op:
        if i == "":
            continue
        spls = i.strip().split("\t")
        nm = spls[0]
        traitls = [] 
        for j in spls[1:]:
            try:
                cur = float(j)
                traitls.append(j)
            except:
                if j == "?" or j == "-":
                    traitls.append(1000000000.0)
                else:
                    print "something is wrong with the data. these traits cannot be read as continuous."
                    print "cannot read starting at line:"
                    print i
                    sys.exit(0)
        
        trait_array = np.array(traitls,dtype=np.double)
        traits[nm]= trait_array
    for node in tree.leaves():
        node.cont_traits = traits[node.label]
    return (ntax,ntraits)
    #return traits

def prune_tip(tree,tip_name):
    for node in tree.iternodes():
        if node.label != tip_name:
            continue
        if len(node.children) == 0: #node is not a SA
            if node.parent == tree:
                par = node.parent
                for c in par.children:
                    if c == node:
                        par.remove_child(c)
                if len(tree.children) == 1:
                    cc = tree.children[0]
                    for c in cc.children:
                        tree.add_child(c)
                    tree.remove_child(cc)
                return tree
            else:
                par = node.parent
                gpar = par.parent
                for c in par.children:
                    if c != node:
                        sib = c
                gpar.add_child(sib)
                sib.length += par.length
                gpar.remove_child(par)
                return tree
        elif len(node.children) > 0:
            print "you haven't yet implemented continuous trait likelhoods for SA trees!!"
            return tree

def prune_node(tree,node):
        if len(node.children) == 0: #node is not a SA
            if node.parent == tree:
                par = node.parent
                for c in par.children:
                    if c == node:
                        par.remove_child(c)
                if len(tree.children) == 1:
                    cc = tree.children[0]
                    for c in cc.children:
                        tree.add_child(c)
                    tree.remove_child(cc)
                return tree
            else:
                par = node.parent
                gpar = par.parent
                for c in par.children:
                    if c != node:
                        sib = c
                gpar.add_child(sib)
                sib.length += par.length
                gpar.remove_child(par)
                return tree
        elif len(node.children) > 0:
            print "you haven't yet implemented continuous trait likelhoods for SA trees!!"
            return tree


def assign_sigsq(tree,sigsq=[1.0]):
    for i in tree.iternodes():
        i.sigsq = sigsq[0]
    return tree
