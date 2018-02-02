#from libcpp cimport bool
from copy import deepcopy
import numpy as np

cdef:
    unsigned int PREORDER = 0
    unsigned int POSTORDER = 1
    #unsigned int BRANCHLENGTH = 0
    #unsiINTERNODES = 1

cdef class Node:
    cdef public dict data
    cdef public bint isroot,istip,recalculate,marked,visited
    cdef public double height,length,old_length,upper,lower,charst,sigsq
    cdef public unsigned int nchildren, number,num_occurrences,rate_class
    cdef public str label
    cdef public object parent
    cdef public list children
    cdef public double[:] cont_traits
    cdef public double contrast_length

    def __init__(self):
        self.data = {}
        self.isroot = False
        self.istip = False
        self.label = ""
        self.length = 0.
        self.old_length = 0.
        self.parent = None
        self.children = []
        self.nchildren = 0
        #self.comment = None
        self.charst = 0.
        self.sigsq = 0.
        self.rate_class = 0
        self.height = 0.0
        self.number = 0
        self.upper = 0.0
        self.lower = 0.0
        self.num_occurrences = 0
        self.cont_traits = np.array([],dtype=np.double)
        self.contrast_length = 0.0
        self.recalculate = True
        self.marked = False
        self.visited = False

    def get_newick_repr(self,showbl=False,show_rate=False):
        ret = ""
        for i in range(len(self.children)):
            if i == 0:
                ret += "("
            ret += self.children[i].get_newick_repr(showbl,show_rate)
            if i == len(self.children)-1:
                ret += ")"
            else:
                ret += ","
        if self.label != None:
            ret += self.label
        if showbl == True:
            ret += ":" + str(self.length)
        if show_rate ==True:
            ret += ":" + str(self.sigsq)
        return ret

    def add_child(self, object child):
        assert child not in self.children
        self.children.append(child)
        child.parent = self
        self.nchildren += 1

    def remove_child(self, object child):
        assert child in self.children
        self.children.remove(child)
        child.parent = None
        self.nchildren -= 1


    def prune_from_node(self):
        for i in self.descendants("POSTORDER"):
            if len(self.children) == 0:
                self.prune()

    def leaves(self):
        return [ n for n in self.iternodes() if n.istip ]

    def iternodes(self, unsigned int order=PREORDER):#, v=None):
        cdef:
            object child, d

        if order == PREORDER:
            yield self
        #print [i.label for i in self.children]
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order == POSTORDER:
            yield self


    def descendants(self, order=PREORDER, v=None):
        if v is None:
            v = []
        #assert order in ("PREORDER", "POSTORDER")
        for child in self.children:
            if order == PREORDER:
                v.append(child)
            else:
                v.insert(0, child)
            if child.children:
                child.descendants(order, v)
        return v

    cpdef bint update_brlens_all(self, double[:] mv):
        cdef:
            object node
            unsigned int count
            double min_bound = 0.000001
            double i,curlen
            bint fixed_tip = True

        #check variables to see that they are >0
        for i in mv:
            if i < min_bound:
                return True

        count = 0
        for node in self.iternodes():
            if node == self:
                continue
            elif fixed_tip == True and node.istip == True: #return LARGE if first tip not fixed
                fixed_tip = False
                continue
            curlen = mv[count]
            node.length = curlen
            node.old_length = curlen
            count += 1
        return False

    cpdef bint update_brlens_single(self, double[:] mv, unsigned int node_num):
        cdef:
            object node
            unsigned int count

        count = 0
        if mv[0] < 0.000001:
            return True
        for node in self.iternodes():
            if count == node_num:
                node.length = mv[0]
                node.old_length= mv[0]
            count += 1
        return False

    cpdef bint update_child_brlens(self, double[:] mv):
        cdef:
            object node
            unsigned int count
            double i

        for i in mv:
            if i < 0.00001:
                return True
        count = 0
        for node in self.children:
            node.length = mv[count]
            node.old_length= mv[count]
            count += 1
        return False

    def find_descendant(self, label):
        if label == self.label:
            return self
        else:
            for child in self.children:
                n = child.find_descendant(label)
                if n:
                    return n
        return None

    def prune(self):
        p = self.parent
        if p:
            p.remove_child(self)
        return p

    def tip_labels(self):
        labs = []
        for i in self.leaves():
            labs.append(i.label)
        return labs

    def nnodes(self, type="internal"):
        n = 0
        if type == "internal":
            for i in self.iternodes():
                if i.istip or i == self:
                    continue
                n += 1
        elif type == "all":
            for i in self.iternodes():
                n+=1
        elif type == "tips":
            for i in self.iternodes():
                if i.istip:
                    n+=1
        return n

    #this changes the trifurcating 'root' of an unrooted tree to self
    cpdef object reroot(self,oldroot):
        cdef:
            list v
            object newpar,newroot,node
            unsigned int i

        if oldroot == self:
            return self

        self.isroot = True
        oldroot.isroot = False
        v = []
        node = self
        #newroot = self

        while 1:
            v.append(node)
            if not node.parent:
                break
            node = node.parent
        v.reverse()

        newpar = None
        for i in range(len(v[:-1])):
            newpar = v[i+1]
            curnode = v[i]
            curnode.remove_child(newpar)
            newpar.add_child(curnode)
            curnode.length = newpar.length
        self.length = 0.0
        return self

    cpdef object unroot_tree(self, object oldroot):
        cdef:
            object s

        s = [i for i in self.children if i!=oldroot][0]
        self.remove_child(oldroot)
        self.remove_child(s)
        oldroot.add_child(s)
        s.length+=oldroot.length
        oldroot.length = 0.
        return oldroot

    cpdef tuple root_tree(self):
        cdef:
            object newroot,c0

        newroot = Node()
        c0 = self.children[2]
        self.remove_child(c0)
        newroot.add_child(c0)
        newroot.add_child(self)
        c0.length = c0.length/2
        self.length = c0.length
        return (newroot,self)


    """
    # this returns all possible NNIs for a single bifurcating node with bifurcating children
    # tree should probably be deep copied before using this
    def nni_set(self):
        if len(self.children) != 2 or len(self.descendants()) < 3:
            print "this only works on bifurcating selfs that parent multiple subtrees (ie. does not lead to only terminal edges)"
            return None

        subtrees = []

        for child in self.children:
            if child.istip == False:
                assert len(child.children) == 2
                for sub in child.children:
                    subtrees.append(sub)

        subtrees += [i for i in self.children if i.istip] #add terminal subtree child --> 'c' in (a,b),c))
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
                nni_trees.append(self.get_newick_repr())

        nni_trees = list(set(nni_trees)) #remove duplicates
        #print len(nni_trees)
        return nni_trees

    """
