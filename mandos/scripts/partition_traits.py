import sys

if len(sys.argv) != 2:
    print "usage: "+sys.argv[0]+" <alignment to be partitioned>"
    sys.exit(0)

fl = open(sys.argv[1],"r")

firstline = fl.readline()

seqs = {}
sites = {}
start = False
taxon_ls = []
for i in fl:
    spls= i.strip().split()
    if len(spls) > 2:
        print "alignments must be PHYLIP format"
        sys.exit(0)
    tax = spls[0]
    taxon_ls.append(tax)
    seq = spls[1]
    for j in range(len(seq)):
        if start == True:
            sites[j]+=seq[j]
        elif start == False:
            sites[j] = [seq[j]]
            if j == len(seq)-1:
                start = True


partitions = [[],[],[],[],[],[],[]]
for i in sites.keys():
    seen = []
    for j in sites[i]:
        if j not in seen and j != "-":
            seen.append(j)
    nstates = len(seen)
    dic = {}
    for j in range(len(seen)):
        dic[seen[j]] = j
    new_seq = ""
    for j in sites[i]:
        if j != "-":
            new_char = dic[j]
        else:
            new_char = j
        new_seq += str(new_char)   
    partitions[nstates].append(new_seq)

#print partitions[1:]
seqs = {}
for i in taxon_ls:
    seqs[i] = ""

part_fl = open(sys.argv[1]+".models","w")
last=1
total_len = 0
for i in range(len(partitions[1:])):
    if len(partitions[i]) == 0:
        continue
    part_fl.write("MK, part"+str(i)+" = "+str(last)+"-"+str( len(partitions[i])+total_len)+"\n")
    last = len(partitions[i])+1
    total_len+=len(partitions[i])
    for j in partitions[i]:
        for k in range(len(j)):
            seqs[taxon_ls[k]] += j[k]
part_fl.close()

print str(len(taxon_ls))+"\t"+str(len(seq))
for i in seqs.keys():
    #a = True
    print i+"\t"+seqs[i]
