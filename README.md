# MANDOS

**M**aking and **AN**alysing **D**endrograms from **O**ccurrences in **S**tratigraphy

# Setting up

Should be able to just run:

        python setup.py build_ext --inplace
        python setup.py build
        python setup.py install

# Reading data into MANDOS

## Trees
can read in trees from a newick file like:

        tree = mandos.tree_utils2.read_tree("examples/cetaceans/cetacean.ML.tre")

## Stratigraphic ranges

Stratigraphic range files should be tab-separated text arranged like:

        Taxon   FAD LAD num_occurrences

These can be read in like:

        ranges = mandos.tree_utils2.read_strat("examples/cetaceans/cetacean.strat")

## Character data

Discrete character matrices should be in phylip format and can be read in like:

        traits = mandos.tree_utils2.read_phylip_file("examples/cetaceans/cetacean.phy")

Handling of continuous traits is a bit of a mess right now, but will be sorted soon.

## Partitioning by state space

It is probably sensible to specify a separate substitution matrix for each state space represented in the data (eg., binary vs. 3-state). You can do this by writing a RAxML-style partitions file:

        part0,1-10
        part1,11-20

There should be a script in the scripts/ folder that can sort traits and generate this file. You can now read in these partitions like:

        sitels = mandos.tree_utils2.read_partition_file("examples/cetaceans/cetacean.phy.models")


# A couple of simple analyses

## Scaling a tree to stratigraphic ranges

Obviously, need to import

        import mandos

After reading in your tree and range data, we want to match the ranges to the tree object and initialise the heights and branch lengths:

        mandos.tree_utils2.match_strat(tree,ranges)
        mandos.tree_utils2.init_heights_strat(tree)

Note that the height initialisation currently doesn't work well with trees that have "singleton" or named ancestral nodes, so it is currently best to read in a fully bifurcating tree and collapse the nodes in MANDOS manually. This can be done like:

        mandos.tree_utils2.make_ancestor(tree, *tax_name*)

the *tax_name* can be either a single taxon or a comma-separated series

We can then optimize the node heights and calculate the stratigraphic likelihood using the preservation model of Huelsenbeck and Rannala 1997:

        opt = mandos.stratoML.optim_lambda_heights(tree,ranges)

## Calculating discrete character likelihood

After having read in character data and partitions (if applicable), we can calculate their likelihood along a tree and optimize branch lengths if we wish.

        morpholike = mandos.tree_likelihood_calculator.calc_mk_like(sitels,tree,seqs,False)

if we want to optimize branch lengths:

        morpholike = mandos.tree_likelihood_calculator.calc_mk_like(sitels,tree,seqs,True)

Branch length optimization is currently really slow. This should improve as I continue the move to Cython and optimize things a bit.

### Morphological likelihood on sampled ancestor trees

Can also calculate likelihood on an SA tree (although may want to optimize branch lengths if comparing to bifurcating, since involves a topological rearrangement)

        mandos.tree_utils2.make_ancestor(tree, *tax_name*)
        morpholike = mandos.tree_likelihood_calculator.calc_mk_like(sitels,tree,seqs,True)

If you are comparing between anagenetic and cladogenetic arrangements, you can prune out the relevant subtree to speed up optimization

        subtree = mandos.tree_utils2.prune_SA(tree,["tax1","tax2","tax3"])

The first argument should be the tree, and the second should be the list of taxa contained in the relevant subtree.


