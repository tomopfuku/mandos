# MANDOS

**M**aking and **AN**alysing **D**endrograms from **O**ccurrences of **S**tratigraphy

(it is a stretch, but I like Tolkien, so deal with it)


# Setting up

Pull the repo from GitHub in a sensible directory:

        git clone https://github.com/carolinetomo/mandos.git

After navigating to the main directory, should be able to just run (Linux or Mac):

        python setup.py build_ext --inplace
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

If a taxon occurs at only one site/locality/whatever, you can set the FAD and LAD to the same, and num_occurrences to 1.

## Character data

Discrete character matrices should be in phylip format and can be read in like:

        traits = mandos.tree_utils2.read_phylip_file("examples/cetaceans/cetacean.phy")

Handling of continuous traits is a bit of a mess right now, but will be sorted soon.

## Partitioning by state space

It is probably sensible to specify a separate substitution matrix for each state space represented in the data (eg., binary vs. 3-state). You can do this by writing a RAxML-style partitions file:

        MODEL, NSTATES, START-STOP

for example:
        
        MULTI, 2 = 1-10
        MULTI, 3 = 11-20

Where the model should follow RAxML nomenclature (usually "MULTI" for discrete morphological data), and the NSTATES should be the 'state-space' of the partiton (e.g., "2" for binary characters). There should be a script in the scripts/ folder that can sort traits and generate this file. You can now read in these partitions like:

        sitels = mandos.tree_utils2.read_partition_file("examples/cetaceans/cetacean.phy.models")


# A couple of simple analyses

## Scaling a tree to stratigraphic ranges

Obviously, need to import

        import mandos

After reading in your tree and range data, we want to match the ranges to the tree object and scale the branch lengths to the ranges provided in the file:

        mandos.tree_utils2.match_strat(tree,ranges)
        mandos.tree_utils2.init_heights_strat(tree)

We can then optimize the node heights and calculate the stratigraphic likelihood using the preservation model of Huelsenbeck and Rannala (1997):

        opt = mandos.stratoML.optim_lambda_heights(tree,ranges)

### Stratigraphic likelihood on sampled ancestor trees

We might want to calculate the stratigraphic likelhood on a tree that containes direct ancestor-descendant relationships. Note that the height initialisation currently doesn't work well with these types of trees, so it is currently best to read in a fully bifurcating tree and collapse the nodes in MANDOS manually. This can be done like:

        mandos.tree_utils2.make_ancestor(tree, *tax_name*)

the *tax_name* can be either a single taxon or a comma-separated series (entered as a string, not a python list). We can then just calculate the likelihood whilst optimizing heights like above.


## Calculating discrete character likelihood

After having read in character data and partitions (if applicable), we can calculate their likelihood along a tree and optimize branch lengths if we wish.

        morpholike = mandos.tree_likelihood_calculator.calc_mk_like(sitels,tree,seqs,False)

if we want to optimize branch lengths:

        morpholike = mandos.tree_likelihood_calculator.calc_mk_like(sitels,tree,seqs,True)

Branch length optimization is currently really slow. This should improve as I continue the move to Cython and optimize things a bit. Branch length optimization is also a difficult task, and so the optimizers don't always converge upon the true global optimum. If using this feature for comparative tests, you would be well-advised to do a couple runs and take the best to make sure your likelihood is really as good as it can be (it often isnt!). 

For good measure, here are the arguments:

        calc_mk_like(sitels,tree,seqs,optimize=False,random_start = True,numstates = None)

You always need to specify a partitions file, but if you want to analyse all of the traits under a single transition matrix, you can just specify the maximum number of states present in the matrix under the numstates argument and place all of the traits into a single partition. 

### Morphological likelihood on sampled ancestor trees

We can also calculate the character likelihood on a tree containing sampled (direct) ancestors. You may want to optimize branch lengths if comparing to bifurcating, since this involves a topological rearrangement. There are more efficient ways of doing this than are currently implemented here, but can run the following for now:

        mandos.tree_utils2.make_ancestor(tree, *tax_name*)
        morpholike = mandos.tree_likelihood_calculator.calc_mk_like(sitels,tree,seqs,True)

If you are comparing between anagenetic and cladogenetic arrangements, you can prune out the relevant subtree to speed up the test

        subtree = mandos.tree_utils2.prune_SA(tree,["tax1","tax2","tax3"])

The first argument should be the tree, and the second should be the list of taxa contained in the relevant subtree.

Can then calculate the likelihood and optimize branch lengths:

        morpholike = mandos.tree_likelihood_calculator.calc_mk_like(sitels,subtree,seqs,True)

When weighing between cladogenetic and anagenetic trees, we need to use something like AIC to account for the difference in parameters:

        print mandos.tree_utils2.tree_AIC(subtree,morpholike,1)

the last argument should be the number of parameters (_k_). _k_ = 1 for each partition in the analyses when using Mk or single rate Brownian motion. 

![alt text](examples/likelihood.png "likelihood calculation on a single character")

## Combining character and stratigraphic data

We might want to evaluate trees using both stratigraphic and character data.

When combining, the number of parameters (_k_) should be
        
        num_nodes+num_tips + num_branch_lengths + morph_parameters + strat_parameters

Note that this is _kind_ of like tip-dating, but the dates only reflect observed stratigraphic ranges, since the ML under the implemented model node heights move as close to the FADs and LADs as is possible given a particular tree shape.


