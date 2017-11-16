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

Coming soon!

# A couple simple analyses

## Scaling a tree to stratigraphic ranges

Obviously, need to import

        import mandos

After reading in your tree and range data, we want to match the ranges to the tree object and initialise the heights and branch lengths:

        mandos.tree_utils2.match_strat(tree,ranges)
        mandos.tree_utils2.init_heights_strat(tree)

Note that the height initialisation currently doesn't work well with trees that have "singleton" or named ancestral nodes, so it is currently best to read in a fully bifurcating tree and collapse the nodes in MANDOS manually. This can be done like:

        mandos.tree_utils2.make_ancestor(tree, *tax_name*)

the *tax_name* can be either a single taxon or a comma-separated series

We can then optimize the node heights using the preservation model of Huelsenbeck and Rannala 1997:

        opt = mandos.stratoML.optim_lambda_heights(tree,ranges)


