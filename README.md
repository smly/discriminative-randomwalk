dwalk
======

dwalk -- an implementation of discriminative random walk ([Callut+, ECML/PKDD 2008])

# Overview

D-walk is a novel technique to tackle semi-supervised classification
 problems in large graphs, proposed in [Callut et al. ECML/PKDD 2008].
D-walks use betweenness measure based on passage times during random walks
 of bounded lengths. Forward and backward recurrences are derived to
 efficiently compute the passage times.

# Format of Input Data

(to be appear)

# Install

    $ ./configure
    $ make

# Usage

    $ ./dwalk -i inputfile -l labelfile -o output_prefix -w bounded_length

## Example

    $ ./dwalk -i input_filename -l label_filename -o output_prefix -w 5
    Number of nodes:           9                                                             
    Number of labels:          2
    Number of labeled nodes:   3
    Number of unlabeled nodes: 6
    calculate alpha variables: .....done
    calculate beta variables:  .....done
    calc bounded dwalks betweeness

    $ cat input_filename
    2:1 3:1
    1:1 3:1
    1:1 2:1 4:1
    3:1 5:1 8:1
    4:1 6:1 7:1
    5:1 7:1
    5:1 6:1
    4:1 9:1
    8:1

    $ cat label_filename
    .
    1
    .
    .
    .
    2
    .
    .
    2

    $ cat output_prefix.weight
    1.000000000 0.000000000
    1.000000000 0.000000000
    1.000000000 0.000000000
    0.230733215 0.769266785
    0.000000000 1.000000000
    0.000000000 1.000000000
    0.000000000 1.000000000
    0.000000000 1.000000000
    0.000000000 1.000000000

    $ cat output_prefix.result
    1
    1
    1
    2
    2
    2
    2
    2
    2

# Reference

Jerome Callut, Kevin Francoisse, Marco Saerens, Pierre Dupont, "Semi-supervised Classification from Discriminative Random Walks", ECML/PKDD 2008.

# Author

  Kohei Ozaki <kohei-o@is.naist.jp>
  Nara Institute of Science and Technology, 
  Graduate School of Information Science, 
  Computational Linguistics Laboratory 

# License

BSD3 License
