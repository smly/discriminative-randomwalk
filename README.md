README
======

# build

    $ autoreconf -ivf
    $ ./configure
    $ make

# usage

## example

    $ ./src/dwalk -i input_filename -l label_filename -o output_prefix
    Number of nodes:           9
    Number of labels:          2
    Number of labeled nodes:   2
    Number of unlabeled nodes: 7
    calculate alpha variables: ....done
    calculate beta variables:  ....done
    prob_y[1] = 0.333333
    prob_y[2] = 0.666667
    [> betweenness:
    0.218 0.152 0.244 0.031 0.001 0.000 0.000 0.002 0.000 
    0.000 0.000 0.017 0.802 1.898 1.034 1.453 0.870 0.375 
    [> normalized:
    1.000 1.000 0.933 0.037 0.001 0.000 0.000 0.002 0.000 
    0.000 0.000 0.067 0.963 0.999 1.000 1.000 0.998 1.000 
    [> classfied result using a MAP:
    1 1 1 2 2 2 2 2 2 
    calculate alpha variables: ..........done
    calculate beta variables:  ..........done
    prob_y[1] = 0.333
    prob_y[2] = 0.667
    [> betweenness:
    0.269 0.246 0.350 0.122 0.036 0.014 0.014 0.049 0.020 
    0.065 0.065 0.151 0.385 0.569 0.361 0.352 0.310 0.150 
    [> normalized:
    0.804 0.790 0.698 0.241 0.060 0.037 0.038 0.136 0.120 
    0.196 0.210 0.302 0.759 0.940 0.963 0.962 0.864 0.880 
    [> classfied result using a MAP:
    1 1 1 2 2 2 2 2 2

## graph data

    $ cat tests/dataset1/sample.in
    2:1 3:1
    1:1 3:1
    1:1 2:1 4:1
    3:1 5:1 8:1
    4:1 6:1 7:1
    5:1 7:1
    5:1 6:1
    4:1 9:1
    8:1

## label data

    $ cat tests/dataset1/sample.label 
    .
    1
    .
    .
    .
    2
    .
    .
    2

# test

    $ autoreconf -ivf
    $ ./configure --enable-devel
    $ make check
    $ make memcheck
