# AttPrior
AttPrior - detecting attractors in a synchronous Boolean network using prior information

# Description
This is a program for detecting attractors in a synchronous Boolean network using prior information on a target (desired) attractor.
Prior information is given as the probability of 0 or 1 at each node,
where it is assumed that prior information on one global state in a periodic attractor is given when the target is a periodic attractor.
The algorithm examines states (among 2^n possible states) from higher probability ones to lower probability ones, from each of which the algorithm
traverses the trajectory of the state transition diagram until reaching to an attractor.
The algorithm stops repetition once the number of examined states exceeds a given limit.

# Installation
gcc -o attprior -lm attprior.c

# Usage
attprior BoolNetFileName PriorInfoFileName

# Input files
Two input files are needed as below.

1. BoolNet file
Description of the Boolean network in the BoolNet format.
Parentheses must not be omitted.

2. Prior information file
For each node, the following information should be provided.

- Node name, 0 or 1, probability.

Note that at most 3 different values can be specified as probabilities.
If the number of nodes is large, the probabilities for most nodes should be 1.

# Outputs
AttPrior outputs some information to stdout. However, please ignore it.
The results are given in the file "found_attractors.txt".
The meaning of this file can be easily understood.

# Example Files
toy_model.boolnet:  BoolNetFile

toy_model_initVals.txt: PriorInfoFile

# Example
$ ./attprior toy_model.boolnet toy_model_initVals.txt

Then, you can see from "found_attractors.txt"  that
this Boolean network has two periodic attractors with periods 3 and 9.

# Parameters
The code includes many parameters.

Since memory spaces are not dynamically allocated, parameter values should be appropriately adjusted depending on the size of a network, and so on, especially when using large scale networks (e.g., the number of nodes is > 100).

Basically, parameter names explain their meanings.

If you want to adjust parameters, please understand details of the code.
(Since there are almost no comment lines, I do not recommend it.)

# Reference
The manuscript related to this program is under preparation.
