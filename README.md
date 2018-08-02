# HMG
Heterogeneous Multiphasic Growth (HMG) model for the simulation of bacterial populations growth in a well-mixed environment

## Introduction

This project was part of a bioengenerring thesis on the exploration of different growth stategies for the accurate simulation of bacterial population in a wider range of growth conditions that the classic Cooper-Helmstetter cell cycle model of growth allows for. The core model is written in ANSI C as libHMG. It contains a mechanistic description of the cell cycle with, among others, the critical mass, eclipse period and others.... as described by the CH model. 

For a more detailed explanation of the particularities of this agent based (or individual based) bacterial model please refer to the following open access paper: 

[Du Lac, Melchior, et al. "Predicting the dynamics and heterogeneity of genomic DNA content within bacterial populations across variable growth regimes." ACS synthetic biology (2016).](https://doi.org/10.1021/acssynbio.5b00217)

Or to the following thesis:

TO BE ADDED

Although the paper and thesis have very particular research purposes, the core model may be used for other purposes, since it mainly follows the classic CH model as described in:

[Abner, Kristo, et al. "Single-cell model of prokaryotic cell cycle." Journal of theoretical biology 341 (2014): 78-87. APA](https://doi.org/10.1016/j.jtbi.2013.09.035)	

## Getting Started

### Prerequisites

To compile libHMG the following packages are required:

* [GSL](https://www.gnu.org/software/gsl/) - To run ODE models in parralell to the cell cycle mechanistic model
* [DEAP](https://github.com/deap) - Distributed Evolutionary Algorithms in Python
* [Numpy](http://www.numpy.org/) - Python scientific packages
* [Scipy](https://www.scipy.org/) - Python scientific packages 
* [Sobol seq](https://pypi.python.org/pypi/sobol_seq/0.1.2) - Sobol sequence generator in python

### Compiling

To compile the libHMG model alone with the example in the main() function in model.c, the following command must be used (assuming that you are compiling the model under a UNIX environment):

```
cd libHMG
gcc -o test model.c cellPopulation.c cell.c utility.c inputModel.c -lm -lgsl -lgslcblas
```

where -lgsl and -lgslcblas are flags for the GSL package and -lm is the flag for the math C packages

To compile libHMG to be ran with the python wrapper functions, it must be compiled as follows:

```
cd libHMG
gcc -shared -o abm.so -fPIC pyConnect.c model.c cellPopulation.c cell.c utility.c inputModel.c -lm -lgsl -lgslcblas
```

For a more complete description of libHMG please refer to the doc in the libHMG folder.

## Running the tests

TODO

## Built With

* [polynomial](https://github.com/wafo-project/pywafo/blob/master/wafo/polynomial.py) - From the pyWafo project

## Authors

Melchior du Lac

## License

GPLv3

## Acknowledgments

* Professor Declan G. Bates
* Associate Professor Joshua N. Leonard
* Dr Andrew H. Scarpelli
* Dr Andrew K. D. Younger
