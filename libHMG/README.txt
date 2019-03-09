docker build -t test .
gcc -shared -o abm.so -fPIC pyConnect.c model.c cellPopulation.c cell.c utility.c inputModel.c -lm -lgsl -lgslcblas
