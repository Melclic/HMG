/**
* @file cellPopulation.h
*
* Cell Population header file that grows the population
*
* @version 1.0
* @author Melchior du Lac
*
*/

#ifndef CELL_POPULATION_H
#define CELL_POPULATION_H

#include "cell.h"
#include "utility.h"
#include "inputModel.h"

typedef struct CELL_POPULATION
{
	int maxCells; /**< max number of dead and alive cells available according to the dedicated RAM*/
	int numCells; /**< Number of cells that are alive*/
	int indexArray; /**< Shows the end of the used array*/
	int freeIndex; /**< Index of a position of a dead cell*/
	Cell * cellArray; /**< Array of cells*/
	float Vi; /**< Mass at initiation*/
	float Vi_plasmid; /**< Mass at initiation*/
	float cNoise; /**< C Gaussian noise*/
	float dNoise; /**< D Gaussion noise*/
	float ViNoise; /**< Critical mass Gaussian noise*/
	float VaNoise; /**< Volume growth Gaussian noise*/
	float chanceInit; /**< Probability term that opens oriC and replication forks*/
	float divNoise; /**< Volume symmetry at division mean Gaussian noise*/
	float divRatio; /**< Ratio between volume symmetry at division mean Gaussian noise (default==0.5)*/
	float partRatio; /**< Chromosome symmetry at division mean Gaussian noise (default==0.5)*/
	float partNoise; /**< Chromosome symmetry at division standard deviation Gaussian noise*/
	float chromDeg; /**< Probability term that determines if the chromosome experiences fatal DNA damage that degreades it. True if >= random number generated from a Gaussian distribution of mean 0 and standard deviation of 1.*/
	float repForkDeg; /**< Probability term that determines if the replicating strand experiences DNA damage that leads to its collapse. True if >= random number generated from a Gaussian distribution of mean 0 and standard deiation of 1.*/
	int numFrozenCells; /**< Count of the number of cells that have reached maximal chromosome number (64)*/
	int numAnucleateCells; /**< Count of the number of cells that have all their chromosomes degraded*/
	int numModelParams; /**< Number of params for the modelInitialsParams*/
	int numModelSpecies; /**< Number of species for the modelInitialSpecies array*/
	int numModelGenes; /**< Number of species for the modelInitialSpecies array*/
	float tau;
	float C1; /**< functional form of the replication rate (C)*/
	float C2; /**< functional form of the replication rate (C)*/
	float C3; /**< functional form of the replication rate (C)*/
	float D1; /**< functional form of the segregation rate (D)*/
	float D2; /**< functional form of the segregation rate (D)*/
	float D3; /**< functional form of the segregation rate (D)*/
        float * totalVolumes; /**< Total population volume changes time series (based on dt)*/
        int lenTotalV; /**< Length of the totalVolumes array*/
	//ODE model stuff
	//double modelInitialParams[8]; /**< Initial params of the model*/
	//double modelInitialSpecies[15]; /**< Initial species of the model*/
	//float modelGeneLocations[3]; /**< Location on the chromosome (between 0.0 and 100.0, expressed in %) of the gene in the model*/
	//int modelGeneLRPos[3]; /**< Position of the gene on either the left or right hand side of the chromosome; where 0 -> right, 1 -> left*/
	//int modelGeneParamsLocations[3]; /**< Location of the parameteris in the model (note not species) that describes the copy number of a gene in particular */
	double modelInitialParams[8]; /**< Initial params of the model*/
	double modelInitialSpecies[15]; /**< Initial species of the model*/
	float modelGeneLocations[3]; /**< Location on the chromosome (between 0.0 and 100.0, expressed in %) of the gene in the model*/
	int modelGeneLRPos[3]; /**< Position of the gene on either the left or right hand side of the chromosome; where 0 -> right, 1 -> left*/
	int modelGeneParamsLocations[3]; /**< Location of the parameteris in the model (note not species) that describes the copy number of a gene in particular */
	//double modelInitialParams[NUM_MODELPARAMS]; /**< Initial params of the model*/
	//double modelInitialSpecies[NUM_MODELSPECIES]; /**< Initial species of the model*/
	//float modelGeneLocations[NUM_MODELGENES]; /**< Location on the chromosome (between 0.0 and 100.0, expressed in %) of the gene in the model*/
	//int modelGeneLRPos[NUM_MODELGENES]; /**< Position of the gene on either the left or right hand side of the chromosome; where 0 -> right, 1 -> left*/
	//int modelGeneParamsLocations[NUM_MODELGENES]; /**< Location of the parameteris in the model (note not species) that describes the copy number of a gene in particular */
} CellPopulation;

int growCells(CellPopulation * cellPopulation,
                float dt,
		float volumeAdd,
		bool isDrugTreat);

void cellDebug(Cell * cellArray, int index);

#endif
