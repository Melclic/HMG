/*
 *   libHMG - Individual based model of the a bacterial population
 *   Copyright (C) 2017  Melchior du Lac
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software Foundation,
 *   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

/**
 * @file model.c
 *
 * Model file with all the methods that initiate and run the model.
 * @version 1.0
 * @author Melchior du Lac
 *
 * 26/10/16
 * -> replaced the exit terms with the standard macro of EXIT_FAILURE and EXIT_SUCCESS
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#include "cell.h"
#include "cellPopulation.h"
#include "model.h"
#include "utility.h"
#include "inputModel.h"


//FILE *f; <---- tmp need to find a better way to output this instead of hardcoding it

//#################################################################################
//################################## ORGANISE ####################################
//#################################################################################

/**
* @brief Re-organise the cellArray so as to speed up the looping process
*
* Re-organise the array of cell by forward looping through it and finding an empty positions.
* If empty position is found then reverse loop to find the last position where a cell is not dead and swap their positions.
* Set the indexArray and freeIndex since the array has been reorganised
*
* @param model Model object
* @return Boolean of success of function
*/
int organiseCellArray(Model * model)
{
	int i, y;
	//loop through all the cells in the sample (could make this more eficient by using model->cellPopulation->numCells)
	for(i=0; i<model->cellPopulation->maxCells; i++)
	{
		//if a dead cell is found loop backwards and find a cell that is not dead and
		//replace the dead cell at poisition i with the live cell at position y
		if(model->cellPopulation->cellArray[i].isDead==true)
		{
			for(y=model->cellPopulation->maxCells-1; y>=0; y--)
			{
				//if the next dead cell position does not find a live cell in front of itself return 0
				if(i>=y)
				{
					model->cellPopulation->indexArray = i-1;
					model->cellPopulation->freeIndex = i;
					return 0;
				}
				//if you find a cell that is not dead in front of position i copy the cell[y] to cell[i] and delete cell[y]
				if(model->cellPopulation->cellArray[y].isDead==false)	
				{
					model->cellPopulation->cellArray[i] = model->cellPopulation->cellArray[y];
					constructCell(model->cellPopulation->cellArray, y);					
					break;
				}
			}		
		}
	}
	//this implies that there are no dead cells
	return 1;
}

//#################################################################################
//################################## SETTERS ####################################
//#################################################################################

/**
* @brief Initiate a model object and all the downstream structures and cellArray.
*
* Initiate a model object and all the downstream structures and cellArray. This includes all the dynamic memory allocation that is only freed when cleanModel() function is called. Initiate and return a model and seed the random start 
*
* @param maxCells Maximal number of cells that simulator may have
* @return Model object
*/
Model * initModel(int maxCells)
{
	Model * model = (Model *)(calloc(1, sizeof(Model)));
	model->cellPopulation = (CellPopulation *)(calloc(1, sizeof(CellPopulation)));
	model->cellPopulation->maxCells = maxCells;
	//############# cell array ##############
	model->cellPopulation->cellArray = (Cell*)(calloc(model->cellPopulation->maxCells, sizeof(Cell)));
	//check that the assignement of the dynamic array is successfull
	if(model->cellPopulation->cellArray==NULL)
	{
		printf("ERROR (initModel): Could not allocate bytes for cellArray\n");
		exit(EXIT_FAILURE);
	}

	//##################### seed random ###############
	if(seedRandom())
	{
		printf("ERROR (initModel): seedRandom has failed.\n");	
		exit(EXIT_FAILURE);
	}

	return model;
}

/**
* @brief Public function to pass the input parameters to the model.
*
* Public function to pass the input parameters to the model. These input parameters are by convention parameters that are constant throughout the simulation.
* We start at this point with a single cell that has not yet been incoluted. The reasoning behind this layout is the future application of a user defined number of cells to inoculate the model. Need to implement the random initiation.   
*
* @param model Model object 
* @param tau Doubling rate of the population (only applicable to exponetially growing populations)
* @param cNoise Gaussian noise standard deviation for a population associated with replication time (C)
* @param dNoise Gaussian noise standard deviation for a population associated with segregation time (D)
* @param Vi Volume at initiation. Also called critical mass
* @param Vi_plasmid Plasmid critical mass (experimental)
* @param ViNoise Gaussian noise standard deviation for a population associated with initiation volume (Vi)
* @param VaNoise Gaussian noise standard deviation for a population associated with cellular growth
* @param chanceInit Probability term that the competent replication fork, or origin of replication opens (if > than random number from Gaussian distribution with mean 0 and standard deviation of 1)
* @param divNoise Gaussian noise standard deviation for a population associated with division assymetry (WT = 0.5, i.e. perfect distribution of mass between mother and daughter cell)
* @param divRatio Gaussian noise mean asymmetry between mother and daughter cell at division (associated with divNoise)
* @param partRatio Gaussian noise mean partition noise for the distribution of chromosomes between mother and dauther cell at division
* @param partNoise Gaussian noise standard deviation for the distribution of chromosomes between mother and daughter cell at division
* @param chromDeg Probability term that a chromosome experiences DNA damage (if > than random number from Gaussian distribution with mean 0 and standard deviation of 1)
* @param repForkDeg Probability term that given a chromosome that experiences DNA damage, that this damage leads to either whole chromosome degradation, or replicating strand degradation
* @param C1 One-phase exponential function replication time first term
* @param C2 One-phase exponential function replication time second term
* @param C3 One-phase exponential function replication time third term
* @param D1 One-phase exponential function segregation time first term
* @param D2 One-phase exponential function segregation time second term
* @param D3 One-phase exponential function segregation time third term
* @param modelInitialParams GSL model initial parameters 
* @param modelInitialSpecies GSL model initial species
* @param modelGeneLocations Relative position of the genes on the chromosome for the GSL model
* @param modelGeneParamsLocations Position of the parameter that dictates the copy number influence on the expression of the GSL model
* @param modelGeneLRPos Position of the genes on either the left hand (0) or right hand (1) side of the chromosome 
* @param dt time step in min-1
*
* @return Error handling integer
*/
int setModel(Model * model, 
		float tau, 
		float cNoise, 
		float dNoise, 
		float Vi, 
		float Vi_plasmid, 
		float ViNoise, 
		float VaNoise, 
		float chanceInit, 
		float divNoise, 
		float divRatio, 
		float partRatio, 
		float partNoise,
		float chromDeg, 
		float repForkDeg, 
		float C1,
		float C2,
		float C3,
		float D1,
		float D2,
		float D3,
		double * modelInitialParams,
		double * modelInitialSpecies,
		float * modelGeneLocations,
		int * modelGeneParamsLocations,
		int * modelGeneLRPos,
		float dt)
{
	//make for a tmp file to write to
	/*
	char fileName[17];
	sprintf(fileName, "phases_su_%d.csv", (int)modelGeneLocations[0]);
        f = fopen(fileName, "w");
        if(f==NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
	fprintf(f, "time,numCells,Ga,stdGa,phase1,phase2,phase3,meanGene1,stdGene1,meanGene2,stdGene2,meanGene3,stdGene3\n");
	*/

	model->dt = dt;
	model->t = 0.0;
	model->stop = 0;
	
	//TODO: kind of stupid, should be able to tell him how many cells i want to inoculate with random number
	model->cellPopulation->numCells = 1;
	model->cellPopulation->freeIndex = 1;
	model->cellPopulation->indexArray = 1;
	model->cellPopulation->tau = tau;
	model->cellPopulation->C1 = C1;
	model->cellPopulation->C2 = C2;
	model->cellPopulation->C3 = C3;
	model->cellPopulation->D1 = D1;
	model->cellPopulation->D2 = D2;
	model->cellPopulation->D3 = D3;
	model->cellPopulation->dNoise = dNoise;
	model->cellPopulation->cNoise = cNoise;
	model->cellPopulation->Vi = Vi;
	model->cellPopulation->Vi_plasmid = Vi_plasmid;
	model->cellPopulation->ViNoise = ViNoise;
	model->cellPopulation->VaNoise = VaNoise;
	model->cellPopulation->chanceInit = chanceInit;
	model->cellPopulation->divNoise = divNoise;
	model->cellPopulation->divRatio = divRatio;
	model->cellPopulation->partNoise = partNoise;
	model->cellPopulation->partRatio = partRatio;	
	model->cellPopulation->chromDeg = chromDeg;
	model->cellPopulation->repForkDeg = repForkDeg;	
	model->cellPopulation->numFrozenCells = 0;	
	memcpy(model->cellPopulation->modelInitialParams, modelInitialParams, sizeof(double)*8); //TODO: don't make this to be hardcoded
	memcpy(model->cellPopulation->modelInitialSpecies, modelInitialSpecies, sizeof(double)*15); //TODO: same
	memcpy(model->cellPopulation->modelGeneLocations, modelGeneLocations, sizeof(float)*3); //TODO: same
	memcpy(model->cellPopulation->modelGeneParamsLocations, modelGeneParamsLocations, sizeof(int)*3); //TODO: same 
	memcpy(model->cellPopulation->modelGeneLRPos, modelGeneLRPos, sizeof(int)*3); //TODO: same
	/*
	memcpy(model->cellPopulation->modelInitialParams, modelInitialParams, sizeof(double)*(NUM_MODELPARAMS+1));
	memcpy(model->cellPopulation->modelInitialSpecies, modelInitialSpecies, sizeof(double)*(NUM_MODELSPECIES+1));
	memcpy(model->cellPopulation->modelGeneLocations, modelGeneLocations, sizeof(float)*(NUM_MODELGENES+1));
	memcpy(model->cellPopulation->modelGeneParamsLocations, modelGeneParamsLocations, sizeof(int)*(NUM_MODELGENES+1));
	memcpy(model->cellPopulation->modelGeneLRPos, modelGeneLRPos, sizeof(int)*(NUM_MODELGENES+1));
	*/
	return 0;
}

/**
* @brief Initialise all cells in cellArray and initiate the first cell
*
* Initiated all the cells in cellArray and flag them as being dead. Loop through the first set of cells and initiate the first cells, numCells, in the array (set as 1 at the moment). All cells are intitiated to have an age of 0.
* TODO: instead of directly accessing the cell.h function, pass through the cellPopulation.h
*
* @param model Model object
* @return Error handling integer
*/
int inoculateModel(Model * model)
{
	int i;
	for(i=0; i<model->cellPopulation->maxCells; i++)
	{
		constructCell(model->cellPopulation->cellArray, i);
	}

	for(i=0; i<model->cellPopulation->numCells; i++)
	{
		//if the input C time is input in minutes
		//TODO: add error handling from the initialiseCell
		if(model->cellPopulation->C2==-1.0)
		{
			initialiseCell(model->cellPopulation->cellArray,
				    i,
				    model->cellPopulation->tau,
				    model->cellPopulation->C1,
				    6646.0*model->cellPopulation->C1/4639221.0, //This is based on the size of the ColE1 in base pairs against the size of the chromosome in base pairs for a bacterial chromosome
				    model->cellPopulation->cNoise,
				    model->cellPopulation->D1,
				    model->cellPopulation->dNoise,
				    model->cellPopulation->Vi,
				    model->cellPopulation->Vi_plasmid,
				    model->cellPopulation->ViNoise,
				    model->cellPopulation->Vi/2.0, //TODO: why is this /2 ????
				    model->cellPopulation->VaNoise,
				    model->cellPopulation->modelInitialSpecies,
				    model->cellPopulation->modelInitialParams,
				    0.0);
		}
		//if the input C time is in its functional form
		else
		{
			initialiseCell(model->cellPopulation->cellArray,
				    i,
				    model->cellPopulation->tau,
				    model->cellPopulation->C1*(1.0+(model->cellPopulation->C2*exp(-model->cellPopulation->C3/(model->cellPopulation->tau/60.0)))),
				    6646.0*(model->cellPopulation->C1*(1.0+(model->cellPopulation->C2*exp(-model->cellPopulation->C3/(model->cellPopulation->tau/60.0)))))/4639221.0, //This is based on the size of the ColE1 in base pairs against the size of the chromosome in base pairs for a bacterial chromosome
				    model->cellPopulation->cNoise,
				    model->cellPopulation->D1*(1.0+(model->cellPopulation->D2*exp(-model->cellPopulation->D3/(model->cellPopulation->tau/60.0)))),
				    model->cellPopulation->dNoise,
				    model->cellPopulation->Vi,
				    model->cellPopulation->Vi_plasmid,
				    model->cellPopulation->ViNoise,
				    model->cellPopulation->Vi/2.0,
				    model->cellPopulation->VaNoise,
				    model->cellPopulation->modelInitialSpecies,
				    model->cellPopulation->modelInitialParams,
				    0.0);
		}
	}
	return 0;
}

/**
* @brief Clean the model for the purpose of freeing the memory correctly
*
* Clean the model. Call constructCell() on each cell in the array to free them, and free the cellArray.
*
* @param model Model object
* @return Error handling integer
*/
int cleanModel(Model * model)
{
	//fclose(f); close the close
	/*
	for(int i=0; i<model->cellPopulation->maxCells; i++)
	{
		gsl_odeiv2_driver_free(model->cellPopulation->cellArray[i].driver);
		constructCell(model->cellPopulation->cellArray, i);
	}	
	free(model->cellPopulation->cellArray);
	*/
	model->cellPopulation->cellArray = NULL;

	model->cellPopulation->maxCells = 0;	
	model->cellPopulation->numCells = 0;	
	model->cellPopulation->indexArray = 0;
	model->cellPopulation->freeIndex = 0;	
	model->cellPopulation->Vi = 0.0;	
	model->cellPopulation->Vi_plasmid = 0.0;	
	model->cellPopulation->cNoise = 0.0;	
	model->cellPopulation->dNoise = 0.0;	
	model->cellPopulation->ViNoise = 0.0;	
	model->cellPopulation->VaNoise = 0.0;	
	model->cellPopulation->chanceInit = 0.0;	
	model->cellPopulation->divNoise = 0.0;	
	model->cellPopulation->divRatio = 0.0;
	model->cellPopulation->partRatio = 0.0;
	model->cellPopulation->partNoise = 0.0;
	model->cellPopulation->chromDeg = 0.0;
	model->cellPopulation->repForkDeg = 0.0;	
	model->cellPopulation->numFrozenCells = 0;	
	model->cellPopulation->numAnucleateCells = 0;

	model->cellPopulation->C1 = 0.0;
	model->cellPopulation->C2 = 0.0;
	model->cellPopulation->C3 = 0.0;
	model->cellPopulation->D1 = 0.0;
	model->cellPopulation->D2 = 0.0;
	model->cellPopulation->D3 = 0.0;

	free(model->cellPopulation->totalVolumes);
        model->cellPopulation->lenTotalV = 0;

	//ODE model parameters
	model->cellPopulation->numModelParams = 0;
        model->cellPopulation->numModelSpecies = 0;
	model->cellPopulation->numModelGenes = 0;

	free(model->cellPopulation);
	model->cellPopulation = NULL;

        model->dt = 0.0;
        model->t = 0.0;
        model->stop = 0;
	
	free(model);
	model = NULL;
	return 0;
}

/**
* @brief Randomly delete live cells to restrict the population
*
* Clean all the cells that have been frozen (Warning: this could skew the results).
* Thereafter randomly delete live cells until the numCells==targetNumCells.  
*
* @param model Model object
* @param targetNumCells Number of cells to be left
* @return Error handling integer
*/
int randomRestrictNumCells(Model * model, int targetNumCells)
{
        int i, stop;
	int startNumCells = getRealNumCells(model);

	//first clean the cells that are isFrozen	
	//WARNING: This could skew the results
	if(model->cellPopulation->numFrozenCells!=0)
	{
		for(i=0; i<model->cellPopulation->maxCells; i++)
		{
			if(model->cellPopulation->cellArray[i].isFrozen==true)
			{
				constructCell(model->cellPopulation->cellArray, i);
				startNumCells--;
			}
		}
		model->cellPopulation->numFrozenCells = 0;
	}

	//randomly select cells and if they are not dead kill them until there are the desired amount
	while(startNumCells>targetNumCells)
	{
		float rnd = getRandomZeroOne();
		int rndIndex = (int)round(rnd*(model->cellPopulation->maxCells));
		if(0<=rndIndex && rndIndex<(model->cellPopulation->maxCells-1))
		{
			if(model->cellPopulation->cellArray[rndIndex].isDead==false)
			{
				constructCell(model->cellPopulation->cellArray, rndIndex);
				startNumCells--;
			}
		}
	}

	stop = organiseCellArray(model);
	if(stop==1)
	{
		printf("WARNING (randomRestrictNumCells): organiseCellArray has failed\n");
		return 1;
	}
	
	//model->cellPopulation->numCells = getRealNumCells(model);
	model->cellPopulation->numCells = startNumCells;
	return 0;
}

/**
* @brief Recalculate the volume addition of the model
*
* Given that we input the total volume calculated from OD, we need to normalise it to the volume of the simulated population
*
* @param model Model object
* @param totalVolumes Array of the total volume of the population, recalculated or not
* @param lenTotalV Length of the totalVolumes array
* @param volPos Position of the simulation in the totalVolumes array
* @return Error handling integer
**/
int recalculateVolumeAddition(Model * model, double * totalVolumes, int lenTotalV, int volPos)
{
	double sumVolume = getTotalVolume(model);
	double toDiv = totalVolumes[volPos]/sumVolume;	
	for(int i=0; i<lenTotalV; i++)
	{
		totalVolumes[i] = totalVolumes[i]/toDiv;
	}
}

//###############################################################################
//############################## Drug Treatment #################################
//###############################################################################

/**
* @brief Randomly block the segregation of cells in the population
*
* To emulate the drug treatment of cephalexin and rifampicin, we block the segregation (D), and the initiation of new origin of replication by setting these two to unreachable values. This is done stochastically (random Gaussian) to emulate that the drugs do not reach cells uniformally.
*
* @param model Model object
* @param drugNoise Probability based on Gaussian equation that the cell experiences the drug treatment
* @return Error handling integer
*/
int blockUnrep(Model * model, float drugNoise)
{
	int i;
	for(i=0; i<model->cellPopulation->maxCells; i++)
	{
		if(model->cellPopulation->cellArray[i].isDead==false && normalDistRandn(0.0,1.0)>drugNoise)
		{
			model->cellPopulation->cellArray[i].D = -1.0;
			model->cellPopulation->cellArray[i].Vi = 9999999999999.0;
			model->cellPopulation->cellArray[i].Vi_plasmid = 9999999999999.0;
		}
	}
	return 0;
}

/**
* @brief Check if any of the cells in the population is actively replicating their chromosomes
*
* Loop through all the cells and all the chromosomes and check if the replication forks are open or not. If any one of them are open return 1, if none are return 0. 
*
* @param model Model object
* @return Integer to see if a replication fork is open (1) or not (0)
*/
int checkIfAnyRep(Model * model)
{
	int i, y;
	for(i=0; i<model->cellPopulation->maxCells; i++)
	{
		for(y=0; y<model->cellPopulation->cellArray[i].numChrom; y++)
		{
			if(model->cellPopulation->cellArray[i].chromArray[y].replicationTimers[0][0][0]>0.0 &&
			model->cellPopulation->cellArray[i].chromArray[y].replicationTimers[0][0][1]>0.0)
			{
				return 1;		
			}
		}
	}		
	return 0;
}

/**
* @brief Run the drug treatment
*
* Loop through all the cells and induce the drug treatment, emulating cephalexin and rifampicin drugs (blockUnrep()), until all the cells in the simulation do not have a single replication fork open.
*
* @param model Model object
* @param maximalExecTime Timer that determines the maximal time permissive
* @param targetCellCount Restrict the number of cells max
* @param drugNoise Probability parameter that a cell experiences the drug treatment
* @return Error handling integer
*/
int runDrugTreatment(Model * model, float maximalExecTime, int targetCellCount, float drugNoise)
{
	model->cellPopulation->Vi = 9999999999999.0;
	model->cellPopulation->Vi_plasmid = 9999999999999.0;
	//model->cellPopulation->D = -1.0;
	clock_t begin = clock();
	float sumDT = 0.0;
	while(checkIfAnyRep(model)==1 && model->stop==0)
	{
		blockUnrep(model, drugNoise);
		model->stop = growCells(model->cellPopulation, model->dt, -2.0, 1.0);
                if(model->stop==2)
                {
                        randomRestrictNumCells(model, targetCellCount);
                        model->stop = 0;
                }
		if((float)(clock() - begin) / CLOCKS_PER_SEC >= maximalExecTime)
		{
			printf("ERROR (runExponential): Exceeded excecution time (500.0 seconds)\n");
			return 1;
		}
		sumDT += model->dt;
		//printf("\r##### dt: %.2f #####", sumDT);
	}		
	//printf("\n");
	//printf("numAnucleateCells: %d\n", model->cellPopulation->numAnucleateCells);
	return model->stop;
}

//#####################################################################################
//##################################### GETTERS #######################################
//#####################################################################################

//#################################### Number of Cells ################################

/**
* @brief Getter of the number of cells in the population
*
* Getter of number of cells from the cellPopulation structure
*
* @param model Model object
* @return Number of cells in the population
*/
int getNumCells(Model * model)
{
	return model->cellPopulation->numCells;
}

/**
* @brief Getter of the number of cells in the population
*
* Getter of number of cells by looping all the cellArray and counting the number of cells that are still alive
*
* @param model Model object
* @return Number of cells in the simulation
*/
int getRealNumCells(Model * model)
{
	int realNumCells = 0;
	int i;
	for(i=0; i<model->cellPopulation->maxCells; i++)
	{
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			realNumCells += 1;
		}
	}
	return realNumCells;
}

//##################################### Segregation Start ##########################

/**
* @brief Getter for the mean start of the segregation
*
* Getter for the mean start time for the segregation timer. Loops through all the cells and returns the mean segregation start time.
* 
* @param model Model object
* @return float Mean segregation time
*/
float getMeanSegStart(Model * model)
{
	int i;
	int count = 0;
	float segStart = 0.0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false && 
			model->cellPopulation->cellArray[i].init_a!=-1.0)
		{
			segStart += model->cellPopulation->cellArray[i].init_a;
			count += 1;
		}
        }
	return segStart/count;
}

/**
* @brief Getter for the standard deviation of the start of segregation
*
* Returns the standard deviation of the start of segregation by looping through all the cells in the simulation and calling getMeanSegStart()
* 
* @param model Model object
* @return stndard deviation for the start of the segregation
*/
float getStdSegStart(Model * model)
{
        float deviation = 0.0;
	float meanSegStart = getMeanSegStart(model);
        int numCells = 0;
        int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
                if(model->cellPopulation->cellArray[i].isDead==false &&
			model->cellPopulation->cellArray[i].init_a!=-1.0)
                {
			deviation += pow((model->cellPopulation->cellArray[i].init_a-meanSegStart),2);
                        numCells += 1;
                }
        }
        return sqrt(deviation/numCells);
}

//################################### doubling time #########################

/**
* @brief Getter for the previous doubling time of the cell
*
* Getter for the previous doubling time in the simulator. Loops through all the cells and returns their DNA content.
* 
* @param model Model object
* @return float Mean doubling time
*/
float getMeanDoublingTime(Model * model)
{
	int i;
	int count = 0;
	float doubTime = 0.0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false && 
			model->cellPopulation->cellArray[i].prev_sum_a!=-1.0)
		{
			doubTime += model->cellPopulation->cellArray[i].prev_sum_a;
			count += 1;
		}
        }
	return doubTime/count;
}

/**
* @brief Getter of oriC standard deviation of the population
*
* Returns the oriC standard deviation of the population in the simulation
* 
* @param model Model object
* @return stndard deviation of doubling time
*/
float getStdDoublingTime(Model * model)
{
	float meanDoublingTime = getMeanDoublingTime(model);
        float deviation = 0.0;
        int numCells = 0;
        int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
                if(model->cellPopulation->cellArray[i].isDead==false &&
			model->cellPopulation->cellArray[i].prev_sum_a!=-1.0)
                {
			deviation += pow((model->cellPopulation->cellArray[i].prev_sum_a-meanDoublingTime),2);
                        numCells += 1;
                }
        }
        return sqrt(deviation/numCells);
}

//############################################### number of oriC ######################

/**
* @brief Getter of the mean number of oriC in the simulator
*
* Getter of the mean number of oriC in the simulator. Loops through all the cells and returns their DNA content.
* 
* @param model Model object
* @return float Mean oriC
*/
float getMeanOriC(Model * model)
{
	int i, y;
	int count = 0;
	float meanOriC = 0.0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			for(y=0; y<MAX_CHROM; y++)
			{
				meanOriC += model->cellPopulation->cellArray[i].chromArray[y].oriC;
			}
			count += 1;
		}
        }
	return meanOriC/count;
}

/**
* @brief Getter of oriC standard deviation of the population
*
* Returns the oriC standard deviation of the population in the simulation. Calls getMeanOriC().
* 
* @param model Model object
* @return float standard deviation of oriC
*/
float getStdOriC(Model * model)
{
	float meanOriC = getMeanOriC(model);
	float deviation = 0.0;
	int numCells = 0;
	int i, y;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			for(y=0; y<MAX_CHROM; y++)
			{
				if(model->cellPopulation->cellArray[i].chromArray[y].oriC!=0)
				{
					deviation += pow((model->cellPopulation->cellArray[i].chromArray[y].oriC-meanOriC),2);
				}
			}
			numCells += 1;
		}
        }
	return sqrt(deviation/numCells);
}

//#################################### chromosome number ###############################

/**
* @brief Getter of the mean number of oriC in the simulator
*
* Getter of the mean number of oriC in the simulator. Loops through all the cells and returns their DNA content.
* 
* @param model Model object
* @return Mean oriC
*/
float getMeanChrom(Model * model)
{
	int i, y;
	int count = 0;
	float meanChrom = 0.0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			meanChrom += (float)model->cellPopulation->cellArray[i].numChrom;
			count += 1;
		}
        }
	return meanChrom/count;
}

/**
* @brief Getter of oriC standard deviation of the population
*
* Returns the oriC standard deviation of the population in the simulation. Calls getMeanChrom()
* 
* @param model Model object
* @return standard deviation of oriC
*/
float getStdChrom(Model * model)
{
	float meanChrom = getMeanChrom(model);
	float deviation = 0.0;
	int numCells = 0;
	int i, y;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			if(model->cellPopulation->cellArray[i].numChrom!=0)
			{
				deviation += pow((model->cellPopulation->cellArray[i].numChrom-meanChrom),2);
			}
			numCells += 1;
		}
        }
	return sqrt(deviation/numCells);
}

//################################### Ga ################################

/**
* @brief Getter of the DNA content of all the individual cells in the simulator
*
* Getter of the DNA content of the simulator. Loops through all the cells and returns their DNA content.
* 
* @param model Model object
* @param DNAContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int getDistGa(Model * model, float * Ga)
{
	int i, y;
	int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{	
			Ga[count] = model->cellPopulation->cellArray[i].Ga;
			count += 1;
		}
        }
	return 0;
}

/**
* @brief Getter of the Mean DNA content of the simulator
*
* Getter of the Mean DNA content of the simulator. Loops through all the cells and returns the mean DNA content of the population
* 
* @param model Model object
* @return Mean DNA content
*/
float getMeanGa(Model * model)
{
	int i;
	int count = 0;
	float meanDNAContent = 0.0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			meanDNAContent += model->cellPopulation->cellArray[i].Ga;
			count += 1;
		}
        }
	return meanDNAContent/count;
}

/**
* @brief Getter of Chromosome DNA content standard deviation of the population
*
* Returns the chromosome DNA content standard deviation of the population in the simulation. Calls the getMeanGa() function
* 
* @param model Model object
* @return Genetic content standard deviation
*/
float getStdGa(Model * model)
{
	float meanGa = getMeanGa(model);
	float deviation = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			deviation += pow((model->cellPopulation->cellArray[i].Ga-meanGa),2);
			numCells += 1;
		}
        }
	return sqrt(deviation/numCells);
}

//################################### Tau #################################

/**
* @brief Getter of the growth rate (tau) distribution of the population
*
* Getter of the growth rate distribution of the population
* 
* @param model Model object
* @param tauContent 1D array of size maxCells
* @return Error handling integer
*/
int getDistTau(Model * model, float * tauContent)
{
	int i, y;
	int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{	
			tauContent[count] = model->cellPopulation->cellArray[i].tau;
			count += 1;
		}
        }
	return 0;
}

/**
* @brief Getter of mean tau of the population
*
* Returns the average doubling time of the population in the simulation
* 
* @param model Model object
* @return Mean doubling time
*/
float getMeanTau(Model * model)
{
	float totalTau = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			totalTau += model->cellPopulation->cellArray[i].tau;
			numCells += 1;
		}
        }
	return totalTau/numCells;
}

/**
* @brief Getter of tau standard deviation of the population
*
* Returns the tau standard deviation of the population in the simulation. Calls getMeanTau()
* 
* @param model Model object
* @return Tau standard deviation
*/
float getStdTau(Model * model)
{
	float TauAv = getMeanTau(model);
	float deviation = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			deviation += pow((model->cellPopulation->cellArray[i].tau-TauAv),2);
			numCells += 1;
		}
        }
	return sqrt(deviation/numCells);
}

//################################# Pa ########################

/**
* @brief Getter of the plasmid DNA content of the simulator
*
* Getter of the plasmid DNA content of the simulator. Loops through all the cells and returns their plasmid DNA content.
* 
* @param model Model object
* @param Pa Empty 1D array of size maxCells
* @return Error handling integer
*/
int getPa(Model * model, float * Pa)
{
	int i;
	int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			Pa[count] = model->cellPopulation->cellArray[i].Pa;
			count += 1;
		}
        }
	return 0;
}

/**
* @brief Getter of mean plasmid DNA of the population
*
* Returns the average plasmid DNA of the population
* 
* @param model Model object
* @return Mean plasmid DNA content
*/
float getMeanPa(Model * model)
{
	float totalPlasmidDNA = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			totalPlasmidDNA += model->cellPopulation->cellArray[i].Pa;
			numCells += 1;
		}
        }
	return totalPlasmidDNA/numCells;
}

/**
* @brief Getter of Plasmid DNA content standard deviation of the population
*
* Returns the plasmid DNA content standard deviation of the population in the simulation. Call the getMeanPa() function
* 
* @param model Model object
* @return Standard deviation plasmid DNA content
*/
float getStdPa(Model * model)
{
	float meanPa = getMeanPa(model);
	float deviation = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			deviation += pow((model->cellPopulation->cellArray[i].Pa-meanPa),2);
			numCells += 1;
		}
        }
	return sqrt(deviation/numCells);
}

//################################## V #############################

/**
 * @brief Getter of sum of population volume
 *
 * Loops through all members of the population and appends their volume
 *  
 * @param model Model object
 * @return Total volume
 */
double getTotalVolume(Model * model)
{
	double totalVolumes = 0.0;
	int i;
	for(i=0; i<model->cellPopulation->maxCells; i++)
	{
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
		totalVolumes += model->cellPopulation->cellArray[i].Va;
		}
	}
	return totalVolumes;
}

/**
* @brief Getter of mean volume of the population
*
* Returns the average volume of the population in the simulation. Loops through all the cells in the population
* 
* @param model Model object
* @return Mean volume
*/
float getMeanVa(Model * model)
{
	float totalVolumes = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			totalVolumes += model->cellPopulation->cellArray[i].Va;
			numCells += 1;
		}
        }
	return totalVolumes/numCells;
}

/**
* @brief Getter of mean volume of the population
*
* Returns the average volume of the population in the simulation. Loops through all the cells in the simulation. Calls the getMeanVa() function
* 
* @param model Model object
* @return Mean volume
*/
float getStdVa(Model * model)
{
	float Vav = getMeanVa(model);
	float deviation = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			deviation += pow((model->cellPopulation->cellArray[i].Va-Vav),2);
			numCells += 1;
		}
        }
	return sqrt(deviation/numCells);
}

/**
* @brief Getter of mean volume of the population
*
* Returns the average volume of the population in the simulation
* 
* @param model Model object
* @return Mean volume
*/
float getTotalV(Model * model)
{
        float totalVolumes = 0.0;
        int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
                if(model->cellPopulation->cellArray[i].isDead==false)
                {
                        totalVolumes += model->cellPopulation->cellArray[i].Va;
                }
        }
        return totalVolumes;
}

/**
* @brief Getter of the volume distribution of the population
*
* Getter of the volume content of the simulator. Loops through all the cells and returns their volume.
* 
* @param model Model object
* @param Va Empty 1D array of size maxCells
* @return Error handling integer
*/
int getDistVa(Model * model, float * Va)
{
	int i, y;
	int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{	
			Va[count] = model->cellPopulation->cellArray[i].Va;
			count += 1;
		}
        }
	return 0;
}

//################################ age ###############################

/**
* @brief Getter of mean age of the population
*
* Returns the average age of the population in the simulation. Loops through all the cells in the population
* 
* @param model Model object
* @return Mean volume
*/
float getMeanAge(Model * model)
{
	float totalVolumes = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			totalVolumes += model->cellPopulation->cellArray[i].a;
			numCells += 1;
		}
        }
	return totalVolumes/numCells;
}

/**
* @brief Getter of standard deviation age of the population
*
* Returns the standard deviation age of the population in the simulation. Loops through all the cells in the simulation. Calls the getMeanA() function
* 
* @param model Model object
* @return Age standard deviation
*/
float getStdAge(Model * model)
{
	float Aav = getMeanAge(model);
	float deviation = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			deviation += pow((model->cellPopulation->cellArray[i].a-Aav),2);
			numCells += 1;
		}
        }
	return sqrt(deviation/numCells);
}

/**
* @brief Getter of the age distribution of the population
*
* Getter of the age of the simulator. Loops through all the cells and returns their age.
* 
* @param model Model object
* @param age Empty 1D array of size maxCells
* @return Error handling integer
*/
int getDistAge(Model * model, float * age)
{
	int i, y;
	int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{	
			age[count] = model->cellPopulation->cellArray[i].a;
			count += 1;
		}
        }
	return 0;
}

//######################################## Replication (C) time ########################
 
/**
* @brief Getter of mean population replication time (C)
*
* Loops through all members of the population and calculates mean replication time (C)
* 
* @param model Model object
* @return Mean C
*/
float getMeanC(Model * model)
{
	float meanC = 0.0;
	int i;
	int numCells = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			meanC += model->cellPopulation->cellArray[i].C;
			numCells += 1;
		}
        }
	return meanC/(float)numCells;	
}

/**
* @brief Getter of standard deviation replication time (C) of the population
*
* Returns the replication time (C) standard deviation of the population in the simulation. Loops through all the cells in the simulation. Calls the getMeanC() function
* 
* @param model Model object
* @return C standard deviation
*/
float getStdC(Model * model)
{
	float Cav = getMeanC(model);
	float deviation = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			deviation += pow((model->cellPopulation->cellArray[i].C-Cav),2);
			numCells += 1;
		}
        }
	return sqrt(deviation/numCells);
}

/**
* @brief Getter of the replication time (C) content of all the individual cells in the simulator
*
* Getter of the replication time (C) of the simulator. Loops through all the cells and returns their DNA content.
* 
* @param model Model object
* @param CContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int getDistC(Model * model, float * CContent)
{
	int i, y;
	int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{	
			CContent[count] = model->cellPopulation->cellArray[i].C;
			count += 1;
		}
        }
	return 0;
}

//######################################## Segregation (D) time ########################

/**
* @brief Getter of mean population segregation time (D)
*
* Loops through all members of the population and calculates mean segregation time (D)
* 
* @param model Model object
* @return Mean D
*/
float getMeanD(Model * model)
{
	float meanD = 0.0;
	int i;
	int numCells = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			meanD += model->cellPopulation->cellArray[i].D;
			numCells += 1;
		}
        }
	return meanD/(float)numCells;	
}

/**
* @brief Getter of standard deviation segregation time (D) of the population
*
* Returns the segregation time (D) standard deviation of the population in the simulation. Loops through all the cells in the simulation. Calls the getMeanD() function
* 
* @param model Model object
* @return D standard deviation
*/
float getStdD(Model * model)
{
	float Dav = getMeanD(model);
	float deviation = 0.0;
	int numCells = 0;
	int i;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			deviation += pow((model->cellPopulation->cellArray[i].D-Dav),2);
			numCells += 1;
		}
        }
	return sqrt(deviation/numCells);
}



/**
* @brief Getter of mean population segregation time (D)
*
* Loops through all members of the population and calculates mean segregation time (D)
* 
* @param model Model object
* @return Mean D
*/
float meanPopD(Model * model)
{
	float meanD = 0.0;
	int i;
	int numCells = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{
			meanD += model->cellPopulation->cellArray[i].D;
			numCells += 1;
		}
        }
	return meanD/(float)numCells;	
}

/**
* @brief Getter of the D time content of all the individual cells in the simulator
*
* Getter of the D time of the simulator. Loops through all the cells and returns their DNA content.
* 
* @param model Model object
* @param DContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int getDistD(Model * model, float * DContent)
{
	int i, y;
	int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
		if(model->cellPopulation->cellArray[i].isDead==false)
		{	
			DContent[count] = model->cellPopulation->cellArray[i].D;
			count += 1;
		}
        }
	return 0;
}

//############################## Prev parameters #####################

/** 
 * @brief Getter for the population distribution previous complete duplication time
 *
 * Returns a 1D array of all the previous age of the cells in the population (from birth to division)
 *
 * @param model Model object
 * @param dist_a Empty 1D array of size maxCells 
 * @return Error handling integer
 */
int getDistPrev_a(Model * model, float * dist_a)
{       
        int i, y; 
        int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {       
                if(model->cellPopulation->cellArray[i].isDead==false)
                {       
                        dist_a[count] = model->cellPopulation->cellArray[i].prev_sum_a;
                        count += 1;
                }
        }
        return 0;
}

/** 
 * @brief Getter for the population distribution previous volume at birth
 *
 * Returns a 1D array of all the previous volume at birth of the cells in the population
 *
 * @param model Model object
 * @param dist_Vb Empty 1D array of size maxCells 
 * @return Error handling integer
 */
int getDistPrev_Vb(Model * model, float * dist_Vb)
{       
        int i, y; 
        int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {       
                if(model->cellPopulation->cellArray[i].isDead==false)
                {       
                        dist_Vb[count] = model->cellPopulation->cellArray[i].prev_newbornVol;
                        count += 1;
                }
        }
        return 0;
}

/** 
 * @brief Getter for the population distribution previous volume at division
 *
 * Returns a 1D array of all the previous volume at division of the cells in the population 
 *
 * @param model Model object
 * @param dist_Vb Empty 1D array of size maxCells 
 * @return Error handling integer
 */
int getDistPrev_Vd(Model * model, float * dist_Vd)
{
        int i, y;
        int count = 0;
        for(i=0; i<model->cellPopulation->maxCells; i++)
        {
                if(model->cellPopulation->cellArray[i].isDead==false)
                {
                        dist_Vd[count] = model->cellPopulation->cellArray[i].prev_divisionVol;
                        count += 1;
                }
        }
        return 0;
}

//############################# Model ############################

/**
* @brief Getter of the mean phase of the population for the implementation of the repressilator in GSL
*
* Getter for the mean phases of the population when implementing the repressilator ODE model with the GSL simulator
*
* @param model Model Object
* @return Array float of the three genes percentage
**/
float * getMeanPhase(Model * model)
{
	static float ret[3];
	float phase1 = 0.0;
	float phase2 = 0.0;
	float phase3 = 0.0;
	int numCells = 0;
	for(int i=0; i<model->cellPopulation->maxCells; i++)
        {
                if(model->cellPopulation->cellArray[i].isDead==false)
                {
			double A = model->cellPopulation->cellArray[i].modelSpecies[0];
			double B = model->cellPopulation->cellArray[i].modelSpecies[5];
			double C = model->cellPopulation->cellArray[i].modelSpecies[10];
			//case where A > B and A > C
			if(A>=B && A>=C)
			{
				phase1 += 1.0;
			}
			else if(B>=A && B>=C)
			{
				phase2 += 1.0;
			}
			else if(C>=A && C>=B)
			{
				phase3 += 1.0;
			}
			numCells += 1;	
		}
	}
	ret[0] = phase1/(float)numCells;
	ret[1] = phase2/(float)numCells;
	ret[2] = phase3/(float)numCells;
	return ret;
}

/**
* @brief Getter of the phase of a single cell for the implementation of the repressilator in GSL
*
* Getter for the phase of a single cell when implementing the repressilator ODE model with the GSL simulator
*
* @param model Model Object
* @return Array float of the three genes
**/
float * getSinglePhase(Model * model, int cellNum)
{
	static float ret[3];
	double A = model->cellPopulation->cellArray[cellNum].modelSpecies[0];
	double B = model->cellPopulation->cellArray[cellNum].modelSpecies[5];
	double C = model->cellPopulation->cellArray[cellNum].modelSpecies[10];
	float phase1 = 0.0;
	float phase2 = 0.0;
	float phase3 = 0.0;
	if(A>=B && A>=C)
	{
		phase1 += 1.0;
	}
	else if(B>=A && B>=C)
	{
		phase2 += 1.0;
	}
	else if(C>=A && C>=B)
	{
		phase3 += 1.0;
	}
	ret[0] = phase1;
	ret[1] = phase2;
	ret[2] = phase3;
	return ret;
}

/**
* @brief Getter for the mean concentration of a population 
* 
* Returns the mean concentration of one of the species of the model within the population
*
* @param model Model object
* @param speciesNum Location in the array of the species of interest
* @return Mean Species
**/
//TODO: The true concentration is actually related to volume of the cell. Not sure how to implement this at this point in time (--> remove the dilution of the molecule from the degredation term and apply it to the simulation)
double getMeanModelSpecies(Model * model, int speciesNum)
{
	double meanSpecies = 0.0;
	int numCells = 0;
	//check that the range on the speciesNum is correct
	//if(speciesNum>=NUM_MODELSPECIES)
	if(speciesNum>=15) //TODO: this is temporary, implement NUM_MODELSPECIES
	{
		printf("WARNING (getMeanModelSpecies): sepciesNum out of range\n");
		return -1.0;
	}
	else
	{
		for(int i=0; i<model->cellPopulation->maxCells; i++)
		{
			if(model->cellPopulation->cellArray[i].isDead==false)
			{
				meanSpecies += model->cellPopulation->cellArray[i].modelSpecies[speciesNum];
				numCells += 1;
			}
		}
	}
	//printf("meanSpecies[%d]: %f/%f\n", speciesNum, meanSpecies, (double)numCells);
	return meanSpecies/(double)numCells;
}

/**
*@ brief Getter for a single concentration of a cell in the population 
* 
* Returns the concentration of one of the species of a single cell in the population
*
* @param model Model object
* @param speciesNum Location in the array of the species of interest
* @param speciesNum Location in the cellArray of the cell of interest
* @return Species amount
**/
double getSingleModelSpecies(Model * model, int cellNum, int speciesNum)
{
	//check that the species location is valid
	//if(speciesNum>=NUM_MODELSPECIES)
	if(speciesNum>=15)
	{
		printf("WARNING(getSingleModelSpecies): speciesNum is out of range\n");
		return -1.0;
	}
	//check that the cell of interest is alive
	if(model->cellPopulation->cellArray[cellNum].isDead==true)
	{
		printf("WARNING(getSingleModelSpecies): cellArray[%d].isDead The cell is dead\n", cellNum);
		return -1.0;
	}
	return model->cellPopulation->cellArray[cellNum].modelSpecies[speciesNum];
}

//####################################### Exposed Genes #################################

/**
* @brief Get all individual cells chromosomal genes that are exposed based on the position in the terminal
*
* Getter of all the cells exposed genes based on their location along the chromosome. Given a gene at position X, where 0.0 is located at the oriC and 1.0 is located at the terC, calculate the gene copy number. Practically this involves looping through all cells and all their chromosomes (in which case already the copy number is +1), and determining if the replication forks are past the input position (and if that is the case, once again that is +1). Because pair of replication fork opening and progression is stochastic, either left or right part of the replication fork must be specified where LR==0 is left hand side of the chromosome and LR==1 is the right hand side and LR==2 is both. 
* 
* @param model Model object
* @param numExposed Number of chromosomes that are exposed (must be the same structure as cell.h replicationTimers. i.e. [7][64][2])
* @param percLoc Location of the gene of interest (must be between >0.0 and <1.0)
* @param LR Because pair of replication fork progression is stochastic, the user may choose the side of the replication fork where the gene of interest is loca
ted. LR location on replicative side of the gene on the chromosome. 0--> Left, 1 --> Right, 2--> both.
* @return Mean number of genes 
*/

//TODO: check this mess underneath
//return the number of genes exposed based on the percentage location on the chromosome
//numExposed -> 2d array: 1rst order is the number of cells, 2nd nuber of max chrom, values are the normalised location (0.0->1.0) of the gene of interest
//NOTE: the array must be initialised to 0
//percLoc -> 1/2 the chrom == 100%. if the gene is on the left hand side of the chromosome 
//int LR location on replicative side of the gene on the chromosome. 0--> Left, 1 --> Right, 2--> both
int getNumExposedGenes(Model * model, int ** numExposed, float percLoc, int LR)
{
        int index;
        int i;
        int y;
        int z;
	int count = 0;
	//int geneCount[MAX_CHROM][getRealNumCells(model)] = {0};
        for(index=0; index<model->cellPopulation->maxCells; index++)
        {
                if(model->cellPopulation->cellArray[index].isDead==false)
                {
			// ############### DEBUG ###########
			/*
			printf("NumChrom: %d\n", model->cellPopulation->cellArray[index].numChrom);
			printf("Ga: %f\n", model->cellPopulation->cellArray[index].Ga);
			printf("C: %f\n", model->cellPopulation->cellArray[index].C);
			for(i=0; i<model->cellPopulation->cellArray[index].numChrom; i++)
			{
				for(y=0; y<log2(model->cellPopulation->cellArray[index].chromArray[i].potentialOriC); y++)
				{
					printf("[");
					for(z=0; z<pow(2,y); z++)
					{
						printf("[");
						for(j=0; j<2; j++)
						{
							printf("%f (%f),", model->cellPopulation->cellArray[index].chromArray[i].replicationTimers[y][z][j], model->cellPopulation->cellArray[index].chromArray[i].replicationTimers[y][z][j]/model->cellPopulation->cellArray[index].C);
						}               
						printf("]");
					}
					printf("]\n");
				}
				printf("]\n");
			}
			*/
			// #####################################
			int tmpPerCell = 0;
                        for(i=0; i<model->cellPopulation->cellArray[index].numChrom; i++)
                        {
                                if(LR==2)
                                {
                                        numExposed[count][i] += 2;
                                }
                                else
                                {
                                        numExposed[count][i] += 1;
                                }
				//for each potential pair or timers
				for(y=0; y<log2(model->cellPopulation->cellArray[index].chromArray[i].potentialOriC); y++)
				{
					for(z=0; z<pow(2,y); z++)
					{
						//0.5 of the
						if(LR==0 || LR==2)
						{
							if(model->cellPopulation->cellArray[index].chromArray[i].replicationTimers[y][z][0]>=percLoc*100.0)
							{
								numExposed[count][i] += 1;                      
							}
						}
						if(LR==1 || LR==2)
						{
							if(model->cellPopulation->cellArray[index].chromArray[i].replicationTimers[y][z][1]>=percLoc*100.0)
							{
								//printf("ADDED\n");
								numExposed[count][i] += 1;                      
							}

						}
					}
				}
                        }
                        count += 1;
                }
        }
	return 0;
}

/**
 * @brief Mean gene copy number of a chromosomal gene based on its relative position
 *
 * This function returns the mean copy number of a gene based on its relative position on the chromosome (where oriC==1.0 and terC==0.0). Because the two replication forks are stochastic in nature, the LR parameter determines if we are considering the left or right hand side of a chromosome. 
 *
 * @param model Model object
 * @param percLoc Location of the gene on the chromosome
 * @param LR Either left (0) or right (1) of the pair of replication forks
 */
float getMeanGeneCount(Model * model, float percLoc, int LR)
{
        int index;
        int i;
        int y;
        int z;
        int count = 0;
	int tmpSampleCount = 0;
	//int geneCount[MAX_CHROM][getRealNumCells(model)] = {0};
        for(index=0; index<model->cellPopulation->maxCells; index++)
        {
                if(model->cellPopulation->cellArray[index].isDead==false)
                {
			int tmpPerCell = 0;
                        for(i=0; i<model->cellPopulation->cellArray[index].numChrom; i++)
                        {
                                if(LR==2)
                                {
                                        tmpSampleCount += 2;
                                }
                                else
                                {
                                        tmpSampleCount += 1;
                                }
				//for each potential pair or timers
				for(y=0; y<log2(model->cellPopulation->cellArray[index].chromArray[i].potentialOriC); y++)
				{
					for(z=0; z<pow(2,y); z++)
					{
						//0.5 of the
						if(LR==0 || LR==2)
						{
							if(model->cellPopulation->cellArray[index].chromArray[i].replicationTimers[y][z][0]>=percLoc*100.0)
							{
								tmpSampleCount += 1;                      
							}
						}
						if(LR==1 || LR==2)
						{
							if(model->cellPopulation->cellArray[index].chromArray[i].replicationTimers[y][z][1]>=percLoc*100.0)
							{
								tmpSampleCount += 1;                      
							}

						}
					}
				}
                        }
                        count += 1;
        		//int a;
        		//scanf("%d", &a);
                }
		//printf("totalSampleExposed: %d\n", tmpSampleCount);
        }
	//printf("##############################\n");
	//printf("tmpSampleCount: %d\n", tmpSampleCount);
	//printf("numCells: %d\n", count);
	//printf("meanCount: %f\n", (float)tmpSampleCount/(float)count);
	//printf("##############################\n");
	return (float)tmpSampleCount/(float)count;
}

/**
 * @brief Standard deviation of gene copy number of a chromosomal gene based on its relative position
 *
 * This function returns the standard deviation for the copy number of a gene based on its relative position on the chromosome (where oriC==1.0 and terC==0.0). Because the two replication forks are stochastic in nature, the LR parameter determines if we are considering the left or right hand side of a chromosome. 
 *
 * @param model Model object
 * @param percLoc Location of the gene on the chromosome
 * @param LR Either left (0) or right (1) of the pair of replication forks
 */
float getStdGeneCount(Model * model, float percLoc, int LR)
{
	float meanGeneCount = getMeanGeneCount(model, percLoc, LR);
        int index;
        int i;
        int y;
        int z;
        int count = 0;
	float deviation = 0.0;
	int tmpGeneCount = 0;

        for(index=0; index<model->cellPopulation->maxCells; index++)
        {
                if(model->cellPopulation->cellArray[index].isDead==false)
                {
                        for(i=0; i<model->cellPopulation->cellArray[index].numChrom; i++)
                        {
                                if(LR==2)
                                {
                                        tmpGeneCount += 2;
                                }
                                else
                                {
                                        tmpGeneCount += 1;
                                }
				//for each potential pair or timers
				for(y=0; y<log2(model->cellPopulation->cellArray[index].chromArray[i].potentialOriC); y++)
				{
					for(z=0; z<pow(2,y); z++)
					{
						//0.5 of the
						if(LR==0 || LR==2)
						{
							if(model->cellPopulation->cellArray[index].chromArray[i].replicationTimers[y][z][0]>=percLoc*100.0)
							{
								tmpGeneCount += 1;                      
							}
						}
						if(LR==1 || LR==2)
						{
							if(model->cellPopulation->cellArray[index].chromArray[i].replicationTimers[y][z][1]>=percLoc*100.0)
							{
								tmpGeneCount += 1;                      
							}

						}
					}
				}
                        }	
                        deviation += pow(((float)tmpGeneCount-meanGeneCount),2);
			tmpGeneCount = 0;
                        count += 1;
                }
        }
        return sqrt(deviation/count);
}



//###############################################################################
//################################## GROWTH #####################################
//###############################################################################

/**
* @brief Run model assuming exponential growth until it reaches targetCellCount with functional forms of C and D 
*
* Run the model, growing the cells assuming Malthusian growth theory. This contains a timer that terminates the function if it exceeds a pre-determined time (500.0 seconds). The input C and D are functional forms, and are calculated from the mu[0].
*
* @param model Model object
* @param maximalExecTime Maximal execution time of the function. Avoids infinite loops
* @param targetCellCount Restrict the number of cells max. It is assumed that each step implies a step of dt. Thus 
* @param totalVolumes Array of the total population volume of the population, extracted from the OD
* @param lenTotalV Length of the totalVolumes array
*
* @return Boolean type parameter if the drug treatment is succesfull or not
*/
//Injection reporting the gene copy numbers every 5 minutes
int runInjection(Model * model, 
			float maximalExecTime, 
			int restrictCellNumber,
                	double * totalVolumes,
                	int lenTotalV)
{	
        float * meanPhases;
	/*
        FILE *f = fopen("phases_inj_01.csv", "w");
        if(f==NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
	*/
	/*
	float meanGene = 0.0;
	model->t = 0.0;
        FILE *f = fopen("copyNumber.csv", "w");
        if(f==NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
	*/
	/**
        FILE *f = fopen("population_270417.csv", "w");
        if(f==NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
        FILE *f1 = fopen("cell0_270417.csv", "w");
        if(f1==NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
	**/
	//fprintf(f, "time,numCells,indVolumeAdd,Ga,stdGa,oriC,stdOriC,chromNum,stdChromNum,meanV,stdV,totalV,divVol,stdDivVol\n");
	//fprintf(f1, "time,Ga,oriC,chromNum,tau,V,divVol\n");
	//fprintf(f, "time,numCells,indVolumeAdd,Ga,stdGa,oriC,stdOriC,chromNum,stdChromNum,meanV,stdV,totalV\n");
	//fprintf(f1, "time,Ga,oriC,chromNum,tau,V,divVol\n");
	//fprintf(f, "time,cellConc,Ga,stdGa,meanGeneCount_1_1,stdGeneCount_1_1,meanGeneCount_1_2,stdGeneCount_1_2,meanGeneCount_1_3,stdGeneCount_1_3,meanGeneCount_2_1,stdGeneCount_2_1,meanGeneCount_2_2,stdGeneCount_2_2,meanGeneCount_2_3,stdGeneCount_2_3,meanGeneCount_3_1,stdGeneCount_3_1,meanGeneCount_3_2,stdGeneCount_3_2,meanGeneCount_3_3,stdGeneCount_3_3\n");
	//fprintf(f, "time,numCells,Ga,stdGa,phase1,phase2,phase3,meanGene1,stdGene1,meanGene2,stdGene2,meanGene3,stdGene3\n");
	/*
	meanPhases = getMeanPhase(model);
	fprintf(f, "%f,", model->t);
	fprintf(f, "%d,", model->cellPopulation->numCells);
	fprintf(f, "%f,", (double)getMeanGa(model));
	fprintf(f, "%f,", (double)getStdGa(model, getMeanGa(model)));
	fprintf(f, "%f,", (double)meanPhases[0]);
	fprintf(f, "%f,", (double)meanPhases[1]);
	fprintf(f, "%f,", (double)meanPhases[2]);
	float meanGene1 = getMeanGeneCount(model, model->cellPopulation->modelGeneLocations[0]/100.0, 1);
	float meanGene2 = getMeanGeneCount(model, model->cellPopulation->modelGeneLocations[1]/100.0, 1);
	float meanGene3 = getMeanGeneCount(model, model->cellPopulation->modelGeneLocations[2]/100.0, 1);
	fprintf(f, "%f,", (double)meanGene1);
	fprintf(f, "%f,", (double)getStdGeneCount(model, model->cellPopulation->modelGeneLocations[0]/100.0, 1, meanGene1));
	fprintf(f, "%f,", (double)meanGene2);
	fprintf(f, "%f,", (double)getStdGeneCount(model, model->cellPopulation->modelGeneLocations[1]/100.0, 1, meanGene2));
	fprintf(f, "%f,", (double)meanGene3);
	fprintf(f, "%f", (double)getStdGeneCount(model, model->cellPopulation->modelGeneLocations[2]/100.0, 1, meanGene3));
	fprintf(f, "\n");
	*/

	int count = 0;
	recalculateVolumeAddition(model, totalVolumes, lenTotalV, count);
	clock_t begin = clock();
	while(model->stop==0 && (count+1)<lenTotalV)
	{
		model->t += model->dt;
		float indVolumeAdd = (totalVolumes[count+1]-totalVolumes[count])/model->cellPopulation->numCells;
		if(indVolumeAdd<0.0)
		{
			indVolumeAdd = 0.0;
		}	
		model->stop = growCells(model->cellPopulation, 
					model->dt, 
					indVolumeAdd,
					//(totalVolumes[count+1]-totalVolumes[count])/model->cellPopulation->numCells,
					0);
		if(model->stop==2)
		{
			printf("WARNING: Restrcting the number of cells\n");
			randomRestrictNumCells(model, restrictCellNumber);
			recalculateVolumeAddition(model, totalVolumes, lenTotalV, count);
			model->stop = 0;
		}
		if(model->stop==1)
		{
			printf("ERROR: model->stop==1 (not sure why this is)\n");
			break;
		}
		count += 1;
		
		if((float)(clock()-begin)/CLOCKS_PER_SEC>=maximalExecTime)
		{
			printf("ERROR (runInjection): Exceeded excecution time (500.0 seconds)\n");
			return 1;
		}
		//printf("\r##### numCells: %d #####", model->cellPopulation->numCells);
		//printf("\r##### numCells: %d (%.2f) #####", model->cellPopulation->numCells, (float)model->t);
		/*
		//PRINT
                fprintf(f, "%f,", model->t);
                fprintf(f, "%d,", model->cellPopulation->numCells);
		fprintf(f, "%f,", indVolumeAdd);
                fprintf(f, "%f,", (double)getMeanGa(model));
                fprintf(f, "%f,", (double)getStdGa(model, getMeanGa(model)));
                fprintf(f, "%f,", (double)getMeanOriC(model));
                fprintf(f, "%f,", (double)getStdOriC(model, getMeanOriC(model)));
                fprintf(f, "%f,", (double)getMeanChrom(model));
                fprintf(f, "%f,", (double)getStdChrom(model, getMeanChrom(model)));
                fprintf(f, "%f,", (double)getMeanVa(model));
                fprintf(f, "%f,", (double)getStdVa(model, getMeanVa(model)));
                fprintf(f, "%f", (double)getTotalV(model));
                //fprintf(f, "%f,", (double)getMeanDivVol(model));
                //fprintf(f, "%f", (double)getStdDivVol(model, getMeanDivVol(model)));
                fprintf(f, "\n");
                //fprintf(f1, "[");
                if(model->cellPopulation->cellArray[0].isDead==false)
                {
                        fprintf(f1, "%f,", model->t);
                        fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].Ga);
                        fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].totalPotentialOriC);
                        fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].numChrom);
                        fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].tau);
                        fprintf(f1, "%f", (double)model->cellPopulation->cellArray[0].Va);
                        //fprintf(f1, "%f", (double)model->cellPopulation->cellArray[0].divVol);
                        fprintf(f1, "\n");
                }
                else
                {
                        printf("WARNING: Cell 0 is dead\n");
                }
		*/
		
		//Print every 10 min
		/*
		if(fmod(roundf(100*model->t)/100,5.0)==0)
		{
			//printf("##### time: %f (%f) --> %f ##### \n", model->t, roundf(100*model->t)/100, fmod(roundf(100*model->t)/100, 10.0));
			fprintf(f, "%f,", model->t);
                	fprintf(f, "%d,", model->cellPopulation->numCells);
                	fprintf(f, "%f,", (double)getMeanGa(model));
                	fprintf(f, "%f,", (double)getStdGa(model, getMeanGa(model)));

			meanGene = getMeanGeneCount(model, 0.1, 0);
			fprintf(f, "%f,", meanGene);
			fprintf(f, "%f,", getStdGeneCount(model, 0.1, 0, meanGene));
			meanGene = getMeanGeneCount(model, 0.5, 0);
			fprintf(f, "%f,", meanGene);
			fprintf(f, "%f,", getStdGeneCount(model, 0.5, 0, meanGene));
			meanGene = getMeanGeneCount(model, 0.9, 0);
			fprintf(f, "%f,", meanGene);
			fprintf(f, "%f,", getStdGeneCount(model, 0.9, 0, meanGene));

			meanGene = getMeanGeneCount(model, 0.1, 0);
			fprintf(f, "%f,", meanGene);
			fprintf(f, "%f,", getStdGeneCount(model, 0.1, 0, meanGene));
			meanGene = getMeanGeneCount(model, 0.3, 0);
			fprintf(f, "%f,", meanGene);
			fprintf(f, "%f,", getStdGeneCount(model, 0.3, 0, meanGene));
			meanGene = getMeanGeneCount(model, 0.4, 0);
			fprintf(f, "%f,", meanGene);
			fprintf(f, "%f,", getStdGeneCount(model, 0.4, 0, meanGene));

			meanGene = getMeanGeneCount(model, 0.7, 0);
			fprintf(f, "%f,", meanGene);
			fprintf(f, "%f,", getStdGeneCount(model, 0.7, 0, meanGene));
			meanGene = getMeanGeneCount(model, 0.8, 0);
			fprintf(f, "%f,", meanGene);
			fprintf(f, "%f,", getStdGeneCount(model, 0.8, 0, meanGene));
			meanGene = getMeanGeneCount(model, 0.9, 0);
			fprintf(f, "%f,", meanGene);
			fprintf(f, "%f", getStdGeneCount(model, 0.9, 0, meanGene));

			fprintf(f, "\n");
		}
		*/
		/*
		if(fmod(roundf(100*model->t)/100,1.0)==0)
		{
			meanPhases = getMeanPhase(model);
			fprintf(f, "%f,", model->t);
			fprintf(f, "%d,", model->cellPopulation->numCells);
			fprintf(f, "%f,", (double)getMeanGa(model));
			fprintf(f, "%f,", (double)getStdGa(model, getMeanGa(model)));
			fprintf(f, "%f,", (double)meanPhases[0]);
			fprintf(f, "%f,", (double)meanPhases[1]);
			fprintf(f, "%f,", (double)meanPhases[2]);
			float meanGene1 = getMeanGeneCount(model, model->cellPopulation->modelGeneLocations[0]/100.0, 1);
			float meanGene2 = getMeanGeneCount(model, model->cellPopulation->modelGeneLocations[1]/100.0, 1);
			float meanGene3 = getMeanGeneCount(model, model->cellPopulation->modelGeneLocations[2]/100.0, 1);
			fprintf(f, "%f,", (double)meanGene1);
			fprintf(f, "%f,", (double)getStdGeneCount(model, model->cellPopulation->modelGeneLocations[0]/100.0, 1, meanGene1));
			fprintf(f, "%f,", (double)meanGene2);
			fprintf(f, "%f,", (double)getStdGeneCount(model, model->cellPopulation->modelGeneLocations[1]/100.0, 1, meanGene2));
			fprintf(f, "%f,", (double)meanGene3);
			fprintf(f, "%f", (double)getStdGeneCount(model, model->cellPopulation->modelGeneLocations[2]/100.0, 1, meanGene3));
                        fprintf(f, "\n");

			//printf("%f,", model->t);
			//printf("%d,", model->cellPopulation->numCells);
			//printf("%f,", (double)getMeanGa(model));
			//printf("%f,", (double)getStdGa(model, getMeanGa(model)));
			//printf("%f,", (double)meanPhases[0]);
			//printf("%f,", (double)meanPhases[1]);
			//printf("%f", (double)meanPhases[2]);
                        //printf("\n");
		}
		*/
	}
	//fclose(f);
	//printf("\n");
	//printf("Clock: %f\n", (float)(clock() - begin) / CLOCKS_PER_SEC);
	//printf("numAnucleateCells: %d\n", model->cellPopulation->numAnucleateCells);
	return model->stop;
}

/**
 * @brief Run the model assuming exponential growth
 *
 * Run the model assuming Malthusian growth starting from a single cell to the targetCellCount. Maximal execution time is made to avoid infinite loop. 
 *
 * @param model Model object
 * @param maximalExecTime Maximal execution time of the simulation
 * @param targetCellCount Target cell number for the simulation
 * @return Error handling integer
 */
int runExpo(Model * model, float maximalExecTime, int targetCellCount)
{	
	/*
        float * meanPhases;
        FILE *f1 = fopen("cell0.csv", "w");
        if(f1==NULL)
        {
                printf("Error opening file!\n");
                exit(1);
        }
	*/
	//fprintf(f, "time,numCells,Ga,stdGa,oriC,stdOriC,chromNum,stdChromNum,meanV,stdV,totalV,divVol,stdDivVol\n");
	//fprintf(f1, "time,Ga,chromNum,oriC,V,p1,p2,p3\n");
	//fprintf(f, "time,numCells,Ga,stdGa,oriC,stdOriC,chromNum,stdChromNum,meanV,stdV,totalV\n");
	//fprintf(f1, "time,Ga,oriC,chromNum,tau,V\n");
	//fprintf(f, "time,numCells,Ga,stdGa,phase1,phase2,phase3,meanGene1,stdGene1,meanGene2,stdGene2,meanGene3,stdGene3\n");
	//printf("time,numCells,Ga,stdGa,phase1,phase2,phase3\n");

	int count = 0;
	clock_t begin = clock();
	while(model->stop==0 && model->cellPopulation->numCells<targetCellCount)
	//while(model->stop==0 && model->t<1000.0)
	{
		model->t += model->dt;
		model->stop = growCells(model->cellPopulation, 
					model->dt, 
					-1.0,
					0);

		if(model->stop==2)
		{
			printf("WARNING: Restrcting the number of cells\n");
			randomRestrictNumCells(model, 2000);
			model->stop = 0;
		}
		if(model->stop==1)
		{
			printf("ERROR: model->stop==1 (not sure why this is)\n");
			break;
		}
		count += 1;
		
		if((float)(clock()-begin)/CLOCKS_PER_SEC>=maximalExecTime)
		{
			printf("ERROR (runInjection): Exceeded excecution time (500.0 seconds)\n");
			return 1;
		}
		//printf("\r##### numCells: %d #####", model->cellPopulation->numCells);
		//printf("\r##### simTime: %.2f #####", (float)model->t);
		//printf("\r##### numCells: %d (%.2f) #####", model->cellPopulation->numCells, (float)model->t);
		/*
		//PRINT
                fprintf(f, "%f,", model->t);
                fprintf(f, "%d,", model->cellPopulation->numCells);
                fprintf(f, "%f,", (double)getMeanGa(model));
                fprintf(f, "%f,", (double)getStdGa(model, getMeanGa(model)));
                fprintf(f, "%f,", (double)getMeanOriC(model));
                fprintf(f, "%f,", (double)getStdOriC(model, getMeanOriC(model)));
                fprintf(f, "%f,", (double)getMeanChrom(model));
                fprintf(f, "%f,", (double)getStdChrom(model, getMeanChrom(model)));
                fprintf(f, "%f,", (double)getMeanVa(model));
                fprintf(f, "%f,", (double)getStdVa(model, getMeanVa(model)));
                fprintf(f, "%f", (double)getTotalV(model));
                //fprintf(f, "%f,", (double)getMeanDivVol(model));
                //fprintf(f, "%f", (double)getStdDivVol(model, getMeanDivVol(model)));
                fprintf(f, "\n");
                //fprintf(f1, "[");
		*/
		/*
		if(model->cellPopulation->cellArray[0].isDead==false)
		{
			if(fmod(roundf(100*model->t)/100,1.0)==0)
			{
				fprintf(f1, "%f,", model->t);
				fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].Ga);
				fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].numChrom);
				fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].totalPotentialOriC);
				fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].Va);
				fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].modelSpecies[1]);
				fprintf(f1, "%f,", (double)model->cellPopulation->cellArray[0].modelSpecies[6]);
				fprintf(f1, "%f", (double)model->cellPopulation->cellArray[0].modelSpecies[11]);
				//fprintf(f1, "%f", (double)model->cellPopulation->cellArray[0].divVol);
				fprintf(f1, "\n");
			}
		}
                else
                {
                        printf("WARNING: Cell 0 is dead\n");
                }
		*/
		/*
		*/
		/*
		if(fmod(roundf(100*model->t)/100,5.0)==0)
		{
			meanPhases = getMeanPhase(model);
			fprintf(f, "%f,", model->t);
			fprintf(f, "%d,", model->cellPopulation->numCells);
			fprintf(f, "%f,", (double)getMeanGa(model));
			fprintf(f, "%f,", (double)getStdGa(model, getMeanGa(model)));
			fprintf(f, "%f,", (double)meanPhases[0]);
			fprintf(f, "%f,", (double)meanPhases[1]);
			fprintf(f, "%f,", (double)meanPhases[2]);
			float meanGene1 = getMeanGeneCount(model, model->cellPopulation->modelGeneLocations[0]/100.0, 1);
			float meanGene2 = getMeanGeneCount(model, model->cellPopulation->modelGeneLocations[1]/100.0, 1);
			float meanGene3 = getMeanGeneCount(model, model->cellPopulation->modelGeneLocations[2]/100.0, 1);
			fprintf(f, "%f,", (double)meanGene1);
			fprintf(f, "%f,", (double)getStdGeneCount(model, model->cellPopulation->modelGeneLocations[0]/100.0, 1, meanGene1));
			fprintf(f, "%f,", (double)meanGene2);
			fprintf(f, "%f,", (double)getStdGeneCount(model, model->cellPopulation->modelGeneLocations[1]/100.0, 1, meanGene2));
			fprintf(f, "%f,", (double)meanGene3);
			fprintf(f, "%f", (double)getStdGeneCount(model, model->cellPopulation->modelGeneLocations[2]/100.0, 1, meanGene3));
                        fprintf(f, "\n");

			//printf("%f,", model->t);
			//printf("%d,", model->cellPopulation->numCells);
			//printf("%f,", (double)getMeanGa(model));
			//printf("%f,", (double)getStdGa(model, getMeanGa(model)));
			//printf("%f,", (double)meanPhases[0]);
			//printf("%f,", (double)meanPhases[1]);
			//printf("%f", (double)meanPhases[2]);
                        //printf("\n");
		}
		*/
	}
	//fclose(f);
	//printf("\n");
	//printf("Clock: %f\n", (float)(clock() - begin) / CLOCKS_PER_SEC);
	//printf("numAnucleateCells: %d\n", model->cellPopulation->numAnucleateCells);
	return model->stop;
}

/**
 * @brief Run the model where with this function a single time step is taken
 *
 * Every time this function is called it takes a single time step as defined with dt. Designed to be externally controlled (such as the python wrapper)
 *
 * @param model Model object
 * @return Error handling integer
 */
int oneTimeStep(Model * model)
{	
	model->t += model->dt;
	model->stop = growCells(model->cellPopulation, model->dt, -1.0, 0);	
	if(model->stop==1)
	{
		printf("ERROR: model->stop==1 (not sure why this is)\n");
		return 1;
	}
	return 0;
}

//##########################################################################################################
//############################################# MAIN #######################################################
//##########################################################################################################
/**
 * @brief Main function
 *
 * Main function 
 *
 * @return Error handling integer
 */
int main()
{
	float Vi = 0.9;
	float Vi_plasmid = 9999999.0;
	float ViNoise = 5.0;
	float VaNoise = 10.0;
	float cNoise = 5.0;
	float dNoise = 5.0;
	float chanceInit = 1.75;
	float divNoise = 10.0;
	float divRatio = 0.5;
	float partRatio = 0.5;
	float partNoise = -1.0;
	float chromDeg = 1000.0;
	float repForkDeg = 1000.0;
	float dt = 0.01;
	int maxCells = 20000;
	//float maxExecTime = 9999999.0;
	float maxExecTime = 1500.0;
	//int targetCellCount = 5000;
	int targetCellCount = 19500;
	int restrictCellNumber = 5000;
	float drugNoise = 3.3;
	float tau = 60.0;
	float maximalExecTime = 500.0;

	/*
	float C1 = 43.0;
	float C2 = -1.0;
	float C3 = -1.0;
	float D1 = 23.0;
	float D2 = -1.0;
	float D3 = -1.0;
	*/
	float C1 = 43.2;
	float C2 = 4.86;
	float C3 = 4.5;
	float D1 = 23.0;
	float D2 = 4.86;
	float D3 = 4.5;

	float mu = log(2.0)/tau;
	
	//float db_h = (log(2)/(tau/60.0));
	//float trad_C = 43.2*(1.0+(4.86*exp(-4.5/db_h)));
	//float trad_D = 23.0*(1.0+(4.86*exp(-4.5/db_h)));
	
	//float C = 40.0;
	//float D = 20.0;


	//###### Disable the model ######
	//double modelInitialParams[8] = {6.0, 1.0, 5.0, 100.0, 5.0, 100.0, 10.0, 20.0};
	//Note that here we start with 10.0 copy numbers of each promoter
	//double modelInitialSpecies[15] = {10.0, 0.0, 2.0, 8.0, 0.0, 
	//				  10.0, 0.0, 8.0, 2.0, 0.0, 
	//				  10.0, 0.0, 5.0, 5.0, 0.0};
	//float modelGeneLocations[3] = {0.1*100.0, 0.5*100.0, 0.9*100.0};
	
	//int modelGeneParamsLocations[3] = {2, 7, 12};
	//int modelGeneLRPos[3] = {1, 1, 1};

	double modelInitialParams[8] = {6.0, 1.0, 5.0, 100.0, 5.0, 100.0, 10.0, 20.0};
	double modelInitialSpecies[15] = {10.0, 0.0, 2.0, 8.0, 0.0,
                                          10.0, 0.0, 8.0, 2.0, 0.0,
                                          10.0, 0.0, 5.0, 5.0, 0.0};
	//float modelGeneLocations[3] = {0.9*100.0, 0.9*100.0, 0.9*100.0};
	float modelGeneLocations[3] = {0.1*100.0, 0.5*100.0, 0.9*100.0};
	int modelGeneParamsLocations[3] = {2, 7, 12};
        int modelGeneLRPos[3] = {1, 1, 1};

	printf("InitModel\n");
	Model * model = initModel(maxCells);
	printf("SetModel\n");
	setModel(model,
                tau,
                cNoise,
                dNoise,
                Vi,
                Vi_plasmid,
                ViNoise,
                VaNoise,
                chanceInit,
                divNoise,
                divRatio,
                partRatio,
                partNoise,
                chromDeg,
                repForkDeg,
                C1,
                C2,
                C3,
                D1,
                D2,
                D3,
                modelInitialParams,
                modelInitialSpecies,
                modelGeneLocations,
                modelGeneParamsLocations,
                modelGeneLRPos,
                dt);
	printf("InoculateModel\n");
	inoculateModel(model);

	//printf("runInjection\n");
	/*	
	int lenTotalV = 100001;
	double od[100001] = {0.0};
	od[0] = 0.001;
	for(int i=0; i<100000; i++)
	{
		od[i+1] = od[i]*(1.0+(log(2.0)/tau)*dt);
	}
	double totalVolumes[100001] = {0.0};
	for(int i=0; i<100001; i++)
	{
		totalVolumes[i] = (3.6*od[i])*pow(10.0, 9.0);
	}

	runInjection(model,
                        maximalExecTime,
                        restrictCellNumber,
                        totalVolumes,
                        lenTotalV);
	*/
	
	printf("runExpo\n");
	runExpo(model, maxExecTime, targetCellCount);
	
	float oriCav = getMeanOriC(model);
	float oriCstd = getStdVa(model);
	float Vav = getMeanVa(model);
	float Vstd = getStdVa(model);
	float TauAv = getMeanTau(model);
	float TauStd = getStdTau(model);
	float GaAv = getMeanGa(model);
	float GaStd = getStdGa(model);
	float PaAv = getMeanPa(model);
	float PaStd = getStdPa(model);

	printf("V: %f (%f)\n", Vav, Vstd);
	printf("OriC: %f (%f)\n", oriCav, oriCstd);
	printf("Tau: %f (%f)\n", TauAv, TauStd);
	printf("Ga: %f (%f)\n", GaAv, GaStd);
	printf("Pa: %f (%f)\n", PaAv, PaStd);
	printf("#################################################\n");

	cleanModel(model);
	return 0;
}
