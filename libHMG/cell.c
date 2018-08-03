/*
 *   libHMG - Individual based model of a bacterial population
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
* @file cell.c
*
* Contains all the functions that are related to the growth of the cell and chromosome. Includes the replication of the cell, its division into mother and daughter cell. Handles the growth of the chromosomes by setting growing all the replication forks and duplicating them once they have finiched replicating themselves. This also includes the implementation of recA related degradation of full chromosomes and the replicating strand if the chromosme experiences DNA damage. Includes an implementation of the eclipse period; i.e, inhibits the start of a new replication cycle until a minimal distance is reached between two replicating chromosomes 
*
* @version 1.0
* @author Melchior du Lac
*
*/

//TODO: consider seperating the functions related to the chromsome and the growth into different files for clarity
//TODO: finish impementation of one and only once replication of chromosomes per cell cycle

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
/*
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
*/

#include "cell.h"
#include "utility.h"
#include "inputModel.h"
#include "libsbmlsim/libsbmlsim/libsbmlsim.h"

//TODO: make all of these parameters in the dynamic by allocating the calloc and not setting this to be static like this
//int MAX_CHROM = 32; //64;
//int LOG2_MAX_REP_TIMERS = 6; //7; //where there can be 64*2 (128) max replication forks on a single chromosome (overkill?)
//int MAX_PLASMID = 2000;
//int NUM_MODELPARAMS = 8;
//int NUM_MODELSPECIES = 15;
//int NUM_MODELGENES = 3;

//################################ PRIVATE FUNCTIONS #############################
/**
 * @brief Debug a single cell
 *
 * Private method used to print the parameters of a cell in the array
 *
 * @param cellArray cellPopulation object
 * @param index Location of the cell of interest in the array
 * @return void
 */
void cDebug(Cell * cellArray, int index)
{
        if(cellArray[index].isDead==false)
        {
                printf("############# Cell: %d #################\n", index);
                printf("cellArray[%d].tau -> %f\n", index, cellArray[index].tau);
                printf("cellArray[%d].Va -> %f\n", index, cellArray[index].Va);
                printf("cellArray[%d].a -> %f\n", index, cellArray[index].a);
                printf("cellArray[%d].prev_sum_a -> %f\n", index, cellArray[index].prev_sum_a);
                printf("cellArray[%d].init_a -> %f\n", index, cellArray[index].init_a);
                printf("cellArray[%d].Vi -> %f\n", index, cellArray[index].Vi);
                printf("cellArray[%d].Vi_plasmid -> %f\n", index, cellArray[index].Vi_plasmid);
                printf("cellArray[%d].numChrom -> %d\n", index, cellArray[index].numChrom);
                printf("cellArray[%d].numPlasmid -> %d\n", index, cellArray[index].numPlasmid);
                printf("cellArray[%d].totalPotentialOriC -> %d\n", index, cellArray[index].totalPotentialOriC);
                printf("cellArray[%d].totalPotentialOriC_plasmid -> %d\n", index, cellArray[index].totalPotentialOriC_plasmid);
                printf("cellArray[%d].totalRepForks -> %d\n", index, cellArray[index].totalRepForks);
                printf("cellArray[%d].injectionDeviation -> %f\n", index, cellArray[index].injectionDeviation);
                printf("cellArray[%d].Ga -> %f\n", index, cellArray[index].Ga);
                printf("cellArray[%d].Pa -> %f\n", index, cellArray[index].Pa);
                printf("cellArray[%d].segregationTimer -> %f\n", index, cellArray[index].segregationTimer);
                printf("cellArray[%d].segregationTimer2 -> %f\n", index, cellArray[index].segregationTimer2);
                //printf("cellArray[%d].divVol -> %f\n", index, cellArray[index].divVol);
                printf("cellArray[%d].C -> %f\n", index, cellArray[index].C);
                printf("cellArray[%d].C_plasmid -> %f\n", index, cellArray[index].C_plasmid);
                printf("cellArray[%d].D -> %f\n", index, cellArray[index].D);
                int i;
                int y;
                int z;
                int j;
		//CHROMOSOME
                for(i=0; i<MAX_CHROM; i++)
                {
			if(cellArray[index].chromArray[i].oriC!=0)
			{
				printf("cellArray[%d].chromArray[%d]\n", index, i);
				printf("\toriC -> %d\n", cellArray[index].chromArray[i].oriC);
				printf("\trepFroks -> %d\n", cellArray[index].chromArray[i].repForks);
				printf("\tpotentialOriC -> %d\n", cellArray[index].chromArray[i].potentialOriC);
				if(cellArray[index].chromArray[i].potentialOriC!=0)
				{
					printf("\tActive Log2(RepForks) -> %f\n", log2(cellArray[index].chromArray[i].potentialOriC));
				}
				printf("\treplicationTimers-> \n[");
				for(y=0; y<log2(cellArray[index].chromArray[i].potentialOriC); y++) //<-- usually 5
				//for(y=0; y<7; y++) //<-- usually 5
				{
					printf("[");
					for(z=0; z<pow(2,y); z++)
					{
						printf("[");
						for(j=0; j<2; j++)
						{
							//printf("%.2f (%.2f),",cellArray[index].chromArray[i].replicationTimers[y][z][j], cellArray[index].chromArray[i].Ctimes[y][z][j]);
							printf("%.2f,",cellArray[index].chromArray[i].replicationTimers[y][z][j]);
						}               
						printf("]");
					}
					printf("]\n");
				}
				printf("]\n");
			}
                }
		//PLASMID
		for(i=0; i<MAX_PLASMID; i++)
		{
			if(cellArray[index].plasmidArray[i].oriC!=0)
			{
				printf("cellArray[%d].plasmidArray[%d]\n", index, i);
				printf("\toriC -> %d\n", cellArray[index].plasmidArray[i].oriC);
				printf("\treplicationTimers -> [%.2f][%.2f]\n", 
					cellArray[index].plasmidArray[i].replicationTimers[0], 
					cellArray[index].plasmidArray[i].replicationTimers[1]);
			}
		}
		//MODEL
		/*
		printf("cellArray[%d].modelPrevSumGenes -> [", index);
		for(i=0; i<NUM_MODELGENES; i++)
		{
			printf("%d, ", cellArray[index].modelPrevSumGenes[i]);
		}
		printf("]\n");

		printf("cellArray[%d].modelSumGenes -> [", index);
		for(i=0; i<NUM_MODELGENES; i++)
		{
			printf("%d, ", cellArray[index].modelSumGenes[i]);
		}
		printf("]\n");

		printf("cellArray[%d].modelParams -> [", index);
		for(i=0; i<NUM_MODELPARAMS; i++)
		{
			printf("%f, ", cellArray[index].modelParams[i]);
		}
		printf("]\n");

		printf("cellArray[%d].modelSpecies -> [", index);
		for(i=0; i<NUM_MODELSPECIES; i++)
		{
			printf("%f, ", cellArray[index].modelSpecies[i]);
		}
		printf("]\n");
		*/
                //printf("\tAddress-> [");
                //for(i=0; i<cellArray[index].numChrom; i++)
                //{
                //      printf("%p,", cellArray[index].chromArray[i].replicationTimers);
                //}
                //printf("]\n");
        }
        else
        {
                printf("Cell is dead\n");
        }
        //int a;
        //scanf("%d", &a);
}

/**
* @brief Delete the chromosome from cell's chromArray
*
* By setting the oriC to 0 the cell is flagged as being non-existing. Loop through all the possible replication forks and set them to 0.0 to make sure the DNA is not counted.
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @param i Chromosome of interest in the chromArray
* @return Error handling integer
*/
int clearChromosome(Cell * cellArray, int index, int i)
{
	cellArray[index].chromArray[i].oriC = 0;
	cellArray[index].chromArray[i].potentialOriC = 0;
	cellArray[index].chromArray[i].repForks = 0;
	int y, z;
	for(y=0; y<LOG2_MAX_REP_TIMERS; y++)
	{
		for(z=0; z<pow(2,y); z++)
		{
			cellArray[index].chromArray[i].replicationTimers[y][z][0] = 0.0;
			cellArray[index].chromArray[i].replicationTimers[y][z][1] = 0.0;
		}
	}
	return 0;
}

/**
* @brief Delete plasmid from cell's chromArray
*
* By setting the oriC to 0 the cell is flagged as being non-existing. Loop through all the possible replication forks and set them to 0.0 to make sure the DNA is not counted.
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @param i Chromosome of interest in the chromArray
* @return Error handling integer
*/
int clearPlasmid(Cell * cellArray, int index, int i)
{
	cellArray[index].plasmidArray[i].oriC = 0;
	cellArray[index].plasmidArray[i].replicationTimers[0] = 0.0;
	cellArray[index].plasmidArray[i].replicationTimers[1] = 0.0;
	return 0;
}

/**
* @brief Reorganise the chromosome of chromArray
*
* Loop forward through all the chromosome in the chromArray and if an empty position is detected reverse loop through it and find the next chromosome and swap the poisitions
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @return Error handling integer
*/
//TODO: remove this function, and implement the queue
int reorganiseChromosomes(Cell * cellArray, int index)
{
	int i, y;
	//faster but dependent on the reliability of numChrom
	//for(i=0; i<cellArray[index].numChrom; i++)
	//slower but more accurate
	for(i=0; i<MAX_CHROM; i++)
	{
		if(cellArray[index].chromArray[i].oriC==0)
		{
			//for(y=MAX_CHROM; y>i; y--)
			for(y=MAX_CHROM-1; y>=0; y--)
			{
				//every position < i are assumed to be filled when this condition is met
				if(i>=y)
				{
					return 0;
				}
				if(cellArray[index].chromArray[y].oriC!=0)
				{
					memcpy(cellArray[index].chromArray[i].replicationTimers, 
					       cellArray[index].chromArray[y].replicationTimers, 
					       sizeof(float)*(LOG2_MAX_REP_TIMERS)*((int)pow(2, LOG2_MAX_REP_TIMERS-1))*2);
					cellArray[index].chromArray[i].oriC = cellArray[index].chromArray[y].oriC;
					cellArray[index].chromArray[i].repForks = cellArray[index].chromArray[y].repForks;
					cellArray[index].chromArray[i].potentialOriC = (int)cellArray[index].chromArray[y].oriC;
					clearChromosome(cellArray, index, y);	
					break;
				}
			}
		}
	}
	//should never really happen unless whole chromArray is filled
	return 1;
}

/**
* @brief Reorganise the plasmids of plasmidArray
*
* Loop forward through all the plasmids in the plasmidArray and if an empty position is detected reverse loop through it and find the next plasmid and swap the poisitions
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @return Error handling integer
*/
int reorganisePlasmids(Cell * cellArray, int index)
{
	int i, y;
	for(i=0; i<MAX_PLASMID; i++)
	{
		if(cellArray[index].plasmidArray[i].oriC==0)
		{
			for(y=MAX_PLASMID-1; y>=0; y--)
			{
				if(i>=y)
				{
					return 0;
				}
				if(cellArray[index].plasmidArray[y].oriC!=0)
				{
					cellArray[index].plasmidArray[i].replicationTimers[0] = cellArray[index].plasmidArray[y].replicationTimers[0]; 
					cellArray[index].plasmidArray[i].replicationTimers[1] = cellArray[index].plasmidArray[y].replicationTimers[1]; 
					cellArray[index].plasmidArray[i].oriC = cellArray[index].plasmidArray[y].oriC;
					clearPlasmid(cellArray, index, y);
					break;
				}
			}
		}
	}
}

//########################### SETTERS #########################
/**
* @brief Construct the cell object by setting eveything to 0
*
* Set/reset the cell to its default (dead) state. Includes resetting all of its attributes to 0 and looping through the chromArray to delete any replication forks that are open
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @return Error handling integer
*/
int constructCell(Cell * cellArray, int index)
{
	cellArray[index].isDead = true;
	cellArray[index].isNewlyRep = false;
	cellArray[index].isInitiated = false;
	cellArray[index].isRep = false;
	cellArray[index].tau = 0.0;
	cellArray[index].C = 0.0;
	cellArray[index].C_plasmid = 0.0;
	cellArray[index].D = 0.0;
	cellArray[index].Vi = -1.0;
	cellArray[index].Vi_plasmid = 0.0;
	cellArray[index].totalPotentialOriC = 0;
	cellArray[index].totalPotentialOriC_plasmid = 0;
	cellArray[index].totalRepForks = 0;
	cellArray[index].a = 0.0;
	cellArray[index].prev_sum_a = -1.0;
	cellArray[index].init_a = -1.0;
	cellArray[index].Va = 0.0;
	cellArray[index].Ga = 0.0;
	cellArray[index].Pa = 0.0;
	cellArray[index].injectionDeviation = 0.0;
	cellArray[index].prev_sum_a = 0.0;
	cellArray[index].prev_newbornVol= 0.0;
	cellArray[index].prev_divisionVol = 0.0;

	//Chromosomes
	cellArray[index].numChrom = 0;
	cellArray[index].isFrozen = false;
	int i;
	for(i=0; i<MAX_CHROM; i++)
	{
		//TODO: change the name to initialiseChromosome
		clearChromosome(cellArray, index, i);
	}
	for(i=0; i<MAX_PLASMID; i++)
	{
		//TODO: change the name to initialiseChromosome
		clearPlasmid(cellArray, index, i);
	}
	cellArray[index].segregationTimer = 0.0;
	cellArray[index].segregationTimer2 = 0.0;
	//cellArray[index].divVol = 0.0;
	//############## Model parameters ###########
	/*
	for(i=0; i<NUM_MODELSPECIES; i++)
	{
		cellArray[index].modelSpecies[i] = 0.0;
	}
	for(i=0; i<NUM_MODELPARAMS; i++)
	{
		cellArray[index].modelParams[i] = 0.0;
	}
	for(i=0; i<NUM_MODELGENES; i++)
	{
		cellArray[index].modelSumGenes[i] = 0;
	}
	for(i=0; i<NUM_MODELGENES; i++)
	{
		cellArray[index].modelPrevSumGenes[i] = 0;
	}
	*/

	//############### libsbmlsim #############	

	if(sbml_results!=NULL)
	{
		free_myResult(cellArray[index].sbml_results);
	}	
	cellArray[index].sbml_err_num = -1;
	cellArray[index].sbml_ato = 0.0;
	cellArray[index].sbml_rtol = 0.0;
	cellArray[index].sbml_facmax = 0.0;

	return 0;
}

/**
* @brief Intialise the cell containing a single chromosome
*
* From a dead cell, we construct the cell initialise the global population parameters (such as C) with noise from their corresponding Gaussian parameters (for example cNoise). Initialise the first chromosome in the chromArray of the cell to have a cell with a DNA chromosome equivalent of 1.0.
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @param tau Exponential doubling rate
* @param C Replication time
* @param C_plasmid Plasmid replication time
* @param cNoise Replication Gaussian standard deviation
* @param D Segregation time
* @param dNoise Segregation Gaussian standard deviation
* @param Vi Critical volume
* @param Vi_plasmid Plasmid critical volume
* @param ViNoise Critical volume Gaussian standard deviation
* @param Va Volume at time a
* @param VaNoise Volume at time a Gaussian standard deviation
* @param modelInitialSpecies GSL model species input
* @param modelInitialParams GSL model parameters input
* @param a Age of the cell (not used. Future application involves if this value >0.0 calculate the DNA content and segragation)
* @return Error handling integer
*/
int initialiseCell(Cell * cellArray,
		    int index,
		    float tau,
		    float C,
		    float C_plasmid,
		    float cNoise,
		    float D,
		    float dNoise,
		    float Vi,
		    float Vi_plasmid,
		    float ViNoise,
		    float Va,
		    float VaNoise,
		    double * modelInitialSpecies,
		    double * modelInitialParams,
		    float a)
{
	cellArray[index].tau = normalDistRandn(tau, (tau*VaNoise/100.0));
	cellArray[index].isDead = false;
	cellArray[index].isNewlyRep = false;
	cellArray[index].isInitiated = false;
	cellArray[index].isRep = false;
	cellArray[index].C = normalDistRandn(C, (C*cNoise/100.0));
	cellArray[index].C_plasmid = normalDistRandn(C_plasmid, (C_plasmid*cNoise/100.0)); //TODO: nNoise is assumed to be the same as the other
	cellArray[index].D = normalDistRandn(D, (D*dNoise/100.0));
	cellArray[index].Vi = normalDistRandn(Vi, (Vi*ViNoise/100.0));
	cellArray[index].Vi_plasmid = normalDistRandn(Vi_plasmid, (Vi_plasmid*ViNoise/100.0)); //TODO: plasmid noise for initiaition is assumed to be the same

	cellArray[index].prev_sum_a = 0.0;
	cellArray[index].prev_newbornVol= 0.0;
	cellArray[index].prev_divisionVol = 0.0;
	cellArray[index].a = a;
	cellArray[index].prev_sum_a = -1.0;
	cellArray[index].init_a = -1.0;

	cellArray[index].Va = Va;
	cellArray[index].Ga = 1.0;
	cellArray[index].Pa = 1.0;
	cellArray[index].segregationTimer = 0.0;
	cellArray[index].segregationTimer2 = 0.0;
	//cellArray[index].divVol = normalDistRandn(adderDivVol(tau), (adderDivVol(tau)*ViNoise)/100.0);
	cellArray[index].totalPotentialOriC = 1;
	cellArray[index].totalPotentialOriC_plasmid = 1;
	cellArray[index].totalRepForks = 0;
	cellArray[index].numChrom = 1;
	cellArray[index].numPlasmid = 1;
	//TODO: only if using the injection growth model
	cellArray[index].injectionDeviation = normalDistRandn(1.0, (VaNoise/100.0));

	cellArray[index].totalRepForks = 0;
	cellArray[index].chromArray[0].oriC = 1;
	cellArray[index].chromArray[0].potentialOriC = 1;	
	cellArray[index].chromArray[0].repForks = 0;	

	cellArray[index].plasmidArray[0].oriC = 1;

	cellArray[index].isFrozen = false;

	//################################# GSL #########################o
	
	//copy the initial parameters to cell
	/*
	memcpy(cellArray[index].modelSpecies, modelInitialSpecies, sizeof(double)*NUM_MODELSPECIES);
	memcpy(cellArray[index].modelParams, modelInitialParams, sizeof(double)*NUM_MODELPARAMS);
	*/
	//memcpy(cellArray[index].modelSpecies, modelInitialSpecies, sizeof(double)*(NUM_MODELSPECIES+1));
	//memcpy(cellArray[index].modelParams, modelInitialParams, sizeof(double)*(NUM_MODELPARAMS+1));
	//set the gene copy numbers from the model to 1.0
	/*
	for(int i=0; i<NUM_MODELGENES; i++)
	{
		cellArray[index].modelSumGenes[i] = 1.0;
		cellArray[index].modelPrevSumGenes[i] = 1.0;
	}
	*/

	// set the gsl ODE solver 
	/*
	cellArray[index].sys.function = model_bennett_repressilator; //TODO: pass this from model.c
	cellArray[index].sys.jacobian = NULL;
	cellArray[index].sys.dimension = NUM_MODELSPECIES;
	cellArray[index].sys.params = &cellArray[index].modelParams;
	
        cellArray[index].driver = gsl_odeiv2_driver_alloc_y_new(&cellArray[index].sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0); //TODO: make the last three inputs part of the model
	*/

	//######################## libsbml model #################
	
	  if (cellArray[index].d == NULL)
	{
	    //return create_myResult_with_errorCode(Unknown);
		printf("ERROR: The SBML document has not been loaded");
		return 1;
	}
	  err_num = SBMLDocument_getNumErrors(cellArray[index].d);
	  if (err_num > 0) 
	  {
	    const XMLError_t *err = (const XMLError_t *)SBMLDocument_getError(cellArray[index].d, 0);
	   if (XMLError_isError(err) || XMLError_isFatal(err)) {
	      XMLErrorCode_t errcode = XMLError_getErrorId(err);
	      switch (errcode) {
		case XMLFileUnreadable:
		  cellArray[index].sbml_results = create_myResult_with_errorCode(FileNotFound);
		  break;
		case XMLFileUnwritable:
		case XMLFileOperationError:
		case XMLNetworkAccessError:
		  cellArray[index].sbml_results = create_myResult_with_errorCode(SBMLOperationFailed);
		  break;
		case InternalXMLParserError:
		case UnrecognizedXMLParserCode:
		case XMLTranscoderError:
		  cellArray[index].sbml_results = create_myResult_with_errorCode(InternalParserError);
		  break;
		case XMLOutOfMemory:
		  cellArray[index].sbml_results = create_myResult_with_errorCode(OutOfMemory);
		  break;
		case XMLUnknownError:
		  cellArray[index].sbml_results = create_myResult_with_errorCode(Unknown);
		  break;
		default:
		  cellArray[index].sbml_results = create_myResult_with_errorCode(InvalidSBML);
		  break;
	      }
	      SBMLDocument_free(d);
	      return 1;
	    }
	  }
	  cellArray[index].m = SBMLDocument_getModel(d);

	return 0;
}

//############################## CELLULAR SUB-FUNCTIONS ############################
//##################################### PRIVATE ####################################

/**
* @brief Divide the cell once replication time (D) is reached
*
* Divide the cell into mother and daughter cell. This includes seperating the chromosome symmetrically, albeit with Gaussian noise, and same with volume distribution. At this point we update the population level parameters (example C, D etc...). Reorganise the chromosomes.
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @param newIndex Empty position in the cellArray for daughter cell
* @param tau Exponential doubling rate
* @param C Replication time
* @param cNoise Replication Gaussian standard deviation
* @param D Segregation time
* @param dNoise Segregation Gaussian standard deviation
* @param Vi Critical volume
* @param ViNoise Critical volume Gaussian standard deviation
* @param divRatio Division mean volume ratio distribution between mother and daughter cell (default=0.5)
* @param divNoise Division standard deviation volume ration distribution between mother and daughter cell 
* @param partRatio Chromosome partition mean ration between mother and daughter cell (default=0.5)
* @param partNoise Chromosome partition standard deviation between mother and daughter cell
* @param modelGeneParamsLocations Location of the parameters in the GSL parameter array that determines the gene number
* @param modelGeneLocations Location of the genes on the chromsome
* @param modelGeneLRPos Determine if the gene is on the left or right hand side of the oriC
* @param dt Time step
* @return Error handling integer
*/
int divideCell(Cell * cellArray,
		int index,
		int newIndex,
		float tau,
		float C,
		float C_plasmid,
		float cNoise,
		float D,
		float dNoise,
		float Vi,
		float Vi_plasmid,
		float ViNoise,
		float VaNoise,
		float divRatio,
		float divNoise,
		float partRatio,
		float partNoise,
		int * modelGeneParamsLocations,
		float * modelGeneLocations,
		int * modelGeneLRPos,
		float dt)
{
	//printf("################## Dividing Cell %d -> %d ##########################", index, newIndex);
	//cDebug(cellArray, index);
	//#####################################################


	int i, y, z;
	//Here we stochastically assign the daughter cell a new Volume by taking half from the mother cell
	//Note: the ratio is also used to transfer the model parameters
	float newVolume = cellArray[index].Va*normalDistRandn(divRatio,((divNoise*divRatio)/100.0));
	int motherTotalOriC = cellArray[index].totalPotentialOriC;
	int daughterTotalOriC = 0;

	//DAUGHTER CELL	-- initialise
	//Use the function instead of defining these manually. The initial parameters for the model in this case need to be overwritten with half of mother cell
	//create 2 tmp species arrays where the mother and dauther concentrations are halfed --> then delete them from the heap to avoid memory leackage
	double tmpMotherSpecies[15] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //FILL THIS with 0.0 when the size is determined
	double tmpDaughterSpecies[15] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	
	//the amount of each species concentration that is trasnferred to the daughter cell is determined by the volume transfer (i.e. newVolume)
	//float percVolTrans = newVolume*100.0/cellArray[index].Va;
	/*
	if(index==0)
	{
		printf("################## DIVIDING ####################\n");
		printf("cellArray[0].modelSpecies: [");
		for(i=0; i<NUM_MODELSPECIES; i++)
		{
			printf("%f,", cellArray[index].modelSpecies[i]);
		}
		printf("]\n");

		for(i=0; i<NUM_MODELSPECIES; i++)
		{
			printf("%f*(%f*100.0/%f)/100.0\n", cellArray[index].modelSpecies[i], newVolume, cellArray[index].Va);
		}
	}
	*/
	for(i=0; i<NUM_MODELSPECIES; i++)
	{
		tmpDaughterSpecies[i] = cellArray[index].modelSpecies[i]*(newVolume*100.0/cellArray[index].Va)/100.0;
		tmpMotherSpecies[i] = cellArray[index].modelSpecies[i]-tmpDaughterSpecies[i];
	}
	/*
	if(index==0)
	{
		printf("tmpDaughterSpecies: [");
		for(i=0; i<NUM_MODELSPECIES; i++)
		{
			printf("%f,", tmpDaughterSpecies[i]);
		}
		printf("]\n");
	}	
	*/

	//RESET THEM TO 0
	for(i=0; i<NUM_MODELGENES; i++)
	{
		cellArray[newIndex].modelSumGenes[i] = 0;
		cellArray[newIndex].modelPrevSumGenes[i] = 0;
	}
	//STUPID maybe pass the pointer then
	//float passVol = cellArray[index].Va-newVolume;
	//double * passModelParams = cellArray[index].modelParams;
	//float zero = 0.0;
	initialiseCell(cellArray, newIndex, tau, C, C_plasmid, cNoise, D, dNoise, Vi, Vi_plasmid, ViNoise, cellArray[index].Va-newVolume, VaNoise, tmpDaughterSpecies, cellArray[index].modelParams, 0.0); 
	//parameters to overwrite from cell initialisation
        cellArray[newIndex].Ga = 0.0;
        cellArray[newIndex].Pa = 0.0;
        cellArray[newIndex].totalPotentialOriC = 0;
        cellArray[newIndex].totalPotentialOriC_plasmid = 0;
        cellArray[newIndex].totalRepForks = 0;
        cellArray[newIndex].numChrom = 0;
        cellArray[newIndex].numPlasmid = 0;
        cellArray[newIndex].totalRepForks = 0;
        cellArray[newIndex].chromArray[0].oriC = 0;
        cellArray[newIndex].chromArray[0].potentialOriC = 0;
        cellArray[newIndex].chromArray[0].repForks = 0;
        cellArray[newIndex].plasmidArray[0].oriC = 0;
	cellArray[newIndex].prev_divisionVol = cellArray[index].Va;
	cellArray[newIndex].prev_sum_a = cellArray[index].a;
        cellArray[newIndex].prev_newbornVol = cellArray[newIndex].Va;
	
	//randomly select half the chromosomes to segregate to daughter
	float tmpGa = 0.0;

	//control the number of chromosomes that are to be passed to the daughter stochasticilly
	int totalChrom = cellArray[index].numChrom;
	int transferChromNum = 0;	

	//######## calculate the ratio of chromosomes that are passed to the daughter cell ####
	// partNoise can be manipulated to make this more noisy
	if(partNoise<=0.0)
	{
		transferChromNum = floor(totalChrom/2);
	}
	else
	{
		transferChromNum = round(normalDistRandn(partRatio, ((partRatio*partNoise)/100.0))*cellArray[index].numChrom);	
		//transferChromNum = round(normalDistRandn(partRatio, partNoise)*cellArray[index].numChrom);	
		//normalDistRandn(C_plasmid, (C_plasmid*cNoise/100.0));
		if(transferChromNum<=0)
		{
			transferChromNum = 0;
		}
		else if(transferChromNum>totalChrom)
		{
			transferChromNum = totalChrom;
		}	
	}	

	//printf("\t\tDEBUG (divideCell): Beginning loop\n");
	//################ TRANSFER HALF OF THE CHROMOSOMES ##################
	while(cellArray[newIndex].numChrom<transferChromNum)
	{
		//printf("\t\t\tDEBUG (divideCell): cellArray[newIndex].numChrom(%d;%d)<transferChromNum(%d)\n", cellArray[newIndex].numChrom, cellArray[index].numChrom, transferChromNum);
		//printf("cellArray[newIndex].numChrom(%d)<transferChromNum(%d)\n", cellArray[newIndex].numChrom, transferChromNum);
		//WARNING--> not reorganising the mother cell chromosome array when it comes out
		int rnd = rand()%totalChrom; //this looks wrong
		//int rnd = rand()%MAX_CHROM;

		//printf("\t\t\trnd: %d\n", rnd);
		if(cellArray[index].chromArray[rnd].oriC!=0)
		{
			memcpy(cellArray[newIndex].chromArray[cellArray[newIndex].numChrom].replicationTimers, 
			       cellArray[index].chromArray[rnd].replicationTimers, 
			       sizeof(float)*(LOG2_MAX_REP_TIMERS)*((int)pow(2,LOG2_MAX_REP_TIMERS-1))*2);
			cellArray[newIndex].chromArray[cellArray[newIndex].numChrom].oriC = cellArray[index].chromArray[rnd].oriC;
			cellArray[newIndex].chromArray[cellArray[newIndex].numChrom].potentialOriC = (int)cellArray[index].chromArray[rnd].oriC;
        		cellArray[newIndex].chromArray[cellArray[newIndex].numChrom].repForks = cellArray[index].chromArray[rnd].repForks;

			//need to do this to determine the Ga content that has been passed
			for(i=0; i<log2(cellArray[index].chromArray[rnd].potentialOriC); i++)
			{
				for(y=0; y<pow(2,i); y++)
				{
					//################ Recalculate the tmpGa ##################
					if(cellArray[index].chromArray[rnd].replicationTimers[i][y][0]!=0.0)
                                        {
						//tmpGa += (float)(cellArray[index].chromArray[i].replicationTimers[i][y][0]/100.0)/2.0;
						tmpGa += (float)(cellArray[index].chromArray[rnd].replicationTimers[i][y][0]/100.0)/2.0;
                                        }
                                        if(cellArray[index].chromArray[rnd].replicationTimers[i][y][1]!=0.0)
                                        {
						//tmpGa += (float)(cellArray[index].chromArray[i].replicationTimers[i][y][1]/100.0)/2.0;
						tmpGa += (float)(cellArray[index].chromArray[rnd].replicationTimers[i][y][1]/100.0)/2.0;
                                        }
					//################# Recalculate the different modelSumGenes  #################
					for(int gN=0; gN<NUM_MODELGENES; gN++)
					{
						if(modelGeneLRPos[gN]==0 &&
							cellArray[index].chromArray[rnd].replicationTimers[i][y][0]>modelGeneLocations[gN])
						{
							cellArray[newIndex].modelSumGenes[gN] += 1;
							cellArray[newIndex].modelPrevSumGenes[gN] += 1;
						}
						if(modelGeneLRPos[gN]==1 &&
							cellArray[index].chromArray[rnd].replicationTimers[i][y][1]>modelGeneLocations[gN])
						{
							cellArray[newIndex].modelSumGenes[gN] += 1;
							cellArray[newIndex].modelPrevSumGenes[gN] += 1;
						}
						//###################################################
					}
				}
			}

			cellArray[newIndex].numChrom += 1;
			daughterTotalOriC += cellArray[index].chromArray[rnd].oriC;
			cellArray[newIndex].totalPotentialOriC += cellArray[index].chromArray[rnd].potentialOriC;
			cellArray[newIndex].totalRepForks += cellArray[index].chromArray[rnd].repForks;	

			cellArray[index].numChrom -= 1;
			motherTotalOriC -= cellArray[index].chromArray[rnd].oriC;
			cellArray[index].totalPotentialOriC -= cellArray[index].chromArray[rnd].potentialOriC;
			cellArray[index].totalRepForks -= cellArray[index].chromArray[rnd].repForks;	

			clearChromosome(cellArray, index, rnd);
		}	
	}
	//printf("\t\tDEBUG (divideCell): Ending loop\n");
	
	//############# TRANSFER HALF OF THE PLASMIDS ####
	// unlike the chromosomes, colE1 does not have active paritioning systems, and thus is random
	float tmpPa = 0.0;
	int totalPlasmidCount = cellArray[index].numPlasmid;
	int transferPlasmidNum = floor(cellArray[index].numPlasmid/2);
	int newIndex_plasmidOriC = 0;
	int index_plasmidOriC = cellArray[index].numPlasmid;
	//while(cellArray[newIndex].numPlasmid<transferPlasmidNum)
	if(false)
	{
		printf("For the moment should not be in here\n");
		int rnd = rand()%totalPlasmidCount;
		if(cellArray[index].plasmidArray[rnd].oriC!=0)
		{	
			if(cellArray[index].plasmidArray[rnd].replicationTimers[0]!=0.0)
			{
				tmpPa += (float)(cellArray[index].plasmidArray[rnd].replicationTimers[0]/100.0)/2.0;		
			}
			if(cellArray[index].plasmidArray[rnd].replicationTimers[1]!=0.0)
			{
				tmpPa += (float)(cellArray[index].plasmidArray[rnd].replicationTimers[1]/100.0)/2.0;
			}
			
			cellArray[newIndex].plasmidArray[cellArray[newIndex].numPlasmid].replicationTimers[0] = cellArray[index].plasmidArray[rnd].replicationTimers[0];
			cellArray[newIndex].plasmidArray[cellArray[newIndex].numPlasmid].replicationTimers[1] = cellArray[index].plasmidArray[rnd].replicationTimers[1];
			cellArray[newIndex].plasmidArray[cellArray[newIndex].numPlasmid].oriC = cellArray[index].plasmidArray[rnd].oriC;	

			cellArray[newIndex].numPlasmid += 1;
			cellArray[index].numPlasmid -= 1;
			newIndex_plasmidOriC += cellArray[index].plasmidArray[rnd].oriC;	
			index_plasmidOriC -= cellArray[index].plasmidArray[rnd].oriC;	
		
			clearPlasmid(cellArray, index, rnd);
		}	
	}	

	//####### Recalculate the different copy number
	//TO overwite: 
	//modelSumGenes <- needs to be recalculated
	//modelParams needs to be updated after modelSumGenes is calculated

	cellArray[newIndex].totalPotentialOriC_plasmid = newIndex_plasmidOriC;
	cellArray[index].totalPotentialOriC_plasmid = index_plasmidOriC;

	//when we divide we do not want the potential oriC to be still there
	cellArray[newIndex].totalPotentialOriC = daughterTotalOriC;

	cellArray[newIndex].Ga = cellArray[newIndex].numChrom+tmpGa;
	cellArray[newIndex].Pa = cellArray[newIndex].numPlasmid+tmpPa;
	//TODO: check the accuracy of this stupid thing
	for(int gN=0; gN<NUM_MODELGENES; gN++)
	{
		//cellArray[newIndex].modelSumGenes[gN] += cellArray[newIndex].numChrom;
		//cellArray[newIndex].modelPrevSumGenes[gN] += cellArray[newIndex].numChrom;
		cellArray[index].modelPrevSumGenes[gN] = cellArray[index].modelSumGenes[gN]-cellArray[newIndex].modelSumGenes[gN];
		cellArray[index].modelSumGenes[gN] = cellArray[index].modelSumGenes[gN]-cellArray[newIndex].modelSumGenes[gN];
		//update the sumGenes to the modelParam
		//cellArray[newIndex].modelParams[modelGeneParamsLocations[gN]] = cellArray[newIndex].modelSumGenes[gN];
		//cellArray[index].modelParams[modelGeneParamsLocations[gN]] = cellArray[index].modelSumGenes[gN];
	}

	cellArray[newIndex].segregationTimer = cellArray[index].segregationTimer2;
	cellArray[newIndex].segregationTimer2 = 0.0; 

	//PARENT CELL
	cellArray[index].Ga = cellArray[index].Ga-cellArray[newIndex].Ga;
	cellArray[index].Pa = cellArray[index].Pa-cellArray[newIndex].Pa;
	cellArray[index].isDead = false;
	cellArray[index].isInitiated = false;
	cellArray[index].isRep = false;
	cellArray[index].prev_divisionVol = cellArray[index].Va;
	cellArray[index].Va = newVolume;	
	cellArray[index].prev_sum_a = cellArray[index].a;
	cellArray[index].a = 0.0;
	cellArray[index].C = normalDistRandn(C, (C*cNoise/100.0));
	cellArray[index].C_plasmid = normalDistRandn(C_plasmid, (C_plasmid*cNoise/100.0));
	cellArray[index].D = normalDistRandn(D, (D*dNoise/100.0));
	cellArray[index].Vi = normalDistRandn(Vi, (Vi*ViNoise/100.0));
	cellArray[index].Vi_plasmid = normalDistRandn(Vi_plasmid, (Vi_plasmid*ViNoise/100.0));
	cellArray[index].tau = normalDistRandn(tau, (tau*VaNoise/100.0));
        cellArray[index].injectionDeviation = normalDistRandn(1.0, (VaNoise/100.0));
	memcpy(cellArray[index].modelSpecies, tmpMotherSpecies, sizeof(double)*NUM_MODELSPECIES);
	//memcpy(cellArray[index].modelSpecies, tmpMotherSpecies, sizeof(double)*(NUM_MODELSPECIES+1));
        cellArray[index].prev_newbornVol = newVolume;

	//reorganise the parent cell chromosome 
	//loop through the numChrom and if empty stop found take the bottom chromosome and move it
	reorganiseChromosomes(cellArray, index);
	reorganisePlasmids(cellArray, index);

	cellArray[index].segregationTimer = cellArray[index].segregationTimer2;	
	cellArray[index].segregationTimer2 = 0.0;	
	//cellArray[index].divVol = adderDivVol(tau);

	//when we divide we do not want the potential oriC to be still there
	cellArray[index].totalPotentialOriC = motherTotalOriC;	
	/*
	printf("##############################################################\n");
	cDebug(cellArray, index);
	cDebug(cellArray, newIndex);
	printf("##############################################################\n");
	int a;
        scanf("%d", &a);
	*/
	/*
	if(index==0)
	{
		printf("cellArray[0].modelSpecies: [");
		for(i=0; i<NUM_MODELSPECIES; i++)
		{
			printf("%f,", cellArray[index].modelSpecies[i]);
		}
		printf("]\n");
	}	
	*/

	//cDebug(cellArray, index);
	//cDebug(cellArray, newIndex);
	
        //int a;
        //scanf("%d", &a);	

	//printf("################## Mother Cell ##########################");
	//cDebug(cellArray, index);
	//#####################################################
	//printf("################## Daughter Cell ##########################");
	//cDebug(cellArray, newIndex);
	return 0;
}

/**
* @brief Grow chromosome and divide the chromosome if the chromosome finishes replicating itself
*
* Grow the chromosome if the replication fork is opened. Check if the origin of replication should be open if potentialOriC>numOriC of the chromosome. This also stands for the replication forks, where we would open stochastically the next replication forks if potentialRepForks>numRepForks. If any of the chromosome have completed replication of its chromosomes, then the segragation timer of the chromoosme is initiated.  
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @param dt Time step
* @param chanceInit Probability term that determines if the origin of replication or the replication fork that should be open, opens. If larger than random number pooled from Gaussian distribution with mean 0 and standard deviation of 1.
* @param C Replication time
* @param cNoise Replication Gaussian standard deviation
* @param isDrugTreat Boolean type parameter that determines if the cell is experiencing drug treatment. If true then stop the initiation of new replication forks, but complete the ones that are currently open.
* @param chromDeg Probability parameter that determines if the chromosomes experiences catastrophic DNA damage that RecA cannot rescue. If larger than random number pooled from Gaussian distribution with mean 0 and standard deviation of 1.
* @param repForkDeg Probability parameter that determines if the replicating strand experiences catastrophic DNA damage that RecA cannot rescue. If larger than random number pooled from Gaussian distribution with mean 0 and standard deviation of 1.
* @param modelGeneLocations Location of the genes on the chromsome
* @param modelGeneLRPos Determine if the gene is on the left or right hand side of the oriC
* @param numCells Total number of cells in the simulation
*
* @return Error handling integer
*/

//TODO: seperate growing and duplicating the chromosome
int growChromosomes(Cell * cellArray, 
			int index, 
			float dt, 
			float chanceInit, 
			float C, 
			float cNoise, 
			bool isDrugTreat, 
			float chromDeg, 
			float repForkDeg, 
			float * modelGeneLocations, 
			int * modelGeneLRPos,
			int numCells)
{
	bool isSegTimer = false;
	bool isSegTimer2 = false;

	/*
	//What does tis do?????
	if(cellArray[index].numChrom>=2)
	{	
		cellArray[index].init_a = cellArray[index].a;
		//cellArray[index].segregationTimer += dt;
		isSegTimer = true;
	}
	*/

	//ONLY ONE check
	int i, y, z;
	float tmpGa = 0.0;
	//TODO: find a way where they get triggered without resetting everything to 0
	//set the copy number of the gene to be the same as the copy number of the chromosom, since each copycontains one copy number of the gene
	//here we remove the copy number of the gene of interest since each non replicating chromosome contains
	for(i=0; i<NUM_MODELGENES; i++)
	{
		cellArray[index].modelSumGenes[i] = cellArray[index].numChrom;
	}
	int sumOriC = 0;
	int sumPotentialOriC = 0;	

	if(cellArray[index].numChrom<MAX_CHROM)
	{
		//increment the replication timers
		for(i=0; i<cellArray[index].numChrom; i++)
		{
			//############# SANITY CHECK ####################
                        if(log2(cellArray[index].chromArray[i].potentialOriC)>LOG2_MAX_REP_TIMERS) //<--4
                        {
                                printf("\nWARNING (CheckCell): log2(cellArray[%d].chromArray[%d].potentialOriC) (%f)>%d\n",
                                        index, i, log2(cellArray[index].chromArray[i].potentialOriC), LOG2_MAX_REP_TIMERS);
                                printf("C: %f, D: %f, tau: %f\n", cellArray[index].C, cellArray[index].D, cellArray[index].tau);
                                //cellDebug(cellArray, index);
                                return 1;
                        }


			//############## ECLIPSE PERIOD ############
			//Here we determine if the eclipse period, that is the minimal time requried between two new replication events, has been elapsed
			bool isPastEclipse = true;
			if(cellArray[index].chromArray[i].oriC>1)
			{
				int prevOriC;
				int prevOriC_two;
				float eclipsePerc = (float)((dt*100.0/((float)cellArray[index].C*0.6))*100.0); //TODO: need to move this as paramerer at the cell initiation to avoid calculation at every time step
				if(cellArray[index].chromArray[i].oriC==2)
				{
					prevOriC = 1;
					prevOriC_two = 0;		
				}
				else
				{	
					//int prev = (int)(findPowerFloor(cellArray[index].chromArray[i].potentialOriC-0.5));
					prevOriC = (int)(log2(findPowerFloor(cellArray[index].chromArray[i].potentialOriC-0.5)));
					prevOriC_two = (int)(log2(findPowerFloor(prevOriC)));
				}
				/* debug for first cell
				if(index==0)
				{
					printf("EclipsePerc: %f\n", (float)((dt*100.0/((float)cellArray[index].C*0.6))*100.0));
					printf("potentialOriC: %d\n", cellArray[index].chromArray[i].potentialOriC);
					//printf("prev: %d\n", prev);
					printf("prevOriC: %d\n", prevOriC);
					printf("prevOriC_two: %d\n", prevOriC_two);
				}
				*/
				int z_p;
				//prevOriC = prev;
				for(z_p=0; z_p<prevOriC; z_p++)
				{
					/*
					if(index==0)
					{
						printf("[%d][%d][0] -> %f\n", prevOriC_two, z_p, cellArray[index].chromArray[i].replicationTimers[prevOriC_two][z_p][0]);
						printf("[%d][%d][1] -> %f\n", prevOriC_two, z_p, cellArray[index].chromArray[i].replicationTimers[prevOriC_two][z_p][1]);
					}
					*/
					if(cellArray[index].chromArray[i].replicationTimers[prevOriC_two][z_p][0]<eclipsePerc &&
						cellArray[index].chromArray[i].replicationTimers[prevOriC_two][z_p][1]<eclipsePerc)
					{
						isPastEclipse = false;
					}
				}
			}
			/*
			if(index==0)
			{
				printf("isPastEclipse: %d\n", isPastEclipse);
			}
			*/
			//############# DEGRADE CHROMOSOME ##################
			//case that the whole chromosome degrades 					
			if(normalDistRandn(0.0,1.0)>chromDeg &&
				numCells>10 &&
				cellArray[index].chromArray[i].oriC!=0 &&
				cellArray[index].chromArray[i].replicationTimers[0][0][0]>0.0 &&
				cellArray[index].chromArray[i].replicationTimers[0][0][1]>0.0)
			{
				if(cellArray[index].numChrom==0)
				{
					return 2;
				}
				//calculate chromosome genetic content and remove it
				float chromGa = 0.0;
				for(y=0; y<log2(findPowerCeil(cellArray[index].chromArray[i].potentialOriC)); y++)
				{
					for(z=0; z<pow(2,y); z++)
					{
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]>0.0)
                                                {
                                                        chromGa += (float)(cellArray[index].chromArray[i].replicationTimers[y][z][0]/100.0)/2.0;
                                                }
                                                if(cellArray[index].chromArray[i].replicationTimers[y][z][1]>0.0)
                                                {
                                                        chromGa += (float)(cellArray[index].chromArray[i].replicationTimers[y][z][1]/100.0)/2.0;
                                                }
					}	
				}
				cellArray[index].totalPotentialOriC -= cellArray[index].chromArray[i].potentialOriC;
				cellArray[index].totalRepForks -= cellArray[index].chromArray[i].repForks;
				clearChromosome(cellArray, index, i);
				cellArray[index].numChrom -= 1;	
				cellArray[index].Ga -= chromGa+1.0;
				reorganiseChromosomes(cellArray, index); //<-- need this to make sure that we can loop more effectively
				continue;
			}				

			//############ GROW AND DUPLICATE THE CHROMOSOME ##############
			//if the timer at position [0][0] exceeds the replication time duplicate the chromosome
			//twist, must have PAIR of replication timers be done to replicate the chromosome
			if(cellArray[index].chromArray[i].replicationTimers[0][0][0]>=100.0 && 
		    	  cellArray[index].chromArray[i].replicationTimers[0][0][1]>=100.0)
			{
				cellArray[index].chromArray[i].replicationTimers[0][0][0] = 0.0;
				cellArray[index].chromArray[i].replicationTimers[0][0][1] = 0.0;	

				/*
				//TODO: include this in another loop of all chromosomes and set it by number
				// NOTE: with chromosome degradation that wouldnt work anymore
				int emptyChromPos = -1;
				//find the position of the next empty chromosome
				for(y=0; y<MAX_CHROM; y++)
				{
					if(cellArray[index].chromArray[y].oriC==0 && emptyChromPos==-1)
					{
						emptyChromPos = y;
						break;
					}
				}
				if(emptyChromPos==-1)
				{
					//printf("WARNING (ReplicateChromosome): Trying to replicate the chromosome when MAX_CHROM (%d) is reached (%d)\n", 
					//	MAX_CHROM, cellArray[index].numChrom);
					return 1;
				}
				*/
				int emptyChromPos = cellArray[index].numChrom; // <-- this is to replace what is above for faster execution
				//make check
				if(cellArray[index].chromArray[emptyChromPos].oriC!=0)
				{
					emptyChromPos = -1;
					printf("WARNING: The empty position used is indeed not empty\n");
					for(y=0; y<MAX_CHROM; y++)
					{
						if(cellArray[index].chromArray[y].oriC==0 && emptyChromPos==-1)
						{
							emptyChromPos = y;
							break;
						}
					}
					if(emptyChromPos==-1) //i.e. did not find an empty chromosome position
					{
						printf("WARNING (ReplicateChromosome): Trying to replicate the chromosome when MAX_CHROM (%d) is reached (%d)\n", 
							MAX_CHROM, cellArray[index].numChrom);
						return 1;
					}
				}
				int motherRepForks = 0;
				int daughterRepForks = 0;
				int motherOriC = 1;
				int daughterOriC = 1;
				
				//for(y=1; y<log2(cellArray[index].chromArray[i].potentialOriC); y++)	
				//if we cannot assume that potential oriC contains every power ceiling then when we loop for updating 
				//the timers; then we find the power ceil then do the same as before	
				for(y=1; y<log2(findPowerCeil(cellArray[index].chromArray[i].potentialOriC)); y++)
				{
					//daughter chromsome -> second half
					for(z=(int)pow(2,y)/2.0; z<(int)pow(2.0,y); z++)
					{
						cellArray[index].chromArray[emptyChromPos].replicationTimers[y-1][z-(int)pow(2,y-1)][0] = 
							cellArray[index].chromArray[i].replicationTimers[y][z][0];
						cellArray[index].chromArray[emptyChromPos].replicationTimers[y-1][z-(int)pow(2,y-1)][1] = 
							cellArray[index].chromArray[i].replicationTimers[y][z][1];
						//ADJUST ASYNC
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]>0.0)
						{
							daughterRepForks += 1;
						}
						if(cellArray[index].chromArray[i].replicationTimers[y][z][1]>0.0)
						{
							daughterRepForks += 1;
						}
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]>0.0 || 
							cellArray[index].chromArray[i].replicationTimers[y][z][1]>0.0)
						{
							daughterOriC += 1;
						}
						//set the mother replication timers in these positions to 0 since they have been moved to the daughter cell
						cellArray[index].chromArray[i].replicationTimers[y][z][0] = 0.0;
						cellArray[index].chromArray[i].replicationTimers[y][z][1] = 0.0;
					}
					//mother chromosome --> first half
					for(z=0; z<(int)pow(2.0,y)/2.0; z++)
					{
						cellArray[index].chromArray[i].replicationTimers[y-1][z][0] = 
							cellArray[index].chromArray[i].replicationTimers[y][z][0];
						cellArray[index].chromArray[i].replicationTimers[y-1][z][1] = 
							cellArray[index].chromArray[i].replicationTimers[y][z][1];
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]>0.0)
						{
							motherRepForks += 1;
						}
						if(cellArray[index].chromArray[i].replicationTimers[y][z][1]>0.0)
						{
							motherRepForks += 1;
						}
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]>0.0 || 
							cellArray[index].chromArray[i].replicationTimers[y][z][1]>0.0)
						{
							motherOriC += 1;
						}
						//we set the last chromosome array location to 0 since it has moved up and will not be overwritten
						if(y==(int)log2(cellArray[index].chromArray[i].potentialOriC)-1)
						{
							cellArray[index].chromArray[i].replicationTimers[y][z][0] = 0.0;
							cellArray[index].chromArray[i].replicationTimers[y][z][1] = 0.0;
						}
					}
				}
				//adjust chromosome variables
				cellArray[index].chromArray[i].repForks = motherRepForks;
				cellArray[index].chromArray[i].oriC = motherOriC;
				//we assume that the potentialOriC that was previously reached will are cancelled
				//ADJUST ASYNC
				//cellArray[index].chromArray[i].potentialOriC = findPowerCeil(cellArray[index].chromArray[i].oriC);
				cellArray[index].chromArray[i].potentialOriC = (int)cellArray[index].chromArray[i].oriC;

				cellArray[index].chromArray[emptyChromPos].repForks = daughterRepForks;
				cellArray[index].chromArray[emptyChromPos].oriC = daughterOriC;
				//cellArray[index].chromArray[emptyChromPos].potentialOriC = findPowerCeil(cellArray[index].chromArray[emptyChromPos].oriC);
				cellArray[index].chromArray[emptyChromPos].potentialOriC = (int)cellArray[index].chromArray[emptyChromPos].oriC;

				//the the cell looses 2 replication forks and gains a chromsome
				cellArray[index].totalRepForks -= 2;
				cellArray[index].numChrom += 1;

				//this is our own implementation of the initiation of segregation -- in the Keasling et al. paper, it initiates when all chromosomes have initiated. In our implementation it starts when the cell finishes the replication of a new chromosome 
				//TODO: Test this is with mass distribution --> considering the individual mass distribution this is WRONG
				if(cellArray[index].segregationTimer==0.0 && isSegTimer==false)
				{
					//record when the segregation timer starts
				
					cellArray[index].init_a = cellArray[index].a;
					cellArray[index].segregationTimer += dt;
					isSegTimer = true;
				}	
				else if(cellArray[index].segregationTimer2==0.0 && isSegTimer2==false)
				{
					cellArray[index].segregationTimer2 += dt;
					isSegTimer2 = true;
				}	
			}

			//############ GROW ##############
			if(cellArray[index].chromArray[i].potentialOriC!=0)
			{	
				int chromNumOriC = 1; //<<--- this is stupid but there is somewhere in this program where the oriC and sum of oriC gets fucked
				int chromNumPotentialOriC = 1;
				//int emptyLoc[2] = {-1, -1}; //<-- when there are mutliple pairs of replication forks, checks and see if the location is empty to be initiated

				//for(y=0; y<log2(cellArray[index].chromArray[i].potentialOriC); y++)
				//because we are not sure that the potentialOriC contains power ceiling values we must find it to loop through it
				for(y=0; y<log2(findPowerCeil(cellArray[index].chromArray[i].potentialOriC)); y++)
				//for(y=0; y<7; y++)
				{
					for(z=0; z<pow(2,y); z++)
					{
						//INCREMENT
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]>0.0 && 
						cellArray[index].chromArray[i].replicationTimers[y][z][0]<=100.0)
						{
							cellArray[index].chromArray[i].replicationTimers[y][z][0] += 
								(float)(dt*100.0/(float)cellArray[index].C);
						}
						if(cellArray[index].chromArray[i].replicationTimers[y][z][1]>0.0 && 
						cellArray[index].chromArray[i].replicationTimers[y][z][1]<=100.0)
						{
							cellArray[index].chromArray[i].replicationTimers[y][z][1] += 
								(float)(dt*100.0/cellArray[index].C);
						}
						
						//CHANCE THAT THE REPLICATION FORKS COLLAPSE AND THE DAUGHTER REPLICATING STRAND DEGRADES (i.e. reset timer)
						//NOTE: we are considering that the last replicating replication froks are the only ones that can be collapsing
						// alternative creates segementation fault
					
						if(normalDistRandn(0.0,1.0)>repForkDeg &&
						   y==(int)(log2(findPowerCeil(cellArray[index].chromArray[i].potentialOriC)))-1 &&
						   cellArray[index].chromArray[i].replicationTimers[y][z][0]>0.0 &&
						   cellArray[index].chromArray[i].replicationTimers[y][z][1]>0.0)
						{
							//printf("\tWARNING: Degrading Replication fork\n");
							cellArray[index].chromArray[i].replicationTimers[y][z][0] = 0.0;
							cellArray[index].chromArray[i].replicationTimers[y][z][1] = 0.0;

							cellArray[index].chromArray[i].oriC -= 1;	
							cellArray[index].chromArray[i].potentialOriC -= 1;
							cellArray[index].chromArray[i].repForks -= 2;

							cellArray[index].totalPotentialOriC -= 1;
							cellArray[index].totalRepForks -= 2;
						}						

						//CORRECT FOR EXCESS
						//need to correct for the timer than exceeds 100.0%
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]>100.0)
						{
							cellArray[index].chromArray[i].replicationTimers[y][z][0] = 100.0;
						}
						if(cellArray[index].chromArray[i].replicationTimers[y][z][1]>100.0)
						{
							cellArray[index].chromArray[i].replicationTimers[y][z][1] = 100.0;
						}

						//OPEN UNINITIATED
						//ADJUST ASYNC
						//we consider that potentialOriC is not necessarily power ceiling 
						//TODO: Make sure that this is done correctly
						//NOTE: we know that because of the structure 'selected', the ones that are first are more likely to open
						// consider the outcome of such a method
						//printf("isDrugTreat: %d\n", isDrugTreat);
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]==0.0 && 
							cellArray[index].chromArray[i].replicationTimers[y][z][1]==0.0 && 
							cellArray[index].chromArray[i].oriC<cellArray[index].chromArray[i].potentialOriC &&
							isPastEclipse==true &&
							isDrugTreat==false)
						{
							/*
							int selected_y = -1;
							int selected_z = -1;
							if(emptyLoc[0]==-1 && emptyLoc[1]==-1)
							{
								emptyLoc[0] = y;
								emptyLoc[1] = z;
								selected_y = y;
								selected_z = z;
							}
							else
							{
								selected_y = emptyLoc[0];
								selected_z = emptyLoc[1];	
							}
							*/
							if(normalDistRandn(0.0,1.0)>chanceInit)
							{
								/*
								cellArray[index].chromArray[i].replicationTimers[selected_y][selected_z][0] += 
									(float)(dt*100.0/(float)cellArray[index].C);
								cellArray[index].chromArray[i].replicationTimers[selected_y][selected_z][1] += 
									(float)(dt*100.0/(float)cellArray[index].C);
								*/
		
								cellArray[index].chromArray[i].replicationTimers[y][z][0] += 
									(float)(dt*100.0/(float)cellArray[index].C);
								cellArray[index].chromArray[i].replicationTimers[y][z][1] += 
									(float)(dt*100.0/(float)cellArray[index].C);
							
								//this is only valid if it is asynchronous
								cellArray[index].chromArray[i].repForks += 2;
								cellArray[index].totalRepForks += 2;
								//we do not no this here because the creation of the pair of replication forks could be asynchronous 
								//cellArray[index].chromArray[i].oriC += 1;
								//emptyLoc[0] = -1;
								//emptyLoc[1] = -1;
							}	
						}

						//CALCULATE THE GENETIC CONTENT: note that this is only the genetic content of the replicating strands
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]>0.0)
						{
							tmpGa += (float)(cellArray[index].chromArray[i].replicationTimers[y][z][0]/100.0)/2.0;
						}
						if(cellArray[index].chromArray[i].replicationTimers[y][z][1]>0.0)
						{
							tmpGa += (float)(cellArray[index].chromArray[i].replicationTimers[y][z][1]/100.0)/2.0;
						}
						//DETERMINE THE COPY NUMBER AS DETERMINED BY modelGeneLocations where their locations are modelGeneParamsLocations in modelParams
						//because the timers only represent the replicating strands, all chromosomes contain at least one copy of the gene (==numChrom)
						//and if the repilication timer is greater than the relative position of the gene then we add a copy number
						for(int gN=0; gN<NUM_MODELGENES; gN++)
						{
							if(modelGeneLRPos[gN]==0 && cellArray[index].chromArray[i].replicationTimers[y][z][0]>modelGeneLocations[gN])
							{
								cellArray[index].modelSumGenes[gN] += 1;
							}
							if(modelGeneLRPos[gN]==1 && cellArray[index].chromArray[i].replicationTimers[y][z][1]>modelGeneLocations[gN])
							{
								cellArray[index].modelSumGenes[gN] += 1;
							}
						}
						//Note: as stated above, this is here only if replication forks initiation is asynchronous. Not the case here.
						if(cellArray[index].chromArray[i].replicationTimers[y][z][0]>0.0 && 
							cellArray[index].chromArray[i].replicationTimers[y][z][1]>0.0)
						{
							chromNumOriC += 1;				
						}
					}
				}
				
				// as stated abouve this is a useless excercise and should be remove TODO.
				cellArray[index].chromArray[i].oriC = chromNumOriC;
				sumOriC += chromNumOriC;
				sumPotentialOriC += cellArray[index].chromArray[i].potentialOriC;
			}	
		}
		//update the amount of DNA per cell
		cellArray[index].Ga = (float)cellArray[index].numChrom+tmpGa;
		cellArray[index].totalPotentialOriC = sumPotentialOriC;

		//need to reorganise the chromosome if indeed it is deleted
		reorganiseChromosomes(cellArray, index);
	}
	else
	{
		printf("Warning: (cell %d) reached (%d) MAX_CHROM (%d)\n", index, cellArray[index].numChrom, MAX_CHROM);
		return 1;
	}

	//increment the segregation timer only if not newly segregated
	if(cellArray[index].segregationTimer>0.0 && isSegTimer==false)
	{
		cellArray[index].segregationTimer += dt;
	}
	if(cellArray[index].segregationTimer2>0.0 && isSegTimer2==false)
	{
		cellArray[index].segregationTimer2 += dt;
	}
	// in Keasling paper states forks>=chrom. We can interpret this as either the pair of forks or the individual replication forks
	// --> warning may cause the cells to initiate segregation timer right after division.... unsure if this actually happens
	//TODO: Test this with mass distribution
	/*
	if((int)(cellArray[index].totalRepForks/2.0)>=cellArray[index].numChrom &&
		cellArray[index].numChrom>=2 && 
		cellArray[index].segregationTimer==0.0)
	{
		cellArray[index].segregationTimer += dt;
	}
	*/
	return 0;
}

/**
* @brief Grow plasmids and duplicate it once it finises replication
*
* Grow the plasmids and replicate once the timer finishes. We do not implement the parallel replication 
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @param dt Time step
*
* @return Error handling integer
*/
//TODO: Need to add the same plasmids copy number for the model as implemented for the plasmids
int growPlasmids(Cell * cellArray, int index, float dt)
{
	int i, y;
	if(cellArray[index].numPlasmid<MAX_PLASMID)
	{
		//increment the replication timers
		for(i=0; i<cellArray[index].numPlasmid; i++)
		{
			if(cellArray[index].plasmidArray[i].oriC!=0)
			{
				//############# DUPLICATE ##################
				if(cellArray[index].plasmidArray[i].replicationTimers[0]==100.0 &&
					cellArray[index].plasmidArray[i].replicationTimers[1]==100.0)
				{
					//find the position of the next empty chromosome
					int emptyPlasmidPos = cellArray[index].numPlasmid;
					//This is a backup method that should be removed if indeed the above works well
					if(cellArray[index].plasmidArray[emptyPlasmidPos].oriC!=0)
					{
						printf("WARNING: The empty plasmid position indeed not empty\n");
						emptyPlasmidPos = -1;
						for(y=0; y<MAX_PLASMID; y++)
						{
							if(cellArray[index].plasmidArray[y].oriC==0 && emptyPlasmidPos==-1)
							{
								emptyPlasmidPos = y;
								break;
							}
						}
						if(emptyPlasmidPos==-1)
						{
							return 1;
						}
					}

					cellArray[index].plasmidArray[i].replicationTimers[0] = 0.0;
					cellArray[index].plasmidArray[i].replicationTimers[1] = 0.0;
					cellArray[index].plasmidArray[i].oriC = 1;				
					
					cellArray[index].plasmidArray[emptyPlasmidPos].replicationTimers[0] = 0.0;				
					cellArray[index].plasmidArray[emptyPlasmidPos].replicationTimers[1] = 0.0;				
					cellArray[index].plasmidArray[emptyPlasmidPos].oriC = 1;				

					cellArray[index].numPlasmid += 1;				
				}
				
				//################ GROW #######################
				if(cellArray[index].plasmidArray[i].replicationTimers[0]>0.0 && 
					cellArray[index].plasmidArray[i].replicationTimers[1]>0.0)
				{
					cellArray[index].plasmidArray[i].replicationTimers[0] += 
						(float)(dt*100.0/(float)cellArray[index].C_plasmid);
					cellArray[index].plasmidArray[i].replicationTimers[1] += 
						(float)(dt*100.0/(float)cellArray[index].C_plasmid);
				}
				//################ INITIATE ######################
				if(cellArray[index].plasmidArray[i].oriC==2)
				{
					cellArray[index].plasmidArray[i].replicationTimers[0] += 
                                                (float)(dt*100.0/(float)cellArray[index].C_plasmid);
                                        cellArray[index].plasmidArray[i].replicationTimers[1] += 
                                                (float)(dt*100.0/(float)cellArray[index].C_plasmid);
				}

				//CORRECT EXCESS
				if(cellArray[index].plasmidArray[i].replicationTimers[0]>100.0)
				{
					cellArray[index].plasmidArray[i].replicationTimers[0] = 100.0;
				}
				if(cellArray[index].plasmidArray[i].replicationTimers[1]>100.0)
				{
					cellArray[index].plasmidArray[i].replicationTimers[1] = 100.0;
				}
				//CALCULATE THE GENETIC CONTENT
				float tmpPa = 0.0;
				if(cellArray[index].plasmidArray[i].replicationTimers[0]>0.0)
				{
					tmpPa += (float)(cellArray[index].plasmidArray[i].replicationTimers[0]/100.0)/2.0;
				}
				if(cellArray[index].plasmidArray[i].replicationTimers[1]>0.0)
				{
					tmpPa += (float)(cellArray[index].plasmidArray[i].replicationTimers[1]/100.0)/2.0;
				}
				cellArray[index].Pa = (float)cellArray[index].numPlasmid+tmpPa;
			}
		}
	}
	else
	{
		printf("ERROR: The cell has reached maximal plasmid copy number\n");
		return 1; 
	}
	return 0;
}

//#################################### CELLULAR FUNCTIONS #############################
//######################################## PUBLIC #####################################

/**
* @brief Grow the cell either assuming exponential growth or injection growth
*
* Public function to grow the cells in the simulation assuming exponential or injection growth. Opens OriC that should be opened and indice cell division if segregation timer has been reached. 
*
* @param cellArray Array containg all the cells
* @param index Position of the cell of interest in the cellArray
* @param newIndex Position of an empty position in the cellArray for potential new cell
* @param tau Exponential doubling rate
* @param cNoise Replication Gaussian standard deviation
* @param dNoise Segregation Gaussian standard deviation
* @param Vi Critical volume
* @param Vi_plasmid Plasmid critical mass
* @param ViNoise Critical volume Gaussian standard deviation
* @param divRatio Division mean volume ratio distribution between mother and daughter cell (default=0.5)
* @param divNoise Division standard deviation volume ration distribution between mother and daughter cell 
* @param partRatio Chromosome partition mean ration between mother and daughter cell (default=0.5)
* @param partNoise Chromosome partition standard deviation between mother and daughter cell
* @param chromDeg Probability parameter that determines if the chromosomes experiences catastrophic DNA damage that RecA cannot rescue. If larger than random number pooled from Gaussian distribution with mean 0 and standard deviation of 1.
* @param repForkDeg Probability parameter that determines if the replicating strand experiences catastrophic DNA damage that RecA cannot rescue. If larger than random number pooled from Gaussian distribution with mean 0 and standard deviation of 1.
* @param numCells Total number of cells in the simulation
* @param C1 One-phase exponential function replication time first term
* @param C2 One-phase exponential function replication time second term
* @param C3 One-phase exponential function replication time third term
* @param D1 One-phase exponential function segregation time first term
* @param D2 One-phase exponential function segregation time second term
* @param D3 One-phase exponential function segregation time third term
* @param dt Time step
* @param volumeAdd If injection growth populaion mean volumetric growth
* @param modelGeneParamsLocations Location of the parameters in the GSL parameter array that determines the gene number
* @param modelGeneLocations Location of the genes on the chromsome
* @param modelGeneLRPos Determine if the gene is on the left or right hand side of the oriC
* @param isDrugTreat Boolean type parameter that determines if the cell is experiencing drug treatment. If true then stop the initiation of new replication forks, but complete the ones that are currently open.
*
* @return Error handling integer
*/
//float trad_C = 43.2*(1.0+(4.86*exp(-4.5/db_h)));
//float trad_D = (18.5+4.5)*(1.0+(exp(-4.5/db_h)));
int growCell(Cell * cellArray,
		int index,
		int newIndex,
		float tau,
		float cNoise,
		float dNoise,
		float Vi,
		float Vi_plasmid,
		float ViNoise,
		float VaNoise,
		float chanceInit,
		float divRatio,
		float divNoise,
		float partRatio,
		float partNoise,
		float chromDeg,
		float repForkDeg,
		int numCells,
		float C1,
		float C2,
		float C3,
		float D1,
		float D2,
		float D3,
		float dt,
		float volumeAdd,
		int * modelGeneParamsLocations,
		float * modelGeneLocations,
		int * modelGeneLRPos,
		bool isDrugTreat)
{
	int stopFlag = 0;
	float mu = -1.0;
	if(volumeAdd==-1.0)
	{
		mu = log(2.0)/cellArray[index].tau;
	}
	else
	{
		//we calculate what is the growth rate based on the mass addition rate of the HMG to use for other functions that are growth rate dependent
		//TODO: See if it is better to set the replication rate at every time step or it is better to set it once a division rate
		mu = -((cellArray[index].Va-(cellArray[index].Va+volumeAdd))/(cellArray[index].Va*dt));
		//float tau = -((cellArray[index].Va-(cellArray[index].Va+volumeAdd))-cellArray[index].Va)/(cellArray[index].Va*dt);
		tau = log(2.0)/mu;
	}

	float C = 0.0;
	float D = 0.0;
	if(C2==-1.0)
	{
		C = C1;
		D = D1;		
	}
	else
	{
		//keasling et al. 
		C = C1*(1.0+(C2*exp(-C3/(tau/60.0)))); //convert per minute and not per hour
		D = D1*(1.0+(D2*exp(-D3/(tau/60.0))));
		//printf("Calculated D: %f\n", D);
		//WARNING: CHANGE THIS ONLY FOR TESTING
		//cellArray[index].C = C;
		//cellArray[index].D = D;
		//cellArray[index].tau = tau;
	}
	//float C_plasmid = 6646.0*C/4639221.0; //6646pb = size of ColE1, 4639221bp = size of chromosome;
	float C_plasmid = 99200.0*C/4639221.0; //99200pb = size of F plasmid
	
	int i, y, z;
	cellArray[index].a += dt;

	//############################## MODEL TIME STEP ####################
	//TODO: seems like we need to add cellArray[index].a to the gsl driver and they append the time step
	//Update the copy number of the genes
	//TODO: change modelGeneParamsLocations to modelGeneLocations 
	for(i=0; i<NUM_MODELGENES; i++)
	{
		//This is the original one where each promoter is 10.0
		cellArray[index].modelSpecies[modelGeneParamsLocations[i]] += (cellArray[index].modelSumGenes[i]-cellArray[index].modelPrevSumGenes[i])*10.0; //WARNING: *10.0 is only because of bennett et al. sum of one 
		cellArray[index].modelPrevSumGenes[i] = cellArray[index].modelSumGenes[i]; 
	}	

	if(NUM_MODELGENES!=0)
	{
		double d_a = (double)cellArray[index].a;
		double d_a_dt = (double)cellArray[index].a+(double)dt;
		int status = gsl_odeiv2_driver_apply(cellArray[index].driver, &d_a, d_a_dt, cellArray[index].modelSpecies);
		if(status!=GSL_SUCCESS)
		{
			printf("ERROR (growCell): GSL error, return value=%d\n", status);
			//cDebug(cellArray, index);
			return 2;
		}
	}

	if(volumeAdd==-1.0)
	{
		cellArray[index].Va = cellArray[index].Va*(1.0+mu*dt);
	}
	else
	{
		cellArray[index].Va = cellArray[index].Va+volumeAdd*cellArray[index].injectionDeviation;
		//cellArray[index].Va = cellArray[index].Va+volumeAdd;
	}

	if(cellArray[index].Va>=(cellArray[index].Vi*cellArray[index].totalPotentialOriC))
	{
                for(i=0; i<cellArray[index].numChrom; i++)
                {
                        if(cellArray[index].chromArray[i].oriC!=0 && 
				(int)log2(findPowerCeil(cellArray[index].chromArray[i].oriC+1))<LOG2_MAX_REP_TIMERS) //<-- make sure the chromosome is present
                        {
				cellArray[index].chromArray[i].potentialOriC = (int)cellArray[index].chromArray[i].oriC*2;
				cellArray[index].isInitiated = true;
			}
		}
                cellArray[index].isRep = true;
	}

	//at this point in time we will not implement overalapping rounds of replication for plasmids for two resaons. The first, since different plasmids have a replication time that is so fast, applying this would be practivally useless
	//if(cellArray[index].Va>=cellArray[index].Vi_plasmid*cellArray[index].numPlasmid)
	if(cellArray[index].Va>=cellArray[index].Vi_plasmid*cellArray[index].totalPotentialOriC_plasmid)
	{
		//terrible implementation but it works. TODO: make random
		for(i=0; i<cellArray[index].numPlasmid; i++)
		{
			if(cellArray[index].plasmidArray[i].oriC==1)
			{
				cellArray[index].plasmidArray[i].oriC = 2;
				cellArray[index].totalPotentialOriC_plasmid += 1;
				break;
			}
		}
	}

	//printf("\tDEBUG (growCell): Opening new potentialOriC\n");
	
	//clock_t begin = clock();

	//printf("\tDEBUG (growCell): Growing chromosomes\n");
	stopFlag = growChromosomes(cellArray, index, dt, chanceInit, C, cNoise, isDrugTreat, chromDeg, repForkDeg, modelGeneLocations, modelGeneLRPos, numCells);
	if(stopFlag!=0)
	{
		return stopFlag;
	}
	stopFlag = growPlasmids(cellArray, index, dt);
	if(stopFlag!=0)
	{
		return stopFlag;
	}

	//########################## cell division ############################
	//second condition for the drug treatment	

	/*
	//ADDER MODEL 
	if(cellArray[index].Vadded>=Vdelta)
	{

	}
	*/	

	/*
	//SIZER MODEL
	if(cellArray[index].Va>=Vtarget)
	{

	}
	*/	

	//############################################# CLASSIC ############################################
	//printf("\tDEBUG (growCell): Dividing cell\n");
	if(cellArray[index].segregationTimer>=cellArray[index].D && cellArray[index].D>0.0)
	//if(cellArray[index].Va>=cellArray[index].divVol && cellArray[index].numChrom>=2)
	//if(cellArray[index].Va>=cellArray[index].divVol)
	{
		/*
		if(cellArray[index].numChrom==1)
		{
			printf("WARNING: Dividing cell only contains a single chromosome\n");
			float newVolume = cellArray[index].Va*normalDistRandn(divRatio,((divNoise*divRatio)/100.0));
			cellArray[index].isDead = false;
			cellArray[index].isInitiated = false;
			cellArray[index].isRep = false;
			cellArray[index].Va = newVolume;
			cellArray[index].a = 0.0;
			cellArray[index].C = normalDistRandn(C, (C*cNoise/100.0));
			cellArray[index].C_plasmid = normalDistRandn(C_plasmid, (C_plasmid*cNoise/100.0));
			cellArray[index].D = normalDistRandn(D, (D*dNoise/100.0));
			cellArray[index].Vi = normalDistRandn(Vi, (Vi*ViNoise/100.0));
			cellArray[index].Vi_plasmid = normalDistRandn(Vi_plasmid, (Vi_plasmid*ViNoise/100.0));
			cellArray[index].tau = tau;
			cellArray[index].segregationTimer = 0.0;
			cellArray[index].segregationTimer2 = 0.0;
			//cellArray[index].divVol = adderDivVol(tau);
			return 3;
		}
		*/
		//else if(cellArray[index].numChrom==0)
		if(cellArray[index].numChrom==0)
		{
			printf("WARNING: Dividing cell contains no chromosomes\n");
			return 2;	
		}
		else if(cellArray[index].numChrom>=2)
		{
			stopFlag = divideCell(cellArray, index, newIndex, tau, C, C_plasmid, cNoise, D, dNoise, Vi, Vi_plasmid, ViNoise, VaNoise, divRatio, divNoise, partRatio, partNoise, modelGeneParamsLocations, modelGeneLocations, modelGeneLRPos, dt);
		}
	}
	//printf("\tDEBUG (growCell): Done dividing cell\n");
	if(stopFlag==3)
	{
		return 3;
	}
	else if(stopFlag==2)
	{
		return 2;
	}

	//################## libsibmlsim ###################
	//move this to the 
	  //dt -> must be brought directly
	  //print_interval = 0;
	  //print_amount = 0;
	  //check that dt, dt is a single time step
	  //cellArray[index].rtn = simulateSBMLModel(cellArray[index].m, dt, dt, print_interval, print_amount, cellArray[index].method, cellArray[index].use_lazy_method, cellArray[index].atol, cellArray[index].rtol, cellArray[index].facmax);
	  //TODO: time the execution time of this function to see if scalable. If not, need to go deeper and set the parameters that contained in this function
	cellArray[index].sbml_results = simulateSBMLModel(cellArray[index].sbml_model, dt, dt, 0, 0, cellArray[index].sbml_simulation_method, cellArray[index].sbml_use_lazy_method, cellArray[index].sbml_atol, cellArray[index].sbml_rtol, cellArray[index].sbml_facmax);
	if(cellArray[index].sbml_results==NULL)
	{
	    cellArray[index].sbml_results = create_myResult_with_errorCode(SimulationFailed);
	    return 1; //<-- not sure about this
	}


	//################# DEBUG: print the cell 0 #############################	
	/*
	if(index==0)
	{
		cDebug(cellArray, index);
	}
	*/
	//cDebug(cellArray, index);

	return 0;
}
