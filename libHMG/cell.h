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
* @file cell.h
*
* Header file of cell
*
* @version 1.0
* @author Melchior du Lac
*
* Object structure of the cell containing the population parameters, such as replication rate (C), that have had Gaussian noise induced to them.
* 26/10/16
	-> Overall tried to reduce the size of the model in the RAM
	-> Removed cellAge
	-> Removed mu, redundant with tau
	-> Included the boolean datatype (stdbool.h) and smaller integer types (inttypes.h)
*/

/*
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
*/
#include "inputModel.h"
#include "libsbmlsim/libsbmlsim/libsbmlsim.h"


#ifndef CELL_H
#define CELL_H

/*
Changed the dynamic memory allocation of the replication timers, Chrom Obj, and chromArray to stack memory for speed
*/

//TODO: change this to dynamic allocation if this is going to be used for a wide range of different 
//extern int MAX_CHROM; //<- 32
//extern int MAX_PLASMID; //<- 2000
//extern int LOG2_MAX_REP_TIMERS; //<- 6
//extern int NUM_MODELPARAMS;
//extern int NUM_MODELSPECIES;
//extern int NUM_MODELGENES;

#define MAX_CHROM 32
#define LOG2_MAX_REP_TIMERS 6
#define MAX_PLASMID 2000

/**
* @brief Structure of the plasmid
*
* Object containing the information for the plasmids (ColE1). Assume that it cannot contain overlapping rounds of replication (timer only 2d)
*
*/
typedef struct PLASMID
{
	float replicationTimers[2]; /**< Replication timers for the pair of replication forks*/
	//TODO: because this can either be 1 or 2, change to either boolean type or define new type for this particular case to save memory
	int oriC; /**< Number of origin of replications (should be either 1 or 2)*/
} Plasmid;

/**
* @brief Structure of the chromosome
*
* Object structure of the chromosme that is contained in the chromArray array 
*
*/
typedef struct CHROM
{
	//chromosome contains 2d array with first dimension being the first replication fork and the rest are the subsequent power of 2 ones.
	//WARNING: if the extern parameters above change this must change as well --> cannot put it as parameters because not variable
	//TODO: figure out why the second term is 64 and not some other number
	//TODO: check to see if the 64 number threshold is ever reached
	//float replicationTimers[7][64][2]; 
	//float replicationTimers[6][32][2]; 
	float replicationTimers[LOG2_MAX_REP_TIMERS][MAX_CHROM][2]; /**< 3D array containing all the potential origin of replications. First dimension (default size = 7) represents the the natural log of the maximal number of replication forks (note that this is arbitrary). Second dimension is the total number of replication forks. The last dimension is the pair of replication forks since the two may be stocastic.*/
	//WARNING: Rememember that this is a waste of memory. Only a better solution if memory is not an issue and speed is more important	
	int oriC; /**< Number of oriC on the chromosome*/
	int repForks; /**< Number of replication forks on the chromosome*/
	int potentialOriC; /** Number of oriC that should be on the chromosome*/
} Chrom;

/**
* @brief Structure of the cell
* 
* Object structure of the cell containing all the parameters of the cell. Includes the GSL objects for the simulation
*
*/
typedef struct CELL
{
	//previous state of the cell
	float prev_sum_a; /**< Previous complete duplication time (to test if the individual cells doubling rate)*/
	float prev_newbornVol; /**< Previous newborn volume of the cell */
	float prev_divisionVol; /**< Previous division volume of the cell*/
	// ######## Sizer model
	float Vadded; /**< Current added size of the model*/
	float Vdelta; /**< Total mass to be added*/
	float Vtarget; /**< Sizer model target*/
	//######## Adder model
	float D; /**< Timer segregation time*/
	float segregationTimer; /**< Segregation timer*/
	float segregationTimer2; /**< Segregation timer to be inherited */
	//TODO: perhaps use this to increase the likelyhood of the cell dying
	//float cellAge; /**< Total age of the cell (cumulative of independent cell cycles). TODO: comment out since only usefull if we are to implement age dependent cell death*/
	float a; /**< Cell internal clock is B+C+D. It it reset at every cell division, mother like duaghter*/
	float init_a; /**< previous complete duplication time (to test if the individual cells doubling rate)*/
	float Va; /**< Mass of the cell at time a*/
	//float V0; /**< Initial mass of the cell at division -- unused*/
	float Vi; /**< Mass at intiation. Also called initiation mass or critical mass*/
	float Vi_plasmid; /**< Mass at intiation. Also called initiation mass or critical mass*/
	float Ga; /**< DNA content at time a*/
	float Pa; /**< DNA content for a plasmid at time a*/
	//float G0; /**< DNA content at end of segregation -- unused*/
	//int totalOriC; /**< Total number of OriC in the cell. TODO: remove parameter to skim the program*/
	int totalPotentialOriC; /**< Total number of oriC that should be opened. Important for the opening of new replication forks*/
	int totalPotentialOriC_plasmid; /**< Total number of oriC that should be opened. Important for the opening of new replication forks*/
	int totalRepForks; /**< Total number of replication forks. TODO: These two are essentially the same (check). Remove one*/
	int numChrom; /**< Total number of chromosomes*/
	int numPlasmid; /**< Total number of chromosomes*/
	float injectionDeviation; /**< Gaussian deviation from the injection growth method*/
	//float divVol;
 	Chrom chromArray[MAX_CHROM]; /**< Array of struct for chromosome (default size = 64). TODO: make this dynamic*/
 	Plasmid plasmidArray[MAX_PLASMID]; /**< Array of struct for Plasmid object TODO: make this a dynamic array that must be initialised by the constructers */
	bool isDead; /**< Boolean type flag that determines if the cell is dead or not*/
	bool isNewlyRep; /**< Boolean type flag that determines if the cell has recently replicated its chromosomes to avoid float replication*/
	/** Check if the cell in question has initiated replication*/
	bool isInitiated; /**< Boolean type parameter that flags if the cell has initiated its chromosomes in its lifetime. Note redundant?*/
	bool isRep; /**< Boolean type parameter that flags if the cell is currently replicating itself. Note: redundant?*/
	float C; /**< Chromosme replicating time*/
	float C_plasmid; /**< Plasmid replication time*/
	float tau; /**< Cell exponential doubling rate*/
	//float mu; /**< Cell exponential instantaneous growth rate. Note redundant?*/
	bool isFrozen; /**< Boolean type parameter that flags if the cell is frozen due to reaching its maximal chromosome number. Note: Redundant?*/
	// array of the return parameters of the model function
	// TODO: do not hardcode the size of these
	//#################### GSL ########################
	
	//double modelParams[NUM_MODELPARAMS+1]; /**< Array of parameters of the model*/ //TODO: have a check that the size of the model is correct
	//double modelSpecies[NUM_MODELSPECIES+1]; /**< Model parameters*/ //TODO: this is a stupid way of implementing the storage of the model parameters, since the user must know the location of the paramters that represent the copy number of the gene in question. Better to have it as an object.
	//int modelSumGenes[NUM_MODELGENES+1]; /**< Number of genes for the cell based on the cellPopulation modelGeneLocations*/
	//int modelPrevSumGenes[NUM_MODELGENES+1]; /**< When calculating the copy number amount, there is a need to determine if the copy number has changed for the cell*/
	//double modelParams[8]; /**< Array of parameters of the model*/ //TODO: have a check that the size of the model is correct
	//double modelSpecies[15]; /**< Model parameters*/ //TODO: this is a stupid way of implementing the storage of the model parameters, since the user must know the location of the paramters that represent the copy number of the gene in question. Better to have it as an object.
	//int modelSumGenes[3]; /**< Number of genes for the cell based on the cellPopulation modelGeneLocations*/
	//int modelPrevSumGenes[3]; /**< When calculating the copy number amount, there is a need to determine if the copy number has changed for the cell*/

	//gsl_odeiv2_system sys; /**< GSL system object*/
	//gsl_odeiv2_driver * driver; /**< GSL driver object*/
	
	//################### libsbmlsim#######################
	
	myResult *sbml_results;
	unsigned int sbml_err_num;
	double sbml_ato;// = 0.0;
	double sbml_rtol;// = 0.0;
	double sbml_facmax;// = 0.0;
	
} Cell;

int constructCell(Cell * cellArray, int index);

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
                    float a);

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
                bool isDrugTreat);

#endif
