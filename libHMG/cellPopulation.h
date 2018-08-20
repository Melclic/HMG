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
* @file cellPopulation.h
*
* Cell Population header file that grows the population. Note that the main reason we are enforcing this structure 
* so that if one needs to study the interactions between two different cell types that have different cell cycle dynamics then one is able to specufy different array space in the allocated memory and cell cycle, noise parameters and SBML model
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

#include "../libsbmlsim/libsbmlsim.h"

/**
 * @brief Cell population object
 *
 * Cell population object containing all the parameters of the population of cells. Includes all the default values and population level parameters before they are assigned to the individual cells in the model with typically Gaussian noise assigned to them. Also contains the array of cells in the nmodel.
 *
 */
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
    
    //int numModelParams; /**< Number of params for the modelInitialsParams*/
    //int numModelSpecies; /**< Number of species for the modelInitialSpecies array*/
    //int numModelGenes; /**< Number of species for the modelInitialSpecies array*/
    
    float tau; /**< cell doubling time (min) */
    float C1; /**< functional form of the replication rate (C)*/
    float C2; /**< functional form of the replication rate (C)*/
    float C3; /**< functional form of the replication rate (C)*/
    float D1; /**< functional form of the segregation rate (D)*/
    float D2; /**< functional form of the segregation rate (D)*/
    float D3; /**< functional form of the segregation rate (D)*/
    float * totalVolumes; /**< Total population volume changes time series (based on dt)*/
    int lenTotalV; /**< Length of the totalVolumes array*/

    //######################## GSL ODE model stuff ################
    
    //double modelInitialParams[8]; /**< Initial params of the model*/
    //double modelInitialSpecies[15]; /**< Initial species of the model*/
    //float modelGeneLocations[3]; /**< Location on the chromosome (between 0.0 and 100.0, expressed in %) of the gene in the model*/
    //int modelGeneLRPos[3]; /**< Position of the gene on either the left or right hand side of the chromosome; where 0 -> right, 1 -> left*/
    //int modelGeneParamsLocations[3]; /**< Location of the parameteris in the model (note not species) that describes the copy number of a gene in particular */
    //double modelInitialParams[8]; /**< Initial params of the model*/
    //double modelInitialSpecies[15]; /**< Initial species of the model*/
    //float modelGeneLocations[3]; /**< Location on the chromosome (between 0.0 and 100.0, expressed in %) of the gene in the model*/
    //int modelGeneLRPos[3]; /**< Position of the gene on either the left or right hand side of the chromosome; where 0 -> right, 1 -> left*/
    //int modelGeneParamsLocations[3]; /**< Location of the parameteris in the model (note not species) that describes the copy number of a gene in particular */
    //double modelInitialParams[NUM_MODELPARAMS]; /**< Initial params of the model*/
    //double modelInitialSpecies[NUM_MODELSPECIES]; /**< Initial species of the model*/
    //float modelGeneLocations[NUM_MODELGENES]; /**< Location on the chromosome (between 0.0 and 100.0, expressed in %) of the gene in the model*/
    //int modelGeneLRPos[NUM_MODELGENES]; /**< Position of the gene on either the left or right hand side of the chromosome; where 0 -> right, 1 -> left*/
    //int modelGeneParamsLocations[NUM_MODELGENES]; /**< Location of the parameteris in the model (note not species) that describes the copy number of a gene in particular */
    
    //####################### libsbmlsim ##############    
    SBMLDocument_t* sbml_document;
    Model_t* sbml_model;
    int sbml_simulation_method;
    boolean sbml_use_lazy_method;
 
    /* num of variables whose quantity is not a constant */
    unsigned int sbml_num_of_all_var_species;
    unsigned int sbml_num_of_all_var_parameters;
    unsigned int sbml_num_of_all_var_compartments;
    unsigned int sbml_num_of_all_var_species_reference;

    /* num of SBase objects */
    unsigned int sbml_num_of_species;
    unsigned int sbml_num_of_parameters;
    unsigned int sbml_num_of_compartments;
    unsigned int sbml_num_of_reactions;
    unsigned int sbml_num_of_rules;
    unsigned int sbml_num_of_events;
    unsigned int sbml_num_of_initialAssignments;
 
    /* num of variables (which is NOT changed by assignment nor algebraic rule) */
    unsigned int sbml_num_of_var_species;
    unsigned int sbml_num_of_var_parameters;
    unsigned int sbml_num_of_var_compartments;
    unsigned int sbml_num_of_var_species_reference;
    
} CellPopulation;

int growCells(CellPopulation * cellPopulation,
                float dt,
                float volumeAdd,
                bool isDrugTreat);

void cellDebug(Cell * cellArray, int index);

#endif
