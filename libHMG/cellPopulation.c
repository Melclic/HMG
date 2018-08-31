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
* @file cellPopulation.c
*
* Cell Population main file that grows the population
*
* @version 1.0
* @author Melchior du Lac
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>

#include "cellPopulation.h"
#include "utility.h"

#include "../libsbmlsim/libsbmlsim.h"

/*
int nestedFindBlockedRep(CellPopulation * cellPopulation, int i)
{
    int y, j, z;
    for(y=0; y<cellPopulation->cellArray[i].numChrom; y++)
    {
        for(j=0; j<cellPopulation->cellArray[i].chromArray[y].potentialOriC; j++)
        {
            for(z=0; z<pow(2,y); z++)
            {
                if(cellPopulation->cellArray[i].chromArray[y].replicationTimers[j][z][0]==-2.0 &&
                cellPopulation->cellArray[i].chromArray[y].replicationTimers[j][z][1]==-2.0)
                {
                    return i;
                }
            }

        }
    }
    return -1;  
}
*/

/**
* @brief Grow cell population 
*
* Loop through all the cellArray and grow the cells using ueither injection or exponential growth. Deal with the different return flags from growCell. 0 -> no errors, 1 -> freeze the cell since it has reached the maximal number of chromosome permitted, 2 -> delete the cell because of an error, 3 -> the cell has does not have any more chromosomes. Deal with the case that the maximal number of cells has been reached by calling the restrictNumCells functions to randomly delete cells in the model. 
*
* @param cellPopulation CellPopulation object
* @param dt Time step
* @param volumeAdd In the case growth is by injection this paramter is >0, else the cells are grown assuming exponential (Malthusian) growth using model->cellPopulation->tau
* @param isDrugTreat Boolean type integer that determines if the model is growing under drug treatment conditions
* @return Error handling integer
*/
int growCells(CellPopulation * cellPopulation, float dt, float volumeAdd, bool isDrugTreat)
{
    int stopFlag = 0;
    int i, y;
    //Need to keep track of two different dynamically allocated arrays
    //for(i=0; i<cellPopulation->maxCells; i++)
//clock_t begin = clock();
    for(i=0; i<cellPopulation->indexArray; i++)
    {
        if(cellPopulation->cellArray[i].isDead==false && cellPopulation->cellArray[i].isNewlyRep==false && cellPopulation->cellArray[i].isFrozen==false)
        {
            //clock_t begin = clock();
            stopFlag = growCell(cellPopulation->cellArray,
                                i,
                                cellPopulation->freeIndex,
                                cellPopulation->tau,
                                cellPopulation->cNoise,
                                cellPopulation->dNoise,
                                cellPopulation->Vi,
                                cellPopulation->Vi_plasmid,
                                cellPopulation->ViNoise,
                                cellPopulation->VaNoise,
                                cellPopulation->chanceInit,
                                cellPopulation->divRatio,
                                cellPopulation->divNoise,
                                cellPopulation->partRatio,
                                cellPopulation->partNoise,
                                cellPopulation->chromDeg, 
                                cellPopulation->repForkDeg,
                                cellPopulation->numCells,
                                cellPopulation->C1,
                                cellPopulation->C2,
                                cellPopulation->C3,
                                cellPopulation->D1,
                                cellPopulation->D2,
                                cellPopulation->D3,
                                dt,
                                volumeAdd,
                                cellPopulation->modelGeneParamsLocations,
                                cellPopulation->modelGeneLocations,
                                cellPopulation->modelGeneLRPos,
                                isDrugTreat);
            //printf("\tgrowCell: %Lf\n", (long float)(clock()-begin));
            if(stopFlag==1)
            {
                //printf("WARNING (growCells): Setting cell %d back to factory settings\n", i);
                printf("WARNING (growCells): Freezing cell %d \n", i);
                cellPopulation->cellArray[i].isFrozen = true;
                cellPopulation->numFrozenCells += 1;
                //constructCell(cellPopulation->cellArray, i);
                //cellPopulation->numCells--;
            }
            else if(stopFlag==2)
            {
                //printf("WARNING (growCells): Anucleate cell (2) %d \n", i);
                constructCell(cellPopulation->cellArray, i);
                cellPopulation->numCells--;
                cellPopulation->numAnucleateCells += 1;
                cellPopulation->freeIndex = i;
                //return 1;
            }
            else if(stopFlag==3)
            {
                //printf("WARNING (growCells): Anucleate cell (3) %d \n", i);
                cellPopulation->numAnucleateCells += 1;
                //cellPopulation->numCells--;
            }
            //find a new index space if the space has been taken by a newly divided cell
            //avoids having to loop throught the whole array everytime you update a cell
            if(cellPopulation->cellArray[cellPopulation->freeIndex].isDead==false)
            {
                //Although we miss one time step using this method it is better than the rest
                if(cellPopulation->indexArray>=cellPopulation->maxCells-5 || cellPopulation->numCells>=cellPopulation->maxCells-5)
                {
                    printf("WARNING (growCells): maxCells (%d) or indexArray (%d) limit has been reached (%d)\n", cellPopulation->maxCells, cellPopulation->indexArray, cellPopulation->numCells);
                    return 2;
                }

                if(cellPopulation->freeIndex==cellPopulation->indexArray)
                {
                    cellPopulation->indexArray += 1;
                }
                cellPopulation->numCells += 1;
                cellPopulation->freeIndex = cellPopulation->indexArray;
            }
        }
        //if there is a cell that is detected to be dead (because it has been randonmly killed), then use that posisition for newly grown cells
        else if(cellPopulation->cellArray[i].isDead==true)
        {
            if(cellPopulation->cellArray[cellPopulation->freeIndex].isDead==false)
            {
                cellPopulation->freeIndex = i;
            }
            else if(cellPopulation->freeIndex>i)
            {
                cellPopulation->freeIndex = i;
            }           
        }
        //this makes sure that if the cells has been newly replicated, it is not updated by growCell, this makes sure it is reset so the next time step it can be grown
        else if(cellPopulation->cellArray[i].isNewlyRep==true)
        {
            cellPopulation->cellArray[i].isNewlyRep = false;
        }
    }
    //DEBUG
    /*
    if((long float)(clock()-begin)/CLOCKS_PER_SEC>0.1)
    {
        cellDebug(cellPopulation->cellArray, i);        
        printf("\tgrowCells: %Lf\n", (long float)(clock()-begin)/CLOCKS_PER_SEC);
    }
    */
    //stop the simulation if there are more than 25% of the cells that are isFrozen
    //if((cellPopulation->numFrozenCells>=cellPopulation->numCells-10 && cellPopulation->numCells>50) || cellPopulation->numFrozenCells==cellPopulation->numCells)
    if((cellPopulation->numFrozenCells>=((int)cellPopulation->numCells/4) && cellPopulation->numCells>50) || (cellPopulation->numFrozenCells==cellPopulation->numCells && cellPopulation->numCells>50))
    {
        printf("ERROR (growCells): Number of isFrozen cells exceeds threshold (%d)\n", cellPopulation->numFrozenCells);
        return 1;
    }
    if(cellPopulation->numCells<=0)
    {
        printf("ERROR (growCells): There are no cells!\n");
        return 1;
    }
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
int initCellPopulation(CellPopulation *cellPopulation)
{
    int i;
    for(i=0; i<cellPopulation->maxCells; i++)
    {
        constructCell(cellPopulation->cellArray, i);
    }

    for(i=0; i<cellPopulation->numCells; i++)
    {
        //if the input C time is input in minutes
        //TODO: add error handling from the initialiseCell
        if(cellPopulation->C2==-1.0)
        {
            initialiseCell(cellPopulation->cellArray,
                    i,
                    cellPopulation->tau,
                    cellPopulation->C1,
                    6646.0*cellPopulation->C1/4639221.0, //This is based on the size of the ColE1 in base pairs against the size of the chromosome in base pairs for a bacterial chromosome
                    cellPopulation->cNoise,
                    cellPopulation->D1,
                    cellPopulation->dNoise,
                    cellPopulation->Vi,
                    cellPopulation->Vi_plasmid,
                    cellPopulation->ViNoise,
                    cellPopulation->Vi/2.0, //TODO: why is this /2 ????
                    cellPopulation->VaNoise,
                    cellPopulation->modelInitialSpecies,
                    cellPopulation->modelInitialParams,
                    0.0);
        }
        //if the input C time is in its functional form
        else
        {
            initialiseCell(cellPopulation->cellArray,
                    i,
                    cellPopulation->tau,
                    cellPopulation->C1*(1.0+(cellPopulation->C2*exp(-cellPopulation->C3/(cellPopulation->tau/60.0)))),
                    //TODO change this so that it is not specific to ColE1 and is valid for minichromosomes and other plasmids
                    6646.0*(cellPopulation->C1*(1.0+(cellPopulation->C2*exp(-cellPopulation->C3/(cellPopulation->tau/60.0)))))/4639221.0, //This is based on the size of the ColE1 in base pairs against the size of the chromosome in base pairs for a bacterial chromosome
                    cellPopulation->cNoise,
                    cellPopulation->D1*(1.0+(cellPopulation->D2*exp(-cellPopulation->D3/(cellPopulation->tau/60.0)))),
                    cellPopulation->dNoise,
                    cellPopulation->Vi,
                    cellPopulation->Vi_plasmid,
                    cellPopulation->ViNoise,
                    cellPopulation->Vi/2.0,
                    cellPopulation->VaNoise,
                    cellPopulation->modelInitialSpecies,
                    cellPopulation->modelInitialParams,
                    0.0);
        }
    }
    return 0;
}

//TODO: add a clean cells functions to make cleaning of the model easier
int setCellPopulation(cellPopulation* cellPopulation 
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
                    float D3)
{
    //TODO: kind of stupid, should be able to tell him how many cells i want to inoculate with random number
    cellPopulation->numCells = 1;
    cellPopulation->freeIndex = 1;
    cellPopulation->indexArray = 1;
    cellPopulation->tau = tau;
    cellPopulation->C1 = C1;
    cellPopulation->C2 = C2;
    cellPopulation->C3 = C3;
    cellPopulation->D1 = D1;
    cellPopulation->D2 = D2;
    cellPopulation->D3 = D3;
    cellPopulation->dNoise = dNoise;
    cellPopulation->cNoise = cNoise;
    cellPopulation->Vi = Vi;
    cellPopulation->Vi_plasmid = Vi_plasmid;
    cellPopulation->ViNoise = ViNoise;
    cellPopulation->VaNoise = VaNoise;
    cellPopulation->chanceInit = chanceInit;
    cellPopulation->divNoise = divNoise;
    cellPopulation->divRatio = divRatio;
    cellPopulation->partNoise = partNoise;
    cellPopulation->partRatio = partRatio;   
    cellPopulation->chromDeg = chromDeg;
    cellPopulation->repForkDeg = repForkDeg; 
    cellPopulation->numFrozenCells = 0;  
    return 0;
}

int setSBMLmodel(CellPopulation * cellPopulation, const char* sbml_file)
{ 
    cellPopulation->method = MTHD_RUNGE_KUTTA;
    cellPopulation->boolean use_lazy_method = false;
    cellPopulation->print_interval = 1;
    cellPopulation->print_amount = 0;
    cellPopulation->atol = ABSOLUTE_ERROR_TOLERANCE;
    cellPopulation->rtol = RELATIVE_ERROR_TOLERANCE;
    cellPopulation->facmax = DEFAULT_FACMAX;
	
	cellPopulation->sbml_mem = allocated_memory_create();
	cellPopulation->sbml_cp_AST = copied_AST_create();

    cellPopulation->sbml_document = readSBMLFromFile(sbml_file);
    if(cellPopulation->sbml_document==NULL)
    {
        //return create_myResult_with_errorCode(Unknown);
        return 1;
    }
    unsigned int err_num = SBMLDocument_getNumErrors(cellPopulation->sbml_document);
    printf("err_num: %d\n", err_num);
    if(err_num>0)
    {
        const XMLError_t *err = (const XMLError_t *)SBMLDocument_getError(cellPopulation->sbml_document, 0);
        if(XMLError_isError(err) || XMLError_isFatal(err))
        {
            XMLErrorCode_t errcode = XMLError_getErrorId(err);
            switch(errcode)
            {
                case XMLFileUnreadable:
                    //myResult rtn = create_myResult_with_errorCode(FileNotFound);
                    break;
                case XMLFileUnwritable:
                case XMLFileOperationError:
                case XMLNetworkAccessError:
                    //myResult rtn = create_myResult_with_errorCode(SBMLOperationFailed);
                    break;
                case InternalXMLParserError:
                case UnrecognizedXMLParserCode:
                case XMLTranscoderError:
                    //myResult rtn = create_myResult_with_errorCode(InternalParserError);
                    break;
                case XMLOutOfMemory:
                    //myResult rtn = create_myResult_with_errorCode(OutOfMemory);
                    break;
                case XMLUnknownError:
                    //myResult rtn = create_myResult_with_errorCode(Unknown);
                    break;
                default:
                    //myResult rtn = create_myResult_with_errorCode(InvalidSBML);
                    break;
            } 
            SBMLDocument_free(cellPopulation->sbml_document);
            return 1;
        }
    }

    cellPopulation->sbml_model = SBMLDocument_getModel(sbml_document);
    cellPopulation->sbml_species = (mySpecies**)malloc(sizeof(mySpecies*) * Model_getNumSpecies(cellPopulation->sbml_model));
    cellPopulation->sbml_parameters = (myParameter**)malloc(sizeof(myParameter*) * Model_getNumParameters(cellPopulation->sbml_model));
    cellPopulation->sbml_compartments = (myCompartment**)malloc(sizeof(mySpecies*) * Model_getNumCompartments(cellPopulation->sbml_model));
    cellPopulation->sbml_reactions = (myReaction**)malloc(sizeof(myReaction*) * Model_getNumReactions(cellPopulation->sbml_model));
    cellPopulation->sbml_rules = (myRule**)malloc(sizeof(myRule*) * Model_getNumRules(cellPopulation->sbml_model));
    cellPopulation->sbml_events = (myEvent**)malloc(sizeof(myEvent*) * Model_getNumEvents(cellPopulation->sbml_model));
    cellPopulation->sbml_initialAssignement = (myInitialAssignment**)malloc(sizeof(myInitialAssignment*) * Model_getNumInitialAssignments(cellPopulation->sbml_model));
	cellPopulation->sbml_algebraicEquations = NULL;

    /*
    char *method_name;
    switch(method) {
    case MTHD_RUNGE_KUTTA: //  Runge-Kutta 
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
    case MTHD_BACKWARD_EULER: //  Backward-Euler 
      method_name = MTHD_NAME_BACKWARD_EULER;
      break;
    case MTHD_CRANK_NICOLSON: //  Crank-Nicolson (Adams-Moulton 2) 
      method_name = MTHD_NAME_CRANK_NICOLSON;
      break;
    case MTHD_ADAMS_MOULTON_3: //  Adams-Moulton 3 
      method_name = MTHD_NAME_ADAMS_MOULTON_3;
      break;
    case MTHD_ADAMS_MOULTON_4: //  Adams-Moulton 4 
      method_name = MTHD_NAME_ADAMS_MOULTON_4;
      break;
    case MTHD_BACKWARD_DIFFERENCE_2: //  Backward-Difference 2 
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_2;
      break;
    case MTHD_BACKWARD_DIFFERENCE_3: //  Backward-Difference 3 
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_3;
      break;
    case MTHD_BACKWARD_DIFFERENCE_4: //  Backward-Difference 4 
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_4;
      break;
    case MTHD_EULER: //  Euler (Adams-Bashforth) 
      method_name = MTHD_NAME_EULER;
      break;
    case MTHD_ADAMS_BASHFORTH_2: //  Adams-Bashforth 2 
      method_name = MTHD_NAME_ADAMS_BASHFORTH_2;
      break;
    case MTHD_ADAMS_BASHFORTH_3: //  Adams-Bashforth 3 
      method_name = MTHD_NAME_ADAMS_BASHFORTH_3;
      break;
    case MTHD_ADAMS_BASHFORTH_4: //  Adams-Bashforth 4 
      method_name = MTHD_NAME_ADAMS_BASHFORTH_4;
      break;
    // Variable Step Size 
    case MTHD_RUNGE_KUTTA_FEHLBERG_5: //  Runge-Kutta-Fehlberg 
      method_name = MTHD_NAME_RUNGE_KUTTA_FEHLBERG_5;
      is_variable_step = true;
      break;
    case MTHD_CASH_KARP: //  Cash-Karp
      method_name = MTHD_NAME_CASH_KARP;
      is_variable_step = true;
      break;
    default:
      cellPopulation->method = MTHD_RUNGE_KUTTA;
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
    }
    */
    cellPopulation->sbml_order = cellPopulation->method / 10;
    cellPopulation->sbml_is_explicit = cellPopulation->method % 10;
    //question remains as to the position of this parameter? Should it be here
    //  --> might be better to put this inside the cell so we can modify it
	int tmp_time = 0; //useless parameter that is used to set the sbml object 
    create_mySBML_objects(false, 
							cellPopulation->sbml_model, 
							cellPopulation->sbml_species, 
							cellPopulation->sbml_parameters, 
							cellPopulation->sbml_compartments, 
							cellPopulation->sbml_reactions,
							cellPopulation->sbml_rules, 
							cellPopulation->sbml_events,
      						cellPopulation->sbml_initialAssignement, 
							cellPopulation->sbml_algebraicEquations,
							cellPopulation->sbml_timeVarAssign,
      						dt, 
							dt,
							&tmp_time, 
							cellPopulation->sbml_mem,
							cellPopulation->cp_AST,
							1);
}


int cleanCellPopulation(CellPopulation *cellPopulation)
{
    //fclose(f); close the file
    for(int i=0; i<cellPopulation->maxCells; i++)
    {
        //gsl_odeiv2_driver_free(cellPopulation->cellArray[i].driver);
        constructCell(cellPopulation->cellArray, i);
    }   
    free(cellPopulation->cellArray);
    cellPopulation->cellArray = NULL;

    cellPopulation->maxCells = 0;    
    cellPopulation->numCells = 0;    
    cellPopulation->indexArray = 0;
    cellPopulation->freeIndex = 0;   
    cellPopulation->Vi = 0.0;    
    cellPopulation->Vi_plasmid = 0.0;    
    cellPopulation->cNoise = 0.0;    
    cellPopulation->dNoise = 0.0;    
    cellPopulation->ViNoise = 0.0;   
    cellPopulation->VaNoise = 0.0;   
    cellPopulation->chanceInit = 0.0;    
    cellPopulation->divNoise = 0.0;  
    cellPopulation->divRatio = 0.0;
    cellPopulation->partRatio = 0.0;
    cellPopulation->partNoise = 0.0;
    cellPopulation->chromDeg = 0.0;
    cellPopulation->repForkDeg = 0.0;    
    cellPopulation->numFrozenCells = 0;  
    cellPopulation->numAnucleateCells = 0;

    cellPopulation->C1 = 0.0;
    cellPopulation->C2 = 0.0;
    cellPopulation->C3 = 0.0;
    cellPopulation->D1 = 0.0;
    cellPopulation->D2 = 0.0;
    cellPopulation->D3 = 0.0;

    free(cellPopulation->totalVolumes);
    cellPopulation->lenTotalV = 0;
}

