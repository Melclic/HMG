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
int growCells(CellPopulation * cellPopulation, float dt, float volumeAdd, bool isDrugTreat, int noderank)
{
        int stopFlag = 0;
        //Need to keep track of two different dynamically allocated arrays
        //for(i=0; i<cellPopulation->maxCells; i++)
	//clock_t begin = clock();
        //for(i=0; i<cellPopulation->indexArray; i++)
	
	for(int i=cellPopulation->nodeStartEnd[noderank-1][0]; i<cellPopulation->nodeStartEnd[noderank-1][1]; i++)
        {
                if(cellPopulation->cellArray[i].isDead==false && 
			cellPopulation->cellArray[i].isNewlyRep==false && 
			cellPopulation->cellArray[i].isFrozen==false &&
			cellPopulation->cellArray[i].isSimAge==false)
                {
			//clock_t begin = clock();
                        stopFlag = growCell(cellPopulation->cellArray,
                                i,
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
				//return 1;
			}
			else if(stopFlag==3)
			{
                                //printf("WARNING (growCells): Anucleate cell (3) %d \n", i);
				cellPopulation->numAnucleateCells += 1;
				//cellPopulation->numCells--;
			}
			else if(stopFlag==4)
			{
				printf("WARNING (growCells): Did not receive valid new cell location\n");
				return 1;
			}

			/*
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
			*/
			cellPopulation->cellArray[index].simAge += dt;
                }
                //if there is a cell that is detected to be dead (because it has been randonmly killed), then use that posisition for newly grown cells
                else if(cellPopulation->cellArray[i].isDead==true)
                {
			if(cellPopulation->cellArray[cellPopulation->nodeFreeCell[noderank-1]].isDead==false)
			{
                        	cellPopulation->nodeFreeCell[noderank-1] = i;
			}
			else if(cellPopulation->nodeFreeCell[noderank-1]>i)
			{
				cellPopulation->nodeFreeCell[noderank-1] = i;
			}			
                }
                //this makes sure that if the cells has been newly replicated, it is not updated by growCell, this makes sure it is reset so the next time step it can be grown
                else if(cellPopulation->cellArray[i].isNewlyRep==true)
                {
                        cellPopulation->cellArray[i].isNewlyRep = false;
                }
		else if(cellPopulation->cellArray[i].isSimAge==true && cellPopulation->cellArray[i].isCounted==false)
		{
			cellPopulation->nodeFinishedCells[noderank-1] += 1;
			cellPopulation->cellArray[i].isCounted==true;	
		}
		//only for MPI
		if(cellPopulation->cellArray[index].isSimAge==false && 
			cellPopulation->cellArray[index].simAge>=cellPopulation->simAge)
		{
			cellPopulation->cellArray[index].isSimAge = true;
		}
	}

	//DEBUG
	/*
	if((long float)(clock()-begin)/CLOCKS_PER_SEC>0.1)
	{
		cellDllArray[indes].simAge += dt;
        if(cellArray[indes].simAge>=simAge)
        {
                cellArray[index].isSimAge = true;
        }
(cellPopulation->cellArray, i);		
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
