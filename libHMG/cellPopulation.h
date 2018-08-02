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
	Cell * cellArray; /**< Array of cells*/

	int * nodeNumCells; /**< Because the shared array is seperated between each node, we keep track of the number of alive cells in each node as to make sure new cells are added proportionally between each node*/ 
	int * nodeFreeCell; /**< Given the number of nodes set the position of the free cell*/
	int * nodeFinisedCells; /**< Number of finished cells per node*/
	int ** nodeStartEnd; /**< Gives the start and the end location of cellArray to each node*/
	float simAge; /**< Global simulated time of the dead and alive cells*/
	bool isSimAge;

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
	float tau; /**< cell doubling time (min) */
	float C1; /**< functional form of the replication rate (C)*/
	float C2; /**< functional form of the replication rate (C)*/
	float C3; /**< functional form of the replication rate (C)*/
	float D1; /**< functional form of the segregation rate (D)*/
	float D2; /**< functional form of the segregation rate (D)*/
	float D3; /**< functional form of the segregation rate (D)*/
        float * totalVolumes; /**< Total population volume changes time series (based on dt)*/
        int lenTotalV; /**< Length of the totalVolumes array*/
} CellPopulation;

int growCells(CellPopulation * cellPopulation,
                float dt,
		float volumeAdd,
		bool isDrugTreat);

void cellDebug(Cell * cellArray, int index);

#endif
