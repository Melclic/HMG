/*
 *   libHMG - Individual based model a bacterial population
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
* @file model.h
*
* Model header file with all the methods that initiate and run the model.
*
* @version 1.0
* @author Melchior du Lac
*
*/

#ifndef MODEL_H
#define MODEL_H

#include <stdint.h>
#include <stdio.h>

#include "cellPopulation.h"

/**
* @brief Structure containing the model description
*
* Structure of the model that contains the time step, time of the simulation, the volumetric changes in the population
* and the structure of the cellPopulation description.
*/
typedef struct MODEL
{
	float dt; /**< Time step */
	float t; /**< Current time of the simulation*/
	int stop; /**< Boolean controlling if to terminate the simulation*/
	int nodeSize; /**< MPI parameter determining the number of nodes that are initiated*/
	CellPopulation * cellPopulation; /**< Struct of population*/
} Model;

//####### GETTERS ##########

int getNumCells(Model * model);

int getRealNumCells(Model * model);

double getTotalVolume(Model * model);

//int getDNAContent(Model * model, float * DNAContent);

int getPlasmidDNAContent(Model * model, float * DNAContent);

//######## ESSENTIAL #########

int inoculateModel(Model * model);

int cleanModel(Model * model);

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
                float dt);

int runInjection(Model * model,
		float maximalExecTime,
		int restrictCellNumber,
		double * totalVolumes,
		int lenTotalV);

int runExpo(Model * model, float maximalExecTime, int targetCellCount);

int runExpoCD(Model * model, float maximalExecTime, int targetCellCount, float C, float D);

int runDrugTreatment(Model * model, float maximalExecTime, int targetCellCount, float drugNoise);

Model * initModel(int maxCells);

int randomRestrictNumCells(Model * model, int targetNumCells);

int getNumExposedGenes(Model * model, int ** numExposed, float percLoc, int LR);

float meanPopD(Model * model);

float meanPopC(Model * model);

float getMeanVa(Model * model);

int getDistVa(Model * model, float * Va);

int getDistGa(Model * model, float * Ga);

int getDistAge(Model * model, float * age);

double getMeanModelSpecies(Model * model, int speciesNum);

double getSingleModelSpecies(Model * model, int speciesNum, int cellNum);

int oneTimeStep(Model * model);

int getDistTau(Model * model, float * tauContent);

int getDistC(Model * model, float * CContent);

int getDistD(Model * model, float * DContent);

int getDistPrev_Vd(Model * model, float * dist_Vd);

int getDistPrev_Vb(Model * model, float * dist_Vb);

int getDistPrev_a(Model * model, float * dist_a);

float getStdVa(Model * model);

float getMeanTau(Model * model);

float getStdTau(Model * model);

#endif
