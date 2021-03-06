/*
 * libHMG - Individual based model of a bacterial population
 * Copyright (C) 2017  Melchior du Lac
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 */


/**
* @file pyConnect.h
*
* Header file for pyConnect
*
* @version 1.0
* @author Melchior du Lac
*/

#ifndef PYCONNECT_H
#define PYCONNECT_H

#include <stdint.h>
#include <stdio.h>

int pyCall_initModel(float cNoise,
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
                        int maxCells,
                        float C1,
                        float C2,
                        float C3,
                        float D1,
                        float D2,
                        float D3,
                        float tau,
                        double * modelInitialParams,
                        double * modelInitialSpecies,
                        float * modelGeneLocations,
                        int * modelGeneParamsLocations,
                        int * modelGeneLRPos,
                        float dt);

int pyCall_runInjection(float maximalExecTime,
                        int restrictCellNumber
                        double * totalVolumes,
                        int lenTotalV);

int pyCall_runExpo(float maximalExecTime, int targetCellCount, float tau)

float pyCall_meanPopC();

float pyCall_meanPopD();

int pyCall_restrictNumCells(int targetNumCells);

float pyCall_getNumExposedGenes(int ** numExposed, float percLoc, int LR);

float pyCall_getTotalVolume();

int pyCall_getRealNumCells();

int pyCall_getDistGa(float * DNAContent);

int pyCall_getDistVa(float * VContent);

int pyCall_cleanModel();

float pyCall_getMeanV();

double pyCall_getMeanModelSpecies(int speciesNum);

double pyCall_getSingleModelSpecies(int speciesNum, int cellNum);

int pyCall_getDistAge(float * age);

int pyCall_oneTimeStep();

int pyCall_getNumCells();

float pyCall_getCellVa(int cellNum);

float pyCall_getCellGa(int cellNum);

int pyCall_getDistTau(float * tauContent);

float pyCall_getMeanVa();

float pyCall_getStdVa(float Vav);

float pyCall_getMeanTau();

float pyCall_getStdTau(float tauAv);

int pyCall_getDistPrev_a(float * dist_a);

int pyCall_getDistPrev_Vb(float * dist_Vb);

int pyCall_getDistPrev_Vd(float * dist_Vd);

int pyCall_runDrugTreatment(float maximalExecTime, int targetCellCount, float drugNoise);

int pyCall_getDistC(float * CContent);

int pyCall_getDistD(float * DContent);

int pyCall_getNumAnucleateCells();

#endif

