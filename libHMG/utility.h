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
* @file utility.h
*
* Contains a random assortment of functions that are used throughout the program, such as random number generator and calculation functions.
*
* @version 1.0
* @author Melchior du Lac
*
*/

#ifndef UTILITY_H
#define UTILITY_H

int seedRandom();

double getRandomZeroOne();

double normalDistRandn(double mu, double sigma);

int globalParameters(double mu, double * globalParam);

double * general_parameters_geussCD(double mu, double C1, double C2, double C3, double D1, double D2, double D3);

double * cell_parameters(double tau, double C, double D, double a, double Vi);

double * retRepTimers(double C, double D, double tau, double a, int maxRepTimers);

int findPowerFloor(int input);

int findPowerCeil(int input);

float adderDivVol(float tau);

#endif
