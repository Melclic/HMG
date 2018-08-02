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
* @file utility.c
*
* Contains a random assortment of functions that are used throughout the program, such as random number generators.
*
* @version 1.0
* @author Melchior du Lac
*
*/

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "utility.h"
#define M_PI 3.14159265358979323846


/**
* @brief Calculate the next upper power of 2 value
*
* Calculates the next integer power of 2 ceiling given an integer input
*
* @param inuput Input value
* @return Output value
*/
int findPowerCeil(int input)
{
	return (int)pow(2,ceil(log(input)/log(2))) ;
}

/**
* @brief Calculate the previous upper power of 2 value
*
* Calculates the previous integer power of 2 ceiling given an integer input
*
* @param inuput Input value
* @return Output value
*/
int findPowerFloor(int input)
{
	return (int)pow(2,floor(log(input)/log(2))) ;
}

/**
* @brief Initiate the random number generator.
*
* Random function of C needs to be intitialised with a seed. i.e. the same seed will return the same sequence of values. This is why this is seeded with the time at execution.
*
* @return Error handling integer
*/
int seedRandom()
{
	srand(time(NULL));
	return 0;
}

/**
* @brief Random number generator
*
* Generate a random number between (0>=1)
*
* @return Output value
*/
double getRandomZeroOne()
{
	return (double)(rand())/(double)(RAND_MAX);
}

/**
* @brief Random number generator based on a Gaussian distribution
*
* Given a Gaussian distribution with mean (mu) and standard deviation (sigma), pool a random number. Taken from http://phoxis.org/2013/05/04/generating-random-numbers-from-normal-distribution-in-c/
*
* @param mu Gaussian mean
* @param sigma Gaussian standard deviation
* @return Output value
*/
double normalDistRandn(double mu, double sigma)
{
	double U1, U2, W, mult;
	static double X1, X2;
	static int call = 0;

	if (call == 1)
	{
		call = !call;
		return (mu + sigma * (double) X2);
	}
	do
	{
		U1 = -1 + ((double) rand () / RAND_MAX) * 2;
		U2 = -1 + ((double) rand () / RAND_MAX) * 2;
		W = pow (U1, 2) + pow (U2, 2);
	}
	while (W >= 1 || W == 0);

	mult = sqrt ((-2 * log (W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;

	return (mu + sigma * (double) X1);
}
