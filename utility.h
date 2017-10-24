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
