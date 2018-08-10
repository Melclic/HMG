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
* @brief EXPERIMENTAL: size at division calculator based on the adder model of growth
*
* EXPERIMENTAL: Based on the adder model, we use the doubling rate of the cell and calculate the volume of the cell at division and use it as a rule to dictate how the cells are growing
*
* @param tau Doubling rate
* @return Volume at division
**/
float adderDivVol(float tau)
{
    //float r = (0.41*pow(2.0, 0.36*(60.0/tau)))/2.0;
    float w = (0.41*pow(2.0, 0.36*(60.0/tau)));
    float divL = 21.74*exp(-0.09895*tau)+4.001; //<- based on our fit of the authors data
    //return (4.0/3.0*(M_PI*pow(r, 2.0))+(2.0*M_PI*r*(divL-2.0*r)));
    return ((M_PI/4.0)*(pow(w, 2.0)))*(divL-(w/3.0));
}

//int const POWER_ARRAY[9] = {0,2,4,8,16,32,64,128,256}; //UNUSED
/*
int findPreviousPowerFloor(int input)
{
    return (int)log2(2,floor(log(input)/log(2)));
}
*/

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

/**
* @brief UNUSED: Exponential distribution population parameters based of growth rate (UNSUSED)
*
* UNUSED: Return parameters that are shared with evevry member of a population (based on specific growth rate). Note that this should be used for EXPONENTIAL growth only.
*
* @param mu Intantaneous growth rate
* @param globalParam Empty array (size 3) of the return parameters (C, D, Vi)
* @return Error handling integer
*/
int globalParameters(double mu, double * globalParam)
{
    double tau = log(2.0)/(mu);
    //Michelsen
    double C = 231.4*exp(-319.7*(1.0/tau))+40.0;
    double D = 80.89*exp(-121.3*(1.0/tau))+20.0;
    //Keasling
    //double C = 43.2*(1.0+4.86*exp(-4.5/(1.0/(60.0/tau))));
    //double D = (18.5+4.5)*(1.0+4.86*exp(-4.5/(1.0/(60.0/tau))));
    //double Vi = 3.148;
    //Global Fit
    //double C = 114.3*exp(-209.6*(1.0/tau))+40.0;  
    //double D = 240.8*exp(-371.4*(1.0/tau))+20.0;
    //NewFit
    //double D = 47.546413109294441;
    //double C = 80.73321259640186+71.36907278660817*exp(-234.32901046825492*mu);
    //double Vi = 2.1763333020923747;

    //Combined fit
    //double C = 73.22*(1.0+4.86*exp(-4.5/(1.0/(60.0/tau))));
    //double D = 39.03*(1.0+4.86*exp(-4.5/(1.0/(60.0/tau))));   
    //double Vi = 2.1763333020923747;

    //user Vi fit
    //double C = 69.48*exp(-1223.0*mu)+83.34;
    //double D = 56.1*exp(-269.8*mu)+54.61;
    //double Vi = 4.424*exp(-113.5*mu);

    //not every dt update!
    //double C = -2636.0*mu+96.77;
    //double C = 41.86*exp(-164.6*mu)+61.01;
    //double D = 54.98;
    //double Vi = 2.948;
        
    //GOOD ONE
    //double C = -1534*mu+97.16;
    //double D = 59.79;
    //double D = 55.0;
    double Vi = 0.9;

    //last attempt
    //double C = 153.15*exp(-(log(2.0)/0.0008)*mu)+(-1534*mu+97.16);
    //double D = 79.8*exp(-(log(2.0)/0.0009)*mu)+59.79;
    //double C = 166.65*exp(-(log(2.0)/0.0015)*mu)+86.5;
    //double D = 79.8*exp(-(log(2.0)/0.0015)*mu)+59.79;
    //double Vi = 3.148;
    
    //lsdldlsd
    //double C = 280*exp(-(log(2)/0.0015)*mu)+77.0;
    //double C = 190*exp(-(log(2)/0.0015)*mu)+60.0;
    //double D = 140.0*exp(-(log(2.0)/0.0013)*mu)+59.79;
    //double D = 100.0*exp(-(log(2.0)/0.0012)*mu)+45;
    //double D = 59.79;
    //double Vi = 3.148;
    //double Vi = 5.0;

    //##################### Volume at Initiation ######################
    //Volkmer fit
    //double Vi = 610.7-8315.1*mu; <--- note that the authors report this to be volume, but simulation evidence sugests that it is volume
    //double Vi = 610.7;
    //###### Would et al. fit
    //double k = floor((C+D)/tau);
    //double ai = 1.0+k-(C+D)/tau;  
    //double Vav = -0.53*pow((tau/60.0),2.0)+2.41*(tau/60.0)+1.62; //<--- original Volkmer et al.
    //double Vav = 153.4*mu+1.334; //<--- fit from Volkmer
    //double Vi = (Vav*pow(2.0,ai))/(2.0*log(2.0)*pow(2.0,k)); //<-- Wold et al.

    ///double Vi = 1.1257414881562995;

    globalParam[0] = C;
    globalParam[1] = D;
    globalParam[2] = Vi;
    return 1;
}

/**
* @brief UNUSED: Calculate the individual cells parameters such as DNA content and replication fork timers based assuming Matlhusian growth (UNUSED)
*
* UNUSED: Given the replication time (C), the segregation time (D), the critical volume (Vi), the doubling rate of the cell (tau) and the age of the cell (a); initiate the cell by calculating the number of replication forks, the resulting DNA contant, and segregation timer. Based on the calculation from Abner et al.
*
* @param tau Doubling rate
* @param C Replication time
* @param D Segregation time
* @param a Cell age (0<=a<=1.0)
* @param Vi Critical volume (default==0.9)
* @return Cellular parameters (size 7): 0: Volume, 1: G0, 2: Ga, 3: Segregation time, 4: Number of oriC, 5: number of terC (or chromosomes).
*/
//TODO: initiate the population assuming exponential distribution holds, instead of starting with a single cell. Could also be a way to test if spread on the population method is valid.
//Remember that k and l relate to the number of termination and oriC only at ai
//Assumes exponential growth
double * cell_parameters(double tau, double C, double D, double a, double Vi)
{
    //double volume = pow(2.0,(C+D+a)/tau-1.0);
        
    //my own work WARNING Vass at time a==0
    //double volume = (30170.0*((log(2.0)/tau))+736.0)/(1.0+(log(2.0)/tau)*0.244);

    double segregation_timer = 0.0;
    //DNA
    double l = floor(D/tau);
    double k = floor((C+D)/tau);

    //double volume = pow(2.0,(C+D+a)/tau-1.0)*(Vi/pow(2.0,k));
    double volume = pow(2.0,(C+D+(a*tau))/tau-1.0);
    //WARNING: ONLY IF WE HAVE real Vi
    //what is the relative cell age at Mi
    double ai = tau*(k+1.0)-C-D; //<-- abner et al
    //double ai = tau*(k+1.0)-C-D; //<-- wold et al
    //double Vu = Vi/pow(2.0,k);
    //volume = volume*Vu;
        

    int terC = pow(2,l);
    int oriC = pow(2,k);

    double G0 = ((pow(2.0,k)-1.0)*(tau/C)) + ((pow(2.0,k)*(fmod((C+D),tau))/C)) -((pow(2.0,l)-1.0)*(tau/C)) - (pow(2.0,l)*((fmod(D,tau))/C));

    double * ain = (double *)malloc(sizeof(double)*((int)k+1));
    ain[(int)k] = tau*(k+1.0)-C-D;
    int i;
    for(i=(int)k-1; i>=0; i--)
    {
        ain[i] = ain[i+1]-tau;
    }

    //Kind of stupid but required to see if there are twice as many terC
    for(i=0; i<((int)k+1); i++)
    {
        if(ain[i]<0.0 && ain[i]+C>0.0)
        {
            if(C-(C+ain[i])+a>=C)
            {
                terC = terC*2;
            }
        }
        else if(ain[i]>=0.0 && a>ain[i])
        {
            if(a-ain[i]>=C)
            {
                terC = terC*2;
            }
        }
    }
    //loop throught all the possible replication forks to calculate how much DNA to add
    double sum_An = 0.0;
    for(i=0; i<((int)k+1); i++)
    {
        double dn = 0.0;
        double bn = 0.0;

        if(ain[i]<0.0)
        {
            dn = 0.0;
        }
        else
        {
            dn = ain[i];
        }

        if(a<ain[i])
        {
            bn = 0.0;
        }
        else if(ain[i]+C>a && a>=ain[i])
        {
            bn = C;
        }
        else if(a>=ain[i]+C)
        {
            if(tau<=(C+D))
            {
                if(ain[i]+C==0.0)
                {
                    bn = 0.0;
                }
                else
                {
                    bn = (a-dn)*C/(ain[i]+C);
                }
            }
            else
            {
                bn = a-dn;
            }
        }
        else
        {
            printf("ERROR: Stupid calculation problem\n");
        }
        if(bn>0.0)
        {
            sum_An += pow(2.0,i)*((a-dn)/bn);
        }
    }
    double Ga = G0+sum_An;
    double at = (1.0-D/tau)*tau;

    //calculate the initiation time and termination time and set eberythign correctly
    if(a>at)
    {
        segregation_timer = a-at;
        if(segregation_timer>=D)
        {
            segregation_timer = 0.0;
            volume= volume/2.0;
            terC = terC/2;
            Ga = Ga/2.0;
            oriC = oriC/2;
        }
    }
    if(a>ai)
    {
        oriC = oriC*2;
    }

    static double ret[6];

    ret[0] = volume;
    ret[1] = G0;
    ret[2] = Ga;
    ret[3] = segregation_timer;
    ret[4] = (double)oriC;
    //just warning why there is a segregation error once if the growth rate return more than MAX_CHROM
    if(terC>8)
    {
        printf("ERROR: Cannot have more than 8 chromosomes\n");
    }
    ret[5] = (double)terC;

    return ret;
}

/**
* @brief UNUSED: Calculate individual cells number of replication timers and their times
*
* UNUSED. Given the replication time (C), the segregation time (D), the critical volume (Vi), the doubling rate of the cell (tau) and the maximal possible replication timers  (this should be determined based on the observable fastest growth rate of the species), we calculate the number of replication forks and their times. 
*
* @param C Replication time
* @param D Segregation time
* @param tau Doubling rate
* @param a Cell age (0<=a<=1.0)
* @param maxRepTimers Ceiling number of replication forks possible
* @return Array representing the times of the replication forks
*/
double * retRepTimers(double C, double D, double tau, double a, int maxRepTimers)
{
    //double l = floor(D/tau);
    double k = floor((C+D)/tau);

    //create the list of all possible initiation times in the life of a cell
    double * replication_timers = malloc(sizeof(double)*maxRepTimers);
    int i;
    for(i=0; i<maxRepTimers; i++)
    {
        replication_timers[i] = 0.0;
    }

    double * ain = (double *)malloc(sizeof(double)*((int)k+1));
    ain[(int)k] = tau*(k+1.0)-C-D;
    for(i=(int)k-1; i>=0; i--)
    {
        ain[i] = ain[i+1]-tau;
    }

    //calculate the inherited timers when initiating the cells at random time a.
    int count = 0;
    //int loc = 0;
    //int y;
    for(i=0; i<((int)k+1); i++)
    {
        double tmpRepTimer = 0.0;
        if(ain[i]<0.0 && ain[i]+C>0.0)
        {
            tmpRepTimer = C-(C+ain[i])+a;
            if(tmpRepTimer>=C)
            { 
               tmpRepTimer = 0.0;
            }
        }
        else if(ain[i]>=0.0 && a>ain[i])
        {
            tmpRepTimer = a-ain[i];
            if(tmpRepTimer>=C)
            {
                tmpRepTimer = 0.0;
            }
        }
        //make the replication fork array
        if(tmpRepTimer!=0.0)
        {
            replication_timers[count] = tmpRepTimer;
        }                
        count++;
        /*
        if(tmpRepTimer!=0.0)
        {
            if(count!=0)
            {
                for(y=0; y<(int)pow(2,count); y++)
                {
                    replication_timers[loc+y] = tmpRepTimer;
                }
                loc = (int)pow(2,count)+1;
                count += 1;
            }
            else
            {
                replication_timers[0] = tmpRepTimer;
                loc = 1;
                count += 1;
            }
        }
        */
    }
    free(ain);
    return replication_timers;
    //return ain;
}
