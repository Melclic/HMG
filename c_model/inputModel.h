/*
 *   libHMG - Individual based model of the a bacterial population
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
* @file input_models.h
*
* Header file of input_models
*
* @version 1.0
* @author Melchior du Lac
*
*/

#ifndef INPUT_MODELS_H
#define INPUT_MODELS_H

//todo make this the way to input parameters into the model function
/*
typedef struct PARAMS
{
	
} Params;
*/

//extern int NUM_MODELPARAMS;
//extern int NUM_MODELSPECIES;
//extern int NUM_MODELGENES;
//WARNING: cannot have 0 sized arrays

//#define NUM_MODELPARAMS 0
//#define NUM_MODELSPECIES 0
//#define NUM_MODELGENES 0

//#define NUM_MODELPARAMS 8
//#define NUM_MODELSPECIES 15
//#define NUM_MODELGENES 3

int model_ori_repressilator(double t, const double y[], double f[], void *params);

int model_bennett_repressilator(double t, const double y[], double f[], void *params);

#endif
