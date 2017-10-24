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
