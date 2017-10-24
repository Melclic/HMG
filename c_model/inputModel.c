#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

//int NUM_MODELPARAMS = 8;
//int NUM_MODELSPECIES = 15;
//int NUM_MODELGENES = 3;
//#define NUM_MODELPARAMS 0
//#define NUM_MODELSPECIES 0
//#define NUM_MODELGENES 0

/**
* @brief Model function of ODE's that describes a system that would be simulated by the user. 
*
* We use the GSL ODE simulation functions. This requires the size of the params, and y to be known outside of the function. The simulation should progress one step size at a time
* TODO: In this situation the user must specify what are the parameters that need to be updated 
* @param t the time step
* @param y array of previous simulation results
* @param f array of the current simulation results
* @param params the constant parameters input into the model
*
**/
int model_ori_repressilator(double t, const double y[], double f[], void *params)
{
	(void)(t);
        double alpha_0 = ((double *)params)[0];
        double alpha = ((double *)params)[1];
        double n = ((double *)params)[2];
        double beta = ((double *)params)[3];

        f[0] = -y[0]+(alpha/(1.0+pow(y[5], n)))+alpha_0;
        f[1] = -y[1]+(alpha/(1.0+pow(y[3], n)))+alpha_0;
        f[2] = -y[2]+(alpha/(1.0+pow(y[4], n)))+alpha_0;
        f[3] = -beta*(y[3]-y[0]);
        f[4] = -beta*(y[4]-y[1]);
        f[5] = -beta*(y[5]-y[2]);

        return GSL_SUCCESS;
}

/*
* @brief Repressilator from Bennett et al. 2, 7 and 12 are the free or unbound promoter sites of the system. These are controlled by the user to be copy number controlled
*/
int model_bennett_repressilator(double t, const double y[], double f[], void *params)
{
        (void)(t);
        double gamma_p = ((double *)params)[0];
        double gamma_m = ((double *)params)[1];
        double kappa_plus = ((double *)params)[2];
        double kappa_minus = ((double *)params)[3];
        double k_plus = ((double *)params)[4];
        double k_minus = ((double *)params)[5];
        double alpha = ((double *)params)[6];
        double sigma = ((double *)params)[7];

        f[0] = -2.0*kappa_plus*pow(y[0], 2.0) + 2.0*kappa_minus*y[1] + sigma*y[4] - gamma_p*y[0];
        f[1] = kappa_plus*pow(y[0], 2.0) - kappa_minus*y[1] - k_plus*y[1]*y[7] + k_minus*y[8];
        f[2] = -k_plus*y[11]*y[2] + k_minus*y[3];
        f[3] = k_plus*y[11]*y[2] - k_minus*y[3];
        f[4] = alpha*y[2] - gamma_m*y[4];

        f[5] = -2.0*kappa_plus*pow(y[5], 2.0) + 2.0*kappa_minus*y[6] + sigma*y[9] - gamma_p*y[5];
        f[6] = kappa_plus*pow(y[5], 2.0) - kappa_minus*y[6] - k_plus*y[6]*y[12] + k_minus*y[13];
        f[7] = -k_plus*y[1]*y[7] + k_minus*y[8];
        f[8] = k_plus*y[1]*y[7] - k_minus*y[8];
        f[9] = alpha*y[7] - gamma_m*y[9];

        f[10] = -2.0*kappa_plus*pow(y[10], 2.0) + 2.0*kappa_minus*y[11] + sigma*y[14] - gamma_p*y[10];
        f[11] = kappa_plus*pow(y[10], 2.0) - kappa_minus*y[11] - k_plus*y[11]*y[2] + k_minus*y[3];
        f[12] = -k_plus*y[6]*y[12] + k_minus*y[13];
        f[13] = k_plus*y[6]*y[12] - k_minus*y[13];
        f[14] = alpha*y[12] - gamma_m*y[14];

        return GSL_SUCCESS;
}

/*
This is the bennett model for the repressilator with incoperation of the plasmid gene copy number. This requires the input of growth rate 
*/
int model_bennett_repressilator_plasmid(double t, const double y[], double f[], void *params)
{
        (void)(t);
        double gamma_p = ((double *)params)[0];
        double gamma_m = ((double *)params)[1];
        double kappa_plus = ((double *)params)[2];
        double kappa_minus = ((double *)params)[3];
        double k_plus = ((double *)params)[4];
        double k_minus = ((double *)params)[5];
        double alpha = ((double *)params)[6];
        double sigma = ((double *)params)[7];

        f[0] = -2.0*kappa_plus*pow(y[0], 2.0) + 2.0*kappa_minus*y[1] + sigma*y[4] - gamma_p*y[0];
        f[1] = kappa_plus*pow(y[0], 2.0) - kappa_minus*y[1] - k_plus*y[1]*y[7] + k_minus*y[8];
        f[2] = -k_plus*y[11]*y[2] + k_minus*y[3];
        f[3] = k_plus*y[11]*y[2] - k_minus*y[3];
        f[4] = alpha*y[2] - gamma_m*y[4];

        f[5] = -2.0*kappa_plus*pow(y[5], 2.0) + 2.0*kappa_minus*y[6] + sigma*y[9] - gamma_p*y[5];
        f[6] = kappa_plus*pow(y[5], 2.0) - kappa_minus*y[6] - k_plus*y[6]*y[12] + k_minus*y[13];
        f[7] = -k_plus*y[1]*y[7] + k_minus*y[8];
        f[8] = k_plus*y[1]*y[7] - k_minus*y[8];
        f[9] = alpha*y[7] - gamma_m*y[9];

        f[10] = -2.0*kappa_plus*pow(y[10], 2.0) + 2.0*kappa_minus*y[11] + sigma*y[14] - gamma_p*y[10];
        f[11] = kappa_plus*pow(y[10], 2.0) - kappa_minus*y[11] - k_plus*y[11]*y[2] + k_minus*y[3];
        f[12] = -k_plus*y[6]*y[12] + k_minus*y[13];
        f[13] = k_plus*y[6]*y[12] - k_minus*y[13];
        f[14] = alpha*y[12] - gamma_m*y[14];

        return GSL_SUCCESS;
}
