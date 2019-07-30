#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//######################### PUBLIC FUNCTIONS #############################

//Remember that k and l relate to the number of termination and oriC only at ai
//Assumes exponential growth
//tau: doubling time in minutes
//C: replication time in minutes
//D: segregation time in minutes
//a: age of the cell between 0.0 and 1.0
//results: array of size 6 to return the results
int cell_parameters(double tau, double C, double D, double a, double *results)
{
        double segregation_timer = 0.0;
        double l = floor(D/tau);
        double k = floor((C+D)/tau);   
        //calculate the age of the cell in minutes according to the doubling rate
        double age = a*tau;
        //double volume = pow(2.0,(C+D+a)/tau-1.0)*(Vi/pow(2.0,k));
        //double volume = pow(2.0,(C+D+(a*tau))/tau-1.0);
        
        double volume = pow(2.0,(C+D+(age))/tau-1.0);
        
        double ai = tau*(k+1.0)-C-D; //<-- abner et al
        //double ai = tau*(k+1.0)-C-D; //<-- wold et al

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
                        //if(C-(C+ain[i])+a>=C)
                        if(C-(C+ain[i])+age>=C)
                        {
                                terC = terC*2;
                        }
                }
                //else if(ain[i]>=0.0 && a>ain[i])
                else if(ain[i]>=0.0 && age>ain[i])
                {
                        //if(a-ain[i]>=C)
                        if(age-ain[i]>=C)
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

                //if(a<ain[i])
                if(age<ain[i])
                {
                        bn = 0.0;
                }
                //else if(ain[i]+C>a && a>=ain[i])
                else if(ain[i]+C>age && age>=ain[i])
                {
                        bn = C;
                }
                //else if(a>=ain[i]+C)
                else if(age>=ain[i]+C)
                {
                        if(tau<=(C+D))
                        {
                                if(ain[i]+C==0.0)
                                {
                                        bn = 0.0;
                                }
                                else
                                {
                                        //bn = (a-dn)*C/(ain[i]+C);
                                        bn = (age-dn)*C/(ain[i]+C);
                                }
                        }
                        else
                        {
                                //bn = a-dn;
                                bn = age-dn;
                        }
                }
                else
                {
                        printf("ERROR: Stupid calculation problem\n");
                }
                if(bn>0.0)
                {
                        //sum_An += pow(2.0,i)*((a-dn)/bn);
                        sum_An += pow(2.0,i)*((age-dn)/bn);
                }
        }

        double Ga = G0+sum_An;
        double at = (1.0-D/tau)*tau;
        //calculate the initiation time and termination time and set eberythign correctly
        //if(a>at)
        if(age>at)
        {
                //segregation_timer = a-at;
                segregation_timer = age-at;
                if(segregation_timer>=D)
                {
                        segregation_timer = 0.0;
                        volume= volume/2.0;
                        terC = terC/2;
                        Ga = Ga/2.0;
                        oriC = oriC/2;
                }
        }
        //if(a>ai)
        if(age>ai)
        {
                oriC = oriC*2;
        }

        free(ain);

        results[0] = volume;
        results[1] = G0;
        results[2] = Ga;
        results[3] = segregation_timer;
        results[4] = (double)oriC;
        results[5] = (double)terC;

        return 0;
}

//to test the script
int main(int argc, const char* argv[])
{	
	double results[6];
	//double *results = cell_parameters(60.0, 40.0, 20.0, 0.5);
	int status = cell_parameters(60.0, 40.0, 20.0, 0.5, results);
	printf("results[0]: %f\n", results[0]);
	printf("results[1]: %f\n", results[1]);
	printf("results[2]: %f\n", results[2]);
	printf("results[3]: %f\n", results[3]);
	printf("results[4]: %f\n", results[4]);
	printf("results[5]: %f\n", results[5]);
}
