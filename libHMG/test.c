/*
 * <!--------------------------------------------------------------------------
 * This file is part of libSBMLSim.  Please visit
 * http://fun.bio.keio.ac.jp/software/libsbmlsim/ for more
 * information about libSBMLSim and its latest version.
 *
 * Copyright (C) 2011-2017 by the Keio University, Yokohama, Japan
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ---------------------------------------------------------------------- -->
*/
#include "libsbmlsim/libsbmlsim.h"
/*
int main(int argc, char *argv[]) {
  myResult *rtn;

  double sim_time = 25;
  double dt = 0.01;
  int print_interval = 10;
  int print_amount = 1;
  int method = MTHD_RUNGE_KUTTA;
  // int method = MTHD_RUNGE_KUTTA_FEHLBERG_5; 
  boolean use_lazy_method = false;

  if (argc < 2) {
    printf("Input SBML file is not specified.\n  Usage: %s sbml.xml\n", argv[0]);
    exit(1);
  }
  rtn = simulateSBMLFromFile(argv[1], sim_time, dt, print_interval, print_amount, method, use_lazy_method);
  if (rtn == NULL) {
    printf("Returned result is NULL\n");
  } else if (myResult_isError(rtn)) {
    printf("ERROR: %s\n", myResult_getErrorMessage(rtn));
  } else {
    //  print_result(rtn); 
    // write_csv(rtn, "cresult.csv"); 
    write_result(rtn, "test.dat");
  }
  if (rtn != NULL)
    free_myResult(rtn);
  return 0;
}
*/


SBMLSIM_EXPORT myResult* main(int argc, char *argv[]) 
{
    printf("Starting reading the file %s \n", argv[1]);


  if (argc < 2) {
    printf("Input SBML file is not specified.\n  Usage: %s sbml.xml\n", argv[0]);
    exit(1);
  }
  
	

  double sim_time = 0.01;
  double dt = 0.01;
  int method = MTHD_RUNGE_KUTTA;
  boolean use_lazy_method = false;
  int print_interval = 1;
  int print_amount = 0;

  myResult *result, *rtn = NULL;
  SBMLDocument_t* d;
  Model_t* m;
  unsigned int err_num;
  double atol = 0.0;
  double rtol = 0.0;
  double facmax = 0.0;
  d = readSBMLFromFile(argv[1]);
  if (d == NULL)
    return create_myResult_with_errorCode(Unknown);
  err_num = SBMLDocument_getNumErrors(d);
  printf("err_num: %d\n", err_num);
  if (err_num > 0) {
    const XMLError_t *err = (const XMLError_t *)SBMLDocument_getError(d, 0);
    if (XMLError_isError(err) || XMLError_isFatal(err)) {
      XMLErrorCode_t errcode = XMLError_getErrorId(err);
      switch (errcode) {
        case XMLFileUnreadable:
          rtn = create_myResult_with_errorCode(FileNotFound);
          break;
        case XMLFileUnwritable:
        case XMLFileOperationError:
        case XMLNetworkAccessError:
          rtn = create_myResult_with_errorCode(SBMLOperationFailed);
          break;
        case InternalXMLParserError:
        case UnrecognizedXMLParserCode:
        case XMLTranscoderError:
          rtn = create_myResult_with_errorCode(InternalParserError);
          break;
        case XMLOutOfMemory:
          rtn = create_myResult_with_errorCode(OutOfMemory);
          break;
        case XMLUnknownError:
          rtn = create_myResult_with_errorCode(Unknown);
          break;
        default:
          rtn = create_myResult_with_errorCode(InvalidSBML);
          break;
      }
      SBMLDocument_free(d);
      return rtn;
    }
  }
    printf("Done reading the file \n");
    printf("Starting reading the document \n");

  m = SBMLDocument_getModel(d);

  //rtn = simulateSBMLModel(m, sim_time, dt, print_interval, print_amount, method, use_lazy_method, atol, rtol, facmax);// BELOW
  
    printf("Here 1\n");

    /*
  if (rtn == NULL)
  {
    printf("rtn detected as NULL \n");
    free_myResult(result);
    return rtn;
  }
  */

    printf("Done reading the document \n");
    printf("Set the parameters for the simulation \n");

  double time = 0;
  int order = 0;
  int is_explicit = 0;
  char *method_name;
  unsigned int num_of_species;
  unsigned int num_of_parameters;
  unsigned int num_of_compartments;
  unsigned int num_of_reactions;
  unsigned int num_of_rules;
  unsigned int num_of_events;
  unsigned int num_of_initialAssignments;
  mySpecies **mySp;
  myParameter **myParam;
  myCompartment **myComp;
  myReaction **myRe;
  myRule **myRu;
  myEvent **myEv;
  myInitialAssignment **myInitAssign;
  // prepare myAlgebraicEquations 
  myAlgebraicEquations *myAlgEq = NULL;
  // prepare timeVariantAssignments 
  timeVariantAssignments *timeVarAssign = NULL;
  // prepare return value 
  // Variables for bifurcation analysis 
  char buf1[256], buf2[256], buf3[256], buf4[256], buf5[256], buf6[256];
  boolean use_bifurcation_analysis = false; // Bifurcation is turned off for this release 
  boolean state_variable_exists = false;
  boolean bifurcation_parameter_exists = false;
  boolean bif_param_is_local = false;
  boolean is_variable_step = false;
  char* sta_var_id = NULL;
  char* bif_param_id = NULL;

  double bif_param_min = 0;
  double bif_param_max = 0;
  double bif_param_stepsize = 0.01;

  // to exclude the transition state 
  double transition_time = 0;

  char* tmp;
  unsigned int i, a, b;

  // for variable stepsize 
  int err_zero_flag = 0;

  allocated_memory *mem;
  copied_AST *cp_AST;

  mem = allocated_memory_create();
  cp_AST = copied_AST_create();

  // Check atol, rtol and facmax, whether it is set to 0.0 
  if (atol == 0.0) {
    atol = ABSOLUTE_ERROR_TOLERANCE;
  }
  if (rtol == 0.0) {
    rtol = RELATIVE_ERROR_TOLERANCE;
  }
  if (facmax == 0.0) {
    facmax = DEFAULT_FACMAX;
  }

  // prepare mySpecies 
  num_of_species = Model_getNumSpecies(m);
  // mySpecies *mySp[num_of_species]; 
  mySp = (mySpecies**)malloc(sizeof(mySpecies*) * num_of_species);
  // prepare myParameters 
  num_of_parameters = Model_getNumParameters(m);
  // myParameter *myParam[num_of_parameters]; 
  myParam = (myParameter**)malloc(sizeof(myParameter*) * num_of_parameters);
  // prepare myCompartments 
  num_of_compartments = Model_getNumCompartments(m);
  // myCompartment *myComp[num_of_compartments]; 
  myComp = (myCompartment**)malloc(sizeof(mySpecies*) * num_of_compartments);
  // prepare myReactions 
  num_of_reactions = Model_getNumReactions(m);
  // myReaction *myRe[num_of_reactions]; 
  myRe = (myReaction**)malloc(sizeof(myReaction*) * num_of_reactions);
  // prepare myRules 
  num_of_rules = Model_getNumRules(m);
  // myRule *myRu[num_of_rules]; 
  myRu = (myRule**)malloc(sizeof(myRule*) * num_of_rules);
  // prepare myEvents 
  num_of_events = Model_getNumEvents(m);
  // myEvent *myEv[num_of_events]; 
  myEv = (myEvent**)malloc(sizeof(myEvent*) * num_of_events);
  // prepare myInitial Assignments 
  num_of_initialAssignments = Model_getNumInitialAssignments(m);
  // myInitialAssignment *myInitAssign[num_of_initialAssignments]; 
  myInitAssign = (myInitialAssignment**)malloc(sizeof(myInitialAssignment*) * num_of_initialAssignments);

  switch(method) {
    case MTHD_RUNGE_KUTTA: //  Runge-Kutta 
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
    case MTHD_BACKWARD_EULER: //  Backward-Euler 
      method_name = MTHD_NAME_BACKWARD_EULER;
      break;
    case MTHD_CRANK_NICOLSON: //  Crank-Nicolson (Adams-Moulton 2) 
      method_name = MTHD_NAME_CRANK_NICOLSON;
      break;
    case MTHD_ADAMS_MOULTON_3: //  Adams-Moulton 3 
      method_name = MTHD_NAME_ADAMS_MOULTON_3;
      break;
    case MTHD_ADAMS_MOULTON_4: //  Adams-Moulton 4 
      method_name = MTHD_NAME_ADAMS_MOULTON_4;
      break;
    case MTHD_BACKWARD_DIFFERENCE_2: //  Backward-Difference 2 
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_2;
      break;
    case MTHD_BACKWARD_DIFFERENCE_3: //  Backward-Difference 3 
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_3;
      break;
    case MTHD_BACKWARD_DIFFERENCE_4: //  Backward-Difference 4 
      method_name = MTHD_NAME_BACKWARD_DIFFERENCE_4;
      break;
    case MTHD_EULER: //  Euler (Adams-Bashforth) 
      method_name = MTHD_NAME_EULER;
      break;
    case MTHD_ADAMS_BASHFORTH_2: //  Adams-Bashforth 2 
      method_name = MTHD_NAME_ADAMS_BASHFORTH_2;
      break;
    case MTHD_ADAMS_BASHFORTH_3: //  Adams-Bashforth 3 
      method_name = MTHD_NAME_ADAMS_BASHFORTH_3;
      break;
    case MTHD_ADAMS_BASHFORTH_4: //  Adams-Bashforth 4 
      method_name = MTHD_NAME_ADAMS_BASHFORTH_4;
      break;
    // Variable Step Size 
    case MTHD_RUNGE_KUTTA_FEHLBERG_5: //  Runge-Kutta-Fehlberg 
      method_name = MTHD_NAME_RUNGE_KUTTA_FEHLBERG_5;
      is_variable_step = true;
      break;
    case MTHD_CASH_KARP: //  Cash-Karp
      method_name = MTHD_NAME_CASH_KARP;
      is_variable_step = true;
      break;
    default:
      method = MTHD_RUNGE_KUTTA;
      method_name = MTHD_NAME_RUNGE_KUTTA;
      break;
  }
  order = method / 10;
  is_explicit = method % 10;
    printf("Done setting parameters for the simulation \n");

    printf("Create the SBML objects \n");
  // create myObjects 
  // free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv,
  //       myInitAssign, myAlgEq, timeVarAssign, mem, cp_AST); 
  create_mySBML_objects(is_variable_step, m, mySp, myParam, myComp, myRe, myRu, myEv,
      myInitAssign, &myAlgEq, &timeVarAssign,
      sim_time, dt, &time, mem, cp_AST, print_interval);

  // create myResult 
  if (is_variable_step) {
    result = create_myResultf(m, mySp, myParam, myComp, sim_time, dt);
  } else {
    result = create_myResult(m, mySp, myParam, myComp, sim_time, dt, print_interval);
  }
    printf("Done creating the SBML objects \n");


  //here we would like to call these once and loop ourselves
  //TEST --> set sim_time==dt


    printf("Starting the loop \n");
  for(int g=0; g<100; g++)
  {
     //   printf("Loop: %d \n", g);
      // simulation 
      if (is_variable_step) {
        if (is_explicit) {
          rtn = simulate_explicitf(m, result, mySp, myParam, myComp, myRe, myRu, myEv,
              myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval,
              &time, order, print_amount, mem, atol, rtol, facmax, cp_AST,
              &err_zero_flag);
        }
      } else {  // Fixed step size 
        if (is_explicit) {
          rtn = simulate_explicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv,
              myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval,
              &time, order, print_amount, mem);
        }else{
          rtn = simulate_implicit(m, result, mySp, myParam, myComp, myRe, myRu, myEv,
              myInitAssign, myAlgEq, timeVarAssign, sim_time, dt, print_interval,
              &time, order, use_lazy_method, print_amount, mem);
        }
      }
      print_result(rtn);
      time = 0;
  }



  // free 
  if (use_bifurcation_analysis == 0 || (use_bifurcation_analysis == 1 && bif_param_is_local == false)) {
    free_mySBML_objects(m, mySp, myParam, myComp, myRe, myRu, myEv,
        myInitAssign, myAlgEq, timeVarAssign, mem, cp_AST);
  }
 
  if (rtn == NULL)
    rtn = create_myResult_with_errorCode(SimulationFailed);
  SBMLDocument_free(d);
  return rtn;
}
