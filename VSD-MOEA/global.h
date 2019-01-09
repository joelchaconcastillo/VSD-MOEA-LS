#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include "random.h"


using namespace std;



//------------- Parameters in test instance ------------------

int     nvar=10,  nobj=2;                    //  the number of variables and objectives
int pops = 100; //population size

double  lowBound = 0,   uppBound = 1;   //  lower and upper bounds of variables
double  vlowBound[100] ,   vuppBound[100];   //  lower and upper bounds of variables

char    strTestInstance[256];
char    currentPATH[1500];
int param_l=5, param_k=5; // the distance and position parameters for the WFG problems..

long int  max_nfes = 10000; //The function evaluation criteria is prefered than generations..
//------------- Parameters in random number ------------------
int     seed    = 177; //Default seed...
long    rnd_uni_init;        

double Initial_lowest_distance_factor=0.2*sqrt(nvar), lowestDistanceFactor; 

//------------- Parameters in VSD-MOEA
double          scale[100];  

int		etax    = 20, 	etam    = 50;   // distribution indexes of crossover and mutation

double  realx=0.9,  realm = 1.0/nvar;    // crossover, mutation, selection probabilities
int run;

#endif
