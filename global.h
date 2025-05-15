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


using namespace std;



#include "random.h"



double  *lbound;
double  *ubound;

int dimension;


// *********************************************************************************

//******** Parameters in random number *********************************************
int     seed    = 177;
long    rnd_uni_init;      
     
int nichsize;


int		max_gen = 500,    //  the maximal number of generations
		max_run = 1,      //  the maximal number of runs
		pops    = 10,    //  the population size
		maxfes;             //  the number of function evluations
//**********************************************************************************
int     nvar;      //  the number of variables
int ID;
double dist_all;
int fes;
int levelgen;
double bestof;
double worstof;
//double similar;
int stepevol;
int *dimvar;
int ndsize;
double Factor;
double averge;
double maxobj;
double minobj;
vector<int> record;
ofstream  fr("p1.txt");
int gp;
#endif