# A New Explicit Penalty Method for Evolutionary Multimodal Optimization

1. Project Overview
This project implements a multimodal optimization algorithm designed to handle multiple optima (global and local) in complex fitness landscapes. It is built upon the CEC2013 benchmark.

2. File Structure

```
   EPM/
   ├── nichtree.cpp       # Main entry point of the program.
   ├── CEC2013.h/cpp      # Interface and implementation for CEC2013 benchmark functions.
   ├── cfunction.h/cpp    # Base and composite function definitions.
   ├── global.h           # Global variables and parameters.
   ├── individual.h       # Definition of the individual class and evolutionary operations.
   ├── niches.h           # Definition of the nich (niche) class for maintaining subpopulations.
   ├── ode.h              # Main evolutionary algorithm logic and population management.
   └── README.md 
``` 

3. Key Functions and Headers
   1. nichtree.cpp
   Main function: Runs the optimization algorithm for multiple test functions (ID 1–20) and evaluates the number of global optima  found.
   Outputs: Results are written to me.txt and result.txt.
  
   2. CEC2013.h/cpp
    Purpose: Provides the interface to the CEC2013 benchmark functions.
    Key Functions:
        evaluate(): Evaluates a solution vector.
        get_lbound(), get_ubound(): Return variable bounds.
        get_fitness_goptima(), get_no_goptima(): Return global optimum fitness and count.
   
   3. cfunction.h/cpp
    Purpose: Defines base functions (e.g., Sphere, Rastrigin) and composite functions (CF1–CF4).
    Key Functions:
        FSphere, FRastrigin, FGriewank, etc.: Basic benchmark functions.
        CF1::evaluate(), CF2::evaluate(), etc.: Composite function evaluations.
   
   4. global.h
    Purpose: Contains global parameters and variables.
    Key Variables:
        lbound, ubound: Variable bounds.
        dimension: Problem dimension.
        maxfes: Maximum function evaluations.
        pops: Population size.

   5. individual.h
    Purpose: Defines the individual class representing a candidate solution.
    Key Methods:
        evaluation(): Evaluates the individual using a CEC2013 function.
        diff(), diff_r2(), etc.: Differential evolution mutation strategies.
        repair(): Ensures variables stay within bounds.

    6. niches.h
     Purpose: Defines the nich class for managing niches (subpopulations).
     Key Methods:
        nevol(): Evolves the niche.
        find_lb(): Finds bounds for the niche.
        select_fitness(): Selects individuals based on fitness and distance.

    7. ode.h
     Purpose: Implements the main evolutionary loop and population management.
     Key Methods:
        run(): Main optimization loop.
        sortfit(), cal_distance(): Sort and compute crowding distance.
        select_elite(): Selects elite individuals.
        iteration(): Generates new offspring.

4. How to Compile and Run
   1. Compile all .cpp files together (e.g., with g++).
   2. Run the executable. Results will be saved in me.txt and result.txt.

5. Output
   1. The algorithm reports the number of global optima found for each function under different accuracy levels (1E-1, 1E-2, ..., 1E-5).


6. License
    This project is for academic research purposes only.