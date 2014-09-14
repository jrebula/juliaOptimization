

#put all the snopt libs somewhere ld can find them: /usr/lib or whatever
# gcc -fpic -Wall -c snset.c
# gcc -fpic -Wall -c snopt.c
# gcc -fpic -Wall -c snerror.c
# gcc -fpic -Wall -c snget.c 
# gfortran -fpic -Wall -c sn_open.f
#gcc -shared -o libcsnopt.so snopt.o snerror.o snset.o sn_open.o snget.o -lsnopt7 -lsnblas -lsnprint7 -lgfortran


module Snopt

# how to load a constant:
#const gsl_brent = unsafe_load(cglobal((:gsl_min_fminimizer_brent,:libgsl), Ptr{Void}))

type SnoptProblem
    snoptProblemReference::Ptr{Void}
    SnoptProblem(s) = new(s)
end

function SnoptProblem(numVariables::Int,
                      numConstraints::Int,
                      numNonzeros::Int,
                      numJacobianNonzeros::Int,
                      numNonlinearConstraints::Int,
                      numNonlinearObjectiveVariables::Int,
                      numNonlinearJacobianVariables::Int,
                      linearObjectiveRow::Int,
                      outputFileName::String,
                      outputStream::String="stdout")
    
    if (length(Base.find_library(["libcsnopt"])) == 0)
        error("make sure libcsnopt.so is somewhere on LD_LIBRARY_PATH (you could put it in /usr/lib if you're root/lazy)!")
    end
    
    snoptProblemReference = c_malloc(10000)
    errorMessage = Ptr{Uint8}
    
    status = ccall((:snInit, "libcsnopt"),
                   Cint,
                   (Ptr{Void},
                    Cint, Cint, Cint, Cint,
                    Cint, Cint, Cint, Cint,
                    Ptr{Uint8}, Ptr{Uint8}),
                   snoptProblemReference,
                   numVariables,
                   numConstraints,
                   numNonzeros,
                   numJacobianNonzeros,
                   numNonlinearConstraints,
                   numNonlinearObjectiveVariables,
                   numNonlinearJacobianVariables,
                   linearObjectiveRow,
                   "manne.out",
                   "stdout")
    
    #const snoptOK = unsafe_load(cglobal((:SN_OK, "libcsnopt"), Ptr{Void}))
    snoptOK = 1;
    
    if status != snoptOK
        ccall((:snGetError, "libcsnopt"),
              Void,
              (Ptr{Void}, Ptr{Uint8}),
              manneReference, errorMessage)
        show(bytestring(errorMessage))
        
        #snGetError( &manne , &error_msg );
        #printf( "%s: error occurred (%d).\n%s", argv[0] , status, error_msg );
        #snDelete( &manne );
        #free(error_msg);
        #exit(EXIT_FAILURE);
    end

    ret = SnoptProblem(snoptProblemReference)
    ret
end







#=
macro wrapSnoptCall(functionToDefine, symbolToCall, outputType, inputTypes,
                    inputsOfNewFunction, allInputs...)
    :(function ($functionToDefine)($inputsOfNewFunction...)
        ccall(($symbolToCall, "libcsnopt"), $outputType, $inputTypes, $(allInputs...))
    end)
end

macroexpand(:(@wrapSnoptCall(getVariables, :snGetVariables, Cdouble, (Ptr{Void},), s::SnoptProblem, s)))

functionsToDefine = ((:getVariablesasdf, :snGetVariables, :Cdouble, :(Ptr{Void},), :(s::SnoptProblem,), (:s,)),
                     (:getVariablesasdf2, :snGetVariables, :Cdouble, :((Ptr{Void},)), :(s::SnoptProblem,), (:s,)));
for i = 1:length(functionsToDefine)
    thisFunctionParams = functionsToDefine[i];
    functionToCreate = thisFunctionParams[1];
    symbolToCall = thisFunctionParams[2]
    outputTpye = thisFunctionParams[3];
    inputTypes = thisFunctionParams[4];
    inputsOfNewFunction = thisFunctionParams[5];
    allInputs = thisFunctionParams[6];
    #show(functionToCreate)
    #@eval show()
    println("asdf")
    eval(:($functionToCreate = ($inputsOfNewFunction) -> ccall((:($symbolToCall), "libcsnopt"), $outputType, $inputTypes, $(allInputs...))))
    println("asdf")
end


f = x -> show(x...)
f(thisFunctionParams)

=#


macro snoptCall(symbolToCall, outputType, inputTypes, allInputs...)
    :(ccall(($symbolToCall, "libcsnopt"), $outputType, $inputTypes, $(allInputs...)))
end


macro snoptGetCall(outputType, symbolToCall)
    :(ccall(($symbolToCall, "libcsnopt"), $outputType, (Ptr{Void},), s.snoptProblemReference))
end

macro snoptGetCallToArray(outputType, symbolToCall, numElements)
    :(pointer_to_array(ccall(($symbolToCall, "libcsnopt"), $outputType, (Ptr{Void},), s.snoptProblemReference), $numElements))
end




macro snoptGetCallIndexed(outputType, symbolToCall)
    :(ccall(($symbolToCall, "libcsnopt"), $outputType, (Ptr{Void}, Cint), s.snoptProblemReference, i))
end



#macroexpand(:(@snoptCall(:snGetVariables, Cdouble, (Ptr{Void},), s)))

#function getNumVariables(s::SnoptProblem)
#    @snoptCall(:snGetNumVariables, Cdouble, (Ptr{Void},), s.snoptProblemReference)
#end

#function getVariables(s::SnoptProblem)
#    @snoptCall(:snGetVariables, Ptr{Cdouble}, (Ptr{Void},), s.snoptProblemReference)
#end



#macro defineSnoptGetFunction(functionToDefine, outputType, symbolToCall)
#    esc(:($functionToDefine = ((s::SnoptProblem),) -> ccall(($symbolToCall, "libcsnopt"), $outputType, (Ptr{Void},), s.snoptProblemReference)))
#end
#macroexpand(:(@defineSnoptGetFunction(:getMajorOptimalityTolerance, Cdouble, :snGetMajorOptimalityTolerance)))



getMajorOptimalityTolerance = (s::SnoptProblem) -> @snoptGetCall(Cdouble, :snGetMajorOptimalityTolerance)
getMajorFeasibilityTolerance = (s::SnoptProblem) -> @snoptGetCall(Cdouble, :snGetMajorFeasibilityTolerance)
getObjectiveAdd = (s::SnoptProblem) -> @snoptGetCall(Cdouble, :snGetObjectiveAdd)
getNumNonlinearConstraints = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetNumNonlinearConstraints)
getNumNonzeros = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetNumNonzeros)
getNumVariables = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetNumVariables)
getNumConstraints = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetNumConstraints)
getNumNonlinearObjectiveVariables = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetNumNonlinearObjectiveVariables)
getNumNonlinearJacobianVariables = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetNumNonlinearJacobianVariables)
getSolveStatus = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetSolveStatus)

getSolveMode = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetSolveMode)
getSuperBasicsLimit = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetSuperBasicsLimit)
getInfinity = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetInfinity)
getMajorIterationLimit = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetMajorIterationLimit)
getMinorIterationLimit = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetMinorIterationLimit)
getMajorOptimalityTolerance = (s::SnoptProblem) -> @snoptGetCall(Cdouble, :snGetMajorOptimalityTolerance)
getMajorFeasibilityTolerance = (s::SnoptProblem) -> @snoptGetCall(Cdouble, :snGetMajorFeasibilityTolerance)
getVerifyLevel = (s::SnoptProblem) -> @snoptGetCall(Cint, :snGetVerifyLevel)


getNonzeroValues = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cdouble}, :snGetNonzeroValues)
getNonzeroRowIndices = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cint}, :snGetNonzeroRowIndices)
getNonzeroColumnPointers = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cint}, :snGetNonzeroColumnPointers)
getVariables = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cdouble}, :snGetVariables)
getConstraints = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cdouble}, :snGetConstraints)
getConstraintUpperBounds = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cdouble}, :snGetConstraintUpperBounds)
getConstraintLowerBounds = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cdouble}, :snGetConstraintLowerBounds)
getVariableUpperBounds = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cdouble}, :snGetVariableUpperBounds)
getVariableLowerBounds = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cdouble}, :snGetVariableLowerBounds)
getMultipliers = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cdouble}, :snGetMultipliers)
getVariableBasisEligibilities = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cint}, :snGetVariableBasisEligibilities)
getConstraintBasisEligibilities = (s::SnoptProblem) -> @snoptGetCall(Ptr{Cint}, :snGetConstraintBasisEligibilities)


arrayFunctions = [
                  :getNonzeroValues,
                  :getNonzeroRowIndices,
                  :getNonzeroColumnPointers,
                  :getVariables,
                  :getConstraints,
                  :getConstraintUpperBounds,
                  :getConstraintLowerBounds,
                  :getVariableUpperBounds,
                  :getVariableLowerBounds,
                  :getMultipliers,
                  :getVariableBasisEligibilities,
                  :getConstraintBasisEligibilities]

# generate fooArray functions from each array function, this converts c pointers to arrays, given the number of elements in the array
for arrayFunction in arrayFunctions
    newFunctionName = arrayFunction + :Array
    #show(newFunctionName)
    expression = :($newFunctionName = (s::SnoptProblem, numberOfElements) -> pointer_to_array($arrayFunction(s), numberOfElements))
    show(expression)
    eval(expression)
    println();
end


#getNonzeroValue = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cdouble, :snGetNonzeroValue)
#getNonzeroRowIndex = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cint, :snGetNonzeroRowIndex)
#getNonzeroColumnPointer = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cint, :snGetNonzeroColumnPointer)
#getVariable = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cdouble, :snGetVariable)
#getConstraint = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cdouble, :snGetConstracint)
#getConstaintUpperBound = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cdouble, :snGetConstaintUpperBound)
#getConstaintLowerBound = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cdouble, :snGetConstaintLowerBound)
#getVariableUpperBound = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cdouble, :snGetVariableUpperBound)
#getVariableLowerBound = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cdouble, :snGetVariableLowerBound)
#getMultiplier = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cdouble, :snGetMultiplier)
#getVariableBasisEligibility = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cint, :snGetVariableBasisEligibility)
#getConstraintBasisEligibility = (s::SnoptProblem, i) -> @snoptGetCallIndexed(Cint, :snGetConstraintBasisEligibility)



function declareIncrementingConstants(names, startNumber=0)
    for i = 0 : (length(names)-1)
        name = names[i+1]
        eval(:(const $name = $(i+startNumber)))
    end
end




#/* boolean */
declareIncrementingConstants((:SN_FALSE, :SN_TRUE))

#/* error codes */
declareIncrementingConstants((
                              :SN_ERROR,
                              :SN_OK,
                              :SN_FILE_ERROR,
                              :SN_PROBLEM_UNINITIALIZED,
                              :SN_PROBLEM_UNSOLVED))
                              

#/* basis eligibility (hs) values */
declareIncrementingConstants((
                              :SN_BASIS_COLD_ELIGIBLE_0,
                              :SN_BASIS_COLD_ELIGIBLE_1,
                              :SN_BASIS_IGNORE,
                              :SN_BASIS_COLD_ELIGIBLE_3,
                              :SN_BASIS_LOWER_BOUND,
                              :SN_BASIS_UPPER_BOUND))

                              
#/* bits for set_flags */
const SN_BIT_MAJOR_FEASIBILITY_TOLERANCE = (1 << 0)
const SN_BIT_MINOR_FEASIBILITY_TOLERANCE = (1 << 1)
const SN_BIT_MAJOR_OPTIMALITY_TOLERANCE  = (1 << 2)
const SN_BIT_MINOR_OPTIMALITY_TOLERANCE  = (1 << 3)
const SN_BIT_SOLVED                      = (1 << 4)

#/* "Normal" result conditions for snSolve() */
declareIncrementingConstants((
                              :SN_SOLUTION_FOUND,
                              :SN_INFEASIBLE    ,
                              :SN_UNBOUNDED     ,
                              :SN_VIOLATION_LIMIT_EXCEEDED,
                              :SN_MINOR_ITERATION_LIMIT_EXCEEDED,
                              :SN_MAJOR_ITERATION_LIMIT_EXCEEDED,
                              :SN_ACCURACY_NOT_ACHIEVED         ,
                              :SN_SUPERBASICS_LIMIT_EXCEEDED        ,
                              :SN_POINT_CANNOT_BE_IMPROVED          ,
                              :SN_CANNOT_SATISFY_GENERAL_CONSTRAINTS,
                              :SN_SINGULAR_BASIS                    ), 1)


end


#= here's the c call the above is trying to reproduce


status =
    snInit
    ( &manne,
     variables,
     constraints,
     nonzeros,
     jacobian_nonzeros,
     nonlinear_constraints,
     nonlinear_objective_variables,
     nonlinear_jacobian_variables,
     linear_objective_row,
     "manne.out",
     "stdout" );

some more variable declarations, these happen at the top of the c file:
const *a_t   = NULL;
const *b_t   = NULL;

const *c_bl  = NULL;
const *v_bl  = NULL;

const *c_bu  = NULL;
const *v_bu  = NULL;

const *v_xs  = NULL;
const *c_xs  = NULL;

const *pi    = NULL;
const *Jcol  = NULL;

cint *v_hs  = NULL;
cint *c_hs  = NULL;

cint *indJ  = NULL;
cint *locJ  = NULL;

cint nT;


if ( status != SN_OK ) {
                        snGetError( &manne , &error_msg );
                        printf( "%s: error occurred (%d).\n%s", argv[0] , status, error_msg );

                        snDelete( &manne );
                        free(error_msg);

                        exit(EXIT_FAILURE);
                        }

    /* Populate the arrays that define problem manne. */
    /* Set some optional parameters.                  */

    manne_init( &manne );

    snSetProblemName        ( &manne, "   Manne            "       );
    snSetOptionString       ( &manne, "Maximize            "       );
    snSetOptionInteger      ( &manne, "Major iterations    ",  400 );
    snSetOptionInteger      ( &manne, "Minor iterations    ", 8000 );
    snSetOptionInteger      ( &manne, "SuperBasics limit   ",  400 );
    snSetOptionString       ( &manne, "Solution  no        "       );

    snSetObjectiveAdd       ( &manne, 0.0 );
    snSetUserFunction       ( &manne, manneUserfun );

    /* Solve the problem */

    status = snoptc ( &manne );

    if ( status == SN_OK )

        switch ( snGetSolveStatus( &manne ) ) {
                                               case SN_SOLUTION_FOUND:
                                               case SN_INFEASIBLE:
                                               case SN_UNBOUNDED:
                                               case SN_VIOLATION_LIMIT_EXCEEDED:
                                               case SN_MINOR_ITERATION_LIMIT_EXCEEDED:
                                               case SN_MAJOR_ITERATION_LIMIT_EXCEEDED:
                                               case SN_ACCURACY_NOT_ACHIEVED:
                                               case SN_SUPERBASICS_LIMIT_EXCEEDED:
                                               case SN_POINT_CANNOT_BE_IMPROVED:
                                               case SN_CANNOT_SATISFY_GENERAL_CONSTRAINTS:
                                               case SN_SINGULAR_BASIS:

                                               printf("normal status -- yippie\n");

                                               fprintf( stdout, "major iterations: %d\n",
                                                       snGetNumMinorIterations( &manne ) );
                                               fprintf( stdout, "minor iterations: %d\n",
                                                       snGetNumMajorIterations( &manne ) );
                                               fprintf( stdout, "objective function evaluations: %d\n",
                                                       snGetNumObjectiveFunctionEvals( &manne ) );
                                               fprintf( stdout, "constraint function evaluations: %d\n",
                                                       snGetNumConstraintFunctionEvals( &manne ) );
                                               fprintf( stdout, "superbasics: %d\n",
                                                       snGetNumSuperBasics( &manne ) );
                                               fprintf( stdout, "degenerate steps: %d\n",
                                                       snGetNumDegenerateSteps( &manne ) );
                                               fprintf( stdout, "fraction degenerate steps: %f\n",
                                                       (float)snGetNumDegenerateSteps( &manne ) /
                                                       (float)snGetNumMinorIterations( &manne ) );
                                               fprintf( stdout, "infeasibilities: %d\n",
                                                       snGetNumInfeasibilities( &manne ) );

                                               fprintf( stdout, "norm of scaled solution: %f\n",
                                                       snGetNormScaledSolution( &manne ) );
                                               fprintf( stdout, "norm of solution: %f\n",
                                                       snGetNormSolution( &manne ) );
                                               fprintf( stdout, "norm of scaled multipliers: %f\n",
                                                       snGetNormScaledMultipliers( &manne ) );
                                               fprintf( stdout, "norm of multipliers: %f\n",
                                                       snGetNormMultipliers( &manne ) );
                                               fprintf( stdout, "penalty: %f\n",
                                                       snGetPenalty( &manne ) );
                                               fprintf( stdout, "objective value: %f\n",
                                                       snGetObjectiveValue( &manne ) );
                                               fprintf( stdout, "linear objective value: %f\n",
                                                       snGetLinearObjectiveValue( &manne ) );
                                               fprintf( stdout, "nonlinear objective value: %f\n",
                                                       snGetNonlinearObjectiveValue( &manne ) );
                                               fprintf( stdout, "sum of infeasibilities: %f\n",
                                                       snGetSumInfeasibilities( &manne ) );

                                               fprintf( stdout, "maximum scaled primal infeasibility: %f\n",
                                                       snGetMaxScaledPrimalInfeasibility( &manne ) );
                                               fprintf( stdout, "index of maximum scaled primal infeasibility: %d\n",
                                                       snGetIMaxScaledPrimalInfeasibility( &manne ) );
                                               fprintf( stdout, "maximum scaled dual infeasibility: %f\n",
                                                       snGetMaxScaledDualInfeasibility( &manne ) );
                                               fprintf( stdout, "index of maximum scaled dual infeasibility: %d\n",
                                                       snGetIMaxScaledDualInfeasibility( &manne ) );
                                               fprintf( stdout, "maximum primal infeasibility: %f\n",
                                                       snGetMaxPrimalInfeasibility( &manne ) );
                                               fprintf( stdout, "index of maximum primal infeasibility: %d\n",
                                                       snGetIMaxPrimalInfeasibility( &manne ) );
                                               fprintf( stdout, "maximum dual infeasibility: %f\n",
                                                       snGetMaxDualInfeasibility( &manne ) );
                                               fprintf( stdout, "index of maximum dual infeasibility: %d\n",
                                                       snGetIMaxDualInfeasibility( &manne ) );

                                               }

    else {

          snGetError( &manne , &error_msg );
          printf( "%s: error occurred (%d).\n%s", argv[0] , status, error_msg );

          manne_delete();
          free(error_msg);
          snDelete( &manne );

          exit(EXIT_FAILURE);

          }

        snDelete( &manne );

        manne_delete();
        exit(EXIT_SUCCESS);

        return 0;
    end












    #void manne_init( snProblem* manne )
    #{
    #  cint     i, j;
    #  cint    *p_locJ;
    #  double   infty = 1e20, dummy = 0.1;

    A = ( xC_0 + xI_0 ) / (xK_0^B);
    gFac = (1.0 + growth)^(1-B);

    b_t = (double*)malloc( sizeof(double) * nT );
    a_t = (double*)malloc( sizeof(double) * nT );

    a_t[0] = A * gFac;
    b_t[0] = Beta;

    for ( i = 1 ; i < nT ; i++ ) {
                                  a_t[i] = a_t[i-1] * gFac;
                                  b_t[i] = b_t[i-1] * Beta;
                                  }

        b_t[nT-1] /= ( 1.0 - Beta );

        v_xs = snGetVariables  ( manne );
        v_hs = snGetVariableBasisEligibilities  ( manne );

        c_xs = snGetConstraints( manne );
        c_hs = snGetConstraintBasisEligibilities( manne );

        v_bu = snGetVariableUpperBounds  ( manne );
        v_bl = snGetVariableLowerBounds  ( manne );

        c_bu = snGetConstraintUpperBounds( manne );
        c_bl = snGetConstraintLowerBounds( manne );
        pi   = snGetMultipliers          ( manne );

        Jcol = snGetNonzeroValues        ( manne );
        indJ = snGetNonzeroRowIndices    ( manne );
        locJ = snGetNonzeroColumnPointers( manne );

        for ( i = 0 ; i < nT ; i++ ) {

                                      /* variables */

                                      v_bl[ i           ] = 3.05;
                                      v_bu[ i           ] = infty;

                                      v_bl[ i + nT      ] = xC_0;
                                      v_bu[ i + nT      ] = infty;

                                      v_bl[ i + nT + nT ] = xI_0;
                                      v_bu[ i + nT + nT ] = infty;

                                      v_xs[ i           ] = xK_0 + i / 10.0;
                                      v_xs[ i + nT      ] = xC_0;
                                      v_xs[ i + nT + nT ] = xI_0;

                                      v_hs[ i           ] = SN_BASIS_COLD_ELIGIBLE_0;
                                      v_hs[ i + nT      ] = SN_BASIS_COLD_ELIGIBLE_0;
                                      v_hs[ i + nT + nT ] = SN_BASIS_COLD_ELIGIBLE_0;

                                      /* constraints */

                                      c_bl[ i      ] = 0.0;
                                      c_bu[ i      ] = infty;

                                      c_bl[ i + nT ] = -infty;
                                      c_bu[ i + nT ] = 0.0;

                                      c_xs[ i      ] = SN_BASIS_COLD_ELIGIBLE_0;
                                      c_xs[ i + nT ] = SN_BASIS_COLD_ELIGIBLE_0;

                                      pi[ i      ] = -1.0;
                                      pi[ i + nT ] =  1.0;
                                      }

            /* exceptional bounds and initial values */

            v_bl[ 0 ] = 3.05;
            v_bu[ 0 ] = 3.05;

            v_bl[ 3 * nT - 3 ] = xI_0;
            v_bu[ 3 * nT - 3 ] = 0.112 * scale;

            v_bl[ 3 * nT - 2 ] = xI_0;
            v_bu[ 3 * nT - 2 ] = 0.114 * scale;

            v_bl[ 3 * nT - 1 ] = xI_0;
            v_bu[ 3 * nT - 1 ] = 0.116 * scale;

            v_xs[ nT - 1 ] = 1000.0;
            v_hs[ nT - 1 ] = SN_BASIS_IGNORE;

            c_bl[ nT - 1 ] =  0.0;
            c_bu[ nT - 1 ] = 10.0;

            c_bl[ nT + nT - 1 ] = -20.0;
            c_bu[ nT + nT - 1 ] =   0.0;

            /* nonzeros */

            p_locJ = locJ;

            /* column 1 */
            p_locJ[0] = 1;

            j = 0 ;   Jcol[j] = dummy;   indJ[j] = 1;
            j = 1 ;   Jcol[j] =  -1.0;   indJ[j] = 1 + nT;

            p_locJ[1] = p_locJ[0] + 2;   p_locJ++;

            /* columns 2 to nT-1 */
            for ( i = 1 ; i < nT - 1 ; i++ ) {

                                              j++;   Jcol[j] = dummy;   indJ[j] = i + 1;
                                              j++;   Jcol[j] =   1.0;   indJ[j] = i +     nT;
                                              j++;   Jcol[j] =  -1.0;   indJ[j] = i + 1 + nT;

                                              p_locJ[1] = p_locJ[0] + 3;  p_locJ++;

                                              }

                /* column nT */
                j++;   Jcol[j] =  dummy;   indJ[j] = nT;
                j++;   Jcol[j] =    1.0;   indJ[j] = 2 * nT - 1;
                j++;   Jcol[j] = growth;   indJ[j] = 2 * nT;

                p_locJ[1] = p_locJ[0] + 3;  p_locJ++;

                /* columns (nT+1) to (2*nT) */
                for ( i = 0 ; i < nT ; i++ ) {

                                              j++;   Jcol[j] = -1.0;   indJ[j] = i + 1;

                                              p_locJ[1] = p_locJ[0] + 1;  p_locJ++;

                                              }

                    /* columns (2*nT+1) to (3*nT) */
                    for ( i = 0 ; i < nT ; i++ ) {

                                                  j++;   Jcol[j] = -1.0;   indJ[j] = i + 1;
                                                  j++;   Jcol[j] = -1.0;   indJ[j] = i + 1 + nT;

                                                  p_locJ[1] = p_locJ[0] + 2;  p_locJ++;

                                                  }

                    end


                    void manne_delete()
                    {
                     free(b_t);
                     free(a_t);

                     }

                    =#


                    #=
                        function manneUserfun(
                                              cint   *mode,
                                              cint   *nnObj,
                                              cint   *nnCon,
                                              cint   *nnJac,
                                              cint   *nnL,
                                              cint   *negCon,
                                              double *x,
                                              double *fObj,
                                              double *gObj,
                                              double *fCon,
                                              double *gCon,
                                              cint   *state,
                                              char   *cu,
                                              cint   *lencu,
                                              double *iu,
                                              cint   *leniu,
                                              double *ru,
                                              cint   *lenru)

                            double  f   = 0.0;
                            double *C   = x + nT;
                            double *G_K = gObj;
                            double *G_C = gObj + nT;
                            double xK;
                            cint   i;

                            /* constraints */
                            for ( i = 0 ; i < nT ; i++ ) {
                                                          xK = x[i];
                                                          fCon[i] = a_t[i] * pow( xK , B );
                                                          gCon[i] = B * fCon[i] / xK;
                                                          }

                                /* objective */
                                for ( i = 0 ; i < nT ; i++ ) {
                                                              f      += b_t[i] * log( C[i] );
                                                              G_K[i]  = 0.0;
                                                              G_C[i]  = b_t[i] / C[i];
                                                              }

                                    *fObj = f;
                                end
                                =#




