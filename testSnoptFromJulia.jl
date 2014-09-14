
# gah, have to manually install glfw: http://www.glfw.org/download.html
#Pkg.add("SymPy")
#using SymPy


include("Snopt.jl")
nT = 10;
scale = nT / 10.0;
numVariables = 3*nT;
numConstraints = 2*nT;
numNonzeros = 6*nT - 1;
numJacobianNonzeros = nT;
numNonlinearConstraints = nT;
numNonlinearJacobianVariables = nT;
numNonlinearObjectiveVariables = 2*nT;
linearObjectiveRow = 0;
s = Snopt.SnoptProblem(numVariables,
                 numConstraints,
                 numNonzeros,
                 numJacobianNonzeros,
                 numNonlinearConstraints,
                 numNonlinearObjectiveVariables,
                 numNonlinearJacobianVariables,
                 linearObjectiveRow,
                 "juliaSnotTest.out")
numVariables = Snopt.getNumVariables(s)
numConstraints = Snopt.getNumConstraints(s)
const inf = 1e20;
const growth = 0.03;
const Beta = 0.95;
const xK_0 = 3.00;
const xC_0 = 0.95;
const xI_0 = 0.05;
const B = 0.25;
scale = 0;
A = 0;
gFac = 0;
dummy = 0.1


#function setupManneExampleProblem(snoptProblem::Ptr{Void})
#{
#  cint     i, j;
#  cint    *p_locJ;
#  double   infty = 1e20, dummy = 0.1;

A = ( xC_0 + xI_0 ) / (xK_0^B);
gFac = (1.0 + growth)^(1-B);

b_t = zeros(nT); #(double*)malloc( sizeof(double) * nT )
a_t = zeros(nT);
a_t[1] = A * gFac;
b_t[1] = Beta;
for i = 2:nT
    a_t[i] = a_t[i-1] * gFac;
    b_t[i] = b_t[i-1] * Beta;
end
b_t[nT] /= ( 1.0 - Beta );


v_xs = Snopt.getVariablesArray(s, numVariables);
v_hs = Snopt.getVariableBasisEligibilitiesArray(s, numVariables);
v_bu = Snopt.getVariableUpperBoundsArray(s, numVariables);
v_bl = Snopt.getVariableLowerBoundsArray(s, numVariables);

c_xs = Snopt.getConstraintsArray(s, numConstraints);
c_hs = Snopt.getConstraintBasisEligibilitiesArray(s, numConstraints);
c_bu = Snopt.getConstraintUpperBoundsArray(s, numConstraints);
c_bl = Snopt.getConstraintLowerBoundsArray(s, numConstraints)

pi = Snopt.getMultipliersArray(s, numConstraints);

Jcol = Snopt.getNonzeroValuesArray(s, 3*nT);
indJ = Snopt.getNonzeroRowIndicesArray(s, 3*nT);
locJ = Snopt.getNonzeroColumnPointersArray(s, 3*nT);

infty = inf + 1e15


for i = 1:nT
    v_bl[ i           ] = 3.05;
    v_bu[ i           ] = infty;
    #    
    v_bl[ i + nT      ] = xC_0;
    v_bu[ i + nT      ] = infty;
    #
    v_bl[ i + nT + nT ] = xI_0;
    v_bu[ i + nT + nT ] = infty;
    #
    v_xs[ i           ] = xK_0 + i / 10.0;
    v_xs[ i + nT      ] = xC_0;
    v_xs[ i + nT + nT ] = xI_0;
    #
    v_hs[ i           ] = Snopt.SN_BASIS_COLD_ELIGIBLE_0;
    v_hs[ i + nT      ] = Snopt.SN_BASIS_COLD_ELIGIBLE_0;
    v_hs[ i + nT + nT ] = Snopt.SN_BASIS_COLD_ELIGIBLE_0;
    #
    # constraints
    #
    c_bl[ i      ] = 0.0;
    c_bu[ i      ] = infty;
    #
    c_bl[ i + nT ] = -infty;
    c_bu[ i + nT ] = 0.0;
    #
    c_xs[ i      ] = Snopt.SN_BASIS_COLD_ELIGIBLE_0;
    c_xs[ i + nT ] = Snopt.SN_BASIS_COLD_ELIGIBLE_0;
    #
    pi[ i      ] = -1.0;
    pi[ i + nT ] =  1.0;
end


#/* exceptional bounds and initial values */

v_bl[ 0 ] = 3.05;
v_bu[ 0 ] = 3.05;

v_bl[ 3 * nT - 3 ] = xI_0;
v_bu[ 3 * nT - 3 ] = 0.112 * scale;

v_bl[ 3 * nT - 2 ] = xI_0;
v_bu[ 3 * nT - 2 ] = 0.114 * scale;

v_bl[ 3 * nT - 1 ] = xI_0;
v_bu[ 3 * nT - 1 ] = 0.116 * scale;

v_xs[ nT - 1 ] = 1000.0;
v_hs[ nT - 1 ] = Snopt.SN_BASIS_IGNORE;

c_bl[ nT - 1 ] =  0.0;
c_bu[ nT - 1 ] = 10.0;

c_bl[ nT + nT - 1 ] = -20.0;
c_bu[ nT + nT - 1 ] =   0.0;

# /* nonzeros */


#/* column 1 */
p_locJ = locJ;
p_locJ[0] = 1;
j = 0 ;   Jcol[j] = dummy;   indJ[j] = 1;
j = 1 ;   Jcol[j] =  -1.0;   indJ[j] = 1 + nT;

p_locJ[1] = p_locJ[0] + 2;   p_locJ++;

#/* columns 2 to nT-1 */
for i = 1:(nT - 2)
    j++;   Jcol[j] = dummy;   indJ[j] = i + 1;
    j++;   Jcol[j] =   1.0;   indJ[j] = i +     nT;
    j++;   Jcol[j] =  -1.0;   indJ[j] = i + 1 + nT;
    p_locJ[1] = p_locJ[0] + 3;  p_locJ++;
end

# /* column nT */
j++;   Jcol[j] =  dummy;   indJ[j] = nT;
j++;   Jcol[j] =    1.0;   indJ[j] = 2 * nT - 1;
j++;   Jcol[j] = growth;   indJ[j] = 2 * nT;

p_locJ[1] = p_locJ[0] + 3;  p_locJ++;

#/* columns (nT+1) to (2*nT) */
for i = 0:(nT-1)               
    j++;   Jcol[j] = -1.0;   indJ[j] = i + 1;
    p_locJ[1] = p_locJ[0] + 1;  p_locJ++;
end

#  /* columns (2*nT+1) to (3*nT) */
for i = 0:(nT-1)
    j++;   Jcol[j] = -1.0;   indJ[j] = i + 1;
    j++;   Jcol[j] = -1.0;   indJ[j] = i + 1 + nT;
    p_locJ[1] = p_locJ[0] + 2;  p_locJ++;
end






show(s)











#ENV["LD_LIBRARY_PATH"] = ENV["LD_LIBRARY_PATH"] * ":" * homedir() * "/workspace/snopt_cpp/lib"
show(ENV["LD_LIBRARY_PATH"])

#include "snopt.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>




#= this is an example of how you can call c from julia... this works
function mycompare{T}(a_::Ptr{T}, b_::Ptr{T})
    a = unsafe_load(a_)
    b = unsafe_load(b_)
    return a < b ? cint(-1) : a > b ? cint(+1) : cint(0)
end
cint(n) = convert(Cint, n)
const mycompare_C = cfunction(mycompare, Cint, (Ptr{Cdouble}, Ptr{Cdouble}))


A = [1.3, -2.7, 4.4, 3.1]
ccall(:qsort, Void, (Ptr{Cdouble}, Csize_t, Csize_t, Ptr{Void}),
      A, length(A), sizeof(eltype(A)), mycompare_C)
show(A)
=#




#function main()

#  cint    i;
#  char   *error_msg;
#  cint    status;
#  cint    variables;
#  cint    constraints;
#  cint    linear_objective_row;
#  cint    nonzeros;
#  cint    jacobian_nonzeros;
#  cint    nonlinear_constraints;
# cint    nonlinear_jacobian_variables;
# cint    nonlinear_objective_variables;

#const gsl_brent = unsafe_load(cglobal((:gsl_min_fminimizer_brent,:libgsl), Ptr{Void}))


#manneReference = c_malloc(10000)
#errorMessage = Ptr{Uint8}


if (false)

    #outputFilename = "manne.out"
    #typeof(outputFilename)
    #argv = [ "a.out", "arg1", "arg2" ]
    
    
    
    #libcsnopt.so 
    
    status = ccall((:snInit, "libcsnopt"), Cint, (Ptr{Void},
                                                  Cint,
                                                  Cint,
                                                  Cint,
                                                  Cint,
                                                  Cint,
                                                  Cint,
                                                  Cint,
                                                  Cint,
                                                  Ptr{Uint8},
                                                  Ptr{Uint8}),
                   manneReference,
                   variables,
                   constraints,
                   nonzeros,
                   jacobian_nonzeros,
                   nonlinear_constraints,
                   nonlinear_objective_variables,
                   nonlinear_jacobian_variables,
                   linear_objective_row,
                   "manne.out",
                   "stdout")


    #const snoptOK = unsafe_load(cglobal((:SN_OK, "libcsnopt"), Ptr{Void}))
    snoptOK = 1;

    
    if ( status != snoptOK )
        #ccall((:snGetError, "libcsnopt"), Void, (Ptr{Void}, Ptr{Uint8}))
        
        #snGetError( &manne , &error_msg );
        #printf( "%s: error occurred (%d).\n%s", argv[0] , status, error_msg );
        
        #snDelete( &manne );
        #free(error_msg);
        
        #exit(EXIT_FAILURE);
    end
    
    
    #jrebula@titian$ gcc -shared -o libcsnopt.so snopt.o snerror.o -lsnopt7 -lsnblas -lsnprint7
    
    #put all the snopt libs somewhere ld can find them: /usr/lib or whatever
    # gcc -fpic -Wall -c snset.c
    # gcc -fpic -Wall -c snopt.c
    # gcc -fpic -Wall -c snerror.c 
    # gfortran -fpic -Wall -c sn_open.f
    #gcc -shared -o libcsnopt.so snopt.o snerror.o snset.o sn_open.o -lsnopt7 -lsnblas -lsnprint7 -lgfortran
    
    
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

  snSetProblemName        ( &manne, "  s           "       );
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












#void manne_init( snProblem*s)
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

  v_xs = snGetVariables  (s);
  v_hs = snGetVariableBasisEligibilities  (s);

  c_xs = snGetConstraints(s);
  c_hs = snGetConstraintBasisEligibilities(s);

  v_bu = snGetVariableUpperBounds  (s);
  v_bl = snGetVariableLowerBounds  (s);

  c_bu = snGetConstraintUpperBounds(s);
  c_bl = snGetConstraintLowerBounds(s);
  pi   = snGetMultipliers          (s);

  Jcol = snGetNonzeroValues        (s);
  indJ = snGetNonzeroRowIndices    (s);
  locJ = snGetNonzeroColumnPointers(s);

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



