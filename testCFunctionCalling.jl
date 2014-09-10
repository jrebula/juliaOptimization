

function mycompare{T}(aIn::Ptr{T}, bIn::Ptr{T})
  a = unsafe_load(aIn)
  b = unsafe_load(bIn)
  println("mycompare($a,$b)")
  return a < b ? cint(-1) : a > b ? cint(+1) : cint(0)
end

cint(n) = convert(Cint, n)

const mycompare_c = cfunction(mycompare, Cint, (Ptr{Cdouble}, Ptr{Cdouble}))


A = [1.3, -2.7, 4.4, 3.1]
#ccall(:qsort, Void, (Ptr{Cdouble}, Csize_t, Csize_t, Ptr{Void}),
#      A, length(A), sizeof(eltype(A)), mycompare_c)

show(A)




#= won't work:
function qsort!{T}(A::Vector{T}, lessthan::Function)
    function mycompare(a_::Ptr{T}, b_::Ptr{T})
        a = unsafe_load(a_)
        b = unsafe_load(b_)
        return lessthan(a, b) ? cint(-1) : cint(+1)
    end
    mycompare_c = cfunction(mycompare, Cint, (Ptr{T}, Ptr{T}))
    ccall(:qsort, Void, (Ptr{T}, Csize_t, Csize_t, Ptr{Void}),
          A, length(A), sizeof(T), mycompare_c)
    A
end
qsort!([1.3, -2.7, 4.4, 3.1], <)

=#

function qsort!_compare{T}(a_::Ptr{T}, b_::Ptr{T}, lessthan_::Ptr{Void})
    a = unsafe_load(a_)
    b = unsafe_load(b_)
    lessthan = unsafe_pointer_to_objref(lessthan_)::Function
    return lessthan(a, b) ? cint(-1) : cint(+1)
end

function myQsort!{T}(A::Vector{T}, lessthan::Function=<)
    compare_c = cfunction(qsort!_compare, Cint, (Ptr{T}, Ptr{T}, Ptr{Void}))
    ccall(:qsort_r, Void, (Ptr{T}, Csize_t, Csize_t, Ptr{Void}, Any),
          A, length(A), sizeof(T), compare_c, lessthan)
    return A
end


#myQsort!(A, (x, y) -> x < y)
myQsort!(A)

show(A)



##
function gsl_function_wrap(x::Cdouble, params::Ptr{Void})
    f = unsafe_pointer_to_objref(params)::Function
    convert(Cdouble, f(x))::Cdouble
end
const gsl_function_wrap_c = cfunction(gsl_function_wrap,
                                      Cdouble, (Cdouble, Ptr{Void}))

type GSL_Function
    func::Ptr{Void}
    params::Any
    GSL_Function(f::Function) = new(gsl_function_wrap_c, f)
end

function gsl_integration_qag(f::Function, a::Real, b::Real, epsrel::Real=1e-12,
                             maxintervals::Integer=10^7)
    s = ccall((:gsl_integration_workspace_alloc,:libgsl), Ptr{Void}, (Csize_t,),
              maxintervals)
    result = Array(Cdouble,1)
    abserr = Array(Cdouble,1)
    ccall((:gsl_integration_qag,:libgsl), Cint,
          (Ptr{GSL_Function}, Cdouble,Cdouble, Cdouble, Csize_t, Cint, Ptr{Void},
           Ptr{Cdouble}, Ptr{Cdouble}),
          &GSL_Function(f), a, b, epsrel, maxintervals, 1, s, result, abserr)
    ccall((:gsl_integration_workspace_free,:libgsl), Void, (Ptr{Void},), s)
    return (result[1], abserr[1])
end



truth = (x) -> x^2

for t = linspace(0, 10, 10)
  (integrated, absErr) = gsl_integration_qag((x) -> 2x, 0, t)
  errVal = integrated - truth(t)
  println("$integrated, $(truth(t)), verVal = $errVal")
end






run(`echo \$PATH`)


type GSL_Minimizer
    m::Ptr{Void} # the gsl_min_fminimizer pointer
    f::Any  # explicit reference to objective, to prevent garbage-collection
    function GSL_Minimizer(t)
       m = ccall((:gsl_min_fminimizer_alloc,:libgsl), Ptr{Void}, (Ptr{Void},), t)
       p = new(m, nothing)
       finalizer(p, p -> ccall((:gsl_min_fminimizer_free,:libgsl),
                               Void, (Ptr{Void},), p.m))
       p
    end
end

function gsl_minimizer_set!(m::GSL_Minimizer, f, x0, xmin, xmax)
    ccall((:gsl_min_fminimizer_set,:libgsl), Cint,
          (Ptr{Void}, Ptr{GSL_Function}, Cdouble, Cdouble, Cdouble),
          m.m, &GSL_Function(f), x0, xmin, xmax)
    m.f = f
    m
end

const gsl_brent = unsafe_load(cglobal((:gsl_min_fminimizer_brent,:libgsl), Ptr{Void}))
GSL_Minimizer() = GSL_Minimizer(gsl_brent)

gsl_minimizer_iterate!(m::GSL_Minimizer) =
    ccall((:gsl_min_fminimizer_iterate,:libgsl), Cint, (Ptr{Void},), m.m)

gsl_minimizer_x(m::GSL_Minimizer) =
    ccall((:gsl_min_fminimizer_x_minimum,:libgsl), Cdouble, (Ptr{Void},), m.m)

gsl_minimizer_f(m::GSL_Minimizer) =
    ccall((:gsl_min_fminimizer_f_minimum,:libgsl), Cdouble, (Ptr{Void},), m.m)

gsl_minimizer_xmin(m::GSL_Minimizer) =
    ccall((:gsl_min_fminimizer_x_lower,:libgsl), Cdouble, (Ptr{Void},), m.m)
gsl_minimizer_xmax(m::GSL_Minimizer) =
    ccall((:gsl_min_fminimizer_x_upper,:libgsl), Cdouble, (Ptr{Void},), m.m)



m = GSL_Minimizer()
gsl_minimizer_set!(m, sin, -1, -3, 1)
while gsl_minimizer_xmax(m) - gsl_minimizer_xmin(m) > 1e-6
    println("iterating at x = $(gsl_minimizer_x(m))")
    gsl_minimizer_iterate!(m)
end

println("found minimum $(gsl_minimizer_f(m)) at x = $(gsl_minimizer_x(m))")


