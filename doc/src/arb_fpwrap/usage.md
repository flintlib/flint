---
title: Using the floating-point wrappers
summary: How to call the arb_fpwrap functions from C, Python and Julia, how to handle their status codes, and what accuracy they guarantee.
order: 1
---

The functions of [](arb_fpwrap) are meant for programs that work with
`double` numbers and want special functions with a guaranteed accuracy,
without using ball arithmetic themselves. Every function has the form

```c
int arb_fpwrap_double_NAME(double * res, double x, ..., int flags);
int arb_fpwrap_cdouble_NAME(complex_double * res, complex_double x, ..., int flags);
```

It writes the result to `res` and returns `FPWRAP_SUCCESS` (0) or, if
it cannot reach the required accuracy (for example at a pole, outside
the real domain of a real function, or when the work limit is reached),
`FPWRAP_UNABLE` (1) with a NaN result. Internally the function is
evaluated in ball arithmetic at 64 bits, then at 128 bits, and so on,
until the enclosure determines the `double` result; a successful result
is therefore never wrong by more than one unit in the last place, and
with `FPWRAP_CORRECT_ROUNDING` it is correctly rounded.

## Calling from C {#c}

```c example name=arb_fpwrap_basic
#include <stdio.h>
#include <flint/arb_fpwrap.h>

int main(void)
{
    double x, y;
    complex_double cx, cy;
    int status;

    x = 2.0;
    status = arb_fpwrap_double_zeta(&y, x, 0);
    printf("zeta(%g) = %.16g  (status %d)\n", x, y, status);

    cx.real = 0.5;
    cx.imag = 123.0;
    status = arb_fpwrap_cdouble_zeta(&cy, cx, 0);
    printf("zeta(%g + %gi) = %.16g + %.16gi  (status %d)\n",
        cx.real, cx.imag, cy.real, cy.imag, status);

    /* a pole: no finite result */
    status = arb_fpwrap_double_zeta(&y, 1.0, 0);
    printf("zeta(1): status %d, result is NaN: %d\n", status, y != y);

    /* cancellation: sin(2^60) needs extra precision internally */
    status = arb_fpwrap_double_sin(&y, 1152921504606846976.0, FPWRAP_CORRECT_ROUNDING);
    printf("sin(2^60) = %.17g  (status %d)\n", y, status);

    /* the 1000th zero of the zeta function */
    status = arb_fpwrap_cdouble_zeta_zero(&cy, 1000, 0);
    printf("rho_1000 = %.16g + %.16gi  (status %d)\n", cy.real, cy.imag, status);

    flint_cleanup();
    return 0;
}
```
```output
zeta(2) = 1.644934066848226  (status 0)
zeta(0.5 + 123i) = 0.006252861175594465 + 0.08206030514520983i  (status 0)
zeta(1): status 1, result is NaN: 1
sin(2^60) = -0.83064921763725463  (status 0)
rho_1000 = 0.5 + 1419.422480945996i  (status 0)
```

The wrappers cache constants and other data as the ball functions do;
`flint_cleanup` frees these caches at the end of the program (this is
only needed to make memory checkers quiet). The functions are
thread-safe.

The `flags` argument combines `FPWRAP_ACCURATE_PARTS` (accurate real
and imaginary parts separately), `FPWRAP_CORRECT_ROUNDING` and a
multiple of `FPWRAP_WORK_LIMIT` (the maximal working precision). For
example, `FPWRAP_CORRECT_ROUNDING + 10 * FPWRAP_WORK_LIMIT` asks for a
correctly rounded result and allows working precisions up to $64 \cdot
2^{10}$ bits.

## Calling from Python {#python}

The functions can be called from Python with `ctypes`; since FLINT
3.0 they are part of `libflint` (in Arb they were in `libarb`). The
following snippet is not tested automatically:

```python
import ctypes
import ctypes.util

libflint = ctypes.CDLL(ctypes.util.find_library('flint'))

class _complex_double(ctypes.Structure):
    _fields_ = [('real', ctypes.c_double),
                ('imag', ctypes.c_double)]

def wrap_double_fun(fun):
    def f(x):
        y = ctypes.c_double()
        if fun(ctypes.byref(y), ctypes.c_double(x), 0):
            raise ValueError(f"unable to evaluate function accurately at {x}")
        return y.value
    return f

def wrap_cdouble_fun(fun):
    fun.argtypes = [ctypes.POINTER(_complex_double), _complex_double, ctypes.c_int]
    def f(x):
        x = complex(x)
        cx = _complex_double(x.real, x.imag)
        cy = _complex_double()
        if fun(ctypes.byref(cy), cx, 0):
            raise ValueError(f"unable to evaluate function accurately at {x}")
        return complex(cy.real, cy.imag)
    return f

zeta = wrap_double_fun(libflint.arb_fpwrap_double_zeta)
czeta = wrap_cdouble_fun(libflint.arb_fpwrap_cdouble_zeta)
print(zeta(2.0))          # 1.6449340668482264
print(czeta(0.5+1e9j))    # (-2.761748029838061-1.6775122409894598j)
print(zeta(1.0))          # pole: raises ValueError
```

The package `python-flint` (see [](guide/interfaces/python)) offers
the ball arithmetic directly.

## Calling from Julia {#julia}

With `ccall`, the wrappers can be called from Julia; `Complex{Float64}`
has the layout of `complex_double` (not tested automatically):

```julia
function zeta(x::Float64)
    cy = Ref{Float64}()
    if Bool(ccall((:arb_fpwrap_double_zeta, :libflint), Cint,
                  (Ptr{Float64}, Float64, Cint), cy, x, 0))
        error("unable to evaluate accurately at ", x)
    end
    return cy[]
end

function zeta(x::Complex{Float64})
    cy = Ref{Complex{Float64}}()
    if Bool(ccall((:arb_fpwrap_cdouble_zeta, :libflint), Cint,
                  (Ptr{Complex{Float64}}, Complex{Float64}, Cint), cy, x, 0))
        error("unable to evaluate accurately at ", x)
    end
    return cy[]
end

println(zeta(2.0))             # 1.6449340668482264
println(zeta(0.5 + 1e9im))     # -2.761748029838061 - 1.6775122409894598im
println(zeta(1.0))             # pole: throws an error
```

## Cost {#cost}

Most calls succeed at the first precision of 64 bits, and the cost is
then that of one evaluation of the ball function at 64 bits, typically
between a few hundred nanoseconds for elementary functions and tens of
microseconds for hypergeometric functions, which is slower than a
floating-point library but without its accuracy problems. Inputs where
the function is ill-conditioned (for example near its zeros) and
correct rounding of values close to a rounding boundary need higher
precisions; each failed attempt doubles the precision, so the total
cost is dominated by the last attempt. Large arguments such as $2^{60}$
in $\sin(2^{60})$ are handled by the ball functions themselves, which
add the precision needed for the argument reduction.
