/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#include <float.h>
#undef ulong
#include "flint.h"
#include "ulong_extras.h"
#include "longlong.h"

double
n_cbrt_estimate(double a)
{ 
    typedef union { 
        slong      uword_val;
#if FLINT64
        double     double_val;
#else
        float      double_val;
#endif
    } uni;

    uni alias;
    ulong n, hi, lo;

#ifdef FLINT64
    const mp_limb_t mul_factor = UWORD(6148914691236517205);
    slong s = UWORD(4607182418800017408);      /* ((1 << 10) - 1) << 52 */
#else
    const mp_limb_t mul_factor = UWORD(1431655765);
    slong s = UWORD(1065353216);               /* ((1 << 7) - 1 << 23)  */
#endif

    alias.double_val = a;
    n = alias.uword_val;
    n -= s;
    umul_ppmm(hi, lo, n, mul_factor);
    n = hi;
    n += s;
    alias.uword_val = n;
    return alias.double_val;
}
