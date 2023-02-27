/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_vec.h"

double
_d_vec_dot_heuristic(const double *vec1, const double *vec2, slong len2,
                     double *err)
{
    double psum = 0, nsum = 0, p, n, d, t;
    int pexp, nexp;
    slong i;

    for (i = 0; i < len2; i++)
    {
        t = vec1[i] * vec2[i];
        if (t >= 0)
            psum += t;
        else
            nsum += t;
    }
    nsum = -nsum;

    if (err != NULL)
    {
        p = frexp(psum, &pexp);
        n = frexp(nsum, &nexp);
        p = ldexp(1.0, pexp - D_BITS);
        n = ldexp(1.0, nexp - D_BITS);
        d = FLINT_MAX(p, n);
        *err = 2 * len2 * d;
    }

    return psum - nsum;
}
