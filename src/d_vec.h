/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef D_VEC_H
#define D_VEC_H

#ifdef D_VEC_INLINES_C
#define D_VEC_INLINE
#else
#define D_VEC_INLINE static inline
#endif

#include "flint.h"
#include "double_extras.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  Memory management  *******************************************************/

double * _d_vec_init(slong len);

void _d_vec_clear(double * vec);

/*  Randomisation  ***********************************************************/

void _d_vec_randtest(double * f, flint_rand_t state,
                        slong len, slong minexp, slong maxexp);

/*  Assignment and basic manipulation  ***************************************/

void _d_vec_set(double * vec1, const double * vec2, slong len2);

void _d_vec_zero(double * vec, slong len);

/*  Comparison  **************************************************************/

int _d_vec_equal(const double * vec1, const double * vec2, slong len);

int _d_vec_approx_equal(const double * vec1, const double * vec2, slong len, double eps);

int _d_vec_is_zero(const double * vec, slong len);

int _d_vec_is_approx_zero(const double * vec, slong len, double eps);

/*  Arithmetic  ****************************************************************/

void _d_vec_add(double * res, const double * vec1, const double * vec2, slong len2);

void _d_vec_sub(double * res, const double * vec1, const double * vec2, slong len2);

D_VEC_INLINE
void _d_vec_mul_2exp(double * res, const double * x, slong len, int e)
{
    slong i;

    if (e >= D_MIN_NORMAL_EXPONENT && e <= D_MAX_NORMAL_EXPONENT)
    {
        double_uint64_u u;
        u.i = ((int64_t) (e + D_EXPONENT_BIAS)) << D_EXPONENT_SHIFT;

        for (i = 0; i < len; i++)
            res[i] = x[i] * u.f;
    }
    else
    {
        for (i = 0; i < len; i++)
            res[i] = ldexp(x[i], e);
    }
}

/*  Dot product and norm  **************************************/

double _d_vec_dot(const double * vec1, const double * vec2, slong len2);

double _d_vec_norm(const double * vec, slong len);

double _d_vec_dot_heuristic(const double * vec1, const double * vec2, slong len2, double * err);

double _d_vec_dot_thrice(const double * vec1, const double * vec2, slong len2, double * err);

#ifdef __cplusplus
}
#endif

#endif
