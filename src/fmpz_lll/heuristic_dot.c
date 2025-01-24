/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_vec.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"

#ifdef __GNUC__
# define ldexp __builtin_ldexp
#else
# include <math.h>
#endif

double
fmpz_lll_heuristic_dot(const double *vec1, const double *vec2, slong len2,
                       const fmpz_mat_t B, slong k, slong j, slong exp_adj)
{
   double sum = _d_vec_dot(vec1, vec2, len2);
   double tmp = _d_vec_norm(vec1, len2);
   double tmp2 = _d_vec_norm(vec2, len2);

   tmp = tmp*tmp2 * ldexp(1.0, -70);
   tmp2 = sum*sum;

   if (tmp2 <= tmp)
    {
        slong exp;
        fmpz_t sp;
        fmpz_init(sp);
        _fmpz_vec_dot(sp, fmpz_mat_row(B, k), fmpz_mat_row(B, j), len2);
        sum = fmpz_get_d_2exp(&exp, sp);
        sum = ldexp(sum, sum - exp_adj);
        fmpz_clear(sp);
    }

    return sum;
}
