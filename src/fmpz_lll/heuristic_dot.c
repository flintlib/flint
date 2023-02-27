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

#include "fmpz_lll.h"

double
fmpz_lll_heuristic_dot(const double *vec1, const double *vec2, slong len2,
                       const fmpz_mat_t B, slong k, slong j, slong exp_adj)
{
   double sum = _d_vec_dot(vec1, vec2, len2);
   double tmp = _d_vec_norm(vec1, len2);
   double tmp2 = _d_vec_norm(vec2, len2);

   tmp = ldexp(tmp*tmp2, -70);
   tmp2 = sum*sum;

   if (tmp2 <= tmp)
    {
        slong exp;
        fmpz_t sp;
        fmpz_init(sp);
        _fmpz_vec_dot(sp, B->rows[k], B->rows[j], len2);
        sum = fmpz_get_d_2exp(&exp, sp);
        sum = ldexp(sum, sum - exp_adj);
        fmpz_clear(sp);
    }

    return sum;
}
