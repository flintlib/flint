/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_lll.h"

double
fmpz_lll_heuristic_dot(const double *vec1, const double *vec2, slong len2,
                       const fmpz_mat_t B, slong k, slong j, slong exp_adj)
{
    double err;
    double sum = _d_vec_dot_heuristic(vec1, vec2, len2, &err);

    if (err > ldexp(1, -D_BITS / 2))
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
