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

#include "d_vec.h"

double
_d_vec_dot_heuristic(const double *vec1, const double *vec2, slong len2,
                     double *err)
{
    double psum = 0, nsum = 0, p, n, d;
    ulong *ps, *ns;
    slong i;

    for (i = 0; i < len2; i++)
    {
        if (vec1[i] * vec2[i] >= 0)
            psum += vec1[i] * vec2[i];
        else
            nsum += -vec1[i] * vec2[i];
    }

    p = psum;
    n = nsum;
    ps = (ulong *) & psum;
    ns = (ulong *) & nsum;
    *ps += 1;
    *ns += 1;
    d = FLINT_MAX(psum - p, nsum - n);
    *err = 2 * len2 * d;

    return psum - nsum;
}
