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

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "d_mat.h"

int
d_mat_is_reduced(const d_mat_t R, double delta, double eta)
{
    double tmp1, tmp2;
    int reduced = 1;
    slong i, d = R->r;

    if (d == 1)
        return 1;

    for (i = 0; (i < d - 1) && (reduced == 1); i++)
    {
        slong j;

        tmp1 = d_mat_entry(R, i + 1, i) * d_mat_entry(R, i + 1, i);
        tmp2 = d_mat_entry(R, i + 1, i + 1) * d_mat_entry(R, i + 1, i + 1);
        tmp1 += tmp2;

        tmp2 = d_mat_entry(R, i, i) * d_mat_entry(R, i, i);
        tmp2 *= delta;

        tmp1 -= tmp2;

        if (tmp1 < 0)
        {
            reduced = 0;
            /* flint_printf("happened at index i = %ld\n", i); */
            break;
        }

        for (j = 0; (j < i) && (reduced == 1); j++)
        {
            tmp2 = d_mat_entry(R, i + 1, i + 1) * eta;
            if (fabs(d_mat_entry(R, j, i)) > fabs(tmp2))
            {
                reduced = 0;
                /* flint_printf
                   ("size reduction problem at index i = %ld, j = %ld\n", i,
                   j); */
                break;
            }
        }
    }

    return reduced;
}
