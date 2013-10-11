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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_randtril)(TEMPLATE(T, mat_t) mat, flint_rand_t state, int unit,
                          const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
        {
            if (j < i)
            {
                TEMPLATE(T, randtest)(TEMPLATE(T, mat_entry)(mat, i, j), state, ctx);
            }
            else if (i == j)
            {
                TEMPLATE(T, randtest)(TEMPLATE(T, mat_entry)(mat, i, j), state, ctx);
                if (unit || TEMPLATE(T, is_zero)(TEMPLATE(T, mat_entry)(mat, i, j), ctx))
                    TEMPLATE(T, one)(TEMPLATE(T, mat_entry)(mat, i, j), ctx);
            }
            else
            {
                TEMPLATE(T, zero)(TEMPLATE(T, mat_entry)(mat, i, j), ctx);
            }
        }
    }
}


#endif
