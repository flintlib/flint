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

    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

slong
TEMPLATE(T, mat_rank) (const TEMPLATE(T, mat_t) A,
                       const TEMPLATE(T, ctx_t) ctx)
{
    slong m, n, rank;
    slong *perm;
    TEMPLATE(T, mat_t) tmp;

    m = A->r;
    n = A->c;

    if (m == 0 || n == 0)
        return 0;

    TEMPLATE(T, mat_init_set) (tmp, A, ctx);
    perm = flint_malloc(sizeof(slong) * m);

    rank = TEMPLATE(T, mat_lu) (perm, tmp, 0, ctx);

    flint_free(perm);
    TEMPLATE(T, mat_clear) (tmp, ctx);
    return rank;
}


#endif
