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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_sub) (TEMPLATE(T, mat_t) res,
                      const TEMPLATE(T, mat_t) mat1,
                      const TEMPLATE(T, mat_t) mat2,
                      const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    if (res->c < 1)
        return;

    for (i = 0; i < res->r; i++)
        _TEMPLATE(T, vec_sub) (res->rows[i], mat1->rows[i], mat2->rows[i],
                               res->c, ctx);
}


#endif
