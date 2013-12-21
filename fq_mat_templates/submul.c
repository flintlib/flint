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
TEMPLATE(T, mat_submul) (TEMPLATE(T, mat_t) D,
                         const TEMPLATE(T, mat_t) C,
                         const TEMPLATE(T, mat_t) A,
                         const TEMPLATE(T, mat_t) B,
                         const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, mat_t) tmp;
    TEMPLATE(T, mat_init) (tmp, A->r, B->c, ctx);
    TEMPLATE(T, mat_mul) (tmp, A, B, ctx);
    TEMPLATE(T, mat_sub) (D, C, tmp, ctx);
    TEMPLATE(T, mat_clear) (tmp, ctx);
}


#endif
