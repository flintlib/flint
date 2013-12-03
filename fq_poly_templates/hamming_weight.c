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

    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

slong
_TEMPLATE(T, poly_hamming_weight) (const TEMPLATE(T, struct) * op, slong len,
                                   const TEMPLATE(T, ctx_t) ctx)
{
    slong i, sum = 0;
    for (i = 0; i < len; i++)
        sum += !TEMPLATE(T, is_zero) (op + i, ctx);

    return sum;
}

slong
TEMPLATE(T, poly_hamming_weight) (const TEMPLATE(T, poly_t) op,
                                  const TEMPLATE(T, ctx_t) ctx)
{

    return _TEMPLATE(T, poly_hamming_weight) (op->coeffs, op->length, ctx);
}


#endif
