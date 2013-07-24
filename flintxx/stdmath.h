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

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/

#ifndef FLINT_CXX_STDMATH_H
#define FLINT_CXX_STDMATH_H

#include "expression.h"

namespace flint {
FLINT_DEFINE_BINOP(pow)
FLINT_DEFINE_BINOP(root)

FLINT_DEFINE_UNOP(abs)
FLINT_DEFINE_UNOP(exp)
FLINT_DEFINE_UNOP(log)
FLINT_DEFINE_UNOP(sqrt)
FLINT_DEFINE_UNOP(inv)

FLINT_DEFINE_UNOP(sqr)

FLINT_DEFINE_UNOP(height)

FLINT_DEFINE_BINOP(mul_2exp)
FLINT_DEFINE_BINOP(fdiv_2exp)
FLINT_DEFINE_BINOP(tdiv)
FLINT_DEFINE_BINOP(tdiv_2exp)
FLINT_DEFINE_BINOP(divexact)

FLINT_DEFINE_BINOP(smod) // "symmetric" %

// poly functions
FLINT_DEFINE_BINOP(reverse)
FLINT_DEFINE_BINOP(resultant)
FLINT_DEFINE_BINOP(poly_bit_pack)
FLINT_DEFINE_BINOP(poly_bit_unpack)
FLINT_DEFINE_BINOP(poly_bit_unpack_unsigned)
FLINT_DEFINE_BINOP(poly_shift_left)
FLINT_DEFINE_BINOP(poly_shift_right)
FLINT_DEFINE_UNOP(poly_2norm)
FLINT_DEFINE_UNOP(content)
FLINT_DEFINE_UNOP(primitive_part)
FLINT_DEFINE_BINOP(inv_series)
FLINT_DEFINE_UNOP(derivative)
FLINT_DEFINE_BINOP(evaluate)
FLINT_DEFINE_BINOP(compose)

// either compose or evaluate. Currently implemented in fmpz_polyxx.h.
FLINT_DEFINE_BINOP(compeval)
} // flint

namespace std {
FLINT_DEFINE_BINOP_HERE(pow)

FLINT_DEFINE_UNOP_HERE(abs)
FLINT_DEFINE_UNOP_HERE(exp)
FLINT_DEFINE_UNOP_HERE(log)
FLINT_DEFINE_UNOP_HERE(sqrt)
}

#endif
