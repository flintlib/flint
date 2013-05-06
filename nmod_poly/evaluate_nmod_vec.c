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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

void
_nmod_poly_evaluate_nmod_vec(mp_ptr ys, mp_srcptr coeffs, long len,
    mp_srcptr xs, long n, nmod_t mod)
{
    if (len < 32)
        _nmod_poly_evaluate_nmod_vec_iter(ys, coeffs, len, xs, n, mod);
    else
        _nmod_poly_evaluate_nmod_vec_fast(ys, coeffs, len, xs, n, mod);
}

void
nmod_poly_evaluate_nmod_vec(mp_ptr ys,
        const nmod_poly_t poly, mp_srcptr xs, long n)
{
    _nmod_poly_evaluate_nmod_vec(ys, poly->coeffs,
                                        poly->length, xs, n, poly->mod);
}
