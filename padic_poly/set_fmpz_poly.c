/*============================================================================

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

    Copyright (C) 2011 Sebastian Pancratz
 
******************************************************************************/

#include "fmpz_poly.h"
#include "padic_poly.h"

void padic_poly_set_fmpz_poly(padic_poly_t f, const fmpz_poly_t g, 
                              const padic_ctx_t ctx)
{
    const long len = g->length;

    padic_poly_fit_length(f, len);
    _padic_poly_set_length(f, len);
    _fmpz_vec_set(f->coeffs, g->coeffs, len);
    f->val = 0;

    padic_poly_canonicalise(f, ctx->p);
    padic_poly_reduce(f, ctx);
}

