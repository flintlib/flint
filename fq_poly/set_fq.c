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

    Copyright (C) 2012 Andres Goens

******************************************************************************/
#include "fq_poly.h"

void
fq_poly_set_fq(fmpz_poly_t poly, fq_ctx_t ctx, const fq_t c)
{
    if( fq_is_zero(c))
        fq_poly_zero(poly);
    else
    {
        fmpz_poly_fit_length(poly, 1);
        flint_calloc(1,sizeof(fq_struct));
        fmpz_set_ui(poly->coeffs, c);
        _fmpz_poly_set_length(poly, 1);
    }
}
