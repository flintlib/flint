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

#include "fq_nmod.h"

void fq_nmod_mul_fmpz(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx)
{
    fmpz_t rx;
    fmpz_init(rx);

    fmpz_mod(rx, x, fq_nmod_ctx_prime(ctx));
    nmod_poly_scalar_mul_nmod(rop, op, fmpz_get_ui(rx));

    fmpz_clear(rx);
}
