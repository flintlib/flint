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

#include "fq_nmod.h"

void
fq_nmod_ctx_randtest_reducible(fq_nmod_ctx_t ctx, flint_rand_t state)
{
    fmpz_t p;
    slong d;
    nmod_poly_t mod;

    fmpz_init(p);
    fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 6), 1));
    d = n_randint(state, 10) + 1;
    nmod_poly_init(mod, fmpz_get_ui(p));
    nmod_poly_randtest_monic(mod, state, d + 1);
    fq_nmod_ctx_init_modulus(ctx, mod, "a");
    nmod_poly_clear(mod);
    fmpz_clear(p);
}
