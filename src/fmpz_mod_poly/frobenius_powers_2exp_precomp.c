/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_frobenius_powers_2exp_precomp(
           fmpz_mod_poly_frobenius_powers_2exp_t pow, const fmpz_mod_poly_t f,
                 const fmpz_mod_poly_t finv, ulong m, const fmpz_mod_ctx_t ctx)
{
    slong i, l = 0;

    if (m == 0)
    {
       pow->len = 0;

       return;
    }

    l = FLINT_CLOG2(m);
    if ((WORD(1) << l) == m)
       l++;

    pow->pow = (fmpz_mod_poly_struct *) flint_malloc(l*sizeof(fmpz_mod_poly_struct));

    for (i = 0; i < l; i++)
       fmpz_mod_poly_init(pow->pow + i, ctx);

    pow->len = l;

   fmpz_mod_poly_powmod_x_fmpz_preinv(pow->pow + 0,
                                      fmpz_mod_ctx_modulus(ctx), f, finv, ctx);

    for (i = 1; i < l; i++)
       fmpz_mod_poly_compose_mod(pow->pow + i, pow->pow + i - 1,
                                                     pow->pow + i - 1, f, ctx);
}
