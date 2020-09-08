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

void
fmpz_mod_poly_frobenius_powers_precomp(fmpz_mod_poly_frobenius_powers_t pow, 
                 const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, ulong m,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i;

    pow->pow = (fmpz_mod_poly_struct *) flint_malloc((m + 1)*sizeof(fmpz_mod_poly_struct));

    for (i = 0; i <= m; i++)
       fmpz_mod_poly_init(pow->pow + i, ctx);

    pow->len = m;

    /* x mod f */
    fmpz_mod_poly_set_coeff_ui(pow->pow + 0, 1, 1, ctx);
    fmpz_mod_poly_set_coeff_ui(pow->pow + 0, 0, 0, ctx);
    _fmpz_mod_poly_set_length(pow->pow + 0, 2);
    if (f->length <= 2)
       fmpz_mod_poly_rem(pow->pow + 0, pow->pow + 0, f, ctx);

    if (m >= 1)
       fmpz_mod_poly_powmod_x_fmpz_preinv(pow->pow + 1,
                                      fmpz_mod_ctx_modulus(ctx), f, finv, ctx);

    for (i = 2; i <= m; i++)
       fmpz_mod_poly_compose_mod(pow->pow + i, pow->pow + i - 1,
                                                         pow->pow + 1, f, ctx);
}
