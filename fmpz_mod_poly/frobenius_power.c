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
fmpz_mod_poly_frobenius_power(fmpz_mod_poly_t res,
                            fmpz_mod_poly_frobenius_powers_2exp_t pow, 
                    const fmpz_mod_poly_t f, ulong m, const fmpz_mod_ctx_t ctx)
{
    slong i = 0;
    ulong bit;
    fmpz_mod_poly_struct * r;
    fmpz_mod_poly_t tr;

    if (res == f)
    {
       fmpz_mod_poly_init(tr, ctx);
       r = tr;
    }
    else
    {
       r = res;
    }

    /* res = x^(p^0) = x */
    if (m == 0)
    {
       fmpz_mod_poly_set_coeff_ui(r, 1, 1, ctx);
       fmpz_mod_poly_set_coeff_ui(r, 0, 0, ctx);
       _fmpz_mod_poly_set_length(r, 2);
       
       /* 
          This is safe wrt impossible inverses, because any zero divisors
          in the leading coefficient of f will have been found in the 
          precomp stage.
       */
       if (f->length <= 2)
          fmpz_mod_poly_rem(r, r, f, ctx);
    }
    else
    {
       /* first nonzero bit */
       while ((m & (WORD(1) << i)) == 0)
          i++;

       /* res = f^(p^(2^i)) */
       fmpz_mod_poly_set(r, pow->pow + i, ctx);
       m ^= (WORD(1) << i);

       while (m != 0)
       {
          i++;
       
          bit = (WORD(1) << i);
          if ((bit & m) != 0)
          {
             fmpz_mod_poly_compose_mod(r, pow->pow + i, r, f, ctx);
             m ^= bit;
          }
       }
    }

    if (res == f)
    {
       fmpz_mod_poly_swap(res, r, ctx);
       fmpz_mod_poly_clear(tr, ctx);
    }
}
