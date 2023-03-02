/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_discriminant(fmpz_t d, const fmpz *poly, slong len, const fmpz_t mod)
{
   fmpz *der = _fmpz_vec_init(len - 1);
   slong dlen = len - 1, exp;
   fmpz_t pow;
   
   _fmpz_mod_poly_derivative(der, poly, len, mod);
   FMPZ_VEC_NORM(der, dlen);

   if (dlen == 0)
   {
      fmpz_zero(d);
   } else
   {
      fmpz_init(pow);

      _fmpz_mod_poly_resultant(d, poly, len, der, dlen, mod);
      exp = len - dlen - 2;
      
      if (exp >= 0)
         fmpz_powm_ui(pow, poly + len - 1, exp, mod);
      else
         fmpz_invmod(pow, poly + len - 1, mod);
      
      fmpz_mul(d, d, pow);
      fmpz_mod(d, d, mod);
      
      if ((len & 3) == 0 || (len & 3) == 3) /* degree is not 0, 1 mod 4 */
         fmpz_negmod(d, d, mod);

      fmpz_clear(pow);
   }
   
   _fmpz_vec_clear(der, len - 1);
}

void fmpz_mod_poly_discriminant(fmpz_t d, const fmpz_mod_poly_t f,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong len = f->length;
    
    if (len <= 1)
        fmpz_zero(d);
    else
       _fmpz_mod_poly_discriminant(d, f->coeffs, len, fmpz_mod_ctx_modulus(ctx));
}

