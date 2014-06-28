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

    Copyright (C) 2013, 2014 William Hart
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

fmpz ** 
_fmpz_mod_poly_powers_precompute(const fmpz * B, slong len, const fmpz_t mod)
{
   slong i;
   fmpz ** powers = flint_malloc(sizeof(fmpz *)*(2*len - 1));
   fmpz_mod_poly_t pow, p;
   fmpz_t c, inv;
   int monic = fmpz_is_one(B + len - 1);

   fmpz_mod_poly_init2(pow, mod, len);
   fmpz_mod_poly_set_ui(pow, 1);
   fmpz_mod_poly_init2(p, mod, len - 1);
   fmpz_init(c);
   fmpz_init(inv);

   if (!monic)
      fmpz_invmod(inv, B + len - 1, mod);

   for (i = 0; i < 2*len - 1; i++)
   {
      powers[i] = _fmpz_vec_init(len - 1);

      if (pow->length == len) /* reduce pow mod B */
      {
         if (monic)
            _fmpz_mod_poly_scalar_mul_fmpz(p->coeffs, B, len - 1, pow->coeffs + len - 1, mod);
         else
         {
            fmpz_mul(c, pow->coeffs + len - 1, inv);
            fmpz_mod(c, c, mod);
            _fmpz_mod_poly_scalar_mul_fmpz(p->coeffs, B, len - 1, c, mod);
         }
         
         _fmpz_mod_poly_set_length(p, len - 1);
         _fmpz_mod_poly_normalise(p);
         fmpz_mod_poly_sub(pow, pow, p);
         _fmpz_mod_poly_set_length(pow, len - 1);
         _fmpz_mod_poly_normalise(pow);
      }

      _fmpz_vec_set(powers[i], pow->coeffs, len - 1);

      fmpz_mod_poly_shift_left(pow, pow, 1);
   }

   fmpz_clear(c);
   fmpz_clear(inv);
   fmpz_mod_poly_clear(pow);
   fmpz_mod_poly_clear(p);

   return powers;
}

void fmpz_mod_poly_powers_precompute(fmpz_mod_poly_powers_precomp_t pinv, 
                                                          fmpz_mod_poly_t poly)
{
    if (poly->length == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_powers_precompute). Division by zero.\n");
        abort();
    }

    pinv->powers = _fmpz_mod_poly_powers_precompute(poly->coeffs, poly->length, &poly->p);
    pinv->len = poly->length;
}
