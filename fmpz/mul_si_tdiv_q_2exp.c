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

    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_mul_si_tdiv_q_2exp(fmpz_t f, const fmpz_t g, len_t x, ulong exp)
{
   fmpz c2 = *g;

   if (x == 0)
   {
       fmpz_zero(f);
       return;
   }
   else if (!COEFF_IS_MPZ(c2)) /* c2 is small */
   {
       mp_limb_t prod[2];
       mp_limb_t uc2;
       mp_limb_t ux;

       if (exp >= 2 * FLINT_BITS)
       {
           fmpz_zero(f);
           return;
       }

       uc2 = FLINT_ABS(c2);
       ux = FLINT_ABS(x);

       umul_ppmm(prod[1], prod[0], uc2, ux);

       if (exp >= FLINT_BITS)
       {
           fmpz_set_ui(f, prod[1] >> (exp - FLINT_BITS));
           if ((c2 ^ x) < 0L)
               fmpz_neg(f, f);
           return;
       }

       if (exp != 0)
       {
           prod[0] = (prod[1] << (FLINT_BITS - exp)) | (prod[0] >> exp);
           prod[1] >>= exp;
       }

       if (!prod[1])
       {
           fmpz_set_ui(f, prod[0]);
           if ((c2 ^ x) < 0L)
               fmpz_neg(f, f);
           return;
       }
       else /* result takes two limbs */
       {
           __mpz_struct *mpz_ptr = _fmpz_promote(f);

           /* two limbs, least significant first, native endian, no
nails, stored in prod */
           mpz_import(mpz_ptr, 2, -1, sizeof(mp_limb_t), 0, 0, prod);
           if ((c2 ^ x) < 0L)
               mpz_neg(mpz_ptr, mpz_ptr);
       }
   }
   else /* c2 is large */
   {
       __mpz_struct *mpz_ptr = _fmpz_promote(f); /* ok without val as
            if aliased both are large */
       mpz_mul_si(mpz_ptr, COEFF_TO_PTR(c2), x);
       mpz_tdiv_q_2exp(mpz_ptr, mpz_ptr, exp);
       _fmpz_demote_val(f);  /* value may be small */
   }
}
