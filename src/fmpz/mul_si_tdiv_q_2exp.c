/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Fredrik Johansson

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

void
fmpz_mul_si_tdiv_q_2exp(fmpz_t f, const fmpz_t g, slong x, ulong exp)
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
           if ((c2 ^ x) < WORD(0))
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
           if ((c2 ^ x) < WORD(0))
               fmpz_neg(f, f);
           return;
       }
       else /* result takes two limbs */
       {
           __mpz_struct * mf = _fmpz_promote(f);

           /* two limbs, least significant first, native endian, no
nails, stored in prod */
           mpz_import(mf, 2, -1, sizeof(mp_limb_t), 0, 0, prod);
           if ((c2 ^ x) < WORD(0))
               mpz_neg(mf, mf);
       }
   }
   else /* c2 is large */
   {
       __mpz_struct * mf = _fmpz_promote(f); /* ok without val as
            if aliased both are large */
       flint_mpz_mul_si(mf, COEFF_TO_PTR(c2), x);
       mpz_tdiv_q_2exp(mf, mf, exp);
       _fmpz_demote_val(f);  /* value may be small */
   }
}
