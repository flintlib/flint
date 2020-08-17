/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int fmpz_is_probabprime_lucas(const fmpz_t n)
{
   fmpz_t A, Q, D, t, m, Vm, Vm1;
   int res = 0;

   if (fmpz_cmp_ui(n, 1) <= 0)
      return 0;

   if (fmpz_is_even(n))
      return fmpz_cmp_ui(n, 2) == 0; /* ensure 2, n coprime */

   if (fmpz_is_square(n))
      return 0;

   fmpz_init(A);
   fmpz_init(Q);
   fmpz_init(D);
   fmpz_init(t);
   fmpz_init(m);
   fmpz_init(Vm);
   fmpz_init(Vm1);

   fmpz_set_si(D, WORD(-3)); /* -3 becomes 5 after first iteration */

   do {
      /* check D = 5, -7, 9, -11 such that (D/n) = -1 */
      do {
         if (fmpz_sgn(D) > 0)
            fmpz_add_ui(D, D, 2);
         else
            fmpz_sub_ui(D, D, 2);

         fmpz_neg(D, D);
      } while (fmpz_jacobi(D, n) != -1); /* this ensures D, n coprime */
   
      fmpz_sub_ui(t, D, 1);
      fmpz_neg(t, t);
      fmpz_tdiv_q_2exp(Q, t, 2);

      fmpz_gcd(t, Q, n); /* require Q, n coprime */
   } while (fmpz_equal(t, n)); 

   if (fmpz_is_one(t)) /* check no factor found */
   {
      fmpz_invmod(A, Q, n); /* A = P^2Q^-1 - 2 mod n, where P = 1 */
      fmpz_sub_ui(A, A, 2);
      if (fmpz_sgn(A) < 0)
         fmpz_add(A, A, n);

      fmpz_add_ui(m, n, 1); /* m = (n - jacobi(D/n))/2 = (n + 1)/2 */
      fmpz_tdiv_q_2exp(m, m, 1);
 
      fmpz_lucas_chain(Vm, Vm1, A, m, n);

      fmpz_mul(Vm, Vm, A);
      fmpz_submul_ui(Vm, Vm1, 2);

      res = fmpz_divisible(Vm, n);
   }

   fmpz_clear(A);
   fmpz_clear(Q);
   fmpz_clear(D);
   fmpz_clear(t);
   fmpz_clear(m);
   fmpz_clear(Vm);
   fmpz_clear(Vm1);

   return res;
}
