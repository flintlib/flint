/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011, 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h> /* qsort */
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

int _compare_poly_lengths(const void * a, const void * b)
{
   const fmpz_poly_struct * p = (fmpz_poly_struct *) a;
   const fmpz_poly_struct * q = (fmpz_poly_struct *) b;

   return p->length - q->length;
}

int fmpz_poly_factor_van_hoeij_check_if_solved(fmpz_mat_t M,
          fmpz_poly_factor_t final_fac, fmpz_poly_factor_t lifted_fac,
                           const fmpz_poly_t f, fmpz_t P, slong exp, fmpz_t lc)
{
   fmpz_poly_factor_t trial_factors;
   fmpz_poly_t prod, q, f_copy;
   fmpz_t temp_lc;
   fmpz_mat_t U;
   nmod_poly_t f2, g2, rem;
   int num_facs, res = 0;
   slong i, j, r;
   slong * part;
 
   r = lifted_fac->num;
   
   part = (slong *) flint_calloc(r, sizeof(slong));

   fmpz_poly_factor_init(trial_factors);
   fmpz_poly_init(prod);
   fmpz_poly_init(q);
   fmpz_poly_init(f_copy);
   fmpz_mat_window_init(U, M, 0, 0, M->r, r);
   fmpz_init(temp_lc);

   nmod_poly_init(f2, 2);
   nmod_poly_init(g2, 2);
   nmod_poly_init(rem, 2);

   if ((num_facs = fmpz_mat_col_partition(part, U, 1)) == 0 || num_facs > r)
      goto cleanup;

   if (num_facs == 1)
   {
      /* f is irreducible */

      fmpz_poly_factor_insert(final_fac, f, exp);
      res = 1;
      goto cleanup;
   }

   /* 
      there is a potential 0-1 basis, so make the potential factors 
   */

   fmpz_set(temp_lc, lc);

   for (i = 1; i <= num_facs; i++)
   {
      fmpz_poly_set_fmpz(prod, temp_lc);

      for (j = 0; j < r; j++)
      {
         if (part[j] == i)
         {
            fmpz_poly_mul(prod, prod, lifted_fac->p + j);
            fmpz_poly_scalar_smod_fmpz(prod, prod, P);
         }
      }

      fmpz_poly_content(temp_lc, prod);
      fmpz_abs(temp_lc, temp_lc);
      fmpz_poly_scalar_divexact_fmpz(prod, prod, temp_lc);

      fmpz_poly_factor_insert(trial_factors, prod, 1);
   }

   /* sort factors by length */

   qsort(trial_factors->p, trial_factors->num, 
                            sizeof(fmpz_poly_struct), _compare_poly_lengths);

   /* trial divide potential factors */

   fmpz_poly_set(f_copy, f);

   for (i = 0; i < trial_factors->num && num_facs > 1; i++)
   {
      /* check if the polynomial divides mod 2 */
      fmpz_poly_get_nmod_poly(f2, f_copy);
      fmpz_poly_get_nmod_poly(g2, trial_factors->p + i);

      nmod_poly_rem(rem, f2, g2);
   
      if (nmod_poly_is_zero(rem) && fmpz_poly_divides(q, f_copy, trial_factors->p + i))
      {
         fmpz_poly_swap(q, f_copy);
         num_facs--;
      } else
         goto cleanup;
   }

   /* if we found all the factors, insert them */

   if (num_facs == 1)
   {
      for (j = 0; j < i; j++)
         fmpz_poly_factor_insert(final_fac, trial_factors->p + j, exp);
      
      fmpz_poly_factor_insert(final_fac, f_copy, exp);

      res = 1; /* we factorised f */
   }

cleanup:

   nmod_poly_clear(f2);
   nmod_poly_clear(g2);
   nmod_poly_clear(rem);
   fmpz_clear(temp_lc);
   fmpz_poly_clear(q);
   fmpz_poly_clear(f_copy);
   fmpz_poly_clear(prod);
   fmpz_poly_factor_clear(trial_factors);
   fmpz_mat_window_clear(U);
   flint_free(part);

   return res;
}
