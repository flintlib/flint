/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011, 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "fmpz_poly.h"

slong _heuristic_van_hoeij_starting_precision(const fmpz_poly_t f, 
                                                            slong r, ulong p)
{
   fmpz_t lead_b, trail_b;
   slong min_b, a;
   
   fmpz_init(lead_b);
   fmpz_init(trail_b);
   
   fmpz_poly_CLD_bound(lead_b, f, f->length - 2);
   fmpz_poly_CLD_bound(trail_b, f, 0);

   min_b = FLINT_MIN(fmpz_bits(lead_b), fmpz_bits(trail_b));

   a = (slong) ((2.5*r + min_b)*log(2) + log(f->length)/2.0)/log(p);

   fmpz_clear(trail_b);
   fmpz_clear(lead_b);

   return a;
}

void fmpz_poly_factor_van_hoeij(fmpz_poly_factor_t final_fac,
         const nmod_poly_factor_t fac, const fmpz_poly_t f, slong exp, ulong p)
{
   fmpz_poly_factor_t lifted_fac; 
   fmpz_mat_t M;
   fmpz_t fp, P, B, lc, bound_sum;
   slong i, r = fac->num;
   slong bit_r = FLINT_MAX(r, 20);
   slong U_exp, a, next_col, num_coeffs, prev_num_coeffs, prev_exp, N, worst_exp;
   ulong sqN;
   fmpz_poly_t * v, * w;
   fmpz_mat_t col, data, U;
   slong * link;
   int hensel_loops;
   fmpz_lll_t fl;

   printf("made it 1\n");

   /* set to identity */
   fmpz_mat_init(M, r, r);

   for (i = 0; i < r; i++)
      fmpz_set_ui(M->rows[i] + i, 1);

   /* we prescale the identity matrix by 2^U_exp */
   U_exp = bit_r;
   
   fmpz_mat_scalar_mul_2exp(M, M, U_exp);
 
   /* compute Mignotte bound */
   fmpz_init(B);
   
   fmpz_poly_factor_mignotte(B, f);
   fmpz_mul_ui(B, B, 2);
   fmpz_add_ui(B, B, 1);
   a = fmpz_clog_ui(B, p);
                
   /* compute heuristic starting precision */
   a = FLINT_MIN(a, _heuristic_van_hoeij_starting_precision(f, r, p));

   /* start Hensel lift */
   fmpz_poly_factor_init(lifted_fac);

   v = flint_malloc((2*r - 2)*sizeof(fmpz_poly_t));
   w = flint_malloc((2*r - 2)*sizeof(fmpz_poly_t));
   link = flint_malloc((2*r - 2)*sizeof(slong));

   for (i = 0; i < 2*r - 2; i++)
   {
      fmpz_poly_init(v[i]);
      fmpz_poly_init(w[i]);
   }

   prev_exp = _fmpz_poly_hensel_start_lift(lifted_fac, link, v, w, f, fac, a);

   /* compute bound */
   fmpz_set_ui(B, r + 1);
   fmpz_mul_2exp(B, B, 2*U_exp);
   fmpz_sqrt(B, B);

   /* compute leading coefficient */
   N = f->length - 1;
   sqN = (ulong) sqrt(N);
   fmpz_init(lc);
   fmpz_set(lc, f->coeffs + N);

   /* main hensel loop */
   hensel_loops = 0;
   fmpz_init(P);
   fmpz_init(fp);
   fmpz_set_ui(fp, p);
   fmpz_pow_ui(P, fp, a);
         
   fmpz_init(bound_sum);
   fmpz_mat_init(col, r, 1);
   fmpz_lll_context_init_default(fl);

   printf("made it 2\n");

   while (!fmpz_poly_factor_van_hoeij_check_if_solved(M, final_fac, lifted_fac, f, P, exp, lc))
   {
      printf("made it 3\n");

      if (hensel_loops < 3 && 3*r > N + 1)
         num_coeffs = r > 200 ? 50 : 30;
      else
         num_coeffs = 10;

      num_coeffs = FLINT_MIN(num_coeffs, (N + 1)/2);
      prev_num_coeffs = 0;

      printf("made it 4\n");

      do {
         fmpz_mat_init(data, r + 1, 2*num_coeffs);
         _fmpz_poly_factor_CLD_mat(data, f, lifted_fac, P, num_coeffs);

         for (next_col = prev_num_coeffs; next_col < 2*num_coeffs - prev_num_coeffs; next_col++)
         {
            /* we alternate taking columns from the right and left */
            slong alt_col, diff = next_col - prev_num_coeffs;
            
            if ((diff % 2) == 0)
               alt_col = prev_num_coeffs + diff/2;
            else
               alt_col = 2*num_coeffs - prev_num_coeffs - (diff + 1)/2;

            fmpz_mul_ui(bound_sum, data->rows[r] + alt_col, sqN);
            worst_exp = fmpz_bits(bound_sum);   
            
            for (i = 0; i < r; i++)
               fmpz_set(col->rows[i], data->rows[i] + alt_col);
            
            fmpz_mat_next_col_van_hoeij(M, P, col, worst_exp, U_exp);

            fmpz_mat_init(U, M->r, M->r);
            fmpz_lll_wrapper_with_removal_knapsack(M, U, B, fl);
            fmpz_mat_clear(U);

            printf("made it 5\n");

            if (fmpz_poly_factor_van_hoeij_check_if_solved(M, final_fac, lifted_fac, f, P, exp, lc))
            {
               fmpz_mat_clear(data);
               goto cleanup;
            }
            printf("made it 6\n");
         }

         prev_num_coeffs = num_coeffs;
         num_coeffs = FLINT_MIN(2*num_coeffs, (N + 1)/2);
         fmpz_mat_clear(data);

         printf("made it 7\n");
      } while (num_coeffs != prev_num_coeffs);
   
      printf("made it 8\n");

      hensel_loops++;

      fmpz_mat_clear(data);
      
      printf("made it 9\n");

      prev_exp = _fmpz_poly_hensel_continue_lift(lifted_fac, link, v, w, f, prev_exp, a, 2*a, fp);

      a = 2*a;
      fmpz_pow_ui(P, fp, a);

      printf("made it 10\n");
   }
   

cleanup:
   printf("made it 11\n");

   fmpz_clear(lc);
   fmpz_clear(fp);
   fmpz_clear(P);
   fmpz_clear(B);
   fmpz_mat_clear(col);
   fmpz_clear(bound_sum);
   fmpz_poly_factor_clear(lifted_fac);
  
   for (i = 0; i < 2*r - 2; i++)
   {
      fmpz_poly_clear(v[i]);
      fmpz_poly_clear(w[i]);
   }
   flint_free(v);
   flint_free(w);
   flint_free(link);

   printf("made it 12\n");
}
