/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011, 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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

/*
   resize the matrix M to num_rows, where num_rows <= M->r
   the algorithm assumes rows have been permuted in memory by LLL
*/
void fmpz_mat_van_hoeij_resize_matrix(fmpz_mat_t M, slong num_rows)
{
   slong i, j;
   fmpz ** empty_rows;
   slong num_empty = 0;
   fmpz * end;
   TMP_INIT;

   if (M->r == num_rows)
      return; /* nothing to be done */

   TMP_START;

   /* space to record empty rows within bounds of new matrix */
   empty_rows = (fmpz **) TMP_ALLOC(M->r*sizeof(fmpz *));

   /* end of new matrix in memory */
   end = M->entries + num_rows*M->c;

   /* clear rows that aren't needed */
   for (i = num_rows; i < M->r; i++)
   {
      _fmpz_vec_zero(M->rows[i], M->c);
      
      /* this row can be repopulated */
      if (M->rows[i] < end)
         empty_rows[num_empty++] = M->rows[i];
   }

   for (i = 0; i < num_rows; i++)
   {
      if (M->rows[i] >= end) /* this row must be moved back to empty spot */
      {
         fmpz * old_row = M->rows[i];
         fmpz * new_row = empty_rows[--num_empty];

         for (j = 0; j < M->c; j++)
            fmpz_swap(old_row + j, new_row + j);

         M->rows[i] = new_row;
      }
   }

   M->r = num_rows;

   TMP_END;   
}

void fmpz_poly_factor_van_hoeij(fmpz_poly_factor_t final_fac,
         const nmod_poly_factor_t fac, const fmpz_poly_t f, slong exp, ulong p)
{
   fmpz_poly_factor_t lifted_fac; 
   fmpz_mat_t M;
   fmpz_t fp, P, B, lc, bound_sum;
   slong i, r = fac->num;
   slong bit_r = FLINT_MAX(r, 20);
   slong U_exp, a, next_col, num_data_cols, num_coeffs;
   slong prev_num_coeffs, prev_exp, N, worst_exp, num_rows;
   ulong sqN;
   fmpz_poly_t * v, * w;
   fmpz_mat_t col, data;
   slong * link;
   int hensel_loops, do_lll;
   fmpz_lll_t fl;

   /* set to identity */
   fmpz_mat_init(M, r, r);

   for (i = 0; i < r; i++)
      fmpz_set_ui(M->rows[i] + i, 1);

   /* we prescale the identity matrix by 2^U_exp */
   U_exp = (slong) FLINT_BIT_COUNT(bit_r);
   
   fmpz_mat_scalar_mul_2exp(M, M, U_exp);
 
   /* compute Mignotte bound */
   fmpz_init(B);
   
   fmpz_poly_factor_mignotte(B, f);
   /*
      bound adjustment, we multiply true factors (which might be
      monic) by the leading coefficient of f in the implementation
      below
   */
   fmpz_mul(B, B, f->coeffs + f->length - 1);
   fmpz_abs(B, B);
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

   while (!fmpz_poly_factor_van_hoeij_check_if_solved(M, final_fac, lifted_fac, f, P, exp, lc))
   {
      if (hensel_loops < 3 && 3*r > N + 1)
         num_coeffs = r > 200 ? 50 : 30;
      else
         num_coeffs = 10;

      num_coeffs = FLINT_MIN(num_coeffs, (N + 1)/2);
      prev_num_coeffs = 0;

      do {
         fmpz_mat_init(data, r + 1, 2*num_coeffs);
         num_data_cols = _fmpz_poly_factor_CLD_mat(data, f, lifted_fac, P, num_coeffs);

         for (next_col = prev_num_coeffs; next_col < num_data_cols - prev_num_coeffs; next_col++)
         {
            /* we alternate taking columns from the right and left */
            slong alt_col, diff = next_col - prev_num_coeffs;
            
            if ((diff % 2) == 0)
               alt_col = prev_num_coeffs + diff/2;
            else
               alt_col = num_data_cols - prev_num_coeffs - (diff + 1)/2;

            fmpz_mul_ui(bound_sum, data->rows[r] + alt_col, sqN);
            worst_exp = fmpz_bits(bound_sum);   
            
            for (i = 0; i < r; i++)
               fmpz_set(col->rows[i], data->rows[i] + alt_col);

            do_lll = fmpz_mat_next_col_van_hoeij(M, P, col, worst_exp, U_exp);

            if (do_lll)
            {
               num_rows = fmpz_lll_wrapper_with_removal_knapsack(M, NULL, B, fl);

               fmpz_mat_van_hoeij_resize_matrix(M, num_rows);

               if (fmpz_poly_factor_van_hoeij_check_if_solved(M, final_fac, lifted_fac, f, P, exp, lc))
               {
                  fmpz_mat_clear(data);
                  goto cleanup;
               }
            }
         }

         prev_num_coeffs = num_coeffs;
         num_coeffs = FLINT_MIN(2*num_coeffs, (N + 1)/2);
         fmpz_mat_clear(data);

      } while (num_coeffs != prev_num_coeffs);
   
      hensel_loops++;
      
      prev_exp = _fmpz_poly_hensel_continue_lift(lifted_fac, link, v, w, f, prev_exp, a, 2*a, fp);

      a = 2*a;
      fmpz_pow_ui(P, fp, a);
   }

cleanup:

   fmpz_clear(lc);
   fmpz_clear(fp);
   fmpz_clear(P);
   fmpz_clear(B);
   fmpz_mat_clear(col);
   fmpz_mat_clear(M);
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
}
