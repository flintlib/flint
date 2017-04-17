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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fq.h"
#include "fq_mat.h"
#include "fmpz_mod_poly.h"
#include "fmpz_sparse.h"

void fmpz_sparse_bp_interp(fmpz_sparse_t res,
    const fmpz_sparse_bp_interp_t evals)
{
  /* Berlekamp-Massey to find the generator polynomial
   * Compute roots of generator poly, along with discrete logs
   * Compute coefficients from transposed Vandermode
   */

  fq_ctx_t Fq;
  fq_mat_t M;
  fq_mat_t b;
  fmpz_mod_poly_t G;
  fmpz_t temp_fmpz;
  fmpz * roots;
  fmpz * expons;
  slong T = evals->length / 2;
  slong t;
  slong i, j;
  slong temp_slong;
  slong * perm;
  const fmpz * w = evals->sample_points + 1;

  fmpz_sparse_zero(res);

  if (T == 0) return;

  fmpz_init(temp_fmpz);
  
  /* build finite field */
  fq_ctx_init(Fq, evals->q, WORD(1), "x");


  /* Build and solve Hankel matrix for Berlekamp-Massey */

  fq_mat_init(M, T, T, Fq);

  for (i=0; i<T; ++i)
  {
    for (j=0; j<T; ++j) 
    {
      fq_set_fmpz(fq_mat_entry(M, i, j), evals->evaluations+(i+j), Fq);
    }
  }

  perm = flint_calloc(T, sizeof(slong));

  temp_slong = fq_mat_lu(perm, M, 1, Fq);
  FLINT_ASSERT (temp_slong);

  fq_mat_init(b, T, 1, Fq);
  for (i=0; i<T; ++i)
  {
    fq_set_fmpz(fq_mat_entry(b, perm[i], 0), evals->evaluations+(i+T), Fq);
  }

  fq_mat_solve_tril(b, M, b, 1, Fq);
  fq_mat_solve_triu(b, M, b, 0, Fq);


  /* compute actual sparsity t */

  t = T;
  for (i=0; i<T && fq_is_zero(fq_mat_entry(b,i,0),Fq); ++i) --t;


  /* create Prony polynomial */

  fmpz_mod_poly_init2(G, evals->q, t+1);
  fmpz_mod_poly_set_coeff_ui(G, t, UWORD(1));

  for (i=0; i<t; ++i) 
  {
    if (fmpz_poly_length(fq_mat_entry(b,i+(T-t),0)) > 0) 
    {
      FLINT_ASSERT(fmpz_poly_length(fq_mat_entry(b,i+(T-t),0)) == WORD(1));
      fmpz_poly_get_coeff_fmpz(temp_fmpz, fq_mat_entry(b,i+(T-t),0), 0);
      fmpz_neg(temp_fmpz, temp_fmpz);
      fmpz_mod_poly_set_coeff_fmpz(G, i, temp_fmpz);
    }
  }

  fq_mat_clear(b, Fq);
  fq_mat_clear(M, Fq);


  /* find roots of Prony polynomial, and their orders */

  roots = _fmpz_vec_init(t);
  expons = _fmpz_vec_init(t);
  fmpz_set_ui(roots+0, UWORD(1));
  fmpz_set_ui(expons+0, UWORD(0));

  i = 0;
  while (1) 
  {
    while (1) 
    {
      fmpz_mod_poly_evaluate_fmpz(temp_fmpz, G, roots+i);
      if (fmpz_is_zero(temp_fmpz)) break;
      fmpz_mul(roots+i, roots+i, w);
      fmpz_add_ui(expons+i, expons+i, UWORD(1));
    }

    if (i == t-1) break;

    ++i;
    fmpz_mul(roots+i, roots+(i-1), w);
    fmpz_add_ui(expons+i, expons+(i-1), UWORD(1));
  }


  /* solve transposed Vandermode to get coeffs */

  fq_mat_init(M, t, t, Fq);

  for (i=0; i<t; ++i)
  {
    fq_set_ui(fq_mat_entry(M,0,i), UWORD(1), Fq);
  }

  for (i=1; i<t; ++i)
  {
    for (j=0; j<t; ++j) 
    {
      fmpz_mul(temp_fmpz, 
          fmpz_poly_get_coeff_ptr(fq_mat_entry(M,i-1,j), 0),
          roots+j);
      fq_set_fmpz(fq_mat_entry(M,i,j), temp_fmpz, Fq);
    }
  }

  temp_slong = fq_mat_lu(perm, M, 1, Fq);
  FLINT_ASSERT (temp_slong);

#ifdef WANT_ASSERT
  for (i=0; i<T; ++i) FLINT_ASSERT(perm[i] == i);
#endif

  fq_mat_init(b, t, 1, Fq);
  for (i=0; i<T; ++i)
  {
    fq_set_fmpz(fq_mat_entry(b,i,0), evals->evaluations+i, Fq);
  }

  fq_mat_solve_tril(b, M, b, 1, Fq);
  fq_mat_solve_triu(b, M, b, 0, Fq);


  /* set coeffs and expons of the actual polynomial */
  
  if (evals->laurent)
  {
    fmpz_fdiv_q_2exp(temp_fmpz, evals->order, UWORD(1));

    for (i=0; i<t; ++i)
    {
      if (fmpz_cmp(expons+i, temp_fmpz) > 0)
      {
        fmpz_sub(expons+i, expons+i, evals->order);
      }
    }
  }

  for (i=0; i<t; ++i)
  {
    fq_struct * coeff = fq_mat_entry(b,i,0);

    if (!fq_is_zero(coeff, Fq))
    {
      FLINT_ASSERT(fmpz_poly_length(coeff) == UWORD(1));
      fmpz_mods(temp_fmpz, fmpz_poly_get_coeff_ptr(coeff,0), evals->q);
      fmpz_sparse_set_coeff(res, temp_fmpz, expons+i);
    }
  }

  fq_mat_clear(b, Fq);
  fq_mat_clear(M, Fq);

  /* clean-up */
  _fmpz_vec_clear(roots, t);
  _fmpz_vec_clear(expons, t);
  fmpz_mod_poly_clear(G);
  flint_free(perm);
  fq_ctx_clear(Fq);
  fmpz_clear(temp_fmpz);
}
