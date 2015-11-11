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

#include "fmpz_sparse.h"

void fmpz_sparse_bp_interp_eval(fmpz_sparse_bp_interp_t res,
    const fmpz_sparse_t poly)
{
  slong i;

  /* TODO can be faster by precomputing evaluations of monomials. */
  for (i=0; i<res->length; ++i)
  {
    fmpz_sparse_evaluate_mod(res->evaluations+i,
        poly, res->sample_points+i, res->q);
  }
}

void fmpz_sparse_bp_interp_mul(fmpz_sparse_bp_interp_t res,
    const fmpz_sparse_t poly)
{
  slong i;
  fmpz_t temp;
  
  fmpz_init2(temp, fmpz_size(res->q));

  for (i=0; i<res->length; ++i)
  {
    fmpz_sparse_evaluate_mod(temp, poly, res->sample_points+i, res->q);
    fmpz_mul(res->evaluations+i, res->evaluations+i, temp);
    fmpz_mod(res->evaluations+i, res->evaluations+i, res->q);
  }

  fmpz_clear(temp);
}


void fmpz_sparse_bp_interp_add(fmpz_sparse_bp_interp_t res,
    const fmpz_t c, const fmpz_sparse_t poly)
{
  slong i;
  fmpz_t temp;
  
  fmpz_init2(temp, fmpz_size(res->q));

  for (i=0; i<res->length; ++i)
  {
    fmpz_sparse_evaluate_mod(temp, poly, res->sample_points+i, res->q);
    fmpz_addmul(res->evaluations+i, c, temp);
    fmpz_mod(res->evaluations+i, res->evaluations+i, res->q);
  }

  fmpz_clear(temp);
}

void fmpz_sparse_bp_interp_pow(fmpz_sparse_bp_interp_t res, ulong pow)
{
  slong i;
  
  for (i=0; i<res->length; ++i)
  {
    fmpz_powm_ui(res->evaluations+i, res->evaluations+i, pow, res->q);
  }
}

