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

    Authored 2015 by Daniel S. Roche and A. Whitman Groves;
    US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"

void fmpz_sparse_randtest(fmpz_sparse_t res, flint_rand_t state, 
    slong terms, const fmpz_t degree, mp_bitcnt_t bits)
{
    slong i, k;
    mp_limb_t j;
    fmpz_t abs;
    fmpz_sparse_zero(res);
    
    fmpz_init(abs);
      
    /* zero terms => 0 */
    if(terms == 0)
      return;

    /* zero degree => random number with 0 exponent */
    if(fmpz_is_zero(degree))
    {
      fmpz_randtest(res->coeffs, state, bits);
      return;
    }

    /*if abs(degree) < terms => terms == degree */
    fmpz_abs(abs, degree);
    if(fmpz_cmp_si(abs, terms) < 0)
    {
      terms = fmpz_get_si(abs);
    }

    if(fmpz_cmp_si(degree, 0) < 0)
      k = 1;
    else
      k = 0;

    _fmpz_sparse_reserve(res, terms);
    for (i=0; i<terms; ++i) {
        fmpz_init(res->coeffs+i);
        fmpz_init(res->expons+i);
        fmpz_randbits(res->coeffs+i, state, bits);
        fmpz_randm(res->expons+i, state, degree);
        if(k)
        {
          j = n_randint(state, 2);
          if(j)
            fmpz_mul_si(res->expons+i, res->expons+i, -1);
        }
    }
    res->length = terms;
    _fmpz_sparse_normalise(res);
}
