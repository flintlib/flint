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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_sparse.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);
    flint_randseed(state, UWORD(3), UWORD(2)); /* TODO remove */

    flint_printf("interp....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
      fmpz_sparse_bp_interp_t f;
      fmpz_sparse_t a, b, c;
      fmpz_t d, h;

      /*random fmpz*/
      fmpz_init(d);
      fmpz_init(h);
      /* fmpz_randtest(d, state, 100); TODO too hard */
      fmpz_randtest(d, state, 6);

      /*random fmpz_sparse*/
      fmpz_sparse_init(a);
      fmpz_sparse_init(b);
      fmpz_sparse_init(c);
      /* fmpz_sparse_randtest(a, state, n_randint(state, 100), d, 200); TODO too hard */
      fmpz_sparse_randtest(a, state, n_randint(state, 4), d, 7);

      /*calculate height of coefficients*/
      fmpz_sparse_height(h, a);
      fmpz_mul_ui(h, h, UWORD(2));
      
      /*initialize interpolation struct and eval random fmpz_sparse*/
      fmpz_sparse_bp_interp_init(f, a->length, h, d, state);
      fmpz_sparse_bp_interp_eval(f, a);
      
      /*b gets result*/
      fmpz_sparse_bp_interp(b, f);

      /*c gets result*/
      fmpz_sparse_bp_interp(c, f);

      /*check that a == b and b == c*/
      result = (fmpz_sparse_equal(a, b) && fmpz_sparse_equal(b,c));
      if (!result)
      {
          flint_printf("FAIL:\n");
          fmpz_sparse_print(a), flint_printf("\n\n");
          fmpz_sparse_print(b), flint_printf("\n\n");
          fmpz_sparse_print(c), flint_printf("\n\n");
          abort();
      }

      fmpz_clear(d);
      fmpz_clear(h);
      fmpz_sparse_clear(a);
      fmpz_sparse_clear(b);
      fmpz_sparse_clear(c);
      fmpz_sparse_bp_interp_clear(f);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
