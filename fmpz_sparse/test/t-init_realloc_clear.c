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

  Authored 2015 by A. Whitman Groves; US Government work in the public domain.

******************************************************************************/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_sparse.h"
#include "ulong_extras.h"

int
main(void)
{
    /**int i;
    FLINT_TEST_INIT(state);

    flint_printf("init/init2/realloc/clear....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a;

        fmpz_sparse_init2(a, n_randint(state, 100));
        fmpz_sparse_clear(a);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a;

        fmpz_sparse_init2(a, n_randint(state, 100));
        fmpz_sparse_realloc(a, n_randint(state, 100));
        fmpz_sparse_clear(a);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a;
        fmpz_t b;

        fmpz_init(b);
        fmpz_randtest_not_zero(b, state, 100);

        fmpz_sparse_init(a);
        fmpz_sparse_randtest(a, state, n_randint(state, 100), b, 200);
        
        fmpz_sparse_clear(a);
        fmpz_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");*/
    return 0;
}
