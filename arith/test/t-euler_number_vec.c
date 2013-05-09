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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_vec.h"


int main()
{
    fmpz * r;
    fmpz_t s, t;
    len_t k, n;

    printf("euler_number_vec....");
    fflush(stdout);

    for (n = 2; n <= 3000; n += (n<100) ? 2 : n/3)
    {
        n += n % 2;
        r = _fmpz_vec_init(n + 1);
        fmpz_init(s);
        fmpz_init(t);

        arith_euler_number_vec(r, n + 1);

        /* sum binomial(n,k) E_k = 0 */
        fmpz_set_ui(t, 1UL);
        for (k = 0; k <= n; k++)
        {
            fmpz_addmul(s, r + k, t);
            fmpz_mul_ui(t, t, n - k);
            fmpz_divexact_ui(t, t, k + 1);
        }

        if (!fmpz_is_zero(s))
        {
            printf("ERROR: sum over 0,...,n = %ld\n", n);
            _fmpz_vec_print(r, n + 1);
            abort();
        }

        fmpz_clear(s);
        fmpz_clear(t);
        _fmpz_vec_clear(r, n + 1);
    }

    mpfr_free_cache();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
