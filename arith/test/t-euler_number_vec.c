/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
    slong k, n;

    FLINT_TEST_INIT(state);

    flint_printf("euler_number_vec....");
    fflush(stdout);

    for (n = 2; n <= 3000; n += (n<100) ? 2 : n/3)
    {
        n += n % 2;
        r = _fmpz_vec_init(n + 1);
        fmpz_init(s);
        fmpz_init(t);

        arith_euler_number_vec(r, n + 1);

        /* sum binomial(n,k) E_k = 0 */
        fmpz_set_ui(t, UWORD(1));
        for (k = 0; k <= n; k++)
        {
            fmpz_addmul(s, r + k, t);
            fmpz_mul_ui(t, t, n - k);
            fmpz_divexact_ui(t, t, k + 1);
        }

        if (!fmpz_is_zero(s))
        {
            flint_printf("ERROR: sum over 0,...,n = %wd\n", n);
            _fmpz_vec_print(r, n + 1);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(s);
        fmpz_clear(t);
        _fmpz_vec_clear(r, n + 1);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
