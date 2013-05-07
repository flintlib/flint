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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

static mp_limb_t
refimpl(long * exp, const fmpz_t x, int bits)
{
    fmpz_t t;
    long xbits;
    mp_limb_t m;

    xbits = fmpz_bits(x);

    fmpz_init(t);
    fmpz_abs(t, x);

    if (xbits >= bits)
        fmpz_tdiv_q_2exp(t, t, xbits - bits);
    else
        fmpz_mul_2exp(t, t, bits - xbits);

    m = fmpz_get_ui(t);
    fmpz_clear(t);

    *exp = xbits - bits;

    return m;
}

int
main(void)
{
    long iter;
    flint_rand_t state;

    printf("abs_lbound_ui_2exp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_t x;
        long bits;
        long exp, yexp;
        mp_limb_t yman, man;

        fmpz_init(x);
        fmpz_randtest_not_zero(x, state, 1 + n_randint(state, 400));

        bits = 1 + n_randint(state, FLINT_BITS - 1);

        yman = refimpl(&yexp, x, bits);
        man = fmpz_abs_lbound_ui_2exp(&exp, x, bits);

        if (FLINT_BIT_COUNT(man) != bits || (man != yman) || (exp != yexp))
        {
            printf("FAIL\n");
            printf("bits = %ld, count = %u\n\n", bits, FLINT_BIT_COUNT(man));
            printf("x = "); fmpz_print(x); printf("\n\n");
            printf("bits(x) = %ld\n\n", fmpz_bits(x));
            printf("man = %lu, exp = %ld\n", man, exp);
            printf("yman = %lu, yexp = %ld\n", yman, yexp);
            abort();
        }

        fmpz_clear(x);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
