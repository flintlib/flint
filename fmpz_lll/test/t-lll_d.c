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

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_lll.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result = 1;
    fmpz_mat_t mat;
    mp_bitcnt_t bits;
    FLINT_TEST_INIT(state);

    flint_printf("lll_d....");
    fflush(stdout);

    /* test using NTRU like matrices */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong q;
        slong r, c;
        double randd, delta, eta;
        d_mat_t Q, R;

        r = 2 * (n_randint(state, 50) + 1);
        c = r;

        fmpz_mat_init(mat, r, c);

        bits = n_randint(state, 20) + 1;
        q = n_randint(state, 200) + 1;

        if (n_randint(state, 1))
            fmpz_mat_randntrulike(mat, state, bits, q);
        else
            fmpz_mat_randntrulike2(mat, state, bits, q);

        randd = d_randtest(state);
        if (randd > 0.5 && n_randint(state, 1))
        {
            delta = 0.25 + (randd - 0.5) * 0.75;
            if (n_randint(state, 1))
                eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
            else
                eta = 0.5 + randd * (sqrt(delta) - 0.5);
        }
        else
        {
            delta = 0.25 + randd * 0.75;
            if (n_randint(state, 1))
                eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
            else
                eta = 0.5 + randd * (sqrt(delta) - 0.5);
        }

        fmpz_lll_d(mat, delta, eta);

        d_mat_init(Q, r, c);
        d_mat_init(R, r, r);
        d_mat_zero(R);

        fmpz_mat_qr_d(Q, R, mat);

        result = d_mat_reduced(R, delta, eta);

        d_mat_clear(Q);
        d_mat_clear(R);

        if (!result)
        {
            flint_printf("FAIL (randntrulike):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("bits = %ld, i = %ld\n", bits, i);
            flint_printf("delta = %g, eta = %g\n", delta, eta);
            abort();
        }

        fmpz_mat_clear(mat);
    }

    /* test using integer relations matrices */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong r, c;
        double randd, delta, eta;
        d_mat_t Q, R;

        r = n_randint(state, 20) + 1;
        c = r + 1;

        fmpz_mat_init(mat, r, c);

        bits = n_randint(state, 200) + 1;

        fmpz_mat_randintrel(mat, state, bits);

        randd = d_randtest(state);
        if (randd > 0.5 && n_randint(state, 1))
        {
            delta = 0.25 + (randd - 0.5) * 0.75;
            if (n_randint(state, 1))
                eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
            else
                eta = 0.5 + randd * (sqrt(delta) - 0.5);
        }
        else
        {
            delta = 0.25 + randd * 0.75;
            if (n_randint(state, 1))
                eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
            else
                eta = 0.5 + randd * (sqrt(delta) - 0.5);
        }

        result = fmpz_lll_d(mat, delta, eta);
        if (result == -1)
        {
            fmpz_mat_clear(mat);
            continue;
        }

        d_mat_init(Q, r, c);
        d_mat_init(R, r, r);
        d_mat_zero(R);

        fmpz_mat_qr_d(Q, R, mat);

        result = d_mat_reduced(R, delta, eta);

        d_mat_clear(Q);
        d_mat_clear(R);

        if (!result)
        {
            flint_printf("FAIL (randintrel):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("bits = %ld, i = %ld\n", bits, i);
            flint_printf("delta = %g, eta = %g\n", delta, eta);
            abort();
        }

        fmpz_mat_clear(mat);
    }

    /* test using ajtai matrices */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong r, c;
        double randd, delta, eta;
        d_mat_t Q, R;

        r = n_randint(state, 10) + 1;
        c = r;

        fmpz_mat_init(mat, r, c);

        bits = n_randint(state, 200) + 1;

        fmpz_mat_randajtai(mat, state, 0.5);

        randd = d_randtest(state);
        if (randd > 0.5 && n_randint(state, 1))
        {
            delta = 0.25 + (randd - 0.5) * 0.75;
            if (n_randint(state, 1))
                eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
            else
                eta = 0.5 + randd * (sqrt(delta) - 0.5);
        }
        else
        {
            delta = 0.25 + randd * 0.75;
            if (n_randint(state, 1))
                eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
            else
                eta = 0.5 + randd * (sqrt(delta) - 0.5);
        }

        fmpz_lll_d(mat, delta, eta);

        d_mat_init(Q, r, c);
        d_mat_init(R, r, r);
        d_mat_zero(R);

        fmpz_mat_qr_d(Q, R, mat);

        result = d_mat_reduced(R, delta, eta);

        d_mat_clear(Q);
        d_mat_clear(R);

        if (!result)
        {
            flint_printf("FAIL (randajtai):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("bits = %ld, i = %ld\n", bits, i);
            flint_printf("delta = %g, eta = %g\n", delta, eta);
            abort();
        }

        fmpz_mat_clear(mat);
    }

    /* test using simultaneous diophantine matrices */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mp_bitcnt_t bits2;
        slong r, c;
        double randd, delta, eta;
        d_mat_t Q, R;

        r = n_randint(state, 200) + 1;
        c = r;

        fmpz_mat_init(mat, r, c);

        bits = n_randint(state, 200) + 1;
        bits2 = n_randint(state, 5) + 1;

        fmpz_mat_randsimdioph(mat, state, bits, bits2);

        randd = d_randtest(state);
        if (randd > 0.5 && n_randint(state, 1))
        {
            delta = 0.25 + (randd - 0.5) * 0.75;
            if (n_randint(state, 1))
                eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
            else
                eta = 0.5 + randd * (sqrt(delta) - 0.5);
        }
        else
        {
            delta = 0.25 + randd * 0.75;
            if (n_randint(state, 1))
                eta = 0.5 + (randd - 0.5) * (sqrt(delta) - 0.5);
            else
                eta = 0.5 + randd * (sqrt(delta) - 0.5);
        }

        fmpz_lll_d(mat, delta, eta);

        d_mat_init(Q, r, c);
        d_mat_init(R, r, r);
        d_mat_zero(R);

        fmpz_mat_qr_d(Q, R, mat);

        result = d_mat_reduced(R, delta, eta);

        d_mat_clear(Q);
        d_mat_clear(R);

        if (!result)
        {
            flint_printf("FAIL (randsimdioph):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("bits = %ld, i = %ld\n", bits, i);
            flint_printf("delta = %g, eta = %g\n", delta, eta);
            abort();
        }

        fmpz_mat_clear(mat);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
