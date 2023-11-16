/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpzi.h"

void
fmpzi_divrem_ref(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)
{
    fmpzi_t t, y_conj;
    fmpz_t v;
    mpz_t ytmp;

    fmpzi_init(t);
    fmpz_init(v);

    /* shallow conjugate */
    *fmpzi_realref(y_conj) = *fmpzi_realref(y);
    if (!COEFF_IS_MPZ(*fmpzi_imagref(y)))
    {
        *fmpzi_imagref(y_conj) = -*fmpzi_imagref(y);
    }
    else
    {
        *ytmp = *COEFF_TO_PTR(*fmpzi_imagref(y));
        mpz_neg(ytmp, ytmp);
        *fmpzi_imagref(y_conj) = PTR_TO_COEFF(ytmp);
    }

    fmpzi_mul(t, x, y_conj);

    fmpz_mul_2exp(fmpzi_realref(t), fmpzi_realref(t), 1);
    fmpz_mul_2exp(fmpzi_imagref(t), fmpzi_imagref(t), 1);

    fmpz_fmma(v, fmpzi_realref(y), fmpzi_realref(y),
                 fmpzi_imagref(y), fmpzi_imagref(y));

    fmpz_add(fmpzi_realref(t), fmpzi_realref(t), v);
    fmpz_add(fmpzi_imagref(t), fmpzi_imagref(t), v);

    fmpz_mul_2exp(v, v, 1);

    fmpz_fdiv_q(fmpzi_realref(q), fmpzi_realref(t), v);
    fmpz_fdiv_q(fmpzi_imagref(q), fmpzi_imagref(t), v);

    if (r != NULL)
    {
        fmpzi_mul(t, q, y);
        fmpzi_sub(r, x, t);
    }

    fmpzi_clear(t);
    fmpz_clear(v);
}

TEST_FUNCTION_START(fmpzi_divrem, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpzi_t x, y, q, r, s, t;
        fmpz_t ny, nr;
        int aliasing;

        fmpzi_init(x);
        fmpzi_init(y);
        fmpzi_init(q);
        fmpzi_init(r);
        fmpzi_init(s);
        fmpzi_init(t);
        fmpz_init(ny);
        fmpz_init(nr);

        fmpzi_randtest(x, state, 2 + n_randint(state, 200));
        do {
            fmpzi_randtest(y, state, 2 + n_randint(state, 200));
        } while (fmpzi_is_zero(y));

        fmpzi_randtest(q, state, 2 + n_randint(state, 200));
        fmpzi_randtest(r, state, 2 + n_randint(state, 200));

        if (n_randint(state, 3) == 0)
        {
            fmpzi_mul(x, y, q);
            if (n_randint(state, 3) == 0)
                fmpzi_add(x, x, r);
        }

        aliasing = n_randint(state, 7);
        switch (aliasing)
        {
            case 0:
                fmpzi_divrem(q, r, x, y);
                break;
            case 1:
                fmpzi_set(q, x);
                fmpzi_divrem(q, r, q, y);
                break;
            case 2:
                fmpzi_set(q, y);
                fmpzi_divrem(q, r, x, q);
                break;
            case 3:
                fmpzi_set(r, x);
                fmpzi_divrem(q, r, r, y);
                break;
            case 4:
                fmpzi_set(r, y);
                fmpzi_divrem(q, r, x, r);
                break;
            case 5:
                fmpzi_set(q, x);
                fmpzi_set(r, y);
                fmpzi_divrem(q, r, q, r);
                break;
            default:
                fmpzi_set(q, y);
                fmpzi_set(r, x);
                fmpzi_divrem(q, r, r, q);
                break;
        }

        fmpzi_mul(s, q, y);
        fmpzi_add(t, s, r);

        fmpzi_norm(ny, y);
        fmpzi_norm(nr, r);
        fmpz_mul_2exp(nr, nr, 1);

        if (!fmpzi_equal(t, x) || !(fmpz_cmp(nr, ny) <= 0))
        {
            flint_printf("FAIL\n");
            flint_printf("aliasing = %d\n", aliasing);
            flint_printf("x = "); fmpzi_print(x); printf("\n");
            flint_printf("y = "); fmpzi_print(y); printf("\n");
            flint_printf("q = "); fmpzi_print(q); printf("\n");
            flint_printf("r = "); fmpzi_print(r); printf("\n");
            flint_printf("s = "); fmpzi_print(s); printf("\n");
            flint_printf("t = "); fmpzi_print(t); printf("\n");
            flint_printf("nr = "); fmpz_print(nr); printf("\n");
            flint_printf("ny = "); fmpz_print(ny); printf("\n");
            flint_abort();
        }

        fmpzi_divrem_ref(s, t, x, y);

        if (!fmpzi_equal(q, s) || !fmpzi_equal(r, t))
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); fmpzi_print(x); printf("\n");
            flint_printf("y = "); fmpzi_print(y); printf("\n");
            flint_printf("q = "); fmpzi_print(q); printf("\n");
            flint_printf("r = "); fmpzi_print(r); printf("\n");
            flint_printf("s = "); fmpzi_print(s); printf("\n");
            flint_printf("t = "); fmpzi_print(t); printf("\n");
            flint_abort();
        }

        fmpzi_clear(x);
        fmpzi_clear(y);
        fmpzi_clear(q);
        fmpzi_clear(r);
        fmpzi_clear(s);
        fmpzi_clear(t);
        fmpz_clear(ny);
        fmpz_clear(nr);
    }

    TEST_FUNCTION_END(state);
}
