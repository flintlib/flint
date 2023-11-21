/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "acb_mat.h"

static void
_acb_mat_set_fmpq_mat_fmpq_mat(acb_mat_t C,
        const fmpq_mat_t rsrc, const fmpq_mat_t isrc, slong prec)
{
    slong i, j;
    for (i = 0; i < acb_mat_nrows(C); i++)
    {
        for (j = 0; j < acb_mat_ncols(C); j++)
        {
            arb_set_fmpq(acb_realref(acb_mat_entry(C, i, j)),
                    fmpq_mat_entry(rsrc, i, j), prec);
            arb_set_fmpq(acb_imagref(acb_mat_entry(C, i, j)),
                    fmpq_mat_entry(isrc, i, j), prec);
        }
    }
}

static void
_acb_mat_conjugate_transpose(acb_mat_t B, const acb_mat_t A)
{
    slong i, j;
    acb_mat_transpose(B, A);
    for (i = 0; i < acb_mat_nrows(B); i++)
        for (j = 0; j < acb_mat_ncols(B); j++)
            acb_conj(acb_mat_entry(B, i, j), acb_mat_entry(B, i, j));
}

static void
_fmpq_mat_sum_of_squares(fmpq_t res, const fmpq_mat_t Q)
{
    slong i, j;
    fmpq_zero(res);
    for (i = 0; i < fmpq_mat_nrows(Q); i++)
    {
        for (j = 0; j < fmpq_mat_ncols(Q); j++)
        {
            fmpq_addmul(res, fmpq_mat_entry(Q, i, j), fmpq_mat_entry(Q, i, j));
        }
    }
}

TEST_FUNCTION_START(acb_mat_frobenius_norm, state)
{
    slong iter;

    /* compare to the exact rational norm */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_mat_t Qr, Qi;
        fmpq_t q;
        acb_mat_t A;
        slong n, qbits, prec;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Qr, n, n);
        fmpq_mat_init(Qi, n, n);
        fmpq_init(q);

        acb_mat_init(A, n, n);

        fmpq_mat_randtest(Qr, state, qbits);
        fmpq_mat_randtest(Qi, state, qbits);

        /* compute the square of the exact rational norm */
        {
            fmpq_t qr, qi;
            fmpq_init(qr);
            fmpq_init(qi);
            _fmpq_mat_sum_of_squares(qr, Qr);
            _fmpq_mat_sum_of_squares(qi, Qi);
            fmpq_add(q, qr, qi);
            fmpq_clear(qr);
            fmpq_clear(qi);
        }

        _acb_mat_set_fmpq_mat_fmpq_mat(A, Qr, Qi, prec);

        /* check that the arb interval contains the exact value */
        {
            arb_t a;
            arb_init(a);

            acb_mat_frobenius_norm(a, A, prec);
            arb_mul(a, a, a, prec);

            if (!arb_contains_fmpq(a, q))
            {
                flint_printf("FAIL (containment, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Qr = \n"); fmpq_mat_print(Qr); flint_printf("\n\n");
                flint_printf("Qi = \n"); fmpq_mat_print(Qi); flint_printf("\n\n");
                flint_printf("frobenius_norm(Q)^2 = \n");
                fmpq_print(q); flint_printf("\n\n");

                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("frobenius_norm(A)^2 = \n");
                arb_printd(a, 15); flint_printf("\n\n");
                flint_printf("frobenius_norm(A)^2 = \n");
                arb_print(a); flint_printf("\n\n");

                flint_abort();
            }

            arb_clear(a);
        }

        /* check that the upper bound is not less than the exact value */
        {
            mag_t b;
            fmpq_t y;

            mag_init(b);
            fmpq_init(y);

            acb_mat_bound_frobenius_norm(b, A);
            mag_mul(b, b, b);
            mag_get_fmpq(y, b);

            if (fmpq_cmp(q, y) > 0)
            {
                flint_printf("FAIL (bound, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Qr = \n"); fmpq_mat_print(Qr); flint_printf("\n\n");
                flint_printf("Qi = \n"); fmpq_mat_print(Qi); flint_printf("\n\n");
                flint_printf("frobenius_norm(Q)^2 = \n");
                fmpq_print(q); flint_printf("\n\n");

                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("bound_frobenius_norm(A)^2 = \n");
                mag_printd(b, 15); flint_printf("\n\n");
                flint_printf("bound_frobenius_norm(A)^2 = \n");
                mag_print(b); flint_printf("\n\n");

                flint_abort();
            }

            mag_clear(b);
            fmpq_clear(y);
        }

        fmpq_mat_clear(Qr);
        fmpq_mat_clear(Qi);
        fmpq_clear(q);
        acb_mat_clear(A);
    }

    /* check trace(A^H A) = frobenius_norm(A)^2 */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, prec;
        acb_mat_t A, AH, AHA;
        acb_t t;

        prec = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        acb_mat_init(A, m, n);
        acb_mat_init(AH, n, m);
        acb_mat_init(AHA, n, n);
        acb_init(t);

        acb_mat_randtest(A, state, 2 + n_randint(state, 100), 10);
        _acb_mat_conjugate_transpose(AH, A);
        acb_mat_mul(AHA, AH, A, prec);
        acb_mat_trace(t, AHA, prec);
        acb_sqrt(t, t, prec);

        /* check the norm bound */
        {
            mag_t low, frobenius;

            mag_init(low);
            acb_get_mag_lower(low, t);

            mag_init(frobenius);
            acb_mat_bound_frobenius_norm(frobenius, A);

            if (mag_cmp(low, frobenius) > 0)
            {
                flint_printf("FAIL (frobenius norm bound)\n", iter);
                flint_printf("m = %wd, n = %wd, prec = %wd\n", m, n, prec);
                flint_printf("\n");

                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");

                flint_printf("lower(sqrt(trace(A^H A))) = \n");
                mag_printd(low, 15); flint_printf("\n\n");

                flint_printf("bound_frobenius_norm(A) = \n");
                mag_printd(frobenius, 15); flint_printf("\n\n");

                flint_abort();
            }

            mag_clear(low);
            mag_clear(frobenius);
        }

        /* check the norm interval */
        {
            arb_t frobenius;

            arb_init(frobenius);
            acb_mat_frobenius_norm(frobenius, A, prec);

            if (!arb_overlaps(acb_realref(t), frobenius))
            {
                flint_printf("FAIL (frobenius norm overlap)\n", iter);
                flint_printf("m = %wd, n = %wd, prec = %wd\n", m, n, prec);
                flint_printf("\n");

                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");

                flint_printf("sqrt(trace(A^H A)) = \n");
                acb_printd(t, 15); flint_printf("\n\n");

                flint_printf("frobenius_norm(A) = \n");
                arb_printd(frobenius, 15); flint_printf("\n\n");

                flint_abort();
            }

            arb_clear(frobenius);
        }

        acb_mat_clear(A);
        acb_mat_clear(AH);
        acb_mat_clear(AHA);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
