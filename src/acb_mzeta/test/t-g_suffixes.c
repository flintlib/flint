/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version. See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_mzeta/impl.h"

TEST_FUNCTION_START(acb_mzeta_G_suffixes, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong i, L, prec;
        int series_ok, high_ok;
        acb_t z;
        acb_ptr a, series, high;

        L = 1 + n_randint(state, 8);
        prec = 40 + n_randint(state, 200);
        a = _acb_vec_init(L);
        series = _acb_vec_init(L + 1);
        high = _acb_vec_init(L + 1);
        acb_init(z);
        acb_one(z);
        acb_mul_2exp_si(z, z, -1);
        if (n_randint(state, 2))
        {
            arb_one(acb_imagref(z));
            arb_mul_2exp_si(acb_imagref(z), acb_imagref(z), -2);
        }

        for (i = 0; i < L; i++)
        {
            if (i + 1 < L && n_randint(state, 4) == 0)
                acb_zero(a + i);
            else
            {
                arb_set_si(acb_realref(a + i), 2 + n_randint(state, 6));
                if (n_randint(state, 2))
                    arb_set_si(acb_imagref(a + i),
                        (slong) n_randint(state, 5) - 2);
                if (n_randint(state, 3) == 0)
                    arb_add_error_2exp_si(acb_realref(a + i), -prec - 40);
            }
        }

        series_ok = _acb_mzeta_G_suffixes(series, a, L, z, prec);
        high_ok = _acb_mzeta_G_suffixes(high, a, L, z, prec + 64);
        if (!series_ok || !high_ok)
        {
            flint_printf("FAIL (algorithm declined, iter %wd)\n", iter);
            flint_abort();
        }
        for (i = 0; i <= L; i++)
        {
            if (!acb_is_finite(series + i) ||
                !acb_overlaps(series + i, high + i))
            {
                flint_printf("FAIL (suffix %wd, iter %wd, prec %wd)\n",
                    i, iter, prec);
                acb_printd(series + i, 20);
                flint_printf("\n");
                acb_printd(high + i, 20);
                flint_printf("\n");
                flint_abort();
            }
        }

        acb_clear(z);
        _acb_vec_clear(a, L);
        _acb_vec_clear(series, L + 1);
        _acb_vec_clear(high, L + 1);
    }

    /* Check empty words and inputs outside the series domain. */
    {
        acb_t a, z;
        acb_ptr out = _acb_vec_init(2);
        acb_init(a);
        acb_init(z);
        acb_one(z);
        if (!_acb_mzeta_G_suffixes(out, a, 0, z, 80) || !acb_is_one(out))
        {
            flint_printf("FAIL (empty word)\n");
            flint_abort();
        }
        acb_zero(a);
        if (_acb_mzeta_G_suffixes(out, a, 1, z, 80))
        {
            flint_printf("FAIL (trailing zero)\n");
            flint_abort();
        }
        acb_one(a);
        if (_acb_mzeta_G_suffixes(out, a, 1, z, 80))
        {
            flint_printf("FAIL (outside series region)\n");
            flint_abort();
        }
        acb_clear(a);
        acb_clear(z);
        _acb_vec_clear(out, 2);
    }

    /* A rotated ratio near the convergence limit exercises disk error
       propagation. A very small argument also exercises depth > N. */
    for (iter = 0; iter < 2; iter++)
    {
        slong i, L = iter ? 12 : 5;
        acb_t z;
        acb_ptr a = _acb_vec_init(L);
        acb_ptr series = _acb_vec_init(L + 1);
        acb_ptr high = _acb_vec_init(L + 1);
        acb_init(z);
        for (i = 0; i < L; i++)
            acb_one(a + i);
        if (iter == 0)
        {
            acb_zero(a);
            acb_zero(a + 3);
            arb_set_si(acb_realref(z), 5);
            arb_set_si(acb_imagref(z), 5);
            acb_mul_2exp_si(z, z, -3);
        }
        else
        {
            acb_one(z);
            acb_mul_2exp_si(z, z, -100);
        }
        if (!_acb_mzeta_G_suffixes(series, a, L, z, 80) ||
            !_acb_mzeta_G_suffixes(high, a, L, z, 144))
        {
            flint_printf("FAIL (boundary or short truncation)\n");
            flint_abort();
        }
        for (i = 0; i <= L; i++)
            if (!acb_overlaps(series + i, high + i))
            {
                flint_printf("FAIL (boundary suffix %wd)\n", i);
                flint_abort();
            }
        acb_clear(z);
        _acb_vec_clear(a, L);
        _acb_vec_clear(series, L + 1);
        _acb_vec_clear(high, L + 1);
    }

    /* Check the direct result against the first suffix for real and
       complex words with both unit and nonunit indices. */
    {
        slong i;
        acb_t z, full;
        acb_ptr a = _acb_vec_init(6);
        acb_ptr high = _acb_vec_init(7);

        acb_init(z);
        acb_init(full);
        acb_one(z);
        acb_set_si(a + 1, 3);
        acb_set_si(a + 2, 4);
        acb_set_si(a + 5, 5);
        for (i = 0; i < 2; i++)
        {
            if (i != 0)
                arb_one(acb_imagref(a + 5));
            if (!_acb_mzeta_G_one(full, a, 6, 144) ||
                !_acb_mzeta_G_suffixes(high, a, 6, z, 256) ||
                acb_rel_accuracy_bits(full) < 128 || !acb_overlaps(full, high))
            {
                flint_printf("FAIL (direct full word, iter %wd)\n", i);
                flint_abort();
            }
        }
        acb_clear(z);
        acb_clear(full);
        _acb_vec_clear(a, 6);
        _acb_vec_clear(high, 7);
    }

    /* Check the tail bound against powers of a logarithm, including
       cutoffs shorter than the depth and ratios near the series limit. */
    {
        slong i, j, L = 12;
        acb_t z, log_value, expected;
        acb_ptr a = _acb_vec_init(L);
        acb_ptr series = _acb_vec_init(L + 1);

        acb_init(z);
        acb_init(log_value);
        acb_init(expected);
        for (i = 0; i < L; i++)
            acb_one(a + i);

        for (j = 0; j < 3; j++)
        {
            if (j == 0)
                acb_zero(z);
            else if (j == 1)
            {
                acb_one(z);
                acb_mul_2exp_si(z, z, -5000);
            }
            else
            {
                acb_set_ui(z, 243);
                acb_mul_2exp_si(z, z, -8);
            }
            acb_neg(log_value, z);
            acb_log1p(log_value, log_value, 256);
            acb_one(expected);

            if (!_acb_mzeta_G_suffixes(series, a, L, z, 128))
            {
                flint_printf("FAIL (repeated letters, iter %wd)\n", j);
                flint_abort();
            }

            for (i = 1; i <= L; i++)
            {
                acb_mul(expected, expected, log_value, 256);
                acb_div_ui(expected, expected, i, 256);
                if (!acb_overlaps(series + L - i, expected) ||
                    (j == 0 && !acb_is_zero(series + L - i)) ||
                    acb_rel_one_accuracy_bits(series + L - i) < 128)
                {
                    flint_printf("FAIL (repeated letters, case %wd, d %wd)\n",
                        j, i);
                    flint_abort();
                }
            }
        }

        acb_clear(z);
        acb_clear(log_value);
        acb_clear(expected);
        _acb_vec_clear(a, L);
        _acb_vec_clear(series, L + 1);
    }

    /* Compare all suffixes at high precision with nonunit indices.
       Cover both short exact ratios and general complex ratios. */
    {
        const slong L = 8;
        acb_t z, full;
        acb_ptr a = _acb_vec_init(L);
        acb_ptr series = _acb_vec_init(L + 1);
        acb_ptr high = _acb_vec_init(L + 1);

        acb_init(z);
        acb_init(full);
        acb_one(z);
        acb_set_si(a + 5, 2);
        acb_set_si(a + 6, 4);
        acb_set_si(a + 7, 8);
        for (iter = 0; iter < 3; iter++)
        {
            if (iter == 1)
            {
                acb_zero(a + 7);
                arb_set_si(acb_imagref(a + 7), 8);
            }
            else if (iter == 2)
            {
                acb_set_si(a + 6, 3);
                arb_one(acb_imagref(a + 6));
                acb_set_si(a + 7, 4);
            }
            if (!_acb_mzeta_G_suffixes(series, a, L, z, 2048) ||
                !_acb_mzeta_G_suffixes(high, a, L, z, 2112) ||
                !_acb_vec_overlaps(series, high, L + 1) ||
                !_acb_mzeta_G_one(full, a, L, 2048) ||
                !acb_overlaps(full, high))
            {
                flint_printf("FAIL (scaled series, iter %wd)\n", iter);
                flint_abort();
            }
        }
        acb_clear(z);
        acb_clear(full);
        _acb_vec_clear(a, L);
        _acb_vec_clear(series, L + 1);
        _acb_vec_clear(high, L + 1);
    }

    /* Depth-one ratios, including every suffix of a zero block. */
    for (iter = 0; iter < 6; iter++)
    {
        slong i, L = iter == 0 ? 1 : 4;
        acb_t z, ratio, ref;
        acb_ptr a = _acb_vec_init(L);
        acb_ptr series = _acb_vec_init(L + 1);
        acb_init(z);
        acb_init(ratio);
        acb_init(ref);
        if (iter == 4)
        {
            /* (2^(FLINT_BITS-1) + 1 + i) / 2^FLINT_BITS. */
            acb_one(z);
            acb_mul_2exp_si(z, z, FLINT_BITS - 1);
            arb_add_ui(acb_realref(z), acb_realref(z), 1, FLINT_BITS + 1);
            arb_one(acb_imagref(z));
            acb_one(a + L - 1);
            acb_mul_2exp_si(a + L - 1, a + L - 1, FLINT_BITS);
            acb_mul_2exp_si(ratio, z, -FLINT_BITS);
        }
        else if (iter == 5)
        {
            /* A non-dyadic ratio near the convergence cutoff. */
            acb_set_si(z, 2);
            arb_set_ui(acb_imagref(z), 2);
            acb_set_si(a + L - 1, 3);
            acb_div(ratio, z, a + L - 1, 1120);
        }
        else
        {
            acb_set_si(z, iter == 0 ? 1 : iter == 1 ? 7 : iter == 2 ? 0 : 6);
            if (iter >= 2)
                arb_set_si(acb_imagref(z), iter == 2 ? -4 : 1);
            acb_set_si(a + L - 1, iter == 0 ? 2 : 8);
            acb_mul_2exp_si(ratio, z, iter == 0 ? -1 : -3);
        }
        if (!_acb_mzeta_G_suffixes(series, a, L, z, 1024))
        {
            flint_printf("FAIL (depth-one dispatch)\n");
            flint_abort();
        }
        for (i = 0; i <= L; i++)
        {
            if (i == L)
                acb_one(ref);
            else
            {
                acb_polylog_si(ref, L - i, ratio, 1088);
                acb_neg(ref, ref);
            }
            if (!acb_overlaps(series + i, ref) ||
                acb_rel_accuracy_bits(series + i) < 1024)
            {
                flint_printf("FAIL (depth-one case %wd, suffix %wd)\n",
                    iter, i);
                flint_abort();
            }
        }
        acb_clear(ref);
        acb_clear(ratio);
        acb_clear(z);
        _acb_vec_clear(a, L);
        _acb_vec_clear(series, L + 1);
    }

    /* Exercise both recurrences at weight 65. */
    {
        const slong L = 65, prec = 512;
        slong j;
        acb_t z, ratio, ref, log;
        acb_ptr a = _acb_vec_init(L);
        acb_ptr res = _acb_vec_init(L + 1);
        acb_init(z);
        acb_init(ratio);
        acb_init(ref);
        acb_init(log);
        acb_one(z);
        acb_set_ui(a + L - 1, 2);
        acb_mul_2exp_si(ratio, z, -1);

        if (!_acb_mzeta_G_suffixes(res, a, L, z, prec) ||
            !acb_is_one(res + L))
            flint_abort();
        for (j = 0; j < L; j++)
        {
            acb_polylog_si(ref, L - j, ratio, 192);
            acb_neg(ref, ref);
            if (!acb_overlaps(res + j, ref) ||
                acb_rel_accuracy_bits(res + j) < 128)
                flint_abort();
        }

        for (j = 0; j < L; j++)
            acb_set_ui(a + j, 2);
        arb_const_log2(acb_realref(log), prec);
        acb_neg(log, log);
        acb_one(ref);
        if (!_acb_mzeta_G_suffixes(res, a, L, z, prec) ||
            !acb_is_one(res + L))
            flint_abort();
        for (j = 1; j <= L; j++)
        {
            acb_mul(ref, ref, log, prec);
            acb_div_ui(ref, ref, j, prec);
            if (!acb_overlaps(res + L - j, ref) ||
                acb_rel_accuracy_bits(res + L - j) < 128)
                flint_abort();
        }

        acb_clear(z);
        acb_clear(ratio);
        acb_clear(ref);
        acb_clear(log);
        _acb_vec_clear(a, L);
        _acb_vec_clear(res, L + 1);
    }

    /* At this precision the repeated-letter series use over 150000 terms. */
    {
        const slong prec = 14500;
        acb_t z, ref, log;
        acb_ptr a = _acb_vec_init(2);
        acb_ptr res = _acb_vec_init(3);
        acb_init(z);
        acb_init(ref);
        acb_init(log);
        acb_set_ui(a, 16);
        acb_set(a + 1, a);
        acb_set_ui(z, 15);
        arb_const_log2(acb_realref(log), prec + 64);
        acb_mul_2exp_si(log, log, 2);
        acb_neg(log, log);
        acb_set(ref, log);

        for (iter = 1; iter <= 2; iter++)
        {
            if (!_acb_mzeta_G_suffixes(res, a, iter, z, prec) ||
                !acb_overlaps(res, ref) || acb_rel_accuracy_bits(res) < prec ||
                !acb_overlaps(res + iter - 1, log) || !acb_is_one(res + iter))
            {
                flint_printf("FAIL (long series at depth %wd)\n", iter);
                flint_abort();
            }
            acb_mul(ref, log, log, prec + 64);
            acb_mul_2exp_si(ref, ref, -1);
        }

        acb_clear(z);
        acb_clear(ref);
        acb_clear(log);
        _acb_vec_clear(a, 2);
        _acb_vec_clear(res, 3);
    }

    /* Blockwise disk propagation controls wrapping for a complex ratio. */
    {
        const slong prec = 256;
        acb_t a, z, endpoint, ref;
        acb_ptr series = _acb_vec_init(2);
        acb_init(a);
        acb_init(z);
        acb_init(endpoint);
        acb_init(ref);
        acb_set_ui(a, 16);
        arb_add_error_2exp_si(acb_realref(a), -240);
        acb_set_ui(z, 12);
        arb_set_ui(acb_imagref(z), 9);
        if (!_acb_mzeta_G_suffixes(series, a, 1, z, prec) ||
            acb_rel_accuracy_bits(series) < 240)
        {
            flint_printf("FAIL (depth-one complex ball)\n");
            flint_abort();
        }
        for (iter = -1; iter <= 1; iter += 2)
        {
            acb_set_si(endpoint, iter);
            acb_mul_2exp_si(endpoint, endpoint, -240);
            acb_add_ui(endpoint, endpoint, 16, prec + 128);
            acb_div(ref, z, endpoint, prec + 128);
            acb_neg(ref, ref);
            acb_log1p(ref, ref, prec + 128);
            if (!acb_contains(series, ref))
            {
                flint_printf("FAIL (depth-one complex ball endpoint)\n");
                flint_abort();
            }
        }
        acb_clear(a);
        acb_clear(z);
        acb_clear(endpoint);
        acb_clear(ref);
        _acb_vec_clear(series, 2);
    }

    /* Return the series enclosure when the input limits the accuracy. */
    {
        acb_t a, z, endpoint, ref;
        acb_ptr x;
        acb_init(a);
        acb_init(z);
        acb_init(endpoint);
        acb_init(ref);
        acb_set_si(a, 2);
        arb_add_error_2exp_si(acb_realref(a), -10);
        acb_one(z);
        acb_mul_2exp_si(z, z, -1);
        x = _acb_vec_init(2);
        if (!_acb_mzeta_G_suffixes(x, a, 1, z, 128))
        {
            flint_printf("FAIL (wide input series)\n");
            flint_abort();
        }
        for (iter = -1; iter <= 1; iter += 2)
        {
            acb_set_si(endpoint, iter);
            acb_mul_2exp_si(endpoint, endpoint, -10);
            acb_add_ui(endpoint, endpoint, 2, 192);
            acb_div(ref, z, endpoint, 192);
            acb_neg(ref, ref);
            acb_log1p(ref, ref, 192);
            if (!acb_contains(x, ref))
            {
                flint_printf("FAIL (wide input endpoint)\n");
                flint_abort();
            }
        }
        _acb_vec_clear(x, 2);
        acb_clear(a);
        acb_clear(z);
        acb_clear(endpoint);
        acb_clear(ref);
    }

    /* The outer ratio is less accurate than the resulting suffix sums. */
    {
        acb_t z;
        acb_ptr a = _acb_vec_init(2);
        acb_ptr series = _acb_vec_init(3);
        acb_ptr high = _acb_vec_init(3);
        acb_init(z);
        acb_set_ui(a, 16);
        acb_set_ui(a + 1, 16);
        arb_add_error_2exp_si(acb_realref(a), -99);
        acb_one(z);
        acb_mul_2exp_si(z, z, -5);
        if (!_acb_mzeta_G_suffixes(series, a, 2, z, 128) ||
            !_acb_mzeta_G_suffixes(high, a, 2, z, 192))
        {
            flint_printf("FAIL (inexact outer ratio)\n");
            flint_abort();
        }
        for (iter = 0; iter <= 2; iter++)
            if (!acb_overlaps(series + iter, high + iter))
            {
                flint_printf("FAIL (inexact outer ratio, suffix %wd)\n", iter);
                flint_abort();
            }
        acb_clear(z);
        _acb_vec_clear(a, 2);
        _acb_vec_clear(series, 3);
        _acb_vec_clear(high, 3);
    }

    TEST_FUNCTION_END(state);
}
