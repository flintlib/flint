/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dft.h"

/* Both of these are defined in t-convol.c and t-dft.c. However, do_f is not the
 * same. */
#define do_f rajtantajtan
typedef void (*do_f) (acb_ptr w, acb_srcptr v, slong len, slong prec);

#ifndef check_vec_eq_prec
#define check_vec_eq_prec check_vec_eq_prec
void
check_vec_eq_prec(acb_srcptr w1, acb_srcptr w2, slong len, slong prec, slong digits, ulong q, char d[], char f1[], char f2[])
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!acb_overlaps(w1 + i, w2 + i))
        {
            flint_printf("FAIL\n\n");
            flint_printf("q = %wu, size = %wd\n", q, len);
            flint_printf("\nDFT %sdiffer from index %wd / %wd \n", d, i, len);
            flint_printf("\n%s =\n", f1);
            acb_vec_printd_index(w1, len, digits);
            flint_printf("\n%s =\n", f2);
            acb_vec_printd_index(w2, len, digits);
            flint_printf("\n\n");
            flint_abort();
        }
        else if (!acb_is_zero(w1+i) && (acb_rel_accuracy_bits(w1 + i) < 30
                || acb_rel_accuracy_bits(w2 + i) < 30))
        {
            flint_printf("FAIL\n\n");
            flint_printf("q = %wu\n", q);
            flint_printf("\nDFT inaccurate from index %wd / %wd \n", i, len);
            flint_printf("\n%s =\n", f1);
            acb_printd(w1 + i, digits);
            flint_printf("\n%s =\n", f2);
            acb_printd(w2 + i, digits);
            flint_printf("\nerrors %wd & %wd [prec = %wd]\n",
                    acb_rel_accuracy_bits(w1 + i),
                    acb_rel_accuracy_bits(w2 + i), prec);
            flint_abort();
        }
    }
}
#endif

TEST_FUNCTION_START(acb_dft, state)
{
    slong k;
    slong prec = 100, digits = 30;
    slong nq = 19;
    ulong q[19] = { 0, 1, 2, 3, 4, 5, 6, 23, 10, 15, 16, 30, 59, 125, 308, 335, 525, 961, 1225};
    slong nr = 5;

    slong f, nf = 5;
    do_f func[5] = { acb_dft_naive, acb_dft_cyc, acb_dft_crt, acb_dft_bluestein, acb_dft };
    char * name[5] = { "naive", "cyc", "crt", "bluestein", "default" };

    /* cyclic dft */
    for (k = 0; k < nq + nr; k++)
    {
        slong i, len, f0;
        acb_ptr v, w1, w2, w3;

        if (k < nq)
            len = q[k];
        else
            len = n_randint(state, 1000);

        v = _acb_vec_init(len);
        w1 = _acb_vec_init(len);
        w2 = _acb_vec_init(len);
        w3 = _acb_vec_init(len);

        for (i = 0; i < len; i++)
            acb_set_si_si(v + i, i, 3 - i);

        /* avoid naive for long transforms */
        f0 = (len > 50);

        for (f = f0; f < nf; f++)
        {

            acb_ptr w = (f == f0) ? w1 : w2;

            if (DFT_VERB)
                flint_printf("\n%s %wu\n", name[f], len);

            /* compute DFT */
            func[f](w, v, len, prec);

            if (len < 500)
            {
                /* check aliasing */
                _acb_vec_set(w3, v, len);
                func[f](w3, w3, len, prec);

                check_vec_eq_prec(w1, w3, len, prec, digits, len, "alias ", name[0], name[f]);
            }

            if (f > f0)
            {
                /* check non aliased */
                check_vec_eq_prec(w1, w2, len, prec, digits, len, "no alias ", name[0], name[f]);
            }
            else
            {
                /* check inverse */
                acb_dft_inverse(w2, w1, len, prec);
                check_vec_eq_prec(v, w2, len, prec, digits, len, "inverse ", "original", "inverse");
            }
        }

        _acb_vec_clear(v, len);
        _acb_vec_clear(w1, len);
        _acb_vec_clear(w2, len);
        _acb_vec_clear(w3, len);
    }

    /* radix2 dft */
    for (k = 0; k < 11; k++)
    {
        slong n = 1 << k, j;
        acb_ptr v, w1, w2;
        v = w2 = _acb_vec_init(n);
        w1 = _acb_vec_init(n);

        for (j = 0; j < n; j++)
            acb_set_si_si(v + j, j, j + 2);

        acb_dft_cyc(w1, v, n, prec);
        acb_dft_rad2_inplace(w2, k, prec);

        check_vec_eq_prec(w1, w2, n, prec, digits, n, "rad2 ", "cyc", "rad2");

        _acb_vec_clear(v, n);
        _acb_vec_clear(w1, n);

    }

    /* multi-threaded radix2 dft */
    for (k = 0; k < 11; k++)
    {
        slong n = 1 << k, j;
        acb_ptr v, w1, w2;
        v = w2 = _acb_vec_init(n);
        w1 = _acb_vec_init(n);

        flint_set_num_threads(k % 5 + 1);

        for (j = 0; j < n; j++)
            acb_set_si_si(v + j, j, j + 2);

        acb_dft_cyc(w1, v, n, prec);
        acb_dft_rad2_inplace_threaded(w2, k, prec);

        check_vec_eq_prec(w1, w2, n, prec, digits, n, "rad2 ", "cyc", "rad2");

        _acb_vec_clear(v, n);
        _acb_vec_clear(w1, n);

    }

    TEST_FUNCTION_END(state);
}
#undef do_f
