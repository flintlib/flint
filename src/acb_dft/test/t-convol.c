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
#define do_f hogusbogusname
typedef void (*do_f) (acb_ptr z, acb_srcptr x, acb_srcptr y, slong len, slong prec);

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

TEST_FUNCTION_START(acb_dft_convol, state)
{
    slong k;
    slong prec = 100, digits = 30;
    slong nq = 17;
    ulong q[17] = { 0, 1, 2, 3, 4, 5, 6, 23, 10, 15, 30, 59, 256, 308, 335, 344, 961};
    ulong nr = 5;

    slong f, nf = 4;
    do_f func[4] = { acb_dft_convol_naive, acb_dft_convol_rad2, acb_dft_convol_dft, acb_dft_convol_mullow };
    char * name[4] = { "naive", "rad2", "dft", "mullow" };

    for (k = 0; k < nq + nr; k++)
    {
        slong i, len;
        acb_ptr z1, z2, x, y;

        if (k < nq)
            len = q[k];
        else
            len = n_randint(state, 2000);

        z1 = _acb_vec_init(len);
        z2 = _acb_vec_init(len);
        x = _acb_vec_init(len);
        y = _acb_vec_init(len);

        for (i = 0; i < len; i++)
        {
#if 1
            acb_set_si(x + i, n_randint(state, 4 * len));
            acb_set_si(y + i, n_randint(state, 4 * len));
#else
            acb_set_si_si(x + i, n_randint(state, 4 * len), n_randint(state, 4 * len));
            acb_set_si_si(y + i, n_randint(state, 4 * len), n_randint(state, 4 * len));
#endif
        }

        for (f = 0; f < nf; f++)
        {

            acb_ptr z = (f == 0) ? z1 : z2;

            func[f](z, x, y, len, prec);

            if (f == 0)
                continue;

            check_vec_eq_prec(z1, z2, len, prec, digits, len, "", name[0], name[f]);

        }

        _acb_vec_clear(x, len);
        _acb_vec_clear(y, len);
        _acb_vec_clear(z1, len);
        _acb_vec_clear(z2, len);
    }

    TEST_FUNCTION_END(state);
}
#undef do_f
