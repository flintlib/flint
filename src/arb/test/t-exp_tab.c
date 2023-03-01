/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong i;

    mpfr_t tabx, expx, y1, y2;
    mpz_t tt;

    flint_printf("exp_tab....");
    fflush(stdout);

    {
        slong prec, bits, num;

        prec = ARB_EXP_TAB1_LIMBS * FLINT_BITS;
        bits = ARB_EXP_TAB1_BITS;
        num = ARB_EXP_TAB1_NUM;

        mpfr_init2(tabx, prec);
        mpfr_init2(expx, prec);
        mpfr_init2(y1, prec);
        mpfr_init2(y2, prec);

        for (i = 0; i < num; i++)
        {
            tt->_mp_d = (mp_ptr) arb_exp_tab1[i];
            tt->_mp_size = prec / FLINT_BITS;
            tt->_mp_alloc = tt->_mp_size;

            while (tt->_mp_size > 0 && tt->_mp_d[tt->_mp_size-1] == 0)
                tt->_mp_size--;

            mpfr_set_z(tabx, tt, MPFR_RNDD);
            mpfr_div_2ui(tabx, tabx, prec, MPFR_RNDD);

            mpfr_set_ui(expx, i, MPFR_RNDD);
            mpfr_div_2ui(expx, expx, bits, MPFR_RNDD);
            mpfr_exp(expx, expx, MPFR_RNDD);

            mpfr_mul_2ui(y1, tabx, prec, MPFR_RNDD);
            mpfr_floor(y1, y1);
            mpfr_div_2ui(y1, y1, prec, MPFR_RNDD);

            mpfr_mul_2ui(y2, expx, prec - 1, MPFR_RNDD);
            mpfr_floor(y2, y2);
            mpfr_div_2ui(y2, y2, prec, MPFR_RNDD);

            if (!mpfr_equal_p(y1, y2))
            {
                flint_printf("FAIL: i = %wd, bits = %wd, prec = %wd\n", i, bits, prec);
                mpfr_printf("y1 = %.1500Rg\n", y1);
                mpfr_printf("y2 = %.1500Rg\n", y2);
                flint_abort();
            }
        }

        mpfr_clear(tabx);
        mpfr_clear(expx);
        mpfr_clear(y1);
        mpfr_clear(y2);
    }

    {
        slong prec, bits, num;

        prec = ARB_EXP_TAB2_LIMBS * FLINT_BITS;
        bits = ARB_EXP_TAB21_BITS;
        num = ARB_EXP_TAB21_NUM;

        mpfr_init2(tabx, prec);
        mpfr_init2(expx, prec);
        mpfr_init2(y1, prec);
        mpfr_init2(y2, prec);

        for (i = 0; i < num; i++)
        {
            tt->_mp_d = (mp_ptr) arb_exp_tab21[i];
            tt->_mp_size = prec / FLINT_BITS;
            tt->_mp_alloc = tt->_mp_size;

            while (tt->_mp_size > 0 && tt->_mp_d[tt->_mp_size-1] == 0)
                tt->_mp_size--;

            mpfr_set_z(tabx, tt, MPFR_RNDD);
            mpfr_div_2ui(tabx, tabx, prec, MPFR_RNDD);

            mpfr_set_ui(expx, i, MPFR_RNDD);
            mpfr_div_2ui(expx, expx, bits, MPFR_RNDD);
            mpfr_exp(expx, expx, MPFR_RNDD);

            mpfr_mul_2ui(y1, tabx, prec, MPFR_RNDD);
            mpfr_floor(y1, y1);
            mpfr_div_2ui(y1, y1, prec, MPFR_RNDD);

            mpfr_mul_2ui(y2, expx, prec - 1, MPFR_RNDD);
            mpfr_floor(y2, y2);
            mpfr_div_2ui(y2, y2, prec, MPFR_RNDD);

            if (!mpfr_equal_p(y1, y2))
            {
                flint_printf("FAIL: i = %wd, bits = %wd, prec = %wd\n", i, bits, prec);
                mpfr_printf("y1 = %.1500Rg\n", y1);
                mpfr_printf("y2 = %.1500Rg\n", y2);
                flint_abort();
            }
        }

        mpfr_clear(tabx);
        mpfr_clear(expx);
        mpfr_clear(y1);
        mpfr_clear(y2);
    }

    {
        slong prec, bits, num;

        prec = ARB_EXP_TAB2_LIMBS * FLINT_BITS;
        bits = ARB_EXP_TAB21_BITS + ARB_EXP_TAB22_BITS;
        num = ARB_EXP_TAB22_NUM;

        mpfr_init2(tabx, prec);
        mpfr_init2(expx, prec);
        mpfr_init2(y1, prec);
        mpfr_init2(y2, prec);

        for (i = 0; i < num; i++)
        {
            tt->_mp_d = (mp_ptr) arb_exp_tab22[i];
            tt->_mp_size = prec / FLINT_BITS;
            tt->_mp_alloc = tt->_mp_size;

            while (tt->_mp_size > 0 && tt->_mp_d[tt->_mp_size-1] == 0)
                tt->_mp_size--;

            mpfr_set_z(tabx, tt, MPFR_RNDD);
            mpfr_div_2ui(tabx, tabx, prec, MPFR_RNDD);

            mpfr_set_ui(expx, i, MPFR_RNDD);
            mpfr_div_2ui(expx, expx, bits, MPFR_RNDD);
            mpfr_exp(expx, expx, MPFR_RNDD);

            mpfr_mul_2ui(y1, tabx, prec, MPFR_RNDD);
            mpfr_floor(y1, y1);
            mpfr_div_2ui(y1, y1, prec, MPFR_RNDD);

            mpfr_mul_2ui(y2, expx, prec - 1, MPFR_RNDD);
            mpfr_floor(y2, y2);
            mpfr_div_2ui(y2, y2, prec, MPFR_RNDD);

            if (!mpfr_equal_p(y1, y2))
            {
                flint_printf("FAIL: i = %wd, bits = %wd, prec = %wd\n", i, bits, prec);
                mpfr_printf("y1 = %.1500Rg\n", y1);
                mpfr_printf("y2 = %.1500Rg\n", y2);
                flint_abort();
            }
        }

        mpfr_clear(tabx);
        mpfr_clear(expx);
        mpfr_clear(y1);
        mpfr_clear(y2);
    }

    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

