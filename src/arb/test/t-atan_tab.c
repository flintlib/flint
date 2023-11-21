/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "arb.h"

TEST_FUNCTION_START(arb_atan_tab, state)
{
    slong i;

    mpfr_t tabx, atanx, y1, y2;
    mpz_t tt;

    {
        slong prec, bits, num;

        prec = ARB_ATAN_TAB1_LIMBS * FLINT_BITS;
        bits = ARB_ATAN_TAB1_BITS;
        num = 1 << ARB_ATAN_TAB1_BITS;

        mpfr_init2(tabx, prec);
        mpfr_init2(atanx, prec);
        mpfr_init2(y1, prec);
        mpfr_init2(y2, prec);

        for (i = 0; i < num; i++)
        {
            tt->_mp_d = (mp_ptr) arb_atan_tab1[i];
            tt->_mp_size = prec / FLINT_BITS;
            tt->_mp_alloc = tt->_mp_size;

            while (tt->_mp_size > 0 && tt->_mp_d[tt->_mp_size-1] == 0)
                tt->_mp_size--;

            mpfr_set_z(tabx, tt, MPFR_RNDD);
            mpfr_div_2ui(tabx, tabx, prec, MPFR_RNDD);

            mpfr_set_ui(atanx, i, MPFR_RNDD);
            mpfr_div_2ui(atanx, atanx, bits, MPFR_RNDD);
            mpfr_atan(atanx, atanx, MPFR_RNDD);

            mpfr_mul_2ui(y1, tabx, prec, MPFR_RNDD);
            mpfr_floor(y1, y1);
            mpfr_div_2ui(y1, y1, prec, MPFR_RNDD);

            mpfr_mul_2ui(y2, atanx, prec, MPFR_RNDD);
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
        mpfr_clear(atanx);
        mpfr_clear(y1);
        mpfr_clear(y2);
    }

    {
        slong prec, bits, num;

        prec = ARB_ATAN_TAB2_LIMBS * FLINT_BITS;
        bits = ARB_ATAN_TAB21_BITS;
        num = 1 << ARB_ATAN_TAB21_BITS;

        mpfr_init2(tabx, prec);
        mpfr_init2(atanx, prec);
        mpfr_init2(y1, prec);
        mpfr_init2(y2, prec);

        for (i = 0; i < num; i++)
        {
            tt->_mp_d = (mp_ptr) arb_atan_tab21[i];
            tt->_mp_size = prec / FLINT_BITS;
            tt->_mp_alloc = tt->_mp_size;

            while (tt->_mp_size > 0 && tt->_mp_d[tt->_mp_size-1] == 0)
                tt->_mp_size--;

            mpfr_set_z(tabx, tt, MPFR_RNDD);
            mpfr_div_2ui(tabx, tabx, prec, MPFR_RNDD);

            mpfr_set_ui(atanx, i, MPFR_RNDD);
            mpfr_div_2ui(atanx, atanx, bits, MPFR_RNDD);
            mpfr_atan(atanx, atanx, MPFR_RNDD);

            mpfr_mul_2ui(y1, tabx, prec, MPFR_RNDD);
            mpfr_floor(y1, y1);
            mpfr_div_2ui(y1, y1, prec, MPFR_RNDD);

            mpfr_mul_2ui(y2, atanx, prec, MPFR_RNDD);
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
        mpfr_clear(atanx);
        mpfr_clear(y1);
        mpfr_clear(y2);
    }

    {
        slong prec, bits, num;

        prec = ARB_ATAN_TAB2_LIMBS * FLINT_BITS;
        bits = ARB_ATAN_TAB21_BITS + ARB_ATAN_TAB22_BITS;
        num = 1 << ARB_ATAN_TAB22_BITS;

        mpfr_init2(tabx, prec);
        mpfr_init2(atanx, prec);
        mpfr_init2(y1, prec);
        mpfr_init2(y2, prec);

        for (i = 0; i < num; i++)
        {
            tt->_mp_d = (mp_ptr) arb_atan_tab22[i];
            tt->_mp_size = prec / FLINT_BITS;
            tt->_mp_alloc = tt->_mp_size;

            while (tt->_mp_size > 0 && tt->_mp_d[tt->_mp_size-1] == 0)
                tt->_mp_size--;

            mpfr_set_z(tabx, tt, MPFR_RNDD);
            mpfr_div_2ui(tabx, tabx, prec, MPFR_RNDD);

            mpfr_set_ui(atanx, i, MPFR_RNDD);
            mpfr_div_2ui(atanx, atanx, bits, MPFR_RNDD);
            mpfr_atan(atanx, atanx, MPFR_RNDD);

            mpfr_mul_2ui(y1, tabx, prec, MPFR_RNDD);
            mpfr_floor(y1, y1);
            mpfr_div_2ui(y1, y1, prec, MPFR_RNDD);

            mpfr_mul_2ui(y2, atanx, prec, MPFR_RNDD);
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
        mpfr_clear(atanx);
        mpfr_clear(y1);
        mpfr_clear(y2);
    }

    TEST_FUNCTION_END(state);
}
