/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <stdlib.h>
#include <string.h>
#include "arb.h"

const char * testdata_floats[] = {
    "0", /* repeated to test empty string later */
    "0",
    "0.0",
    "0.",
    ".0",
    "+0",
    "+0.0",
    "+0.",
    "   +.0  ",
    "-0",
    "-0.0",
    "-0.",
    "-.0",

    " 0e3",
    "0.0e3",
    "0.e3",
    " .0e3",
    "+0e3",
    "+0.0e3",
    "+0.e3",
    "+.0e3",
    "-0e3",
    "-0.0e3",
    "-0.e3",
    "-.0e3",

    "0e+3",
    "0.0e+3",
    "0.e+3",
    "  .0e+3",
    "+0e+3",
    "+0.0e+3",
    "+0.e+3",
    "+.0e+3",
    "-0E+3",
    "-0.0E+3   ",
    "-0.e+3",
    "-.0e+3",

    "0e-3",
    "0.0e-3",
    "0.e-3",
    ".0E-3",
    "+0e-3",
    "+0.0e-3",
    "+0.E-3",
    "+.0e-3",
    "-0e-3",
    "-0.0e-3",
    "-0.e-3",
    "-.0e-3",

    "03.125",
    "+03.125",
    "-03.125",
    "03.12500",
    "+03.12500",
    "-03.12500",
    "03.125e+3",
    "+03.125e+3",
    "  -03.125e+3  ",
    "03.12500e+3",
    "+03.12500e+3",
    "-03.12500e+3",
    "03.125e3",
    "+03.125E3",
    "-03.125e3",
    "03.12500e3",
    "+03.12500E3",
    "-03.12500E3",

    "25000.0e-2",
    "-25000.0e-2",
    "   25000.0000000000000000000e-2",
    "-25000.0000000000000000000e-2",
    " 000025000.0000000000000000000e-2 ",
    "-000025000.0000000000000000000e-2",
    "25000.e-2",
    "-25000.e-2",

    "12345.125",
    "-12345.125",
    "+12345.125",

    "1.0625",
    "1.03125",
    "1.015625",
    "1.0078125",
    "1.00390625",
    "1.001953125",
    "1.0009765625",
    "1.00048828125",
    "1.000244140625",
    "1.0001220703125",
    "1.00006103515625",
    "1.000030517578125",
    "1.0000152587890625",
    "1.00000762939453125",
    "1.000003814697265625",
    "1.0000019073486328125",
    "1.00000095367431640625",
    "1.000000476837158203125",
    "1.0000002384185791015625",
    "1.00000011920928955078125",
    "1.000000059604644775390625",
    "1.0000000298023223876953125",
    "1.00000001490116119384765625",
    "1.000000007450580596923828125",
    "1.0000000037252902984619140625",
    "1.00000000186264514923095703125",
    "1.000000000931322574615478515625",
    "1.0000000004656612873077392578125",
    "1.00000000023283064365386962890625",
    "1.000000000116415321826934814453125",
    "1.0000000000582076609134674072265625",
    "1.00000000002910383045673370361328125",
    "1.000000000014551915228366851806640625",
    "1.0000000000072759576141834259033203125",
    "1.00000000000363797880709171295166015625",
    "1.000000000001818989403545856475830078125",
    "1.0000000000009094947017729282379150390625",
    "1.00000000000045474735088646411895751953125",
    "1.000000000000227373675443232059478759765625",
    "1.0000000000001136868377216160297393798828125",
    "1.00000000000005684341886080801486968994140625",
    "1.000000000000028421709430404007434844970703125",
    "1.0000000000000142108547152020037174224853515625",
    "1.00000000000000710542735760100185871124267578125",
    "1.000000000000003552713678800500929355621337890625",
    "1.0000000000000017763568394002504646778106689453125",
    "1.00000000000000088817841970012523233890533447265625",
    "1.000000000000000444089209850062616169452667236328125",
    "1.0000000000000002220446049250313080847263336181640625",

    "inf",
    "-inf",
    "+inf",
    "Inf",
    "-INF",
    "+Inf",

    "NAN",
    "-NaN",
    "+NAN",

    NULL,
};

const char * testdata_invalid[] = {
    "",
    ".",
    "+.",
    "-.",
    ".e+3",
    "-.e+5",
    "+.e-5",
    "2+3",
    "12 34",
    "150a.25",
    "-e+4",
    "10.25x",
    "10.3.5",
    "125e3.6",
    "125e-3.6",
    "3.14 e+5",
    "3.140 e3",
    "3.14+e5",
    "3.14e+ 5",
    " 3.14e- 5",
    "3.14e+-5",
    "3.14e+/-5",
    ".0.",
    "..",
    ":)",
    "  +/- ",
    " 0 0",
    "  +/- EEE ",
    "-3.5e+x5 +/-",
    "4.7 +/- -3.5e+x5",
    "4.7 +/-",
    NULL,
};

const char * testdata_rounding[] = {
    "1.1",
    "-1.1",
    "1.01",
    "1.001",
    "1.0001",
    "1.00001",
    "1.000001",
    "1.0000001",
    "1.00000001",
    "1.000000001",
    "1.0000000001",
    "1.00000000001",
    "1.000000000001",
    "1.0000000000001",
    "1.00000000000001",
    "1.000000000000001",
    "1.0000000000000001",
    "1.00000000000000001",
    "100000000000000001",
    "1000000000000000001",
    "10000000000000000001",
    NULL,
};

TEST_FUNCTION_START(arb_set_str, state)
{
    arb_t t, u, v;
    double x;
    int error, bracket;
    char tmp[256];
    slong i, j;

    arb_init(t);
    arb_init(u);
    arb_init(v);

    for (i = 0; testdata_floats[i] != NULL; i++)
    {
        arb_const_pi(t, 53);

        error = arb_set_str(t, testdata_floats[i], 53);

        x = strtod(testdata_floats[i], NULL);

        if (x != x)
        {
            arb_indeterminate(u);
        }
        else
        {
            arf_set_d(arb_midref(u), x);
            mag_zero(arb_radref(u));
        }

        if (error != 0 || !arb_equal(t, u))
        {
            flint_printf("FAIL (valid input): %s\n", testdata_floats[i]);
            arb_printd(t, 15); flint_printf("\n");
            arb_printd(u, 15); flint_printf("\n");
            flint_abort();
        }
    }

    for (i = 0; testdata_floats[i] != NULL; i++)
    {
        for (j = 0; testdata_floats[j] != NULL; j++)
        {
            for (bracket = 0; bracket < 2; bracket++)
            {
                arb_const_pi(t, 53);

                bracket = n_randint(state, 2);

                strcpy(tmp, "");

                if (bracket)
                    strcat(tmp, "[");

                /* allow empty string for midpoint */
                strcat(tmp, (i == 0) ? "" : testdata_floats[i]);
                strcat(tmp, "+/-");
                strcat(tmp, testdata_floats[j]);

                if (bracket)
                    strcat(tmp, "]");

                error = arb_set_str(t, tmp, 53);

                x = strtod((i == 0) ? "0" : testdata_floats[i], NULL);

                if (x != x)
                {
                    arb_indeterminate(u);
                }
                else
                {
                    arf_set_d(arb_midref(u), x);
                    mag_zero(arb_radref(u));
                }

                x = strtod(testdata_floats[j], NULL);
                arf_set_d(arb_midref(v), x);
                mag_zero(arb_radref(v));

                arb_abs(v, v);
                arb_add_error(u, v);

                if (error != 0 || !arb_equal(t, u))
                {
                    flint_printf("FAIL (valid input): %s\n", tmp);
                    arb_printd(t, 15); flint_printf("\n");
                    arb_printd(u, 15); flint_printf("\n");
                    flint_abort();
                }
            }
        }
    }

    for (i = 0; testdata_invalid[i] != NULL; i++)
    {
        arb_const_pi(t, 53);

        error = arb_set_str(t, testdata_invalid[i], 53);

        if (error == 0)
        {
            flint_printf("FAIL (invalid input): %s\n", testdata_invalid[i]);
            arb_printd(t, 15); flint_printf("\n");
            flint_abort();
        }
    }

    for (i = 0; testdata_invalid[i] != NULL; i++)
    {
        for (j = 0; testdata_invalid[j] != NULL; j++)
        {
            for (bracket = 0; bracket < 2; bracket++)
            {
                arb_const_pi(t, 53);

                bracket = n_randint(state, 2);

                strcpy(tmp, "");

                if (bracket)
                    strcat(tmp, "[");

                strcat(tmp, testdata_invalid[i]);
                strcat(tmp, "+/-");
                strcat(tmp, testdata_invalid[j]);

                if (bracket)
                    strcat(tmp, "]");

                error = arb_set_str(t, tmp, 53);

                if (error == 0)
                {
                    flint_printf("FAIL (invalid input): %s\n", tmp);
                    arb_printd(t, 15); flint_printf("\n");
                    flint_abort();
                }
            }
        }
    }

    for (i = 0; testdata_rounding[i] != NULL; i++)
    {
        arb_const_pi(t, 53);

        error = arb_set_str(t, testdata_rounding[i], 53);

        /* It is not *guaranteed* that the rounding will be to half-ulp
         * accuracy, just like it is not guaranteed that all inputs will be
         * correctly rounded. But it is very likely for a random string, and
         * it *is* guaranteed if the input is the result of get_str performed
         * with a sufficient number of digits.
         * Thus, we are justified in testing for these things here. */
        x = strtod(testdata_rounding[i], NULL);
        arf_set_d(arb_midref(u), x);
        arf_mag_set_ulp(arb_radref(u), arb_midref(u), 54);

        if (error != 0 || !arf_equal(arb_midref(t), arb_midref(u)) ||
            mag_cmp(arb_radref(t), arb_radref(u)) > 0)
        {
            flint_printf("FAIL (valid input): %s\n", testdata_rounding[i]);
            arb_printd(t, 19); flint_printf("\n");
            arb_printd(u, 19); flint_printf("\n");
            flint_abort();
        }

        error = arb_set_str(t, testdata_rounding[i], 24);

        x = strtof(testdata_rounding[i], NULL);
        arf_set_d(arb_midref(u), x);
        arf_mag_set_ulp(arb_radref(u), arb_midref(u), 25);

        if (error != 0 || !arf_equal(arb_midref(t), arb_midref(u)) ||
            mag_cmp(arb_radref(t), arb_radref(u)) > 0)
        {
            flint_printf("FAIL (valid input): %s\n", testdata_rounding[i]);
            arb_printd(t, 19); flint_printf("\n");
            arb_printd(u, 19); flint_printf("\n");
            flint_abort();
        }
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);

    TEST_FUNCTION_END(state);
}
