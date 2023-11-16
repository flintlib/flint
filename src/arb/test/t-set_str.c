/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);

    TEST_FUNCTION_END(state);
}
