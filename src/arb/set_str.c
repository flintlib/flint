/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <ctype.h>
#include "arb.h"

static int
arb_set_float_str(arb_t res, const char * inp, slong prec)
{
    char * emarker;
    char * buf;
    int error;
    slong i;
    fmpz_t exp;
    fmpz_t man;
    slong num_int, num_frac;
    int after_radix;

    if (inp[0] == '+')
    {
        return arb_set_float_str(res, inp + 1, prec);
    }

    if (inp[0] == '-')
    {
        error = arb_set_float_str(res, inp + 1, prec);
        arb_neg(res, res);
        return error;
    }

    if (strcmp(inp, "inf") == 0)
    {
        arb_pos_inf(res);
        return 0;
    }

    if (strcmp(inp, "nan") == 0)
    {
        arb_indeterminate(res);
        return 0;
    }

    error = 0;
    fmpz_init(exp);
    fmpz_init(man);
    buf = flint_malloc(strlen(inp) + 1);

    emarker = strchr(inp, 'e');

    /* parse exponent (0 by default) */
    if (emarker != NULL)
    {
        /* allow e+42 as well as e42 */
        if (emarker[1] == '+')
        {
            if (!(emarker[2] >= '0' && emarker[2] <= '9'))
                error = 1;
            else
                error = fmpz_set_str(exp, emarker + 2, 10);
        }
        else
        {
            if (emarker[1] == '-' && !(emarker[2] >= '0' && emarker[2] <= '9'))
                error = 1;
            else
                error = fmpz_set_str(exp, emarker + 1, 10);
        }

        if (error)
            goto cleanup;
    }

    /* parse floating-point part */
    {
        num_int = 0;
        num_frac = 0;
        after_radix = 0;

        for (i = 0; inp + i != emarker && inp[i] != '\0'; i++)
        {
            if (inp[i] == '.' && !after_radix)
            {
                after_radix = 1;
            }
            else if (inp[i] >= '0' && inp[i] <= '9')
            {
                buf[num_int + num_frac] = inp[i];

                num_frac += after_radix;
                num_int += !after_radix;
            }
            else
            {
                error = 1;
                goto cleanup;
            }
        }

        buf[num_int + num_frac] = '\0';

        /* put trailing zeros into the exponent */
        while (num_int + num_frac > 1 && buf[num_int + num_frac - 1] == '0')
        {
            buf[num_int + num_frac - 1] = '\0';
            num_frac--;
        }

        fmpz_sub_si(exp, exp, num_frac);

        error = fmpz_set_str(man, buf, 10);
        if (error)
            goto cleanup;
    }

    if (fmpz_is_zero(man))
    {
        arb_zero(res);
    }
    else if (fmpz_is_zero(exp))
    {
        arb_set_round_fmpz(res, man, prec);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_set_ui(t, 10);
        arb_set_fmpz(res, man);

        if (fmpz_sgn(exp) > 0)
        {
            arb_pow_fmpz_binexp(t, t, exp, prec + 4);
            arb_mul(res, res, t, prec);
        }
        else
        {
            fmpz_neg(exp, exp);
            arb_pow_fmpz_binexp(t, t, exp, prec + 4);
            arb_div(res, res, t, prec);
        }

        arb_clear(t);
    }

cleanup:
    fmpz_clear(exp);
    fmpz_clear(man);
    flint_free(buf);

    if (error)
        arb_indeterminate(res);

    return error;
}

int
arb_set_str(arb_t res, const char * inp, slong prec)
{
    char * buf;
    char * split;
    char * first;
    char * last;
    slong i, len;
    int error;

    error = 0;
    len = strlen(inp);
    buf = flint_malloc(len + 1);

    for (i = 0; i <= len; i++)
        buf[i] = tolower(inp[i]);

    split = strstr(buf, "+/-");

    if (split == NULL)
    {
        /* strip whitespace and brackets */
        first = buf;
        while (isspace(first[0]) || first[0] == '[')
            first++;
        last = buf + len;
        while (last - first > 0 && (isspace(last[-1]) || last[-1] == ']'))
            last--;
        last[0] = '\0';

        error = arb_set_float_str(res, first, prec);
    }
    else
    {
        arb_t rad;
        arb_init(rad);

        /* strip whitespace and brackets */
        first = buf;
        while (isspace(first[0]) || first[0] == '[')
            first++;
        last = split;
        while (last - first > 0 && (isspace(last[-1]) || last[-1] == ']'))
            last--;
        last[0] = '\0';

        if (first == last)
            arb_zero(res);
        else
            error = arb_set_float_str(res, first, prec);

        if (!error)
        {
            /* strip whitespace and brackets */
            first = split + 3;
            while (isspace(first[0]) || first[0] == '[')
                first++;
            last = buf + len;
            while (last - first > 0 && (isspace(last[-1]) || last[-1] == ']'))
                last--;
            last[0] = '\0';

            error = arb_set_float_str(rad, first, prec);
            arb_abs(rad, rad);
            arb_add_error(res, rad);
        }

        arb_clear(rad);
    }

    flint_free(buf);
    return error;
}

