/*
    Copyright (C) 2019 Julian Rüth

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>

#include "arf.h"

static void
arf_get_fmpz_2exp_dump(fmpz_t m, fmpz_t e, const arf_t x) {
    if (arf_is_special(x))
    {
        fmpz_zero(m);
        if (arf_is_zero(x)) fmpz_zero(e);
        else if (arf_is_pos_inf(x)) fmpz_set_si(e, -1);
        else if (arf_is_neg_inf(x)) fmpz_set_si(e, -2);
        else if (arf_is_nan(x)) fmpz_set_si(e, -3);
        else
        {
            /* Impossible to happen; all the special values have been treated above. */
            flint_abort();
        }
        return;
    }

    arf_get_fmpz_2exp(m, e, x);
}

char *
arf_dump_str(const arf_t x)
{
    size_t res_len;
    char * res;

    fmpz_t mantissa, exponent;

    fmpz_init(mantissa);
    fmpz_init(exponent);

    arf_get_fmpz_2exp_dump(mantissa, exponent, x);

    res_len = (fmpz_sgn(mantissa) < 0) + fmpz_sizeinbase(mantissa, 16) + 1
        + (fmpz_sgn(exponent) < 0) + fmpz_sizeinbase(exponent, 16);
    res = (char*)flint_malloc(res_len + 1);

    fmpz_get_str(res, 16, mantissa);
    strcat(res, " ");
    fmpz_get_str(res + strlen(res), 16, exponent);

    fmpz_clear(mantissa);
    fmpz_clear(exponent);

    if(strlen(res) != res_len) flint_abort(); /* assert */

    return res;
}
