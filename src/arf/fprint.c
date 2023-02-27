/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2015 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "arf.h"
#include "arb.h"

char * arf_get_str(const arf_t x, slong d)
{
    if (arf_is_special(x))
    {
        char * s = flint_malloc(5);

        if (arf_is_zero(x))
            strcpy(s, "0");
        else if (arf_is_pos_inf(x))
            strcpy(s, "+inf");
        else if (arf_is_neg_inf(x))
            strcpy(s, "-inf");
        else
            strcpy(s, "nan");

        return s;
    }
    else
    {
        arb_t t;
        *arb_midref(t) = *x;
        mag_init(arb_radref(t));  /* no need to free */
        return arb_get_str(t, FLINT_MAX(d, 1), ARB_STR_NO_RADIUS);
    }
}

void
arf_fprint(FILE * file, const arf_t x)
{
    if (arf_is_normal(x))
    {
        fmpz_t man, exp;

        fmpz_init(man);
        fmpz_init(exp);

        arf_get_fmpz_2exp(man, exp, x);

        flint_fprintf(file, "(");
        fmpz_fprint(file, man);
        flint_fprintf(file, " * 2^");
        fmpz_fprint(file, exp);
        flint_fprintf(file, ")");

        fmpz_clear(man);
        fmpz_clear(exp);
    }
    else
    {
        if (arf_is_zero(x)) flint_fprintf(file, "(0)");
        else if (arf_is_pos_inf(x)) flint_fprintf(file, "(+inf)");
        else if (arf_is_neg_inf(x)) flint_fprintf(file, "(-inf)");
        else flint_fprintf(file, "(nan)");
    }
}

void
arf_fprintd(FILE * file, const arf_t x, slong d)
{
    char * s = arf_get_str(x, d);

    fprintf(file, "%s", s);
    flint_free(s);
}
