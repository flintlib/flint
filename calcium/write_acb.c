/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "calcium.h"

static char * arb_get_str2(const arb_t x, slong digits, ulong flags)
{
    char * s;

    if (!arb_is_finite(x) && (flags & ARB_STR_NO_RADIUS))
    {
        s = flint_malloc(4);
        strcpy(s, "???");
        return s;
    }

    s = arb_get_str(x, digits, flags);

    if ((flags & ARB_STR_NO_RADIUS) && (s[0] == '['))
    {
        slong i, j;
        fmpz_t e;

        fmpz_init(e);
        for (i = 1; s[i] != '\0'; i++)
        {
            if (s[i] == 'e')
            {
                for (j = i + 1; s[j] != '\0'; j++)
                {
                    if (s[j] == ']')
                    {
                        s[j] = '\0';
                        break;
                    }
                }

                fmpz_set_str(e, s + i + 1, 10);
                fmpz_add_ui(e, e, 1);
                s[0] = '0';
                s[1] = 'e';
                if (fmpz_sgn(e) < 0)
                {
                    fmpz_get_str(s + 2, 10, e);
                }
                else
                {
                    s[3] = '-';
                    fmpz_get_str(s + 3, 10, e);
                }
                break;

            }
        }
        fmpz_clear(e);
    }

    return s;
}

void
calcium_write_acb(calcium_stream_t out, const acb_t z, slong digits, ulong flags)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        calcium_write_free(out, arb_get_str2(acb_realref(z), digits, flags));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        calcium_write_free(out, arb_get_str2(acb_imagref(z), digits, flags));
        calcium_write(out, "*I");
    }
    else
    {
        calcium_write_free(out, arb_get_str2(acb_realref(z), digits, flags));

        if ((arb_is_exact(acb_imagref(z)) || (flags & ARB_STR_NO_RADIUS))
                && arf_sgn(arb_midref(acb_imagref(z))) < 0)
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, acb_imagref(z));
            calcium_write(out, " - ");
            calcium_write_free(out, arb_get_str2(t, digits, flags));
            arb_clear(t);
        }
        else
        {
            calcium_write(out, " + ");
            calcium_write_free(out, arb_get_str2(acb_imagref(z), digits, flags));
        }

        calcium_write(out, "*I");
    }
}
