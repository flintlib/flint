/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

static slong
bisect(slong needle, const slong * haystack, slong len)
{
    slong a, b, mid;

    a = 0;
    b = len - 1;

    while (a < b)
    {
        mid = (a + b) / 2;

        if (haystack[mid] < needle)
            a = mid + 1;
        else
            b = mid;
    }

    if (a == b && haystack[a] == needle)
    {
        return a;
    }

    return -1;
}

/* write p = 2a with a in P */
static int
write_as_2a(slong * i1, slong * i2, slong p, const slong * P, slong Plen)
{
    if (p % 2 == 0)
    {
        slong i = bisect(p / 2, P, Plen);

        if (i != -1)
        {
            *i1 = *i2 = i;
            return 1;
        }
    }

    return 0;
}

/* write p = a + b with a, b in P */
static int
write_as_a_b(slong * i1, slong * i2, slong p, const slong * P, slong Plen)
{
    slong i, j, pi;

    for (i = 0; i < Plen; i++)
    {
        pi = P[i];

        j = bisect(p - pi, P, Plen);

        if (j != -1)
        {
            *i1 = i;
            *i2 = j;
            return 1;
        }
    }

    return 0;
}

/* write p = 2a + b with a, b in P */
static int
write_as_2a_b(slong * i1, slong * i2, slong p, const slong * P, slong Plen)
{
    slong i, j, pi;

    for (i = 0; i < Plen; i++)
    {
        pi = P[i];

        if (2 * pi > p)
            break;

        j = bisect(p - 2*pi, P, Plen);

        if (j != -1)
        {
            *i1 = i;
            *i2 = j;
            return 1;
        }
    }

    return 0;
}

void
acb_modular_addseq_theta(slong * exponents, slong * aindex, slong * bindex, slong num)
{
    slong i;
    slong c;

    for (i = 0; i < num; i++)
    {
        if (i == 0)
        {
            exponents[i] = 1;
            aindex[i] = 0;
            bindex[i] = 0;
            continue;
        }
        else
        {
            c = ((i + 2) * (i + 2)) / 4;

            exponents[i] = c;

            if (write_as_2a(aindex + i, bindex + i, c, exponents, i))
                continue;

            if (write_as_a_b(aindex + i, bindex + i, c, exponents, i))
                continue;

            if (write_as_2a_b(aindex + i, bindex + i, c, exponents, i))
                continue;

            flint_throw(FLINT_ERROR, "i = %wd, c = %wu: bad addition sequence!\n", i, c);
        }
    }
}

void
acb_modular_addseq_eta(slong * exponents, slong * aindex, slong * bindex, slong num)
{
    slong i;
    slong c;

    for (i = 0; i < num; i++)
    {
        if (i == 0)
        {
            exponents[i] = 1;
            aindex[i] = 0;
            bindex[i] = 0;
            continue;
        }
        else
        {
            c = ((i + 2) / 2) * ((3 * i + 5) / 2) / 2;

            exponents[i] = c;

            if (write_as_2a(aindex + i, bindex + i, c, exponents, i))
                continue;

            if (write_as_a_b(aindex + i, bindex + i, c, exponents, i))
                continue;

            if (write_as_2a_b(aindex + i, bindex + i, c, exponents, i))
                continue;

            flint_throw(FLINT_ERROR, "i = %wd, c = %wu: bad addition sequence!\n", i, c);
        }
    }
}

