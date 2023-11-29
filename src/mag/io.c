/*
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2015 Arb authors
    Copyright (C) 2019 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "mag.h"
#include "arf.h"

/* string *********************************************************************/

static void
mag_set_arf_dump(mag_t x, const arf_t y)
{
    if (arf_is_special(y))
    {
        if (arf_is_zero(y))
        {
            mag_zero(x);
        }
        else if (arf_is_pos_inf(y))
        {
            mag_inf(x);
        }
        else
        {
            /* a mag cannot be negative infinity or NaN */
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);
        }
    }
    else
    {
        fmpz_t mantissa, exponent;
        fmpz_init(mantissa);
        fmpz_init(exponent);

        arf_get_fmpz_2exp(mantissa, exponent, y);

        if (fmpz_cmp_ui(mantissa, 1 << MAG_BITS) >= 0)
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        mag_set_ui(x, fmpz_get_ui(mantissa));

        mag_mul_2exp_fmpz(x, x, exponent);

        fmpz_clear(exponent);
        fmpz_clear(mantissa);
    }
}

int
mag_load_str(mag_t x, const char* data)
{
    int err = 0;
    arf_t y;

    arf_init(y);

    err = arf_load_str(y, data);
    if (err)
    {
        arf_clear(y);
        return err;
    }

    mag_set_arf_dump(x, y);

    arf_clear(y);
    return err;
}

char *
mag_dump_str(const mag_t x)
{
    char * res;
    arf_t y;

    arf_init_set_mag_shallow(y, x);

    res = arf_dump_str(y);
    return res;
}

/* printing *******************************************************************/

void
mag_fprint(FILE * file, const mag_t x)
{
    flint_fprintf(file, "(");
    if (mag_is_zero(x))
        flint_fprintf(file, "0");
    else if (mag_is_inf(x))
        flint_fprintf(file, "inf");
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_sub_ui(t, MAG_EXPREF(x), MAG_BITS);
        flint_fprintf(file, "%wu * 2^", MAG_MAN(x));
        fmpz_fprint(file, t);
        fmpz_clear(t);
    }
    flint_fprintf(file, ")");
}

void
mag_fprintd(FILE * file, const mag_t x, slong d)
{
    arf_t t;
    arf_init(t);
    arf_set_mag(t, x);
    arf_fprintd(file, t, d);
    arf_clear(t);
}

void mag_print(const mag_t x) { mag_fprint(stdout, x); }
void mag_printd(const mag_t x, slong d) { mag_fprintd(stdout, x, d); }

/* file I/O *******************************************************************/

int
mag_load_file(mag_t x, FILE * stream)
{
    int err = 0;
    arf_t y;

    arf_init(y);

    err = arf_load_file(y, stream);

    if (err)
    {
        arf_clear(y);
        return err;
    }

    mag_set_arf_dump(x, y);

    arf_clear(y);
    return err;
}

int
mag_dump_file(FILE * stream, const mag_t x)
{
    int nwrite;
    char* data = mag_dump_str(x);

    nwrite = fputs(data, stream);
    if (nwrite == EOF)
        return nwrite;

    flint_free(data);
    return 0;
}
