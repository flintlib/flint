/*
    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Arb authors
    Copyright (C) 2019 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>
#include "arb.h"
#include "arf.h"
#include "mag.h"

/* strings ********************************************************************/

char *
arb_dump_str(const arb_t x)
{
    char * mid;
    char * mag;
    size_t res_len;
    char * res;

    mid = arf_dump_str(arb_midref(x));
    mag = mag_dump_str(arb_radref(x));

    res_len = strlen(mid) + 1 + strlen(mag);
    res = (char*)flint_malloc(res_len + 1);
    strcpy(res, mid);
    strcat(res, " ");
    strcat(res, mag);

    flint_free(mid);
    flint_free(mag);

    return res;
}

int
arb_load_str(arb_t x, const char* data)
{
    size_t midlen, maglen;
    char * mid;
    char * mag;
    int err = 0;

    const char* split = strchr(data, ' ');
    if (split == NULL)
    {
        return 1;
    }
    split = strchr(split + 1, ' ');
    if (split == NULL)
    {
        return 1;
    }

    midlen = (size_t)(split - data);
    mid = (char*)flint_malloc(midlen + 1);
    strncpy(mid, data, midlen);
    mid[midlen] = '\0';

    maglen = strlen(data) - midlen - 1;
    mag = (char*)flint_malloc(maglen + 1);
    strncpy(mag, split + 1, maglen);
    mag[maglen] = '\0';

    err = arf_load_str(arb_midref(x), mid);
    if (err)
    {
        flint_free(mid);
        flint_free(mag);
        return err;
    }

    err = mag_load_str(arb_radref(x), mag);

    flint_free(mid);
    flint_free(mag);

    return err;
}

/* printing *******************************************************************/

void
arb_fprint(FILE * file, const arb_t x)
{
    arf_fprint(file, arb_midref(x));
    flint_fprintf(file, " +/- ");
    mag_fprint(file, arb_radref(x));
}

void
arb_fprintd(FILE * file, const arb_t x, slong digits)
{
    arf_fprintd(file, arb_midref(x), FLINT_MAX(digits, 1));
    flint_fprintf(file, " +/- ");
    mag_fprintd(file, arb_radref(x), 5);
}

void
arb_fprintn(FILE * file, const arb_t x, slong digits, ulong flags)
{
    char * s = arb_get_str(x, digits, flags);
    flint_fprintf(file, "%s", s);
    flint_free(s);
}

void arb_print(const arb_t x) { arb_fprint(stdout, x); }
void arb_printd(const arb_t x, slong digits) { arb_fprintd(stdout, x, digits); }
void arb_printn(const arb_t x, slong digits, ulong flags) { arb_fprintn(stdout, x, digits, flags); }

void
_arb_vec_printn(arb_srcptr vec, slong len, slong ndigits, ulong flags)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        arb_printn(vec + i, ndigits, flags);
        if (i < len - 1)
            flint_printf(", ");
    }
}

void
_arb_vec_printd(arb_srcptr vec, slong len, slong ndigits)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        arb_printd(vec + i, ndigits);
        if (i < len - 1)
        {
            flint_printf(", ");
        }
    }
    flint_printf("\n");
}


/* file I/O *******************************************************************/

int
arb_dump_file(FILE* stream, const arb_t x)
{
    int nwrite;
    char* data = arb_dump_str(x);

    nwrite = fputs(data, stream);
    if (nwrite == EOF)
        return nwrite;

    flint_free(data);
    return 0;
}

int
arb_load_file(arb_t x, FILE* stream)
{
    int err;

    err = arf_load_file(arb_midref(x), stream);

    if (err) return err;

    err = mag_load_file(arb_radref(x), stream);

    return err;
}
