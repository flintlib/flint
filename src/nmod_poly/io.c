/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Jean-Pierre Flori

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "nmod_poly.h"

/* printing *******************************************************************/

int nmod_poly_fprint(FILE * f, const nmod_poly_t poly)
{
    char * s;
    int r;

    s = nmod_poly_get_str(poly);
    r = fputs(s, f);
    flint_free(s);

    return (r < 0) ? r : 1;
}

int nmod_poly_fprint_pretty(FILE * f, const nmod_poly_t a, const char * x)
{
    size_t r;
    slong i;

    if (a->length == 0)
    {
        r = fputc('0', f);
        r = ((int) r != EOF) ? 1 : EOF;
        return r;
    }
    else if (a->length == 1)
    {
        r = flint_fprintf(f, "%wu", a->coeffs[0]);
        return r;
    }

    i = a->length - 1;
    r = 1;
    if (i == 1)
    {
        switch (a->coeffs[1])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "%s", x);
                break;
            default:
                r = flint_fprintf(f, "%wu*%s", a->coeffs[1], x);
        }
        --i;
    }
    else
    {
        switch (a->coeffs[i])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "%s^%wd", x, i);
                break;
            default:
                r = flint_fprintf(f, "%wu*%s^%wd", a->coeffs[i], x, i);
        }
        --i;
    }
    for (; (r > 0) && (i > 1); --i)
    {
        switch (a->coeffs[i])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "+%s^%wd", x, i);
                break;
            default:
                r = flint_fprintf(f, "+%wu*%s^%wd", a->coeffs[i], x, i);
        }
    }
    if (r > 0 && i == 1)
    {
        switch (a->coeffs[1])
        {
            case UWORD(0):
                break;
            case UWORD(1):
                r = flint_fprintf(f, "+%s", x);
                break;
            default:
                r = flint_fprintf(f, "+%wu*%s", a->coeffs[1], x);
        }
    }
    if (r > 0)
    {
        if (a->coeffs[0] != UWORD(0))
            r = flint_fprintf(f, "+%wu", a->coeffs[0]);
    }

    return (int) r;
}

int nmod_poly_print(const nmod_poly_t a) { return nmod_poly_fprint(stdout, a); }
int nmod_poly_print_pretty(const nmod_poly_t a, const char * x) { return nmod_poly_fprint_pretty(stdout, a, x); }

/* reading ********************************************************************/

int nmod_poly_fread(FILE * f, nmod_poly_t poly)
{
    slong i, length;
    mp_limb_t n;

    if (flint_fscanf(f, "%wd %wu", &length, &n) != 2)
        return 0;

    nmod_poly_clear(poly);
    nmod_poly_init(poly,n);
    nmod_poly_fit_length(poly, length);
    poly->length = length;

    for (i = 0; i < length; i++)
    {
        if (!flint_fscanf(f, "%wu", &poly->coeffs[i]))
        {
            poly->length = i;
            return 0;
        }
    }

    _nmod_poly_normalise(poly);

    return 1;
}

int nmod_poly_read(nmod_poly_t poly) { return nmod_poly_fread(stdin, poly); }
