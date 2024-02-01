/*
    Copyright (C) 2018 Vincent Delecroix
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <errno.h>
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

/* TODO: Remove me */
static void
__fmpq_poly_set_array_mpq(fmpz * poly, fmpz_t den, const mpq_t * a, slong n)
{
    slong i;
    mpz_t d, t;

    flint_mpz_init_set_ui(d, 1);
    mpz_init(t);
    for (i = 0; i < n; i++)
    {
        mpz_lcm(d, d, mpq_denref(a[i]));
    }

    for (i = 0; i < n; i++)
    {
        __mpz_struct *ptr = _fmpz_promote(poly + i);

        mpz_divexact(t, d, mpq_denref(a[i]));
        mpz_mul(ptr, mpq_numref(a[i]), t);
        _fmpz_demote_val(poly + i);
    }

    fmpz_set_mpz(den, d);
    mpz_clear(d);
    mpz_clear(t);
}

int
_fmpq_poly_set_str(fmpz * poly, fmpz_t den, const char * str, slong len)
{
    char * w;
    slong i;
    mpq_t * a;

    if (!len)
        return *str == '\0';
    if (*str == '\0')
        return -1;

    /* Find maximal gap between spaces and allocate w */
    {
        const char * s = str;
        slong max;
        for (max = 0; *s != '\0';)
        {
            slong cur;
            for (s++, cur = 1; *s != ' ' && *s != '\0'; s++, cur++) ;
            if (max < cur)
                max = cur;
        }

        w = (char *) flint_malloc((max + 1) * sizeof(char));
    }

    a = (mpq_t *) flint_malloc(len * sizeof(mpq_t));

    str--;
    for (i = 0; i < len; i++)
    {
        char * v;
        int ans;

        for (str++, v = w; *str != ' ' && *str != '\0';)
            *v++ = *str++;
        *v = '\0';
        mpq_init(a[i]);
        ans = mpq_set_str(a[i], w, 10);

        /* If the format is not correct, clear up and return -1 */
        if (ans || (i + 1 < len && *str == '\0'))
        {
            int j;
            for (j = 0; j <= i; j++)
                mpq_clear(a[j]);
            flint_free(a);
            flint_free(w);
            return -1;
        }
    }

    __fmpq_poly_set_array_mpq(poly, den, (const mpq_t *) a, len);

    for (i = 0; i < len; i++)
        mpq_clear(a[i]);
    flint_free(a);
    flint_free(w);

    if (*str != '\0')
        return -1;
    else
        return 0;
}

int
fmpq_poly_set_str(fmpq_poly_t poly, const char * str)
{
    int ans;
    slong len;
    char * endptr;

    /* Get the length (a positive integer) */
    if (*str < 48 || *str > 57)
    {
        fmpq_poly_zero(poly);
        return -1;
    }
    errno = 0;
    len = strtol(str, &endptr, 10);
    if (errno || len < 0 || (len > 0 && *endptr == '\0'))
    {
        fmpq_poly_zero(poly);
        return -1;
    }
    if (len == 0)
    {
        fmpq_poly_zero(poly);
        if (*(str+1) != '\0')
            return -1;
        return 0;
    }

    /* Check that we have two spaces after the length */
    endptr++;
    if (*endptr != ' ')
    {
        fmpq_poly_zero(poly);
        return -1;
    }
    endptr++;

    /* Now get the coefficients */
    fmpq_poly_fit_length(poly, len);
    ans = _fmpq_poly_set_str(poly->coeffs, poly->den, endptr, len);

    if (ans)
    {
        fmpq_poly_zero(poly);
        return -1;
    }

    _fmpq_poly_set_length(poly, len);
    _fmpq_poly_normalise(poly);
    return 0;
}
