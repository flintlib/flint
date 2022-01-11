/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz_vec.h"
#include "profiler.h"

#define MAXLENGTH 100

typedef struct
{
   slong bits;
}
info_t;

void
_fmpz_vec_content_old(fmpz_t res, const fmpz * vec, slong len)
{
    while (len > 0 && fmpz_is_zero(vec + 0))
    {
        len--;
        vec++;
    }

    while (len > 1 && fmpz_is_zero(vec + len - 1))
        len--;

    if (len <= 2)
    {
        if (len == 0)
            fmpz_zero(res);
        else if (len == 1)
            fmpz_abs(res, vec + 0);
        else
            fmpz_gcd(res, vec + 0, vec + 1);
        return;
    }

    if (fmpz_is_pm1(vec + 0) || fmpz_is_pm1(vec + len - 1))
    {
        fmpz_one(res);
        return;
    }

    fmpz_gcd3(res, vec + 0, vec + 1, vec + len - 1);
    vec += 2;
    len -= 3;

    while (len >= 2 && !fmpz_is_one(res))
    {
        fmpz_gcd3(res, vec + 0, vec + len - 1, res);
        vec++;
        len -= 2;
    }

    if (len != 0 && !fmpz_is_one(res))
        fmpz_gcd(res, res, vec + 0);
}

void
sample_new(void * arg, ulong count)
{
    fmpz_t res, multiplier;
    fmpz * vec;
    slong len;
    int ix;
    info_t * info = (info_t *) arg;
    slong bits = info->bits;

    fmpz_init(res);
    fmpz_init(multiplier);
    vec = _fmpz_vec_init(MAXLENGTH);

    FLINT_TEST_INIT(state);

    for (ix = 0; ix < 10000 * count; ix++)
    {
        len = n_randint(state, MAXLENGTH);
        while (len == 0)
            len = n_randint(state, MAXLENGTH);

        _fmpz_vec_randtest(vec, state, len, bits);
        fmpz_randtest_not_zero(multiplier, state, bits);

        _fmpz_vec_scalar_mul_fmpz(vec, vec, len, multiplier);

        prof_start();
        _fmpz_vec_content(res, vec, len);
        prof_stop();
    }

    fmpz_clear(res);
    fmpz_clear(multiplier);
    _fmpz_vec_clear(vec, MAXLENGTH);

    flint_randclear(state);
}

void
sample_old(void * arg, ulong count)
{
    fmpz_t res, multiplier;
    fmpz * vec;
    slong len;
    int ix;
    info_t * info = (info_t *) arg;
    slong bits = info->bits;

    fmpz_init(res);
    fmpz_init(multiplier);
    vec = _fmpz_vec_init(MAXLENGTH);

    FLINT_TEST_INIT(state);

    for (ix = 0; ix < 10000 * count; ix++)
    {
        len = n_randint(state, MAXLENGTH);
        while (len == 0)
            len = n_randint(state, MAXLENGTH);

        _fmpz_vec_randtest(vec, state, len, bits);
        fmpz_randtest_not_zero(multiplier, state, bits);

        _fmpz_vec_scalar_mul_fmpz(vec, vec, len, multiplier);

        prof_start();
        _fmpz_vec_content_old(res, vec, len);
        prof_stop();
    }

    fmpz_clear(res);
    fmpz_clear(multiplier);
    _fmpz_vec_clear(vec, MAXLENGTH);

    flint_randclear(state);
}

int
main()
{
    double minnew, maxnew, minold, maxold;
    int bits;
    info_t as;

    for (bits = 5; bits <= 150; bits += 5)
    {
        as.bits = bits;

        prof_repeat(&minnew, &maxnew, sample_new, &as);
        prof_repeat(&minold, &maxold, sample_old, &as);
        flint_printf("%3d bits:   (min) %.2fx speedup   (max) %.2fx speedup\n",
                bits, minold / minnew, maxold / maxnew);
    }
}
