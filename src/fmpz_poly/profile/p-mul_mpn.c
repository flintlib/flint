/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "profiler.h"

/* Usage: p-mul_mpn [mode] [maxlen]
   Compares full products: mpn (auto), classical, toom_karatsuba, toom_scalar,
   and the current dispatcher, for a grid of bit sizes and lengths.
   mode: 0 = random signed, 1 = nonnegative, 2 = squaring */

typedef void (* mulfunc)(fmpz *, const fmpz *, slong, const fmpz *, slong);

static void
toom(fmpz * r, const fmpz * a, slong la, const fmpz * b, slong lb)
{
    if (!_fmpz_poly_mul_toom_scalar(r, a, la, b, lb))
        _fmpz_poly_mul_classical(r, a, la, b, lb);
}

static double
tm(mulfunc f, fmpz * r, const fmpz * a, slong la, const fmpz * b, slong lb)
{
    double t, w;
    TIMEIT_START
    f(r, a, la, b, lb);
    TIMEIT_STOP_VALUES(t, w);
    (void) w;
    return t * 1e9;
}

int main(int argc, char * argv[])
{
    slong mode, maxlen, bits, len, i, j;
    slong bl[] = {40, 64, 100, 128, 192, 256, 320, 400, 512, 700, 1000, 1500, 2000, 3000, 4000, 6000};
    slong nb = sizeof(bl) / sizeof(bl[0]);
    FLINT_TEST_INIT(state);

    mode = argc > 1 ? atol(argv[1]) : 0;
    maxlen = argc > 2 ? atol(argv[2]) : 10;

    flint_printf("ns per product (mode %wd); columns: mpn classical toom_karatsuba toom_scalar dispatch  | winner\n", mode);

    for (len = 2; len <= maxlen; len = (len < 8) ? len + 1 : len + len / 2)
    {
        for (i = 0; i < nb; i++)
        {
            fmpz * a, * b, * r;
            double t[5];
            const char * names[5] = {"mpn", "classical", "toomk", "toom", "dispatch"};
            int best = 0;

            bits = bl[i];
            a = _fmpz_vec_init(len);
            b = _fmpz_vec_init(len);
            r = _fmpz_vec_init(2 * len);

            for (j = 0; j < len; j++)
            {
                fmpz_randbits(a + j, state, bits);
                fmpz_randbits(b + j, state, bits);
                if (mode == 1) { fmpz_abs(a + j, a + j); fmpz_abs(b + j, b + j); }
            }

            if (mode == 2) { _fmpz_vec_clear(b, len); b = a; }

            t[0] = tm(_fmpz_poly_mul_mpn, r, a, len, b, len);
            t[1] = tm(_fmpz_poly_mul_classical, r, a, len, b, len);
            t[2] = tm(_fmpz_poly_mul_toom_karatsuba, r, a, len, b, len);
            t[3] = (2 * len - 1 <= 20) ? tm(toom, r, a, len, b, len) : -1;
            t[4] = tm(_fmpz_poly_mul, r, a, len, b, len);

            for (j = 1; j < 4; j++)
                if (t[j] >= 0 && t[j] < t[best])
                    best = j;

            flint_printf("len %2wd bits %5wd: %8.0f %8.0f %8.0f %8.0f %8.0f  | %s (mpn/best %.2f)\n",
                len, bits, t[0], t[1], t[2], t[3], t[4], names[best], t[0] / t[best]);

            _fmpz_vec_clear(a, len);
            if (mode != 2) _fmpz_vec_clear(b, len);
            _fmpz_vec_clear(r, 2 * len);
        }
    }

    flint_rand_clear(state);
    return 0;
}
