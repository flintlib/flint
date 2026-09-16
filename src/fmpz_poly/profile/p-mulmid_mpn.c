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
#include "mpn_extras.h"
#include "fmpz_poly/impl.h"
#include "ulong_extras.h"
#include "profiler.h"

/* Usage: p-mulmid_mpn [bits] [len1] [len2] [mode]
   mode: 0 = random signed, 1 = nonnegative, 2 = squaring (signed), 3 = squaring (nonnegative)
   Prints ns per product. */

static double
time_it(void (*f)(fmpz *, const fmpz *, slong, const fmpz *, slong, slong, slong),
        fmpz * res, const fmpz * a, slong len1, const fmpz * b, slong len2)
{
    double tcpu, twall;
    TIMEIT_START
    f(res, a, len1, b, len2, 0, len1 + len2 - 1);
    TIMEIT_STOP_VALUES(tcpu, twall);
    return tcpu * 1e9;
}

static void
kar_sgn(fmpz * res, const fmpz * a, slong len1, const fmpz * b, slong len2, slong nlo, slong nhi)
{
    _flint_mpn_poly_mulmid_force_method = 1;
    _fmpz_poly_mulmid_mpn(res, a, len1, b, len2, nlo, nhi);
    _flint_mpn_poly_mulmid_force_method = -1;
}

static void
kar_bias(fmpz * res, const fmpz * a, slong len1, const fmpz * b, slong len2, slong nlo, slong nhi)
{
    _flint_mpn_poly_mulmid_force_method = 2;
    _fmpz_poly_mulmid_mpn(res, a, len1, b, len2, nlo, nhi);
    _flint_mpn_poly_mulmid_force_method = -1;
}

int main(int argc, char * argv[])
{
    slong bits, len1, len2, mode, i, cutoff;
    fmpz * a, * b, * r1, * r2;
    FLINT_TEST_INIT(state);

    bits = argc > 1 ? atol(argv[1]) : 100;
    len1 = argc > 2 ? atol(argv[2]) : 32;
    len2 = argc > 3 ? atol(argv[3]) : len1;
    mode = argc > 4 ? atol(argv[4]) : 0;

    a = _fmpz_vec_init(len1);
    b = _fmpz_vec_init(len2);
    r1 = _fmpz_vec_init(len1 + len2 - 1);
    r2 = _fmpz_vec_init(len1 + len2 - 1);

    for (i = 0; i < len1; i++) { fmpz_randbits(a + i, state, bits); if (mode == 1 || mode == 3) fmpz_abs(a + i, a + i); }
    for (i = 0; i < len2; i++) { fmpz_randbits(b + i, state, bits); if (mode == 1 || mode == 3) fmpz_abs(b + i, b + i); }
    if (mode >= 2) { _fmpz_vec_clear(b, len2); len2 = len1; b = a; }

    _fmpz_poly_mulmid_classical(r1, a, len1, b, len2, 0, len1 + len2 - 1);
    _fmpz_poly_mulmid_mpn(r2, a, len1, b, len2, 0, len1 + len2 - 1);
    if (!_fmpz_vec_equal(r1, r2, len1 + len2 - 1)) { flint_printf("FAIL\n"); flint_abort(); }

    flint_printf("bits=%wd len1=%wd len2=%wd mode=%wd :", bits, len1, len2, mode);
    flint_printf("  auto %.0f", time_it(_fmpz_poly_mulmid, r1, a, len1, b, len2));
    flint_printf("  mpn %.0f", time_it(_fmpz_poly_mulmid_mpn, r1, a, len1, b, len2));
    _flint_mpn_poly_mulmid_force_cutoff = WORD_MAX;
    flint_printf("  mpn_cl %.0f", time_it(_fmpz_poly_mulmid_mpn, r1, a, len1, b, len2));
    _flint_mpn_poly_mulmid_force_cutoff = -1;
    if (mode == 0 || mode == 2)
    {
        flint_printf("  kar_sgn %.0f", time_it(kar_sgn, r1, a, len1, b, len2));
        flint_printf("  kar_bias %.0f", time_it(kar_bias, r1, a, len1, b, len2));
    }
    if (argc > 5)
    {
        flint_printf("  cutoffs:");
        for (cutoff = 4; cutoff <= 64; cutoff *= 2)
        {
            _flint_mpn_poly_mulmid_force_cutoff = cutoff;
            flint_printf(" [%wd] %.0f", cutoff, time_it(_fmpz_poly_mulmid_mpn, r1, a, len1, b, len2));
        }
        _flint_mpn_poly_mulmid_force_cutoff = -1;
    }
    flint_printf("\n");

    _fmpz_vec_clear(a, len1);
    if (mode < 2) _fmpz_vec_clear(b, len2);
    _fmpz_vec_clear(r1, len1 + len2 - 1);
    _fmpz_vec_clear(r2, len1 + len2 - 1);
    flint_rand_clear(state);
    return 0;
}
