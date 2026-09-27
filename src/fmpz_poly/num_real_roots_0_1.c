/*
   Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly/impl.h"

slong fmpz_poly_num_real_roots_0_1(const fmpz_poly_t pol)
{
    const fmpz * p = pol->coeffs;
    slong len = pol->length;
    slong i, bits, res, n_exact, n_interval;

    if (len == 0)
        flint_throw(FLINT_ERROR, "zero polynomial in %s\n", __func__);

    /* Remove roots at 0, which are not counted. */
    for (i = 0; i < len && fmpz_is_zero(p + i); i++)
        ;
    p += i;
    len -= i;

    if (len <= 1)
        return 0;

    /* Cheap Descartes test on (0, 1); this is also the first step of VCA,
       and succeeds immediately for most generic polynomials. */
    res = _fmpz_poly_descartes_bound_0_1(p, len, 1);
    if (res <= 1)
        return res;

    /* The primitive Sturm sequence wins when the remainders stay small
       (e.g. classical orthogonal polynomials, sparse polynomials); abort
       as soon as they grow and fall back on VCA. */
    bits = FLINT_ABS(_fmpz_vec_max_bits(p, len));
    res = _fmpz_poly_num_real_roots_sturm_bounded(p, len, 1,
        _fmpz_poly_num_real_roots_sturm_bound(len, bits));
    if (res >= 0)
        return res;

    n_exact = n_interval = 0;
    _fmpz_poly_isolate_real_roots_0_1_vca(NULL, &n_exact, NULL, NULL,
        &n_interval, p, len);
    return n_exact + n_interval;
}
