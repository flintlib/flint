/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"

int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p)
{
    slong n;
    nn_ptr d, x, r;
    int success;
    TMP_INIT;

    if (b == a || b == p)
    {
        int ans;
        fmpz_t t;

        fmpz_init(t);
        ans = fmpz_sqrtmod(t, a, p);
        fmpz_swap(b, t);
        fmpz_clear(t);
        return ans;
    }

    if (fmpz_sgn(p) <= 0)
        return 0;

    fmpz_mod(b, a, p);

    if (fmpz_cmp_ui(b, 1) <= 0)
        return 1;

    if (fmpz_is_even(p))
        return 0;

    /*
        Primality is assumed, not checked, but a modulus that is a perfect
        square has no quadratic nonresidue at all and would send the search
        for one on a long walk; it is cheap to rule out here.
    */
    if (fmpz_is_square(p))
    {
        fmpz_zero(b);
        return 0;
    }

    n = fmpz_size(p);

    TMP_START;

    d = TMP_ALLOC((3 * n) * sizeof(ulong));
    x = d + n;
    r = x + n;

    fmpz_get_ui_array(d, n, p);
    fmpz_get_ui_array(x, n, b);

    /* a failure of the algorithm is reported as no root, as documented */
    success = (flint_mpn_sqrtmod(r, x, d, n) == 1);

    if (success)
        fmpz_set_ui_array(b, r, n);
    else
        fmpz_zero(b);

    TMP_END;

    return success;
}
