/*
    Copyright (C) 2016 Vincent Delecroix
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

void
_qqbar_randtest(qqbar_t res, flint_rand_t state, slong deg, slong bits, int real)
{
    fmpz_poly_t pol;
    slong prec, i, rdeg, r1, r2;
    acb_ptr roots;

    deg = FLINT_MAX(deg, 1);
    bits = FLINT_MAX(bits, 1);

    if ((deg == 1 || n_randint(state, 4) == 0) && real != 2)
    {
        fmpq_t t;
        fmpq_init(t);
        do {
            fmpq_randtest(t, state, bits);
        } while (fmpz_bits(fmpq_numref(t)) > bits ||
                 fmpz_bits(fmpq_denref(t)) > bits);
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
        return;
    }

    fmpz_poly_init(pol);

    do {
        fmpz_poly_randtest_irreducible(pol, state, deg + 1, bits);
        rdeg = fmpz_poly_degree(pol);
        r1 = rdeg;
        r2 = 0;
        if (real)
            fmpz_poly_signature(&r1, &r2, pol);
    }
    while (rdeg < 1 || (real == 1 && r1 < 1) || (real == 2 && r2 < 1));

    if (fmpz_sgn(pol->coeffs + rdeg) < 0)
        fmpz_poly_neg(pol, pol);

    roots = _acb_vec_init(rdeg);
    if (real == 0)
        i = n_randint(state, rdeg);
    else if (real == 1)
        i = n_randint(state, r1);
    else
        i = r1 + n_randint(state, 2 * r2);

    for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
    {
        arb_fmpz_poly_complex_roots(roots, pol, 0, prec);

        if (_qqbar_validate_uniqueness(roots + i, pol, roots + i, 2 * prec))
        {
            fmpz_poly_set(QQBAR_POLY(res), pol);
            acb_set(QQBAR_ENCLOSURE(res), roots + i);
            break;
        }
    }

    _acb_vec_clear(roots, rdeg);
    fmpz_poly_clear(pol);
}

void
qqbar_randtest(qqbar_t res, flint_rand_t state, slong deg, slong bits)
{
    _qqbar_randtest(res, state, deg, bits, 0);
}

void
qqbar_randtest_real(qqbar_t res, flint_rand_t state, slong deg, slong bits)
{
    _qqbar_randtest(res, state, deg, bits, 1);
}

void
qqbar_randtest_nonreal(qqbar_t res, flint_rand_t state, slong deg, slong bits)
{
    if (deg <= 1)
    {
        flint_throw(FLINT_ERROR, "qqbar_randtest_nonreal: must have deg >= 2\n");
    }

    _qqbar_randtest(res, state, deg, bits, 2);
}
