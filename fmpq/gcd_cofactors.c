/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

void _fmpq_gcd_cofactors(
    fmpz_t ng_, fmpz_t dg_,
    fmpz_t abar,
    fmpz_t bbar,
    const fmpz_t na, const fmpz_t da,
    const fmpz_t nb, const fmpz_t db)
{
    fmpz_t ng, dg, nabar, dabar, nbbar, dbbar;
#if FLINT_WANT_ASSERT
    fmpq_t cqt_g;
    int input_is_canonical = _fmpq_is_canonical(na, da) &&
                             _fmpq_is_canonical(nb, db);
#endif

    fmpz_init(ng);

    fmpz_gcd(ng, na, nb);
    if (fmpz_is_zero(ng))
    {
        fmpz_zero(ng_);
        fmpz_one(dg_);
        fmpz_zero(abar);
        fmpz_zero(bbar);
        fmpz_clear(ng);
        return;
    }

#if FLINT_WANT_ASSERT
    fmpq_init(cqt_g);
    _fmpq_gcd(fmpq_numref(cqt_g), fmpq_denref(cqt_g), na, da, nb, db);
#endif

    fmpz_init(dg);
    fmpz_init(nabar);
    fmpz_init(dabar);
    fmpz_init(nbbar);
    fmpz_init(dbbar);

    fmpz_divexact(nabar, na, ng);
    fmpz_divexact(nbbar, nb, ng);

    fmpz_gcd(dg, da, db);
    fmpz_divexact(dabar, da, dg);
    fmpz_divexact(dbbar, db, dg);

    fmpz_mul(abar, nabar, dbbar);
    fmpz_mul(bbar, dabar, nbbar);
    fmpz_mul(dg_, da, dbbar);
    fmpz_swap(ng_, ng);

#if FLINT_WANT_ASSERT
    if (input_is_canonical)
    {
        FLINT_ASSERT(fmpz_equal(fmpq_numref(cqt_g), ng_));
        FLINT_ASSERT(fmpz_equal(fmpq_denref(cqt_g), dg_));
        fmpz_gcd(ng, abar, bbar);
        FLINT_ASSERT(fmpz_is_one(ng));
    }
    fmpq_clear(cqt_g);
#endif

    fmpz_clear(ng);
    fmpz_clear(dg);
    fmpz_clear(nabar);
    fmpz_clear(dabar);
    fmpz_clear(nbbar);
    fmpz_clear(dbbar);
}
