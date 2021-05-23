/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz.h"

void
fmpz_xgcd_minimal(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
{
    mpz_t dd, aa, bb;

    if (fmpz_is_zero(g))
    {
        fmpz_abs(d, f);
        fmpz_set_si(a, fmpz_sgn(f));
        fmpz_zero(b);
        return;
    }
    else if (fmpz_is_zero(f))
    {
        fmpz_abs(d, g);
        fmpz_set_si(b, fmpz_sgn(g));
        fmpz_zero(a);
        return;
    }

	mpz_init(dd);
	mpz_init(aa);
	mpz_init(bb);
	
    if (!COEFF_IS_MPZ(*f) && !COEFF_IS_MPZ(*g))  /* both are small */
    {
        mpz_t ff, gg;
        fmpz_t absf, absg;

        fmpz_init(absf);
        fmpz_init(absg);
        fmpz_abs(absf, f);
        fmpz_abs(absg, g);

        ff->_mp_alloc = 1;
        gg->_mp_alloc = 1;
        ff->_mp_size  = fmpz_sgn(f);
        gg->_mp_size  = fmpz_sgn(g);
        ff->_mp_d = (mp_limb_t *) absf;
        gg->_mp_d = (mp_limb_t *) absg;

        mpz_gcdext(dd, aa, bb, ff, gg);

        fmpz_clear(absf);
        fmpz_clear(absg);
    }
    else if (!COEFF_IS_MPZ(*f))  /* only f is small */
    {
        mpz_t ff;
        fmpz_t absf;

        fmpz_init(absf);
        fmpz_abs(absf, f);

        ff->_mp_alloc = 1;
        ff->_mp_size = fmpz_sgn(f);
        ff->_mp_d = (mp_limb_t *) absf;

        mpz_gcdext(dd, aa, bb, ff, COEFF_TO_PTR(*g));

        fmpz_clear(absf);
    }
    else if (!COEFF_IS_MPZ(*g))  /* only g is small */
    {
        mpz_t gg;
        fmpz_t absg;

        fmpz_init(absg);
        fmpz_abs(absg, g);

        gg->_mp_alloc = 1;
        gg->_mp_size = fmpz_sgn(g);
        gg->_mp_d = (mp_limb_t *) absg;

        mpz_gcdext(dd, aa, bb, COEFF_TO_PTR(*f), gg);

        fmpz_clear(absg);
    }
    else /* both are big */
    {
        mpz_gcdext(dd, aa, bb, COEFF_TO_PTR(*f), COEFF_TO_PTR(*g));
    }

	_fmpz_promote_val(d);
    _fmpz_promote_val(a);
    _fmpz_promote_val(b);

	mpz_swap(COEFF_TO_PTR(*d), dd);
	mpz_swap(COEFF_TO_PTR(*a), aa);
	mpz_swap(COEFF_TO_PTR(*b), bb);

	mpz_clear(dd);
	mpz_clear(aa);
	mpz_clear(bb);

    _fmpz_demote_val(d);
    _fmpz_demote_val(a);
    _fmpz_demote_val(b);
}
