/*
    Copyright (C) 2012 Thomas M. DuBuisson
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"

void fmpz_and(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
        fmpz c1,c2;
        c1 = *g;
        c2 = *h;
        if(!COEFF_IS_MPZ(c1))
        {
            if(!COEFF_IS_MPZ(c2)) /* both inputs are small */
            {
                fmpz_set_si(f, c1 & c2);
            } else /* g is small, h is large */
            {
                mpz_t tmp;
                mpz_ptr mpz3 = _fmpz_promote(f);
                flint_mpz_init_set_si(tmp, c1);
                mpz_and(mpz3, COEFF_TO_PTR(c2), tmp);
                _fmpz_demote_val(f);
                mpz_clear(tmp);
            }
        } else
        {
            if(!COEFF_IS_MPZ(c2)) /* g is large, h is small */
            {
                mpz_t tmp;
                mpz_ptr mpz3 = _fmpz_promote(f);
                flint_mpz_init_set_si(tmp, c2);
                mpz_and(mpz3, COEFF_TO_PTR(c1), tmp);
                _fmpz_demote_val(f);
                mpz_clear(tmp);
            } else /* g and h are large */
            {
                mpz_ptr mpz3 = _fmpz_promote(f);
                mpz_ptr mpz1 = COEFF_TO_PTR(c1);
                mpz_ptr mpz2 = COEFF_TO_PTR(c2);
                mpz_and(mpz3, mpz1, mpz2);
                _fmpz_demote_val(f);
            }
        }
}

void fmpz_complement(fmpz_t r, const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f)) /* f is small */
    {
	slong res = ~(*f);
	fmpz_set_si(r, res);
    } else /* f is big */
    {
        if(r != f) { /* not aliased */
            mpz_ptr ptr, ptr2;
            ptr = _fmpz_promote(r);
            ptr2 = COEFF_TO_PTR(*f);
            mpz_com(ptr, ptr2);
            _fmpz_demote_val(r);
        } else { /* alaised */
            fmpz_t tmp;
            mpz_ptr ptr, ptr2;
            fmpz_init(tmp);
            ptr = _fmpz_promote(tmp);
            ptr2 = COEFF_TO_PTR(*f);
            mpz_com(ptr, ptr2);
            _fmpz_demote_val(tmp);
            fmpz_set(r,tmp);
            fmpz_clear(tmp);
        }
    }
}

void fmpz_or(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
        fmpz c1,c2;
        c1 = *g;
        c2 = *h;
        if(!COEFF_IS_MPZ(c1))
        {
            if(!COEFF_IS_MPZ(c2)) /* both inputs are small */
            {
                fmpz_set_si(f, c1 | c2);
            } else /* g is small, h is large */
            {
                mpz_t tmp;
                mpz_ptr mpz3 = _fmpz_promote(f);
                flint_mpz_init_set_si(tmp,c1);
                mpz_ior(mpz3, COEFF_TO_PTR(c2), tmp);
                _fmpz_demote_val(f);
                mpz_clear(tmp);
            }
        } else
        {
            if(!COEFF_IS_MPZ(c2)) /* g is large, h is small */
            {
                mpz_t tmp;
                mpz_ptr mpz3 = _fmpz_promote(f);
                flint_mpz_init_set_si(tmp,c2);
                mpz_ior(mpz3, COEFF_TO_PTR(c1), tmp);
                _fmpz_demote_val(f);
                mpz_clear(tmp);
            } else /* g and h are large */
            {
                mpz_ptr mpz3 = _fmpz_promote(f);
                mpz_ptr mpz1 = COEFF_TO_PTR(c1);
                mpz_ptr mpz2 = COEFF_TO_PTR(c2);
                mpz_ior(mpz3, mpz1, mpz2);
                _fmpz_demote_val(f);
            }
        }
}

void fmpz_xor(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
        fmpz c1,c2;
        c1 = *g;
        c2 = *h;
        if(!COEFF_IS_MPZ(c1))
        {
            if(!COEFF_IS_MPZ(c2)) /* both inputs are small */
            {
                fmpz_set_si(f, c1 ^ c2);
            } else /* g is small, h is large */
            {
                mpz_t tmp;
                mpz_ptr mpz3 = _fmpz_promote(f);
                flint_mpz_init_set_si(tmp, c1);
                mpz_xor(mpz3, COEFF_TO_PTR(c2), tmp);
                _fmpz_demote_val(f);
                mpz_clear(tmp);
            }
        } else
        {
            if(!COEFF_IS_MPZ(c2)) /* g is large, h is small */
            {
                mpz_t tmp;
                mpz_ptr mpz3 = _fmpz_promote(f);
                flint_mpz_init_set_si(tmp, c2);
                mpz_xor(mpz3, COEFF_TO_PTR(c1), tmp);
                _fmpz_demote_val(f);
                mpz_clear(tmp);
            } else /* g and h are large */
            {
                mpz_ptr mpz3 = _fmpz_promote(f);
                mpz_ptr mpz1 = COEFF_TO_PTR(c1);
                mpz_ptr mpz2 = COEFF_TO_PTR(c2);
                mpz_xor(mpz3, mpz1, mpz2);
                _fmpz_demote_val(f);
            }
        }
}

void fmpz_clrbit(fmpz_t f, ulong i)
{
    if (!COEFF_IS_MPZ(*f))
    {
        if (i < SMALL_FMPZ_BITCOUNT_MAX)
        {
            *f &= ~(WORD(1) << i);
        }
        /* i >= FLINT_BITS  --> nop */
    }
    else
    {
        mpz_ptr ptr = COEFF_TO_PTR(*f);

        mpz_clrbit(ptr, i);
        _fmpz_demote_val(f);
    }
}

void fmpz_combit(fmpz_t f, ulong i)
{
    if (!COEFF_IS_MPZ(*f))
    {
        if (i < SMALL_FMPZ_BITCOUNT_MAX)
        {
            *f ^= (WORD(1) << i);
        }
        else  /* i >= FLINT_BITS */
        {
            mpz_ptr ptr = _fmpz_promote_val(f);
            mpz_combit(ptr, i);
            _fmpz_demote_val(f);
        }
    }
    else
    {
        mpz_ptr ptr = COEFF_TO_PTR(*f);
        mpz_combit(ptr, i);
        _fmpz_demote_val(f);
    }
}

flint_bitcnt_t fmpz_popcnt(const fmpz_t c)
{
    ulong d;
    fmpz c1 = *c;

    if (!COEFF_IS_MPZ(c1))
    {
        if (c1 < 0)
            return 0;

        d = c1;
        return mpn_popcount(&d, 1);
    }
    else
    {
        mpz_ptr t = COEFF_TO_PTR(c1);

        if (mpz_sgn(t) < 0)
            return 0;
        else
            return mpz_popcount(t);
    }
}
