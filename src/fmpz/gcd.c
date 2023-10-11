/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        ulong u1;

        if (c1 == 0)
        {
            fmpz_abs(f, h);
            return;
        }

        u1 = FLINT_ABS(c1);
        if (!COEFF_IS_MPZ(c2))  /* and h is also small */
        {
            ulong u2;

            if (c2 == 0)
            {
                fmpz_abs(f, g);
                return;
            }

            u2 = FLINT_ABS(c2);
            fmpz_set_ui(f, mpn_gcd_1((mp_srcptr) &u2, (mp_size_t) 1, u1));
        }
        else                    /* but h is large */
        {
            __mpz_struct * mpzc2 = COEFF_TO_PTR(c2);
            mp_size_t size = mpzc2->_mp_size;
            /* The sign is stored in the size of an mpz, and gcd_1 only takes
             * positive integers. */
            fmpz_set_ui(f, mpn_gcd_1(mpzc2->_mp_d, FLINT_ABS(size), u1));
        }
    }
    else                        /* g is large */
    {
        if (!COEFF_IS_MPZ(c2))  /* but h is small */
        {
            ulong u2;
            __mpz_struct * mpzc1;
            mp_size_t size;

            if (c2 == 0)
            {
                fmpz_abs(f, g);
                return;
            }

            u2 = FLINT_ABS(c2);
            mpzc1 = COEFF_TO_PTR(c1);
            size = mpzc1->_mp_size;
            fmpz_set_ui(f, mpn_gcd_1(mpzc1->_mp_d, FLINT_ABS(size), u2));
        }
        else
        {
            /* TODO: Change to mpn_gcd in order to save some calculations that
             * have already been already made. */
            __mpz_struct * mf = _fmpz_promote(f);
            mpz_gcd(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);
        }
    }
}

/* Assumes that c fits in fmpz */
static void
_fmpz_gcd3_small(fmpz_t res, const fmpz_t a, const fmpz_t b, ulong c)
{
    if (c <= 1)
    {
        if (c == 1)
            fmpz_one(res);
        else
            fmpz_gcd(res, a, b);
    }
    else
    {
        if (!COEFF_IS_MPZ(*a))
        {
            if (*a != 0)
                c = mpn_gcd_1(&c, 1, FLINT_ABS(*a));

            if (c != 1)
            {
                if (!COEFF_IS_MPZ(*b))
                {
                    if (*b != 0)
                        c = mpn_gcd_1(&c, 1, FLINT_ABS(*b));
                }
                else
                {
                    __mpz_struct * mb = COEFF_TO_PTR(*b);
                    c = mpn_gcd_1(mb->_mp_d, FLINT_ABS(mb->_mp_size), c);
                }
            }
        }
        else
        {
            __mpz_struct * ma = COEFF_TO_PTR(*a);

            if (!COEFF_IS_MPZ(*b))
            {
                if (*b != 0)
                    c = mpn_gcd_1(&c, 1, FLINT_ABS(*b));

                if (c != 1)
                    c = mpn_gcd_1(ma->_mp_d, FLINT_ABS(ma->_mp_size), c);
            }
            else
            {
                c = mpn_gcd_1(ma->_mp_d, FLINT_ABS(ma->_mp_size), c);

                if (c != 1)
                {
                    ma = COEFF_TO_PTR(*b);
                    c = mpn_gcd_1(ma->_mp_d, FLINT_ABS(ma->_mp_size), c);
                }
            }
        }

        if (COEFF_IS_MPZ(*res))
            _fmpz_demote(res);
        *res = c;
    }
}

void
fmpz_gcd3(fmpz_t res, const fmpz_t a, const fmpz_t b, const fmpz_t c)
{
    if (!COEFF_IS_MPZ(*a))
    {
        _fmpz_gcd3_small(res, b, c, FLINT_ABS(*a));
    }
    else if (!COEFF_IS_MPZ(*b))
    {
        _fmpz_gcd3_small(res, a, c, FLINT_ABS(*b));
    }
    else if (!COEFF_IS_MPZ(*c))
    {
        _fmpz_gcd3_small(res, a, b, FLINT_ABS(*c));
    }
    else
    {
        /* Three-way mpz_gcd. */
        __mpz_struct *rp, *ap, *bp, *cp, *tp;
        mp_size_t an, bn, cn, mn;

        /* If res is small, it cannot be aliased with a, b, c, so promoting is fine. */
        rp = _fmpz_promote(res);

        ap = COEFF_TO_PTR(*a);
        bp = COEFF_TO_PTR(*b);
        cp = COEFF_TO_PTR(*c);

        an = FLINT_ABS(ap->_mp_size);
        bn = FLINT_ABS(bp->_mp_size);
        cn = FLINT_ABS(cp->_mp_size);

        /* Select c to be the largest operand; we do the smaller gcd first. */
        mn = FLINT_MAX(FLINT_MAX(an, bn), cn);
        tp = cp;
        if (mn != cn)
        {
            if (mn == an)
            {
                cp = ap;
                ap = tp;
            }
            else
            {
                cp = bp;
                bp = tp;
            }

            cn = mn;
        }

        /* Handle aliasing */
        if (rp == cp)
        {
            mpz_t t;
            TMP_INIT;
            TMP_START;

            /* It would be more efficient to allocate temporary space for
               gcd(a, b), but we can't be sure that mpz_gcd never attempts
               to reallocate the output. */
            t->_mp_d = TMP_ALLOC(sizeof(mp_limb_t) * cn);
            t->_mp_size = t->_mp_alloc = cn;
            flint_mpn_copyi(t->_mp_d, cp->_mp_d, cn);

            mpz_gcd(rp, ap, bp);
            if (mpz_cmpabs_ui(rp, 1) != 0)
                mpz_gcd(rp, rp, t);

            TMP_END;
        }
        else
        {
            mpz_gcd(rp, ap, bp);
            if (mpz_cmpabs_ui(rp, 1) != 0)
                mpz_gcd(rp, rp, cp);
        }

        /* The result may be small */
        _fmpz_demote_val(res);
    }
}

void
fmpz_gcd_ui(fmpz_t res, const fmpz_t a, ulong b)
{
    if (b == 0)
        fmpz_abs(res, a);
    else if (!COEFF_IS_MPZ(*a))
    {
        if (*a != 0)
        {
            _fmpz_demote(res);
            *res = mpn_gcd_1(&b, 1, FLINT_ABS(*a));
        }
        else
            fmpz_set_ui(res, b);
    }
    else
    {
        __mpz_struct * ma = COEFF_TO_PTR(*a);
        fmpz_set_ui(res, mpn_gcd_1(ma->_mp_d, FLINT_ABS(ma->_mp_size), b));
    }
}

void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g)
{
    FLINT_ASSERT(fmpz_cmp(f, g) < 0);

    if (fmpz_is_zero(f))
    {
        fmpz_set(d, g);
        fmpz_set_ui(a, 0);
	return;
    }

    if (!COEFF_IS_MPZ(*g))  /* g is small, hence f is small */
    {
        fmpz ff, gg;
        ff = *f;
        gg = *g;

        _fmpz_demote(d);
        _fmpz_demote(a);

        *d = n_gcdinv((mp_limb_t *) a, ff, gg);
    }
    else  /* g is large */
    {
        mpz_t atemp, dtemp;

	mpz_init(atemp);
	mpz_init(dtemp);

	_fmpz_promote_val(d);
        _fmpz_promote_val(a);

        if (!COEFF_IS_MPZ(*f))  /* f is small */
        {
            mpz_t fptr;

            fptr->_mp_alloc = 1;
            fptr->_mp_size  = 1;
            fptr->_mp_d     = (mp_limb_t *) f;

            mpz_gcdext(dtemp, atemp, NULL,
                       fptr, COEFF_TO_PTR(*g));
        }
        else  /* f is large */
        {
            mpz_gcdext(dtemp, atemp, NULL,
                       COEFF_TO_PTR(*f), COEFF_TO_PTR(*g));
        }

	if (mpz_cmp_ui(atemp, 0) < 0)
           mpz_add(atemp, atemp, COEFF_TO_PTR(*g));

	mpz_swap(COEFF_TO_PTR(*d), dtemp);
	mpz_swap(COEFF_TO_PTR(*a), atemp);

	mpz_clear(atemp);
	mpz_clear(dtemp);

        _fmpz_demote_val(d);
        _fmpz_demote_val(a);
    }
}
