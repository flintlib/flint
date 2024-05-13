/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
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
            fmpz_set_ui(f, mpn_gcd_1((nn_srcptr) &u2, (slong) 1, u1));
        }
        else                    /* but h is large */
        {
            mpz_ptr mpzc2 = COEFF_TO_PTR(c2);
            slong size = mpzc2->_mp_size;
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
            mpz_ptr mpzc1;
            slong size;

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
            mpz_ptr mf = _fmpz_promote(f);
            mpz_gcd(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);
        }
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
        mpz_ptr ma = COEFF_TO_PTR(*a);
        fmpz_set_ui(res, mpn_gcd_1(ma->_mp_d, FLINT_ABS(ma->_mp_size), b));
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
                    mpz_ptr mb = COEFF_TO_PTR(*b);
                    c = mpn_gcd_1(mb->_mp_d, FLINT_ABS(mb->_mp_size), c);
                }
            }
        }
        else
        {
            mpz_ptr ma = COEFF_TO_PTR(*a);

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
        mpz_ptr rp, ap, bp, cp, tp;
        slong an, bn, cn, mn;

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
            t->_mp_d = TMP_ALLOC(sizeof(ulong) * cn);
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
