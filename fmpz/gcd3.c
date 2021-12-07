/*
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

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
