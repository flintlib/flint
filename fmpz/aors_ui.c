/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

static void
_fmpz_add_mpn_1(fmpz_t f, const mp_limb_t * glimbs, mp_size_t gsz, mp_limb_t x);

static void
_fmpz_sub_mpn_1(fmpz_t f, const mp_limb_t * glimbs, mp_size_t gsz, mp_limb_t x);

void
fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    __mpz_struct * mf;
    slong g1 = *g;
    slong f1 = *f;

    if (!COEFF_IS_MPZ(g1))  /* g is small */
    {
        mp_size_t sz = 2;
        if (g1 >= 0)
        {
            {   /* add with jump if carry */
                ulong tmp = g1;
                g1 += x;
                if (((ulong) g1) < tmp)
                    goto carry;
            }
            if (((ulong) g1) <= COEFF_MAX)
            {
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = g1;
                return;
            }
nocarry:    sz = 1; /* No carry, but result does not fit in small fmpz */
carry:      if (COEFF_IS_MPZ(f1))
                mf = COEFF_TO_PTR(f1);
            else
            {
                mf = _fmpz_new_mpz();
                *f = PTR_TO_COEFF(mf);
            }
            mf->_mp_size = sz;
            mf->_mp_d[0] = g1;
            mf->_mp_d[1] = 1; /* Set carry (not used if sz = 1) */
        }
        else /* g < 0 */
        {
            g1 += x;
            if (((slong) x) >= 0 && g1 <= COEFF_MAX)
            {
                /* If x > 0 does not have its top bit set
                 * and COEFF_MIN <= g < 0, we can interpret x + g as a slong.
                 * So if the result in g1 is smaller than COEFF_MAX, it is a
                 * small fmpz. */
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = g1;
                return;
            }
            else
            {
                /* 1) If top bit is set in x, the result is going to be positive
                 *    but will be larger than COEFF_MAX since
                 *
                 *          x + g  >=  (2^63) - (2^62 - 1)  =  2^62 + 1.
                 *
                 *    However, it will be contained in one limb since g < 0.
                 *
                 * 2) If top bit is not set, then result is larger than
                 *    COEFF_MAX, and so it cannot be a small fmpz. However, it
                 *    must fit in one limb since
                 *
                 *          x + g  <=  (2^63 - 1) + (-1)  =  2^63 - 2,
                 *
                 *    which is contained in one limb. */

                goto nocarry;
            }
        }
    }
    else
    {
        __mpz_struct * mg = COEFF_TO_PTR(g1);
        mp_size_t gsz = mg->_mp_size;
        mp_limb_t * glimbs = mg->_mp_d;

        if (gsz > 0)
            _fmpz_add_mpn_1(f, glimbs, gsz, x);
        else
            _fmpz_sub_mpn_1(f, glimbs, gsz, x);
    }
}

/* "Add" two number with same sign. Decide sign from g. */
static void
_fmpz_add_mpn_1(fmpz_t f, const mp_limb_t * glimbs, mp_size_t gsz, mp_limb_t x)
{
    __mpz_struct * mf;
    mp_limb_t * flimbs;
    mp_size_t gabssz = FLINT_ABS(gsz);

    /* Promote f as it is guaranteed to be large */
    if (COEFF_IS_MPZ(*f))
        mf = COEFF_TO_PTR(*f);
    else
    {
        mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
    }
    flimbs = mf->_mp_d;

    if (mf->_mp_alloc < (gabssz + 1)) /* Ensure result fits */
    {
        mp_limb_t * tmp = flimbs;
        flimbs = _mpz_realloc(mf, gabssz + 1);

        /* If f and g are aliased, then we need to change glimbs as well. */
        if (tmp == glimbs)
            glimbs = flimbs;
    }

    /* Use GMP to calculate result */
    flimbs[gabssz] = mpn_add_1(flimbs, glimbs, gabssz, x);

    /* flimbs[gabssz] is the carry from mpn_add_1,
     * and so gabssz + flimbs[gabssz] is valid to determine the size of f */
    mf->_mp_size = gabssz + flimbs[gabssz];
    if (gsz < 0)
    {
        /* g and x has same sign. If g is negative, we negate the result */
        mf->_mp_size = -mf->_mp_size;
    }
}

/* Subtract two limbs (they have different sign) and decide the sign via g. */
static void
_fmpz_sub_mpn_1(fmpz_t f, const mp_limb_t * glimbs, mp_size_t gsz, mp_limb_t x)
{
    __mpz_struct * mf;
    mp_limb_t * flimbs;
    mp_size_t gabssz = FLINT_ABS(gsz);

    /* If size of g is 1, we have a higher probability of the result being
     * small. */
    if (gabssz == 1)
    {
        if (x <= glimbs[0]) /* Result is zero or has the same sign as g */
        {
            x = glimbs[0] - x;
L1:         if (x <= COEFF_MAX) /* Fits in small fmpz */
            {
                if (COEFF_IS_MPZ(*f))
                    _fmpz_clear_mpz(*f);
                *f = (gsz > 0) ? x : -x; /* With consideration of sign */
            }
            else /* Does not fit in small fmpz */
            {
                if (COEFF_IS_MPZ(*f))
                    mf = COEFF_TO_PTR(*f);
                else
                {
                    mf = _fmpz_new_mpz();
                    *f = PTR_TO_COEFF(mf);
                }
                mf->_mp_d[0] = x;
                mf->_mp_size = gsz; /* Sign of f is the same as for g */
            }
        }
        else /* |x| > |g|, which implies f has opposite sign of g */
        {
            /* Set x to the absolute value of |g - x|. By switching sign of
             * gsz, we can reuse the code above. */
            x -= glimbs[0];
            gsz = -gsz;
            goto L1;
        }
        return;
    }

    /* As g has more than one limb, it is a very high probability that result
     * does not fit inside small fmpz. */
    if (COEFF_IS_MPZ(*f))
        mf = COEFF_TO_PTR(*f);
    else
    {
        mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
    }
    flimbs = mf->_mp_d;

    if (gabssz == 2)
    {
        /* Special case. Can result in a small fmpz, but as |g| > |x| the sign
         * cannot change. */
        sub_ddmmss(flimbs[1], flimbs[0], glimbs[1], glimbs[0], 0, x);
        if (flimbs[1] != 0)
        {
            /* Most likely: upper limb not zero, so we just have set the sign
             * of f to g's. */
            mf->_mp_size = gsz;
        }
        else if (flimbs[0] > COEFF_MAX)
        {
            /* Still very likely: Upper limb is zero but lower limb does not
             * fit inside a small fmpz. Sign is the same as for g, but the
             * absolute value of the size is one. */
            mf->_mp_size = (gsz > 0) ? 1 : -1;
        }
        else
        {
            /* Upper limb is zero and lower limb fits inside a small fmpz.
             * Therefore we set f to +/- flimbs[0] and clear the mpz associated
             * to f. */
            slong tmp = flimbs[0]; /* We will clear this mpz, so save first. */
            _fmpz_clear_mpz(*f);
            *f = (gsz > 0) ? tmp : -tmp;
        }
    }
    else
    {
        /* As the absolute value of g's size is larger than 2, the result won't
         * fit inside a small fmpz. */
        if (mf->_mp_alloc < gabssz) /* Ensure result fits */
        {
            /* The allocation size of g is always larger than the absolute value
             * of g. Therefore, if f's allocation size is smaller than g's
             * size, they cannot be aliased. */
            flimbs = _mpz_realloc(mf, gabssz);
        }

        mpn_sub_1(flimbs, glimbs, gabssz, x); /* Subtract via GMP */

        /* If last limb is zero, we have to set f's absolute size to one less
         * than g's. */
        mf->_mp_size = gabssz - (flimbs[gabssz - 1] == 0);
        if (gsz < 0) /* If g is negative, then f is negative as well. */
            mf->_mp_size = -mf->_mp_size;
    }
}

void
fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    __mpz_struct * mf;
    slong g1 = *g;
    slong f1 = *f;

    if (!COEFF_IS_MPZ(g1))  /* g is small */
    {
        mp_size_t sz = -2;
        if (g1 <= 0)
        {
            /* "add" with jump if carry */
            g1 = x - g1; /* g1 = x + |g| */
            if (((ulong) g1) < x)
                goto carry;

            if (((ulong) g1) <= COEFF_MAX)
            {
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = -g1;
                return;
            }
nocarry:    sz = -1; /* No carry, but result is not a small fmpz */
carry:      if (COEFF_IS_MPZ(f1))
                mf = COEFF_TO_PTR(f1);
            else
            {
                mf = _fmpz_new_mpz();
                *f = PTR_TO_COEFF(mf);
            }
            mf->_mp_size = sz;
            mf->_mp_d[0] = g1;
            mf->_mp_d[1] = 1; /* Set carry (not used if sz = -1) */
        }
        else
        {
            g1 = x - g1; /* -(g - x) */
            if (((slong) x) >= 0 && g1 <= COEFF_MAX)
            {
                /* If x > 0 does not have its top bit set
                 * and 0 < g <= COEFF_MAX, we can interpret x - g as a slong.
                 * So if the result in g1 is smaller than COEFF_MAX, it is a
                 * small fmpz. */
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = -g1; /* g - x = -(x - g) */
                return;
            }
            else
            {
                /* 1) If top bit is set in x, the result is going to be negative
                 *    but will be larger than COEFF_MAX since
                 *
                 *          x - g  >=  2^63 - (2^62 - 1)  =  2^62 + 1.
                 *
                 *    However, it will be contained in one limb since g > 0.
                 *
                 * 2) If top bit is not set, then result is smaller than
                 *    COEFF_MIN, and so it cannot be a small fmpz. However, it
                 *    must fit in one limb since
                 *
                 *          x - g  <=  (2^63 - 1) - 1  =  2^63 - 2,
                 *
                 *    which is contained in one limb. */

                goto nocarry;
            }
        }
    }
    else
    {
        __mpz_struct * mg = COEFF_TO_PTR(g1);
        mp_size_t gsz = mg->_mp_size;
        mp_limb_t * glimbs = mg->_mp_d;

        if (gsz > 0)
            _fmpz_sub_mpn_1(f, glimbs, gsz, x);
        else
            _fmpz_add_mpn_1(f, glimbs, gsz, x);
    }
}
