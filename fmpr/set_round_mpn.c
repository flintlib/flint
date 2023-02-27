/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

/* like mpn_scan0, but takes an upper size */
static __inline__ flint_bitcnt_t
mpn_scan0b(mp_srcptr up, mp_size_t size, flint_bitcnt_t from_bit)
{
    mp_limb_t t;
    slong i, c;

    i = from_bit / GMP_NUMB_BITS;
    c = from_bit % FLINT_BITS;
    t = ((~up[i]) >> c) << c;

    while (t == 0)
    {
        i++;
        if (i == size)
            return size * FLINT_BITS;
        else
            t = ~up[i];
    }

    count_trailing_zeros(c, t);
    return (i * FLINT_BITS) + c;
}

slong
_fmpr_set_round_mpn(slong * shift, fmpz_t man, mp_srcptr x, mp_size_t xn, int negative, slong prec, fmpr_rnd_t rnd)
{
    slong bc, val, val_bits, val_limbs, ret;
    int increment;

    /* compute the total bit length of x */
    count_leading_zeros(bc, x[xn - 1]);
    bc = FLINT_BITS - bc;
    bc += (xn - 1) * FLINT_BITS;

    /* already odd */
    if (x[0] & 1)
    {
        /* quick exit */
        if (bc <= prec)
        {
            if (bc <= FLINT_BITS - 2)
            {
                mp_limb_t t = x[0];
                _fmpz_demote(man);
                *man = negative ? -t : t;
            }
            else
                fmpz_set_mpn_large(man, x, xn, negative);

            *shift = 0;
            return FMPR_RESULT_EXACT;
        }
        else
        {
            val_limbs = val_bits = val = 0;
        }
    }
    else
    {
        /* trailing zero bits: val = val_limbs * FLINT_BITS + val_bits */
        val_limbs = 0;
        while (x[val_limbs] == 0)
            val_limbs++;
        count_trailing_zeros(val_bits, x[val_limbs]);
        val = val_bits + (val_limbs * FLINT_BITS);
    }

    /* no rounding necessary; just copy or shift to destination */
    if (bc - val <= prec)
    {
        ret = FMPR_RESULT_EXACT;
        increment = 0;
    }
    else
    {
        /* truncation */
        if (!rounds_up(rnd, negative))
        {
            val = mpn_scan1(x, bc - prec);
            increment = 0;
        }
        /* round to next higher odd mantissa */
        else
        {
            val = mpn_scan0b(x, xn, bc - prec);

            /* can overflow to next power of 2 */
            if (val == bc)
            {
                fmpz_set_si(man, negative ? -1 : 1);
                *shift = bc;
                return prec - 1;
            }

            /* otherwise, we are cutting off at a zero bit, and
               incrementing at that position will not cause carry
               propagation below */
            increment = 1;
        }

        val_limbs = val / FLINT_BITS;
        val_bits = val % FLINT_BITS;
        ret = prec - (bc - val);
    }

    /* the output mantissa is a small fmpz */
    if (bc - val <= FLINT_BITS - 2)
    {
        mp_limb_t h;

        if (val_limbs + 1 == xn || val_bits == 0)
            h = x[val_limbs] >> val_bits;
        else
            h = (x[val_limbs] >> val_bits) | (x[val_limbs + 1] << (FLINT_BITS - val_bits));

        h += increment;
        _fmpz_demote(man);
        *man = negative ? -h : h;
    }
    /* the output mantissa is an mpz */
    else
    {
        mp_ptr dest;
        slong res_limbs, res_alloc;
        __mpz_struct * zptr = _fmpz_promote(man);

        res_limbs = ((bc - val) + FLINT_BITS - 1) / FLINT_BITS;
        /* todo: only allocate the required size, not the size before shifting */
        res_alloc = xn - val_limbs;

        if (zptr->_mp_alloc < res_alloc)
            mpz_realloc2(zptr, res_alloc * FLINT_BITS);

        dest = zptr->_mp_d;

        /* right shift by val */
        if (val_bits == 0)
            flint_mpn_copyi(dest, x + val_limbs, res_limbs);
        else
            mpn_rshift(dest, x + val_limbs, xn - val_limbs, val_bits);

        dest[0] += increment;
        zptr->_mp_size = negative ? -res_limbs : res_limbs;
    }

    *shift = val;
    return ret;
}

