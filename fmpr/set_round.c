/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

/* like mpn_scan0b, but takes an upper size */
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
_fmpr_set_round(fmpz_t rman, fmpz_t rexp,
    const fmpz_t man, const fmpz_t exp, slong prec, fmpr_rnd_t rnd)
{
    if (!COEFF_IS_MPZ(*man))
    {
        slong lead, trail, bc, v, w, shift, ret;
        int negative;

        v = *man;

        count_trailing_zeros(trail, v);
        v >>= trail;
        shift = trail;
        ret = FMPR_RESULT_EXACT;

        /* may need to round */
        if (prec < FLINT_BITS - 2 - trail)
        {
            if (v < 0)
            {
                negative = 1;
                w = -v;
            }
            else
            {
                negative = 0;
                w = v;
            }

            count_leading_zeros(lead, w);
            bc = FLINT_BITS - lead;

            /* round */
            if (prec < bc)
            {
                w = (w >> (bc - prec)) + rounds_up(rnd, negative);
                shift += bc - prec;
                count_trailing_zeros(trail, w);
                w >>= trail;
                shift += trail;
                v = negative ? -w : w;
                ret = trail;

                /* special case: if w overflowed to the next power of two,
                   the error bound must be multiplied by 2 */
                ret -= (trail == prec);
            }
        }

        /* the input is small, so the output must be small too */
        _fmpz_demote(rman);
        *rman = v;
        fmpz_add_ui_inline(rexp, exp, shift);
        return ret;
    }
    else
    {
        slong size, bc, val, val_bits, val_limbs, ret;
        int negative, increment;
        mp_ptr d;
        __mpz_struct * z = COEFF_TO_PTR(*man);

        size = z->_mp_size;
        negative = size < 0;
        size = FLINT_ABS(size);
        d = z->_mp_d;

        /* bit size */
        count_leading_zeros(bc, d[size - 1]);
        bc = FLINT_BITS - bc;
        bc += (size - 1) * FLINT_BITS;

        /* quick exit */
        if (bc <= prec && (d[0] & 1))
        {
            if (rman != man) fmpz_set(rman, man);
            if (rexp != exp) fmpz_set(rexp, exp);
            return FMPR_RESULT_EXACT;
        }

        /* trailing zeros */
        val_limbs = 0;
        while (d[val_limbs] == 0)
            val_limbs++;
        count_trailing_zeros(val_bits, d[val_limbs]);
        val = val_bits + (val_limbs * FLINT_BITS);

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
                val = mpn_scan1(d, bc - prec);
                increment = 0;
            }
            /* round to next higher odd mantissa */
            else
            {
                val = mpn_scan0b(d, size, bc - prec);

                /* can overflow to next power of 2 */
                if (val == bc)
                {
                    fmpz_set_si(rman, negative ? -1 : 1);
                    fmpz_add_ui_inline(rexp, exp, bc);
                    return prec - 1;
                }

                /* otherwise, incrementing will not cause overflow below */
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
            if (val_limbs + 1 == size || val_bits == 0)
                h = d[val_limbs] >> val_bits;
            else
                h = (d[val_limbs] >> val_bits)
                    | (d[val_limbs + 1] << (FLINT_BITS - val_bits));
            h += increment;
            _fmpz_demote(rman);
            *rman = negative ? -h : h;
        }
        /* the output mantissa is an mpz */
        else
        {
            if (rman == man)
            {
                mpz_tdiv_q_2exp(z, z, val);
                if (increment) z->_mp_d[0]++;
            }
            else
            {
                __mpz_struct * w = _fmpz_promote(rman);
                /* must reload pointer, as promoting rman could change it */
                z = COEFF_TO_PTR(*man);
                mpz_tdiv_q_2exp(w, z, val);

                if (increment)
                    w->_mp_d[0]++;
            }
        }

        fmpz_add_ui_inline(rexp, exp, val);
        return ret;
    }
}

