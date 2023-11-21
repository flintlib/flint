/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

static inline void
_arb_get_mag_lower(mag_t z, const arf_t mid, const mag_t rad)
{
    if (arf_is_special(mid) || mag_is_special(rad))
    {
        if (mag_is_zero(rad))
        {
            arf_get_mag_lower(z, mid);
        }
        else if (arf_is_inf(mid) && mag_is_finite(rad))
        {
            mag_inf(z);
        }
        else
        {
            mag_zero(z);
        }
    }
    else
    {
        slong shift, fix;

        shift = _fmpz_sub_small(MAG_EXPREF(mid), MAG_EXPREF(rad));

        /* mid < rad */
        if (shift < 0)
        {
            mag_zero(z);
        }
        else
        {
            mp_limb_t m, xm, rm;

            ARF_GET_TOP_LIMB(xm, mid);
            xm = xm >> (FLINT_BITS - MAG_BITS);

            if (shift <= MAG_BITS)
                rm = (MAG_MAN(rad) >> shift) + 1;
            else
                rm = 1;

            m = xm - rm;

            if (shift > 1)  /* more than one bit cancellation not possible */
            {
                fix = !(m >> (MAG_BITS - 1));
                m <<= fix;
                MAG_MAN(z) = m;
                _fmpz_add_fast(MAG_EXPREF(z), MAG_EXPREF(mid), -fix);
            }
            else if (rm < xm && m > (1 << (MAG_BITS - 4))) /* not too much cancellation */
            {
                fix = MAG_BITS - FLINT_BIT_COUNT(m);
                m <<= fix;
                MAG_MAN(z) = m;
                _fmpz_add_fast(MAG_EXPREF(z), MAG_EXPREF(mid), -fix);
            }
            else
            {
                arf_t t;
                arf_init(t);

                arf_set_mag(t, rad);

                if (arf_sgn(mid) > 0)
                {
                    arf_sub(t, mid, t, MAG_BITS, ARF_RND_DOWN);
                }
                else
                {
                    arf_add(t, mid, t, MAG_BITS, ARF_RND_DOWN);
                    arf_neg(t, t);
                }

                if (arf_sgn(t) <= 0)
                    mag_zero(z);
                else
                    arf_get_mag_lower(z, t);

                arf_clear(t);
            }
        }
    }
}

void
arb_get_mag_lower(mag_t z, const arb_t x)
{
    _arb_get_mag_lower(z, arb_midref(x), arb_radref(x));
}

