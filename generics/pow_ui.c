/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2013 William Hart

******************************************************************************/

#include "generics.h"



static void
pow_binexp_generic(elem_ptr res, elem_srcptr base, ulong exp, const ring_t ring)
{
    ulong bit = ~((~0UL) >> 1);
    elem_ptr v, R, S, T;

    if (exp <= 2)
    {
        if (exp == 0)
            elem_one(res, ring);
        else if (exp == 1)
            elem_set(res, base, ring);
        else if (exp == 2)
            elem_mul(res, base, base, ring);
        return;
    }

    if (elem_is_zero(base, ring))
    {
        elem_zero(res, ring);
        return;
    }

    if (elem_is_one(base, ring))
    {
        elem_one(res, ring);
        return;
    }

    if (res == base)
    {
        elem_ptr tmp;
        ELEM_TMP_INIT(tmp, ring);
        pow_binexp_generic(tmp, base, exp, ring);
        elem_swap(res, tmp, ring);
        ELEM_TMP_CLEAR(tmp, ring);
        return;
    }

    ELEM_TMP_INIT(v, ring);

    while ((bit & exp) == 0UL)
        bit >>= 1;
    bit >>= 1;

    {
        unsigned int swaps = 0U;
        ulong bit2 = bit;
        if ((bit2 & exp))
            swaps = ~swaps;
        while (bit2 >>= 1)
            if ((bit2 & exp) == 0UL)
                swaps = ~swaps;
        
        if (swaps == 0U)
        {
            R = res;
            S = v;
        }
        else
        {
            R = v;
            S = res;
        }
    }

    elem_mul(R, base, base, ring);

    if ((bit & exp))
    {
        elem_mul(S, R, base, ring);
        T = R;
        R = S;
        S = T;
    }
    
    while ((bit >>= 1))
    {
        if ((bit & exp))
        {
            elem_mul(S, R, R, ring);
            elem_mul(R, S, base, ring);
        }
        else
        {
            elem_mul(S, R, R, ring);
            T = R;
            R = S;
            S = T;
        }
    }

    ELEM_TMP_CLEAR(v, ring);
}


void
elem_pow_ui(elem_ptr res, elem_srcptr op, ulong exp, const ring_t ring)
{
    switch (ring->type)
    {
        case TYPE_FMPZ:
            fmpz_pow_ui(res, op, exp);
            break;

        case TYPE_MPZ:
            mpz_pow_ui(res, op, exp);
            break;

        case TYPE_POLY:
            elem_poly_pow_ui(res, op, exp, ring);
            break;

        case TYPE_MOD:
            {
                switch (RING_PARENT(ring)->type)
                {
                    case TYPE_LIMB:
                        *((mp_ptr) res) = n_powmod2_preinv(*((mp_srcptr) op),
                            exp, ring->nmod.n, ring->nmod.ninv);
                        break;

                    case TYPE_FMPZ:
                        fmpz_powm_ui(res, op, exp, RING_MODULUS(ring));
                        break;

                    default:
                        NOT_IMPLEMENTED("pow_ui (mod)", ring);
                }
            }
            break;

        case TYPE_FRAC:
            elem_pow_ui(NUMER(res, ring), NUMER(op, ring), exp, RING_NUMER(ring));
            elem_pow_ui(DENOM(res, ring), DENOM(op, ring), exp, RING_DENOM(ring));
            break;

        default:
            pow_binexp_generic(res, op, exp, ring);
    }
}

