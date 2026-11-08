/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "fmpz.h"

/* the loop of the case FLINT_BITS < bits < 2 FLINT_BITS of
   _nmod_poly_bit_unpack, each coefficient reduced with RED */
#define _BIT_UNPACK2_LOOP(RED)                                              \
    for (i = nlo; i < nhi; i++)                                             \
    {                                                                       \
        if (current_bit == 0)                                               \
        {                                                                   \
            temp_lower = mpn[current_limb++];                               \
            temp_upper = mpn[current_limb] & mask;                          \
                                                                            \
            RED(res[i], temp_upper, temp_lower, mod);                       \
                                                                            \
            current_bit = bits - FLINT_BITS;                                \
        }                                                                   \
        else if (current_bit > double_boundary_limit_bit)                   \
        {                                                                   \
            /* the coeff will be across two limb boundaries */              \
            temp_lower = mpn[current_limb++] >> current_bit;                \
            temp_lower |=                                                   \
                (mpn[current_limb] << (FLINT_BITS - current_bit));          \
                                                                            \
            temp_upper = mpn[current_limb++] >> current_bit;                \
            temp_upper |=                                                   \
                (mpn[current_limb] << (FLINT_BITS - current_bit));          \
            temp_upper &= mask;                                             \
                                                                            \
            RED(res[i], temp_upper, temp_lower, mod);                       \
                                                                            \
            current_bit += bits - 2 * FLINT_BITS;                           \
        }                                                                   \
        else                                                                \
        {                                                                   \
            /* the coeff will be across one limb boundary */                \
            temp_lower =                                                    \
                (mpn[current_limb] >> current_bit) | (mpn[current_limb + 1] \
                                                      << (FLINT_BITS -      \
                                                          current_bit));    \
            current_limb++;                                                 \
                                                                            \
            temp_upper = mpn[current_limb] >> current_bit;                  \
            temp_upper &= mask;                                             \
                                                                            \
            RED(res[i], temp_upper, temp_lower, mod);                       \
                                                                            \
            current_bit += bits - FLINT_BITS;                               \
            if (current_bit == FLINT_BITS)                                  \
            {                                                               \
                current_bit = 0;                                            \
                current_limb++;                                             \
            }                                                               \
        }                                                                   \
    }

/* the loop of the case 2 FLINT_BITS < bits < 3 FLINT_BITS of
   _nmod_poly_bit_unpack, each coefficient reduced with RED */
#define _BIT_UNPACK3_LOOP(RED)                                              \
    for (i = nlo; i < nhi; i++)                                             \
    {                                                                       \
        if (current_bit == 0)                                               \
        {                                                                   \
            temp_lower = mpn[current_limb++];                               \
            temp_upper = mpn[current_limb++];                               \
            temp_upper2 = mpn[current_limb] & mask;                         \
                                                                            \
            RED(res[i], temp_upper2, temp_upper, temp_lower, mod);          \
                                                                            \
            current_bit = bits - 2 * FLINT_BITS;                            \
        }                                                                   \
        else if (current_bit <= double_boundary_limit_bit)                  \
        {                                                                   \
            /* the coeff will be across two limb boundaries */              \
            temp_lower = mpn[current_limb++] >> current_bit;                \
            temp_lower |=                                                   \
                (mpn[current_limb] << (FLINT_BITS - current_bit));          \
                                                                            \
            temp_upper = mpn[current_limb++] >> current_bit;                \
            temp_upper |=                                                   \
                (mpn[current_limb] << (FLINT_BITS - current_bit));          \
                                                                            \
            temp_upper2 = mpn[current_limb] >> current_bit;                 \
            temp_upper2 &= mask;                                            \
                                                                            \
            RED(res[i], temp_upper2, temp_upper, temp_lower, mod);          \
                                                                            \
            current_bit += bits - 2 * FLINT_BITS;                           \
            if (current_bit == FLINT_BITS)                                  \
            {                                                               \
                current_bit = 0;                                            \
                current_limb++;                                             \
            }                                                               \
        }                                                                   \
        else                                                                \
        {                                                                   \
            /* the coeff will be across three limb boundaries */            \
            temp_lower = mpn[current_limb++] >> current_bit;                \
            temp_lower |=                                                   \
                (mpn[current_limb] << (FLINT_BITS - current_bit));          \
                                                                            \
            temp_upper = mpn[current_limb++] >> current_bit;                \
            temp_upper |=                                                   \
                (mpn[current_limb] << (FLINT_BITS - current_bit));          \
                                                                            \
            temp_upper2 = mpn[current_limb++] >> current_bit;               \
            temp_upper2 |=                                                  \
                (mpn[current_limb] << (FLINT_BITS - current_bit));          \
                                                                            \
            temp_upper2 &= mask;                                            \
                                                                            \
            RED(res[i], temp_upper2, temp_upper, temp_lower, mod);          \
                                                                            \
            current_bit += bits - 3 * FLINT_BITS;                           \
        }                                                                   \
    }

/* NMOD_RED3 for a top word that may be >= n */
#define _BIT_UNPACK_RED3_ANY(r, a_hi, a_me, a_lo, mod)                      \
    do {                                                                    \
        ulong t_hi;                                                         \
        NMOD_RED(t_hi, a_hi, mod);                                          \
        NMOD_RED3(r, t_hi, a_me, a_lo, mod);                                \
    } while (0)

/* Assumes len > 0, bits > 0. */
void
_nmod_poly_bit_unpack(nn_ptr res, slong nlo, slong nhi, nn_srcptr mpn, flint_bitcnt_t bits,
                      nmod_t mod)
{
    slong i;
    ulong current_bit = 0, current_limb = 0;
    ulong temp_lower, temp_upper, temp_upper2;

    if (nlo != 0)
    {
        current_limb = (nlo * bits) / FLINT_BITS;
        current_bit = (nlo * bits) % FLINT_BITS;
        res -= nlo;
    }

    if (bits < FLINT_BITS)
    {
        ulong boundary_limit_bit = FLINT_BITS - bits;

        ulong mask = (UWORD(1) << bits) - UWORD(1);

        for (i = nlo; i < nhi; i++)
        {
            if (current_bit > boundary_limit_bit)
            {
                temp_lower = mpn[current_limb++] >> current_bit;
                temp_upper = mpn[current_limb] << (FLINT_BITS - current_bit);

                temp_upper |= temp_lower;
                temp_upper &= mask;

                NMOD_RED(res[i], temp_upper, mod);

                current_bit += bits - FLINT_BITS;
            }
            else
            {
                /* the coeff will fit in the current limb */
                temp_upper = (mpn[current_limb] >> current_bit) & mask;

                NMOD_RED(res[i], temp_upper, mod);

                current_bit += bits;

                if (current_bit == FLINT_BITS)
                {
                    current_bit = 0;
                    current_limb++;
                }
            }
        }
    }
    else if (bits == FLINT_BITS)
    {
        for (i = nlo; i < nhi; i++)
            NMOD_RED(res[i], mpn[i], mod);
    }
    else if (bits == 2 * FLINT_BITS)
    {
        for (i = nlo; i < nhi; i++)
        {
            NMOD2_RED2(res[i], mpn[current_limb + 1], mpn[current_limb], mod);
            current_limb += 2;
        }
    }
    else if (bits < 2 * FLINT_BITS) /* FLINT_BITS < bits < 2*FLINT_BITS */
    {
        ulong double_boundary_limit_bit = 2 * FLINT_BITS - bits;

        ulong mask = (UWORD(1) << (bits - FLINT_BITS)) - UWORD(1);

        /* the coefficients are below 2^bits: NMOD_RED2_BOUNDED applies if
           bits <= nu + FLINT_BITS - 2 and nu < FLINT_BITS, nu the bit length
           of n (then nu >= 3, since bits > FLINT_BITS) */
        if (mod.norm != 0 && bits <= 2 * FLINT_BITS - 2 - mod.norm)
        {
            _BIT_UNPACK2_LOOP(NMOD_RED2_BOUNDED)
        }
        else
        {
            _BIT_UNPACK2_LOOP(NMOD2_RED2)
        }
    }
    else                        /* 2*FLINT_BITS < bits < 3*FLINT_BITS */
    {
        ulong double_boundary_limit_bit = 3 * FLINT_BITS - bits;

        ulong mask = (UWORD(1) << (bits - 2 * FLINT_BITS)) - UWORD(1);

        /* the fields are below 2^bits: NMOD_RED3_BOUNDED applies if
           bits <= nu + 2 FLINT_BITS - 2 and nu < FLINT_BITS, nu the bit length
           of n (then nu >= 3); otherwise NMOD_RED3, which needs the top word
           below n: it is below 2^(bits - 2 FLINT_BITS) <= 2^(nu - 1) <= n if
           bits - 2 FLINT_BITS < nu, and else is reduced first */
        if (mod.norm != 0 && bits <= 3 * FLINT_BITS - 2 - mod.norm)
        {
            _BIT_UNPACK3_LOOP(NMOD_RED3_BOUNDED)
        }
        else if (bits - 2 * FLINT_BITS < NMOD_BITS(mod))
        {
            _BIT_UNPACK3_LOOP(NMOD_RED3)
        }
        else
        {
            _BIT_UNPACK3_LOOP(_BIT_UNPACK_RED3_ANY)
        }
    }
}

#undef _BIT_UNPACK2_LOOP
#undef _BIT_UNPACK3_LOOP
#undef _BIT_UNPACK_RED3_ANY

void
nmod_poly_bit_unpack(nmod_poly_t poly, const fmpz_t f, flint_bitcnt_t bit_size)
{
    slong len;
    mpz_t tmp;

    if (fmpz_sgn(f) < 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_bit_unpack). f < 0.\n");
    }

    if (bit_size == 0 || fmpz_is_zero(f))
    {
        nmod_poly_zero(poly);
        return;
    }

    len = (fmpz_bits(f) + bit_size - 1) / bit_size;

    mpz_init2(tmp, bit_size*len);
    flint_mpn_zero(tmp->_mp_d, tmp->_mp_alloc);
    fmpz_get_mpz(tmp, f);

    nmod_poly_fit_length(poly, len);

    _nmod_poly_bit_unpack(poly->coeffs, 0, len, tmp->_mp_d, bit_size, poly->mod);
    poly->length = len;
    _nmod_poly_normalise(poly);

    mpz_clear(tmp);
}
