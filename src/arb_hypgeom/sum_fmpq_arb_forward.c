/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

void
arb_hypgeom_sum_fmpq_arb_forward(arb_t res, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
{
    arb_t t, u, sden;
    slong i, k, num_max_bits, den_max_bits, bits1, bits2, Nbits, den_prec;

    if (N <= 1)
    {
        if (N == 1)
            arb_one(res);
        else
            arb_zero(res);
        return;
    }

    arb_init(t);
    arb_init(u);
    arb_init(sden);

    Nbits = FLINT_BIT_COUNT(N);
    num_max_bits = 0;
    for (i = 0; i < alen; i++)
    {
        bits1 = fmpz_bits(fmpq_numref(a + i));
        bits2 = fmpz_bits(fmpq_denref(a + i)) + Nbits;
        num_max_bits += FLINT_MAX(bits1, bits2) + 1;
    }

    den_max_bits = 0;
    for (i = 0; i < blen; i++)
    {
        bits1 = fmpz_bits(fmpq_numref(b + i));
        bits2 = fmpz_bits(fmpq_denref(b + i)) + Nbits;
        den_max_bits += FLINT_MAX(bits1, bits2) + 1;
    }

    den_prec = FLINT_MAX(128, prec + 2 * Nbits);

    if (num_max_bits <= FLINT_BITS - 1 && den_max_bits <= FLINT_BITS - 1)
    {
        slong num, den;

        den = 1;
        num = 1;
        for (i = 0; i < alen; i++)
            den *= *fmpq_denref(a + i);
        for (i = 0; i < blen; i++)
            num *= *fmpq_denref(b + i);

        if (reciprocal)
        {
            arb_mul_si(u, z, den, prec);
            arb_set_si(t, num);
            arb_div(u, t, u, prec);
        }
        else
        {
            arb_mul_si(u, z, num, prec);
            arb_div_si(u, u, den, prec);
        }

        arb_one(res);
        arb_one(t);
        arb_one(sden);

        for (k = 1; k < N && !arb_is_zero(t); k++)
        {
            arb_mul(t, t, u, prec);

            if (alen != 0)
            {
                num = *fmpq_denref(a + 0) * (k - 1) + *fmpq_numref(a + 0);

                for (i = 1; i < alen; i++)
                    num *= (*fmpq_denref(a + i) * (k - 1) + *fmpq_numref(a + i));

                arb_mul_si(t, t, num, prec);
            }

            if (blen != 0)
            {
                den = *fmpq_denref(b + 0) * (k - 1) + *fmpq_numref(b + 0);

                for (i = 1; i < blen; i++)
                    den *= (*fmpq_denref(b + i) * (k - 1) + *fmpq_numref(b + i));

                arb_fma_si(res, res, den, t, prec);
                arb_mul_si(sden, sden, den, den_prec);
            }
            else
            {
                arb_add(res, res, t, prec);
            }
        }

        arb_div(res, res, sden, prec);
    }
    else
    {
        fmpz_t num, den, c;

        fmpz_init(num);
        fmpz_init(den);
        fmpz_init(c);

        fmpz_one(den);
        fmpz_one(num);
        for (i = 0; i < alen; i++)
            fmpz_mul(den, den, fmpq_denref(a + i));
        for (i = 0; i < blen; i++)
            fmpz_mul(num, num, fmpq_denref(b + i));

        if (reciprocal)
        {
            arb_mul_fmpz(u, z, den, prec);
            arb_set_fmpz(t, num);
            arb_div(u, t, u, prec);
        }
        else
        {
            arb_mul_fmpz(u, z, num, prec);
            arb_div_fmpz(u, u, den, prec);
        }

        arb_one(res);
        arb_one(t);
        arb_one(sden);

        for (k = 1; k < N && !arb_is_zero(t); k++)
        {
            arb_mul(t, t, u, prec);

            if (alen != 0)
            {
                fmpz_mul_ui(num, fmpq_denref(a + 0), k - 1);
                fmpz_add(num, num, fmpq_numref(a + 0));

                for (i = 1; i < alen; i++)
                {
                    fmpz_mul_ui(c, fmpq_denref(a + i), k - 1);
                    fmpz_add(c, c, fmpq_numref(a + i));
                    fmpz_mul(num, num, c);
                }

                arb_mul_fmpz(t, t, num, prec);
            }

            if (blen != 0)
            {
                fmpz_mul_ui(den, fmpq_denref(b + 0), k - 1);
                fmpz_add(den, den, fmpq_numref(b + 0));

                for (i = 1; i < blen; i++)
                {
                    fmpz_mul_ui(c, fmpq_denref(b + i), k - 1);
                    fmpz_add(c, c, fmpq_numref(b + i));
                    fmpz_mul(den, den, c);
                }

                arb_fma_fmpz(res, res, den, t, prec);
                arb_mul_fmpz(sden, sden, den, den_prec);
            }
            else
            {
                arb_add(res, res, t, prec);
            }
        }

        arb_div(res, res, sden, prec);

        fmpz_clear(num);
        fmpz_clear(den);
        fmpz_clear(c);
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(sden);
}
