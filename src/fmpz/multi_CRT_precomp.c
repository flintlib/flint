/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"


void _fmpz_multi_CRT_precomp(
    fmpz * outputs,
    const fmpz_multi_CRT_t P,
    const fmpz * inputs,
    int sign)
{
    slong i, a, b, c;
    slong len = P->length;
    const fmpz * m = P->moduli;
    const fmpz * mf = P->fracmoduli;
    fmpz * A, * B, * C, * t1, * t2, * t3, * t4;

    t1 = outputs + P->temp1loc;
    t2 = outputs + P->temp2loc;
    t3 = outputs + P->temp3loc;
    t4 = outputs + P->temp4loc;

    FLINT_ASSERT(len < 1 || P->good);

    if (len > 0)
    {
        for (i = P->moduli_count - 1; i > 0; i--)
        {
            if (!fmpz_equal(inputs + 0, inputs + i))
                goto doit;
        }
    }

    _fmpz_smod(outputs + 0, inputs + 0, P->final_modulus, sign, t4);
    return;

doit:

    for (i = 0; i < len; i++)
    {
        a = P->prog[i].a_idx;
        b = P->prog[i].b_idx;
        c = P->prog[i].c_idx;

        A = outputs + a;
        B = outputs + b;
        C = outputs + c;

        if (b < 0)
        {
            b = -b - 1;
            B = t1;

            fmpz_mul(t3, inputs + b, mf + b);
            _fmpz_smod(B, t3, m + b, sign, t4);
        }

        if (c < 0)
        {
            c = -c - 1;
            C = t2;

            fmpz_mul(t3, inputs + c, mf + c);
            _fmpz_smod(C, t3, m + c, sign, t4);
        }

        /* A = B*c_m + C*b_m */
        fmpz_mul(A, B, P->prog[i].c_modulus);
        fmpz_mul(t3, C, P->prog[i].b_modulus);
        fmpz_add(A, A, t3);
    }

    _fmpz_smod(outputs + 0, A, P->final_modulus, sign, t4);
}


void fmpz_multi_CRT_precomp(
    fmpz_t output,
    const fmpz_multi_CRT_t P,
    const fmpz * inputs,
    int sign)
{
    slong i;
    fmpz * out;
    TMP_INIT;

    TMP_START;
    out = TMP_ARRAY_ALLOC(P->localsize, fmpz);
    for (i = 0; i < P->localsize; i++)
        fmpz_init(out + i);

    fmpz_swap(out + 0, output);
    _fmpz_multi_CRT_precomp(out, P, inputs, sign);
    fmpz_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
        fmpz_clear(out + i);

    TMP_END;
}

