/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"


void _fmpz_multi_mod_precomp(
    fmpz * outputs,
    const fmpz_multi_mod_t P,
    const fmpz_t input,
    int sign,
    fmpz * T)
{
    slong i, a, b;
    slong len = P->length;
    _fmpz_multi_mod_instr * instr = P->prog;
    fmpz * t1 = T + P->temp1loc;
    unsigned char * org;
    TMP_INIT;

    TMP_START;

    /*
        Efficiently propogate small inputs without copying:
        ord[i] = 1 means T[i] should be read from input
    */
    org = TMP_ARRAY_ALLOC(P->localsize, unsigned char);

#if FLINT_WANT_ASSERT
    for (i = 0; i < P->localsize; i++)
        org[i] = 2;
#endif

    for (i = 0; i < len; i++)
    {
        a = P->prog[i].in_idx;
        b = P->prog[i].out_idx;

        FLINT_ASSERT(a < 1 || org[a] < 2);

        if (a > 0 && org[a] == 0)
        {
            /* read input from T[a] */

            if (b < 0)
            {
                _fmpz_smod(outputs - b - 1, T + a, instr[i].modulus, sign, t1);
            }
            else
            {
                org[b] = 0;
                fmpz_tdiv_qr(t1, T + b, T + a, instr[i].modulus);
            }
        }
        else
        {
            /* read input from input */

            if (b < 0)
            {
                _fmpz_smod(outputs - b - 1, input, instr[i].modulus, sign, t1);
            }
            else if (fmpz_cmpabs(instr[i].modulus, input) > 0)
            {
                org[b] = 1;
            }
            else
            {
                org[b] = 0;
                fmpz_tdiv_qr(t1, T + b, input, instr[i].modulus);
            }
        }
    }

    TMP_END;
}

void fmpz_multi_mod_precomp(
    fmpz * outputs,
    const fmpz_multi_mod_t P,
    const fmpz_t input,
    int sign)
{
    slong i;
    fmpz * tmp;
    TMP_INIT;

    TMP_START;
    tmp = TMP_ARRAY_ALLOC(P->localsize, fmpz);
    for (i = 0; i < P->localsize; i++)
        fmpz_init(tmp + i);

    _fmpz_multi_mod_precomp(outputs, P, input, sign, tmp);

    for (i = 0; i < P->localsize; i++)
        fmpz_clear(tmp + i);

    TMP_END;
}

