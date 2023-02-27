/*
    Copyright (C) 2019-2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void fmpz_deprecated_multi_crt_init(fmpz_multi_crt_t P)
{
    fmpz_multi_CRT_init(P);
}

void fmpz_deprecated_multi_crt_clear(fmpz_multi_crt_t P)
{
    fmpz_multi_CRT_clear(P);
}

int fmpz_deprecated_multi_crt_precompute(
    fmpz_multi_crt_t P,
    const fmpz * moduli,
    slong len)
{
    return fmpz_multi_CRT_precompute(P, moduli, len);
}

int fmpz_deprecated_multi_crt_precompute_p(
    fmpz_multi_crt_t P,
    const fmpz * const * moduli,
    slong len)
{
    int success;
    slong i;
    fmpz * m = FLINT_ARRAY_ALLOC(len, fmpz);

    for (i = 0; i < len; i++)
        m[i] = *moduli[i];

    success = fmpz_multi_CRT_precompute(P, m, len);

    flint_free(m);

    return success;
}

void fmpz_deprecated_multi_crt_precomp(
    fmpz_t output,
    const fmpz_multi_crt_t P,
    const fmpz * inputs)
{
    fmpz_multi_CRT_precomp(output, P, inputs, 1);
}

void fmpz_deprecated_multi_crt_precomp_p(
    fmpz_t output,
    const fmpz_multi_crt_t P,
    const fmpz * const * inputs)
{
    slong i;
    fmpz * ins = FLINT_ARRAY_ALLOC(P->moduli_count, fmpz);

    for (i = 0; i < P->moduli_count; i++)
        ins[i] = *inputs[i];

    fmpz_multi_CRT_precomp(output, P, ins, 1);

    flint_free(ins);
}

int fmpz_deprecated_multi_crt(
    fmpz_t output,
    const fmpz * moduli,
    const fmpz * values,
    slong len)
{
    return fmpz_multi_CRT(output, moduli, values, len, 1);
}

void _fmpz_deprecated_multi_crt_run(
    fmpz * outputs,
    const fmpz_multi_crt_t P,
    const fmpz * inputs)
{
    _fmpz_multi_CRT_precomp(outputs, P, inputs, 1);
}

slong _fmpz_deprecated_multi_crt_local_size(const fmpz_multi_crt_t CRT)
{
    return CRT->localsize;
}

void _fmpz_deprecated_multi_crt_run_p(
    fmpz * outputs,
    const fmpz_multi_crt_t P,
    const fmpz * const * inputs)
{
    slong i;
    fmpz * ins = FLINT_ARRAY_ALLOC(P->moduli_count, fmpz);

    for (i = 0; i < P->moduli_count; i++)
        ins[i] = *inputs[i];

    _fmpz_multi_CRT_precomp(outputs, P, ins, 1);

    flint_free(ins);
}

