/*
    Copyright (C) 2008-2011 William Hart
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"

static void _fmpz_vec_get_fft_coeff(mp_limb_t ** coeffs_f,
                       const fmpz * coeffs_m, slong l, slong i)
{
    slong size_f = l + 1;
    mp_limb_t * coeff;
    slong size_j, c;
    int signed_c;
    c = coeffs_m[i];
    signed_c = 0;

    if (!COEFF_IS_MPZ(c)) /* coeff is small */
    {
        size_j = 1;

        if (c < 0)
        {
            signed_c = 1;
            c = -c;
            coeff = (mp_limb_t *) &c;
        }
        else
            coeff = (mp_limb_t *) coeffs_m + i;
    }
    else /* coeff is an mpz_t */
    {
        __mpz_struct * mc = COEFF_TO_PTR(c);
        size_j = mc->_mp_size;
        if (size_j < 0)
        {
            signed_c = 1;
            size_j = -size_j;
        }
        coeff = mc->_mp_d;
    }

    if (signed_c) /* write out FFT coefficient, ensuring sign is correct */
    {
        mpn_neg(coeffs_f[i], coeff, size_j);
        flint_mpn_store(coeffs_f[i] + size_j, size_f - size_j, WORD(-1));
    }
    else
    {
        flint_mpn_copyi(coeffs_f[i], coeff, size_j);
        flint_mpn_zero(coeffs_f[i] + size_j, size_f - size_j);
    }
}

typedef struct
{
    mp_limb_t ** coeffs_f;
    const fmpz * coeffs_m;
    slong limbs;
}
work_t;

static void
worker(slong i, work_t * work)
{
    _fmpz_vec_get_fft_coeff(work->coeffs_f, work->coeffs_m, work->limbs, i);
}

void _fmpz_vec_get_fft(mp_limb_t ** coeffs_f,
                       const fmpz * coeffs_m, slong limbs, slong length)
{
    work_t work;
    slong max_threads;

    work.coeffs_f = coeffs_f;
    work.coeffs_m = coeffs_m;
    work.limbs = limbs;

    max_threads = flint_get_num_threads();
    max_threads = FLINT_MIN(max_threads, 1e-5 * limbs * length + 1);

    flint_parallel_do((do_func_t) worker, &work, length, max_threads, FLINT_PARALLEL_UNIFORM);
}
