/*
    Copyright (C) 2008-2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz_vec.h"

static void _fmpz_vec_set_fft_coeff(fmpz * coeffs_m, slong i,
                          const mp_ptr * coeffs_f, slong limbs, slong sign)
{
    slong size;
    mp_limb_t * data;
    __mpz_struct * mcoeffs_m;

    coeffs_m += i;

    if (sign)
    {
        mp_limb_t halflimb = UWORD(1) << (FLINT_BITS - 1);

        {
            mcoeffs_m = _fmpz_promote(coeffs_m);
            data = FLINT_MPZ_REALLOC(mcoeffs_m, limbs);

			if ((coeffs_f[i][limbs - 1] > halflimb) || coeffs_f[i][limbs])
            {
                mpn_neg(data, coeffs_f[i], limbs);
                mpn_add_1(data, data, limbs, WORD(1));
                size = limbs;
                while ((size) && (data[size - 1] == 0)) size--; /* normalise */
                mcoeffs_m->_mp_size = -size;
                if (size >= WORD(-1)) _fmpz_demote_val(coeffs_m); /* coefficient may be small*/
            }
            else
            {
                flint_mpn_copyi(data, coeffs_f[i], limbs);
                size = limbs;
                while ((size) && (data[size - 1] == WORD(0))) size--; /* normalise */
                mcoeffs_m->_mp_size = size;
                if (size <= 1) _fmpz_demote_val(coeffs_m); /* coefficient may be small */
            }
        }
    }
    else
    {
        {
            mcoeffs_m = _fmpz_promote(coeffs_m);
            data = FLINT_MPZ_REALLOC(mcoeffs_m, limbs);
            flint_mpn_copyi(data, coeffs_f[i], limbs);
            size = limbs;
            while ((size) && (data[size - 1] == WORD(0))) size--; /* normalise */
            mcoeffs_m->_mp_size = size;
            if (size <= 1) _fmpz_demote_val(coeffs_m); /* coefficient may be small */
        }
    }
}

typedef struct
{
    fmpz * coeffs_m;
    const mp_ptr * coeffs_f;
    slong limbs;
    int sign;
}
work_t;

static void
worker(slong i, work_t * work)
{
    _fmpz_vec_set_fft_coeff(work->coeffs_m, i, work->coeffs_f, work->limbs, work->sign);
}

void _fmpz_vec_set_fft(fmpz * coeffs_m, slong length,
                          const mp_ptr * coeffs_f, slong limbs, slong sign)
{
    work_t work;
    slong max_threads;

    work.coeffs_m = coeffs_m;
    work.coeffs_f = coeffs_f;
    work.limbs = limbs;
    work.sign = sign;

    max_threads = flint_get_num_threads();
    max_threads = FLINT_MIN(max_threads, 1e-5 * limbs * length + 1);

    flint_parallel_do((do_func_t) worker, &work, length, max_threads, FLINT_PARALLEL_UNIFORM);
}
