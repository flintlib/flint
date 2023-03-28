/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_add(fmpz * res, const fmpz * vec1, const fmpz * vec2, slong len2)
{
    slong i;
    for (i = 0; i < len2; i++)
        fmpz_add(res + i, vec1 + i, vec2 + i);
}

#if defined(__AVX2__) || defined(__AVX512F__)

/* Note that the representation of an mpz in an fmpz is that the MSB is 0 while
 * the second MSB is 1. Hence we can simplify this check to just check if x >
 * COEFF_MAX. This should probably be done throughout the repository, but it has
 * been rejected previously. */
#define IS_MPZ(x) ((x) > COEFF_MAX)
#define TURNED_MPZ(x) ((x) > COEFF_MAX || (x) < COEFF_MIN)

__GMP_DECLSPEC extern void * (*__gmp_reallocate_func) (void *, size_t, size_t);
#define ALLOC(x) ((x)->_mp_alloc)
#define SIZ(x) ((x)->_mp_size)
#define PTR(x) ((x)->_mp_d)

static void FLINT_NOINLINE
fmpz_add_large2(fmpz * rp, const fmpz f_ip)
{
    fmpz f_rp = *rp;
    ulong f_abs;
    mpz_ptr mp;
    slong mp_size;
    mp_ptr mp_d;
    slong mp_abssize;

    if (!IS_MPZ(f_ip))
    {
        /* rp big, ip small */
        f_abs = FLINT_ABS(f_ip);
        mp = COEFF_TO_PTR(f_rp);
        mp_size = SIZ(mp);
        mp_d = PTR(mp);
        mp_abssize = FLINT_ABS(mp_size);

        if ((mp_size ^ f_ip) < 0 && mp_abssize == 1)
        {
            /* Treat case when demotion is possible */
            ulong tmp = *mp_d - f_abs;
            if (tmp <= (ulong) COEFF_MAX)
            {
                _fmpz_clear_mpz(f_rp);
                *rp = (mp_size > 0) ? tmp : -tmp;
            }
            else
                *mp_d = tmp;
            return;
        }
        else if ((mp_size ^ f_ip) < 0)
            goto sub_1;
        else
        {
            if (mp_abssize == ALLOC(mp))
            {
                mp_d = PTR(mp) = __gmp_reallocate_func(mp_d,
                        sizeof(mp_limb_t) * ALLOC(mp),
                        sizeof(mp_limb_t) * (ALLOC(mp) + 1));
                ALLOC(mp)++;
            }
            goto add_1;
        }
add_1:
        {
            /* Carry can happen */
            ulong carry;
            __GMPN_ADD_1(carry, mp_d, mp_d, mp_abssize, f_abs);
            if (carry)
            {
                mp_d[mp_abssize] = carry;
                SIZ(mp) += FLINT_SGN(SIZ(mp));
            }
            return;
        }
sub_1:
        {
            /* Borrow cannot happen, but the most significant limb can become zero */
            ulong FLINT_SET_BUT_UNUSED(borrow_unused);
            __GMPN_SUB_1(borrow_unused, mp_d, mp_d, mp_abssize, f_abs);
            if (mp_d[mp_abssize - 1] == 0)
            {
                SIZ(mp) -= FLINT_SGN(SIZ(mp));
            }
            return;
        }
    }
    else if (!IS_MPZ(f_rp))
    {
        /* rp small, ip big */
        mpz_ptr tmp_mp;
        slong tmp_mp_size;
        mp_ptr tmp_mp_d;
        slong tmp_mp_abssize;
        slong ix;

        f_abs = FLINT_ABS(f_rp);
        tmp_mp = COEFF_TO_PTR(f_ip);
        tmp_mp_size = SIZ(tmp_mp);
        tmp_mp_d = PTR(tmp_mp);
        tmp_mp_abssize = FLINT_ABS(tmp_mp_size);

        if ((f_rp ^ tmp_mp_size) < 0 && tmp_mp_abssize == 1)
        {
            /* Treat case when promotion can be avoided */
            ulong tmp = *tmp_mp_d - f_abs;
            if (tmp <= (ulong) COEFF_MAX)
                *rp = (tmp_mp_size > 0) ? tmp : -tmp;
            else
            {
                mp = _fmpz_new_mpz();
                *rp = PTR_TO_COEFF(mp);
                SIZ(mp) = tmp_mp_size;
                *PTR(mp) = tmp;
            }

            return;
        }

        /* Result won't fit in a small fmpz */
        mp = _fmpz_new_mpz2(tmp_mp_abssize + 1);
        *rp = PTR_TO_COEFF(mp);
        mp_size = SIZ(mp) = tmp_mp_size;
        mp_abssize = tmp_mp_abssize;
        mp_d = PTR(mp);
        for (ix = 0; ix < mp_abssize; ix++)
            mp_d[ix] = tmp_mp_d[ix];

        if ((mp_size ^ f_rp) < 0)
            goto sub_1;
        else
            goto add_1;
    }
    else
    {
        /* NOTE: We could probably optimize for larger fmpz, but then the
         * function body becomes larger as well. We mostly want to optimize for
         * small fmpz. */
        mpz_ptr tmp_mp;

        mp = COEFF_TO_PTR(f_rp);
        tmp_mp = COEFF_TO_PTR(f_ip);

        mpz_add(mp, mp, tmp_mp);
        _fmpz_demote_val(rp); /* FIXME: should be inlined */
    }
}
#undef ALLOC
#undef SIZ
#undef PTR

#include <immintrin.h>

#ifdef __AVX512F__
# define NUM 8
# define simd_type __m512i
# define simd_load_qword _mm512_set1_epi64
# define simd_load_mem(x) _mm512_loadu_epi64(x)
# define simd_store_mem(x, y) _mm512_storeu_epi64((x), (y))
# define simd_cmpgt _mm512_cmpgt_epi64_mask
# define simd_add _mm512_add_epi64
#else
# define NUM 4
# define simd_type __m256i
# define simd_load_qword _mm256_set1_epi64x
# define simd_load_mem(x) _mm256_loadu_si256((simd_type *) (x))
# define simd_store_mem(x, y) _mm256_storeu_si256((simd_type *) (x), (y))
# define simd_cmpgt _mm256_cmpgt_epi64
# define simd_testz _mm256_testz_si256
# define simd_add _mm256_add_epi64
#endif

void
_fmpz_vec_add2(fmpz * rp, const fmpz * ip, slong len)
{
    slong ix;
    simd_type simd_coeff_max, simd_offset, simd_cmp_got_large;
    slong counter = 0;

    for (ix = 0; ix < len % NUM; ix++)
    {
        if (!IS_MPZ(rp[ix]) && !IS_MPZ(ip[ix]))
        {
            rp[ix] += ip[ix];

            if (TURNED_MPZ(rp[ix]))
            {
                mpz_ptr mrp = _fmpz_new_mpz();

                /* At least one limb is allocated when calling _fmpz_new_mpz */
                mrp->_mp_d[0] = FLINT_ABS(rp[ix]);
                mrp->_mp_size = FLINT_SGN(rp[ix]);
                rp[ix] = PTR_TO_COEFF(mrp);
            }
        }
        else if (ip[ix] != 0)
            fmpz_add_large2(rp + ix, ip[ix]);
        counter++;
    }

    rp += len % NUM;
    ip += len % NUM;
    len -= len % NUM;

    simd_coeff_max = simd_load_qword(COEFF_MAX);
    simd_offset = simd_load_qword(WORD_MIN - COEFF_MIN);
    simd_cmp_got_large = simd_load_qword(COEFF_MAX + (WORD_MIN - COEFF_MIN));

    for (ix = 0; ix < len / NUM; ix++)
    {
        simd_type s1, s2;
#ifdef __AVX512F__
#else
        simd_type c1, c2;
#endif
        int zf;

        s1 = simd_load_mem(rp);
        s2 = simd_load_mem(ip);

#ifdef __AVX512F__
        zf = !(simd_cmpgt(s1, simd_coeff_max) && simd_cmpgt(s2, simd_coeff_max));
#else
        c1 = simd_cmpgt(s1, simd_coeff_max);
        c2 = simd_cmpgt(s2, simd_coeff_max);

        zf = simd_testz(c1, c1) && simd_testz(c2, c2);
#endif

        if (zf)
        {
            s1 = simd_add(s1, s2);
            s2 = simd_add(s1, simd_offset);
#ifdef __AVX512F__
            zf = simd_cmpgt(s2, simd_cmp_got_large);
#else
            c1 = simd_cmpgt(s2, simd_cmp_got_large);
            zf = simd_testz(c1, c1);
#endif
            if (zf)
                simd_store_mem(rp, s1);
            else
                goto not_all_small;
        }
        else
        {
not_all_small:
            slong jx;
            for (jx = 0; jx < NUM; jx++)
            {
                if (!IS_MPZ(rp[jx]) && !IS_MPZ(ip[jx]))
                {
                    rp[jx] += ip[jx];

                    if (TURNED_MPZ(rp[jx]))
                    {
                        mpz_ptr mrp = _fmpz_new_mpz();

                        /* At least one limb is allocated when calling _fmpz_new_mpz */
                        mrp->_mp_d[0] = FLINT_ABS(rp[jx]);
                        mrp->_mp_size = FLINT_SGN(rp[jx]);
                        rp[jx] = PTR_TO_COEFF(mrp);
                    }
                }
                else if (ip[jx] != 0)
                    fmpz_add_large2(rp + jx, ip[jx]);
            }
        }

        rp += NUM;
        ip += NUM;
    }
}

#undef NUM
#undef IS_MPZ
#undef simd_type
#undef simd_load_qword
#undef simd_load_mem
#undef simd_store_mem
#undef simd_cmpgt
#undef simd_testz
#undef simd_add
#endif
