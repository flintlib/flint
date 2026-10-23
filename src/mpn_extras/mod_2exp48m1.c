/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#if FLINT_BITS == 64

#if defined(__AVX2__)
#include <immintrin.h>
#endif

/*
    Returns a value congruent to a modulo 2^48 - 1 (the same quantity as
    GMP's internal mpn_mod_34lsub1). Since 2^48 == 1 and 2^64 == 2^16
    modulo 2^48 - 1, limb i contributes with weight 2^(16 (i mod 3)). With
    AVX2, 12 limbs per step are added into three vector accumulators whose
    lanes carry a fixed pattern of weights; each limb is split into its low
    48 and high 16 bits so that no lane overflows for fewer than 2^15
    limbs, which is the period of the folding. This runs at about half the
    time of GMP's assembly from 128 limbs on (2048 limbs: 0.25 vs 0.49 us),
    so it is preferred to mpn_mod_34lsub1 whenever AVX2 is enabled; short
    inputs and the tail use add/adc carry chains for the three weight
    classes. Requires n < 2^47.
*/
mp_limb_t
flint_mpn_mod_2exp48m1(mp_srcptr a, mp_size_t n)
{
    const mp_limb_t M = (UWORD(1) << 48) - 1;
    mp_limb_t c0 = 0, c1 = 0, c2 = 0, r;
    mp_size_t i = 0;

#if FLINT_HAVE_NATIVE_mpn_mod_34lsub1
    /* GMP's assembly is faster for short inputs */
    if (n < 48)
        return mpn_mod_34lsub1(a, n);
#endif

#if defined(__AVX2__)
    if (n >= 48)
    {
        /* lane weight classes: v0 = (0,1,2,0), v1 = (1,2,0,1), v2 = (2,0,1,2) */
        const __m256i mask = _mm256_set1_epi64x(M);
        __m256i A0 = _mm256_setzero_si256(), A1 = A0, A2 = A0;
        mp_limb_t t0[4], t1[4], t2[4];
        mp_size_t chunk;

        while (i + 12 <= n)
        {
            chunk = FLINT_MIN((n - i) / 12, 8000) * 12;
            for ( ; chunk > 0; chunk -= 12, i += 12)
            {
                __m256i v0 = _mm256_loadu_si256((const __m256i *) (a + i));
                __m256i v1 = _mm256_loadu_si256((const __m256i *) (a + i + 4));
                __m256i v2 = _mm256_loadu_si256((const __m256i *) (a + i + 8));
                A0 = _mm256_add_epi64(A0, _mm256_add_epi64(_mm256_and_si256(v0, mask), _mm256_srli_epi64(v0, 48)));
                A1 = _mm256_add_epi64(A1, _mm256_add_epi64(_mm256_and_si256(v1, mask), _mm256_srli_epi64(v1, 48)));
                A2 = _mm256_add_epi64(A2, _mm256_add_epi64(_mm256_and_si256(v2, mask), _mm256_srli_epi64(v2, 48)));
            }
            /* fold each lane back below 2^49 */
            A0 = _mm256_add_epi64(_mm256_and_si256(A0, mask), _mm256_srli_epi64(A0, 48));
            A1 = _mm256_add_epi64(_mm256_and_si256(A1, mask), _mm256_srli_epi64(A1, 48));
            A2 = _mm256_add_epi64(_mm256_and_si256(A2, mask), _mm256_srli_epi64(A2, 48));
        }

        _mm256_storeu_si256((__m256i *) t0, A0);
        _mm256_storeu_si256((__m256i *) t1, A1);
        _mm256_storeu_si256((__m256i *) t2, A2);
        c0 = t0[0] + t0[3] + t1[2] + t2[1];
        c1 = t0[1] + t1[0] + t1[3] + t2[2];
        c2 = t0[2] + t1[1] + t2[0] + t2[3];
        c0 = (c0 & M) + (c0 >> 48);
        c1 = (c1 & M) + (c1 >> 48);
        c2 = (c2 & M) + (c2 >> 48);
    }
#endif

    {
        mp_limb_t h0 = 0, l0 = 0, h1 = 0, l1 = 0, h2 = 0, l2 = 0;

        /* double-limb sums of the remaining limbs in each weight class
           (i is a multiple of 3 here) */
        for ( ; i + 3 <= n; i += 3)
        {
            add_ssaaaa(h0, l0, h0, l0, 0, a[i]);
            add_ssaaaa(h1, l1, h1, l1, 0, a[i + 1]);
            add_ssaaaa(h2, l2, h2, l2, 0, a[i + 2]);
        }
        if (i < n)
            add_ssaaaa(h0, l0, h0, l0, 0, a[i]);
        if (i + 1 < n)
            add_ssaaaa(h1, l1, h1, l1, 0, a[i + 1]);

        /* h 2^64 + l == h 2^16 + (l mod 2^48) + (l div 2^48) with h < n */
        c0 += (l0 & M) + (l0 >> 48) + (h0 << 16);
        c1 += (l1 & M) + (l1 >> 48) + (h1 << 16);
        c2 += (l2 & M) + (l2 >> 48) + (h2 << 16);
        c0 = (c0 & M) + (c0 >> 48);
        c1 = (c1 & M) + (c1 >> 48);
        c2 = (c2 & M) + (c2 >> 48);
    }

    /* c0 + c1 2^16 + c2 2^32 with each c < 2^49: split c1, c2 across the
       48-bit boundary */
    r = c0 + ((c1 << 16) & M) + (c1 >> 32) + ((c2 << 32) & M) + (c2 >> 16);
    return r;
}

#else
/* keep the translation unit non-empty on 32-bit builds */
typedef int _flint_mpn_mod_2exp48m1_unused;
#endif
