/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "arith.h"
#include "thread_support.h"
#include "impl.h"

/* Representation numbers of sums of three squares. */

/* r3(N) for N < X: theta(x)^2 is formed in linear time and multiplied
   once by theta(x) (arith_sum_of_squares_vec) */
uint32_t *
_acb_mzeta_r3_all_table(ulong X)
{
    slong n = FLINT_MAX(X, 1), k;
    fmpz * r = _fmpz_vec_init(n);
    uint32_t * r3 = flint_malloc(sizeof(uint32_t) * n);

    arith_sum_of_squares_vec(r, 3, n);
    for (k = 0; k < n; k++)
        r3[k] = fmpz_get_ui(r + k);
    _fmpz_vec_clear(r, n);
    return r3;
}

/* r3(3 + 8k) for k < K: with A = 2a + 1, (A^2 - 1)/8 = a(a + 1)/2, so
   r3(8k + 3) = 8 [x^k] T(x)^3, T(x) = sum_{a >= 0} x^(a(a+1)/2) (Gauss);
   T(x)^2 directly from the about pi K / 2 pairs, then one multiplication */
uint32_t *
_acb_mzeta_r3_3mod8_table(ulong K)
{
    slong n = FLINT_MAX(K, 1), k;
    fmpz * T = _fmpz_vec_init(n);
    fmpz * T2 = _fmpz_vec_init(n);
    fmpz * T3 = _fmpz_vec_init(n);
    uint32_t * r3 = flint_malloc(sizeof(uint32_t) * n);
    ulong a, b, sa, sb;

    for (a = 0; (sa = a * (a + 1) / 2) < (ulong) n; a++)
    {
        fmpz_one(T + sa);
        for (b = 0; (sb = sa + b * (b + 1) / 2) < (ulong) n; b++)
            fmpz_add_ui(T2 + sb, T2 + sb, 1);
    }
    _fmpz_poly_mullow(T3, T2, n, T, n, n);
    for (k = 0; k < n; k++)
        r3[k] = 8 * fmpz_get_ui(T3 + k);
    _fmpz_vec_clear(T, n);
    _fmpz_vec_clear(T2, n);
    _fmpz_vec_clear(T3, n);
    return r3;
}

/* t2(j) = #{a, b >= 0 : a(a+1)/2 + b(b+1)/2 = j} = sum_{d | 4j+1} chi_4(d)
   (from 8j + 2 = (2a+1)^2 + (2b+1)^2 and r2(n) = 4 sum_{d|n} chi_4(d)).
   The divisors pair up as d, (4j+1)/d with equal chi_4, so it suffices to
   sieve over d <= sqrt(4j+1): + chi_4(d) at j = (d^2 - 1)/4 (4j+1 = d^2)
   and + 2 chi_4(d) at the later j = (d^2 - 1)/4 + i d.  Segmented, the
   segments in parallel; t2 < 2^16 for 4j + 1 < 2^32 (2 bytes per entry). */
#define T2_SEG 65536

typedef struct
{
    uint16_t * t2;
    ulong K;
}
t2_ctx;

static void
t2_worker(slong s, void * vctx)
{
    t2_ctx * T = vctx;
    ulong J0 = (ulong) s * T2_SEG, J1 = FLINT_MIN(T->K, J0 + T2_SEG), d, j;
    int32_t * buf;

    if (J0 >= J1)
        return;

    buf = flint_calloc(J1 - J0, sizeof(int32_t));

    for (d = 1; ; d += 2)
    {
        ulong jsq = (d * d - 1) / 4;
        int32_t chi = (d % 4 == 1) ? 1 : -1;

        if (jsq >= J1)
            break;

        if (jsq >= J0)
        {
            buf[jsq - J0] += chi;
            j = jsq + d;
        }
        else
            j = jsq + ((J0 - jsq + d - 1) / d) * d;

        for ( ; j < J1; j += d)
            buf[j - J0] += 2 * chi;
    }

    for (j = J0; j < J1; j++)
    {
        if (buf[j - J0] < 0 || buf[j - J0] > 65535)
            flint_throw(FLINT_ERROR, "t2 out of range\n");
        T->t2[j] = (uint16_t) buf[j - J0];
    }

    flint_free(buf);
}

uint16_t *
_acb_mzeta_t2_table(ulong K, slong nthreads)
{
    t2_ctx T;
    T.K = K;
    T.t2 = flint_malloc(sizeof(uint16_t) * FLINT_MAX(K, 1));
    flint_parallel_do(t2_worker, &T, (K + T2_SEG - 1) / T2_SEG, nthreads, FLINT_PARALLEL_STRIDED);
    return T.t2;
}

/* r3(8k + 3) = 8 t3(k), k0 <= k < k0 + len, with
   t3(k) = sum_{a : a(a+1)/2 <= k} t2(k - a(a+1)/2): about sqrt(2k) passes
   of len contiguous additions */
void
_acb_mzeta_r3_3mod8_block(uint32_t * out, const uint16_t * t2, ulong k0, ulong len)
{
    ulong a, s, k1 = k0 + len, k;

    memset(out, 0, len * sizeof(uint32_t));
    for (a = 0; (s = a * (a + 1) / 2) < k1; a++)
    {
        ulong lo = FLINT_MAX(k0, s), n = k1 - lo;
        const uint16_t * src = t2 + (lo - s);
        uint32_t * dst = out + (lo - k0);
        for (k = 0; k < n; k++)
            dst[k] += src[k];
    }
    for (k = 0; k < len; k++)
        out[k] *= 8;
}
