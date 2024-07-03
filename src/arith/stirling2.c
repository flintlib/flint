/*
    Copyright (C) 2010, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "arith.h"

#ifdef __GNUC__
# define log __builtin_log
# define exp __builtin_exp
# define frexp __builtin_frexp
#else
# include <math.h>
#endif

/* S(n,k) <= (1/2) binomial(n,k) * k^(n-k) */
static slong
stirling_2_bound_2exp(ulong n, ulong k)
{
    double bnk;
    int exp;
    slong bnk_exp, j;

    /* binomial coefficients */
    bnk = 1.0;
    bnk_exp = 0;

    for (j = 1; j <= k; j++)
    {
        bnk = (bnk * (n + 1 - j)) / j;
        bnk = frexp(bnk, &exp);
        bnk_exp += exp;
    }

    /* (1/2) * binomial(n,k) * k^(n-k) */
    /* add instead of subtract 1 to ensure upper bound in fp */
    return bnk_exp + (n - k) * log(k) * 1.44269504088896 + 1;
}

static slong
stirling_2_bound_2exp_vec(slong * bound, ulong n, slong len)
{
    double bnk;
    int exp;
    slong bnk_exp, k, max, kmax;

    kmax = len - 1;
    bound[0] = 0;

    max = 0;

    /* binomial coefficients */
    bnk = 1.0;
    bnk_exp = 0;

    for (k = 1; k <= kmax; k++)
    {
        /* binomial recurrence */
        bnk = (bnk * (n + 1 - k)) / k;
        bnk = frexp(bnk, &exp);
        bnk_exp += exp;

        /* (1/2) * binomial(n,k) * k^(n-k) */
        bound[k] = bnk_exp + (n-k) * log(k) * 1.44269504088896 + 1;
        max = FLINT_MAX(max, bound[k]);
    }

    return max;
}

#if FLINT_BITS == 32
#define MAX_N_1LIMB 16
#define MAX_N_2LIMB 26
#else
#define MAX_N_1LIMB 26
#define MAX_N_2LIMB 43
#endif

static void
triangular_1(nn_ptr c, slong n, slong klen)
{
    slong m, k;

    c[0] = 0; c[1] = 1; c[2] = 3; c[3] = 1;

    for (m = 4; m <= n; m++)
    {
        if (m < klen)
            c[m] = 1;

        for (k = FLINT_MIN(m, klen) - 1; k >= 2; k--)
            c[k] = c[k] * k + c[k - 1];
    }
}

static void
triangular_2(nn_ptr c, slong n, slong klen)
{
    ulong hi, lo;
    slong m, k;

    triangular_1(c, MAX_N_1LIMB, klen);

    for (k = FLINT_MIN(MAX_N_1LIMB, klen - 1); k >= 0; k--)
    {
        c[2 * k] = c[k];
        c[2 * k + 1] = 0;
    }

    for (m = MAX_N_1LIMB + 1; m <= n; m++)
    {
        if (m < klen)
        {
            c[2 * m] = 1;
            c[2 * m + 1] = 0;
        }

        for (k = FLINT_MIN(m, klen) - 1; k >= 2; k--)
        {
            umul_ppmm(hi, lo, c[2 * k], k);
            hi += c[2 * k + 1] * k;
            add_ssaaaa(c[2 * k + 1], c[2 * k],
                hi, lo, c[2 * (k - 1) + 1], c[2 * (k - 1)]);
        }
    }
}

void
arith_stirling_number_2_vec_triangular(fmpz * row, slong n, slong klen)
{
    ulong c[2 * MAX_N_2LIMB + 2];
    slong m, k;

    if (klen <= 0)
        return;

    if (n >= 1)
    {
        if (n <= MAX_N_1LIMB)
        {
            triangular_1(c, n, klen);
            for (k = 0; k <= FLINT_MIN(klen - 1, n); k++)
                fmpz_set_ui(row + k, c[k]);
        }
        else
        {
            m = FLINT_MIN(n, MAX_N_2LIMB);

            triangular_2(c, m, klen);
            for (k = 0; k <= FLINT_MIN(klen - 1, m); k++)
                fmpz_set_uiui(row + k, c[2 * k + 1], c[2 * k]);

            for (m = MAX_N_2LIMB + 1 ; m <= n; m++)
            {
                if (m < klen)
                    fmpz_one(row + m);

                for (k = FLINT_MIN(m, klen) - 1; k >= 2; k--)
                {
                    fmpz_mul_ui(row + k, row + k, k);
                    fmpz_add(row + k, row + k - 1, row + k);
                }
            }
        }
    }

    for (k = n; k < klen; k++)
        fmpz_set_ui(row + k, k == n);
}

void
arith_stirling_number_2_vec_convolution(fmpz * res, ulong n, slong klen)
{
    slong k, kodd, len;
    ulong e;
    fmpz *t, *u, *v;

    if (klen <= 0)
        return;

    len = FLINT_MIN(klen - 1, n - 1);

    t = _fmpz_vec_init(len + 1);
    u = _fmpz_vec_init(len);
    v = _fmpz_vec_init(len);

    if (n >= 1 && len >= 1)
    {
        fmpz_one(t + len);
        for (k = len - 1; k >= 0; k--)
            fmpz_mul_ui(t + k, t + k + 1, k + 1);

        for (kodd = 1; kodd <= len; kodd += 2)
        {
            fmpz_set_ui(v, kodd);
            fmpz_pow_ui(v, v, n);

            for (k = kodd, e = 0; k <= len; k *= 2, e++)
            {
                /* k^n / k! */
                fmpz_mul(u + k - 1, v, t + k);
                fmpz_mul_2exp(u + k - 1, u + k - 1, e * n);
            }
        }

        for (k = 1; k < len; k += 2)
            fmpz_neg(t + k, t + k);

        _fmpz_poly_mullow(v, t, len, u, len, len);

        fmpz_mul(t, t, t);
        for (k = 0; k < len; k++)
            fmpz_divexact(res + k + 1, v + k, t);
    }

    fmpz_set_ui(res + 0, n == 0);
    for (k = n; k < klen; k++)
        fmpz_set_ui(res + k, n == k);

    _fmpz_vec_clear(t, len + 1);
    _fmpz_vec_clear(u, len);
    _fmpz_vec_clear(v, len);
}

static void
divisor_table(unsigned int * tab, slong len)
{
    slong i, j;

    for (i = 0; i < len; i++)
    {
        tab[2 * i] = 1;
        tab[2 * i + 1] = i;
    }

    for (i = 2; i < len; i++)
    {
        for (j = 2; j <= i && i * j < len; j++)
        {
            tab[2 * i * j] = j;
            tab[2 * i * j + 1] = i;
        }
    }
}

static void
arith_stirling_number_2_nmod_vec(nn_ptr res, const unsigned int * divtab, ulong n, slong len, nmod_t mod)
{
    nn_ptr t, u;
    slong i;
    ulong c;
    TMP_INIT;
    TMP_START;

    t = TMP_ALLOC(len * sizeof(ulong));
    u = TMP_ALLOC(len * sizeof(ulong));

    /* compute inverse factorials */
    t[len - 1] = 1;
    for (i = len - 2; i >= 0; i--)
        t[i] = _nmod_mul_fullword(t[i + 1], i + 1, mod);

    c = nmod_inv(t[0], mod);
    t[0] = 1;
    for (i = 1; i < len; i++)
        t[i] = _nmod_mul_fullword(t[i], c, mod);

    /* compute powers */
    u[0] = nmod_pow_ui(0, n, mod);
    u[1] = nmod_pow_ui(1, n, mod);

    for (i = 2; i < len; i++)
    {
        if (divtab[2 * i] == 1)
            u[i] = nmod_pow_ui(i, n, mod);
        else
            u[i] = _nmod_mul_fullword(u[divtab[2 * i]], u[divtab[2 * i + 1]], mod);
    }

    for (i = 1; i < len; i++)
        u[i] = _nmod_mul_fullword(u[i], t[i], mod);

    for (i = 1; i < len; i += 2)
        t[i] = nmod_neg(t[i], mod);

    _nmod_poly_mullow(res, t, len, u, len, len, mod);

    TMP_END;
}

#define CRT_MAX_RESOLUTION 16

void
arith_stirling_number_2_vec_multi_mod(fmpz * res, ulong n, slong klen)
{
    fmpz_comb_t comb[CRT_MAX_RESOLUTION];
    fmpz_comb_temp_t temp[CRT_MAX_RESOLUTION];
    nn_ptr primes, residues;
    nn_ptr * polys;
    nmod_t mod;
    slong i, j, k, len, num_primes, num_primes_k, resolution;
    slong need_bits, size, prime_bits;
    slong *bounds;
    unsigned int * divtab;
    /* per comb */
    slong * local_len;
    slong * local_num_primes;

    if (klen <= 0)
        return;

    if (n <= 2)
    {
        arith_stirling_number_2_vec_triangular(res, n, klen);
        return;
    }

    if (klen > n + 1)
    {
        _fmpz_vec_zero(res + n + 1, klen - n - 1);
        klen = n + 1;
    }

    len = klen;

    bounds = flint_malloc(sizeof(slong) * len);

    need_bits = stirling_2_bound_2exp_vec(bounds, n, len);
    need_bits = FLINT_MAX(need_bits, 1);

    /* make bounds nonincreasing */
    for (k = len - 2; k >= 0; k--)
        bounds[k] = FLINT_MAX(bounds[k], bounds[k + 1]);

    resolution = FLINT_MAX(1, FLINT_MIN(CRT_MAX_RESOLUTION, n / 16));

    size = need_bits;
    prime_bits = FLINT_BITS - 1;
    num_primes = (size + prime_bits - 1) / prime_bits;

    primes = flint_malloc(num_primes * sizeof(ulong));
    residues = flint_malloc(num_primes * sizeof(ulong));
    polys = flint_malloc(num_primes * sizeof(nn_ptr));
    divtab = flint_malloc(2 * len * sizeof(unsigned int));

    divisor_table(divtab, len);

    local_len = flint_malloc(resolution * sizeof(slong));
    local_num_primes = flint_malloc(resolution * sizeof(slong));

    primes[0] = n_nextprime(UWORD(1) << prime_bits, 0);
    for (k = 1; k < num_primes; k++)
        primes[k] = n_nextprime(primes[k-1], 0);

    for (i = 0; i < resolution; i++)
    {
        local_num_primes[i] = FLINT_MAX(1, num_primes * (i + 1) / resolution);

        fmpz_comb_init(comb[i], primes, local_num_primes[i]);
        fmpz_comb_temp_init(temp[i], comb[i]);
        local_len[i] = len;

        if (i > 0)
        {
            while (local_len[i] > 0 && bounds[local_len[i] - 1] < prime_bits * local_num_primes[i - 1])
                local_len[i]--;
        }
    }

    for (j = 0; j < num_primes; j++)
    {
        i = resolution - 1;
        while (i > 0 && j < local_num_primes[i - 1])
            i--;

        polys[j] = _nmod_vec_init(local_len[i]);
        nmod_init(&mod, primes[j]);
        arith_stirling_number_2_nmod_vec(polys[j], divtab, n, local_len[i], mod);
    }

    for (k = 0; k < len; k++)
    {
        i = 0;
        while (i + 1 < resolution && bounds[k] >= comb[i]->num_primes * prime_bits)
            i++;

        /* Use only as large a comb as needed */
        num_primes_k = comb[i]->num_primes;

        for (j = 0; j < num_primes_k; j++)
            residues[j] = polys[j][k];

        fmpz_multi_CRT_ui(res + k, residues, comb[i], temp[i], 0);
    }

    /* Cleanup */
    for (k = 0; k < num_primes; k++)
        _nmod_vec_clear(polys[k]);

    for (i = 0; i < resolution; i++)
    {
        fmpz_comb_temp_clear(temp[i]);
        fmpz_comb_clear(comb[i]);
    }

    flint_free(primes);
    flint_free(residues);
    flint_free(polys);
    flint_free(divtab);

    flint_free(bounds);

    flint_free(local_len);
    flint_free(local_num_primes);
}

void
arith_stirling_number_2_vec(fmpz * row, ulong n, slong klen)
{
    if (n <= 80)
        arith_stirling_number_2_vec_triangular(row, n, klen);
    else if (klen < n / 2)
        arith_stirling_number_2_vec_convolution(row, n, klen);
    else
        arith_stirling_number_2_vec_multi_mod(row, n, klen);
}

static void
stirling_2_egf(fmpz_t res, ulong n, ulong k)
{
    fmpz * num;
    fmpz * rnum;
    fmpz_t den;
    fmpz_t rden;
    slong i, len;

    if (k >= n || k == 0)
    {
        fmpz_set_ui(res, n == k);
        return;
    }

    len = n - k + 1;

    num = _fmpz_vec_init(len);
    rnum = _fmpz_vec_init(len);
    fmpz_init(den);
    fmpz_init(rden);

    fmpz_one(num + len - 1);
    for (i = len - 2; i >= 0; i--)
        fmpz_mul_ui(num + i, num + i + 1, i + 2);

    fmpz_set(den, num + 0);

    _fmpq_poly_pow_trunc(rnum, rden, num, den, len, k, len);

    fmpz_set_ui(num, k);
    fmpz_add_ui(num, num, 1);
    fmpz_rfac_ui(num, num, n - k);

    fmpz_mul(num, num, rnum + n - k);
    fmpz_divexact(res, num, rden);

    _fmpz_vec_clear(num, len);
    _fmpz_vec_clear(rnum, len);
    fmpz_clear(den);
    fmpz_clear(rden);
}

static void
stirling_2_powsum(fmpz_t s, ulong n, ulong k)
{
    fmpz_t t, u;
    fmpz *b;
    slong i, j, m, max_b;

    max_b = (k + 1) / 2;

    fmpz_init(t);
    fmpz_init(u);
    b = _fmpz_vec_init(max_b + 1);

    fmpz_one(b + 0);
    for (j = 1; j <= max_b; j++)
    {
        fmpz_mul_ui(b + j, b + j - 1, k + 1 - j);
        fmpz_divexact_ui(b  + j, b + j, j);
    }

    fmpz_zero(s);
    for (j = 1; j <= k; j += 2)
    {
        fmpz_ui_pow_ui(u, j, n);

        m = j;
        while (1)  /* Process each m = 2^p * j */
        {
            i = (m <= max_b) ? m : k - m;

            if ((k + m) & 1)
                fmpz_submul(s, b + i, u);
            else
                fmpz_addmul(s, b + i, u);

            m *= 2;
            if (m > k)
                break;
            fmpz_mul_2exp(u, u, n);
        }
    }

    _fmpz_vec_clear(b, max_b + 1);
    fmpz_fac_ui(t, k);
    fmpz_divexact(s, s, t);
    fmpz_clear(t);
    fmpz_clear(u);
}

/* req: k >= 2 */
static ulong
stirling_2_nmod(const unsigned int * divtab, ulong n, ulong k, nmod_t mod)
{
    nn_ptr t, u;
    slong i, bin_len, pow_len;
    ulong s1, s2, bden, bd;
    dot_params_t params;
    TMP_INIT;
    TMP_START;

    pow_len = k + 1;
    bin_len = FLINT_MIN(pow_len, k / 2 + 1);

    t = TMP_ALLOC(bin_len * sizeof(ulong));
    u = TMP_ALLOC(pow_len * sizeof(ulong));

    /* compute binomial coefficients + denominator */
    t[0] = 1;
    for (i = 1; i < bin_len; i++)
        t[i] = _nmod_mul_fullword(t[i - 1], k + 1 - i, mod);

    bd = t[bin_len - 1 - (k + 1) % 2];
    bden = 1;
    for (i = bin_len - 1; i >= 0; i--)
    {
        bden = _nmod_mul_fullword(bden, i + 1, mod);
        t[i] = _nmod_mul_fullword(t[i], bden, mod);
    }

    /* compute powers */
    u[0] = nmod_pow_ui(0, n, mod);
    u[1] = nmod_pow_ui(1, n, mod);

    for (i = 2; i < pow_len; i++)
    {
        if (divtab[2 * i] == 1)
            u[i] = nmod_pow_ui(i, n, mod);
        else
            u[i] = _nmod_mul_fullword(u[divtab[2 * i]], u[divtab[2 * i + 1]], mod);
    }

    for (i = 1; i < bin_len; i += 2)
        t[i] = nmod_neg(t[i], mod);

    params = _nmod_vec_dot_params(bin_len, mod);
    s1 = _nmod_vec_dot(t, u, bin_len, mod, params);

    if (pow_len > bin_len)
    {
        params = _nmod_vec_dot_params(pow_len - bin_len, mod);
        s2 = _nmod_vec_dot_rev(u + bin_len, t + k - pow_len + 1, pow_len - bin_len, mod, params);
        if (k % 2)
            s1 = nmod_sub(s1, s2, mod);
        else
            s1 = nmod_add(s1, s2, mod);
    }

    TMP_END;

    if (k % 2)
        s1 = nmod_neg(s1, mod);

    bden = nmod_mul(bden, bden, mod);
    bden = nmod_mul(bden, bd, mod);
    bden = nmod_inv(bden, mod);
    return nmod_mul(s1, bden, mod);
}

static void
_fmpz_crt_combine(fmpz_t r1r2, fmpz_t m1m2, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2, const fmpz_t m2)
{
    fmpz_invmod(m1m2, m1, m2);
    fmpz_mul(m1m2, m1m2, m1);
    fmpz_sub(r1r2, r2, r1);
    fmpz_mul(r1r2, r1r2, m1m2);
    fmpz_add(r1r2, r1r2, r1);
    fmpz_mul(m1m2, m1, m2);
    fmpz_mod(r1r2, r1r2, m1m2);
}

static void
tree_crt(fmpz_t r, fmpz_t m, nn_srcptr residues, nn_srcptr primes, slong len)
{
    if (len == 0)
    {
        fmpz_zero(r);
        fmpz_one(m);
    }
    else if (len == 1)
    {
        fmpz_set_ui(r, residues[0]);
        fmpz_set_ui(m, primes[0]);
    }
    else
    {
        fmpz_t r1, m1, r2, m2;

        fmpz_init(r1);
        fmpz_init(m1);
        fmpz_init(r2);
        fmpz_init(m2);

        tree_crt(r1, m1, residues, primes, len / 2);
        tree_crt(r2, m2, residues + len / 2, primes + len / 2, len - len / 2);
        _fmpz_crt_combine(r, m, r1, m1, r2, m2);

        fmpz_clear(r1);
        fmpz_clear(m1);
        fmpz_clear(r2);
        fmpz_clear(m2);
    }
}

static void
stirling_2_multi_mod(fmpz_t res, ulong n, ulong k)
{
    fmpz_t tmp;
    nmod_t mod;
    nn_ptr primes, residues;
    slong i, num_primes;
    flint_bitcnt_t size, prime_bits;
    unsigned int * divtab;

    size = stirling_2_bound_2exp(n, k);
    prime_bits = FLINT_BITS - 1;
    num_primes = (size + prime_bits - 1) / prime_bits;

    fmpz_init(tmp);
    primes = flint_malloc(num_primes * sizeof(ulong));
    residues = flint_malloc(num_primes * sizeof(ulong));

    divtab = flint_malloc(2 * sizeof(unsigned int) * (n + 1));
    divisor_table(divtab, n + 1);

    primes[0] = n_nextprime(UWORD(1) << prime_bits, 0);
    for (i = 1; i < num_primes; i++)
        primes[i] = n_nextprime(primes[i - 1], 0);

    for (i = 0; i < num_primes; i++)
    {
        nmod_init(&mod, primes[i]);
        residues[i] = stirling_2_nmod(divtab, n, k, mod);
    }

    tree_crt(res, tmp, residues, primes, num_primes);

    flint_free(primes);
    flint_free(residues);
    flint_free(divtab);
    fmpz_clear(tmp);
}

void
arith_stirling_number_2(fmpz_t res, ulong n, ulong k)
{
    if (k >= n)
    {
        fmpz_set_ui(res, n == k);
    }
    else if (k <= 1)
    {
        fmpz_set_ui(res, k);
    }
    else if (k == n - 1)  /* S(n, n-1) = binomial(n, 2) */
    {
        fmpz_set_ui(res, n);
        fmpz_mul_ui(res, res, n - 1);
        fmpz_tdiv_q_2exp(res, res, 1);
    }
    else if (k == 2) /* S(n,2) = 2^(n-1)-1 */
    {
        fmpz_one(res);
        fmpz_mul_2exp(res, res, n - 1);
        fmpz_sub_ui(res, res, 1);
    }
    else if (n <= MAX_N_1LIMB)
    {
        ulong c[MAX_N_2LIMB + 1];
        triangular_1(c, n, k + 1);
        fmpz_set_ui(res, c[k]);
    }
    else if (n <= MAX_N_2LIMB)
    {
        ulong c[2 * MAX_N_2LIMB + 2];
        triangular_2(c, n, k + 1);
        fmpz_set_uiui(res, c[2 * k + 1], c[2 * k]);
    }
    else
    {
        double low_cutoff, high_cutoff;

        if (n < 200)
        {
            low_cutoff = high_cutoff = 0.9;
        }
        else
        {
            if (n < 3000)
                low_cutoff = 0.95 * exp(-0.00022 * n);
            else
                low_cutoff = 1500 / n;

            low_cutoff = FLINT_MAX(low_cutoff, 0.0002);
            low_cutoff = FLINT_MIN(low_cutoff, 0.8);

            high_cutoff = 0.92 + 0.005 * log(n);
            high_cutoff = FLINT_MIN(high_cutoff, 0.98);
        }

        if (k <= low_cutoff * n)
            stirling_2_powsum(res, n, k);
        else if (k >= high_cutoff * n)
            stirling_2_egf(res, n, k);
        else
            stirling_2_multi_mod(res, n, k);
    }
}
