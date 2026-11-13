/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "thread_support.h"
#include "acb.h"
#include "impl.h"

/* Multiple polylogarithms through Goncharov's iterated integrals

     G(a_1, ..., a_L; z) = int_{z > t_1 > ... > t_L > 0} dt_1/(t_1 - a_1) ... dt_L/(t_L - a_L)

   (along the segment [0, z]; a_L != 0 for convergence at 0), with

     Li_{s_1,...,s_r}(x_1, ..., x_r) = sum_{n_1 > ... > n_r >= 1} x_1^n_1 ... x_r^n_r / (n_1^s_1 ... n_r^s_r)
         = (-1)^r G(0^(s_1-1), 1/y_1, ..., 0^(s_r-1), 1/y_r; 1),  y_i = x_1 ... x_i.

   Write the word in blocks of s_i - 1 zeros followed by b_i != 0,
   for i = 1, ..., d. For a suffix w beginning with e - 1 zeros followed
   by b_i, put r_j = z/b_j and n_{d+1} = 0. For |z| < min_j |b_j|,
   the series at zero is

     G(w; z) = (-1)^(d-i+1) sum_{n_i > ... > n_d >= 1}
                 r_i^(n_i-n_{i+1}) / n_i^e
                 prod_{j=i+1}^d (r_j^(n_j-n_{j+1}) / n_j^(s_j)).

   The degree is n_i. Grouping terms by n_i = n gives the partial sum
   (-1)^(d-i+1) sum_{n=1}^N C_i(n)/n^e, where

     C_d(n) = r_d^n,
     C_i(n) = r_i (C_i(n-1) + C_{i+1}(n-1)/(n-1)^s_{i+1})
                  (i < d, n >= 2).

   For i < d, start with C_i(1) = 0. Suffixes starting in the same
   block share C_i(n). See g_suffixes_scaled for integer scaling and
   g_suffixes_depth_one for the blocked Horner evaluation at depth one.
   See g_series_tail for the common truncation bound. */

typedef struct
{
    slong L, depth, max_index;
    slong * state_index;   /* index of a nonzero letter in reverse order */
    slong * block_end;     /* position of the last letter in each block */
    acb_ptr ratios;        /* z / a_l for nonzero letters */
}
g_series_ctx;

/* Propagate complex states as disks to avoid wrapping under rotation.
   The radius bounds the Euclidean distance to the midpoint. */
typedef struct
{
    acb_t ratio_midpoint, midpoint, input, value;
    mag_t ratio_abs, ratio_error, error, scratch;
}
g_disk_struct;

static void
g_disk_init(g_disk_struct * disk, const acb_t ratio)
{
    acb_init(disk->ratio_midpoint);
    acb_init(disk->midpoint);
    acb_init(disk->input);
    acb_init(disk->value);
    mag_init(disk->ratio_abs);
    mag_init(disk->ratio_error);
    mag_init(disk->error);
    mag_init(disk->scratch);

    acb_get_mid(disk->ratio_midpoint, ratio);
    acb_get_mag(disk->ratio_abs, ratio);
    mag_hypot(disk->ratio_error, arb_radref(acb_realref(ratio)),
        arb_radref(acb_imagref(ratio)));
}

static void
g_disk_clear(g_disk_struct * disk)
{
    acb_clear(disk->ratio_midpoint);
    acb_clear(disk->midpoint);
    acb_clear(disk->input);
    acb_clear(disk->value);
    mag_clear(disk->ratio_abs);
    mag_clear(disk->ratio_error);
    mag_clear(disk->error);
    mag_clear(disk->scratch);
}

static void
g_disk_update_midpoint(g_disk_struct * disk)
{
    mag_hypot(disk->scratch, arb_radref(acb_realref(disk->value)),
        arb_radref(acb_imagref(disk->value)));
    mag_add(disk->error, disk->error, disk->scratch);
    acb_get_mid(disk->midpoint, disk->value);
}

static void
g_disk_get_value(acb_t res, const g_disk_struct * disk)
{
    acb_set(res, disk->midpoint);
    acb_add_error_mag(res, disk->error);
}

static void
g_disk_multiply(g_disk_struct * disk, const acb_t input, slong prec)
{
    /* |r s - r0 s0| <= |r| |s - s0| + |s0| |r - r0|. */
    mag_mul(disk->error, disk->error, disk->ratio_abs);
    if (!mag_is_zero(disk->ratio_error))
    {
        acb_get_mag(disk->scratch, input);
        mag_addmul(disk->error, disk->scratch, disk->ratio_error);
    }
    acb_mul(disk->value, disk->ratio_midpoint, input, prec);
}

static void
g_disk_scaled_step(g_disk_struct * disk,
    const g_disk_struct * next, const fmpz_t scale, const fmpz_t next_scale,
    slong prec)
{
    acb_mul_fmpz(disk->input, disk->midpoint, scale, prec);
    mag_mul_fmpz(disk->error, disk->error, scale);
    if (next != NULL)
    {
        if (fmpz_is_one(next_scale))
        {
            acb_add(disk->input, disk->input, next->midpoint, prec);
            mag_add(disk->error, disk->error, next->error);
        }
        else
        {
            acb_addmul_fmpz(disk->input, next->midpoint, next_scale, prec);
            mag_mul_fmpz(disk->scratch, next->error, next_scale);
            mag_add(disk->error, disk->error, disk->scratch);
        }
    }
    mag_hypot(disk->scratch, arb_radref(acb_realref(disk->input)),
        arb_radref(acb_imagref(disk->input)));
    mag_add(disk->error, disk->error, disk->scratch);
    acb_get_mid(disk->input, disk->input);
    g_disk_multiply(disk, disk->input, prec);
    g_disk_update_midpoint(disk);
    g_disk_get_value(disk->value, disk);
}

/* With C = max_i s_i, let U_i(n) = ((n-1)!)^C C_i(n) and let B_{i,e}(n)
   be the unsigned partial sum through degree n scaled by (n!)^C. Then

     U_i(n+1) = r_i (n^C U_i(n) + n^(C-s_{i+1}) U_{i+1}(n)),
     B_{i,e}(n) = n^C B_{i,e}(n-1) + n^(C-e) U_i(n).

   The innermost update has only the first term, with U_d(1) = r_d.
   Initialize U_i(1) = 0 for i < d and B_{i,e}(0) = 0.
   In a block [lo, hi], scale B_{i,e}(lo-1) by (hi!/(lo-1)!)^C
   once, then add each n^(C-e) U_i(n) weighted by (hi!/n!)^C.
   Divide by (N!)^C and apply the suffix sign at the end.
   Update from the outermost block inwards so the next state still has
   its degree-n value when used by the current block. */
static void
g_suffixes_scaled(acb_ptr res, const g_series_ctx * G, slong N, slong wp,
    int is_real, double qd)
{
    enum
    {
        BLOCK_SIZE = 4
    };
    const slong min_prec = FLINT_MIN(64, wp);
    const double decay = qd == 0.0 ? 0.0 : -log2(qd);
    const double guard = ceil((G->depth - 1) * log2(1.0 + log((double) N)))
        + FLINT_BIT_COUNT(N);
    slong n, lo, hi, k, start, last, next, j, state;
    acb_ptr real_terms = NULL;
    g_disk_struct * disks = NULL;
    arb_t denominator;
    fmpz powers[BLOCK_SIZE] = {0}, products[BLOCK_SIZE + 1] = {0};
    fmpz_t factor;

    arb_init(denominator);
    fmpz_init(factor);
    if (is_real)
    {
        real_terms = _acb_vec_init(G->depth);
        acb_set(real_terms, G->ratios + G->L - 1);
    }
    else
    {
        disks = flint_malloc(G->depth * sizeof(g_disk_struct));
        for (j = 0; j < G->L; j++)
            if (G->state_index[j] >= 0)
                g_disk_init(disks + G->state_index[j], G->ratios + j);
        acb_one(disks[0].midpoint);
        g_disk_multiply(disks, disks[0].midpoint, wp);
        g_disk_update_midpoint(disks);
        g_disk_get_value(disks[0].value, disks);
    }
    _acb_vec_zero(res, G->L);
    acb_one(res + G->L);

    lo = 1;
    hi = 0;
    for (n = 1; n <= N; n++)
    {
        const fmpz * power, * remaining_product;
        /* Unscaled states are bounded by q^n (1 + log N)^(d-1).
           Allow for this factor and the sum of rounding errors before
           reducing state precision. Accumulate all sums at wp. */
        double drop = FLINT_MAX(0.0, floor(n * decay - guard));
        slong state_prec = wp - (slong) FLINT_MIN(wp - min_prec, drop);

        if (n > hi)
        {
            lo = n;
            hi += FLINT_MIN(BLOCK_SIZE, N - hi);
            /* products[k] scales the sum from degree lo+k-1 to hi. */
            fmpz_one(products + hi - lo + 1);
            for (k = hi - lo; k >= 0; k--)
            {
                fmpz_ui_pow_ui(powers + k, lo + k, G->max_index);
                fmpz_mul(products + k, products + k + 1, powers + k);
            }
            for (j = 0; j < G->L; j++)
                acb_mul_fmpz(res + j, res + j, products, wp);
        }
        power = powers + n - lo;
        remaining_product = products + n - lo + 1;

        for (start = 0; start < G->L; start = last + 1)
        {
            acb_srcptr term;

            last = G->block_end[start];
            state = G->state_index[last];
            term = is_real ? real_terms + state : disks[state].value;

            fmpz_ui_pow_ui(factor, n, G->max_index - (last - start + 1));
            fmpz_mul(factor, factor, remaining_product);
            for (j = start; j <= last; j++)
            {
                if (fmpz_is_one(factor))
                    acb_add(res + j, res + j, term, wp);
                else
                    acb_addmul_fmpz(res + j, term, factor, wp);
                if (j < last)
                    fmpz_mul_ui(factor, factor, n);
            }

            if (n < N)
            {
                if (state > 0)
                {
                    next = G->block_end[last + 1];
                    fmpz_ui_pow_ui(factor, n, G->max_index - (next - last));
                }
                if (is_real)
                {
                    acb_ptr current = real_terms + state;

                    acb_mul_fmpz(current, current, power, state_prec);
                    if (state > 0)
                    {
                        if (fmpz_is_one(factor))
                            acb_add(current, current, current - 1, state_prec);
                        else
                            acb_addmul_fmpz(current, current - 1, factor,
                                state_prec);
                    }
                    acb_mul(current, current, G->ratios + last, state_prec);
                }
                else
                {
                    g_disk_scaled_step(disks + state,
                        state > 0 ? disks + state - 1 : NULL,
                        power, factor, state_prec);
                }
            }
        }
    }

    arb_fac_ui(denominator, N, wp);
    arb_pow_ui(denominator, denominator, G->max_index, wp);
    for (j = 0; j < G->L; j++)
    {
        acb_div_arb(res + j, res + j, denominator, wp);
        if ((G->state_index[G->block_end[j]] + 1) & 1)
            acb_neg(res + j, res + j);
    }

    fmpz_clear(factor);
    for (k = 0; k < BLOCK_SIZE; k++)
        fmpz_clear(powers + k);
    for (k = 0; k <= BLOCK_SIZE; k++)
        fmpz_clear(products + k);
    arb_clear(denominator);
    if (is_real)
        _acb_vec_clear(real_terms, G->depth);
    else
    {
        for (state = 0; state < G->depth; state++)
            g_disk_clear(disks + state);
        flint_free(disks);
    }
}

/* For a block [lo, hi], set D = (lo * ... * hi)^m and evaluate
   sum_{n=lo}^{hi} (D/n^m) r^(n-lo) by Horner's method. Combine blocks
   in decreasing degree order, multiplying the accumulator by
   r^(hi-lo+1) before adding the block divided by D. Multiply by -r
   at the end to account for the starting degree and suffix sign. */
static void
g_suffixes_depth_one(acb_ptr res, const g_series_ctx * G, slong N, slong wp,
    int is_real)
{
    const slong block_size = 32;
    acb_srcptr ratio = G->ratios + G->L - 1;
    slong hi, lo, width, j, m, n;
    fmpz_t product, den, coeff;
    acb_t block, block_power;
    arb_t den_ball;
    g_disk_struct disk;

    if (!is_real)
        g_disk_init(&disk, ratio);
    fmpz_init(product);
    fmpz_init(den);
    fmpz_init(coeff);
    acb_init(block);
    acb_init(block_power);
    arb_init(den_ball);
    acb_pow_ui(block_power, ratio, block_size, wp);
    _acb_vec_zero(res, G->L);

    for (hi = N; hi > 0; hi -= block_size)
    {
        width = FLINT_MIN(hi, block_size);
        lo = hi - width + 1;
        if (width < block_size)
            acb_pow_ui(block_power, ratio, width, wp);
        fmpz_rfac_uiui(product, lo, width);
        fmpz_set(den, product);

        for (j = G->L - 1; j >= 0; j--)
        {
            m = G->L - j;
            if (is_real)
                acb_zero(block);
            else
            {
                acb_zero(disk.midpoint);
                mag_zero(disk.error);
            }
            for (n = hi; n >= lo; n--)
            {
                if (m == 1)
                    fmpz_divexact_ui(coeff, den, n);
                else
                {
                    fmpz_ui_pow_ui(coeff, n, m);
                    fmpz_divexact(coeff, den, coeff);
                }
                if (is_real)
                {
                    acb_mul(block, block, ratio, wp);
                    acb_add_fmpz(block, block, coeff, wp);
                }
                else
                {
                    g_disk_multiply(&disk, disk.midpoint, wp);
                    acb_add_fmpz(disk.value, disk.value, coeff, wp);
                    g_disk_update_midpoint(&disk);
                }
            }

            if (!is_real)
                g_disk_get_value(block, &disk);
            arb_set_fmpz(den_ball, den);
            acb_div_arb(block, block, den_ball, wp);
            acb_mul(res + j, res + j, block_power, wp);
            acb_add(res + j, res + j, block, wp);

            if (j > 0)
                fmpz_mul(den, den, product);
        }
    }

    for (j = 0; j < G->L; j++)
    {
        acb_mul(res + j, res + j, ratio, wp);
        acb_neg(res + j, res + j);
    }
    acb_one(res + G->L);

    fmpz_clear(product);
    fmpz_clear(den);
    fmpz_clear(coeff);
    acb_clear(block);
    acb_clear(block_power);
    arb_clear(den_ball);
    if (!is_real)
        g_disk_clear(&disk);
}

static int
g_suffixes_series(acb_ptr res, const g_series_ctx * G, slong N,
    const mag_t tail, slong wp, double qd)
{
    slong l;
    int is_real = _acb_vec_is_real(G->ratios, G->L);
    int ok = 1;

    if (G->depth == 1)
        g_suffixes_depth_one(res, G, N, wp, is_real);
    else
        g_suffixes_scaled(res, G, N, wp, is_real, qd);

    for (l = 0; l < G->L; l++)
    {
        acb_add_error_mag(res + l, tail);
        if (!acb_is_finite(res + l))
            ok = 0;
    }

    return ok;
}

/* For a suffix of depth k, the degree-n terms have absolute sum at most
   q^n H_(n-1)^(k-1) / (n (k-1)!), where q bounds all |z/b_j|. The gaps
   give the factor q^n, and the distinct inner indices give the factorial.
   Write M(H) = max_{0 <= j < d} H^j/j!. The maximum is attained at
   min(d - 1, floor(H)), by the ratio of successive terms. Consequently
   M(1 + log(n))/n is nonincreasing: on each interval with fixed maximizing
   j its logarithmic derivative with respect to log(n) is j/(1 + log(n)) - 1
   <= 0. The maximum is continuous when the maximizing index changes.
   At T = N + 1 the resulting common absolute remainder bound is

     q^T / (1 - q) * M(1 + log(T)) / T.

   Compute this with outward rounding. g_series_cutoff increases N
   until the bound is at most 2^(-prec-8). */
static void
g_series_tail(mag_t tail, const mag_t q, slong d, slong N)
{
    mag_t h, t;
    slong j;
    ulong n = N + 1;

    mag_init(h);
    mag_init(t);
    mag_log_ui(h, n);
    mag_add_ui(h, h, 1);
    /* In this range mag_get_d represents h exactly. The maximizer
       must be taken for the rounded upper bound h used below. */
    j = FLINT_MIN(d - 1, (slong) mag_get_d(h));
    mag_pow_ui(h, h, j);
    mag_rfac_ui(t, j);
    mag_mul(h, h, t);
    mag_div_ui(h, h, n);
    mag_geom_series(tail, q, n);
    mag_mul(tail, tail, h);
    mag_clear(h);
    mag_clear(t);
}

static slong
g_series_cutoff(mag_t tail, const mag_t q, double qd, slong d, slong prec)
{
    const slong min_n = 8;
    const slong max_n = WORD_MAX / 2;
    slong N, i, j;
    double rate, target, n, h;

    if (mag_is_zero(q))
    {
        mag_zero(tail);
        return min_n;
    }

    rate = -log2(qd);
    target = (double) prec + 8 - log2(1 - qd);
    n = FLINT_MAX((double) min_n, target / rate);

    /* Estimate the cutoff in double precision, then certify it in mag. */
    for (i = 0; i < 2; i++)
    {
        h = 1 + log(n);
        j = FLINT_MIN(d - 1, (slong) h);
        n = FLINT_MAX((double) min_n,
            (target + j * log2(h) - lgamma(j + 1.0) / log(2.0) - log2(n))
            / rate);
    }

    if (!(n < (double) max_n))
        return 0;

    N = (slong) ceil(n);
    g_series_tail(tail, q, d, N);
    while (mag_cmp_2exp_si(tail, -prec - 8) > 0)
    {
        if (N > max_n - (1 + N / 16))
            return 0;

        N += 1 + N / 16;
        g_series_tail(tail, q, d, N);
    }

    return N;
}

int
_acb_mzeta_G_suffixes(acb_ptr res, acb_srcptr a, slong L, const acb_t z,
    slong prec)
{
    g_series_ctx G;
    slong l, d, run, N, wp;
    mag_t R, q, tail, t;
    double qd;
    int ok = 1;

    if (L == 0)
    {
        acb_one(res);
        return 1;
    }
    if (acb_is_zero(a + L - 1))
        return 0;

    G.L = L;
    G.state_index = flint_malloc(sizeof(slong) * L);
    G.block_end = flint_malloc(sizeof(slong) * L);
    G.ratios = _acb_vec_init(L);

    mag_init(R);
    mag_init(q);
    mag_init(tail);
    mag_init(t);

    /* nonzero letters (as balls: a letter whose ball contains 0 is not
       allowed), runs of zeros, R = min |a_l| */
    d = 0;
    run = 0;
    G.max_index = 1;
    mag_inf(R);
    for (l = L - 1; l >= 0; l--)
    {
        if (acb_is_zero(a + l))
        {
            G.state_index[l] = -1;
            G.block_end[l] = G.block_end[l + 1];
            run++;
            G.max_index = FLINT_MAX(G.max_index, run + 1);
        }
        else
        {
            if (acb_contains_zero(a + l))
            {
                ok = 0;
                break;
            }
            G.state_index[l] = d++;
            G.block_end[l] = l;
            run = 0;
            acb_get_mag_lower(t, a + l);
            mag_min(R, R, t);
        }
    }
    G.depth = d;

    /* q = |z| / R < 1 */
    if (ok)
    {
        acb_get_mag(q, z);
        mag_div(q, q, R);
        qd = mag_get_d(q);
        ok = (qd < 0.95);
    }

    if (ok)
    {
        N = g_series_cutoff(tail, q, qd, d, prec);
        ok = (N != 0 && !mag_is_inf(tail));
    }

    if (ok)
    {
        wp = prec + 2 * FLINT_BIT_COUNT(N) + 2 * FLINT_BIT_COUNT(d + L) + 16;

        for (l = 0; l < L; l++)
            if (G.state_index[l] >= 0)
                acb_div(G.ratios + l, z, a + l, wp);

        ok = g_suffixes_series(res, &G, N, tail, wp, qd);
    }

    flint_free(G.state_index);
    flint_free(G.block_end);
    _acb_vec_clear(G.ratios, L);
    mag_clear(R);
    mag_clear(q);
    mag_clear(tail);
    mag_clear(t);
    return ok;
}

/* G(a_1, ..., a_L; 1): directly if q = 1/R <= 1/2, else split at 1/2:
   G(a; 1) = sum_j (-1)^j G(1 - a_j, ..., 1 - a_1; 1/2) G(a_{j+1}, ..., a_L; 1/2)
   (t -> 1 - t above 1/2), if the letters a and 1 - a are at distance
   >= 0.75 from 0 (ratio <= 2/3).  Returns 0 if not applicable. */
typedef struct
{
    acb_ptr res[2];
    acb_srcptr a[2];
    slong L, prec;
    int ok[2];
}
gjob;

static void
gworker(slong i, void * vj)
{
    gjob * J = vj;
    acb_t z;
    acb_init(z);
    acb_one(z);
    acb_mul_2exp_si(z, z, -1);
    J->ok[i] = _acb_mzeta_G_suffixes(J->res[i], J->a[i], J->L, z, J->prec);
    acb_clear(z);
}

int
_acb_mzeta_G_one(acb_t res, acb_srcptr a, slong L, slong prec)
{
    slong l, wp;
    acb_ptr b, A, B;
    mag_t t, R1, R2;
    int ok;

    if (L == 0)
    {
        acb_one(res);
        return 1;
    }
    if (acb_is_zero(a + L - 1))
        return 0;

    mag_init(t);
    mag_init(R1);
    mag_init(R2);
    b = _acb_vec_init(L);
    mag_inf(R1);
    mag_inf(R2);
    for (l = 0; l < L; l++)
    {
        acb_sub_ui(b + l, a + l, 1, prec + 16);
        acb_neg(b + l, b + l);          /* 1 - a_l */
        if (!acb_is_zero(a + l))
        {
            acb_get_mag_lower(t, a + l);
            mag_min(R1, R1, t);
        }
        if (!acb_is_zero(b + l))
        {
            acb_get_mag_lower(t, b + l);
            mag_min(R2, R2, t);
        }
    }

    wp = prec + FLINT_BIT_COUNT(L) + 8;

    if (mag_cmp_2exp_si(R1, 1) >= 0)
    {
        /* all letters at distance >= 2: series at 1 directly */
        acb_ptr r = _acb_vec_init(L + 1);
        acb_t z;
        acb_init(z);
        acb_one(z);
        ok = _acb_mzeta_G_suffixes(r, a, L, z, wp);
        if (ok)
            acb_set_round(res, r, prec);
        acb_clear(z);
        _acb_vec_clear(r, L + 1);
    }
    else if (mag_cmp_2exp_si(R1, -1) > 0 && mag_cmp_2exp_si(R2, -1) > 0
        && mag_get_d(R1) >= 0.75 && mag_get_d(R2) >= 0.75 && !acb_is_zero(b + 0))
    {
        /* the reversed reflected word for the upper part */
        acb_ptr ar = _acb_vec_init(L);
        gjob J;
        for (l = 0; l < L; l++)
            acb_set(ar + l, b + L - 1 - l);
        A = _acb_vec_init(L + 1);
        B = _acb_vec_init(L + 1);
        J.res[0] = A; J.res[1] = B;
        J.a[0] = ar; J.a[1] = a;
        J.L = L;
        J.prec = wp;
        flint_parallel_do(gworker, &J, 2, flint_get_num_threads(), FLINT_PARALLEL_UNIFORM);
        ok = J.ok[0] && J.ok[1];
        if (ok)
        {
            /* A[L - j] = G(1 - a_j, ..., 1 - a_1; 1/2) (j letters),
               B[j] = G(a_{j+1}, ..., a_L; 1/2) */
            /* (t -> 1 - t: dt/(t - a) = -ds/(s - (1 - a)), a sign per
               letter of the upper part) */
            acb_zero(res);
            for (l = 0; l <= L; l++)
            {
                if (l & 1)
                    acb_submul(res, A + L - l, B + l, wp);
                else
                    acb_addmul(res, A + L - l, B + l, wp);
            }
            acb_set_round(res, res, prec);
        }
        _acb_vec_clear(ar, L);
        _acb_vec_clear(A, L + 1);
        _acb_vec_clear(B, L + 1);
    }
    else
        ok = 0;

    _acb_vec_clear(b, L);
    mag_clear(t);
    mag_clear(R1);
    mag_clear(R2);
    return ok;
}

/* from the partial products y_i = x_1 ... x_i (exact when possible) */
static void
polylog_y(acb_t res, const ulong * s, acb_srcptr y, slong r, slong prec)
{
    slong i, k = 0, L;
    acb_ptr a;
    int ok = 1;

    for (i = 0; i < r; i++)
        k += s[i];
    a = _acb_vec_init(k);
    L = 0;
    for (i = 0; i < r && ok; i++)
    {
        ulong t;
        if (acb_contains_zero(y + i))
            ok = 0;
        for (t = 1; t < s[i]; t++)
            acb_zero(a + L++);
        acb_inv(a + L++, y + i, prec + 32);
    }
    /* convergence at t = 1 (the letter a_1 must not be 1) is checked by
       the evaluation: the reflected letter 1 - a_1 must be nonzero */
    if (ok)
        ok = _acb_mzeta_G_one(res, a, L, prec);
    if (ok && (r & 1))
        acb_neg(res, res);
    if (!ok)
        acb_indeterminate(res);
    _acb_vec_clear(a, k);
}

void
acb_mzeta_polylog(acb_t res, const ulong * s, acb_srcptr x, slong r, slong prec)
{
    slong i, k = 0;
    acb_ptr y;

    if (r == 0)
    {
        acb_one(res);
        return;
    }
    for (i = 0; i < r; i++)
    {
        if (s[i] == 0 || s[i] > (UWORD(1) << 20))
        {
            acb_indeterminate(res);
            return;
        }
        k += s[i];
    }

    y = _acb_vec_init(r);
    acb_set(y, x);
    for (i = 1; i < r; i++)
        acb_mul(y + i, y + i - 1, x + i, prec + 32);
    polylog_y(res, s, y, r, prec);
    _acb_vec_clear(y, r);
    (void) k;
}

void
acb_mzeta_mzv_colored(acb_t res, const ulong * s, const ulong * e, slong r, ulong N, slong prec)
{
    acb_ptr x;
    slong i;

    if (N == 0)
    {
        acb_indeterminate(res);
        return;
    }
    for (i = 0; i < r; i++)
        if (s[i] == 0 || s[i] > (UWORD(1) << 20))
        {
            acb_indeterminate(res);
            return;
        }
    if (r == 0)
    {
        acb_one(res);
        return;
    }
    /* the partial products y_i = exp(2 pi i (e_1 + ... + e_i)/N), exact
       for the exponents 0, N/4, N/2, 3N/4 */
    x = _acb_vec_init(r);
    {
        ulong sum = 0;
        for (i = 0; i < r; i++)
        {
            fmpq_t q;
            fmpq_init(q);
            sum = (sum + e[i] % N) % N;
            fmpq_set_si(q, 2 * sum, N);
            arb_sin_cos_pi_fmpq(acb_imagref(x + i), acb_realref(x + i), q, prec + 32);
            fmpq_clear(q);
        }
    }
    polylog_y(res, s, x, r, prec);
    _acb_vec_clear(x, r);
}
