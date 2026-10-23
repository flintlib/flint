/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "acb.h"
#include "acb_dirichlet.h"
#include "acb_dirichlet/impl.h"
#include "arb_hypgeom.h"
#include "acb_dft.h"
#include "thread_support.h"
#include "dfloat.h"

FLINT_FORCE_INLINE void
_acb_vec_kronecker_mul(acb_ptr z, acb_srcptr x, acb_srcptr y, slong len, slong prec)
{
    slong k;
    for (k = 0; k < len; k++)
        acb_mul(z + k, x + k, y + k, prec);
}

static void
_acb_dot_arb(acb_t res, const acb_t initial, int subtract,
             acb_srcptr x, slong xstep, arb_srcptr y, slong ystep,
             slong len, slong prec)
{
    arb_ptr a;
    arb_srcptr b, c;
    if (sizeof(acb_struct) != 2*sizeof(arb_struct))
    {
        flint_throw(FLINT_ERROR, "expected sizeof(acb_struct)=%zu "
                "to be twice sizeof(arb_struct)=%zu\n",
                sizeof(acb_struct), sizeof(arb_struct));
    }
    if (initial == NULL)
    {
        flint_throw(FLINT_ERROR, "not implemented for NULL initial value\n");
    }

    a = acb_realref(res);
    b = acb_realref(initial);
    c = acb_realref(x);
    arb_dot(a, b, subtract, c, xstep*2, y, ystep, len, prec);

    a = acb_imagref(res);
    b = acb_imagref(initial);
    c = acb_imagref(x);
    arb_dot(a, b, subtract, c, xstep*2, y, ystep, len, prec);
}

static void
_arb_add_d(arb_t z, const arb_t x, double d, slong prec)
{
    arb_t u;
    arb_init(u);
    arb_set_d(u, d);
    arb_add(z, x, u, prec);
    arb_clear(u);
}

static void
_arb_div_si_si(arb_t z, slong a, slong b, slong prec)
{
    arb_set_si(z, a);
    arb_div_si(z, z, b, prec);
}

static void
_arb_inv_si(arb_t z, slong n, slong prec)
{
    arb_set_si(z, n);
    arb_inv(z, z, prec);
}

static void
platt_g_gamma_term(acb_t out, const arb_t t0, const acb_t t, slong prec)
{
    acb_t z;
    acb_init(z);
    acb_add_arb(z, t, t0, prec);
    acb_mul_onei(z, z);
    acb_mul_2exp_si(z, z, 1);
    acb_add_ui(z, z, 1, prec);
    acb_mul_2exp_si(z, z, -2);
    acb_gamma(out, z, prec);
    acb_clear(z);
}

static void
platt_g_exp_term(acb_t out,
        const arb_t t0, const arb_t h, const acb_t t, slong prec)
{
    arb_t pi;
    acb_t z1, z2;
    arb_init(pi);
    acb_init(z1);
    acb_init(z2);
    arb_const_pi(pi, prec);
    acb_add_arb(z1, t, t0, prec);
    acb_mul_arb(z1, z1, pi, prec);
    acb_mul_2exp_si(z1, z1, -2);
    acb_div_arb(z2, t, h, prec);
    acb_sqr(z2, z2, prec);
    acb_mul_2exp_si(z2, z2, -1);
    acb_sub(out, z1, z2, prec);
    acb_exp(out, out, prec);
    arb_clear(pi);
    acb_clear(z1);
    acb_clear(z2);
}

static void
platt_g_base(acb_t out, const acb_t t, slong prec)
{
    arb_t pi;
    arb_init(pi);
    arb_const_pi(pi, prec);
    acb_mul_arb(out, t, pi, prec);
    acb_mul_onei(out, out);
    acb_mul_2exp_si(out, out, 1);
    acb_neg(out, out);
    arb_clear(pi);
}


/* The gamma-exponential factors: the table of the multi-evaluation
   has the entries coeff_i base_i^k (row k, i < N); only coeff and base
   are stored, rounded to the transform precision tprec (the rows are
   generated where they are transformed). */
typedef struct
{
    acb_ptr coeff, base;
    slong A, B, prec, tprec, chunk;
    arb_srcptr t0, h;
}
platt_g_table_arg;

static void
platt_g_table_worker(slong c, void * arg_ptr)
{
    const platt_g_table_arg * a = arg_ptr;
    slong A = a->A, N = a->A * a->B, prec = a->prec;
    slong i, n, i0 = c * a->chunk, i1 = FLINT_MIN(N, (c + 1) * a->chunk);
    acb_t t, gamma_term, exp_term;

    acb_init(t);
    acb_init(gamma_term);
    acb_init(exp_term);

    for (i = i0; i < i1; i++)
    {
        n = i - N/2;

        acb_set_si(t, n);
        acb_div_si(t, t, A, prec);

        platt_g_base(a->base + i, t, a->tprec);

        platt_g_gamma_term(gamma_term, a->t0, t, prec);
        platt_g_exp_term(exp_term, a->t0, a->h, t, prec);
        acb_mul(a->coeff + i, gamma_term, exp_term, a->tprec);
    }

    acb_clear(t);
    acb_clear(gamma_term);
    acb_clear(exp_term);
}

static void
platt_g_table(acb_ptr coeff, acb_ptr base, slong A, slong B,
        const arb_t t0, const arb_t h, slong prec, slong tprec)
{
    platt_g_table_arg a;
    slong N = A*B;
    a.coeff = coeff;
    a.base = base;
    a.A = A;
    a.B = B;
    a.prec = prec;
    a.tprec = tprec;
    a.t0 = t0;
    a.h = h;
    a.chunk = 256;
    flint_parallel_do(platt_g_table_worker, &a, (N + a.chunk - 1) / a.chunk, 0,
        FLINT_PARALLEL_STRIDED);
}

static void
logjsqrtpi(arb_t out, const fmpz_t j, slong prec)
{
    arb_const_sqrt_pi(out, prec);
    arb_mul_fmpz(out, out, j, prec);
    arb_log(out, out, prec);
}

static void
_acb_add_error_arb_mag(acb_t res, const arb_t x)
{
    mag_t err;
    mag_init(err);
    arb_get_mag(err, x);
    acb_add_error_mag(res, err);
    mag_clear(err);
}

static void
_acb_vec_scalar_add_error_mag(acb_ptr res, slong len, const mag_t err)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        acb_add_error_mag(res + i, err);
    }
}

static void
_acb_vec_scalar_add_error_arb_mag(acb_ptr res, slong len, const arb_t x)
{
    mag_t err;
    mag_init(err);
    arb_get_mag(err, x);
    _acb_vec_scalar_add_error_mag(res, len, err);
    mag_clear(err);
}

/*
 * For each integer m in [0, A*B) find the smallest integer j such that
 * log(j*sqrt(pi))/(2*pi) >= m/B - 1/(2*B)
 */
void get_smk_points(fmpz * res, slong A, slong B)
{
    slong m, N, prec;
    arb_t x, u, v;
    fmpz_t z;

    arb_init(x);
    arb_init(u); /* pi / B */
    arb_init(v); /* 1 / sqrt(pi) */
    fmpz_init(z);

    N = A*B;
    prec = 4;
    arb_indeterminate(u);
    arb_indeterminate(v);
    for (m = 0; m < N; m++)
    {
        while (1)
        {
            arb_set_si(x, 2*m - 1);
            arb_mul(x, x, u, prec);
            arb_exp(x, x, prec);
            arb_mul(x, x, v, prec);
            arb_ceil(x, x, prec);
            if (arb_get_unique_fmpz(z, x))
            {
                fmpz_set(res + m, z);
                break;
            }
            else
            {
                prec *= 2;
                arb_const_pi(u, prec);
                arb_div_si(u, u, B, prec);
                arb_const_sqrt_pi(v, prec);
                arb_inv(v, v, prec);
            }
        }
    }
    arb_clear(x);
    arb_clear(u);
    arb_clear(v);
    fmpz_clear(z);
}

slong platt_get_smk_index(slong B, const fmpz_t j, slong prec)
{
    slong m;
    arb_t pi, x;
    fmpz_t z;

    arb_init(pi);
    arb_init(x);
    fmpz_init(z);

    m = -1;

    while (1)
    {
        arb_const_pi(pi, prec);
        logjsqrtpi(x, j, prec);
        arb_div(x, x, pi, prec);
        arb_mul_2exp_si(x, x, -1);
        arb_mul_si(x, x, B, prec);
        _arb_add_d(x, x, 0.5, prec);
        arb_floor(x, x, prec);

        if (arb_get_unique_fmpz(z, x))
        {
            m = fmpz_get_si(z);
            break;
        }
        else
        {
            prec *= 2;
        }
    }

    arb_clear(pi);
    arb_clear(x);
    fmpz_clear(z);

    return m;
}

typedef struct
{
    slong bmax;
    slong b;
    slong K;
    arb_ptr M; /* (b, K) */
    acb_ptr v; /* (b, ) */
}
smk_block_struct;
typedef smk_block_struct smk_block_t[1];

static void
smk_block_init(smk_block_t p, slong K, slong bmax)
{
    p->bmax = bmax;
    p->b = 0;
    p->K = K;
    p->M = _arb_vec_init(K*bmax);
    p->v = _acb_vec_init(bmax);
}

static void
smk_block_clear(smk_block_t p)
{
    _arb_vec_clear(p->M, p->K * p->bmax);
    _acb_vec_clear(p->v, p->bmax);
}

static int
smk_block_is_full(smk_block_t p)
{
    return p->b == p->bmax;
}

static void
smk_block_reset(smk_block_t p)
{
    p->b = 0;
}

static void
smk_block_increment(smk_block_t p, const acb_t z, arb_srcptr v)
{
    if (smk_block_is_full(p))
    {
        flint_throw(FLINT_ERROR, "trying to increment a full block\n");
    }
    acb_set(p->v + p->b, z);
    _arb_vec_set(p->M + p->K * p->b, v, p->K);
    p->b += 1;
}

static void
smk_block_accumulate(smk_block_t p, acb_ptr res, slong prec)
{
    slong i;
    for (i = 0; i < p->K; i++)
        _acb_dot_arb(res + i, res + i, 0, p->v, 1, p->M + i, p->K, p->b, prec);
}

void
_platt_smk(acb_ptr table, acb_ptr startvec, acb_ptr stopvec,
        const fmpz * smk_points, const arb_t t0, slong A, slong B,
        const fmpz_t jstart, const fmpz_t jstop, slong mstart, slong mstop,
        slong K, slong prec)
{
    fmpz_t j, j_plus_one;
    slong k, m;
    slong N = A * B;
    smk_block_t block;
    acb_ptr accum;
    arb_ptr diff_powers;
    arb_t rpi, logsqrtpi, rsqrtj, um, a, base;
    acb_t z;

    fmpz_init(j);
    fmpz_init(j_plus_one);
    arb_init(rpi);
    arb_init(logsqrtpi);
    arb_init(rsqrtj);
    arb_init(um);
    arb_init(a);
    arb_init(base);
    acb_init(z);
    smk_block_init(block, K, 32);
    diff_powers = _arb_vec_init(K);
    accum = _acb_vec_init(K);

    arb_const_pi(rpi, prec);
    arb_inv(rpi, rpi, prec);
    arb_const_sqrt_pi(logsqrtpi, prec);
    arb_log(logsqrtpi, logsqrtpi, prec);

    m = platt_get_smk_index(B, jstart, prec);
    _arb_div_si_si(um, m, B, prec);

    for (fmpz_set(j, jstart); fmpz_cmp(j, jstop) <= 0; fmpz_add_ui(j, j, 1))
    {
        arb_log_fmpz(a, j, prec);
        arb_add(a, a, logsqrtpi, prec);
        arb_mul(a, a, rpi, prec);

        /* todo arb_rsqrt_fmpz */
        arb_sqrt_fmpz(rsqrtj, j, prec);
        arb_inv(rsqrtj, rsqrtj, prec);

        acb_set_arb(z, t0);
        acb_mul_arb(z, z, a, prec);
        acb_neg(z, z);
        acb_exp_pi_i(z, z, prec);
        acb_mul_arb(z, z, rsqrtj, prec);

        while (m < N - 1 && fmpz_cmp(smk_points + m + 1, j) <= 0)
        {
            m += 1;
            _arb_div_si_si(um, m, B, prec);
        }

        if (m < mstart || m > mstop)
        {
            flint_throw(FLINT_ERROR, "out of bounds error: m = %wd not in [%wd, %wd]\n",
                          m, mstart, mstop);
        }

        arb_mul_2exp_si(base, a, -1);
        arb_sub(base, base, um, prec);

        _arb_vec_set_powers(diff_powers, base, K, prec);
        smk_block_increment(block, z, diff_powers);

        {
            int j_stops, m_increases;

            fmpz_add_ui(j_plus_one, j, 1);
            j_stops = fmpz_equal(j, jstop);
            m_increases = (m < N - 1 &&
                           fmpz_cmp(smk_points + m + 1, j_plus_one) <= 0);
            if (j_stops || m_increases || smk_block_is_full(block))
            {
                smk_block_accumulate(block, accum, prec);
                smk_block_reset(block);
            }
            if (j_stops || m_increases)
            {
                if (startvec && m == mstart)
                {
                    _acb_vec_set(startvec, accum, K);
                }
                else if (stopvec && m == mstop)
                {
                    _acb_vec_set(stopvec, accum, K);
                }
                else
                {
                    for (k = 0; k < K; k++)
                        acb_set(table + N*k + m, accum + k);
                }
                _acb_vec_zero(accum, K);
            }
        }
    }

    fmpz_clear(j);
    fmpz_clear(j_plus_one);
    arb_clear(rpi);
    arb_clear(logsqrtpi);
    arb_clear(rsqrtj);
    arb_clear(um);
    arb_clear(a);
    arb_clear(base);
    acb_clear(z);
    smk_block_clear(block);
    _arb_vec_clear(diff_powers, K);
    _acb_vec_clear(accum, K);
}


/* The K convolutions are summed in the frequency domain, so that a
   single inverse transform is needed (2K + 1 transforms of length 2N
   instead of 3K); the products are accumulated over fixed groups of
   PLATT_CONV_GROUP consecutive k in parallel, and the group sums added
   in order, so that the result does not depend on the number of
   threads. */
#define PLATT_CONV_GROUP 4

/* the row k of the sums over j, into row (N entries) at prec */
typedef void (* platt_S_row_func)(acb_ptr row, slong k, const void * ctx, slong prec);

typedef struct
{
    acb_ptr * partial;
    acb_srcptr coeff, base;
    arb_srcptr A7;                      /* the Lemma A7 bounds, k < K */
    platt_S_row_func S_row;
    const void * S_ctx;
    arb_srcptr inv_fac;
    const acb_dft_pre_struct * pre;     /* length 2N */
    const acb_dft_pre_struct * pre_N;   /* length N */
    arb_srcptr t0, h;
    slong A, B, N, K, sigma, prec, g0;
}
platt_conv_arg;

/* The rows k of the group g: the table row coeff_i base_i^k with the
   truncation bound of Lemma A5, its DFT (after swapping the halves)
   divided by A, the bound of Lemma A7 (precomputed for all k) and the
   division by k! (a multiplication by 1/k!); then its convolution with the row k of S, in the frequency
   domain (zero-padded transforms of length 2N), accumulated over the
   group. */
static void
platt_conv_worker(slong gg, void * arg_ptr)
{
    const platt_conv_arg * a = arg_ptr;
    slong N = a->N, prec = a->prec, i, k, g = a->g0 + gg;
    slong k0 = g * PLATT_CONV_GROUP, k1 = FLINT_MIN(a->K, k0 + PLATT_CONV_GROUP);
    acb_ptr row, pw, fp, gp, acc;
    arb_t err;

    arb_init(err);
    row = _acb_vec_init(N);
    pw = _acb_vec_init(N);
    fp = _acb_vec_init(N*2);
    gp = _acb_vec_init(N*2);
    acc = _acb_vec_init(N*2);

    for (i = 0; i < N; i++)
        acb_pow_ui(pw + i, a->base + i, k0, prec);

    for (k = k0; k < k1; k++)
    {
        if (k > k0)
            _acb_vec_kronecker_mul(pw, pw, a->base, N, prec);
        _acb_vec_kronecker_mul(row, a->coeff, pw, N, prec);

        acb_dirichlet_platt_lemma_A5(err, a->B, a->h, k, prec);
        _acb_vec_scalar_add_error_arb_mag(row, N, err);
        for (i = 0; i < N/2; i++)
            acb_swap(row + i, row + i + N/2);
        acb_dft_precomp(row, row, a->pre_N, prec);
        _acb_vec_scalar_div_ui(row, row, N, (ulong) a->A, prec);
        _acb_vec_scalar_add_error_arb_mag(row, N, a->A7 + k);
        if (k >= 2)
            _acb_vec_scalar_mul_arb(row, row, N, a->inv_fac + k, prec);

        /* the S row, zero-padded and reversed (index i -> 2N - i) */
        a->S_row(gp, k, a->S_ctx, prec);
        _acb_vec_zero(fp, N*2);
        acb_set(fp, gp);
        for (i = 1; i < N; i++)
            acb_set(fp + N*2 - i, gp + i);
        /* the table row, zero-padded */
        _acb_vec_set(gp, row, N);
        _acb_vec_zero(gp + N, N);

        acb_dft_precomp(fp, fp, a->pre, prec);
        acb_dft_precomp(gp, gp, a->pre, prec);
        if (k == k0)
            _acb_vec_kronecker_mul(acc, gp, fp, N*2, prec);
        else
        {
            _acb_vec_kronecker_mul(gp, gp, fp, N*2, prec);
            _acb_vec_add(acc, acc, gp, N*2, prec);
        }
    }

    a->partial[g] = acc;

    arb_clear(err);
    _acb_vec_clear(row, N);
    _acb_vec_clear(pw, N);
    _acb_vec_clear(fp, N*2);
    _acb_vec_clear(gp, N*2);
}

/* The groups run in waves of as many groups as threads, whose partial
   sums are added (in the order of the groups) before the next wave,
   so that at most one partial sum per thread is alive. */
static void
do_convolutions(acb_ptr out_table, platt_conv_arg * a_in)
{
    slong i, g, G, W, g0, g1, N = a_in->N, K = a_in->K, prec = a_in->prec;
    acb_ptr total, padded_out_table;
    acb_dft_pre_t pre;
    platt_conv_arg a = *a_in;

    G = (K + PLATT_CONV_GROUP - 1) / PLATT_CONV_GROUP;
    W = FLINT_MAX(1, flint_get_num_threads());
    total = NULL;
    acb_dft_precomp_init(pre, N*2, prec);

    a.partial = flint_calloc(G, sizeof(acb_ptr));
    a.pre = pre;

    for (g0 = 0; g0 < G; g0 = g1)
    {
        g1 = FLINT_MIN(G, g0 + W);
        a.g0 = g0;
        flint_parallel_do(platt_conv_worker, &a, g1 - g0, 0, FLINT_PARALLEL_STRIDED);
        for (g = g0; g < g1; g++)
        {
            if (total == NULL)
                total = a.partial[g];
            else
            {
                _acb_vec_add(total, total, a.partial[g], N*2, prec);
                _acb_vec_clear(a.partial[g], N*2);
            }
        }
    }

    padded_out_table = _acb_vec_init(N*2);
    acb_dft_inverse_precomp(padded_out_table, total, pre, prec);
    _acb_vec_clear(total, N*2);
    flint_free(a.partial);

    for (i = 0; i <= N/2; i++)
    {
        acb_add(out_table + i, out_table + i, padded_out_table + i, prec);
    }

    _acb_vec_clear(padded_out_table, N*2);
    acb_dft_precomp_clear(pre);
}

static void
remove_gaussian_window(arb_ptr out, slong A, slong B, const arb_t h, slong prec)
{
    slong i, n;
    slong N = A*B;
    arb_t t, x;
    arb_init(t);
    arb_init(x);
    for (i = 0; i < N; i++)
    {
        n = i - N/2;
        arb_set_si(t, n);
        arb_div_si(t, t, A, prec);
        arb_div(x, t, h, prec);
        arb_sqr(x, x, prec);
        arb_mul_2exp_si(x, x, -1);
        arb_exp(x, x, prec);
        arb_mul(out + i, out + i, x, prec);
    }
    arb_clear(t);
    arb_clear(x);
}

typedef struct
{
    acb_srcptr S;
    slong N;
}
platt_S_acb_ctx;

typedef struct
{
    const double * S5;
    slong N;
    acb_t c;
}
platt_S_dd_ctx;

static slong
_platt_transform_prec(const arb_t t0, slong prec)
{
    slong dprec;
    dprec = prec - FLINT_MAX(0, arf_abs_bound_lt_2exp_si(arb_midref(t0))) + 16;
    dprec = FLINT_MAX(dprec, 64);
    return FLINT_MIN(dprec, prec);
}

static void
_platt_S_row_acb(acb_ptr row, slong k, const void * ctx, slong prec)
{
    const platt_S_acb_ctx * c = ctx;
    _acb_vec_set(row, c->S + k * c->N, c->N);
}

static void
_platt_S_row_dd(acb_ptr row, slong k, const void * ctx, slong prec)
{
    const platt_S_dd_ctx * c = ctx;
    slong i;
    for (i = 0; i < c->N; i++)
    {
        const double * e = c->S5 + 5 * (k * c->N + i);
        mag_t r;
        arf_t t;
        mag_init(r);
        arf_init(t);
        arf_set_d(arb_midref(acb_realref(row + i)), e[0]);
        arf_set_d(t, e[1]);
        arf_add(arb_midref(acb_realref(row + i)), arb_midref(acb_realref(row + i)), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_d(arb_midref(acb_imagref(row + i)), e[2]);
        arf_set_d(t, e[3]);
        arf_add(arb_midref(acb_imagref(row + i)), arb_midref(acb_imagref(row + i)), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        mag_set_d(r, e[4]);
        mag_set(arb_radref(acb_realref(row + i)), r);
        mag_set(arb_radref(acb_imagref(row + i)), r);
        acb_mul(row + i, row + i, c->c, prec);
        mag_clear(r);
        arf_clear(t);
    }
}

static void
_platt_multieval_rows(arb_ptr out, platt_S_row_func S_row, const void * S_ctx,
        const arb_t t0, slong A, slong B, const arb_t h, const fmpz_t J,
        slong K, slong sigma, slong prec);

void
_acb_dirichlet_platt_multieval(arb_ptr out, acb_srcptr S_table,
        const arb_t t0, slong A, slong B, const arb_t h, const fmpz_t J,
        slong K, slong sigma, slong prec)
{
    platt_S_acb_ctx c;
    c.S = S_table;
    c.N = A * B;
    _platt_multieval_rows(out, _platt_S_row_acb, &c, t0, A, B, h, J, K, sigma, prec);
}

static void
_platt_multieval_rows(arb_ptr out, platt_S_row_func S_row, const void * S_ctx,
        const arb_t t0, slong A, slong B, const arb_t h, const fmpz_t J,
        slong K, slong sigma, slong prec)
{
    /* The transforms and convolutions work on the gamma factors and
       sums already computed at prec bits, of size about 1, whose
       phases (of size t0 log t0) carry about prec - log2(t0) bits: the
       transforms are done at that precision plus 16 bits (at least 64,
       at most prec); with the heuristic parameters, this gives the
       same grid as prec (2 limbs instead of 3 at the zeros from 1e15
       on). */
    slong dprec;
    slong N = A*B;
    slong i, k;
    acb_ptr coeff, base, out_a, out_b;
    arb_t t, x, k_factorial, err, ratio, c, xi;
    acb_t z;
    acb_dft_pre_t pre_N;

    arb_init(t);
    arb_init(x);
    arb_init(k_factorial);
    arb_init(err);
    arb_init(ratio);
    arb_init(c);
    arb_init(xi);
    acb_init(z);
    coeff = _acb_vec_init(N);
    base = _acb_vec_init(N);
    out_a = _acb_vec_init(N);
    out_b = _acb_vec_init(N);
    dprec = _platt_transform_prec(t0, prec);
    acb_dft_precomp_init(pre_N, N, dprec);

    _arb_inv_si(xi, B, prec);
    arb_mul_2exp_si(xi, xi, -1);

    platt_g_table(coeff, base, A, B, t0, h, prec, dprec);

    {
        platt_conv_arg ca;
        arb_ptr A7, inv_fac = _arb_vec_init(K);

        /* 1/k! */
        arb_one(k_factorial);
        for (k = 0; k < K; k++)
        {
            if (k >= 2)
                arb_mul_ui(k_factorial, k_factorial, (ulong) k, prec);
            arb_inv(inv_fac + k, k_factorial, prec);
        }

        /* the Lemma A7 bounds (an upper bound, and tiny: at 64 bits) */
        A7 = _arb_vec_init(K);
        _acb_dirichlet_platt_lemma_A7_vec(A7, sigma, t0, h, K, A, FLINT_MIN(prec, 64));

        ca.coeff = coeff;
        ca.base = base;
        ca.A7 = A7;
        ca.inv_fac = inv_fac;
        ca.pre_N = pre_N;
        ca.t0 = t0;
        ca.h = h;
        ca.A = A;
        ca.B = B;
        ca.N = N;
        ca.K = K;
        ca.sigma = sigma;
        ca.prec = dprec;
        ca.S_row = S_row;
        ca.S_ctx = S_ctx;
        do_convolutions(out_a, &ca);
        _arb_vec_clear(inv_fac, K);
        _arb_vec_clear(A7, K);
    }

    for (i = 0; i < N/2 + 1; i++)
    {
        arb_set_si(x, i);
        arb_div_si(x, x, B, prec);
        acb_dirichlet_platt_lemma_32(err, h, t0, x, prec);
        _acb_add_error_arb_mag(out_a + i, err);
    }

    acb_dirichlet_platt_lemma_B1(err, sigma, t0, h, J, prec);
    _acb_vec_scalar_add_error_arb_mag(out_a, N/2 + 1, err);

    arb_sqrt_fmpz(c, J, prec);
    arb_mul_2exp_si(c, c, 1);
    arb_sub_ui(c, c, 1, prec);
    acb_dirichlet_platt_lemma_B2(err, K, h, xi, prec);
    arb_mul(err, err, c, prec);
    _acb_vec_scalar_add_error_arb_mag(out_a, N/2 + 1, err);

    for (i = 1; i < N/2; i++)
    {
        acb_conj(out_a + N - i, out_a + i);
    }

    acb_dirichlet_platt_lemma_A9(err, sigma, t0, h, A, prec);
    _acb_vec_scalar_add_error_arb_mag(out_a, N, err);

    acb_dft_inverse_precomp(out_b, out_a, pre_N, dprec);
    _acb_vec_scalar_mul_ui(out_b, out_b, N, (ulong) A, prec);
    for (i = 0; i < N/2; i++)
    {
        acb_swap(out_b + i, out_b + i + N/2);
    }

    acb_dirichlet_platt_lemma_A11(err, t0, h, B, prec);
    _acb_vec_scalar_add_error_arb_mag(out_b, N, err);

    for (i = 0; i < N; i++)
    {
        arb_swap(out + i, acb_realref(out_b + i));
    }

    remove_gaussian_window(out, A, B, h, prec);

    arb_clear(t);
    arb_clear(x);
    arb_clear(k_factorial);
    arb_clear(err);
    arb_clear(ratio);
    arb_clear(c);
    arb_clear(xi);
    acb_clear(z);
    _acb_vec_clear(coeff, N);
    _acb_vec_clear(base, N);
    _acb_vec_clear(out_a, N);
    _acb_vec_clear(out_b, N);
    acb_dft_precomp_clear(pre_N);
}

int acb_dirichlet_platt_use_dfloat = 1;

/* The sums over j with dfloat balls (dfloat module): the terms in
   triple-double balls with phases in quad-double balls and the moments
   in triple-double balls, double-double and double arithmetic give a
   grid about as accurate as the arb version at prec = log2(T) + 106
   (about 75 bits after the binary point: the moment k is amplified by
   about 2^(24 + 8.6 k) with the heuristic parameters), for which the
   default precision of the zeta_zeros example asks; the arb version is
   used above that precision, or if the dfloat version fails. */
static int
_platt_multieval_dfloat(arb_ptr out, const fmpz_t T, slong A, slong B,
        const arb_t h, const fmpz_t J, slong K, slong sigma, slong prec)
{
    slong N = A*B;
    platt_S_dd_ctx c;
    double * S5;
    arb_t t0, w;
    fmpz * smk_points;
    int ok;

    if (!acb_dirichlet_platt_use_dfloat || !fmpz_abs_fits_ui(J) ||
        fmpz_sgn(T) <= 0 || prec > (slong) fmpz_bits(T) + 110)
        return 0;

    smk_points = _fmpz_vec_init(N);
    get_smk_points(smk_points, A, B);
    arb_init(t0);
    arb_set_fmpz(t0, T);
    S5 = flint_malloc(sizeof(double) * 5 * K * N);

    /* the sums as compact balls (5 doubles each), without the factor
       c = exp(-i t0 log sqrt(pi)), which is applied as the rows are
       read */
    ok = _dfloat_platt_smk_dd(S5, smk_points, t0, A, B, fmpz_get_ui(J), K, 4, 0);
    if (ok)
    {
        acb_init(c.c);
        arb_init(w);
        arb_const_sqrt_pi(w, prec + 64);
        arb_log(w, w, prec + 64);
        arb_mul(w, w, t0, prec + 64);
        arb_neg(w, w);
        arb_sin_cos(acb_imagref(c.c), acb_realref(c.c), w, prec);
        c.S5 = S5;
        c.N = N;
        _platt_multieval_rows(out, _platt_S_row_dd, &c, t0, A, B, h, J, K, sigma, prec);
        acb_clear(c.c);
        arb_clear(w);
    }

    arb_clear(t0);
    flint_free(S5);
    _fmpz_vec_clear(smk_points, N);
    return ok;
}

void
acb_dirichlet_platt_multieval(arb_ptr out, const fmpz_t T, slong A, slong B,
        const arb_t h, const fmpz_t J, slong K, slong sigma, slong prec)
{
    if (_platt_multieval_dfloat(out, T, A, B, h, J, K, sigma, prec))
        return;

    if (flint_get_num_threads() > 1)
    {
        acb_dirichlet_platt_multieval_threaded(
                out, T, A, B, h, J, K, sigma, prec);
    }
    else
    {
        slong N = A*B;
        acb_ptr S;
        arb_t t0;
        fmpz_t one;
        fmpz * smk_points;

        smk_points = _fmpz_vec_init(N);
        get_smk_points(smk_points, A, B);

        fmpz_init(one);
        fmpz_one(one);

        arb_init(t0);
        S =  _acb_vec_init(K*N);

        arb_set_fmpz(t0, T);

        _platt_smk(S, NULL, NULL, smk_points,
                   t0, A, B, one, J, 0, N-1, K, prec);

        _acb_dirichlet_platt_multieval(out, S, t0, A, B, h, J, K, sigma, prec);

        arb_clear(t0);
        fmpz_clear(one);
        _acb_vec_clear(S, K*N);
        _fmpz_vec_clear(smk_points, N);
    }
}
