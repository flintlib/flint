/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Jet version of acb_dirichlet_l_fmpq_afe: computes the Taylor expansion
    of L(s+x, chi) in x to length len for rational s, using the same
    approximate functional equation

        Gamma((s+kappa)/2) L(s) = S(s) + w (pi/q)^(s-1/2) conj(S(1-s))
                                  [+ pi^(s/2) / (s (s-1)) if q = 1]

        S(s) = sum_n chi(n) n^(-s) Gamma((s+kappa)/2, pi n^2 / q)

    with every quantity replaced by a power series in x. The incomplete
    gamma function jets with respect to the parameter are computed by the
    jet versions of the binary splitting and bit-burst algorithms in
    arb_hypgeom.

    Precision management follows the scalar code; the only changes are
    that magnitude estimates are taken as the maximum over all
    coefficients and that error bounds are valid for all coefficients.
*/

#include <math.h>
#include "fmpq.h"
#include "ulong_extras.h"
#include "arb_poly.h"
#include "arb_hypgeom.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"
#include "thread_pool.h"
#include "thread_support.h"

#define MAG_PREC 30
#define INV_LOG2 1.4426950408889634074
#define PI 3.1415926535897932385

static double
log_gamma_upper_approx(double a, double z)
{
    if (a < z)
        return (a - 1) * log(z) - z;
    else
        return a * (log(a) - 1);
}

/* log2(max_{0 <= j < len} L^j / j!) for L >= 0 */
static double
jet_max_bits(double L, slong len)
{
    double t = 0.0, m = 0.0;
    slong j;

    L = fabs(L);

    for (j = 1; j < len; j++)
    {
        t += log(L / j) * INV_LOG2;
        m = FLINT_MAX(m, t);
        if (j > L)
            break;
    }

    return m;
}

/* log2(sum_{0 <= j < len} L^j / j!) for L >= 0 */
static double
jet_sum_bits(double L, slong len)
{
    double t = 1.0, s = 1.0;
    slong j;

    L = fabs(L);

    if (len > L + 50)
        return L * INV_LOG2;

    for (j = 1; j < len; j++)
    {
        t = t * L / j;
        s += t;
    }

    return log(s) * INV_LOG2;
}

static void
acb_dirichlet_root_number2(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    acb_dirichlet_root_number(res, G, chi, prec);

    if (dirichlet_char_is_real(G, chi))
        arb_zero(acb_imagref(res));
}

static void
_arf_trunc(arf_t x)
{
    if (arf_sgn(x) < 0)
        arf_ceil(x, x);
    else
        arf_floor(x, x);
}

static void
arb_extract_bits(arb_t t, const arb_t z, slong b)
{
    arb_mul_2exp_si(t, z, b);
    _arf_trunc(arb_midref(t));
    mag_zero(arb_radref(t));
    arb_mul_2exp_si(t, t, -b);
}

/*
    Bound, valid for every coefficient of the power series in y, for

        sum_{n >= N} |chi(n) n^(-sigma-y) Gamma(a + y/2, pi n^2/q)|,
        a = (sigma + parity)/2 = sd2.

    By Cauchy's estimate on |y| = 1, using t >= z >= 1, the coefficients of
    each term are bounded by n^(1-sigma) Gamma(b, z) with b = a + 1/2. If
    z >= max(1, 2(b-1)), then Gamma(b, z) <= K z^(b-1) e^(-z) with K = 2
    (K = 1 if b <= 1), which gives the bound

        K (pi/q)^(b-1) n^parity e^(-pi n^2/q).

    Summing over n >= N with (1+k/N)^parity <= e^(k parity/N) gives

        K (pi/q)^(b-1) N^parity e^(-pi N^2/q) / (1 - e^(-delta)),
        delta = 2 pi N / q - parity / N > 0.
*/
static void
acb_dirichlet_afe_tail_bound_jet(mag_t res, const fmpq_t sd2, slong N, ulong q, int parity)
{
    arb_t z, t, u;
    fmpq_t e;
    mag_t zl, d, m;

    arb_init(z);
    arb_init(t);
    arb_init(u);
    fmpq_init(e);
    mag_init(zl);
    mag_init(d);
    mag_init(m);

    /* z = pi N^2 / q */
    arb_const_pi(z, MAG_PREC);
    arb_mul_ui(z, z, N, MAG_PREC);
    arb_mul_ui(z, z, N, MAG_PREC);
    arb_div_ui(z, z, q, MAG_PREC);
    arb_get_mag_lower(zl, z);

    /* e = b - 1 = sd2 - 1/2 */
    fmpq_set_si(e, 1, 2);
    fmpq_sub(e, sd2, e);

    /* t = 2(b-1) */
    arb_set_fmpq(t, e, MAG_PREC);
    arb_mul_2exp_si(t, t, 1);

    /* delta */
    arb_const_pi(u, MAG_PREC);
    arb_mul_ui(u, u, 2 * N, MAG_PREC);
    arb_div_ui(u, u, q, MAG_PREC);
    if (parity)
    {
        arb_t v;
        arb_init(v);
        arb_set_ui(v, N);
        arb_inv(v, v, MAG_PREC);
        arb_sub(u, u, v, MAG_PREC);
        arb_clear(v);
    }

    if (mag_cmp_2exp_si(zl, 0) < 0 || !arb_is_positive(u) ||
        (fmpq_sgn(e) > 0 && !arb_lt(t, z)))
    {
        mag_inf(res);
    }
    else
    {
        /* 1 + 1/delta */
        arb_get_mag_lower(d, u);
        mag_inv(d, d);
        mag_add_ui(d, d, 1);

        /* (pi/q)^(b-1) */
        arb_const_pi(t, MAG_PREC);
        arb_div_ui(t, t, q, MAG_PREC);
        arb_pow_fmpq(t, t, e, MAG_PREC);
        arb_get_mag(m, t);

        mag_expinv(res, zl);
        mag_mul(res, res, m);
        mag_mul(res, res, d);

        if (parity)
            mag_mul_ui(res, res, N);

        if (fmpq_sgn(e) > 0)
            mag_mul_2exp_si(res, res, 1);
    }

    arb_clear(z);
    arb_clear(t);
    arb_clear(u);
    fmpq_clear(e);
    mag_clear(zl);
    mag_clear(d);
    mag_clear(m);
}

/* The bound for the scalar case (len = 1), as in l_fmpq_afe.c */
static void
acb_dirichlet_afe_tail_bound(mag_t res, const fmpq_t sd2, slong N, ulong q, int parity)
{
    mag_t pi_n2_q, t, u;
    fmpz_t sprime;

    mag_init(pi_n2_q);
    mag_init(t);
    mag_init(u);
    fmpz_init(sprime);

    /* pi_n2_q = pi * N^2 / q (lower bound) */
    mag_const_pi_lower(pi_n2_q);
    mag_mul_ui_lower(pi_n2_q, pi_n2_q, N);
    mag_mul_ui_lower(pi_n2_q, pi_n2_q, N);
    mag_set_ui(t, q);
    mag_div_lower(pi_n2_q, pi_n2_q, t);

    /* upper bound for sd2 */
    fmpz_cdiv_q(sprime, fmpq_numref(sd2), fmpq_denref(sd2));

    /* require pi_n2_q > s' */
    mag_set_fmpz(t, sprime);
    if (fmpz_sgn(sprime) > 0 && mag_cmp(pi_n2_q, t) <= 0)
    {
        mag_inf(res);
    }
    else
    {
        mag_expinv(res, pi_n2_q);

        mag_div_ui(res, res, N);
        if (!parity)
            mag_div_ui(res, res, N);

        /* (1 + q/pi) */
        mag_set_ui(t, q);
        mag_const_pi_lower(u);
        mag_div(t, t, u);
        mag_add_ui(t, t, 1);
        mag_mul(res, res, t);

        /* max(1, 2^s') */
        if (fmpz_sgn(sprime) > 0)
            mag_mul_2exp_fmpz(res, res, sprime);

        /* (pi/q)^(s'-1) */
        fmpz_sub_ui(sprime, sprime, 1);
        if (fmpz_sgn(sprime) >= 0)
        {
            mag_const_pi(t);
            mag_set_ui_lower(u, q);
            mag_div(t, t, u);
            mag_pow_fmpz(t, t, sprime);
        }
        else
        {
            mag_const_pi_lower(t);
            mag_set_ui(u, q);
            mag_div_lower(t, t, u);

            fmpz_neg(sprime, sprime);
            mag_pow_fmpz_lower(t, t, sprime);
            mag_inv(t, t);
        }

        mag_mul(res, res, t);
    }

    mag_clear(pi_n2_q);
    mag_clear(t);
    mag_clear(u);
    fmpz_clear(sprime);
}

/* res = 1/Gamma(a + x/2), a rational */
static void
_arb_rgamma_fmpq_jet_half(arb_ptr res, const fmpq_t a, slong len, slong prec)
{
    arb_ptr t;
    slong j;

    t = _arb_vec_init(len);

    if (len == 1)
    {
        /* the value: no cancellation problem */
        if (fmpz_is_one(fmpq_denref(a)) && fmpz_sgn(fmpq_numref(a)) <= 0)
        {
            arb_zero(res);
        }
        else
        {
            arb_gamma_fmpq(t, a, prec);
            arb_inv(res, t, prec);
        }
    }
    else if (fmpq_cmp_ui(a, 1) < 0)
    {
        /* 1/Gamma(a+e) = (a+e) (a+1+e) ... (a+k-1+e) / Gamma(a+k+e) with
           1 <= a+k < 2. Inverting the series of Gamma(a+e) directly would
           lose about len log2(1/dist(a, poles)) bits to cancellation. */
        slong k, i;
        fmpq_t f;
        fmpz_t kk;
        arb_t c;

        fmpq_init(f);
        fmpz_init(kk);
        arb_init(c);

        /* k = ceil(1 - a) */
        fmpq_sub_ui(f, a, 1);
        fmpq_neg(f, f);
        fmpz_cdiv_q(kk, fmpq_numref(f), fmpq_denref(f));
        k = fmpz_get_si(kk);
        fmpq_add_fmpz(f, a, kk);

        _arb_hypgeom_gamma_fmpq_jet(t, f, len, prec);
        _arb_poly_inv_series(res, t, len, len, prec);

        for (i = 0; i < k; i++)
        {
            /* res <- res * (a + i + e) */
            fmpq_add_si(f, a, i);
            arb_set_fmpq(c, f, prec);

            for (j = len - 1; j >= 1; j--)
            {
                arb_mul(res + j, res + j, c, prec);
                arb_add(res + j, res + j, res + j - 1, prec);
            }
            arb_mul(res, res, c, prec);
        }

        fmpq_clear(f);
        fmpz_clear(kk);
        arb_clear(c);
    }
    else
    {
        _arb_hypgeom_gamma_fmpq_jet(t, a, len, prec);
        _arb_poly_inv_series(res, t, len, len, prec);
    }

    for (j = 1; j < len; j++)
        arb_mul_2exp_si(res + j, res + j, -j);

    _arb_vec_clear(t, len);
}

/* res = c^x as power series, given log(c) */
static void
_arb_exp_jet(arb_ptr res, const arb_t logc, slong len, slong prec)
{
    slong j;

    arb_one(res);
    for (j = 1; j < len; j++)
    {
        arb_mul(res + j, res + j - 1, logc, prec);
        arb_div_ui(res + j, res + j, j, prec);
    }
}

static void
_acb_vec_mullow_arb_vec(acb_ptr res, acb_srcptr x, arb_srcptr y, slong len, slong prec)
{
    arb_ptr a, b;
    slong j;

    a = _arb_vec_init(len);
    b = _arb_vec_init(len);

    for (j = 0; j < len; j++)
        arb_set(a + j, acb_realref(x + j));
    _arb_poly_mullow(b, a, len, y, len, len, prec);
    for (j = 0; j < len; j++)
        arb_set(a + j, acb_imagref(x + j));
    for (j = 0; j < len; j++)
        arb_swap(acb_realref(res + j), b + j);
    _arb_poly_mullow(b, a, len, y, len, len, prec);
    for (j = 0; j < len; j++)
        arb_swap(acb_imagref(res + j), b + j);

    _arb_vec_clear(a, len);
    _arb_vec_clear(b, len);
}

/* Given exact dyadic numbers 0 < z0, z1 with z1 = z0 (1 + small) and an
   enclosure of log(z0), sets res to an enclosure of log(z1), using
   log(z1) = log(z0) + 2 atanh((z1 - z0) / (z1 + z0)). This is much cheaper
   than a full logarithm since the atanh argument is a small rational
   number of low height. */
static void
_arb_log_update(arb_t res, const arb_t logz0, const arb_t z0, const arb_t z1, slong prec)
{
    fmpz_t A, B;
    slong e;

    if (!arb_is_exact(z0) || !arb_is_exact(z1) || arf_is_zero(arb_midref(z0)))
    {
        arb_log(res, z1, prec);
        return;
    }

    fmpz_init(A);
    fmpz_init(B);

    e = FLINT_MIN(arf_abs_bound_lt_2exp_si(arb_midref(z0)) - arf_bits(arb_midref(z0)),
                  arf_abs_bound_lt_2exp_si(arb_midref(z1)) - arf_bits(arb_midref(z1)));

    if (!arf_get_fmpz_fixed_si(A, arb_midref(z1), e) && !arf_get_fmpz_fixed_si(B, arb_midref(z0), e))
    {
        arb_t t;
        arb_init(t);
        fmpz_sub(A, A, B);
        fmpz_mul_2exp(B, B, 1);
        fmpz_add(B, B, A);
        /* A = z1 - z0, B = z1 + z0 */
        arb_atan_frac_bsplit(t, A, B, 1, prec);
        arb_mul_2exp_si(t, t, 1);
        arb_add(res, logz0, t, prec);
        arb_clear(t);
    }
    else
    {
        arb_log(res, z1, prec);
    }

    fmpz_clear(A);
    fmpz_clear(B);
}

/*
    The sums S(s + y) = sum_n chi(n) n^(-s-y) Gamma((s+parity)/2 + y/2,
    pi n^2/q). The terms are independent. A sum is evaluated in phases:

    1. (serial) the number of terms N, from the tail bound;
    2. (parallel) planning of each term: the character value, the tolerance
       for the incomplete gamma function, the starting point z0 of the
       bit-burst evaluation and whether the asymptotic series can be used
       there (if not, the term needs Gamma((s+parity)/2 + eps));
    3. (serial) Gamma((s+parity)/2 + eps), shared by the terms that need
       it, to the smallest tolerance of those terms;
    4. (parallel) evaluation of the terms, which are added to the sum in
       order of increasing n as soon as all earlier terms are done, so that
       the result does not depend on the number of threads;
    5. (serial) adding the tail bound.

    The work in phases 2 and 4 is handed out dynamically to the threads,
    in order of increasing n, which roughly means decreasing cost. Terms
    that finish before an earlier term wait in a small ring buffer; a
    thread waits before starting a new term when the buffer is full, so
    the memory is O(threads * len * prec) rather than O(N * len * prec).
*/

#define AFE_TERM_ZERO 0         /* chi(n) = 0 */
#define AFE_TERM_ASYMP 1        /* asymptotic series at z0 */
#define AFE_TERM_LOWER 2        /* Gamma(a) - lower series at z0 */
#define AFE_TERM_SINGULAR 3     /* series at z0 for a nonpositive integer */

typedef struct
{
    fmpq_t s;
    fmpq_t s2;                   /* (s + parity) / 2 */
    const dirichlet_group_struct * G;
    const dirichlet_char_struct * chi;
    mag_t abs_tol;
    slong len;
    slong prec;
    int parity;
    int gamma_singular;
    int is_real;
    ulong q;
    slong N;                     /* terms 1 <= n < N */
    mag_t tail;                  /* bound for the tail n >= N */
    double gamma_mag;            /* log2 of max |Gamma((s+parity)/2+eps)| coefficient */
    double abs_tol_mag;          /* log2(abs_tol) */
    double s_d, s2_d;            /* s, s2 as doubles */

    /* per term, indexed by n - 1 */
    int * method;
    ulong * expo;                /* chi(n) = exp(2 pi i expo / G->expo) */
    arb_ptr z0;
    slong * NN;                  /* number of terms of the series at z0 */
    slong * wp;                  /* working precision sans cancellation */
    slong * wp2;                 /* working precision with cancellation */
    mag_ptr AE;                  /* its truncation error */
    mag_ptr tol;                 /* tolerance for the incomplete gamma function */

    /* Gamma(s2 + eps), or Gamma(1 + eps) of length len + 1 if gamma_singular
       (only if needed by some term) */
    arb_ptr Ga;
    slong Ga_len;

    acb_ptr sum;                 /* sum of the terms added so far (len) */
}
afe_sum_struct;

static void
afe_sum_init(afe_sum_struct * S, const fmpq_t s, const dirichlet_group_t G,
    const dirichlet_char_t chi, const mag_t abs_tol, slong len, slong prec)
{
    slong n;

    fmpq_init(S->s);
    fmpq_init(S->s2);
    mag_init(S->abs_tol);
    mag_init(S->tail);

    fmpq_set(S->s, s);
    mag_set(S->abs_tol, abs_tol);
    S->G = G;
    S->chi = chi;
    S->len = len;
    S->prec = prec;

    if (G == NULL)
    {
        S->parity = 0;
        S->q = 1;
        S->is_real = 1;
    }
    else
    {
        S->parity = dirichlet_parity_char(G, chi);
        S->q = G->q;
        S->is_real = dirichlet_char_is_real(G, chi);
    }

    /* s2 = (s+parity)/2 */
    fmpq_add_ui(S->s2, s, S->parity);
    fmpq_div_2exp(S->s2, S->s2, 1);

    S->gamma_singular = (fmpz_is_one(fmpq_denref(S->s2)) && fmpz_sgn(fmpq_numref(S->s2)) <= 0);

    S->abs_tol_mag = mag_get_d_log2_approx(abs_tol);
    S->s_d = fmpq_get_d(s);
    S->s2_d = fmpq_get_d(S->s2);

    /* number of terms */
    for (n = 1; ; n++)
    {
        if (len == 1)
            acb_dirichlet_afe_tail_bound(S->tail, S->s2, n, S->q, S->parity);
        else
            acb_dirichlet_afe_tail_bound_jet(S->tail, S->s2, n, S->q, S->parity);

        if (mag_cmp(S->tail, abs_tol) < 0)
            break;
    }

    S->N = n;

    /* estimate of the magnitude of Gamma((s+parity)/2+eps), all
       coefficients (only the exponents are used) */
    S->gamma_mag = 0.0;
    if (!S->gamma_singular && S->N > 1 && len == 1)
    {
        /* log2 |Gamma(a)| + 1 (as ARF_EXP of the value) */
        double lg = lgamma(fmpq_get_d(S->s2)) * INV_LOG2;
        S->gamma_mag = (lg > -1e15 && lg < 1e15) ? floor(lg) + 1 : 0.0;
    }
    else if (!S->gamma_singular && S->N > 1)
    {
        arb_ptr t;
        slong j;

        t = _arb_vec_init(len);
        _arb_hypgeom_gamma_fmpq_jet(t, S->s2, len, 64);

        S->gamma_mag = -WORD_MAX;
        for (j = 0; j < len; j++)
            if (!arb_is_zero(t + j))
                S->gamma_mag = FLINT_MAX(S->gamma_mag, ARF_EXP(arb_midref(t + j)));

        _arb_vec_clear(t, len);
    }

    n = FLINT_MAX(S->N - 1, 1);
    S->method = flint_malloc(sizeof(int) * n);
    S->expo = flint_malloc(sizeof(ulong) * n);
    S->z0 = _arb_vec_init(n);
    S->NN = flint_malloc(sizeof(slong) * 3 * n);
    S->wp = S->NN + n;
    S->wp2 = S->NN + 2 * n;
    S->AE = _mag_vec_init(n);
    S->tol = _mag_vec_init(n);
    S->sum = _acb_vec_init(len);
    S->Ga = NULL;
    S->Ga_len = 0;
}

static void
afe_sum_clear(afe_sum_struct * S)
{
    slong n = FLINT_MAX(S->N - 1, 1);

    fmpq_clear(S->s);
    fmpq_clear(S->s2);
    mag_clear(S->abs_tol);
    mag_clear(S->tail);
    flint_free(S->method);
    flint_free(S->expo);
    _arb_vec_clear(S->z0, n);
    flint_free(S->NN);
    _mag_vec_clear(S->AE, n);
    _mag_vec_clear(S->tol, n);
    _acb_vec_clear(S->sum, S->len);
    if (S->Ga != NULL)
        _arb_vec_clear(S->Ga, S->Ga_len);
}

/* log2 of the tolerance for the incomplete gamma function in term n:
   we want Gamma(a,z) n^-s with abs_tol, i.e. Gamma(a,z) with
   abs_tol * n^s, divided by the growth of the coefficients of n^-y */
static double
afe_term_tol_mag(const afe_sum_struct * S, slong n)
{
    return S->abs_tol_mag + S->s_d * log(n) * INV_LOG2
        - ((S->len > 1) ? jet_sum_bits(log(n), S->len) : 0.0);
}

/* working precisions for term n: wp sans cancellation, wp2 with
   cancellation (series at 0) */
static void
afe_term_precs(slong * wp_out, slong * wp2_out, const afe_sum_struct * S, slong n)
{
    slong wp, wp2, len = S->len;
    double gammainc_mag, ns_mag, n_extra, z_extra, aa, zz, abs_tol_mag;

    abs_tol_mag = S->abs_tol_mag;
    aa = S->s2_d;
    zz = PI * n * n / S->q;

    /* growth of the coefficients of z^eps and n^(-y) */
    z_extra = (len > 1) ? jet_max_bits(log(zz), len) : 0.0;
    n_extra = (len > 1) ? jet_sum_bits(log(n), len) : 0.0;

    /* Gamma((s+parity)/2, z)  (want lower bound, to estimate cancellation) */
    gammainc_mag = log_gamma_upper_approx(aa, zz) * INV_LOG2 + z_extra;
    /* n^-s */
    ns_mag = -S->s_d * log(n) * INV_LOG2;

    /* wp = Precision needed sans cancellation. */
    wp = gammainc_mag + ns_mag + n_extra - abs_tol_mag + 5;
    wp = FLINT_MAX(wp, 30);

    /* wp2 = Precision needed with cancellation. */
    if (S->gamma_singular)
    {
        /* Max term is roughly n^(-s) * z^((s+parity)/2) * exp(z) */
        wp2 = -abs_tol_mag + ns_mag + n_extra + aa * log(zz) * INV_LOG2 + zz * INV_LOG2 + z_extra + 5;
        wp2 = FLINT_MAX(wp2, 30);
    }
    else
    {
        wp2 = FLINT_MAX(S->gamma_mag, gammainc_mag) + ns_mag + n_extra - abs_tol_mag + 5;
        wp2 = FLINT_MAX(wp2, 30);
    }

    *wp_out = wp;
    *wp2_out = wp2;
}

/* phase 2: plan term n */
static void
afe_sum_plan(afe_sum_struct * S, slong n)
{
    slong i = n - 1, wp, wp2, zprec;
    arb_t z;

    if (S->G == NULL)
        S->expo[i] = 0;
    else
        S->expo[i] = dirichlet_chi(S->G, S->chi, n);

    if (S->expo[i] == DIRICHLET_CHI_NULL)
    {
        S->method[i] = AFE_TERM_ZERO;
        return;
    }

    mag_set_ui_2exp_si(S->tol + i, 1, afe_term_tol_mag(S, n));

    /* starting point: pi n^2 / q truncated to 32 fractional bits, from an
       approximation with at most 128 bits (when the working precision is
       lower, z0 has fewer bits, which makes the series at z0 cheaper) */
    afe_term_precs(&wp, &wp2, S, n);
    S->wp[i] = wp;
    S->wp2[i] = wp2;
    zprec = FLINT_MIN(wp2, 128);
    arb_init(z);
    arb_const_pi(z, zprec);
    arb_mul_ui(z, z, n, zprec);
    arb_mul_ui(z, z, n, zprec);
    arb_div_ui(z, z, S->q, zprec);
    arb_extract_bits(S->z0 + i, z, 32);
    arb_clear(z);

    /* Can we use the asymptotic series? */
    S->NN[i] = _arb_hypgeom_gamma_upper_fmpq_inf_choose_N_jet(S->AE + i, S->s2, S->z0 + i, S->tol + i, S->len);

    if (S->NN[i] != -1)
        S->method[i] = AFE_TERM_ASYMP;
    else if (S->gamma_singular)
        S->method[i] = AFE_TERM_SINGULAR;
    else
        S->method[i] = AFE_TERM_LOWER;
}

/* phase 3: the shared gamma function jet */
static void
afe_sum_gamma(afe_sum_struct * S)
{
    slong i, j, gprec, len = S->len;
    int need = 0;
    mag_t tol, t;

    mag_init(tol);
    mag_init(t);
    mag_inf(tol);

    for (i = 0; i < S->N - 1; i++)
    {
        /* the scalar singular series does not use the gamma function */
        if (S->method[i] == AFE_TERM_LOWER ||
            (S->method[i] == AFE_TERM_SINGULAR && len > 1))
        {
            need = 1;
            mag_min(tol, tol, S->tol + i);
        }
    }

    if (need)
    {
        fmpq_t one;
        fmpq_init(one);
        fmpq_one(one);

        S->Ga_len = S->gamma_singular ? len + 1 : len;
        S->Ga = _arb_vec_init(S->Ga_len);

        gprec = S->prec * 1.05 + 30;

        while (1)
        {
            _arb_hypgeom_gamma_fmpq_jet(S->Ga, S->gamma_singular ? (const fmpq *) one : S->s2,
                S->Ga_len, gprec);

            mag_zero(t);
            for (j = 0; j < S->Ga_len; j++)
                mag_max(t, t, arb_radref(S->Ga + j));

            /* the cap is a safeguard; the terms then lose accuracy but
               remain correct */
            if (mag_cmp(t, tol) <= 0 || gprec > 64 * S->prec + 1000)
                break;

            gprec *= 2;
        }

        fmpq_clear(one);
    }

    mag_clear(tol);
    mag_clear(t);
}

/* phase 4: compute term n into res (length len); uses only local storage
   apart from reading the shared data, so terms can be computed in parallel */
static void
afe_sum_term(acb_ptr res, afe_sum_struct * S, slong n)
{
    slong i = n - 1, j, start_bits, bits, wp, wp2;
    mag_t AE;
    arb_t ns, t, z, z0, z1, expmz0, logn, logz0, logz1;
    arb_ptr Gz0, Gz1, P, T;
    acb_t c;
    slong len = S->len, prec = S->prec;
    ulong q = S->q;
    const fmpq * s2 = S->s2;
    const mag_struct * abs_tol_gamma = S->tol + i;

    if (S->method[i] == AFE_TERM_ZERO)
    {
        _acb_vec_zero(res, len);
        return;
    }

    acb_init(c);
    mag_init(AE);
    arb_init(ns);
    arb_init(t);
    arb_init(z);
    arb_init(z0);
    arb_init(z1);
    arb_init(expmz0);
    arb_init(logn);
    arb_init(logz0);
    arb_init(logz1);
    Gz0 = _arb_vec_init(4 * len);
    Gz1 = Gz0 + len;
    P = Gz0 + 2 * len;
    T = Gz0 + 3 * len;

    /* local precision (from the plan) */
    wp = S->wp[i];
    wp2 = S->wp2[i];

    /* chi(n) */
    if (S->G == NULL)
    {
        acb_one(c);
    }
    else
    {
        fmpq_t e;
        fmpq_init(e);
        fmpq_set_si(e, 2 * S->expo[i], S->G->expo);
        arb_sin_cos_pi_fmpq(acb_imagref(c), acb_realref(c), e, wp);
        fmpq_clear(e);
    }

    arb_const_pi(z, wp2);
    arb_mul_ui(z, z, n, wp2);
    arb_mul_ui(z, z, n, wp2);
    arb_div_ui(z, z, q, wp2);

    start_bits = 32;
    arb_set(z0, S->z0 + i);

    /* log(z0) is needed for the jet prefactors z0^eps */
    if (len > 1)
        arb_log(logz0, z0, wp2);

    if (S->method[i] == AFE_TERM_ASYMP)
    {
        _arb_hypgeom_gamma_upper_fmpq_inf_bsplit_jet(Gz0, s2, z0, (len > 1) ? logz0 : NULL, S->NN[i], len, wp);
        for (j = 0; j < len; j++)
            arb_add_error_mag(Gz0 + j, S->AE + i);
    }
    else
    {
        arb_ptr Ga = S->Ga;
        slong NN;

        if (S->method[i] == AFE_TERM_SINGULAR)
        {
            slong nn = *fmpq_numref(s2);

            if (COEFF_IS_MPZ(nn))
            {
                _arb_vec_indeterminate(Gz0, len);
            }
            else
            {
                nn = -nn;
                NN = _arb_hypgeom_gamma_upper_singular_si_choose_N_jet(AE, nn, z0, abs_tol_gamma, len);
                _arb_hypgeom_gamma_upper_singular_si_bsplit_jet(Gz0, nn, z0, (len > 1) ? logz0 : NULL, NN, (len > 1) ? Ga : NULL, len, wp2);
                for (j = 0; j < len; j++)
                    arb_add_error_mag(Gz0 + j, AE);
            }
        }
        else
        {
            NN = _arb_hypgeom_gamma_lower_fmpq_0_choose_N_jet(AE, s2, z0, abs_tol_gamma, len);
            _arb_hypgeom_gamma_lower_fmpq_0_bsplit_jet(Gz0, s2, z0, (len > 1) ? logz0 : NULL, NN, len, wp2);
            for (j = 0; j < len; j++)
                arb_add_error_mag(Gz0 + j, AE);

            _arb_vec_sub(Gz0, Ga, Gz0, len, wp);
        }
    }

    arb_neg(expmz0, z0);
    arb_exp(expmz0, expmz0, wp);

    /* Bit-burst steps */
    for (bits = start_bits * 2; bits < wp / 8; bits *= 2)
    {
        arb_extract_bits(z1, z, bits);
        _arb_gamma_upper_fmpq_step_bsplit_jet(Gz1, s2, z0, (len > 1) ? logz0 : NULL, z1, Gz0, expmz0, abs_tol_gamma, len, wp);
        arb_sub(t, z0, z1, wp);
        arb_exp(t, t, wp);
        arb_mul(expmz0, expmz0, t, wp);
        _arb_vec_swap(Gz0, Gz1, len);
        if (len > 1)
        {
            _arb_log_update(logz1, logz0, z0, z1, wp);
            arb_swap(logz0, logz1);
        }
        arb_set(z0, z1);
    }

    /* Final step, including error bound */
    _arb_gamma_upper_fmpq_step_bsplit_jet(Gz1, s2, z0, (len > 1) ? logz0 : NULL, z, Gz0, expmz0, abs_tol_gamma, len, wp);
    _arb_vec_swap(Gz0, Gz1, len);

    /* Convert from eps = y/2 to y */
    for (j = 1; j < len; j++)
        arb_mul_2exp_si(Gz0 + j, Gz0 + j, -j);

    /* Multiply by prefactor n^(-s-y) */
    arb_set_ui(ns, n);
    arb_pow_fmpq(ns, ns, S->s, wp);
    if (len > 1)
    {
        arb_inv(ns, ns, wp);
        arb_log_ui(logn, n, wp);
        arb_neg(logn, logn);
        _arb_exp_jet(P, logn, len, wp);
        _arb_vec_scalar_mul(P, P, len, ns, wp);
        _arb_poly_mullow(T, Gz0, len, P, len, len, wp);
    }
    else
    {
        arb_div(T, Gz0, ns, wp);
    }

    for (j = 0; j < len; j++)
        acb_mul_arb(res + j, c, T + j, prec);

    mag_clear(AE);
    arb_clear(ns);
    arb_clear(t);
    arb_clear(z);
    arb_clear(z0);
    arb_clear(z1);
    arb_clear(expmz0);
    arb_clear(logn);
    arb_clear(logz0);
    arb_clear(logz1);
    _arb_vec_clear(Gz0, 4 * len);
    acb_clear(c);
}

/* phase 5: res = sum of the terms, plus the tail bound */
static void
afe_sum_finish(acb_ptr res, const afe_sum_struct * S)
{
    slong j, len = S->len;

    _acb_vec_set(res, S->sum, len);

    for (j = 0; j < len; j++)
    {
        if (S->is_real)
            arb_add_error_mag(acb_realref(res + j), S->tail);
        else
            acb_add_error_mag(res + j, S->tail);
    }
}

/*
    Parallel work queue over the terms of one or two sums. Work item i is
    term n = i / nsums + 1 of sum i % nsums. In the evaluation phase, item i
    is computed into slot i % nslots of a ring buffer, and the items are
    added to their sums in the order of i (hence each sum in the order of
    n) as soon as all earlier items are done. Items are handed out in
    increasing order, so the items in progress or waiting to be added are
    lo <= i < next, and a thread waits before taking item next when
    next - lo = nslots.
*/
typedef struct
{
    afe_sum_struct * S[2];
    slong nsums;
    slong num;                   /* number of work items */
    slong next;                  /* next item to hand out */
    int phase;                   /* 0: plan, 1: evaluate */
    int threaded;

    /* evaluation phase */
    slong len;
    slong nslots;
    acb_ptr slots;               /* nslots * len */
    char * done;                 /* nslots */
    slong lo;                    /* first item not yet added */
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
    pthread_cond_t cond;
#endif
}
afe_work_struct;

/* add the finished items lo, lo + 1, ... to the sums (with the lock held) */
static void
afe_work_flush(afe_work_struct * W)
{
    while (W->lo < W->next && W->done[W->lo % W->nslots])
    {
        slong i = W->lo;
        afe_sum_struct * S = W->S[i % W->nsums];
        slong n = i / W->nsums + 1;

        if (n < S->N)
            _acb_vec_add(S->sum, S->sum, W->slots + (i % W->nslots) * W->len, W->len, S->prec);

        W->done[i % W->nslots] = 0;
        W->lo++;
    }
}

static void
afe_worker(void * arg)
{
    afe_work_struct * W = (afe_work_struct *) arg;
    afe_sum_struct * S;
    slong i, n;

    while (1)
    {
#if FLINT_USES_PTHREAD
        if (W->threaded)
        {
            pthread_mutex_lock(&W->mutex);
            if (W->phase == 1)
                while (W->next < W->num && W->next - W->lo >= W->nslots)
                    pthread_cond_wait(&W->cond, &W->mutex);
            i = W->next++;
            pthread_mutex_unlock(&W->mutex);
        }
        else
#endif
        {
            i = W->next++;
        }
        if (i >= W->num)
            break;

        n = i / W->nsums + 1;
        S = W->S[i % W->nsums];

        if (W->phase == 0)
        {
            if (n < S->N)
                afe_sum_plan(S, n);
        }
        else
        {
            if (n < S->N)
                afe_sum_term(W->slots + (i % W->nslots) * W->len, S, n);

#if FLINT_USES_PTHREAD
            if (W->threaded)
            {
                pthread_mutex_lock(&W->mutex);
                W->done[i % W->nslots] = 1;
                afe_work_flush(W);
                pthread_cond_broadcast(&W->cond);
                pthread_mutex_unlock(&W->mutex);
            }
            else
#endif
            {
                W->done[i % W->nslots] = 1;
                afe_work_flush(W);
            }
        }
    }
}

/*
    Runs one phase. Threads are used only when there is enough work:
    the evaluation phase is skipped for threading at low precision, where
    the cost of waking the threads exceeds the cost of the terms (the plan
    phase is cheap and only threaded at high precision).
*/
static void
afe_run_phase(afe_work_struct * W, int phase, slong thread_limit)
{
    slong i, num_workers;
    thread_pool_handle * handles;

    W->phase = phase;
    W->next = 0;
    W->lo = 0;

    thread_limit = FLINT_MIN(thread_limit, W->num);

    W->threaded = 0;

    if (thread_limit <= 1)
    {
        if (phase == 1)
        {
            W->nslots = 1;
            W->slots = _acb_vec_init(W->len);
            W->done = flint_calloc(1, sizeof(char));
        }

        afe_worker(W);
    }
    else
    {
        num_workers = flint_request_threads(&handles, thread_limit);
        W->threaded = (num_workers > 0);

        if (phase == 1)
        {
            /* a few slots per thread, so that threads rarely wait */
            W->nslots = 4 * (num_workers + 1);
            W->slots = _acb_vec_init(W->nslots * W->len);
            W->done = flint_calloc(W->nslots, sizeof(char));
        }

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, afe_worker, W);

        afe_worker(W);

        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        flint_give_back_threads(handles, num_workers);
    }

    if (phase == 1)
    {
        FLINT_ASSERT(W->lo == W->num);
        _acb_vec_clear(W->slots, W->nslots * W->len);
        flint_free(W->done);
    }
}

static void
afe_sums_compute(afe_sum_struct ** S, slong nsums)
{
    afe_work_struct W;
    slong i, maxN, threads, maxprec;

    maxN = 0;
    maxprec = 0;
    for (i = 0; i < nsums; i++)
    {
        W.S[i] = S[i];
        maxN = FLINT_MAX(maxN, S[i]->N);
        maxprec = FLINT_MAX(maxprec, S[i]->prec);
    }

    W.nsums = nsums;
    W.num = (maxN - 1) * nsums;
    W.len = S[0]->len;

    threads = flint_get_num_threads();

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&W.mutex, NULL);
    pthread_cond_init(&W.cond, NULL);
#endif

    afe_run_phase(&W, 0, (maxprec * S[0]->len >= 20000) ? threads : 1);

    for (i = 0; i < nsums; i++)
        afe_sum_gamma(S[i]);

    afe_run_phase(&W, 1, (maxprec * S[0]->len >= 1000) ? threads : 1);

#if FLINT_USES_PTHREAD
    pthread_cond_destroy(&W.cond);
    pthread_mutex_destroy(&W.mutex);
#endif
}

/* chi must be primitive */
static void
_acb_dirichlet_l_fmpq_afe_jet_primitive(acb_ptr res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec)
{
    arb_t t, u;
    acb_t w;
    acb_ptr S1, S2;
    arb_ptr R, V, W;
    fmpq_t s2, half;
    mag_t tol1, tol2;
    double ds, m1, m2, m2pre;
    slong origprec, prec1, prec2, j;
    ulong q;
    int parity, s_is_half, s_is_one;

    if (len <= 0)
        return;

    q = (G == NULL) ? 1 : G->q;
    parity = (G == NULL) ? 0 : dirichlet_parity_char(G, chi);

    if (q != 1)
        deflate = 0;

    s_is_one = fmpq_is_one(s);
    s_is_half = fmpz_is_one(fmpq_numref(s)) && fmpz_equal_ui(fmpq_denref(s), 2);

    if (q == 1 && s_is_one && !deflate)
    {
        _acb_vec_indeterminate(res, len);
        return;
    }


    origprec = prec;
    prec = prec * 1.001 + 2 * FLINT_BIT_COUNT(q) + ((len > 1) ? 2 * FLINT_BIT_COUNT(len) : 0);

    S1 = _acb_vec_init(len);
    S2 = _acb_vec_init(len);
    R = _arb_vec_init(len + 1);
    V = _arb_vec_init(len + 1);
    W = _arb_vec_init(len + 1);
    acb_init(w);
    arb_init(t);
    arb_init(u);
    fmpq_init(s2);
    fmpq_init(half);
    mag_init(tol1);
    mag_init(tol2);

    fmpq_set_si(half, 1, 2);

    ds = fmpq_get_d(s);

    m1 = log_gamma_upper_approx(0.5 * (ds + parity), PI / q) * INV_LOG2;
    m2 = log_gamma_upper_approx(0.5 * (1.0 - ds + parity), PI / q) * INV_LOG2;

    m2pre = (ds - 0.5) * log(PI / q) * INV_LOG2;

    mag_one(tol1);
    mag_mul_2exp_si(tol1, tol1, FLINT_MAX(m1, m2 + m2pre) - prec);

    mag_mul_2exp_si(tol2, tol1, -m2pre);

    prec1 = prec - (FLINT_MAX(m1, m2 + m2pre) - m1);
    prec1 = FLINT_MAX(prec1, 32);

    prec2 = prec - (FLINT_MAX(m1, m2 + m2pre) - (m2 + m2pre));
    prec2 = FLINT_MAX(prec2, 32);

    /* S1 = S(s + x) and S2 = S(1 - s - x), with the terms of both sums
       computed together (in parallel when threads are available);
       S(1/2 + y) = S1 when s = 1/2 */
    {
        afe_sum_struct sum1, sum2;
        afe_sum_struct * sums[2];

        afe_sum_init(&sum1, s, G, chi, tol1, len, prec1);
        sums[0] = &sum1;

        if (!s_is_half)
        {
            fmpq_sub_ui(s2, s, 1);
            fmpq_neg(s2, s2);
            afe_sum_init(&sum2, s2, G, chi, tol2, len, prec2);
            sums[1] = &sum2;
        }

        afe_sums_compute(sums, s_is_half ? 1 : 2);

        afe_sum_finish(S1, &sum1);
        afe_sum_clear(&sum1);

        if (s_is_half)
        {
            _acb_vec_set(S2, S1, len);
        }
        else
        {
            afe_sum_finish(S2, &sum2);
            afe_sum_clear(&sum2);
        }
    }

    for (j = 0; j < len; j++)
    {
        acb_conj(S2 + j, S2 + j);
        if (j & 1)
            acb_neg(S2 + j, S2 + j);
    }

    /* times w (pi/q)^(s-1/2+x) */
    arb_const_pi(t, prec);
    arb_div_ui(t, t, q, prec);
    arb_log(u, t, prec);
    _arb_exp_jet(R, u, len, prec);
    fmpq_sub(s2, s, half);
    arb_pow_fmpq(t, t, s2, prec);
    _arb_vec_scalar_mul(R, R, len, t, prec);
    _acb_vec_mullow_arb_vec(S2, S2, R, len, prec);

    if (q != 1)
    {
        acb_dirichlet_root_number2(w, G, chi, prec);
        _acb_vec_scalar_mul(S2, S2, len, w, prec);
    }

    _acb_vec_add(S1, S1, S2, len, prec);

    /* divide by Gamma((s+parity)/2 + x/2) */
    fmpq_add_ui(s2, s, parity);
    fmpq_div_2exp(s2, s2, 1);
    _arb_rgamma_fmpq_jet_half(V, s2, len, prec);
    _acb_vec_mullow_arb_vec(res, S1, V, len, prec);

    /* Add the polar part pi^((s+x)/2) / ((s+x) (s+x-1) Gamma((s+x)/2))
       = f(x) / (s+x-1),  f(x) = pi^((s+x)/2) / (2 Gamma((s+x)/2 + 1)),
       where the factor 1/(s+x) has been absorbed into the gamma function
       (dividing the series of 1/Gamma((s+x)/2) by s+x would lose about
       len log2(1/|s|) bits). With deflation, we compute
       (f(x) - 1) / (s+x-1) which is regular at s+x = 1 since f = 1 there. */
    if (q == 1)
    {
        slong wp, plen;

        wp = prec;

        /* the series division by (s-1+x) amplifies errors by |s-1|^(-len)
           when deflating; compensate (this part is cheap) */
        if (deflate && !s_is_one)
        {
            fmpq_sub_ui(s2, s, 1);
            wp += 10 + (slong) (len * FLINT_MAX(0.0, -log(fabs(fmpq_get_d(s2))) * INV_LOG2));
        }

        /* one more coefficient is needed for the deflation at s = 1 */
        plen = s_is_one ? len + 1 : len;

        /* R = pi^((s+x)/2) / 2 */
        arb_const_pi(t, wp);
        arb_log(u, t, wp);
        arb_mul_2exp_si(u, u, -1);
        _arb_exp_jet(R, u, plen, wp);
        fmpq_div_2exp(s2, s, 1);
        arb_pow_fmpq(t, t, s2, wp);
        arb_mul_2exp_si(t, t, -1);
        _arb_vec_scalar_mul(R, R, plen, t, wp);

        /* V = 1/Gamma(s/2 + 1 + x/2) */
        fmpq_add_ui(s2, s2, 1);
        _arb_rgamma_fmpq_jet_half(V, s2, plen, wp);

        /* W = f */
        _arb_poly_mullow(W, R, plen, V, plen, plen, wp);

        if (s_is_one)
        {
            /* deflated: (f(x) - 1) / x; the constant term is exactly 1 */
            for (j = 0; j < len; j++)
                arb_swap(W + j, W + j + 1);
        }
        else
        {
            if (deflate)
                arb_sub_ui(W, W, 1, wp);

            fmpq_sub_ui(s2, s, 1);
            arb_set_fmpq(R, s2, wp);
            arb_one(R + 1);
            _arb_poly_div_series(V, W, len, R, 2, len, wp);
            _arb_vec_swap(W, V, len);
        }

        for (j = 0; j < len; j++)
            arb_add(acb_realref(res + j), acb_realref(res + j), W + j, prec);
    }

    _acb_vec_set_round(res, res, len, origprec);

    _acb_vec_clear(S1, len);
    _acb_vec_clear(S2, len);
    _arb_vec_clear(R, len + 1);
    _arb_vec_clear(V, len + 1);
    _arb_vec_clear(W, len + 1);
    acb_clear(w);
    arb_clear(t);
    arb_clear(u);
    fmpq_clear(s2);
    fmpq_clear(half);
    mag_clear(tol1);
    mag_clear(tol2);
}

/*
    An imprimitive character chi modulo q is reduced to the primitive
    character chi0 inducing it, using
    L(s+x, chi) = L(s+x, chi0) prod_{p | q} (1 - chi0(p) p^(-s-x)).
    Deflation only has an effect for q = 1; with a principal character
    modulo q > 1 the output is indeterminate when deflating.
*/
void
acb_dirichlet_l_fmpq_afe_jet(acb_ptr res, const fmpq_t s, const dirichlet_group_t G,
    const dirichlet_char_t chi, int deflate, slong len, slong prec)
{
    ulong q, cond;

    if (len <= 0)
        return;

    q = (G == NULL) ? 1 : G->q;
    cond = (q == 1) ? 1 : dirichlet_conductor_char(G, chi);

    if (cond == q)
    {
        _acb_dirichlet_l_fmpq_afe_jet_primitive(res, s, G, chi, deflate, len, prec);
    }
    else if (deflate && cond == 1)
    {
        _acb_vec_indeterminate(res, len);
    }
    else
    {
        dirichlet_group_t G0;
        dirichlet_char_t chi0;
        n_factor_t fac;
        acb_ptr E, F, T;
        arb_t t, u;
        acb_t c;
        fmpq_t ms;
        slong k, j, wp;

        wp = prec + 10;

        if (cond != 1)
        {
            dirichlet_subgroup_init(G0, G, cond);
            dirichlet_char_init(chi0, G0);
            dirichlet_char_lower(chi0, G0, chi, G);
        }

        E = _acb_vec_init(len);
        F = _acb_vec_init(len);
        T = _acb_vec_init(len);
        arb_init(t);
        arb_init(u);
        acb_init(c);
        fmpq_init(ms);

        if (cond == 1)
            _acb_dirichlet_l_fmpq_afe_jet_primitive(T, s, NULL, NULL, 0, len, wp);
        else
            _acb_dirichlet_l_fmpq_afe_jet_primitive(T, s, G0, chi0, 0, len, wp);

        /* Euler factors 1 - chi0(p) p^(-s) exp(-x log(p)) */
        acb_one(E);
        fmpq_neg(ms, s);
        n_factor_init(&fac);
        n_factor(&fac, q, 1);

        for (k = 0; k < fac.num && acb_is_finite(T); k++)
        {
            ulong p = fac.p[k];

            if (cond % p == 0)
                continue;

            if (cond == 1)
                acb_one(c);
            else
                acb_dirichlet_chi(c, G0, chi0, p, wp);

            arb_set_ui(t, p);
            arb_pow_fmpq(t, t, ms, wp);
            acb_mul_arb(c, c, t, wp);

            arb_log_ui(u, p, wp);
            arb_neg(u, u);

            /* F = -c exp(-x log p) + 1 */
            acb_neg(F, c);
            for (j = 1; j < len; j++)
            {
                acb_mul_arb(F + j, F + j - 1, u, wp);
                acb_div_ui(F + j, F + j, j, wp);
            }
            acb_add_ui(F, F, 1, wp);

            _acb_poly_mullow(res, E, len, F, len, len, wp);
            _acb_vec_swap(E, res, len);
        }

        _acb_poly_mullow(res, T, len, E, len, len, wp);

        if (dirichlet_char_is_real(G, chi))
            for (j = 0; j < len; j++)
                arb_zero(acb_imagref(res + j));

        _acb_vec_set_round(res, res, len, prec);

        _acb_vec_clear(E, len);
        _acb_vec_clear(F, len);
        _acb_vec_clear(T, len);
        arb_clear(t);
        arb_clear(u);
        acb_clear(c);
        fmpq_clear(ms);

        if (cond != 1)
        {
            dirichlet_char_clear(chi0);
            dirichlet_group_clear(G0);
        }
    }
}
