/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "dfloat.h"
#include "impl.h"

/* Tables of F(t) = e^t t^-a Gamma(a, t) = int_0^oo (1 + v)^(a-1) e^(-tv) dv
   (real a) for the double-precision lattice sums, which then need only
   one exponential e^-t per point (computed with the vector kernels of
   dfloat) and a short polynomial evaluation.

   F satisfies t F' = (t - a) F - 1, so its Taylor coefficients at t0 obey
       t0 (k + 1) c_{k+1} = (t0 - a - k) c_k + c_{k-1} - [k = 0].
   The values F(t0) at the interval centers are propagated downward from
   one evaluation at the top (_acb_mzeta_gamma_upper_scaled): downward is
   the stable direction, since the homogeneous solution e^t t^-a decays
   relative to F.  Interval arithmetic through the recurrence would
   grossly overestimate the effect of the error r0 of F(t0); instead the
   coefficients are computed from the exact midpoint and r0 is propagated
   exactly: it becomes r0 h(t)/h(t0), h = e^t t^-a.

   Coefficient bounds: |c_k| <= tau^-(k+1), tau = t0 - max(0, a - 1) (from
   |(1 + v)^(a-1)| <= e^(max(0, a-1) v)), used on intervals of half-width
   tau/16 when max(0, a - 1) <= t0/4; otherwise the Cauchy bound
   |c_k| <= B (t0/2)^-k with B = sup_{Re t >= t0/2} |F| <= e^(t0/2)
   (t0/2)^-a Gamma(a) (a > 0), on intervals of half-width t0/32.

   The evaluation (_acb_mzeta_ftab_eval) uses Horner's rule in plain
   double (n = 1, with fma) or double-double (n = 2) arithmetic with the
   a priori bound sum_k (k + 1) |c_k| rho^k times the relative error of a
   Horner step. */

#define FT_R 16.0
#define FT_RC 32.0

/* the bucket of a positive double: its top 16 bits (exponent and 4 bits
   of mantissa), increasing with t; and the largest double of a bucket */
static slong
ft_key(double t)
{
    uint64_t u;
    memcpy(&u, &t, sizeof(u));
    return (slong) (u >> 48);
}

static double
ft_key_top(slong b)
{
    uint64_t u = ((uint64_t) b << 48) | ((UINT64_C(1) << 48) - 1);
    double t;
    memcpy(&t, &u, sizeof(t));
    return t;
}

static double
ft_rho(double t, double s, int * regime)
{
    if (s <= 0.25 * t)
    {
        *regime = 0;
        return (t - s) / FT_R;
    }
    *regime = 1;
    return t / FT_RC;
}

/* Taylor coefficients c_0..c_K of F at t0 from c0 = F(t0) */
static void
ft_coeffs(arb_ptr c, slong K, const arb_t c0, const arb_t t0, double a, slong prec)
{
    slong k;
    arb_t u, v;

    arb_init(u);
    arb_init(v);
    arb_set(c, c0);
    for (k = 0; k < K; k++)
    {
        /* a + k exactly (a may have low-order bits) */
        arb_set_d(u, a);
        arb_add_ui(u, u, k, prec);
        arb_sub(u, t0, u, prec);
        arb_mul(v, u, c + k, prec);
        if (k > 0)
            arb_add(v, v, c + k - 1, prec);
        else
            arb_sub_ui(v, v, 1, prec);
        arb_mul_ui(u, t0, k + 1, prec);
        arb_div(c + k + 1, v, u, prec);
    }
    arb_clear(u);
    arb_clear(v);
}

void
_acb_mzeta_ftab_init(_acb_mzeta_ftab_struct * T, double a, int n, double tmin, double tmax)
{
    slong prec = 64 * n + 48, K, K2, k, j, alloc = 32, i;
    arb_t t0, f, h, hp, tmp;
    arb_ptr c;
    slong bprec = 64;
    double t0d, s;
    int regime, unusable;

    K = (n == 1) ? 15 : 28;
    T->K1 = (n == 1) ? K + 1 : 14;
    K2 = (n == 1) ? 32 : 60;

    T->a = a;
    T->n = n;
    T->K = K;
    T->tmin_req = tmin;
    T->tmax_req = tmax;
    /* the callers check the range; a NaN or out-of-range request gives a
       small table that is marked unusable (all points get radius inf)
       instead of looping */
    unusable = !(tmin > 0.0 && tmin <= 1e15 && tmax <= 1e15 && fabs(a) <= 1e6);
    if (unusable)
    {
        a = T->a = 0.5;
        tmin = 1.0;
        tmax = 2.0;
    }
    tmax = FLINT_MAX(tmax, tmin * 1.5);
    s = FLINT_MAX(0.0, a - 1.0);
    T->t0 = flint_malloc(alloc * sizeof(double));
    T->rho = flint_malloc(alloc * sizeof(double));
    T->E = flint_malloc(alloc * sizeof(double));
    T->D = flint_malloc(alloc * sizeof(double));
    T->c = flint_malloc(alloc * (K + 1) * n * sizeof(double));

    arb_init(t0); arb_init(f); arb_init(h); arb_init(hp); arb_init(tmp);
    c = _arb_vec_init(K2 + 1);

    /* F(tmax) = e^tmax G(tmax) */
    t0d = tmax;
    arb_set_d(t0, t0d);
    {
        acb_t aa, g;
        acb_init(aa);
        acb_init(g);
        acb_set_d(aa, a);
        _acb_mzeta_gamma_upper_scaled(g, aa, t0, NULL, prec);
        arb_exp(tmp, t0, prec);
        arb_mul(f, acb_realref(g), tmp, prec);
        acb_clear(aa);
        acb_clear(g);
    }

    j = 0;
    for (;;)
    {
        double rho0 = ft_rho(t0d, s, &regime), rho = 1.02 * rho0;
        mag_t m;
        arb_t tau, Rho, Bc, Mh, Mt, E, r0, trunc, D, z;

        mag_init(m);
        arb_init(Mt);
        arb_init(tau); arb_init(Rho); arb_init(Bc); arb_init(Mh); arb_init(E);
        arb_init(r0); arb_init(trunc); arb_init(D); arb_init(z);

        if (j == alloc)
        {
            alloc *= 2;
            T->t0 = flint_realloc(T->t0, alloc * sizeof(double));
            T->rho = flint_realloc(T->rho, alloc * sizeof(double));
            T->E = flint_realloc(T->E, alloc * sizeof(double));
            T->D = flint_realloc(T->D, alloc * sizeof(double));
            T->c = flint_realloc(T->c, alloc * (K + 1) * n * sizeof(double));
        }

        /* the error bounds are computed in ball arithmetic from the
           (exact) doubles t0, rho, a and the coefficients */
        arb_set_d(Rho, rho);
        if (regime == 0)
        {
            /* tau = t0 - max(0, a - 1) */
            arb_set_d(tau, a);
            arb_sub_ui(tau, tau, 1, bprec);
            if (!arb_is_positive(tau))
                arb_zero(tau);
            arb_sub(tau, t0, tau, bprec);
        }
        else
        {
            /* B = e^rc rc^-a Gamma(a) >= sup_{Re t >= rc} |F|, rc = t0/2 */
            arb_mul_2exp_si(tau, t0, -1);
            arb_exp(Bc, tau, bprec);
            arb_set_d(z, -a);
            arb_pow(z, tau, z, bprec);
            arb_mul(Bc, Bc, z, bprec);
            arb_set_d(z, a);
            arb_gamma(z, z, bprec);
            arb_mul(Bc, Bc, z, bprec);
        }

        /* the coefficients from the exact midpoint of F(t0); its radius r0
           enters through the homogeneous solution: r0 h(t)/h(t0) with
           h(t)/h(t0) <= amp on the interval */
        arf_set_mag(arb_midref(r0), arb_radref(f));
        arb_set_arf(tmp, arb_midref(f));
        ft_coeffs(c, K2, tmp, t0, a, prec);

        arb_zero(Mh);
        arb_zero(E);
        for (k = 0; k <= K; k++)
        {
            double ck0 = 0.0;
            arb_set(tmp, c + k);
            for (i = 0; i < n; i++)
            {
                double ci = arf_get_d(arb_midref(tmp), ARF_RND_NEAR);
                if (i == 0)
                    ck0 = ci;
                arb_set_d(z, ci);
                arb_sub(tmp, tmp, z, prec);
                T->c[(j * (K + 1) + k) * n + i] = ci;
            }
            /* sum |c_k - stored| rho^k and sum (k+1) |c_k| rho^k */
            arb_pow_ui(z, Rho, k, bprec);
            arb_get_mag(m, tmp);
            arb_set_d(trunc, 0.0);
            arf_set_mag(arb_midref(trunc), m);
            arb_addmul(E, trunc, z, bprec);
            arb_set_d(trunc, fabs(ck0));
            arb_mul_ui(trunc, trunc, k + 1, bprec);
            arb_addmul(Mh, trunc, z, bprec);
            if (k >= T->K1)
                arb_addmul(Mt, trunc, z, bprec);
        }

        /* truncation after degree K, and the derivative bound D */
        if (regime == 0)
        {
            arb_div(z, Rho, tau, bprec);
            arb_pow_ui(trunc, z, K + 1, bprec);
            arb_sub(z, tau, Rho, bprec);
            arb_div(trunc, trunc, z, bprec);
            arb_sqr(D, z, bprec);
            arb_inv(D, D, bprec);
        }
        else
        {
            arb_div(z, Rho, tau, bprec);              /* q = rho / rc */
            arb_pow_ui(trunc, z, K + 1, bprec);
            arb_mul(trunc, trunc, Bc, bprec);
            arb_sub_ui(z, z, 1, bprec);
            arb_neg(z, z);
            arb_div(trunc, trunc, z, bprec);
            arb_sub(z, tau, Rho, bprec);              /* rc - rho */
            arb_sqr(z, z, bprec);
            arb_mul(D, Bc, tau, bprec);
            arb_div(D, D, z, bprec);
        }
        arb_add(E, E, trunc, bprec);

        /* Horner rounding: 2^-52 (n = 1: fma steps) resp. 2^-102 (n = 2:
           double-double product and sum, 16 u^2) times sum (k+1) |c_k| rho^k;
           for n = 2 the terms of degree >= K1 are summed in plain double
           with the high parts of c_k and h: the dropped low parts (u |c_k|,
           (k - K1) u |h|^k) and the rounding ((k - K1 + 1) u) give at most
           2^-51 sum_{k >= K1} (k + 1) |c_k| rho^k */
        arb_mul_2exp_si(Mh, Mh, (n == 1) ? -52 : -102);
        arb_add(E, E, Mh, bprec);
        if (n == 2)
        {
            arb_mul_2exp_si(Mt, Mt, -51);
            arb_add(E, E, Mt, bprec);
        }

        /* r0 amp: the polynomial is the degree-K truncation of the
           perturbed series, off by r0 trunc_K(h(t0 + y)/h(t0)),
           h(t) = e^t t^-a; e^rho (1 - rho/t0)^-|a| majorizes
           e^y (1 + y/t0)^-a for either sign of a, since
           |binom(-a, k)| <= (|a|)_k / k! */
        arb_exp(trunc, Rho, bprec);
        arb_sub(z, t0, Rho, bprec);
        arb_div(z, t0, z, bprec);
        {
            arb_t ea;
            arb_init(ea);
            arb_set_d(ea, fabs(a));
            arb_pow(z, z, ea, bprec);
            arb_clear(ea);
        }
        arb_mul(trunc, trunc, z, bprec);
        arb_addmul(E, r0, trunc, bprec);

        arb_get_mag(m, E);
        T->E[j] = mag_get_d(m);
        arb_get_mag(m, D);
        T->D[j] = mag_get_d(m);
        T->t0[j] = t0d;
        T->rho[j] = rho;
        j++;

        if (!arb_is_finite(E) || !arb_is_finite(D) || !(T->E[j - 1] < 1.0))
        {
            /* should not happen; make the table unusable rather than wrong */
            T->E[j - 1] = INFINITY;
        }

        if (t0d - rho <= tmin)
        {
            mag_clear(m);
            arb_clear(Mt);
            arb_clear(tau); arb_clear(Rho); arb_clear(Bc); arb_clear(Mh); arb_clear(E);
            arb_clear(r0); arb_clear(trunc); arb_clear(D); arb_clear(z);
            break;
        }

        /* the next (lower) center: t1 + rho0(t1) = t0 - rho0(t0) */
        {
            double t1, hd;
            int reg1;
            t1 = (t0d - rho0 + s / FT_R) / (1.0 + 1.0 / FT_R);
            ft_rho(t1, s, &reg1);
            if (reg1 == 1)
                t1 = (t0d - rho0) / (1.0 + 1.0 / FT_RC);
            hd = t1 - t0d;
            arb_set_d(h, hd);
            arb_zero(f);
            arb_one(hp);
            for (k = 0; k <= K2; k++)
            {
                arb_addmul(f, c + k, hp, prec);
                arb_mul(hp, hp, h, prec);
            }
            /* truncation: |h|/tau resp. |h|/rc is at most 1/8 by the
               choice of t1; checked */
            arb_abs(h, h);
            arb_div(z, h, tau, bprec);
            arb_sub_ui(tmp, z, 1, bprec);
            if (!arb_is_negative(tmp))
                arb_indeterminate(z);
            arb_pow_ui(trunc, z, K2 + 1, bprec);
            arb_sub_ui(z, z, 1, bprec);
            arb_neg(z, z);
            arb_div(trunc, trunc, z, bprec);
            if (regime == 0)
                arb_div(trunc, trunc, tau, bprec);
            else
                arb_mul(trunc, trunc, Bc, bprec);
            /* the propagated r0: r0 times the degree-K2 truncation of
               g(y) = h(t0 + y)/h(t0) = e^y (1 + y/t0)^-a at y = hd, which
               is g(hd) = e^hd (t1/t0)^-a (< 1 going down: the stable
               direction) plus the tail of g after degree K2; the
               coefficients of g are majorized by those of
               e^x (1 - x/t0)^-|a|, whose maximum on |x| = t0/2 is
               M = e^(t0/2) 2^|a|, so the tail is at most
               M q^(K2+1) / (1 - q), q = 2 |hd| / t0 <= 1/4 */
            {
                arb_t ea, q;
                arb_init(ea);
                arb_init(q);
                /* g(hd) */
                arb_set_d(z, t1);
                arb_div(z, z, t0, bprec);
                arb_set_d(ea, -a);
                arb_pow(z, z, ea, bprec);
                arb_set_d(D, hd);
                arb_exp(D, D, bprec);
                arb_mul(D, D, z, bprec);
                /* the tail */
                arb_mul_2exp_si(q, h, 1);
                arb_div(q, q, t0, bprec);
                arb_sub_ui(tmp, q, 1, bprec);
                if (!arb_is_negative(tmp))
                    arb_indeterminate(q);
                arb_pow_ui(z, q, K2 + 1, bprec);
                arb_sub_ui(q, q, 1, bprec);
                arb_neg(q, q);
                arb_div(z, z, q, bprec);
                arb_mul_2exp_si(q, t0, -1);
                arb_exp(q, q, bprec);
                arb_mul(z, z, q, bprec);
                arb_set_d(ea, fabs(a));
                arb_set_ui(q, 2);
                arb_pow(q, q, ea, bprec);
                arb_mul(z, z, q, bprec);
                arb_add(D, D, z, bprec);
                arb_clear(ea);
                arb_clear(q);
            }
            arb_addmul(trunc, r0, D, bprec);
            arb_get_mag(m, trunc);
            arb_add_error_mag(f, m);
            t0d = t1;
            arb_set_d(t0, t0d);
        }

        mag_clear(m);
        arb_clear(Mt);
        arb_clear(tau); arb_clear(Rho); arb_clear(Bc); arb_clear(Mh); arb_clear(E);
        arb_clear(r0); arb_clear(trunc); arb_clear(D); arb_clear(z);
    }

    T->len = j;
    if (unusable)
        for (k = 0; k < j; k++)
            T->E[k] = INFINITY;
    T->thi = T->t0[0] + T->rho[0];
    T->tlo = T->t0[j - 1] - T->rho[j - 1];

    /* lookup table: for each bucket of t (exponent and the top 4 bits of
       the mantissa, so 16 buckets per octave), the first interval (largest
       center) whose lower end lies below the top of the bucket; the
       evaluation then only walks forward a step or two */
    {
        slong b, jj = 0;
        T->kmin = ft_key(T->tlo);
        T->kmax = ft_key(T->thi);
        T->lut = flint_malloc((T->kmax - T->kmin + 1) * sizeof(slong));
        for (b = T->kmax; b >= T->kmin; b--)
        {
            double top = ft_key_top(b);
            while (jj < T->len - 1 && T->t0[jj] - T->rho[jj] > top)
                jj++;
            T->lut[b - T->kmin] = jj;
        }
    }

    arb_clear(t0); arb_clear(f); arb_clear(h); arb_clear(hp); arb_clear(tmp);
    _arb_vec_clear(c, K2 + 1);
}

void
_acb_mzeta_ftab_clear(_acb_mzeta_ftab_struct * T)
{
    flint_free(T->t0);
    flint_free(T->rho);
    flint_free(T->E);
    flint_free(T->D);
    flint_free(T->c);
    flint_free(T->lut);
}

static void
ftab_clear_obj(void * obj)
{
    _acb_mzeta_ftab_clear(obj);
}

typedef struct
{
    double a, tmin, tmax;
    int n;
}
ftab_key;

static int
ftab_match(const void * obj, const void * key)
{
    const _acb_mzeta_ftab_struct * T = obj;
    const ftab_key * k = key;
    return T->a == k->a && T->n == k->n && T->tmin_req <= k->tmin && T->tmax_req >= k->tmax;
}

const _acb_mzeta_ftab_struct *
_acb_mzeta_ftab_get(acb_mzeta_cache_t cache, double a, int n, double tmin, double tmax)
{
    _acb_mzeta_ftab_struct * T;
    ftab_key key;

    key.a = a;
    key.n = n;
    key.tmin = tmin;
    key.tmax = tmax;
    T = _acb_mzeta_cache_find(cache, _ACB_MZETA_CACHE_FTAB, ftab_match, &key);
    if (T == NULL)
    {
        T = flint_malloc(sizeof(_acb_mzeta_ftab_struct));
        _acb_mzeta_ftab_init(T, a, n, tmin, tmax);
        _acb_mzeta_cache_insert(cache, _ACB_MZETA_CACHE_FTAB, T, ftab_clear_obj);
    }
    return T;
}

/* double-double arithmetic for the n = 2 kernel (fma for every product
   whose rounding matters, so that contraction cannot change the result) */
#define DD_MUL(zh, zl, xh, xl, yh, yl) \
    do { \
        double _ph, _pl; \
        DFLOAT_TWO_PROD(_ph, _pl, (xh), (yh)); \
        _pl = fma((xh), (yl), _pl); \
        _pl = fma((xl), (yh), _pl); \
        DFLOAT_FAST_TWO_SUM((zh), (zl), _ph, _pl); \
    } while (0)

#define DD_ADD(zh, zl, xh, xl, yh, yl) \
    do { \
        double _s, _e, _t, _f; \
        DFLOAT_TWO_SUM(_s, _e, (xh), (yh)); \
        DFLOAT_TWO_SUM(_t, _f, (xl), (yl)); \
        _e += _t; \
        DFLOAT_FAST_TWO_SUM(_s, _e, _s, _e); \
        _e += _f; \
        DFLOAT_FAST_TWO_SUM((zh), (zl), _s, _e); \
    } while (0)

#define FT_ILP 4

void
_acb_mzeta_ftab_eval(double * res, const double * t, slong len, const _acb_mzeta_ftab_struct * T)
{
    int n = T->n;
    slong K = T->K, K1 = T->K1, i, j = 0, k, g, stride = n + 1;
    slong jj[FT_ILP], ii[FT_ILP];
    double hh[FT_ILP], hl[FT_ILP], hr[FT_ILP];
    int cnt;

    /* pass over the points, collecting groups of FT_ILP valid points whose
       independent Horner chains are interleaved (the evaluation of one
       chain is bound by the latency of its fma / double-double steps) */
    cnt = 0;
    for (i = 0; i <= len; i++)
    {
        if (i < len)
        {
            const double * ti = t + i * stride;
            double * ri = res + i * stride;
            double tm = ti[0], h, l = 0.0, r;

            if (!(tm >= T->tlo && tm <= T->thi && ti[n] < 1e-3))
            {
                for (k = 0; k < n; k++)
                    ri[k] = 0.0;
                ri[n] = INFINITY;
                continue;
            }

            /* centers decrease with j; start from the lookup table */
            j = T->lut[ft_key(tm) - T->kmin];
            while (j < T->len - 1 && tm < T->t0[j] - T->rho[j])
                j++;

            if (n == 1)
            {
                h = tm - T->t0[j];
                /* the rounding of h (at most ulp(h)/2) goes into the radius */
                r = ti[1] + fabs(h) * 0x1p-52;
            }
            else
            {
                /* h = (t_hi - t0) + t_lo as a double-double */
                double s2, e2;
                DFLOAT_TWO_SUM(h, l, tm, -T->t0[j]);
                DFLOAT_TWO_SUM(s2, e2, h, ti[1]);
                e2 += l;
                DFLOAT_FAST_TWO_SUM(h, l, s2, e2);
                r = ti[2] + fabs(h) * 0x1p-103;
            }

            if (!(fabs(h) + r <= T->rho[j] * (1.0 - 0x1p-50)))
            {
                for (k = 0; k < n; k++)
                    ri[k] = 0.0;
                ri[n] = INFINITY;
                continue;
            }

            ii[cnt] = i;
            jj[cnt] = j;
            hh[cnt] = h;
            hl[cnt] = l;
            hr[cnt] = r;
            cnt++;
            if (cnt < FT_ILP)
                continue;
        }

        if (cnt == 0)
            continue;

        /* pad the group with copies of its first point */
        for (g = cnt; g < FT_ILP; g++)
        {
            ii[g] = -1;
            jj[g] = jj[0];
            hh[g] = hh[0];
            hl[g] = hl[0];
            hr[g] = hr[0];
        }

        if (n == 1)
        {
            double p[FT_ILP];
            const double * c[FT_ILP];
            for (g = 0; g < FT_ILP; g++)
            {
                c[g] = T->c + jj[g] * (K + 1);
                p[g] = c[g][K];
            }
            for (k = K - 1; k >= 0; k--)
                for (g = 0; g < FT_ILP; g++)
                    p[g] = fma(p[g], hh[g], c[g][k]);
            for (g = 0; g < cnt; g++)
            {
                double * ri = res + ii[g] * stride;
                ri[0] = p[g];
                ri[1] = _dfloat_rad_finish(T->E[jj[g]] + _dfloat_rad_mul(T->D[jj[g]], hr[g]));
                if (!isfinite(p[g]))
                {
                    ri[0] = 0.0;
                    ri[1] = INFINITY;
                }
            }
        }
        else
        {
            double ph[FT_ILP], pl[FT_ILP];
            const double * c[FT_ILP];
            for (g = 0; g < FT_ILP; g++)
            {
                c[g] = T->c + jj[g] * (K + 1) * 2;
                ph[g] = c[g][2 * K];
                pl[g] = 0.0;
            }
            /* the terms of degree >= K1 are below 2^-56 of the value:
               plain double Horner with the high parts (the error bound of
               the table includes this) */
            for (k = K - 1; k >= K1; k--)
                for (g = 0; g < FT_ILP; g++)
                    ph[g] = fma(ph[g], hh[g], c[g][2 * k]);
            for (k = K1 - 1; k >= 0; k--)
                for (g = 0; g < FT_ILP; g++)
                {
                    DD_MUL(ph[g], pl[g], ph[g], pl[g], hh[g], hl[g]);
                    DD_ADD(ph[g], pl[g], ph[g], pl[g], c[g][2 * k], c[g][2 * k + 1]);
                }
            for (g = 0; g < cnt; g++)
            {
                double * ri = res + ii[g] * stride;
                ri[0] = ph[g];
                ri[1] = pl[g];
                ri[2] = _dfloat_rad_finish(T->E[jj[g]] + _dfloat_rad_mul(T->D[jj[g]], hr[g]));
                if (!isfinite(ph[g]) || !isfinite(pl[g]))
                {
                    ri[0] = ri[1] = 0.0;
                    ri[2] = INFINITY;
                }
            }
        }

        cnt = 0;
    }
}
