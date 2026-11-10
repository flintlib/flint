/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_mat.h"
#include "mag.h"
#include "mp_real.h"
#include "impl.h"

/* The dual terms of the theta splitting of Khinchin's constant
   (const_khinchin_afe.c):

       J_m = sqrt(pi) nu_(3/2)(c),  c = pi^2 (m + 1/2)^2,

   with the moments (u0 = 4^-a, E = e^(-c/u0), E(x) = sum x^n/(n n!))

       nu_t = int_0^u0 u^-t e^(-c/u) E(4u) du,
       Y_t  = int_0^u0 u^-t e^(-c/u) e^(4u) du,
       C_t  = int_0^u0 u^-t e^(-c/u) du,

   for which integration by parts gives

       c nu_(t+1) = E g1 u0^(1-t) - (1-t) nu_t - Y_t + C_t,
       c Y_(t+1)  = E e4 u0^(1-t) - (1-t) Y_t - 4 Y_(t-1),
       c C_(t+1)  = E u0^(1-t) - (1-t) C_t,

   g1 = E(4 u0), e4 = e^(4 u0), and -d/dc shifts t by one.  With the
   scaled state X_j = c^j (nu, Y, Y_(-1)/c, C) at t = 3/2 + j, J_m is the
   value at c of a holonomic function of order four whose only finite
   singularity is c = 0.  It is evaluated by the bit-burst algorithm:

   STARTING VALUES at a 64-bit dyadic c0 near c, from the shift series:
   with D_k = int_0^u0 u^(k - 1/2) e^(-c0/u) du = alpha_k Gb + beta_k E,
   Gb = D_(-1) = c0^(-1/2) Gamma(1/2, c0/u0) and
   (k + 1/2) D_k + c0 D_(k-1) = E u0^(k + 1/2), b_n = 4^n/n!,

       P_n = b_n alpha_(n-1),  R_n = b_n beta_(n-1),  T_n = b_n u0^(n+1/2),
       P_(n+1) = -4 c0 P_n / ((n+1)(n+1/2)),
       R_(n+1) = 4 (T_n - c0 R_n) / ((n+1)(n+1/2)),  T_(n+1) = 4 u0 T_n/(n+1),

   and seven sums (nu = Gb s0 + E s1, Y = Gb s2 + E s3,
   Y_(-1) = -c0 Gb s4 + E (s5 - c0 s6)) by one matrix binary splitting.
   Each D_k cancels by up to e^(c0/u0) between its two parts, but the
   weighted sums only by about e^(4 sqrt(c0)), the size of
   sum |b_n alpha_(n-1)| = sum (16 c0)^n/(2n)!, so that many guard bits
   suffice.

   TAYLOR STAGES from c0 to c1 = c0 + delta (c_j with 64 4^j fractional
   bits): Z_t(c0 + delta) = sum_l (-delta)^l / l! Z_(t+l)(c0), the
   coefficients from the recurrences above, as a matrix binary splitting
   with the forcing E c^j u0^-j (three constants g1, e4, 1) as the
   scalar-geometric forcing of _mp_real_linrec.  Short stages (late
   ones, with few terms of huge bit size) are summed directly.

   ERROR BOUNDS.  On (0, u0], 0 < E(4u) <= e^(4u) - 1 and 1 <= e^(4u),
   so 0 < nu_t, C_t <= Y_t <= e^(4 u0) C_t with
   C_s(c0) = c0^(1-s) Gamma(s - 1, c0/u0).  The Taylor remainder after L
   terms is bounded through Gamma(s - 1, x) <= Gamma(s - 1) by
   2 e^(4 u0) c0^(-1/2) rho^L / (1 - rho), rho = |delta|/c0, and through
   Gamma(a, x) <= x^(a-1) e^-x / (1 - (a-1)/x) for 1 <= a < 1 + x by a
   bound carrying e^(-c0/u0), so that late stages need only the
   significant bits of J_m (the smaller bound is used).  The starting
   values' tails are at most twice their first terms.  Finally the last
   stage evaluates J at a c_l within 2^-B of c, and |J'| is bounded on
   [c - 1, inf). */

#define NX 4    /* state */
#define NG 3    /* forcing constants */
#define IV 3    /* starting values: state */
#define IW 7    /* starting values: sums */

/* Taylor stages ************************************************************/

typedef struct
{
    fmpz_t C, D, dn, dd;    /* c0 = C/D (D = 2^e), delta = dn/dd (dd = 2^B) */
    ulong u0d;
    slong ns;
    int sel[NX];
}
_step_struct;

/* leaf l: r M_l with r = -delta / (c0 (l+1)), as the integer matrix
   -dn M' over q = 2 dd C (l+1), M' = 2 D M:
     nu:  [(2l+1)D, -2D, 0, 2D | 2D, 0, 0]
     Y:   [0, (2l+1)D, -8C, 0  | 0, 2D, 0]
     Y-1: [0, 2D, 0, 0         | 0, 0, 0]
     C:   [0, 0, 0, (2l+1)D    | 0, 0, 2D],
   forcing factor -delta/((l+1) u0): phi = -2 C dn u0d; the sums take
   the state before the step */
static void
_taylor_leaf(fmpz_mat_t T, fmpz_t phi, fmpz_t q, slong l, const void * vctx)
{
    const _step_struct * st = (const _step_struct *) vctx;
    fmpz_t a, b;
    slong i;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_mul(b, st->D, st->dn);
    fmpz_mul_si(a, b, -(2 * l + 1));
    fmpz_mul_si(b, b, -2);
    fmpz_set(fmpz_mat_entry(T, 0, 0), a);
    fmpz_neg(fmpz_mat_entry(T, 0, 1), b);
    fmpz_set(fmpz_mat_entry(T, 0, 3), b);
    fmpz_set(fmpz_mat_entry(T, 0, NX + 0), b);
    fmpz_set(fmpz_mat_entry(T, 1, 1), a);
    fmpz_mul(fmpz_mat_entry(T, 1, 2), st->C, st->dn);
    fmpz_mul_2exp(fmpz_mat_entry(T, 1, 2), fmpz_mat_entry(T, 1, 2), 3);
    fmpz_set(fmpz_mat_entry(T, 1, NX + 1), b);
    fmpz_set(fmpz_mat_entry(T, 2, 1), b);
    fmpz_set(fmpz_mat_entry(T, 3, 3), a);
    fmpz_set(fmpz_mat_entry(T, 3, NX + 2), b);
    fmpz_mul(phi, st->C, st->dn);
    fmpz_mul_ui(phi, phi, st->u0d);
    fmpz_mul_si(phi, phi, -2);
    fmpz_mul(q, st->dd, st->C);
    fmpz_mul_ui(q, q, 2 * (ulong) (l + 1));
    for (i = 0; i < st->ns; i++)
        fmpz_set(fmpz_mat_entry(T, NX + i, st->sel[i]), q);
    fmpz_clear(a);
    fmpz_clear(b);
}

/* ln2 of |x| for a positive fmpz (an estimate) */
static double
_fmpz_log2(const fmpz_t x)
{
    return fmpz_dlog(x) * 1.4426950408889634;
}

/* The Taylor remainder bound (see the top of the file) for L terms of
   the step from c0 = C 2^-e by delta = dn 2^-B. */
static void
_taylor_tail(mag_t res, const _step_struct * st, slong e, slong B, slong L,
    const _mp_real_khinchin_afe_struct * ctx)
{
    mag_t c0lo, rho, t, u, e4u, w;
    fmpz_t lm, z;
    slong a = ctx->a;

    mag_init(c0lo); mag_init(rho); mag_init(t); mag_init(u); mag_init(e4u); mag_init(w);
    fmpz_init(lm);
    fmpz_init(z);

    mag_set_fmpz_lower(c0lo, st->C);
    mag_mul_2exp_si(c0lo, c0lo, -e);
    fmpz_abs(z, st->dn);
    mag_set_fmpz(rho, z);
    mag_mul_2exp_si(rho, rho, -B);
    mag_div(rho, rho, c0lo);                        /* rho = |delta|/c0 */

    /* e^(4 u0) */
    mag_set_ui_2exp_si(e4u, 1, 2 - 2 * a);
    mag_exp(e4u, e4u);

    /* (a) 2 e^(4 u0) c0^(-1/2) rho^L / (1 - rho) */
    mag_one(t);
    mag_sub_lower(t, t, rho);
    if (mag_is_zero(t))
    {
        mag_inf(res);
        goto cleanup;
    }
    mag_pow_ui(res, rho, L);
    mag_div(res, res, t);
    mag_rsqrt(u, c0lo);
    mag_mul(res, res, u);
    mag_mul(res, res, e4u);
    mag_mul_2exp_si(res, res, 1);

    /* (b) for x = c0/u0 > 4, lm = floor(x/2) >= L and |delta|/u0 <= (L+1)/2:
       2 [4 e^(4 u0) e^-x (|delta|/u0)^L / L! (u0^(1/2)/c0 + u0^(3/2)/c0^2)
          + e^(4 u0) c0^(-1/2) rho^(lm+1) / (1 - rho)] */
    fmpz_mul_2exp(lm, st->C, 2 * a);
    fmpz_fdiv_q_2exp(lm, lm, e + 1);                /* floor(x/2) */
    if (fmpz_cmp_ui(lm, 2) > 0 && fmpz_cmp_si(lm, L) >= 0 && fmpz_abs_fits_ui(lm))
    {
        /* |dn| 4^a 2 <= (L + 1) 2^B */
        fmpz_t l1, r1;
        fmpz_init(l1);
        fmpz_init(r1);
        fmpz_mul_2exp(l1, z, 2 * a + 1);
        fmpz_set_ui(r1, L + 1);
        fmpz_mul_2exp(r1, r1, B);
        if (fmpz_cmp(l1, r1) <= 0)
        {
            mag_t x, v, s;
            mag_init(x); mag_init(v); mag_init(s);

            /* 4 e^(4u0) e^-x (|delta| 4^a)^L / L! (2^-a/c0 + 2^-3a/c0^2) */
            mag_set(x, c0lo);
            mag_mul_2exp_si(x, x, 2 * a);
            mag_expinv(v, x);
            mag_set_fmpz(s, z);
            mag_mul_2exp_si(s, s, 2 * a - B);
            mag_pow_ui(s, s, L);
            mag_mul(v, v, s);
            mag_rfac_ui(s, L);
            mag_mul(v, v, s);
            mag_one(s);
            mag_mul_2exp_si(s, s, -a);
            mag_div(s, s, c0lo);
            mag_set(x, s);
            mag_mul_2exp_si(x, x, -2 * a);
            mag_div(x, x, c0lo);
            mag_add(s, s, x);
            mag_mul(v, v, s);
            mag_mul(v, v, e4u);
            mag_mul_2exp_si(v, v, 2);
            /* + e^(4u0) c0^(-1/2) rho^(lm+1) / (1 - rho) */
            mag_pow_ui(s, rho, fmpz_get_ui(lm) + 1);
            mag_div(s, s, t);
            mag_mul(s, s, u);
            mag_mul(s, s, e4u);
            mag_add(v, v, s);
            mag_mul_2exp_si(v, v, 1);
            if (mag_cmp(v, res) < 0)
                mag_swap(res, v);
            mag_clear(x); mag_clear(v); mag_clear(s);
        }
        fmpz_clear(l1);
        fmpz_clear(r1);
    }

cleanup:
    mag_clear(c0lo); mag_clear(rho); mag_clear(t); mag_clear(u); mag_clear(e4u); mag_clear(w);
    fmpz_clear(lm);
    fmpz_clear(z);
}

/* the number of Taylor terms (doubles: an estimate; the remainder is
   bounded rigorously by _taylor_tail) */
static slong
_taylor_terms(const _step_struct * st, slong e, slong B, slong abs_bits, double u0)
{
    double lc = _fmpz_log2(st->C) - e;
    double ld = (double) fmpz_bits(st->dn) - B;
    double lrho = ld - lc;
    slong L, L2;

    L = (slong) ((abs_bits + 2.5 * FLINT_MAX(lc, 0) + 48) / (-lrho)) + 2;

    /* the e^(-c0/u0)-aware bound usually needs far fewer terms when
       c0/u0 is large */
    {
        double x = pow(2.0, lc) / u0, lu = log2(u0);
        double lm = floor(x / 2), lim = -(double) abs_bits - 10;
        double cst = 3 + 4 * u0 * 1.4426950408889634 - x * 1.4426950408889634
            + FLINT_MAX(0.5 * lu - lc, 1.5 * lu - 2 * lc);
        double tl = 1 + 4 * u0 * 1.4426950408889634 - 0.5 * lc + (lm + 1) * lrho;

        if (tl < lim && lm > 2)
        {
            for (L2 = 2; L2 <= lm && L2 < L; L2++)
            {
                double hd = cst + L2 * (ld - lu) - lgamma(L2 + 1.0) * 1.4426950408889634;
                if (hd < lim && pow(2.0, ld - lu) <= (L2 + 1) / 2.0)
                    break;
            }
            if (L2 < L)
                L = L2;
        }
    }
    return FLINT_MAX(L, 2);
}

/* below this many terms a stage is summed directly */
#define TAYLOR_DIRECT_TERMS 100

/* One step from c0 = C 2^-e to c1 = C1 2^-B on the state X and E =
   e^(-c0/u0); last: only nu is needed.  wp: relative working precision
   (bits) of the state; abs_bits: target of the remainder. */
static void
_taylor_step(mp_real_ptr X, mp_real_t E, const fmpz_t C, slong e, const fmpz_t C1,
    slong B, slong wp, slong abs_bits, int last,
    const _mp_real_khinchin_afe_struct * ctx)
{
    _step_struct st;
    double u0 = ldexp(1.0, -2 * ctx->a), lc;
    slong L, i, W, n = mp_real_prec_bits(wp);
    mp_real_struct Acc[NX];
    mp_real_t F;
    mag_t tail;

    fmpz_init_set(st.C, C);
    fmpz_init(st.D);
    fmpz_one(st.D);
    fmpz_mul_2exp(st.D, st.D, e);
    fmpz_init(st.dn);
    fmpz_mul_2exp(st.dn, C, B - e);
    fmpz_sub(st.dn, C1, st.dn);
    fmpz_init(st.dd);
    fmpz_one(st.dd);
    fmpz_mul_2exp(st.dd, st.dd, B);
    st.u0d = ctx->u0d;
    if (last)
    {
        st.ns = 1;
        st.sel[0] = 0;
    }
    else
    {
        st.ns = NX;
        for (i = 0; i < NX; i++)
            st.sel[i] = (int) i;
    }

    for (i = 0; i < NX; i++)
        mp_real_init(Acc + i);
    mp_real_init(F);
    mag_init(tail);

    lc = _fmpz_log2(C) - e;
    L = _taylor_terms(&st, e, B, abs_bits, u0);
    _taylor_tail(tail, &st, e, B, L, ctx);

    /* F_0 = E u0^(1 - t0) = E / sqrt(u0) */
    mp_real_mul_2exp_si(F, E, ctx->a);

    if (L <= TAYLOR_DIRECT_TERMS)
    {
        /* few terms: directly on the state,
             X_(l+1) = r_l (M_l X_l + f_l (g1, e4, 0, 1)),  f_(l+1) = phi_l f_l,
             r_l = -delta/(c0 (l+1)),  phi_l = -delta/((l+1) u0),
             M_l = [(2l+1)/2, -1, 0, 1; 0, (2l+1)/2, -4 c0, 0; 0, 1, 0, 0;
                    0, 0, 0, (2l+1)/2],
           term l + 1 about |delta| max(1/u0, (l + 3/2)/c0)/(l + 1) times
           term l, so its working precision is lowered by the log of the
           product (the radii are propagated in any case) */
        mp_real_struct x[NX], y[NX];
        mp_real_t f, r0, ph0, rl, c0, t, h;
        double dlog = 0.0, c0d = pow(2.0, lc), ldel = (double) fmpz_bits(st.dn) - B;
        slong l;

        for (i = 0; i < NX; i++)
        {
            mp_real_init(x + i);
            mp_real_init(y + i);
            mp_real_set(x + i, X + i);
        }
        mp_real_init(f); mp_real_init(r0); mp_real_init(ph0); mp_real_init(rl);
        mp_real_init(c0); mp_real_init(t); mp_real_init(h);

        mp_real_set_fmpz(c0, C);
        mp_real_mul_2exp_si(c0, c0, -e);
        mp_real_set_fmpz(ph0, st.dn);
        mp_real_mul_2exp_si(ph0, ph0, -B);
        mp_real_neg(ph0, ph0);                         /* -delta */
        mp_real_div(r0, ph0, c0, n);                   /* -delta/c0 */
        mp_real_mul_2exp_si(ph0, ph0, 2 * ctx->a);     /* -delta/u0 */
        mp_real_set(f, F);

        for (l = 0; l < L; l++)
        {
            slong wl;

            dlog += ldel + log2(FLINT_MAX(1.0 / u0, (l + 1.5) / c0d)) - log2(l + 1.0);
            wl = mp_real_prec_bits(FLINT_MAX(64, wp + (slong) FLINT_MIN(0.0, dlog) + 128));
            for (i = 0; i < NX; i++)
                mp_real_add(Acc + i, Acc + i, x + i, n);
            if (l == L - 1)
                break;
            mp_real_set_ui(h, 2 * l + 1);
            mp_real_mul_2exp_si(h, h, -1);
            /* y = M_l x + f (g1, e4, 0, 1) */
            mp_real_mul(y + 0, x + 0, h, wl);
            mp_real_sub(y + 0, y + 0, x + 1, wl);
            mp_real_add(y + 0, y + 0, x + 3, wl);
            mp_real_mul(t, f, ctx->g1, wl);
            mp_real_add(y + 0, y + 0, t, wl);
            mp_real_mul(y + 1, x + 1, h, wl);
            mp_real_mul(t, x + 2, c0, wl);
            mp_real_mul_2exp_si(t, t, 2);
            mp_real_sub(y + 1, y + 1, t, wl);
            mp_real_mul(t, f, ctx->e4, wl);
            mp_real_add(y + 1, y + 1, t, wl);
            mp_real_set(y + 2, x + 1);
            mp_real_mul(y + 3, x + 3, h, wl);
            mp_real_add(y + 3, y + 3, f, wl);
            mp_real_div_ui(rl, r0, l + 1, wl);
            for (i = 0; i < NX; i++)
                mp_real_mul(x + i, y + i, rl, wl);
            mp_real_mul(f, f, ph0, wl);
            mp_real_div_ui(f, f, l + 1, wl);
        }

        for (i = 0; i < NX; i++)
        {
            mp_real_clear(x + i);
            mp_real_clear(y + i);
        }
        mp_real_clear(f); mp_real_clear(r0); mp_real_clear(ph0); mp_real_clear(rl);
        mp_real_clear(c0); mp_real_clear(t); mp_real_clear(h);
    }
    else
    {
        _mp_real_linrec_struct rec;
        mp_real_struct ghat[NG], Y[NX], Xs[NX];

        rec.d = NX;
        rec.s = st.ns;
        rec.g = NG;
        rec.leaf = _taylor_leaf;
        rec.bottom = NULL;
        rec.bottom_len = 0;
        rec.ctx = &st;

        for (i = 0; i < NG; i++)
            mp_real_init(ghat + i);
        for (i = 0; i < NX; i++)
        {
            mp_real_init(Y + i);
            mp_real_init(Xs + i);
        }
        mp_real_set(ghat + 0, ctx->g1);
        mp_real_set(ghat + 1, ctx->e4);
        mp_real_set_ui(ghat + 2, 1);

        /* the state from X, the sums (those selected) into Y */
        W = wp + 64 + 4 * (slong) FLINT_MAX(lc, 0);
        for (i = 0; i < NX; i++)
        {
            mp_real_set(Xs + i, X + i);
            mp_real_zero(Y + i);
        }
        _mp_real_linrec_eval(Xs, Y, F, &rec, 0, L, ghat, 1, mp_real_prec_bits(W));
        for (i = 0; i < st.ns; i++)
            mp_real_swap(Acc + st.sel[i], Y + i);

        for (i = 0; i < NX; i++)
        {
            mp_real_clear(Y + i);
            mp_real_clear(Xs + i);
        }
        for (i = 0; i < NG; i++)
            mp_real_clear(ghat + i);
    }

    for (i = 0; i < st.ns; i++)
        _mp_real_add_error_mag(Acc + st.sel[i], tail);

    if (!last)
    {
        mp_real_t r, t;

        mp_real_init(r);
        mp_real_init(t);
        /* the lagged component is scaled by 1/c: times c0 / c1 */
        mp_real_set_fmpz(r, C);
        mp_real_mul_2exp_si(r, r, B - e);
        mp_real_set_fmpz(t, C1);
        mp_real_div(r, r, t, n);
        mp_real_swap(X + 0, Acc + 0);
        mp_real_swap(X + 1, Acc + 1);
        mp_real_mul(X + 2, Acc + 2, r, n);
        mp_real_swap(X + 3, Acc + 3);
        /* E = e^(-c1/u0) = E e^(-delta/u0) */
        mp_real_set_fmpz(t, st.dn);
        mp_real_mul_2exp_si(t, t, 2 * ctx->a - B);
        mp_real_neg(t, t);
        mp_real_exp_bits(t, t, FLINT_BITS * n);
        mp_real_mul(E, E, t, n);
        mp_real_clear(r);
        mp_real_clear(t);
    }
    else
        mp_real_swap(X + 0, Acc + 0);

    for (i = 0; i < NX; i++)
        mp_real_clear(Acc + i);
    mp_real_clear(F);
    mag_clear(tail);
    fmpz_clear(st.C);
    fmpz_clear(st.D);
    fmpz_clear(st.dn);
    fmpz_clear(st.dd);
}

/* starting values **********************************************************/

typedef struct
{
    const fmpz * C;
    const fmpz * D;
    ulong u0d;
}
_init_struct;

/* rows 0..2: (P, R, T) -> next; rows 3..9: the sums s0..s6, all over
   the common q = (n+1)(2n+1) n' D u0d, n' = max(n, 1):
   the transition 1/((n+1)(n+1/2)) = 2/((n+1)(2n+1)), c0 = C/D, the
   weights 1/n, 1 and 2/(2n+1) (the last divides the transition's
   denominator already) */
static void
_init_leaf(fmpz_mat_t T, fmpz_t phi, fmpz_t q, slong n, const void * vctx)
{
    const _init_struct * c = (const _init_struct *) vctx;
    ulong nn = (n == 0) ? 1 : n;
    fmpz_t x;
    slong i, j;

    fmpz_init(x);
    fmpz_set_ui(q, n + 1);
    fmpz_mul_ui(q, q, 2 * (ulong) n + 1);
    fmpz_mul_ui(q, q, nn);
    fmpz_mul(q, q, c->D);

    /* P_(n+1) = -4 c P_n 2/((n+1)(2n+1)): -8 C nn */
    fmpz_mul_ui(fmpz_mat_entry(T, 0, 0), c->C, 8 * nn);
    fmpz_neg(fmpz_mat_entry(T, 0, 0), fmpz_mat_entry(T, 0, 0));
    /* R_(n+1) = 4 (T_n - c R_n) 2/((n+1)(2n+1)) */
    fmpz_mul_ui(fmpz_mat_entry(T, 1, 2), c->D, 8 * nn);
    fmpz_set(fmpz_mat_entry(T, 1, 1), fmpz_mat_entry(T, 0, 0));
    /* T_(n+1) = 4 T_n / (n+1): 4 (2n+1) nn D */
    fmpz_mul_ui(fmpz_mat_entry(T, 2, 2), c->D, 4 * nn);
    fmpz_mul_ui(fmpz_mat_entry(T, 2, 2), fmpz_mat_entry(T, 2, 2), 2 * (ulong) n + 1);
    /* sums */
    if (n >= 1)
    {
        fmpz_divexact_ui(x, q, n);
        fmpz_set(fmpz_mat_entry(T, IV + 0, 0), x);
        fmpz_set(fmpz_mat_entry(T, IV + 1, 1), x);
    }
    fmpz_set(fmpz_mat_entry(T, IV + 2, 0), q);
    fmpz_set(fmpz_mat_entry(T, IV + 3, 1), q);
    fmpz_divexact_ui(x, q, 2 * (ulong) n + 1);
    fmpz_mul_2exp(x, x, 1);
    fmpz_set(fmpz_mat_entry(T, IV + 4, 0), x);
    fmpz_set(fmpz_mat_entry(T, IV + 5, 2), x);
    fmpz_set(fmpz_mat_entry(T, IV + 6, 1), x);
    /* u0 = 1/u0d: T_(n+1) = 4 u0 T_n/(n+1), so everything but that entry
       is scaled by u0d, and q too */
    for (i = 0; i < IV + IW; i++)
        for (j = 0; j < IV; j++)
            if (!(i == 2 && j == 2))
                fmpz_mul_ui(fmpz_mat_entry(T, i, j), fmpz_mat_entry(T, i, j), c->u0d);
    fmpz_mul_ui(q, q, c->u0d);
    fmpz_clear(x);
}

/* the number of terms: b_n |alpha| ~ (4c)^n / (n!)^2 and b_n u0^n
   ~ (4 u0)^n / n! below 2^-wp (an estimate; the tail is bounded
   rigorously) */
static slong
_init_terms(double cd, double u0, slong wp)
{
    double ly = 0.0, lz = 0.0;
    slong N;

    for (N = 1; ; N++)
    {
        ly += log2(4 * cd + 1) - 2 * log2(N + 0.5);
        lz += 2 + log2(u0) - log2((double) N);
        if (N > 10 && ly + 10 < -(double) wp - 20 && lz + 10 < -(double) wp - 20)
            break;
    }
    N++;
    /* the conditions of the tail bound */
    while ((2.0 * N + 1) * (2.0 * N + 2) < 32 * cd + 32 || N + 1 < 8 * u0 + 1)
        N++;
    return N;
}

/* the tail bound after N terms, for every component:
     4 Gb |P_N| + 2 E T_N (1 / (u0 (N - 1/2)) + 1/c),
   |P_N| = (16 c)^N / (2N)!,  T_N = 4^N u0^(N+1/2) / N! */
static void
_init_tail(mag_t res, const fmpz_t C, slong e, const mp_real_t Gb,
    const mp_real_t E, slong N, slong a)
{
    mag_t c, cl, P, T, t, s;

    mag_init(c); mag_init(cl); mag_init(P); mag_init(T); mag_init(t); mag_init(s);
    mag_set_fmpz(c, C);
    mag_mul_2exp_si(c, c, -e);
    mag_set_fmpz_lower(cl, C);
    mag_mul_2exp_si(cl, cl, -e);

    /* the conditions of the bound, (2N + 1)(2N + 2) >= 32 c + 32 (the
       ratios of the |P_n| at most 1/2) and N >= 2 (those of the T_n,
       4 u0 / (n + 1) <= 1/2 for u0 <= 1/4), checked exactly (_init_terms
       chooses N in doubles) */
    mag_mul_ui(t, c, 32);
    mag_add_ui(t, t, 32);
    mag_set_ui_lower(s, 2 * (ulong) N + 1);
    mag_mul_ui_lower(s, s, 2 * (ulong) N + 2);
    if (N < 2 || mag_cmp(s, t) < 0)
    {
        mag_inf(res);
        goto cleanup;
    }

    mag_mul_ui(P, c, 16);
    mag_pow_ui(P, P, N);
    mag_rfac_ui(t, 2 * N);
    mag_mul(P, P, t);
    _mp_real_get_mag(t, Gb);
    mag_mul(s, P, t);
    mag_mul_2exp_si(s, s, 2);                       /* 4 Gb P_N */

    mag_set_ui_2exp_si(T, 1, 2 * N - 2 * a * N - a);    /* 4^N u0^(N + 1/2) */
    mag_rfac_ui(t, N);
    mag_mul(T, T, t);
    mag_mul_2exp_si(t, T, 2 * a + 1);
    mag_div_ui(t, t, 2 * N - 1);                    /* T_N / (u0 (N - 1/2)) */
    mag_div(P, T, cl);
    mag_add(t, t, P);                               /* + T_N / c */
    _mp_real_get_mag(P, E);
    mag_mul(t, t, P);
    mag_mul_2exp_si(t, t, 1);
    mag_add(res, s, t);

cleanup:
    mag_clear(c); mag_clear(cl); mag_clear(P); mag_clear(T); mag_clear(t); mag_clear(s);
}

/* the limbs of x > 0 (a small value is stored in *buf) */
static void
_fmpz_limbs(nn_srcptr * d, slong * len, ulong * buf, const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
    {
        *buf = (ulong) *x;
        *d = buf;
        *len = 1;
    }
    else
    {
        *d = COEFF_TO_PTR(*x)->_mp_d;
        *len = COEFF_TO_PTR(*x)->_mp_size;
    }
}

/* the state X at c0 = C 2^-e and E = e^(-c0/u0), to relative accuracy
   about 2^-wp */
static void
_initial_values(mp_real_ptr X, mp_real_t E, const fmpz_t C, slong e, slong wp,
    const _mp_real_khinchin_afe_struct * ctx)
{
    _mp_real_linrec_struct rec;
    _init_struct ic;
    fmpz_t D;
    fmpq_t y;
    mp_real_struct V0[IV], s[IW];
    mp_real_t Gb, c, t;
    mag_t tb;
    double cd, u0 = ldexp(1.0, -2 * ctx->a);
    slong N, W, i, n;

    fmpz_init(D);
    fmpz_one(D);
    fmpz_mul_2exp(D, D, e);
    fmpq_init(y);
    mp_real_init(Gb);
    mp_real_init(c);
    mp_real_init(t);
    mag_init(tb);
    for (i = 0; i < IV; i++)
        mp_real_init(V0 + i);
    for (i = 0; i < IW; i++)
        mp_real_init(s + i);

    cd = ldexp(fmpz_get_d(C), -e);
    /* the guard bits for the cancellation in the weighted sums */
    W = wp + (slong) (6.0 * sqrt(cd)) + 64;
    n = mp_real_prec_bits(W);
    N = _init_terms(cd, u0, W);

    ic.C = C;
    ic.D = D;
    ic.u0d = ctx->u0d;
    rec.d = IV;
    rec.s = IW;
    rec.g = 0;
    rec.leaf = _init_leaf;
    rec.bottom = NULL;
    rec.bottom_len = 0;
    rec.ctx = &ic;
    /* V_0 = (P_0, R_0, T_0) = (1, 0, sqrt(u0)) */
    mp_real_set_ui(V0 + 0, 1);
    mp_real_set_ui(V0 + 2, 1);
    mp_real_mul_2exp_si(V0 + 2, V0 + 2, -ctx->a);
    for (i = 0; i < IW; i++)
        mp_real_zero(s + i);
    _mp_real_linrec_eval(V0, s, NULL, &rec, 0, N, NULL, 1, mp_real_prec_bits(W + 64));

    /* c0, y = c0/u0 (exact), E = e^-y, Gb = c0^(-1/2) Gamma(1/2, y) */
    mp_real_set_fmpz(c, C);
    mp_real_mul_2exp_si(c, c, -e);
    fmpz_mul_ui(fmpq_numref(y), C, ctx->u0d);
    fmpz_set(fmpq_denref(y), D);
    fmpq_canonicalise(y);
    mp_real_set_fmpz(t, fmpq_numref(y));
    mp_real_mul_2exp_si(t, t, -(slong) fmpz_val2(fmpq_denref(y)));
    mp_real_neg(t, t);
    mp_real_exp_bits(E, t, FLINT_BITS * n);
    {
        nn_srcptr pd, qd;
        slong pn, qn;
        ulong pb, qb;
        _fmpz_limbs(&pd, &pn, &pb, fmpq_numref(y));
        _fmpz_limbs(&qd, &qn, &qb, fmpq_denref(y));
        mp_real_gamma_upper_half_frac(Gb, pd, pn, qd, qn, n);
    }
    mp_real_rsqrt(t, c, n);
    mp_real_mul(Gb, Gb, t, n);

    /* nu = Gb s0 + E s1, Y = Gb s2 + E s3, Y_(-1)/c = (-c Gb s4 + E (s5 - c s6))/c,
       C = Gb */
    mp_real_mul(X + 0, Gb, s + 0, n);
    mp_real_mul(t, E, s + 1, n);
    mp_real_add(X + 0, X + 0, t, n);
    mp_real_mul(X + 1, Gb, s + 2, n);
    mp_real_mul(t, E, s + 3, n);
    mp_real_add(X + 1, X + 1, t, n);
    mp_real_mul(X + 2, Gb, s + 4, n);
    mp_real_neg(X + 2, X + 2);
    mp_real_div(t, s + 5, c, n);
    mp_real_sub(t, t, s + 6, n);
    mp_real_mul(t, t, E, n);
    mp_real_add(X + 2, X + 2, t, n);
    mp_real_set(X + 3, Gb);

    _init_tail(tb, C, e, Gb, E, N, ctx->a);
    for (i = 0; i < NX; i++)
        _mp_real_add_error_mag(X + i, tb);

    for (i = 0; i < IV; i++)
        mp_real_clear(V0 + i);
    for (i = 0; i < IW; i++)
        mp_real_clear(s + i);
    fmpz_clear(D);
    fmpq_clear(y);
    mp_real_clear(Gb);
    mp_real_clear(c);
    mp_real_clear(t);
    mag_clear(tb);
}

/* the driver ****************************************************************/

/* bits of the starting point and the growth factor of the bit-burst */
#define DUAL_B0 64
#define DUAL_GROW 4

void
_mp_real_khinchin_afe_dual(mp_real_t res, const mp_real_t cstar,
    const _mp_real_khinchin_afe_struct * ctx)
{
    mp_real_struct X[NX];
    mp_real_t E, s;
    fmpz_t C, C1;
    mag_t mm, err;
    slong prec = ctx->prec, abs_bits = prec + 10, wp, e, B, Bfin, i;
    double cd;

    for (i = 0; i < NX; i++)
        mp_real_init(X + i);
    mp_real_init(E);
    mp_real_init(s);
    fmpz_init(C);
    fmpz_init(C1);
    mag_init(mm);
    mag_init(err);

    {
        fmpz_t z;
        fmpz_init(z);
        _mp_real_get_fmpz_fixed(z, err, cstar, 0);
        cd = fmpz_get_d(z);
        fmpz_clear(z);
    }

    /* J_m is about e^(-c/u0): the relative precision needed */
    wp = prec - (slong) (cd * ctx->u0d * 1.4426950408889634) + 64;
    wp = FLINT_MAX(wp, 64);

    /* the final bits of c: |J'| carries the factor e^(-c/u0) (below) */
    Bfin = FLINT_MIN(wp + 128, prec + 64);
    Bfin = FLINT_MAX(Bfin, 2 * DUAL_B0);

    e = DUAL_B0;
    mag_zero(err);
    _mp_real_get_fmpz_fixed(C, err, cstar, -e);     /* c0 <= c */
    _initial_values(X, E, C, e, wp + 64, ctx);

    for (;;)
    {
        int last;

        B = e * DUAL_GROW;
        last = (B >= Bfin);
        if (last)
            B = Bfin;
        mag_zero(err);
        _mp_real_get_fmpz_fixed(C1, err, cstar, -B);
        _taylor_step(X, E, C, e, C1, B, wp + 32, abs_bits, last, ctx);
        fmpz_swap(C, C1);
        e = B;
        if (last)
            break;
    }

    mp_real_const_pi4(s, mp_real_prec_bits(wp), 1);
    mp_real_mul_2exp_si(s, s, 2);
    mp_real_sqrt(s, s, mp_real_prec_bits(wp));
    mp_real_mul(res, X + 0, s, mp_real_prec_bits(wp));

    /* J was evaluated at c_l with |c - c_l| <= delta = 2^-e + rad(c), and
       for xi >= c - delta, |J'(xi)| = sqrt(pi) int_0^u0 u^(-5/2) e^(-xi/u) E(4u) du
       (E(4u) <= E(4 u0) = g1) is at most
       (a) sqrt(pi) g1 Gamma(3/2) xi^(-3/2) = (pi/2) g1 xi^(-3/2), and
       (b) 2 sqrt(pi) g1 u0^(-1/2) e^(-xi/u0) / xi
           (Gamma(a, x) <= x^(a-1) e^-x / (1 - (a-1)/x), xi/u0 >= 1);
       the smaller one is used.  Taking xi only delta below c (rather than
       a fixed distance) keeps (b) within a small factor of |J|, which
       matters since the last step has only about wp bits relative to
       |J|: a fixed distance 1 lost about 1.44/u0 bits (370 at u0 = 1/256). */
    {
        mag_t xi, d, d2, g1, t;

        mag_init(xi); mag_init(d); mag_init(d2); mag_init(g1); mag_init(t);
        mag_set_ui_2exp_si(mm, 1, -e);
        _mp_real_get_rad_mag(t, cstar);
        mag_add(mm, mm, t);
        _mp_real_get_mag_lower(xi, cstar);
        mag_sub_lower(xi, xi, mm);
        _mp_real_get_mag(g1, ctx->g1);
        if (mag_is_zero(xi))
            mag_inf(d);
        else
        {
            /* (a) */
            mag_rsqrt(d, xi);
            mag_pow_ui(d, d, 3);
            mag_mul(d, d, g1);
            mag_const_pi(t);
            mag_mul(d, d, t);
            mag_mul_2exp_si(d, d, -1);
            /* (b) */
            mag_mul_2exp_si(t, xi, 2 * ctx->a);
            mag_expinv(d2, t);
            mag_div(d2, d2, xi);
            mag_mul_2exp_si(d2, d2, ctx->a + 1);
            mag_mul(d2, d2, g1);
            mag_const_pi(t);
            mag_sqrt(t, t);
            mag_mul(d2, d2, t);
            mag_min(d, d, d2);
        }
        mag_mul(mm, mm, d);
        _mp_real_add_error_mag(res, mm);
        mag_clear(xi); mag_clear(d); mag_clear(d2); mag_clear(g1); mag_clear(t);
    }

    for (i = 0; i < NX; i++)
        mp_real_clear(X + i);
    mp_real_clear(E);
    mp_real_clear(s);
    fmpz_clear(C);
    fmpz_clear(C1);
    mag_clear(mm);
    mag_clear(err);
}
