/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_poly.h"
#include "arb.h"
#include "acb.h"
#include "acf.h"
#include "impl.h"

/*
    Near-boundary values of the elementary functions.

    A part V of f(x) that the compositions leave undecided lies within a
    tiny distance of a rounding boundary beta.  For a transcendental f
    this happens (beyond the precisions that the retries reach) only near
    an exact point z0 of f, where an expansion with exactly known
    coefficients has a dyadic leading term: with h = x - z0 exact and
    tiny,

        f(x) = C + u s sum_{j = j0}^{N} q_j (pi^k H)^j + O(|H|^(N + 1)),

    with H = v h for a unit v (or H = 1 / x at infinity), q_j rational
    (fmpq_poly series), j0 = -1 for a simple pole, a constant C in
    Q(i) pi^c, a unit u and a dyadic factor s (sqrt(w h) at a square-root
    branch point, when it is dyadic).  The sign of V - beta is then that
    of a sum of exact monomials in the parts of H (and s) with integer
    coefficients and powers of pi, after multiplying by a positive common
    denominator (pi^k |H|^2 at a pole, |x|^(2N) at infinity): exact
    within each power of pi (_acf_pt_sign), with pi as a ball across
    them, decided when the sum excludes zero by more than the bound on
    the omitted terms.
*/

/* term lists: values times pi^k */
typedef struct
{
    arf_struct * t;
    slong * k;
    slong n, alloc;
}
_acf_pt_struct;

typedef _acf_pt_struct _acf_pt_t[1];

static void
_acf_pt_init(_acf_pt_t T)
{
    T->t = NULL;
    T->k = NULL;
    T->n = T->alloc = 0;
}

static void
_acf_pt_clear(_acf_pt_t T)
{
    slong i;
    for (i = 0; i < T->alloc; i++)
        arf_clear(T->t + i);
    flint_free(T->t);
    flint_free(T->k);
}

/* T += c x y pi^k (x, y may be NULL) */
static void
_acf_pt_push(_acf_pt_t T, const arf_t c, slong k, arf_srcptr x, arf_srcptr y)
{
    arf_ptr t;

    if (T->n == T->alloc)
    {
        slong i, alloc = FLINT_MAX(16, 2 * T->alloc);
        T->t = flint_realloc(T->t, alloc * sizeof(arf_struct));
        T->k = flint_realloc(T->k, alloc * sizeof(slong));
        for (i = T->alloc; i < alloc; i++)
            arf_init(T->t + i);
        T->alloc = alloc;
    }

    t = T->t + T->n;
    arf_set(t, c);
    if (x != NULL) arf_mul(t, t, x, ARF_PREC_EXACT, ARF_RND_DOWN);
    if (y != NULL) arf_mul(t, t, y, ARF_PREC_EXACT, ARF_RND_DOWN);
    if (!arf_is_zero(t))
    {
        T->k[T->n] = k;
        T->n++;
    }
}

/* R += c X Y (sign: +-1) */
static void
_acf_pt_addmul(_acf_pt_t R, const _acf_pt_t X, const _acf_pt_t Y, int sign)
{
    slong i, j;
    arf_t c;

    arf_init(c);
    for (i = 0; i < X->n; i++)
        for (j = 0; j < Y->n; j++)
        {
            arf_set_si(c, sign);
            _acf_pt_push(R, c, X->k[i] + Y->k[j], X->t + i, Y->t + j);
        }
    arf_clear(c);
}

/* R += c X pi^k */
static void
_acf_pt_addscal(_acf_pt_t R, const _acf_pt_t X, const arf_t c, slong k)
{
    slong i;
    for (i = 0; i < X->n; i++)
        _acf_pt_push(R, c, X->k[i] + k, X->t + i, NULL);
}

/* complex term lists */
typedef struct
{
    _acf_pt_t re, im;
}
_acf_cpt_struct;

typedef _acf_cpt_struct _acf_cpt_t[1];

static void _acf_cpt_init(_acf_cpt_t X) { _acf_pt_init(X->re); _acf_pt_init(X->im); }
static void _acf_cpt_clear(_acf_cpt_t X) { _acf_pt_clear(X->re); _acf_pt_clear(X->im); }
static void _acf_cpt_zero(_acf_cpt_t X) { X->re->n = 0; X->im->n = 0; }

static void
_acf_cpt_swap(_acf_cpt_t X, _acf_cpt_t Y)
{
    _acf_cpt_struct t = *X;
    *X = *Y;
    *Y = t;
}

/* R += X Y */
static void
_acf_cpt_addmul(_acf_cpt_t R, const _acf_cpt_t X, const _acf_cpt_t Y)
{
    _acf_pt_addmul(R->re, X->re, Y->re, 1);
    _acf_pt_addmul(R->re, X->im, Y->im, -1);
    _acf_pt_addmul(R->im, X->re, Y->im, 1);
    _acf_pt_addmul(R->im, X->im, Y->re, 1);
}

/* R += c i^u X pi^k */
static void
_acf_cpt_addscal(_acf_cpt_t R, const _acf_cpt_t X, const arf_t c, slong k, int u)
{
    arf_t nc;
    arf_init(nc);
    arf_neg(nc, c);
    u &= 3;
    if (u == 0) { _acf_pt_addscal(R->re, X->re, c, k);  _acf_pt_addscal(R->im, X->im, c, k); }
    if (u == 1) { _acf_pt_addscal(R->re, X->im, nc, k); _acf_pt_addscal(R->im, X->re, c, k); }
    if (u == 2) { _acf_pt_addscal(R->re, X->re, nc, k); _acf_pt_addscal(R->im, X->im, nc, k); }
    if (u == 3) { _acf_pt_addscal(R->re, X->im, c, k);  _acf_pt_addscal(R->im, X->re, nc, k); }
    arf_clear(nc);
}

/* X = a + bi */
static void
_acf_cpt_set(_acf_cpt_t X, const arf_t a, const arf_t b)
{
    arf_t one;
    arf_init(one);
    arf_one(one);
    _acf_cpt_zero(X);
    _acf_pt_push(X->re, one, 0, a, NULL);
    _acf_pt_push(X->im, one, 0, b, NULL);
    arf_clear(one);
}

static int
_acf_cmp_pt_desc(const void * x, const void * y)
{
    arf_srcptr u = *(arf_srcptr const *) x, v = *(arf_srcptr const *) y;
    return -fmpz_cmp(ARF_EXPREF(u), ARF_EXPREF(v));
}

/* the sign of the sum of the exact terms t[0], ..., t[n - 1] (as in
   _acf_sign_sum), and in S that sum, exactly or within less than half
   of itself */
static int
_acf_sign_sum_ptr(arf_t S, arf_srcptr * t, slong n)
{
    fmpz_t e;
    slong i;

    arf_zero(S);
    if (n == 0)
        return 0;

    qsort(t, n, sizeof(arf_srcptr), _acf_cmp_pt_desc);
    fmpz_init(e);
    for (i = 0; i < n; i++)
    {
        if (!arf_is_zero(S))
        {
            fmpz_add_ui(e, ARF_EXPREF(t[i]), FLINT_BIT_COUNT(n - i) + 2);
            if (fmpz_cmp(e, ARF_EXPREF(S)) <= 0)
                break;
        }
        arf_add(S, S, t[i], ARF_PREC_EXACT, ARF_RND_DOWN);
    }
    fmpz_clear(e);
    return arf_sgn(S);
}

/* The sum of the n exact terms t as a ball at precision w: the correctly
   rounded sum (arf_sum, for any spread of the exponents) within an ulp,
   so that the ball has the sign of the sum, or is an exact zero. */
static void
_acf_sum_ball(arb_t res, arf_srcptr t, slong n, slong w)
{
    if (arf_sum(arb_midref(res), t, n, w, ARF_RND_DOWN))
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), w);
    else
        mag_zero(arb_radref(res));
}

/* The sign of the sum of the terms of P (with their powers of pi), when
   its absolute value exceeds rem: 1 and *side, or 0 (undecided).  Each
   power of pi gathers a sum G_k (all k >= 0), exact for a single power,
   and otherwise as a ball at precision w with the sign of G_k, entering
   the sum of the pi^k G_k in ball arithmetic. */
static int
_acf_pt_sign(int * side, const _acf_pt_t P, const mag_t rem)
{
    slong i, j, kmin, kmax, w;
    arf_srcptr * t;
    arf_t S;
    int ok = 0;

    if (P->n == 0)
        return 0;

    kmin = kmax = P->k[0];
    for (i = 1; i < P->n; i++)
    {
        kmin = FLINT_MIN(kmin, P->k[i]);
        kmax = FLINT_MAX(kmax, P->k[i]);
    }

    t = flint_malloc(P->n * sizeof(arf_srcptr));
    arf_init(S);

    if (kmin == kmax)
    {
        /* one power of pi: exact */
        int s;
        mag_t m;
        for (i = 0; i < P->n; i++)
            t[i] = P->t + i;
        s = _acf_sign_sum_ptr(S, t, P->n);
        mag_init(m);
        arf_get_mag_lower(m, S);
        mag_mul_2exp_si(m, m, -1);
        /* (pi^kmin >= 1) */
        FLINT_ASSERT(kmin >= 0);
        if (s != 0 && mag_cmp(m, rem) > 0)
        {
            *side = s;
            ok = 1;
        }
        mag_clear(m);
    }
    else
    {
        arb_t tot, g, pik, pi;
        arf_struct * tj = flint_malloc(P->n * sizeof(arf_struct));
        arb_init(tot);
        arb_init(g);
        arb_init(pik);
        arb_init(pi);

        for (w = 128; w <= (WORD(1) << 16) && !ok; w *= 2)
        {
            arb_zero(tot);
            arb_const_pi(pi, w);
            for (j = kmin; j <= kmax; j++)
            {
                slong n = 0;

                /* (shallow copies of the terms with pi^j, read only) */
                for (i = 0; i < P->n; i++)
                    if (P->k[i] == j)
                        tj[n++] = P->t[i];
                if (n == 0)
                    continue;
                _acf_sum_ball(g, tj, n, w);
                if (arb_is_zero(g))
                    continue;
                arb_pow_ui(pik, pi, j, w);
                arb_addmul(tot, g, pik, w);
            }
            if (arb_is_nonzero(tot))
            {
                mag_t m;
                mag_init(m);
                arb_get_mag_lower(m, tot);
                if (mag_cmp(m, rem) > 0)
                {
                    *side = arb_sgn_nonzero(tot);
                    ok = 1;
                }
                mag_clear(m);
                if (!ok)
                    break;      /* (the remainder too large: undecided) */
            }
        }

        arb_clear(tot);
        arb_clear(g);
        arb_clear(pik);
        arb_clear(pi);
        flint_free(tj);
    }

    arf_clear(S);
    flint_free(t);
    return ok;
}

/* A local expansion of f at an exact point (see above). */
typedef struct
{
    const fmpq_poly_struct * Q;     /* q_{j0}, q_{j0 + 1}, ... */
    slong j0;                       /* 0, or -1 for a simple pole */
    int kpi;                        /* the variable pi^kpi H */
    int unit;                       /* the series times i^unit */
    const acf_struct * s;           /* a dyadic factor, or NULL */
    fmpq_t cre, cim;                /* the constant (cre + cim i) pi^cpi */
    int cpi;
    double M, rho;                  /* |q_j| <= M rho^j */
}
_acf_local_struct;

/* sgn(Part(f(x)) - pi^bpi beta) from the expansion L at H (exact; inf:
   H = 1 / x, x given), or 0 when undecided */
static int
_acf_local_side(int * side, const _acf_local_struct * L, const acf_t H, int inf, int part, const arf_t beta, int bpi)
{
    _acf_cpt_t P, Hj, T, Hc, A, D;
    slong N, j, len;
    fmpz_t den;
    arf_t c, one;
    mag_t hm, rem, t;
    int ok = 0;

    len = fmpq_poly_length(L->Q);
    N = L->j0 + len - 1;            /* the last power */
    if (N < 0)
        return 0;

    mag_init(hm);
    mag_init(rem);
    mag_init(t);

    /* |kappa H| (inf: 1 / |x|) */
    {
        mag_t u;
        mag_init(u);
        arf_get_mag(hm, acf_realref(H));
        arf_get_mag(u, acf_imagref(H));
        mag_hypot(hm, hm, u);
        if (inf)
        {
            arf_get_mag_lower(hm, acf_realref(H));
            arf_get_mag_lower(u, acf_imagref(H));
            mag_max(hm, hm, u);
            mag_inv(hm, hm);
        }
        if (L->kpi)
            mag_mul_ui(hm, hm, 4);
        mag_mul_ui(hm, hm, (ulong) (L->rho * 16 + 1));
        mag_mul_2exp_si(hm, hm, -4);
        mag_clear(u);
    }
    /* rho kappa |H| <= 2^-8 */
    if (mag_cmp_2exp_si(hm, -8) > 0)
    {
        mag_clear(hm);
        mag_clear(rem);
        mag_clear(t);
        return 0;
    }

    _acf_cpt_init(P);
    _acf_cpt_init(Hj);
    _acf_cpt_init(T);
    _acf_cpt_init(Hc);
    _acf_cpt_init(A);
    _acf_cpt_init(D);
    fmpz_init(den);
    arf_init(c);
    arf_init(one);
    arf_one(one);

    /* the common denominator of the coefficients and the constant */
    fmpz_lcm(den, L->Q->den, fmpq_denref(L->cre));
    fmpz_lcm(den, den, fmpq_denref(L->cim));

    /* the positive multiplier D: 1, pi^kpi |H|^2 (pole) or |x|^(2N) (inf) */
    _acf_cpt_zero(D);
    if (L->j0 < 0)
    {
        _acf_pt_push(D->re, one, L->kpi, acf_realref(H), acf_realref(H));
        _acf_pt_push(D->re, one, L->kpi, acf_imagref(H), acf_imagref(H));
    }
    else if (inf)
    {
        /* (a^2 + b^2)^N */
        _acf_cpt_t X2;
        _acf_cpt_init(X2);
        _acf_cpt_zero(X2);
        _acf_pt_push(X2->re, one, 0, acf_realref(H), acf_realref(H));
        _acf_pt_push(X2->re, one, 0, acf_imagref(H), acf_imagref(H));
        _acf_pt_push(D->re, one, 0, NULL, NULL);
        for (j = 0; j < N; j++)
        {
            _acf_cpt_zero(T);
            _acf_pt_addmul(T->re, D->re, X2->re, 1);
            _acf_cpt_swap(T, D);
        }
        _acf_cpt_clear(X2);
    }
    else
        _acf_pt_push(D->re, one, 0, NULL, NULL);

    /* Hc: H, or conj(x) at infinity */
    _acf_cpt_set(Hc, acf_realref(H), acf_imagref(H));
    if (inf)
    {
        arf_t nb;
        arf_init(nb);
        arf_neg(nb, acf_imagref(H));
        _acf_cpt_set(Hc, acf_realref(H), nb);
        arf_clear(nb);
    }

    /* A = sum_j den q_j pi^(kpi j) H^j, multiplied by D */
    _acf_cpt_zero(A);
    _acf_cpt_zero(Hj);
    _acf_pt_push(Hj->re, one, 0, NULL, NULL);       /* H^0 */
    for (j = 0; j <= N; j++)
    {
        slong i = j - L->j0;      /* the index of q_j in Q */

        if (j > 0)
        {
            _acf_cpt_zero(T);
            _acf_cpt_addmul(T, Hj, Hc);
            _acf_cpt_swap(T, Hj);
        }
        if (i >= 0 && i < len && !fmpz_is_zero(L->Q->coeffs + i))
        {
            /* c = den q_j (an integer) */
            fmpz_t q;
            fmpz_init(q);
            fmpz_divexact(q, den, L->Q->den);
            fmpz_mul(q, q, L->Q->coeffs + i);
            arf_set_fmpz(c, q);
            fmpz_clear(q);

            _acf_cpt_zero(T);
            if (L->j0 < 0)
            {
                /* H^j pi^(kpi j) times pi^kpi |H|^2 */
                _acf_cpt_addmul(T, Hj, D);
                _acf_cpt_addscal(A, T, c, L->kpi * j, 0);
            }
            else if (inf)
            {
                /* conj(x)^j |x|^(2 (N - j)) */
                _acf_cpt_t E, X2;
                slong m;
                _acf_cpt_init(E);
                _acf_cpt_init(X2);
                _acf_pt_push(E->re, one, 0, NULL, NULL);
                _acf_pt_push(X2->re, one, 0, acf_realref(H), acf_realref(H));
                _acf_pt_push(X2->re, one, 0, acf_imagref(H), acf_imagref(H));
                for (m = 0; m < N - j; m++)
                {
                    _acf_cpt_zero(T);
                    _acf_pt_addmul(T->re, E->re, X2->re, 1);
                    _acf_cpt_swap(T, E);
                }
                _acf_cpt_zero(T);
                _acf_cpt_addmul(T, Hj, E);
                _acf_cpt_addscal(A, T, c, L->kpi * j, 0);
                _acf_cpt_clear(E);
                _acf_cpt_clear(X2);
            }
            else
                _acf_cpt_addscal(A, Hj, c, L->kpi * j, 0);
        }
    }
    if (L->j0 < 0 && !fmpz_is_zero(L->Q->coeffs))
    {
        /* q_{-1} / (pi^kpi H) times pi^kpi |H|^2 = q_{-1} conj(H) */
        fmpz_t q;
        arf_t nb;
        fmpz_init(q);
        arf_init(nb);
        fmpz_divexact(q, den, L->Q->den);
        fmpz_mul(q, q, L->Q->coeffs);
        arf_set_fmpz(c, q);
        _acf_pt_push(A->re, c, 0, acf_realref(H), NULL);
        arf_neg(nb, acf_imagref(H));
        _acf_pt_push(A->im, c, 0, nb, NULL);
        fmpz_clear(q);
        arf_clear(nb);
    }

    /* P = i^unit s A, plus the constant times den D, minus beta */
    _acf_cpt_zero(P);
    if (L->s != NULL)
    {
        _acf_cpt_set(T, acf_realref(L->s), acf_imagref(L->s));
        _acf_cpt_zero(Hj);
        _acf_cpt_addmul(Hj, T, A);
        _acf_cpt_addscal(P, Hj, one, 0, L->unit);
    }
    else
        _acf_cpt_addscal(P, A, one, 0, L->unit);

    {
        fmpz_t q;
        fmpz_init(q);
        fmpz_divexact(q, den, fmpq_denref(L->cre));
        fmpz_mul(q, q, fmpq_numref(L->cre));
        arf_set_fmpz(c, q);
        _acf_pt_addscal(P->re, D->re, c, L->cpi);
        fmpz_divexact(q, den, fmpq_denref(L->cim));
        fmpz_mul(q, q, fmpq_numref(L->cim));
        arf_set_fmpz(c, q);
        _acf_pt_addscal(P->im, D->re, c, L->cpi);
        fmpz_clear(q);
    }

    arf_set_fmpz(c, den);
    arf_mul(c, c, beta, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_neg(c, c);
    _acf_pt_addscal(part ? P->im : P->re, D->re, c, bpi);

    /* the omitted terms: den |s| D 2 M (rho kappa |H|)^(N + 1), with D
       bounded by (|H|^2 pi^kpi) resp. |x|^(2N) */
    mag_pow_ui(rem, hm, N + 1);
    mag_mul_ui(rem, rem, (ulong) (2 * L->M + 1));
    mag_set_fmpz(t, den);
    mag_mul(rem, rem, t);
    if (L->s != NULL)
    {
        mag_t u;
        mag_init(u);
        arf_get_mag(t, acf_realref(L->s));
        arf_get_mag(u, acf_imagref(L->s));
        mag_hypot(t, t, u);
        mag_mul(rem, rem, t);
        mag_clear(u);
    }
    if (L->j0 < 0 || inf)
    {
        mag_t u;
        mag_init(u);
        arf_get_mag(t, acf_realref(H));
        arf_get_mag(u, acf_imagref(H));
        mag_hypot(t, t, u);
        if (L->j0 < 0)
        {
            mag_mul(t, t, t);
            if (L->kpi)
                mag_mul_ui(t, t, 4);
        }
        else
            mag_pow_ui(t, t, 2 * N);
        mag_mul(rem, rem, t);
        mag_clear(u);
    }
    if (bpi > 0 || L->cpi > 0 || L->kpi > 0)
    {
        /* (the powers of pi among the omitted terms) */
        mag_mul_2exp_si(rem, rem, 2 * (N + 2));
    }

    ok = _acf_pt_sign(side, part ? P->im : P->re, rem);

    _acf_cpt_clear(P);
    _acf_cpt_clear(Hj);
    _acf_cpt_clear(T);
    _acf_cpt_clear(Hc);
    _acf_cpt_clear(A);
    _acf_cpt_clear(D);
    fmpz_clear(den);
    arf_clear(c);
    arf_clear(one);
    mag_clear(hm);
    mag_clear(rem);
    mag_clear(t);
    return ok;
}

/* the series of the expansions, q_0 .. q_{n-1} in t */
enum
{
    S_EXP, S_EXPM1, S_SIN, S_COS, S_TAN, S_SEC, S_TCOT, S_TCSC, S_LOG1P,
    S_ATAN, S_ASIN, S_SINC, S_TANQ, S_COTQ, S_ACOSB
};

static void
_acf_series(fmpq_poly_t Q, int kind, slong n)
{
    fmpq_poly_t t, u;

    fmpq_poly_init(t);
    fmpq_poly_init(u);
    fmpq_poly_set_coeff_si(t, 1, 1);        /* t */

    switch (kind)
    {
        case S_EXP: fmpq_poly_exp_series(Q, t, n); break;
        case S_EXPM1: fmpq_poly_exp_series(Q, t, n); fmpq_poly_set_coeff_si(Q, 0, 0); break;
        case S_SIN: fmpq_poly_sin_series(Q, t, n); break;
        case S_COS: fmpq_poly_cos_series(Q, t, n); break;
        case S_TAN: fmpq_poly_tan_series(Q, t, n); break;
        case S_SEC: fmpq_poly_cos_series(u, t, n); fmpq_poly_inv_series(Q, u, n); break;
        case S_TCOT:
            /* t cos t / sin t = cos t / (sin t / t) */
            fmpq_poly_sin_series(u, t, n + 1);
            fmpq_poly_shift_right(u, u, 1);
            fmpq_poly_cos_series(Q, t, n);
            fmpq_poly_div_series(Q, Q, u, n);
            break;
        case S_TCSC:
            fmpq_poly_sin_series(u, t, n + 1);
            fmpq_poly_shift_right(u, u, 1);
            fmpq_poly_inv_series(Q, u, n);
            break;
        case S_LOG1P: fmpq_poly_set_coeff_si(t, 0, 1); fmpq_poly_log_series(Q, t, n); break;
        case S_ATAN: fmpq_poly_atan_series(Q, t, n); break;
        case S_ASIN: fmpq_poly_asin_series(Q, t, n); break;
        case S_SINC:
            fmpq_poly_sin_series(Q, t, n + 1);
            fmpq_poly_shift_right(Q, Q, 1);
            break;
        case S_TANQ:
        case S_COTQ:
            /* (1 + T) / (1 - T) resp. its inverse, T = tan t */
            fmpq_poly_tan_series(u, t, n);
            fmpq_poly_one(Q);
            fmpq_poly_add(Q, Q, u);
            fmpq_poly_one(t);
            fmpq_poly_sub(t, t, u);
            if (kind == S_TANQ)
                fmpq_poly_div_series(Q, Q, t, n);
            else
                fmpq_poly_div_series(Q, t, Q, n);
            break;
        case S_ACOSB:
            {
                /* asin(sqrt(v)) / sqrt(v): the odd series halved in
                   degree */
                slong i;
                fmpq_poly_asin_series(u, t, 2 * n);
                fmpq_poly_zero(Q);
                for (i = 0; i < n; i++)
                {
                    fmpq_t c;
                    fmpq_init(c);
                    fmpq_poly_get_coeff_fmpq(c, u, 2 * i + 1);
                    fmpq_poly_set_coeff_fmpq(Q, i, c);
                    fmpq_clear(c);
                }
            }
            break;
        default:
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    }

    fmpq_poly_clear(t);
    fmpq_poly_clear(u);
}

#define ACF_NEAR_TERMS 10

/* the nearest integer to 2^e Re(x), when it fits */
static int
_acf_near_index(slong * n, const arf_t a, slong e)
{
    arf_t t;
    fmpz_t z;
    int ok;

    arf_init(t);
    fmpz_init(z);
    arf_mul_2exp_si(t, a, e);
    arf_get_fmpz(z, t, ARF_RND_NEAR);
    ok = fmpz_bits(z) < FLINT_BITS - 4;
    if (ok)
        *n = fmpz_get_si(z);
    arf_clear(t);
    fmpz_clear(z);
    return ok;
}

/* h = x - z0 (z0 = n 2^-e + m i) rotated by i^v: into H (exactly) */
static void
_acf_near_offset(acf_t H, const acf_t x, slong n, slong e, slong m, int v)
{
    arf_t t;
    arf_init(t);
    arf_set_si_2exp_si(t, n, -e);
    arf_sub(acf_realref(H), acf_realref(x), t, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_set_si(t, m);
    arf_sub(acf_imagref(H), acf_imagref(x), t, ARF_PREC_EXACT, ARF_RND_DOWN);
    v &= 3;
    while (v-- > 0)
    {
        /* i (p + qi) = -q + pi */
        arf_swap(acf_realref(H), acf_imagref(H));
        arf_neg(acf_realref(H), acf_realref(H));
    }
    arf_clear(t);
}

/* Q(c t): the coefficients scaled by c^j */
static void
_acf_series_scale(fmpq_poly_t Q, slong cnum, ulong cden)
{
    fmpq_t p, q, c;
    slong i;
    fmpq_init(p);
    fmpq_init(q);
    fmpq_init(c);
    fmpq_set_si(c, cnum, cden);
    fmpq_one(p);
    for (i = 0; i < fmpq_poly_length(Q); i++)
    {
        fmpq_poly_get_coeff_fmpq(q, Q, i);
        fmpq_mul(q, q, p);
        fmpq_poly_set_coeff_fmpq(Q, i, q);
        fmpq_mul(p, p, c);
    }
    fmpq_clear(p);
    fmpq_clear(q);
    fmpq_clear(c);
}

/* the squared distance |x - z0|^2 as a mag, z0 = n + m i */
static void
_acf_near_dist(mag_t d, const acf_t x, slong n, slong m)
{
    arf_t t;
    mag_t u;
    arf_init(t);
    mag_init(u);
    arf_sub_si(t, acf_realref(x), n, MAG_BITS, ARF_RND_UP);
    arf_get_mag(d, t);
    arf_sub_si(t, acf_imagref(x), m, MAG_BITS, ARF_RND_UP);
    arf_get_mag(u, t);
    mag_hypot(d, d, u);
    arf_clear(t);
    mag_clear(u);
}

/* the square-root factor s = sqrt(c h) at a branch point, when dyadic */
static int
_acf_near_branch(acf_t s, const acf_t h, slong c)
{
    acf_t t;
    int ok;
    acf_init(t);
    arf_mul_si(acf_realref(t), acf_realref(h), c, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul_si(acf_imagref(t), acf_imagref(h), c, ARF_PREC_EXACT, ARF_RND_DOWN);
    ok = !arf_is_zero(acf_realref(t)) && !arf_is_zero(acf_imagref(t)) && _acf_sqrt_exact(s, t);
    acf_clear(t);
    return ok;
}

/* the expansion of the base function g at the exact point nearest to x:
   L (with Q), H and *inf; 0 when there is none */
static int
_acf_near_base(_acf_local_struct * L, fmpq_poly_t Q, acf_t H, acf_t S, int * inf, int g, const acf_t x)
{
    slong n;
    int kind = -1, unit = 0, pole = 0, kpi = 0;

    *inf = 0;
    L->s = NULL;
    L->M = 1;
    L->rho = 1;

    switch (g)
    {
        case ACF_NEAR_EXP:
        case ACF_NEAR_SIN:
        case ACF_NEAR_COS:
        case ACF_NEAR_TAN:
        case ACF_NEAR_COT:
        case ACF_NEAR_SEC:
        case ACF_NEAR_CSC:
        case ACF_NEAR_EXPM1:
        case ACF_NEAR_SINC:
        case ACF_NEAR_SINC_PI:
            acf_set(H, x);
            kind = (g == ACF_NEAR_EXP) ? S_EXP : (g == ACF_NEAR_SIN) ? S_SIN : (g == ACF_NEAR_COS) ? S_COS :
                   (g == ACF_NEAR_TAN) ? S_TAN : (g == ACF_NEAR_COT) ? S_TCOT : (g == ACF_NEAR_SEC) ? S_SEC :
                   (g == ACF_NEAR_CSC) ? S_TCSC : (g == ACF_NEAR_EXPM1) ? S_EXPM1 : S_SINC;
            pole = (g == ACF_NEAR_COT || g == ACF_NEAR_CSC);
            kpi = (g == ACF_NEAR_SINC_PI);
            L->M = (g == ACF_NEAR_TAN || g == ACF_NEAR_COT || g == ACF_NEAR_SEC || g == ACF_NEAR_CSC) ? 2 : 1;
            break;

        case ACF_NEAR_EXP_PI_I:
            /* exp(i pi (n/2 + h)) = i^n exp(pi (i h)) */
            if (!_acf_near_index(&n, acf_realref(x), 1))
                break;
            _acf_near_offset(H, x, n, 1, 0, 1);
            kind = S_EXP;
            kpi = 1;
            unit = (int) (((n % 4) + 4) % 4);
            break;

        case ACF_NEAR_SIN_PI:
        case ACF_NEAR_COS_PI:
            /* sin(pi (n/2 + h)) = (-1)^(n/2) sin(pi h) (n even),
               (-1)^((n-1)/2) cos(pi h) (n odd); cos(pi (n/2 + h)) =
               (-1)^(n/2) cos(pi h) (n even), -(-1)^((n-1)/2) sin(pi h) */
            if (!_acf_near_index(&n, acf_realref(x), 1))
                break;
            _acf_near_offset(H, x, n, 1, 0, 0);
            kpi = 1;
            if (n % 2 == 0)
            {
                kind = (g == ACF_NEAR_SIN_PI) ? S_SIN : S_COS;
                unit = ((n / 2) % 2 != 0) ? 2 : 0;
            }
            else
            {
                kind = (g == ACF_NEAR_SIN_PI) ? S_COS : S_SIN;
                unit = ((((n - 1) / 2) % 2 != 0) ^ (g == ACF_NEAR_COS_PI)) ? 2 : 0;
            }
            break;

        case ACF_NEAR_TAN_PI:
        case ACF_NEAR_COT_PI:
            {
                /* at k/4: tan(pi h), (1 + T) / (1 - T), -cot(pi h),
                   -(1 - T) / (1 + T) for k = 0, 1, 2, 3 mod 4; cot
                   likewise */
                static const int tk[4] = { S_TAN, S_TANQ, S_TCOT, S_COTQ };
                static const int ck[4] = { S_TCOT, S_COTQ, S_TAN, S_TANQ };
                int m;
                if (!_acf_near_index(&n, acf_realref(x), 2))
                    break;
                _acf_near_offset(H, x, n, 2, 0, 0);
                kpi = 1;
                L->M = 2;
                L->rho = 2;
                m = (int) (((n % 4) + 4) % 4);
                kind = (g == ACF_NEAR_TAN_PI) ? tk[m] : ck[m];
                unit = (m >= 2) ? 2 : 0;
                pole = (kind == S_TCOT);
            }
            break;

        case ACF_NEAR_SEC_PI:
        case ACF_NEAR_CSC_PI:
            /* sec(pi (n/2 + h)) = (-1)^(n/2) sec(pi h), -(-1)^((n-1)/2)
               csc(pi h); csc(pi (n/2 + h)) = (-1)^(n/2) csc(pi h),
               (-1)^((n-1)/2) sec(pi h) */
            if (!_acf_near_index(&n, acf_realref(x), 1))
                break;
            _acf_near_offset(H, x, n, 1, 0, 0);
            kpi = 1;
            L->M = 2;
            kind = ((n % 2 == 0) ^ (g == ACF_NEAR_CSC_PI)) ? S_SEC : S_TCSC;
            pole = (kind == S_TCSC);
            if (n % 2 == 0)
                unit = ((n / 2) % 2 != 0) ? 2 : 0;
            else
                unit = ((((n - 1) / 2) % 2 != 0) ^ (g == ACF_NEAR_SEC_PI)) ? 2 : 0;
            break;

        case ACF_NEAR_LOG:
            {
                /* log(z0 + h) = log(z0) + log1p(conj(z0) h) at the nearest
                   z0 = 1, i, -1, -i, log(z0) = theta pi i (theta = 0,
                   1/2, +-1 by the side of the cut, -1/2) */
                arf_srcptr a = acf_realref(x), b = acf_imagref(x);
                int q;
                if (arf_cmpabs(a, b) >= 0)
                    q = (arf_sgn(a) > 0) ? 0 : 2;
                else
                    q = (arf_sgn(b) > 0) ? 1 : 3;
                _acf_near_offset(H, x, (q == 0) ? 1 : (q == 2) ? -1 : 0, 0, (q == 1) ? 1 : (q == 3) ? -1 : 0, 4 - q);
                kind = S_LOG1P;
                if (q == 1)
                    fmpq_set_si(L->cim, 1, 2);
                else if (q == 3)
                    fmpq_set_si(L->cim, -1, 2);
                else if (q == 2)
                    fmpq_set_si(L->cim, (arf_sgn(b) < 0) ? -1 : 1, 1);
                L->cpi = 1;
            }
            break;

        case ACF_NEAR_ASIN:
        case ACF_NEAR_ACOS:
            {
                /* at 0: asin(h), acos(h) = pi/2 - asin(h); at +-1:
                   acos(1 + h) = s G(-h/2) with s = sqrt(-2h), acos(-1 + h)
                   = pi - s G(h/2) with s = sqrt(2h), G(v) = asin(sqrt(v)) /
                   sqrt(v), asin = pi/2 - acos */
                mag_t d0, d1;
                int z;
                mag_init(d0);
                mag_init(d1);
                _acf_near_dist(d0, x, 0, 0);
                _acf_near_dist(d1, x, (arf_sgn(acf_realref(x)) >= 0) ? 1 : -1, 0);
                z = (mag_cmp(d0, d1) <= 0) ? 0 : (arf_sgn(acf_realref(x)) >= 0) ? 1 : -1;
                mag_clear(d0);
                mag_clear(d1);
                _acf_near_offset(H, x, z, 0, 0, 0);
                if (z == 0)
                {
                    kind = S_ASIN;
                    if (g == ACF_NEAR_ACOS)
                    {
                        unit = 2;
                        fmpq_set_si(L->cre, 1, 2);
                        L->cpi = 1;
                    }
                }
                else
                {
                    if (!_acf_near_branch(S, H, -2 * z))
                        break;
                    L->s = S;
                    kind = S_ACOSB;
                    _acf_series(Q, kind, ACF_NEAR_TERMS);
                    _acf_series_scale(Q, -z, 2);
                    L->rho = 1;
                    /* acos: s G (z = 1), pi - s G (z = -1); asin = pi/2 - acos */
                    if (g == ACF_NEAR_ACOS)
                    {
                        unit = (z == 1) ? 0 : 2;
                        fmpq_set_si(L->cre, (z == 1) ? 0 : 1, 1);
                    }
                    else
                    {
                        unit = (z == 1) ? 2 : 0;
                        fmpq_set_si(L->cre, (z == 1) ? 1 : -1, 2);
                    }
                    L->cpi = 1;
                    L->Q = Q;
                    L->j0 = 0;
                    L->kpi = 0;
                    L->unit = unit;
                    return 1;
                }
            }
            break;

        case ACF_NEAR_ATAN:
        case ACF_NEAR_ACOT:
            {
                /* atan: at 0 atan(h), at infinity sigma pi/2 - atan(1/x);
                   acot: at infinity atan(1/x), at 0 sigma pi/2 - atan(h);
                   sigma = sgn(Re x) (the axes are handled by the callers) */
                int large = (arf_cmpabs_2exp_si(acf_realref(x), 0) > 0 || arf_cmpabs_2exp_si(acf_imagref(x), 0) > 0);
                int sigma = (arf_sgn(acf_realref(x)) >= 0) ? 1 : -1;
                acf_set(H, x);
                *inf = large;
                kind = S_ATAN;
                if (large ^ (g == ACF_NEAR_ACOT))
                {
                    unit = 2;
                    fmpq_set_si(L->cre, sigma, 2);
                    L->cpi = 1;
                }
            }
            break;

        case ACF_NEAR_ASEC:
        case ACF_NEAR_ACSC:
            /* at infinity: asec(x) = pi/2 - asin(1/x), acsc(x) = asin(1/x) */
            acf_set(H, x);
            *inf = 1;
            kind = S_ASIN;
            if (g == ACF_NEAR_ASEC)
            {
                unit = 2;
                fmpq_set_si(L->cre, 1, 2);
                L->cpi = 1;
            }
            break;

        default:
            break;
    }

    if (kind < 0)
        return 0;

    _acf_series(Q, kind, ACF_NEAR_TERMS);
    L->Q = Q;
    L->j0 = pole ? -1 : 0;
    L->kpi = kpi;
    L->unit = unit;
    return 1;
}

/* sgn(Part(f(x)) - pi^bpi beta) for the function fn of acf (the part:
   0 real, 1 imaginary) from its expansion at the exact point nearest to
   x; returns 1 when decided (in *side), 0 otherwise.  The hyperbolic
   inverse functions are rotations of the circular ones (asinh x = -i
   asin(ix), atanh x = -i atan(ix), acsch x = -i acsc(-ix), acoth x = -i
   acot(-ix), acosh x = sigma i acos(x) and asech x = -sigma i asec(x)
   with sigma = sgn(Im x)), the pi versions compare f with pi beta, and
   log1p(x) = log(1 + x). */
int
_acf_near_side(int * side, int fn, int part, const acf_t x, const arf_t beta, int bpi)
{
    _acf_local_struct L;
    fmpq_poly_t Q;
    acf_t H, S, y;
    int ok = 0, inf, g = fn, rot = 0, r;

    fmpq_poly_init(Q);
    acf_init(H);
    acf_init(S);
    acf_init(y);
    fmpq_init(L.cre);
    fmpq_init(L.cim);
    L.cpi = 0;
    acf_set(y, x);

    switch (fn)
    {
        case ACF_NEAR_ASINH: g = ACF_NEAR_ASIN; rot = 1; break;
        case ACF_NEAR_ATANH: g = ACF_NEAR_ATAN; rot = 1; break;
        case ACF_NEAR_ACSCH: g = ACF_NEAR_ACSC; rot = 3; break;
        case ACF_NEAR_ACOTH: g = ACF_NEAR_ACOT; rot = 3; break;
        case ACF_NEAR_ASIN_PI: g = ACF_NEAR_ASIN; bpi++; break;
        case ACF_NEAR_ACOS_PI: g = ACF_NEAR_ACOS; bpi++; break;
        case ACF_NEAR_ATAN_PI: g = ACF_NEAR_ATAN; bpi++; break;
        case ACF_NEAR_ACOSH: g = ACF_NEAR_ACOS; break;
        case ACF_NEAR_ASECH: g = ACF_NEAR_ASEC; break;
        case ACF_NEAR_LOG1P:
            g = ACF_NEAR_LOG;
            arf_add_ui(acf_realref(y), acf_realref(y), 1, ARF_PREC_EXACT, ARF_RND_DOWN);
            break;
        default: break;
    }

    /* y = i^rot x */
    for (r = 0; r < rot; r++)
    {
        arf_swap(acf_realref(y), acf_imagref(y));
        arf_neg(acf_realref(y), acf_realref(y));
    }

    if (_acf_near_base(&L, Q, H, S, &inf, g, y))
    {
        /* the result times i^u: -i for the rotations, sigma i for acosh
           and -sigma i for asech */
        int u = 0;
        if (rot != 0)
            u = 3;
        else if (fn == ACF_NEAR_ACOSH || fn == ACF_NEAR_ASECH)
        {
            int sigma = arf_sgn(acf_imagref(x)) * ((fn == ACF_NEAR_ASECH) ? -1 : 1);
            u = (sigma > 0) ? 1 : 3;
        }
        if (u != 0 && !(arf_sgn(acf_imagref(x)) == 0 && (fn == ACF_NEAR_ACOSH || fn == ACF_NEAR_ASECH)))
        {
            /* i^u (cre + cim i) */
            fmpq_t t;
            fmpq_init(t);
            for (r = 0; r < u; r++)
            {
                fmpq_set(t, L.cre);
                fmpq_neg(L.cre, L.cim);
                fmpq_set(L.cim, t);
            }
            fmpq_clear(t);
            L.unit = (L.unit + u) & 3;
            ok = _acf_local_side(side, &L, H, inf, part, beta, bpi);
        }
        else if (u == 0)
            ok = _acf_local_side(side, &L, H, inf, part, beta, bpi);
    }

    fmpq_poly_clear(Q);
    acf_clear(H);
    acf_clear(S);
    acf_clear(y);
    fmpq_clear(L.cre);
    fmpq_clear(L.cim);
    return ok;
}

/* The sum of the terms of T (all with k = 0) as a ball at precision w
   with the sign of the sum (or an exact zero) */
static void
_acf_pt_get_arb(arb_t res, const _acf_pt_t T, slong w)
{
    arb_zero(res);
    if (T->n != 0)
        _acf_sum_ball(res, T->t, T->n, w);
}

/* the rounding of a ball r to prec bits in the strict mode rnd: the
   flag, or ARF_ELEM_UNCERTAIN when its endpoints round differently */
static int
_acf_ball_round(arf_t res, const arb_t r, slong prec, arf_rnd_t rnd)
{
    arf_t lo, hi;
    slong wp;
    int s1, s2, st;

    if (!arb_is_finite(r))
        return ARF_ELEM_UNCERTAIN;
    if (arb_is_exact(r))
        return arf_set_round(res, arb_midref(r), prec, rnd);
    if (arb_contains_zero(r))
        return ARF_ELEM_UNCERTAIN;

    wp = FLINT_MAX(prec, (slong) arf_bits(arb_midref(r))) + 64;
    arf_init(lo);
    arf_init(hi);
    arb_get_lbound_arf(lo, r, wp);
    arb_get_ubound_arf(hi, r, wp);
    s1 = arf_set_round(lo, lo, prec, rnd);
    s2 = arf_set_round(hi, hi, prec, rnd);
    (void) s1; (void) s2;
    if (arf_equal(lo, hi))
    {
        arf_swap(res, lo);
        st = 1;
    }
    else
        st = ARF_ELEM_UNCERTAIN;
    arf_clear(lo);
    arf_clear(hi);
    return st;
}

/*
    x^y near its exact cases: x = z0 (1 + h) with z0 = 2^K i^u and h
    tiny, where x^y = z0^y sum_j binom(y, j) h^j, the sum exact as term
    lists (multiplied by N!) and z0^y exact when it is dyadic (y real with
    K y and u y integers), otherwise a ball whose part is exactly zero
    when its angle K log(2) Im(y) + pi u Re(y) / 2 is a multiple of pi/2;
    and y tiny, where x^y = 1 + expm1(y log x) with expm1(t) = expm1(t_r)
    cos(t_i) - 2 sin(t_i / 2)^2 + i exp(t_r) sin(t_i) in ball arithmetic
    (whose parts have their own relative accuracy).  A part with a dyadic
    leading term near a rounding boundary is rounded from the exact side
    of that boundary, the others from balls at increasing precision.
    Returns 1 with the result in res (and the flags), 0 when neither
    applies or the precision limit is reached.
*/
int
_acf_pow_near(acf_t res, int * status, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x), yr = acf_realref(y), yi = acf_imagref(y);
    arf_rnd_t r1 = (rnd >= ARF_RND_FAST) ? arf_rnd_relaxed_to_strict(rnd) : rnd;
    slong N = ACF_NEAR_TERMS, K, w, j;
    int u, ok = 0, s[2], i;
    acf_t h, z, t1;
    arf_t m0, R[2];
    mag_t hm, ym;

    if (arf_is_zero(a) && arf_is_zero(b))
        return 0;

    acf_init(h);
    acf_init(z);
    acf_init(t1);
    arf_init(m0);
    arf_init(R[0]);
    arf_init(R[1]);
    mag_init(hm);
    mag_init(ym);

    /* |y| + 1 */
    arf_get_mag(ym, yr);
    arf_get_mag(hm, yi);
    mag_hypot(ym, ym, hm);
    mag_add_ui(ym, ym, 1);

    /* z0 = +-2^K or +-i 2^K: the larger part rounded to one bit */
    u = (arf_cmpabs(a, b) >= 0) ? ((arf_sgn(a) > 0) ? 0 : 2) : ((arf_sgn(b) > 0) ? 1 : 3);
    arf_abs(m0, (u % 2 == 0) ? a : b);
    arf_set_round(m0, m0, 1, ARF_RND_NEAR);
    if (!COEFF_IS_MPZ(ARF_EXP(m0)))
    {
        K = ARF_EXP(m0) - 1;
        /* h = x / z0 - 1 (exactly) */
        acf_set(h, x);
        for (i = 0; i < 4 - u; i++)
        {
            arf_swap(acf_realref(h), acf_imagref(h));
            arf_neg(acf_realref(h), acf_realref(h));
        }
        arf_mul_2exp_si(acf_realref(h), acf_realref(h), -K);
        arf_mul_2exp_si(acf_imagref(h), acf_imagref(h), -K);
        arf_sub_ui(acf_realref(h), acf_realref(h), 1, ARF_PREC_EXACT, ARF_RND_DOWN);
        /* (from here u counts the signed quarter turns of the principal
           argument of z0: -1 for -i, and -2 for -1 below the cut, where
           the principal log x = log z0 + log(1 + h) takes arg z0 = -pi) */
        if (u == 3)
            u = -1;
        else if (u == 2 && arf_sgn(b) < 0)
            u = -2;
        arf_get_mag(hm, acf_realref(h));
        {
            mag_t v;
            mag_init(v);
            arf_get_mag(v, acf_imagref(h));
            mag_hypot(hm, hm, v);
            mag_clear(v);
        }
        mag_mul(hm, hm, ym);        /* rho |h| */

        if (mag_cmp_2exp_si(hm, -prec - 32) <= 0 && !mag_is_zero(hm))
        {
            _acf_cpt_t E, Hc, Hj, T, C;
            arf_t one, c;
            fmpz_t L, q;
            mag_t rem;
            int c0exact, zpart = -1;

            _acf_cpt_init(E);
            _acf_cpt_init(Hc);
            _acf_cpt_init(Hj);
            _acf_cpt_init(T);
            _acf_cpt_init(C);
            arf_init(one);
            arf_init(c);
            fmpz_init(L);
            fmpz_init(q);
            mag_init(rem);
            arf_one(one);

            /* E = N! sum_j binom(y, j) h^j: coefficients N!/j! prod_{m<j} (y - m) */
            fmpz_fac_ui(L, N);
            _acf_cpt_set(Hc, acf_realref(h), acf_imagref(h));
            _acf_cpt_zero(Hj);
            _acf_pt_push(Hj->re, one, 0, NULL, NULL);
            arf_one(acf_realref(z));
            arf_zero(acf_imagref(z));         /* prod (y - m) */
            for (j = 0; j <= N; j++)
            {
                if (j > 0)
                {
                    acf_set(t1, y);
                    arf_sub_si(acf_realref(t1), acf_realref(t1), j - 1, ARF_PREC_EXACT, ARF_RND_DOWN);
                    acf_mul(z, z, t1, ARF_PREC_EXACT, ARF_RND_DOWN);
                    _acf_cpt_zero(T);
                    _acf_cpt_addmul(T, Hj, Hc);
                    _acf_cpt_swap(T, Hj);
                }
                fmpz_fac_ui(q, j);
                fmpz_divexact(q, L, q);
                arf_set_fmpz(c, q);
                _acf_cpt_set(C, acf_realref(z), acf_imagref(z));
                _acf_cpt_zero(T);
                _acf_cpt_addmul(T, C, Hj);
                _acf_cpt_addscal(E, T, c, 0, 0);
            }

            /* the omitted terms: 2 (rho |h|)^(N + 1), times |z0^y| below
               (none for an integer 0 <= y <= N: the binomial sum ends) */
            mag_pow_ui(rem, hm, N + 1);
            mag_mul_2exp_si(rem, rem, 1);
            if (arf_is_zero(yi) && arf_is_int(yr) && arf_sgn(yr) >= 0 && arf_cmp_si(yr, N) <= 0)
                mag_zero(rem);

            /* z0^y: exact when y is real with K y and u y integers */
            c0exact = 0;
            if (arf_is_zero(yi))
            {
                arf_t ky, uy;
                arf_init(ky);
                arf_init(uy);
                arf_mul_si(ky, yr, K, ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_mul_si(uy, yr, u, ARF_PREC_EXACT, ARF_RND_DOWN);
                if (arf_is_int(ky) && arf_is_int(uy) && arf_cmpabs_2exp_si(ky, FLINT_BITS - 4) < 0)
                {
                    /* 2^(K y) i^(u y) */
                    fmpz_t e;
                    slong m;
                    fmpz_init(e);
                    arf_get_fmpz(e, ky, ARF_RND_DOWN);
                    m = arf_get_si(uy, ARF_RND_DOWN);
                    m = ((m % 4) + 4) % 4;
                    arf_one(acf_realref(z));
                    arf_zero(acf_imagref(z));
                    arf_mul_2exp_fmpz(acf_realref(z), acf_realref(z), e);
                    while (m-- > 0)
                    {
                        arf_swap(acf_realref(z), acf_imagref(z));
                        arf_neg(acf_realref(z), acf_realref(z));
                    }
                    c0exact = 1;
                    fmpz_clear(e);
                }
                arf_clear(ky);
                arf_clear(uy);
            }
            if (!c0exact && (K == 0 || arf_is_zero(yi)))
            {
                /* the angle pi u Re(y) / 2: a part exactly zero when u y_r
                   is an integer (zpart: that part) */
                arf_t uy;
                arf_init(uy);
                arf_mul_si(uy, yr, u, ARF_PREC_EXACT, ARF_RND_DOWN);
                if (arf_is_int(uy) && arf_cmpabs_2exp_si(uy, FLINT_BITS - 4) < 0)
                {
                    slong m = arf_get_si(uy, ARF_RND_DOWN);
                    zpart = (int) ((((m % 2) + 2) % 2) == 0);   /* even: imaginary part zero */
                }
                arf_clear(uy);
            }

            for (w = prec + 64; w <= 4 * prec + 65536 && !ok; w *= 2)
            {
                arb_t ER, EI, c0r, c0i, v, tt;
                acb_t c0;
                arb_init(ER);
                arb_init(EI);
                arb_init(c0r);
                arb_init(c0i);
                arb_init(v);
                arb_init(tt);
                acb_init(c0);

                _acf_pt_get_arb(ER, E->re, w);
                _acf_pt_get_arb(EI, E->im, w);

                if (c0exact)
                {
                    arb_set_arf(c0r, acf_realref(z));
                    arb_set_arf(c0i, acf_imagref(z));
                }
                else
                {
                    /* exp(y (K log 2 + i pi u / 2)) */
                    acb_t l, yy;
                    acb_init(l);
                    acb_init(yy);
                    arb_const_log2(acb_realref(l), w);
                    arb_mul_si(acb_realref(l), acb_realref(l), K, w);
                    arb_const_pi(acb_imagref(l), w);
                    arb_mul_si(acb_imagref(l), acb_imagref(l), u, w);
                    arb_mul_2exp_si(acb_imagref(l), acb_imagref(l), -1);
                    arb_set_arf(acb_realref(yy), yr);
                    arb_set_arf(acb_imagref(yy), yi);
                    acb_mul(l, l, yy, w);
                    if (zpart >= 0)
                    {
                        /* |c0| i^m exactly in direction */
                        arb_exp(tt, acb_realref(l), w);
                        acb_zero(c0);
                        arb_set(zpart ? acb_realref(c0) : acb_imagref(c0), tt);
                        {
                            arf_t uy;
                            slong m;
                            arf_init(uy);
                            arf_mul_si(uy, yr, u, ARF_PREC_EXACT, ARF_RND_DOWN);
                            m = arf_get_si(uy, ARF_RND_DOWN);
                            m = ((m % 4) + 4) % 4;
                            if (m == 2 || m == 3)
                                acb_neg(c0, c0);
                            arf_clear(uy);
                        }
                    }
                    else
                        acb_exp(c0, l, w);
                    arb_set(c0r, acb_realref(c0));
                    arb_set(c0i, acb_imagref(c0));
                    acb_clear(l);
                    acb_clear(yy);
                }

                ok = 1;
                for (i = 0; i < 2 && ok; i++)
                {
                    /* Re: c0r ER - c0i EI, Im: c0r EI + c0i ER (over N!) */
                    mag_t cm;
                    arb_zero(v);
                    if (!arb_is_zero(c0r))
                        arb_addmul(v, c0r, i ? EI : ER, w);
                    if (!arb_is_zero(c0i))
                    {
                        arb_mul(tt, c0i, i ? ER : EI, w);
                        if (i)
                            arb_add(v, v, tt, w);
                        else
                            arb_sub(v, v, tt, w);
                    }
                    arb_div_fmpz(v, v, L, w);
                    mag_init(cm);
                    {
                        mag_t cn;
                        mag_init(cn);
                        arb_get_mag(cm, c0r);
                        arb_get_mag(cn, c0i);
                        mag_hypot(cm, cm, cn);
                        mag_clear(cn);
                    }
                    mag_mul(cm, cm, rem);
                    arb_add_error_mag(v, cm);
                    mag_clear(cm);

                    s[i] = _acf_ball_round(R[i], v, prec, r1);
                    if (s[i] == ARF_ELEM_UNCERTAIN && c0exact && arb_rel_accuracy_bits(v) > prec + 40)
                    {
                        /* the exact side of the boundary near the value:
                           Part(z0^y E) - N! beta, with the omitted terms */
                        _acf_cpt_t P;
                        arf_t beta;
                        mag_t rr;
                        int side;
                        _acf_cpt_init(P);
                        arf_init(beta);
                        mag_init(rr);
                        _acf_round_boundary(beta, arb_midref(v), prec, r1);
                        _acf_cpt_set(C, acf_realref(z), acf_imagref(z));
                        _acf_cpt_zero(P);
                        _acf_cpt_addmul(P, C, E);
                        arf_set_fmpz(c, L);
                        arf_mul(c, c, beta, ARF_PREC_EXACT, ARF_RND_DOWN);
                        arf_neg(c, c);
                        _acf_pt_push(i ? P->im : P->re, c, 0, NULL, NULL);
                        {
                            mag_t rv;
                            mag_init(rv);
                            arf_get_mag(rr, acf_realref(z));
                            arf_get_mag(rv, acf_imagref(z));
                            mag_hypot(rr, rr, rv);
                            mag_clear(rv);
                        }
                        mag_mul(rr, rr, rem);
                        mag_mul_fmpz(rr, rr, L);
                        if (_acf_pt_sign(&side, i ? P->im : P->re, rr))
                            s[i] = _acf_round_side(R[i], beta, side, prec, r1);
                        _acf_cpt_clear(P);
                        arf_clear(beta);
                        mag_clear(rr);
                    }
                    if (s[i] == ARF_ELEM_UNCERTAIN)
                        ok = 0;
                }

                arb_clear(ER);
                arb_clear(EI);
                arb_clear(c0r);
                arb_clear(c0i);
                arb_clear(v);
                arb_clear(tt);
                acb_clear(c0);
            }

            _acf_cpt_clear(E);
            _acf_cpt_clear(Hc);
            _acf_cpt_clear(Hj);
            _acf_cpt_clear(T);
            _acf_cpt_clear(C);
            arf_clear(one);
            arf_clear(c);
            fmpz_clear(L);
            fmpz_clear(q);
            mag_clear(rem);
        }
    }

    if (!ok)
    {
        /* y tiny: 1 + expm1(y log x) */
        mag_t lm;
        mag_init(lm);
        arf_get_mag(hm, yr);
        arf_get_mag(lm, yi);
        mag_hypot(hm, hm, lm);
        mag_clear(lm);
        if (!mag_is_zero(hm) && mag_cmp_2exp_si(hm, -prec - 64) <= 0)
        {
            for (w = prec + 64; w <= 4 * prec + 65536 && !ok; w *= 2)
            {
                acb_t l, yy;
                arb_t er, ei, v, tt;
                acb_init(l);
                acb_init(yy);
                arb_init(er);
                arb_init(ei);
                arb_init(v);
                arb_init(tt);

                arb_set_arf(acb_realref(l), a);
                arb_set_arf(acb_imagref(l), b);
                acb_log(l, l, w);
                arb_set_arf(acb_realref(yy), yr);
                arb_set_arf(acb_imagref(yy), yi);
                acb_mul(l, l, yy, w);       /* t = y log x */

                if (acb_is_finite(l) && arb_is_nonzero(acb_realref(l)) && arb_is_nonzero(acb_imagref(l))
                    && mag_cmp_2exp_si(arb_radref(acb_realref(l)), -16) < 0)
                {
                    /* Re: expm1(t_r) cos(t_i) - 2 sin(t_i / 2)^2 */
                    arb_expm1(er, acb_realref(l), w);
                    arb_cos(tt, acb_imagref(l), w);
                    arb_mul(er, er, tt, w);
                    arb_mul_2exp_si(tt, acb_imagref(l), -1);
                    arb_sin(tt, tt, w);
                    arb_sqr(tt, tt, w);
                    arb_mul_2exp_si(tt, tt, 1);
                    arb_sub(er, er, tt, w);
                    /* Im: exp(t_r) sin(t_i) */
                    arb_exp(ei, acb_realref(l), w);
                    arb_sin(tt, acb_imagref(l), w);
                    arb_mul(ei, ei, tt, w);

                    ok = 1;
                    /* the real part 1 + er: from the sign of er when it is
                       below the rounding scale */
                    if (arb_is_nonzero(er) && mag_cmp_2exp_si(arb_radref(er), 0) < 0 && arf_cmpabs_2exp_si(arb_midref(er), -prec - 8) < 0
                        && (arb_is_positive(er) || arb_is_negative(er)))
                    {
                        arf_one(m0);
                        s[0] = _acf_round_side(R[0], m0, arb_sgn_nonzero(er), prec, r1);
                    }
                    else
                    {
                        arb_add_ui(v, er, 1, w);
                        s[0] = _acf_ball_round(R[0], v, prec, r1);
                    }
                    s[1] = _acf_ball_round(R[1], ei, prec, r1);
                    if (s[0] == ARF_ELEM_UNCERTAIN || s[1] == ARF_ELEM_UNCERTAIN)
                        ok = 0;
                }

                acb_clear(l);
                acb_clear(yy);
                arb_clear(er);
                arb_clear(ei);
                arb_clear(v);
                arb_clear(tt);
            }
        }
    }

    if (ok)
    {
        arf_swap(acf_realref(res), R[0]);
        arf_swap(acf_imagref(res), R[1]);
        *status = s[0] | (s[1] << 1);
    }

    acf_clear(h);
    acf_clear(z);
    acf_clear(t1);
    arf_clear(m0);
    arf_clear(R[0]);
    arf_clear(R[1]);
    mag_clear(hm);
    mag_clear(ym);
    return ok;
}
