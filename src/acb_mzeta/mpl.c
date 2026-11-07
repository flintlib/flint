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

   Power series at 0: with F_l(t) = G(a_l, ..., a_L; t) = sum_n f_l(n) t^n,
   f_l(n) = f_{l+1}(n)/n if a_l = 0, and f_l(n) = -S_l(n-1)/(n a_l) with
   S_l(n) = S_l(n-1)/a_l + f_{l+1}(n) otherwise.  Scaled by z^n
   (g = f z^n, T = S z^n), for n >= 1:

     g_l(n) = g_{l+1}(n)/n                        (a_l = 0)
     g_l(n) = -z T_l(n-1) / (n a_l)               (a_l != 0)
     T_l(n) = (z/a_l) T_l(n-1) + g_{l+1}(n),   g_{L+1}(n) = 0,
     A_l(n) = A_l(n-1) + g_l(n)                   (A_l(N) -> F_l(z)),

   where g_{l+1}(n) is expanded down to the next nonzero letter, so that
   X(n) = M(n) X(n-1) with the T_l (a_l != 0) and the A_l, lower
   triangular in the order of decreasing l; X(0) has T_L = 1.  The rows
   are multiplied by n^C (C the longest run of zeros plus 1) and the
   product is computed by binary splitting in ball arithmetic (exact while
   the entries are exact and short).

   Truncation: the f_l are majorized by the coefficients for the letters
   |a_l| (positive series), and with R = min |a_l| over the nonzero
   letters, rho = R N/(N+1): F~(rho) <= log(R/(R - rho))^d = log(N+1)^d
   (an integral over a letter 0 does not increase the bound of a series
   with nonnegative coefficients vanishing at 0, one over a letter b
   multiplies it by at most log(R/(R - rho))), so
   |f(n) z^n| <= log(N+1)^d (|z|/rho)^n and the tail beyond N is at most
   log(N+1)^d e (1 + 1/N) q^(N+1) / (1 - q (1 + 1/N)), q = |z|/R. */

typedef struct
{
    acb_ptr e;
    slong n;
}
cmat_struct;

#define CE(M, i, j) ((M)->e + (i) * (M)->n + (j))

static void
cmat_init(cmat_struct * M, slong n)
{
    M->n = n;
    M->e = _acb_vec_init(n * n);
}

static void
cmat_clear(cmat_struct * M)
{
    _acb_vec_clear(M->e, M->n * M->n);
}

static void
cmat_mul(cmat_struct * C, const cmat_struct * A, const cmat_struct * B, slong prec)
{
    slong n = A->n, i, j, l;
    for (i = 0; i < n; i++)
        for (j = 0; j <= i; j++)
        {
            acb_ptr c = CE(C, i, j);
            acb_zero(c);
            for (l = j; l <= i; l++)
                if (!acb_is_zero(CE(A, i, l)) && !acb_is_zero(CE(B, l, j)))
                    acb_addmul(c, CE(A, i, l), CE(B, l, j), prec);
        }
}

typedef struct
{
    slong L, d, dim, C;
    slong * tidx;       /* state index of T_l (a_l != 0), else -1 */
    slong * next;       /* the first l' >= l with a_{l'} != 0 */
    acb_ptr zb;         /* z / a_l for nonzero letters */
    acb_ptr mzb;        /* -z / a_l */
}
gctx;

/* state indices: T_l for nonzero letters in order of decreasing l, then
   the A_l in order of decreasing l */
static slong
aidx(const gctx * G, slong l)
{
    return G->d + (G->L - 1 - l);
}

static void
leaf(cmat_struct * Q, const gctx * G, ulong n, slong prec)
{
    slong L = G->L, l;
    acb_t t;
    fmpz_t p;

    acb_init(t);
    fmpz_init(p);
    _acb_vec_zero(Q->e, Q->n * Q->n);

    for (l = 0; l < L; l++)
    {
        /* g_l(n) = (1/n^c) (-z/(n a_l')) T_l'(n-1), l' = next(l), c = l' - l;
           times n^C */
        slong lp = G->next[l], c = lp - l + 1, r = aidx(G, l);
        fmpz_set_ui(p, n);
        fmpz_pow_ui(p, p, G->C - c);
        acb_mul_fmpz(t, G->mzb + lp, p, prec);
        acb_set(CE(Q, r, G->tidx[lp]), t);
        /* A_l(n-1) n^C */
        fmpz_set_ui(p, n);
        fmpz_pow_ui(p, p, G->C);
        acb_set_fmpz(CE(Q, r, r), p);

        if (G->tidx[l] >= 0)
        {
            /* T_l(n) = (z/a_l) T_l(n-1) + g_{l+1}(n) */
            slong rt = G->tidx[l];
            acb_mul_fmpz(CE(Q, rt, rt), G->zb + l, p, prec);
            if (l + 1 < L)
            {
                slong lq = G->next[l + 1], cq = lq - (l + 1) + 1;
                fmpz_set_ui(p, n);
                fmpz_pow_ui(p, p, G->C - cq);
                acb_mul_fmpz(CE(Q, rt, G->tidx[lq]), G->mzb + lq, p, prec);
            }
        }
    }

    acb_clear(t);
    fmpz_clear(p);
}

static void
bsplit(cmat_struct * P, fmpz_t Den, const gctx * G, ulong a, ulong b, slong prec)
{
    if (b - a == 1)
    {
        leaf(P, G, b, prec);
        fmpz_set_ui(Den, b);
        fmpz_pow_ui(Den, Den, G->C);
    }
    else
    {
        ulong m = a + (b - a) / 2;
        cmat_struct P1, P2;
        fmpz_t D2;
        cmat_init(&P1, G->dim);
        cmat_init(&P2, G->dim);
        fmpz_init(D2);
        bsplit(&P1, Den, G, a, m, prec);
        bsplit(&P2, D2, G, m, b, prec);
        cmat_mul(P, &P2, &P1, prec);
        fmpz_mul(Den, Den, D2);
        cmat_clear(&P1);
        cmat_clear(&P2);
        fmpz_clear(D2);
    }
}

/* the first column of the product, the lowest segment as a vector */
static void
bsplit_col(acb_ptr v, fmpz_t Den, const gctx * G, ulong a, ulong b, slong col, slong prec)
{
    slong n = G->dim, i, l;
    if (b - a <= 4)
    {
        cmat_struct P;
        cmat_init(&P, n);
        bsplit(&P, Den, G, a, b, prec);
        for (i = 0; i < n; i++)
            acb_set(v + i, CE(&P, i, col));
        cmat_clear(&P);
    }
    else
    {
        ulong m = a + (b - a) / 2;
        cmat_struct P2;
        acb_ptr v1;
        fmpz_t D2;
        cmat_init(&P2, n);
        v1 = _acb_vec_init(n);
        fmpz_init(D2);
        bsplit_col(v1, Den, G, a, m, col, prec);
        bsplit(&P2, D2, G, m, b, prec);
        for (i = 0; i < n; i++)
        {
            acb_zero(v + i);
            for (l = 0; l <= i; l++)
                if (!acb_is_zero(CE(&P2, i, l)) && !acb_is_zero(v1 + l))
                    acb_addmul(v + i, CE(&P2, i, l), v1 + l, prec);
        }
        fmpz_mul(Den, Den, D2);
        cmat_clear(&P2);
        _acb_vec_clear(v1, n);
        fmpz_clear(D2);
    }
}

int
_acb_mzeta_G_suffixes(acb_ptr res, acb_srcptr a, slong L, const acb_t z, slong prec)
{
    gctx G;
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
    G.tidx = flint_malloc(sizeof(slong) * L);
    G.next = flint_malloc(sizeof(slong) * L);
    G.zb = _acb_vec_init(L);
    G.mzb = _acb_vec_init(L);

    mag_init(R);
    mag_init(q);
    mag_init(tail);
    mag_init(t);

    /* nonzero letters (as balls: a letter whose ball contains 0 is not
       allowed), runs of zeros, R = min |a_l| */
    d = 0;
    run = 0;
    G.C = 1;
    mag_inf(R);
    for (l = L - 1; l >= 0; l--)
    {
        if (acb_is_zero(a + l))
        {
            G.tidx[l] = -1;
            G.next[l] = G.next[l + 1];
            run++;
            G.C = FLINT_MAX(G.C, run + 1);
        }
        else
        {
            if (acb_contains_zero(a + l))
            {
                ok = 0;
                break;
            }
            G.tidx[l] = d++;
            G.next[l] = l;
            run = 0;
            acb_get_mag_lower(t, a + l);
            mag_min(R, R, t);
        }
    }
    G.d = d;
    G.dim = d + L;

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
        /* N: log(N+1)^d e q^(N+1) / (1 - q (1 + 1/N)) below 2^-(prec + 4) */
        N = 8;
        while ((N + 1) * -log2(qd) - d * log2(log(N + 1.0) + 1e-300) - 2.0 < prec + 8)
            N += 1 + N / 16;

        /* rigorous bound */
        {
            arb_t u, v, w;
            arb_init(u);
            arb_init(v);
            arb_init(w);
            arb_set_ui(u, N + 1);
            arb_log(u, u, 30);
            arb_pow_ui(u, u, d, 30);
            arf_set_mag(arb_midref(v), q);
            mag_zero(arb_radref(v));
            arb_pow_ui(w, v, N + 1, 30);
            arb_mul(u, u, w, 30);
            /* (1 + 1/N)^(N+1) <= e (N+1)/N */
            arb_const_e(w, 30);
            arb_mul(u, u, w, 30);
            arb_mul_ui(u, u, N + 1, 30);
            arb_div_ui(u, u, N, 30);
            arb_mul_ui(w, v, N + 1, 30);
            arb_div_ui(w, w, N, 30);
            arb_sub_ui(w, w, 1, 30);
            arb_neg(w, w);
            if (arb_is_positive(w))
            {
                arb_div(u, u, w, 30);
                arb_get_mag(tail, u);
            }
            else
                mag_inf(tail);
            arb_clear(u);
            arb_clear(v);
            arb_clear(w);
        }

        wp = prec + 2 * FLINT_BIT_COUNT(N) + 2 * FLINT_BIT_COUNT(G.dim) + 16;

        for (l = 0; l < L; l++)
            if (G.tidx[l] >= 0)
            {
                acb_div(G.zb + l, z, a + l, wp);
                acb_neg(G.mzb + l, G.zb + l);
            }

        {
            acb_ptr v = _acb_vec_init(G.dim);
            fmpz_t Den;
            arb_t Dr;
            fmpz_init(Den);
            arb_init(Dr);
            /* X(0) = e_{T_{L-1}}: T of the last letter is the first state */
            bsplit_col(v, Den, &G, 0, N, G.tidx[L - 1], wp);
            arb_set_fmpz(Dr, Den);
            for (l = 0; l < L; l++)
            {
                acb_div_arb(res + l, v + aidx(&G, l), Dr, wp);
                acb_add_error_mag(res + l, tail);
            }
            acb_one(res + L);
            _acb_vec_clear(v, G.dim);
            fmpz_clear(Den);
            arb_clear(Dr);
        }
    }

    flint_free(G.tidx);
    flint_free(G.next);
    _acb_vec_clear(G.zb, L);
    _acb_vec_clear(G.mzb, L);
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
            acb_set_round(res, r + 0, prec);
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
