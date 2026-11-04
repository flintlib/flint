/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mpn_extras.h"
#include "impl.h"

/* Multiple zeta values by Akhilesh's double tails [Akh2017].

   For a binary word w = e_1 ... e_k (letters x0 = dt/t, x1 = dt/(1-t))
   the double tails are

     z(w)_{m,n} = int_{1 > t_1 > ... > t_k > 0} (1 - t_1)^m w(t) t_k^n,

   z(w)_{0,0} the iterated integral (the MZV if w is admissible).  They
   satisfy (integrating out t_1 or t_k):

     z(x1 u)_{m,n} = z(u)_{m,n} / m,            z(u x0)_{m,n} = z(u)_{m,n} / n,
     z(x0 u)_{m,n} - z(x0 u)_{m+1,n} = z(u)_{m+1,n} / (m + 1),
     z(u x1)_{m,n} - z(u x1)_{m,n+1} = z(u)_{m,n+1} / (n + 1),
     z()_{m,n} = m! n! / (m + n)!,

   so for u' = x0 u x1, stepping along the diagonal,

     z(u')_{m,n} = z(u')_{m+1,n+1} + z(u x1)_{m+1,n} / (m + 1) + z(x0 u)_{m+1,n+1} / (n + 1).

   The subword w[i..j) only occurs with m - n = d_i, the number of x0 in
   w[0..i), so level n of the recursion holds one value
   V_n = z(u)_{n + d, n} per state (u, d) = (w[i..j), d_i); the levels are
   computed from n = N down to 0 (level N + 1 taken as zero), the states of
   a level by increasing length.  Several words are evaluated together,
   sharing the states they have in common (the merged indices of star
   values, all words of a weight).  At level 0 only the states reachable
   from the words at level 0 are needed.

   Arithmetic: fixed point with B fraction bits (one integer limb on top),
   truncating divisions, so that every computed value is a lower bound
   (all quantities are nonnegative) and its error is counted in ulps
   (exact integer bookkeeping).  The empty word values are propagated as
   4^n m! n! / (m + n)!, whose step ratio is about 1, and shifted down.

   Truncation: z(s)_{m,n} <= max_t (1 - t)^(m - a) t^(n - b) for a subword
   with a letters x1 and b letters x0 (bound x1 by 1/(1 - t_1), x0 by 1/t_k,
   and the integrand by its value at t_1 = t_k; the simplex has volume
   <= 1).  If G bounds the dropped level N + 1 values, the error at level n
   >= 1 is at most G (n + 2)/n times that at level n + 1 (the diagonal step
   contributes (1 + 1/(n+1)) G_{n+1} plus 1/(m+1) <= 1/(n+1) of a shorter
   word's error at the same level; the other rules divide by m, n >= 1),
   and at level 0 at most 6 times that at level 1, so in total
   3 (N + 1) (N + 2) G. */

/* the word of a composition: x0^(s_1 - 1) x1 ... x0^(s_r - 1) x1 */
slong
_acb_mzeta_mzv_word(unsigned char * w, const ulong * s, slong r)
{
    slong i, k = 0;
    ulong t;
    for (i = 0; i < r; i++)
    {
        for (t = 1; t < s[i]; t++)
            w[k++] = 0;
        w[k++] = 1;
    }
    return k;
}

/* res = floor(x / d), error count err -> ceil(err / d) + 1 */
static inline void
fdiv(nn_ptr res, nn_srcptr x, slong len, ulong d)
{
    mpn_divrem_1(res, 0, x, len, d);
}

static inline ulong
err_div(ulong err, ulong d)
{
    return (err + d - 1) / d + 1;
}

/* mag bound for max_t (1 - t)^a t^b = a^a b^b / (a + b)^(a + b), a, b >= 0 */
static void
mag_beta_max(mag_t res, ulong a, ulong b)
{
    arb_t t, u;
    if (a == 0 || b == 0)
    {
        mag_one(res);
        return;
    }
    arb_init(t);
    arb_init(u);
    /* exp(a log(a/(a+b)) + b log(b/(a+b))) */
    arb_set_ui(t, a);
    arb_div_ui(t, t, a + b, 30);
    arb_log(t, t, 30);
    arb_mul_ui(t, t, a, 30);
    arb_set_ui(u, b);
    arb_div_ui(u, u, a + b, 30);
    arb_log(u, u, 30);
    arb_mul_ui(u, u, b, 30);
    arb_add(t, t, u, 30);
    arb_exp(t, t, 30);
    arb_get_mag(res, t);
    arb_clear(t);
    arb_clear(u);
}

#define KDIV 4

/* a state: the subword u = word[i..i+L) of some input word, with offset d */
typedef struct
{
    const unsigned char * u;
    slong L;
    slong d;
    slong c1, c2;           /* children (see below), -1 if none */
    int type;               /* 0 empty, 1 x1 first, 2 x0 last, 3 diagonal */
    int need0;              /* needed at level 0 */
}
mzv_state;

static int
state_cmp(const void * a, const void * b)
{
    const mzv_state * x = a, * y = b;
    if (x->L != y->L)
        return (x->L < y->L) ? -1 : 1;
    if (x->d != y->d)
        return (x->d < y->d) ? -1 : 1;
    return memcmp(x->u, y->u, x->L);
}

static slong
state_find(const mzv_state * S, slong ns, const unsigned char * u, slong L, slong d)
{
    mzv_state key;
    const mzv_state * r;
    key.u = u;
    key.L = L;
    key.d = d;
    r = bsearch(&key, S, ns, sizeof(mzv_state), state_cmp);
    if (r == NULL)
        flint_throw(FLINT_ERROR, "mzv_akhilesh: missing state\n");
    return r - S;
}

/* truncation bound for dropping level N + 1: 3 (N + 1)(N + 2) max U over
   the nonempty states */
static void
trunc_bound(mag_t res, const mzv_state * S, slong ns, slong N)
{
    slong s, t;
    mag_t u;
    mag_init(u);
    mag_zero(res);
    for (s = 0; s < ns; s++)
    {
        slong a = 0, b = 0, m = N + 1 + S[s].d, n = N + 1;
        if (S[s].L == 0)
            continue;
        for (t = 0; t < S[s].L; t++)
        {
            if (S[s].u[t])
                a++;
            else
                b++;
        }
        if (m < a || n < b)
        {
            mag_inf(res);
            break;
        }
        mag_beta_max(u, m - a, n - b);
        mag_max(res, res, u);
    }
    mag_mul_ui(res, res, N + 1);
    mag_mul_ui(res, res, N + 2);
    mag_mul_ui(res, res, 3);
    mag_clear(u);
}

/* e = floor of a lower bound for Ehat(n, d) = 4^n (n+d)! n! / (2n+d)!
   in fixed point with B fraction bits (len + 1 limbs), and its error in
   ulps */
static void
ehat_init(nn_ptr e, ulong * err, slong n, slong d, slong B, slong len)
{
    arb_t t, u;
    fmpz_t z;
    mag_t r;

    arb_init(t);
    arb_init(u);
    fmpz_init(z);
    mag_init(r);

    arb_fac_ui(t, n + d, B + 64);
    arb_fac_ui(u, n, B + 64);
    arb_mul(t, t, u, B + 64);
    arb_fac_ui(u, 2 * n + d, B + 64);
    arb_div(t, t, u, B + 64);
    arb_mul_2exp_si(t, t, 2 * n);
    /* floor of the lower bound; error = ceil(2 rad) + 1 ulps */
    arb_get_lbound_arf(arb_midref(u), t, B + 64);
    arf_get_fmpz_fixed_si(z, arb_midref(u), -B);
    if (fmpz_sgn(z) < 0)
        fmpz_zero(z);
    fmpz_get_ui_array(e, len + 1, z);
    mag_mul_2exp_si(r, arb_radref(t), B + 1);
    *err = (ulong) mag_get_d(r) + 3;

    arb_clear(t);
    arb_clear(u);
    fmpz_clear(z);
    mag_clear(r);
}


void
_acb_mzeta_mzv_akhilesh_vec(arb_ptr res, const unsigned char * const * words,
    const slong * lens, slong nw, slong prec)
{
    slong N, B, nl, len, la, i, j, n, ndel, ns, kmax, s, nj;
    slong jobs[KDIV];
    nn_ptr dbuf[KDIV];
    ulong dd[KDIV], dinv[KDIV];
    mzv_state * S;
    slong * root;
    nn_ptr V, W, E, tmp;
    ulong * eV, * eW, * eE, emax;
    mag_t tb;

    kmax = 0;
    ns = 0;
    for (i = 0; i < nw; i++)
    {
        kmax = FLINT_MAX(kmax, lens[i]);
        ns += (lens[i] + 1) * (lens[i] + 2) / 2;
    }

    /* all states of all words, deduplicated */
    S = flint_malloc(sizeof(mzv_state) * FLINT_MAX(ns, 1));
    ns = 0;
    for (i = 0; i < nw; i++)
    {
        const unsigned char * w = words[i];
        slong k = lens[i], a, b, d = 0;
        for (a = 0; a <= k; a++)
        {
            for (b = a; b <= k; b++)
            {
                S[ns].u = w + a;
                S[ns].L = b - a;
                S[ns].d = d;
                S[ns].need0 = 0;
                ns++;
            }
            if (a < k && w[a] == 0)
                d++;
        }
    }
    qsort(S, ns, sizeof(mzv_state), state_cmp);
    j = 0;
    for (i = 0; i < ns; i++)
        if (j == 0 || state_cmp(S + i, S + j - 1) != 0)
            S[j++] = S[i];
    ns = j;

    /* children: x1 first: (u[1..], d) same level;  x0 last: (u[..-1], d)
       same level;  diagonal: c1 = (u[1..], d + 1) same level,
       c2 = (u[..-1], d) previous level */
    for (s = 0; s < ns; s++)
    {
        const unsigned char * u = S[s].u;
        slong L = S[s].L, d = S[s].d;
        S[s].c1 = S[s].c2 = -1;
        if (L == 0)
            S[s].type = 0;
        else if (u[0] == 1)
        {
            S[s].type = 1;
            S[s].c1 = state_find(S, ns, u + 1, L - 1, d);
        }
        else if (u[L - 1] == 0)
        {
            S[s].type = 2;
            S[s].c1 = state_find(S, ns, u, L - 1, d);
        }
        else
        {
            S[s].type = 3;
            S[s].c1 = state_find(S, ns, u + 1, L - 1, d + 1);
            S[s].c2 = state_find(S, ns, u, L - 1, d);
        }
    }

    /* the roots, and the states needed at level 0 (sorted by length, so
       a backward pass propagates) */
    root = flint_malloc(sizeof(slong) * FLINT_MAX(nw, 1));
    for (i = 0; i < nw; i++)
    {
        root[i] = state_find(S, ns, words[i], lens[i], 0);
        S[root[i]].need0 = 1;
    }
    for (s = ns - 1; s >= 0; s--)
    {
        if (S[s].need0 && S[s].type != 0)
        {
            if (S[s].type == 2)
                flint_throw(FLINT_ERROR, "mzv_akhilesh: n = 0 in a local rule\n");
            S[S[s].c1].need0 = 1;
        }
    }

    ndel = 0;
    for (s = 0; s < ns; s++)
        ndel = FLINT_MAX(ndel, S[s].d + 1);

    mag_init(tb);

    /* N: the truncation bound ~ N^2 4^-(N - k) below 2^-(prec + 4) */
    N = prec / 2 + kmax + 8;
    while (1)
    {
        trunc_bound(tb, S, ns, N);
        if (mag_cmp_2exp_si(tb, -prec - 4) <= 0)
            break;
        N += 4 + N / 64;
    }

    /* fraction bits: the ulp counts grow like O(k N) */
    B = prec + 2 * FLINT_BIT_COUNT(N) + FLINT_BIT_COUNT(kmax) + 16;
    nl = (B + FLINT_BITS - 1) / FLINT_BITS;
    B = nl * FLINT_BITS;
    len = nl + 1;          /* one integer limb */

    V = flint_calloc(ns * len, sizeof(ulong));
    W = flint_calloc(ns * len, sizeof(ulong));
    E = flint_calloc(ndel * (len + 1), sizeof(ulong));
    tmp = flint_malloc(sizeof(ulong) * (len + 1) * KDIV);
    for (i = 0; i < KDIV; i++)
        dbuf[i] = tmp + i * (len + 1);
    eV = flint_calloc(ns, sizeof(ulong));
    eW = flint_calloc(ns, sizeof(ulong));
    eE = flint_calloc(ndel, sizeof(ulong));

#define VAL(A, s) ((A) + (s) * len)

    /* scaled empty-word values Ehat(n, d) = 4^n (n+d)! n! / (2n+d)! at
       n = N, as fixed point lower bounds (len limbs, entry has len + 1 for
       the products) */
    for (i = 0; i < ndel; i++)
        ehat_init(E + i * (len + 1), eE + i, N, i, B, len);

    /* Ehat(n, d) grows by about 2^d from n = N down to 0 (from about
       sqrt(pi N) 2^-d to 1), and its error in ulps with it: when the count
       exceeds what the other values accumulate (O(N) ulps), the value is
       recomputed from the factorials */
    emax = 4 * (ulong) N + 64;

    for (n = N; n >= 0; n--)
    {
        /* Ehat(n, d) from Ehat(n + 1, d): times (2n+d+2)(2n+d+1) / (4 (n+d+1)(n+1)) */
        if (n < N)
        {
            for (i = 0; i < ndel; i++)
            {
                nn_ptr e = E + i * (len + 1);
                ulong num = (ulong) (2 * n + i + 2) * (ulong) (2 * n + i + 1);
                ulong den = 4 * (ulong) (n + i + 1) * (ulong) (n + 1);
                ulong hi = mpn_mul_1(e, e, len, num);
                e[len] = hi;
                mpn_divrem_1(e, 0, e, len + 1, den);
                /* error: ceil(err num / den) + 1, exactly */
                {
                    ulong pr[2], q[2];
                    umul_ppmm(pr[1], pr[0], eE[i], num);
                    mpn_divrem_1(q, 0, pr, 2, den);
                    if (q[1] != 0 || q[0] > emax)
                        ehat_init(e, eE + i, n, i, B, len);
                    else
                        eE[i] = q[0] + 2;
                }
            }
        }

        /* tapering: every value at level n >= kmax is below
           4^-(n-kmax) (2n+1) (the bound above for nonempty subwords, and
           1/binom(2n,n) <= (2n+1)/4^n for the empty word), so only the
           low la limbs can be nonzero; la does not decrease with n
           decreasing, so the higher limbs of reused slots stay zero */
        if (n >= kmax)
        {
            slong z = 2 * (n - kmax) - FLINT_BIT_COUNT(2 * n + 1) - 2;
            la = (B - z + FLINT_BITS - 1) / FLINT_BITS + 1;
            la = FLINT_MAX(la, 1);
            la = FLINT_MIN(la, len);
        }
        else
            la = len;

        /* The divisions dominate (mpn_divrem_1 is latency bound): the two
           of a diagonal step are merged into one by (m+1)(n+1), the
           dividends are built with mul_1 / addmul_1 with the normalizing
           shift in the multipliers, and up to KDIV independent divisions
           (states of the same length, whose children are shorter) run
           interleaved with precomputed inverses. */
        nj = 0;
        for (s = 0; s <= ns; s++)
        {
            slong m;
            nn_ptr v;

            /* flush the pending divisions at a change of length or when full */
            if (nj > 0 && (s == ns || nj == KDIV || S[s].L != S[jobs[0]].L))
            {
                ulong R[KDIV], q;
                slong t;
                for (t = 0; t < KDIV; t++)
                    R[t] = 0;
                for (i = la; i >= 0; i--)
                    for (t = 0; t < nj; t++)
                    {
                        udiv_qrnnd_preinv(q, R[t], R[t], dbuf[t][i], dd[t], dinv[t]);
                        if (i < la)
                            VAL(V, jobs[t])[i] = q;
                    }
                for (t = 0; t < nj; t++)
                    if (S[jobs[t]].type == 3)
                        mpn_add_n(VAL(V, jobs[t]), VAL(V, jobs[t]), VAL(W, jobs[t]), la);
                nj = 0;
            }

            if (s == ns)
                break;

            if (n == 0 && !S[s].need0)
                continue;

            m = n + S[s].d;
            v = VAL(V, s);

            if (S[s].type == 0)
            {
                /* E(n, d) = Ehat(n, d) / 4^n */
                nn_srcptr e = E + S[s].d * (len + 1);
                slong sh = 2 * n, ls = sh / FLINT_BITS, bs = sh % FLINT_BITS;
                if (ls >= len)
                {
                    flint_mpn_zero(v, len);
                    eV[s] = 1;
                }
                else
                {
                    if (bs)
                        mpn_rshift(v, e + ls, len - ls, bs);
                    else
                        flint_mpn_copyi(v, e + ls, len - ls);
                    /* the limb above (e[len]) is zero for Ehat < 2^64 */
                    flint_mpn_zero(v + len - ls, ls);
                    eV[s] = (eE[S[s].d] >> FLINT_MIN(sh, FLINT_BITS - 1)) + 2;
                }
            }
            else
            {
                /* dividend in la + 1 limbs, normalized divisor */
                nn_ptr x = dbuf[nj];
                ulong d, sh;
                if (S[s].type == 1 || S[s].type == 2)
                {
                    d = (S[s].type == 1) ? (ulong) m : (ulong) n;
                    sh = flint_clz(d);
                    if (sh)
                        x[la] = mpn_lshift(x, VAL(V, S[s].c1), la, sh);
                    else
                    {
                        flint_mpn_copyi(x, VAL(V, S[s].c1), la);
                        x[la] = 0;
                    }
                    eV[s] = err_div(eV[S[s].c1], d);
                }
                else
                {
                    /* (V_n(c1) (n+1) + V_{n+1}(c2) (m+1)) / ((m+1)(n+1)),
                       then + V_{n+1}(s); one rounding: the error is at most
                       e1/(m+1) + e2/(n+1) + 1 ulps */
                    d = (ulong) (m + 1) * (ulong) (n + 1);
                    sh = flint_clz(d);
                    x[la] = mpn_mul_1(x, VAL(V, S[s].c1), la, (ulong) (n + 1) << sh);
                    x[la] += mpn_addmul_1(x, VAL(W, S[s].c2), la, (ulong) (m + 1) << sh);
                    eV[s] = eW[s] + (eV[S[s].c1] + m) / (m + 1) + (eW[S[s].c2] + n) / (n + 1) + 1;
                }
                dd[nj] = d << sh;
                dinv[nj] = n_preinvert_limb_prenorm(dd[nj]);
                jobs[nj] = s;
                nj++;
            }
        }

        {
            nn_ptr t = V; V = W; W = t;
        }
        {
            ulong * t = eV; eV = eW; eW = t;
        }
    }

    /* the results: lower bounds, plus the ulp errors and the truncation */
    for (i = 0; i < nw; i++)
    {
        nn_srcptr v = VAL(W, root[i]);
        slong sl = len;
        mag_t e;
        arf_t h;
        arb_ptr r = res + i;

        if (lens[i] == 0)
        {
            arb_one(r);
            continue;
        }

        mag_init(e);
        arf_init(h);
        while (sl > 0 && v[sl - 1] == 0)
            sl--;
        if (sl == 0)
            arf_zero(arb_midref(r));
        else
            arf_set_mpn(arb_midref(r), v, sl, 0);
        arf_mul_2exp_si(arb_midref(r), arb_midref(r), -B);
        /* value in [v, v + err ulps + trunc]: midpoint v + h, radius h */
        mag_set_ui(e, eW[root[i]]);
        mag_mul_2exp_si(e, e, -B);
        mag_add(e, e, tb);
        mag_mul_2exp_si(e, e, -1);
        arf_set_mag(h, e);
        arf_add(arb_midref(r), arb_midref(r), h, ARF_PREC_EXACT, ARF_RND_DOWN);
        mag_set(arb_radref(r), e);
        arb_set_round(r, r, prec);
        mag_clear(e);
        arf_clear(h);
    }

    flint_free(S);
    flint_free(root);
    flint_free(V);
    flint_free(W);
    flint_free(E);
    flint_free(tmp);
    flint_free(eV);
    flint_free(eW);
    flint_free(eE);
    mag_clear(tb);
}

void
_acb_mzeta_mzv_akhilesh(arb_t res, const unsigned char * w, slong k, slong prec)
{
    if (k == 0)
    {
        arb_one(res);
        return;
    }
    if (w[0] != 0 || w[k - 1] != 1)
    {
        /* not admissible */
        arb_indeterminate(res);
        return;
    }
    _acb_mzeta_mzv_akhilesh_vec(res, &w, &k, 1, prec);
}
