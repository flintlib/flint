/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "qfb.h"
#include "fmpq.h"
#include "arb.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "ecpp.h"

extern int ecpp_verbose;
#ifdef ECPP_PROFILE
#include <time.h>
static double _tnow(void)
{
    struct timespec t;
    timespec_get(&t, TIME_UTC);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}
#endif
#ifdef ECPP_PROFILE
/* accumulated over a proof: class group series, values (arb), decomposition (arb), descent (mod n) */
double ecpp_tower_prof[4];
slong ecpp_tower_count, ecpp_tower_hsum, ecpp_tower_retries;
#endif

/*
    A root of the Hilbert class polynomial H_D modulo a prime n that splits
    completely in the Hilbert class field H, through the class field tower
    (Enge, Morain: fast decomposition of polynomials with known Galois
    group; the representation follows Enge's mpfrcx and CM libraries).

    The Galois group of H/K is the class group Cl. A composition series
    1 = S_0 < S_1 < ... < S_L = Cl with prime quotients p_i = |S_i / S_{i-1}|
    gives a tower K = L_L < ... < L_1 < L_0 = H of fields L_i = fixed field
    of S_i, with [L_{i-1} : L_i] = p_i. The j-values are ordered so that
    consecutive blocks of size |S_i| are cosets of S_i. From the top: the
    relative minimal polynomial U over L_1 of j has degree p_1 and
    coefficients in L_1; one of them, y, generates L_1 over K, with minimal
    polynomial V of degree h / p_1 and integer coefficients (V is fixed by
    complex conjugation). The other coefficients c of U are written in
    Hecke form c V'(y) = W_c(y) with W_c the Lagrange numerator
    sum_k c_k prod_{l != k} (X - y_l), which has integer coefficients as
    well: no denominators appear. Then the conjugates of y take the role of
    the roots at the next level down, until the bottom polynomial over K of
    degree p_L. Modulo n every polynomial splits completely, so a root of
    the bottom polynomial gives an embedding of L_{L-1}, the Hecke forms
    evaluated at it give the next relative polynomial, and so on up to a
    root of U, which is a root of H_D. Every root finding is in a prime
    degree p_i instead of h.

    Returns 1 and sets j on success, 0 on failure (a generator could not
    be found or the coefficients not identified; the caller falls back).
*/

/* reduced form index lookup: forms sorted by (a, b) into a table of
   (a, b, index), binary search */
typedef struct { slong a, b, idx; } _form_key;

static int
_form_key_cmp(const void * x, const void * y)
{
    const _form_key * p = (const _form_key *) x, * q = (const _form_key *) y;
    if (p->a != q->a) return (p->a < q->a) ? -1 : 1;
    if (p->b != q->b) return (p->b < q->b) ? -1 : 1;
    return 0;
}

static _form_key *
_form_table(const qfb * forms, slong h)
{
    _form_key * t = flint_malloc(h * sizeof(_form_key));
    slong i;
    for (i = 0; i < h; i++)
    {
        t[i].a = fmpz_get_si(forms[i].a);
        t[i].b = fmpz_get_si(forms[i].b);
        t[i].idx = i;
    }
    qsort(t, h, sizeof(_form_key), _form_key_cmp);
    return t;
}

static slong
_form_index(const _form_key * t, slong h, const qfb_t f)
{
    _form_key key, * r;
    key.a = fmpz_get_si(f->a);
    key.b = fmpz_get_si(f->b);
    r = bsearch(&key, t, h, sizeof(_form_key), _form_key_cmp);
    return (r == NULL) ? -1 : r->idx;
}

/*
    Orders the classes as described: perm[k] is the index of the k-th
    class, degs[0..levels-1] the relative degrees from the bottom
    (degs[levels-1] = |S_1|). Returns levels, or 0 on failure.
*/
static slong
_classgroup_series(slong * perm, slong * degs, const qfb * forms, slong h, slong D)
{
    fmpz_t Df, L;
    qfb_t x, y, t;
    slong * cur;        /* elements of the current subgroup, in order */
    slong size = 1, levels = 0, i, k, idx;
    char * in;
    int ok = 1;

    _form_key * table = _form_table(forms, h);

    fmpz_init_set_si(Df, D);
    fmpz_init(L);
    fmpz_set_si(L, -D);
    fmpz_root(L, L, 4);
    qfb_init(x); qfb_init(y); qfb_init(t);

    /* the principal form */
    for (i = 0; i < h; i++)
        if (fmpz_is_one(forms[i].a))
            break;
    if (i == h)
    {
        ok = 0;
        goto cleanup;
    }
    cur = perm;
    cur[0] = i;
    in = flint_calloc(h, sizeof(char));
    in[i] = 1;

    while (size < h && ok)
    {
        slong ord, p, cnt;

        /*
            An element x outside the subgroup of prime order p modulo it,
            with p the largest prime factor of the order of x modulo the
            subgroup, among the first few candidates: the large primes are
            then the top levels of the tower (where the Kummer descent
            applies) and 2 and 3 the bottom (radicals).
        */
        {
            slong cand, ncand = 0, best_p = 0, best_idx = -1, best_ord = 0;
            for (cand = 0; cand < h && ncand < 16 && ok; cand++)
            {
                slong q, rest, pc;
                if (in[cand])
                    continue;
                ncand++;
                qfb_set(x, (qfb *) forms + cand);
                qfb_set(t, x);
                ord = 1;
                while (_form_index(table, h, t) < 0 || !in[_form_index(table, h, t)])
                {
                    qfb_nucomp(t, t, x, Df, L);
                    qfb_reduce(t, t, Df);
                    ord++;
                    if (ord > h)
                    {
                        ok = 0;
                        break;
                    }
                }
                if (!ok)
                    break;
                for (q = 2, rest = ord, pc = 1; rest > 1; q++)
                    if (rest % q == 0)
                    {
                        pc = q;
                        while (rest % q == 0)
                            rest /= q;
                    }
                if (pc > best_p)
                {
                    best_p = pc;
                    best_idx = cand;
                    best_ord = ord;
                }
            }
            if (!ok || best_idx < 0)
            {
                ok = 0;
                break;
            }
            idx = best_idx;
            p = best_p;
            ord = best_ord;
            qfb_set(x, (qfb *) forms + idx);
        }
        if (ord != p)
        {
            /* x^(ord/p) has order p modulo the subgroup */
            qfb_pow_ui(x, x, Df, ord / p);
            qfb_reduce(x, x, Df);
        }

        /* new subgroup: S, xS, x^2 S, ..., x^(p-1) S */
        for (k = 1; k < p; k++)
        {
            for (i = 0; i < size; i++)
            {
                slong pos;
                qfb_nucomp(t, x, (qfb *) forms + cur[(k - 1) * size + i], Df, L);
                qfb_reduce(t, t, Df);
                pos = _form_index(table, h, t);
                if (pos < 0 || in[pos])
                {
                    ok = 0;
                    break;
                }
                cur[k * size + i] = pos;
                in[pos] = 1;
            }
            if (!ok)
                break;
        }
        cnt = size * p;
        degs[levels++] = p;
        size = cnt;
    }
    flint_free(in);

    /* degs were filled from the top (S_1 first): reverse to bottom-first */
    for (i = 0; i < levels / 2; i++)
    {
        slong tmp = degs[i];
        degs[i] = degs[levels - 1 - i];
        degs[levels - 1 - i] = tmp;
    }

cleanup:
    flint_free(table);
    fmpz_clear(Df);
    fmpz_clear(L);
    qfb_clear(x); qfb_clear(y); qfb_clear(t);
    return ok ? levels : 0;
}

/*
    Subproduct tree for P = prod_k (X - y_k) and the Lagrange numerators
    W_j = sum_k c_{j,k} prod_{l != k} (X - y_l), j < m: for a split
    S = S1 u S2, P = P1 P2 and W = W1 P2 + W2 P1. Products only, so the
    error bounds stay tight (synthetic division by X - y_k loses a number
    of bits proportional to the height at each step).
*/
static void
_hecke_tree(acb_poly_t P, acb_poly_struct * W, slong m, acb_srcptr ys,
                        const acb_struct * const * cs, slong lo, slong hi, slong prec)
{
    slong j;

    if (hi - lo == 1)
    {
        acb_t t;
        acb_init(t);
        acb_poly_zero(P);
        acb_poly_set_coeff_si(P, 1, 1);
        acb_neg(t, ys + lo);
        acb_poly_set_coeff_acb(P, 0, t);
        for (j = 0; j < m; j++)
        {
            acb_poly_zero(W + j);
            acb_poly_set_coeff_acb(W + j, 0, cs[j] + lo);
        }
        acb_clear(t);
    }
    else
    {
        slong mid = lo + (hi - lo) / 2;
        acb_poly_t P1, P2, T;
        acb_poly_struct * W1 = flint_malloc(m * sizeof(acb_poly_struct));
        acb_poly_struct * W2 = flint_malloc(m * sizeof(acb_poly_struct));
        acb_poly_init(P1); acb_poly_init(P2); acb_poly_init(T);
        for (j = 0; j < m; j++)
        {
            acb_poly_init(W1 + j);
            acb_poly_init(W2 + j);
        }
        _hecke_tree(P1, W1, m, ys, cs, lo, mid, prec);
        _hecke_tree(P2, W2, m, ys, cs, mid, hi, prec);
        acb_poly_mul(P, P1, P2, prec);
        for (j = 0; j < m; j++)
        {
            acb_poly_mul(W + j, W1 + j, P2, prec);
            acb_poly_mul(T, W2 + j, P1, prec);
            acb_poly_add(W + j, W + j, T, prec);
        }
        for (j = 0; j < m; j++)
        {
            acb_poly_clear(W1 + j);
            acb_poly_clear(W2 + j);
        }
        flint_free(W1); flint_free(W2);
        acb_poly_clear(P1); acb_poly_clear(P2); acb_poly_clear(T);
    }
}

/*
    Kummer data for a level of prime degree p >= 5. The extension
    L_{i-1} = L_i(x) is cyclic, generated by sigma (the class group element
    of order p), with conjugates x_k = sigma^k(x). For a primitive p-th
    root of unity zeta the Lagrange resolvent theta_zeta = sum_k zeta^k x_k
    satisfies sigma(theta_zeta) = zeta^{-1} theta_zeta, so
    A = theta_zeta^p and, for 2 <= s <= p-1, beta_s = theta_{zeta^s}
    theta_zeta^{p-s} are algebraic integers in L_i(zeta) = L_i + L_i zeta
    + ... + L_i zeta^{p-2}. Their coordinates (times p^{p-1}, to clear the
    denominators of the power basis of zeta) are algebraic integers of
    L_i, hence have Hecke forms like the coefficients of the relative
    polynomial; they are recovered from the values at the p-1 primitive
    roots by solving a small linear system. Modulo n with p | n-1 the
    level is then descended with one p-th root: theta = A^{1/p},
    theta_s = beta_s theta^s / A, and x = (tr + sum_s theta_s)/p where
    tr = theta_1 is the trace, instead of a root of a degree-p polynomial.
    The polynomial is kept for the other n.

    Layout of W[lev]: [0, p) coefficients of the relative polynomial,
    [p] V', then p-1 coordinates of A, then (p-2)(p-1) coordinates of
    b_2, ..., b_{p-1}.
*/
static slong
_level_size(slong p)
{
    if (p < 5)
        return p + 1;
    return p + 1 + 2 * ((p - 1) + (p - 2) * (p - 1));
}

/* the Kummer numerators have coefficients in K = Q(sqrt D), stored as two
   integer polynomials: coefficient = (u + v sqrt D) / 2 */
#define KUMMER_A(p, j) ((p) + 1 + 2 * (j))
#define KUMMER_B(p, s, j) ((p) + 1 + 2 * (((p) - 1) + ((s) - 2) * ((p) - 1) + (j)))

/*
    Given the p roots x_0..x_{p-1} of a block, the values, at the coset,
    of the p-1 coordinates of A and of the b_s (s = 2..p-1): out has
    (p-1) + (p-2)(p-1) entries. Minv is the inverse of the (p-1) x (p-1)
    matrix [zeta^{s j}], s = 1..p-1, j = 0..p-2.
*/
static void
_kummer_values(acb_ptr out, acb_srcptr x, slong p, const acb_mat_t Minv, slong prec)
{
    acb_ptr theta = _acb_vec_init(p), vals = _acb_vec_init(p - 1), zeta = _acb_vec_init(p);
    acb_t t;
    slong u, k, s, j;

    acb_init(t);
    /* zeta[u] = exp(2 pi i u / p) */
    _acb_vec_unit_roots(zeta, p, p, prec);
    /* theta_u = sum_k zeta^{u k} x_k */
    for (u = 0; u < p; u++)
    {
        acb_zero(theta + u);
        for (k = 0; k < p; k++)
            acb_addmul(theta + u, zeta + (u * k) % p, x + k, prec);
    }
    /* A at the primitive roots: theta_u^p, u = 1..p-1 -> coordinates */
    for (u = 1; u < p; u++)
        acb_pow_ui(vals + u - 1, theta + u, p, prec);
    for (j = 0; j < p - 1; j++)
    {
        acb_zero(out + j);
        for (u = 1; u < p; u++)
            acb_addmul(out + j, acb_mat_entry(Minv, j, u - 1), vals + u - 1, prec);
    }
    /* scaled by p^{p-1}: the coordinates in the power basis of zeta have
       denominators dividing a power of p */
    acb_set_si(t, p);
    acb_pow_ui(t, t, p - 1, prec);
    for (j = 0; j < p - 1; j++)
        acb_mul(out + j, out + j, t, prec);
    /* beta_s = theta_{zeta^s} theta_zeta^{p-s} (an algebraic integer, unlike
       theta_{zeta^s} / theta_zeta^s) at the primitive root zeta^u */
    for (s = 2; s < p; s++)
    {
        for (u = 1; u < p; u++)
        {
            acb_pow_ui(vals + u - 1, theta + u, p - s, prec);
            acb_mul(vals + u - 1, vals + u - 1, theta + ((u * s) % p), prec);
        }
        for (j = 0; j < p - 1; j++)
        {
            acb_ptr o = out + (p - 1) + (s - 2) * (p - 1) + j;
            acb_zero(o);
            for (u = 1; u < p; u++)
                acb_addmul(o, acb_mat_entry(Minv, j, u - 1), vals + u - 1, prec);
            acb_mul(o, o, t, prec);
        }
    }
    acb_clear(t);
    _acb_vec_clear(theta, p);
    _acb_vec_clear(vals, p - 1);
    _acb_vec_clear(zeta, p);
}

/* w^p = z modulo n = 1 mod p (Adleman-Manders-Miller), z a p-th power */
static int
_rootmod_p(fmpz_t w, const fmpz_t z, slong p, flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_t t, e, y, g, G, tmp, pw, m, Gi, cur, pw1;
    slong sv = 0, k, j, tries, digit;
    int result = 0;

    if (fmpz_is_zero(z))
    {
        fmpz_zero(w);
        return 1;
    }
    fmpz_init(t); fmpz_init(e); fmpz_init(y); fmpz_init(g); fmpz_init(G);
    fmpz_init(tmp); fmpz_init(pw); fmpz_init(m); fmpz_init(Gi); fmpz_init(cur); fmpz_init(pw1);

    fmpz_sub_ui(t, n, 1);
    while (fmpz_divisible_si(t, p))
    {
        fmpz_divexact_si(t, t, p);
        sv++;
    }
    if (sv == 0)
        goto cleanup;
    fmpz_set_ui(e, p);
    fmpz_invmod(e, e, t);
    fmpz_powm(w, z, e, n);
    fmpz_powm(y, z, t, n);
    if (fmpz_is_one(y))
    {
        result = 1;
        goto cleanup;
    }
    fmpz_mul_ui(tmp, e, p);
    fmpz_sub_ui(tmp, tmp, 1);
    fmpz_divexact(tmp, tmp, t);
    fmpz_powm(y, y, tmp, n);            /* y^k in the p-Sylow subgroup, a p-th power there */

    for (tries = 0; tries < 100; tries++)
    {
        fmpz_randm(g, state, n);
        fmpz_powm(G, g, t, n);
        fmpz_set(pw, G);
        for (j = 0; j < sv - 1; j++)
            fmpz_mod_pow_ui(pw, pw, p, ctx);
        if (!fmpz_is_one(pw) && !fmpz_is_zero(pw))
            break;
    }
    if (tries == 100)
        goto cleanup;
    /* pw1 = G^{p^{sv-1}}: an element of order p */
    fmpz_set(pw1, pw);
    fmpz_mod_inv(Gi, G, ctx);
    fmpz_zero(m);
    for (k = 0; k < sv; k++)
    {
        fmpz_powm(cur, Gi, m, n);
        fmpz_mod_mul(cur, cur, y, ctx);
        for (j = 0; j < sv - 1 - k; j++)
            fmpz_mod_pow_ui(cur, cur, p, ctx);
        /* cur = pw1^digit */
        fmpz_one(tmp);
        for (digit = 0; digit < p; digit++)
        {
            if (fmpz_equal(cur, tmp))
                break;
            fmpz_mod_mul(tmp, tmp, pw1, ctx);
        }
        if (digit == p)
            break;
        fmpz_set_ui(tmp, 1);
        for (j = 0; j < k; j++)
            fmpz_mul_ui(tmp, tmp, p);
        fmpz_addmul_ui(m, tmp, digit);
    }
    if (k == sv && fmpz_divisible_si(m, p))
    {
        fmpz_divexact_si(m, m, p);
        fmpz_powm(cur, Gi, m, n);
        fmpz_mod_mul(w, w, cur, ctx);
        fmpz_mod_pow_ui(cur, w, p, ctx);
        result = fmpz_equal(cur, z);
    }

cleanup:
    fmpz_clear(t); fmpz_clear(e); fmpz_clear(y); fmpz_clear(g); fmpz_clear(G);
    fmpz_clear(tmp); fmpz_clear(pw); fmpz_clear(m); fmpz_clear(Gi); fmpz_clear(cur); fmpz_clear(pw1);
    return result;
}

/*
    Descent through a Kummer level: r the current root (value of the
    generator y of L_i), deninv = 1 / V'(r); x is set to a root of the
    relative polynomial U (given for verification). Returns 1 on success.
*/
static void _eval_K(fmpz_t c, const fmpz_poly_struct * W, slong idx, const fmpz_t r,
                        const fmpz_t sqrtD, fmpz_mod_poly_t Wm, const fmpz_mod_ctx_t ctx);

/*
    Arithmetic in F_{n^2} = F_n[s] / (s^2 - c), c a non-residue, for the
    Kummer descent when zeta_p lies in F_{n^2} but not in F_n
    (n = -1 mod p): elements (a, b) = a + b s.
*/
static void
_f2_mul(fmpz_t ra, fmpz_t rb, const fmpz_t a1, const fmpz_t b1,
        const fmpz_t a2, const fmpz_t b2, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    fmpz_t t, u;
    fmpz_init(t); fmpz_init(u);
    fmpz_mod_mul(t, b1, b2, ctx);
    fmpz_mod_mul(t, t, c, ctx);
    fmpz_mod_addmul(t, t, a1, a2, ctx);         /* a1 a2 + c b1 b2 */
    fmpz_mod_mul(u, a1, b2, ctx);
    fmpz_mod_addmul(u, u, a2, b1, ctx);         /* a1 b2 + a2 b1 */
    fmpz_swap(ra, t);
    fmpz_swap(rb, u);
    fmpz_clear(t); fmpz_clear(u);
}

static void
_f2_pow(fmpz_t ra, fmpz_t rb, const fmpz_t a, const fmpz_t b, const fmpz_t e,
                                        const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    fmpz_t xa, xb;
    slong i;
    fmpz_init(xa); fmpz_init(xb);
    fmpz_one(xa);
    fmpz_zero(xb);
    for (i = fmpz_bits(e) - 1; i >= 0; i--)
    {
        _f2_mul(xa, xb, xa, xb, xa, xb, c, ctx);
        if (fmpz_tstbit(e, i))
            _f2_mul(xa, xb, xa, xb, a, b, c, ctx);
    }
    fmpz_swap(ra, xa);
    fmpz_swap(rb, xb);
    fmpz_clear(xa); fmpz_clear(xb);
}

static void
_f2_inv(fmpz_t ra, fmpz_t rb, const fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                        const fmpz_mod_ctx_t ctx)
{
    fmpz_t nrm, t;
    fmpz_init(nrm); fmpz_init(t);
    fmpz_mod_mul(nrm, b, b, ctx);
    fmpz_mod_mul(nrm, nrm, c, ctx);
    fmpz_mod_mul(t, a, a, ctx);
    fmpz_mod_sub(nrm, t, nrm, ctx);             /* a^2 - c b^2 */
    fmpz_mod_inv(nrm, nrm, ctx);
    fmpz_mod_mul(t, a, nrm, ctx);
    fmpz_mod_mul(rb, b, nrm, ctx);
    fmpz_mod_neg(rb, rb, ctx);
    fmpz_swap(ra, t);
    fmpz_clear(nrm); fmpz_clear(t);
}

/*
    The Kummer descent in F_{n^2} for n = -1 mod p, p exactly dividing
    n + 1: a p-th root of a p-th power A is A^e with p e = 1 mod
    (n^2 - 1) / p, and the result x = (tr + sum theta_s) / p lies in F_n.
*/
static int
_kummer_descent_f2(fmpz_t x, const fmpz_poly_struct * W, slong p, const fmpz_t r,
        const fmpz_t deninv, const fmpz_t sqrtD, const fmpz_mod_poly_t U,
        flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_mod_poly_t Wm;
    fmpz_t c, e, t, za, zb, pa, pb, Aa, Ab, tha, thb, ta, tb, ba, bb, ia, ib, xa, xb, scale;
    slong j, tries;
    int ok = 0;

    fmpz_mod_poly_init(Wm, ctx);
    fmpz_init(c); fmpz_init(e); fmpz_init(t); fmpz_init(za); fmpz_init(zb);
    fmpz_init(pa); fmpz_init(pb); fmpz_init(Aa); fmpz_init(Ab); fmpz_init(tha); fmpz_init(thb);
    fmpz_init(ta); fmpz_init(tb); fmpz_init(ba); fmpz_init(bb); fmpz_init(ia); fmpz_init(ib);
    fmpz_init(xa); fmpz_init(xb); fmpz_init(scale);

    /* p must divide n + 1 exactly once */
    fmpz_add_ui(t, n, 1);
    if (!fmpz_divisible_si(t, p))
        goto cleanup;
    fmpz_divexact_si(t, t, p);
    if (fmpz_divisible_si(t, p))
        goto cleanup;

    /* a non-residue c and a primitive p-th root of unity zeta = g^((n^2-1)/p) */
    for (tries = 0; tries < 100; tries++)
    {
        fmpz_randm(c, state, n);
        if (fmpz_jacobi(c, n) == -1)
            break;
    }
    if (tries == 100)
        goto cleanup;
    fmpz_mul(e, n, n);
    fmpz_sub_ui(e, e, 1);
    fmpz_divexact_si(e, e, p);
    for (tries = 0; tries < 100; tries++)
    {
        fmpz_randm(ta, state, n);
        fmpz_randm(tb, state, n);
        _f2_pow(za, zb, ta, tb, e, c, ctx);
        if (!(fmpz_is_one(za) && fmpz_is_zero(zb)) && !(fmpz_is_zero(za) && fmpz_is_zero(zb)))
            break;
    }
    if (tries == 100)
        goto cleanup;

    /* 1 / p^{p-1} */
    fmpz_set_ui(scale, p);
    fmpz_mod_pow_ui(scale, scale, p - 1, ctx);
    fmpz_mod_inv(scale, scale, ctx);

    /* A = sum_j A_j zeta^j / p^{p-1} */
    fmpz_zero(Aa); fmpz_zero(Ab);
    fmpz_one(pa); fmpz_zero(pb);
    for (j = 0; j < p - 1; j++)
    {
        _eval_K(t, W, KUMMER_A(p, j), r, sqrtD, Wm, ctx);
        fmpz_mod_mul(t, t, deninv, ctx);
        fmpz_mod_mul(t, t, scale, ctx);
        fmpz_mod_addmul(Aa, Aa, t, pa, ctx);
        fmpz_mod_addmul(Ab, Ab, t, pb, ctx);
        _f2_mul(pa, pb, pa, pb, za, zb, c, ctx);
    }
    if (fmpz_is_zero(Aa) && fmpz_is_zero(Ab))
        goto cleanup;
    /* theta = A^e2 with p e2 = 1 mod (n^2-1)/p */
    fmpz_set_ui(t, p);
    fmpz_invmod(t, t, e);
    _f2_pow(tha, thb, Aa, Ab, t, c, ctx);
    /* check theta^p = A */
    fmpz_set_ui(t, p);
    _f2_pow(ta, tb, tha, thb, t, c, ctx);
    if (!fmpz_equal(ta, Aa) || !fmpz_equal(tb, Ab))
        goto cleanup;
    _f2_inv(ia, ib, Aa, Ab, c, ctx);            /* 1 / A */

    /* x = (tr + theta + sum_{s>=2} beta_s theta^s / A) / p */
    fmpz_mod_poly_set_fmpz_poly(Wm, W + (p - 1), ctx);
    fmpz_mod_poly_evaluate_fmpz(xa, Wm, r, ctx);
    fmpz_mod_mul(xa, xa, deninv, ctx);
    fmpz_mod_neg(xa, xa, ctx);
    fmpz_zero(xb);
    fmpz_mod_add(xa, xa, tha, ctx);
    fmpz_mod_add(xb, xb, thb, ctx);
    fmpz_set(ta, tha); fmpz_set(tb, thb);       /* theta^s */
    {
        slong sidx;
        for (sidx = 2; sidx < p; sidx++)
        {
            _f2_mul(ta, tb, ta, tb, tha, thb, c, ctx);
            fmpz_zero(ba); fmpz_zero(bb);
            fmpz_one(pa); fmpz_zero(pb);
            for (j = 0; j < p - 1; j++)
            {
                _eval_K(t, W, KUMMER_B(p, sidx, j), r, sqrtD, Wm, ctx);
                fmpz_mod_mul(t, t, deninv, ctx);
                fmpz_mod_mul(t, t, scale, ctx);
                fmpz_mod_addmul(ba, ba, t, pa, ctx);
                fmpz_mod_addmul(bb, bb, t, pb, ctx);
                _f2_mul(pa, pb, pa, pb, za, zb, c, ctx);
            }
            _f2_mul(ba, bb, ba, bb, ta, tb, c, ctx);
            _f2_mul(ba, bb, ba, bb, ia, ib, c, ctx);
            fmpz_mod_add(xa, xa, ba, ctx);
            fmpz_mod_add(xb, xb, bb, ctx);
        }
    }
    if (!fmpz_is_zero(xb))
        goto cleanup;                           /* not in F_n: wrong data */
    fmpz_set_ui(t, p);
    fmpz_mod_inv(t, t, ctx);
    fmpz_mod_mul(x, xa, t, ctx);
    fmpz_mod_poly_evaluate_fmpz(t, U, x, ctx);
    ok = fmpz_is_zero(t);

cleanup:
    fmpz_mod_poly_clear(Wm, ctx);
    fmpz_clear(c); fmpz_clear(e); fmpz_clear(t); fmpz_clear(za); fmpz_clear(zb);
    fmpz_clear(pa); fmpz_clear(pb); fmpz_clear(Aa); fmpz_clear(Ab); fmpz_clear(tha); fmpz_clear(thb);
    fmpz_clear(ta); fmpz_clear(tb); fmpz_clear(ba); fmpz_clear(bb); fmpz_clear(ia); fmpz_clear(ib);
    fmpz_clear(xa); fmpz_clear(xb); fmpz_clear(scale);
    return ok;
}

/* (u(r) + v(r) sqrtD) / 2 for a K-valued Hecke numerator */
static void
_eval_K(fmpz_t c, const fmpz_poly_struct * W, slong idx, const fmpz_t r,
                        const fmpz_t sqrtD, fmpz_mod_poly_t Wm, const fmpz_mod_ctx_t ctx)
{
    fmpz_t t, h;
    fmpz_init(t); fmpz_init(h);
    fmpz_mod_poly_set_fmpz_poly(Wm, W + idx, ctx);
    fmpz_mod_poly_evaluate_fmpz(c, Wm, r, ctx);
    fmpz_mod_poly_set_fmpz_poly(Wm, W + idx + 1, ctx);
    fmpz_mod_poly_evaluate_fmpz(t, Wm, r, ctx);
    fmpz_mod_mul(t, t, sqrtD, ctx);
    fmpz_mod_add(c, c, t, ctx);
    fmpz_set_ui(h, 2);
    fmpz_mod_inv(h, h, ctx);
    fmpz_mod_mul(c, c, h, ctx);
    fmpz_clear(t); fmpz_clear(h);
}

static int
_kummer_descent(fmpz_t x, const fmpz_poly_struct * W, slong p, const fmpz_t r,
        const fmpz_t deninv, const fmpz_t sqrtD, const fmpz_mod_poly_t U,
        flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    fmpz_mod_poly_t Wm;
    fmpz_t zeta, zpow, A, theta, ths, b, c, t, e, g;
    slong j, s, tries;
    int ok = 0;

    fmpz_mod_poly_init(Wm, ctx);
    fmpz_init(zeta); fmpz_init(zpow); fmpz_init(A); fmpz_init(theta); fmpz_init(ths);
    fmpz_init(b); fmpz_init(c); fmpz_init(t); fmpz_init(e); fmpz_init(g);

    /* a primitive p-th root of unity */
    fmpz_sub_ui(e, n, 1);
    fmpz_divexact_si(e, e, p);
    for (tries = 0; tries < 100; tries++)
    {
        fmpz_randm(g, state, n);
        fmpz_powm(zeta, g, e, n);
        if (!fmpz_is_one(zeta) && !fmpz_is_zero(zeta))
            break;
    }
    if (tries == 100)
        goto cleanup;

    /* A = sum_j A_j zeta^j / p^{p-1} */
    fmpz_zero(A);
    fmpz_one(zpow);
    for (j = 0; j < p - 1; j++)
    {
        _eval_K(c, W, KUMMER_A(p, j), r, sqrtD, Wm, ctx);
        fmpz_mod_mul(c, c, deninv, ctx);
        fmpz_mod_mul(c, c, zpow, ctx);
        fmpz_mod_add(A, A, c, ctx);
        fmpz_mod_mul(zpow, zpow, zeta, ctx);
    }
    fmpz_set_ui(t, p);
    fmpz_mod_pow_ui(t, t, p - 1, ctx);
    fmpz_mod_inv(t, t, ctx);            /* 1 / p^{p-1} */
    fmpz_mod_mul(A, A, t, ctx);
    if (fmpz_is_zero(A) || !_rootmod_p(theta, A, p, state, ctx))
        goto cleanup;
    fmpz_mod_inv(g, A, ctx);            /* 1 / A */

    /* x = (tr + theta + sum_{s >= 2} b_s theta^s) / p, tr = -c_{p-1} */
    fmpz_mod_poly_set_fmpz_poly(Wm, W + (p - 1), ctx);
    fmpz_mod_poly_evaluate_fmpz(x, Wm, r, ctx);
    fmpz_mod_mul(x, x, deninv, ctx);
    fmpz_mod_neg(x, x, ctx);
    fmpz_mod_add(x, x, theta, ctx);
    fmpz_set(ths, theta);
    for (s = 2; s < p; s++)
    {
        fmpz_mod_mul(ths, ths, theta, ctx);     /* theta^s */
        fmpz_zero(b);
        fmpz_one(zpow);
        for (j = 0; j < p - 1; j++)
        {
            _eval_K(c, W, KUMMER_B(p, s, j), r, sqrtD, Wm, ctx);
            fmpz_mod_mul(c, c, deninv, ctx);
            fmpz_mod_mul(c, c, zpow, ctx);
            fmpz_mod_add(b, b, c, ctx);
            fmpz_mod_mul(zpow, zpow, zeta, ctx);
        }
        fmpz_mod_mul(b, b, t, ctx);         /* / p^{p-1} */
        fmpz_mod_mul(b, b, ths, ctx);       /* beta_s theta^s */
        fmpz_mod_mul(b, b, g, ctx);         /* / A */
        fmpz_mod_add(x, x, b, ctx);
    }
    fmpz_set_ui(c, p);
    fmpz_mod_inv(c, c, ctx);
    fmpz_mod_mul(x, x, c, ctx);

    fmpz_mod_poly_evaluate_fmpz(t, U, x, ctx);
    ok = fmpz_is_zero(t);

cleanup:
    fmpz_mod_poly_clear(Wm, ctx);
    fmpz_clear(zeta); fmpz_clear(zpow); fmpz_clear(A); fmpz_clear(theta); fmpz_clear(ths);
    fmpz_clear(b); fmpz_clear(c); fmpz_clear(t); fmpz_clear(e); fmpz_clear(g);
    return ok;
}

/* rounds an acb polynomial with coefficients (u + v sqrt D) / 2 in O_K into
   the integer polynomials u, v; sqrtD = i sqrt|D|; 0 if not identified */
static int
_acb_poly_get_K_poly(fmpz_poly_t u, fmpz_poly_t v, const acb_poly_t p, const acb_t sqrtD, slong prec)
{
    slong i, len = acb_poly_length(p);
    fmpz_t c;
    arb_t t;
    int ok = 1;

    fmpz_init(c);
    arb_init(t);
    fmpz_poly_zero(u);
    fmpz_poly_zero(v);
    for (i = 0; i < len && ok; i++)
    {
        const acb_struct * z = acb_poly_get_coeff_ptr((acb_poly_struct *) p, i);
        arb_mul_2exp_si(t, acb_realref(z), 1);
        if (!arb_get_unique_fmpz(c, t))
            ok = 0;
        else
        {
            fmpz_poly_set_coeff_fmpz(u, i, c);
            arb_div(t, acb_imagref(z), acb_imagref(sqrtD), prec);
            arb_mul_2exp_si(t, t, 1);
            if (!arb_get_unique_fmpz(c, t))
                ok = 0;
            else
                fmpz_poly_set_coeff_fmpz(v, i, c);
        }
    }
    fmpz_clear(c);
    arb_clear(t);
    return ok;
}

/* rounds an acb polynomial with integer coefficients; 0 if not identified */
static int
_acb_poly_get_fmpz_poly(fmpz_poly_t r, const acb_poly_t p)
{
    slong i, len = acb_poly_length(p);
    fmpz_t c;
    int ok = 1;

    fmpz_init(c);
    fmpz_poly_zero(r);
    for (i = 0; i < len && ok; i++)
    {
        const acb_struct * z = acb_poly_get_coeff_ptr((acb_poly_struct *) p, i);
        if (!arb_contains_zero(acb_imagref(z)) || !arb_get_unique_fmpz(c, acb_realref(z)))
            ok = 0;
        else
            fmpz_poly_set_coeff_fmpz(r, i, c);
    }
    fmpz_clear(c);
    return ok;
}

/*
    Weber's function f as class invariant, for fundamental D = -4m
    (m = 1, 2, 3, 5, 6, 7 mod 8). The theorems on which values are class
    invariants are those of Weber, Yui and Zagier (Math. Comp. 66, 1997)
    and Schertz (J. Theor. Nombres Bordeaux 14, 2002, whose N-systems fix
    the form representatives); the normalisations below (the power of f,
    the factors of sqrt 2, the cube when 3 | D and the sign (2/a)) are the
    ones of Enge's CM library, used here as the reference for which
    normalisation yields a polynomial with integer coefficients. The value is taken at tau = (-b + sqrt(D)) / (2a)
    for a form of the class in a 48-system: gcd(a, 48) = 1 and b = 0 mod
    96. With f = zeta_48^{-1} eta((tau+1)/2) / eta(tau) and
    f1 = eta(tau/2) / eta(tau):
        m = 1 mod 8: g = f^2 / sqrt 2      m = 5 mod 8: g = f^4 / 2
        m = 3 mod 8: g = f                 m = 7 mod 8: g = f / sqrt 2
        m = 2, 6 mod 8: g = f1^2 / sqrt 2
    then g = g^3 if 3 | D, and g = -g if (2/a) = -1 (except m = 3, 5). The
    class polynomial of g has integer coefficients and height a factor
    72 / e smaller than that of H_D, e in {1, 2, 4} times 3 if 3 | D.
    Back from a root g:
        m = 1, 2, 6: F = (8 g^{6/k})^2,  m = 5: F = 64 g^{6/k},
        m = 3: F = (g^{6/k})^4,  m = 7: F = (8 g^{6/k})^4
    (k = 3 if 3 | D, else 1), and j = (F - 16)^3 / F for m odd,
    j = (F + 16)^3 / F for m even (F = f^24 resp. f1^24).
*/

/*
    The residue class computations on D below use |D| = -D (a positive
    number) and shifts and masks rather than the division and remainder of
    a negative slong: the 32-bit build with gcc 15.2 on Alpine took the
    Weber path for D = -4047 with the latter (D % 4 != 0 being false),
    which no sanitizer or other compiler reproduces.
*/
static int
_weber_ok(slong D)
{
    ulong m;
    if (D >= 0 || (((ulong) (-D)) & 3) != 0)
        return 0;
    m = (((ulong) (-D)) >> 2) & 7;
    return (m == 1 || m == 2 || m == 3 || m == 5 || m == 6 || m == 7);
}

/* the height factor 72 / e */
static double
_weber_factor(slong D)
{
    ulong m8 = (((ulong) (-D)) >> 2) & 7;
    slong e = (m8 == 5) ? 4 : (m8 == 3 || m8 == 7) ? 1 : 2;
    if (((ulong) (-D)) % 3 == 0)
        e *= 3;
    return 72.0 / e;
}

/*
    An equivalent form (a', b', c') with gcd(a', N) = 1 and b' = b0 mod 2N
    (an N-system entry in the sense of Schertz), for a primitive form
    (a, b, c) of discriminant D with b0 = D mod 2. The first coefficient
    of the form transformed by [[x, u], [y, v]] in SL_2(Z) is Q(x, y): a
    pair (x, y) with gcd(Q(x, y), N) = 1 exists among small coprime pairs
    since the form is primitive, and is completed to a matrix with the
    extended Euclidean algorithm; a translation by k then moves b to
    b + 2 a' k = b0 mod 2N, i.e. k = (b0 - b)/2 * a'^{-1} mod N.
*/
static void
_nsystem_entry(slong * a, slong * b, slong N, slong b0, slong D)
{
    slong c = (*b * *b - D) / (4 * *a);
    slong x, y, ap, bp, u, v, k, r, s;
    ulong g;

    /* Q(x, y) coprime to N, |x|, |y| small: try (1,0), (0,1), then a grid */
    ap = 0; x = 1; y = 0;
    for (r = 0; r <= 8 && n_gcd(FLINT_ABS(ap), N) != 1; r++)
    {
        for (s = -r; s <= r && (r == 0 || n_gcd(FLINT_ABS(ap), N) != 1); s++)
        {
            x = (r == 0) ? 1 : r;
            y = (r == 0) ? 0 : s;
            if (n_gcd(FLINT_ABS(x), FLINT_ABS(y)) != 1)
                continue;
            ap = *a * x * x + *b * x * y + c * y * y;
            if (n_gcd(FLINT_ABS(ap), N) == 1)
                break;
            x = s; y = r;   /* the pair (s, r) as well */
            if (n_gcd(FLINT_ABS(x), FLINT_ABS(y)) != 1)
                continue;
            ap = *a * x * x + *b * x * y + c * y * y;
        }
    }
    if (n_gcd(FLINT_ABS(ap), N) != 1)
        return;     /* not reached for primitive forms; leave the form alone */

    /* complete (x, y) to [[x, u], [y, v]] with x v - y u = 1 */
    {
        fmpz_t gg, uu, vv, xx, yy;
        fmpz_init(gg); fmpz_init(uu); fmpz_init(vv);
        fmpz_init_set_si(xx, x); fmpz_init_set_si(yy, y);
        fmpz_xgcd(gg, vv, uu, xx, yy);      /* vv x + uu y = 1 */
        v = fmpz_get_si(vv);
        u = -fmpz_get_si(uu);               /* x v - y u = x v + y (-u) */
        fmpz_clear(gg); fmpz_clear(uu); fmpz_clear(vv); fmpz_clear(xx); fmpz_clear(yy);
    }
    /* b' = 2 a x u + b (x v + y u) + 2 c y v */
    bp = 2 * *a * x * u + *b * (x * v + y * u) + 2 * c * y * v;

    /* translate: k = (b0 - b')/2 * a'^{-1} mod N */
    k = (b0 - bp) / 2;
    k %= N;
    if (k < 0)
        k += N;
    g = n_invmod(n_mod2_preinv(FLINT_ABS(ap) % N, N, n_preinvert_limb(N)), N);
    if (ap < 0)
        g = (N - g) % N;
    k = (slong) n_mulmod2_preinv((ulong) k, g, N, n_preinvert_limb(N));
    *a = ap;
    *b = bp + 2 * ap * k;
}

/* g at the class of the form (a, b) */
static void
_weber_value(acb_t g, slong a0, slong b0, slong D, slong prec)
{
    slong a = a0, b = b0, m = (slong) (((ulong) (-D)) >> 2);
    acb_t tau, t, e1, e2;
    arb_t s2;

    _nsystem_entry(&a, &b, 48, 0, D);

    acb_init(tau); acb_init(t); acb_init(e1); acb_init(e2);
    arb_init(s2);
    arb_set_si(acb_realref(tau), -b);
    arb_set_si(acb_imagref(tau), -D);
    arb_sqrt(acb_imagref(tau), acb_imagref(tau), prec);
    acb_div_si(tau, tau, 2 * a, prec);
    arb_sqrt_ui(s2, 2, prec);

    acb_modular_eta(e2, tau, prec);
    if (m % 8 == 2 || m % 8 == 6)
    {
        /* f1 = eta(tau/2) / eta(tau) */
        acb_mul_2exp_si(t, tau, -1);
        acb_modular_eta(e1, t, prec);
        acb_div(g, e1, e2, prec);
        acb_mul(g, g, g, prec);
        acb_mul_arb(g, g, s2, prec);
        acb_mul_2exp_si(g, g, -1);
    }
    else
    {
        /* f = zeta_48^{-1} eta((tau+1)/2) / eta(tau) */
        acb_add_ui(t, tau, 1, prec);
        acb_mul_2exp_si(t, t, -1);
        acb_modular_eta(e1, t, prec);
        acb_div(g, e1, e2, prec);
        {
            /* zeta_48^{-1} = exp(-2 pi i / 48) = cos(pi/24) - i sin(pi/24) */
            fmpq_t x;
            fmpq_init(x);
            fmpq_set_si(x, -1, 24);
            arb_sin_cos_pi_fmpq(acb_imagref(t), acb_realref(t), x, prec);
            fmpq_clear(x);
        }
        acb_mul(g, g, t, prec);
        if (m % 8 == 5)
        {
            acb_mul(g, g, g, prec);
            acb_mul(g, g, g, prec);
            acb_mul_2exp_si(g, g, -1);
        }
        else if (m % 8 == 3)
        {
            /* g = f */
        }
        else if (m % 8 == 7)
        {
            acb_mul_arb(g, g, s2, prec);
            acb_mul_2exp_si(g, g, -1);
        }
        else
        {
            acb_mul(g, g, g, prec);
            acb_mul_arb(g, g, s2, prec);
            acb_mul_2exp_si(g, g, -1);
        }
    }
    if (((ulong) (-D)) % 3 == 0)
    {
        acb_mul(t, g, g, prec);
        acb_mul(g, g, t, prec);
    }
    if (m % 8 != 3 && m % 8 != 5 && n_jacobi_unsigned(2, FLINT_ABS(a)) == -1)
        acb_neg(g, g);

    acb_clear(tau); acb_clear(t); acb_clear(e1); acb_clear(e2);
    arb_clear(s2);
}

/* j from a root g of the class polynomial of the Weber invariant, mod n */
static void
_weber_to_j(fmpz_t j, const fmpz_t g, slong D, const fmpz_mod_ctx_t ctx)
{
    slong m = (slong) (((ulong) (-D)) >> 2);
    fmpz_t F, t;
    fmpz_init(F); fmpz_init(t);
    fmpz_mod_pow_ui(F, g, ((((ulong) (-D)) % 3) == 0) ? 2 : 6, ctx);
    if (m % 8 == 5)
        fmpz_mod_mul_ui(F, F, 64, ctx);
    else if (m % 8 == 3)
        fmpz_mod_pow_ui(F, F, 4, ctx);
    else if (m % 8 == 7)
    {
        fmpz_mod_mul_ui(F, F, 8, ctx);
        fmpz_mod_pow_ui(F, F, 4, ctx);
    }
    else
    {
        fmpz_mod_mul_ui(F, F, 8, ctx);
        fmpz_mod_mul(F, F, F, ctx);
    }
    if (m % 2 == 1)
        fmpz_mod_sub_ui(t, F, 16, ctx);
    else
        fmpz_mod_add_ui(t, F, 16, ctx);
    fmpz_mod_pow_ui(t, t, 3, ctx);
    fmpz_mod_inv(F, F, ctx);
    fmpz_mod_mul(j, t, F, ctx);
    fmpz_clear(F); fmpz_clear(t);
}

int
ecpp_class_poly_tower(fmpz_t j, slong D, flint_rand_t state, const fmpz_mod_ctx_t ctx)
{
    int flags = 0;
    if (fmpz_bits(fmpz_mod_ctx_modulus(ctx)) >= ECPP_KUMMER_BITS)
        flags |= ECPP_TOWER_KUMMER;
    return _ecpp_class_poly_tower(j, D, flags, state, ctx);
}

int
_ecpp_class_poly_tower(fmpz_t j, slong D, int flags, flint_rand_t state,
                                                        const fmpz_mod_ctx_t ctx)
{
    qfb * forms;
    slong h, levels, i, k, lev, prec;
    slong * perm, * degs;
    fmpz_poly_struct ** W;      /* W[lev][0..d] : level lev's Hecke forms, W[0][0] bottom */
    int * kummer;               /* level has the Kummer data */
    acb_ptr roots;
    arb_t sqrtD;
    acb_t z;
    double lgh;
    int success = 0, kummer_retries = 0;
    int weber = _weber_ok(D) && !(flags & ECPP_TOWER_J);
    /* the invariant is only a class invariant for D = 0 mod 4 */
    if (weber && (((ulong) (-D)) & 3) != 0)
        weber = 0;
    /* the Kummer data (and its retries) only pay when a powering of a
       degree-p polynomial is expensive, i.e. for large n; the cost model
       in disc.c prices the descents accordingly */
    int want_kummer = (flags & ECPP_TOWER_KUMMER) != 0;

    h = qfb_reduced_forms(&forms, D);
    if (h <= 0)
        return 0;
    if (h == 1)
    {
        fmpz_poly_t H;
        fmpz_poly_init(H);
        acb_modular_hilbert_class_poly(H, D);
        fmpz_mod_set_fmpz(j, H->coeffs + 0, ctx);
        fmpz_mod_neg(j, j, ctx);
        fmpz_poly_clear(H);
        qfb_array_clear(&forms, h);
        return 1;
    }

    perm = flint_malloc(h * sizeof(slong));
    degs = flint_malloc((FLINT_BIT_COUNT(h) + 2) * sizeof(slong));
#ifdef ECPP_PROFILE
    {
        double t0 = _tnow();
        levels = _classgroup_series(perm, degs, forms, h, D);
        ecpp_tower_prof[0] += _tnow() - t0;
        ecpp_tower_count++; ecpp_tower_hsum += h;
    }
#else
    levels = _classgroup_series(perm, degs, forms, h, D);
#endif
    if (levels == 0)
        goto cleanup_forms;

    W = flint_malloc(levels * sizeof(fmpz_poly_struct *));
    kummer = flint_calloc(levels, sizeof(int));
    for (lev = 0; lev < levels; lev++)
    {
        slong nd = (lev == 0) ? 1 : _level_size(degs[lev]);
        W[lev] = flint_malloc(nd * sizeof(fmpz_poly_struct));
        for (i = 0; i < nd; i++)
            fmpz_poly_init(W[lev] + i);
        if (lev > 0)
            kummer[lev] = 0;
    }

    /*
        Precision. The Lagrange numerators of the Hecke forms are sums of
        terms of the size of V times the coefficients, with heavy
        cancellation, and this happens at every level: in practice about
        six times the height of H_D is needed (with j as the invariant;
        class invariants of smaller height would reduce this in proportion,
        which is what makes class numbers in the hundreds affordable in
        Enge's CM library).
    */
    lgh = 0.0;
    for (i = 0; i < h; i++)
        lgh += 1.0 / fmpz_get_d(forms[i].a);
    prec = 3.141593 * sqrt((double) -D) * lgh * 1.442696 * 1.5 / (weber ? _weber_factor(D) : 1.0) + 64 + 4 * h;

    roots = _acb_vec_init(h);
    arb_init(sqrtD);
    acb_init(z);

    for (;;)
    {
        slong num = h;      /* number of roots at the current level */
        int kummer_failed = 0;
        acb_ptr ys = _acb_vec_init(h);
        acb_poly_t U, V, T, lin;
        acb_poly_init(U); acb_poly_init(V); acb_poly_init(T); acb_poly_init(lin);

#ifdef ECPP_PROFILE
        double tj0 = _tnow();
#endif
        /* j-values in class group order */
        arb_set_si(sqrtD, -D);
        arb_sqrt(sqrtD, sqrtD, prec);
        for (i = 0; i < h; i++)
        {
            const qfb * f = forms + perm[i];
            if (weber)
                _weber_value(roots + i, fmpz_get_si(f->a), fmpz_get_si(f->b), D, prec);
            else
            {
                arb_set_fmpz(acb_realref(z), f->b);
                arb_neg(acb_realref(z), acb_realref(z));
                arb_set(acb_imagref(z), sqrtD);
                acb_div_fmpz(z, z, f->a, prec);
                acb_mul_2exp_si(z, z, -1);
                acb_modular_j(roots + i, z, prec);
            }
        }

        if (ecpp_verbose)
            flint_printf("ecpp: tower for D = %wd, h = %wd, %wd levels, precision %wd (%s)\n",
                    D, h, levels, prec, weber ? "Weber" : "j");
#ifdef ECPP_PROFILE
        ecpp_tower_prof[1] += _tnow() - tj0;
        tj0 = _tnow();
#endif
        success = 1;
        for (lev = levels - 1; lev >= 1 && success; lev--)
        {
            slong m = degs[lev], nn = num / m, jj, found = -1;
            acb_poly_struct * Us = flint_malloc(nn * sizeof(acb_poly_struct));

            for (k = 0; k < nn; k++)
            {
                acb_poly_init(Us + k);
                acb_poly_product_roots(Us + k, roots + m * k, m, prec);
            }
            /* a generator among the coefficients, the trace first */
            for (jj = m - 1; jj >= 0 && found < 0; jj--)
            {
                int distinct = 1;
                for (k = 0; k < nn && distinct; k++)
                    acb_poly_get_coeff_acb(ys + k, Us + k, jj);
                for (k = 1; k < nn && distinct; k++)
                {
                    /* equal up to rounding? compare with a margin */
                    acb_sub(z, ys + 0, ys + k, prec);
                    if (acb_contains_zero(z))
                        distinct = 0;
                }
                if (distinct)
                    found = jj;
            }
            if (found < 0)
            {
                success = 0;
            }
            else
            {
                /* V = prod (X - y_k) and the Hecke numerators of the
                   coefficients c_{jj,k} of the U_k, by a subproduct tree */
                /* Kummer vectors: their numerators are larger than the
                   polynomial ones, so only where the numerics are cheap */
                slong nk = (m >= 5 && (weber || h <= 64) && want_kummer) ? (m - 1) + (m - 2) * (m - 1) : 0;
                slong nv = m + nk;
                acb_ptr * cs = flint_malloc(nv * sizeof(acb_ptr));
                acb_poly_struct * Ws = flint_malloc(nv * sizeof(acb_poly_struct));
                for (jj = 0; jj < nv; jj++)
                {
                    cs[jj] = _acb_vec_init(nn);
                    acb_poly_init(Ws + jj);
                }
                for (jj = 0; jj < m; jj++)
                    for (k = 0; k < nn; k++)
                        acb_poly_get_coeff_acb(cs[jj] + k, Us + k, jj);
                if (nk > 0)
                {
                    acb_mat_t Mz, Minv;
                    acb_ptr vals = _acb_vec_init(nk);
                    acb_t t;
                    slong u;
                    acb_ptr zm = _acb_vec_init(m);
                    acb_init(t);
                    acb_mat_init(Mz, m - 1, m - 1);
                    acb_mat_init(Minv, m - 1, m - 1);
                    /* zm[k] = exp(2 pi i k / m) */
                    _acb_vec_unit_roots(zm, m, m, prec);
                    for (u = 1; u < m; u++)
                        for (jj = 0; jj < m - 1; jj++)
                            acb_set(acb_mat_entry(Mz, u - 1, jj), zm + (u * jj) % m);
                    _acb_vec_clear(zm, m);
                    acb_mat_inv(Minv, Mz, prec);
                    for (k = 0; k < nn; k++)
                    {
                        _kummer_values(vals, roots + m * k, m, Minv, prec);
                        for (jj = 0; jj < nk; jj++)
                            acb_set(cs[m + jj] + k, vals + jj);
                    }
                    acb_clear(t);
                    _acb_vec_clear(vals, nk);
                    acb_mat_clear(Mz);
                    acb_mat_clear(Minv);
                }
                _hecke_tree(V, Ws, nv, ys, (const acb_struct * const *) cs, 0, nn, prec);
                if (!_acb_poly_get_fmpz_poly(W[lev] + m, V))   /* temporarily V */
                    success = 0;
                for (jj = 0; jj < m && success; jj++)
                    if (!_acb_poly_get_fmpz_poly(W[lev] + jj, Ws + jj))
                        success = 0;
                /* the Kummer data may fail to be integral (e.g. zeta in L_i):
                   then only the polynomial is kept */
                kummer[lev] = (nk > 0);
                if (kummer[lev])
                {
                    acb_t sD;
                    acb_init(sD);
                    arb_zero(acb_realref(sD));
                    arb_set(acb_imagref(sD), sqrtD);
                    for (jj = 0; jj < nk && kummer[lev]; jj++)
                        if (!_acb_poly_get_K_poly(W[lev] + m + 1 + 2 * jj, W[lev] + m + 2 + 2 * jj,
                                                    Ws + m + jj, sD, prec))
                        {
                            kummer[lev] = 0;
                            kummer_failed = 1;
                            if (ecpp_verbose)
                                flint_printf("ecpp: Kummer data not identified at degree %wd (precision %wd)\n", m, prec);
                        }
                    acb_clear(sD);
                }
                for (jj = 0; jj < nv; jj++)
                {
                    _acb_vec_clear(cs[jj], nn);
                    acb_poly_clear(Ws + jj);
                }
                flint_free(cs);
                flint_free(Ws);
                if (success)
                {
                    /* W[lev][m] = V' */
                    fmpz_poly_derivative(W[lev] + m, W[lev] + m);
                    /* the conjugates of y are the roots of the next level */
                    for (k = 0; k < nn; k++)
                        acb_set(roots + k, ys + k);
                    num = nn;
                }
            }
            for (k = 0; k < nn; k++)
                acb_poly_clear(Us + k);
            flint_free(Us);
        }
        if (success)
        {
            acb_poly_product_roots(V, roots, num, prec);
            if (!_acb_poly_get_fmpz_poly(W[0] + 0, V))
                success = 0;
        }

        acb_poly_clear(U); acb_poly_clear(V); acb_poly_clear(T); acb_poly_clear(lin);
        _acb_vec_clear(ys, h);
#ifdef ECPP_PROFILE
        ecpp_tower_prof[2] += _tnow() - tj0;
        if (!success || kummer_failed) ecpp_tower_retries++;
#endif

        /* the Kummer numerators are larger than the polynomial ones: a
           couple of retries at higher precision before giving them up */
        if (success && (!kummer_failed || kummer_retries >= 1))
            break;
        if (success)
        {
            kummer_retries++;
            prec = 2 * prec + 64;
            continue;
        }
        if (prec > 40 * (3.141593 * sqrt((double) -D) * lgh * 1.442696 + 1000))
            break;
        prec = prec * 3 / 2 + 64;
    }

    /* descent modulo n */
#ifdef ECPP_PROFILE
    {
    double td0 = _tnow();
#endif
    if (success)
    {
        fmpz_mod_poly_t f;
        fmpz_t r, den, c, sqrtDn;
        fmpz_mod_poly_init(f, ctx);
        fmpz_init(r); fmpz_init(den); fmpz_init(c);
        fmpz_init(sqrtDn);
        /* sqrt(D) mod n, needed by the Kummer data (either sign works
           for one of the two conjugate data sets; both are tried) */
        {
            slong lev2;
            int need = 0;
            for (lev2 = 1; lev2 < levels; lev2++)
                need |= kummer[lev2];
            if (need)
            {
                fmpz_set_si(sqrtDn, D);
                fmpz_mod_set_fmpz(sqrtDn, sqrtDn, ctx);
                if (!fmpz_sqrtmod(sqrtDn, sqrtDn, fmpz_mod_ctx_modulus(ctx)))
                    for (lev2 = 1; lev2 < levels; lev2++)
                        kummer[lev2] = 0;
            }
        }

        fmpz_mod_poly_set_fmpz_poly(f, W[0] + 0, ctx);
        success = ecpp_poly_root(r, f, state, ctx);
        for (lev = 1; lev < levels && success; lev++)
        {
            slong m = degs[lev], jj;
            ulong nm1;
            fmpz_mod_poly_t Wm;
            fmpz_mod_poly_init(Wm, ctx);
            fmpz_mod_poly_set_fmpz_poly(Wm, W[lev] + m, ctx);
            fmpz_mod_poly_evaluate_fmpz(den, Wm, r, ctx);
            if (fmpz_is_zero(den))
                success = 0;
            else
            {
                fmpz_mod_inv(den, den, ctx);
                fmpz_mod_poly_zero(f, ctx);
                fmpz_mod_poly_set_coeff_ui(f, m, 1, ctx);
                for (jj = 0; jj < m; jj++)
                {
                    fmpz_mod_poly_set_fmpz_poly(Wm, W[lev] + jj, ctx);
                    fmpz_mod_poly_evaluate_fmpz(c, Wm, r, ctx);
                    fmpz_mod_mul(c, c, den, ctx);
                    fmpz_mod_poly_set_coeff_fmpz(f, jj, c, ctx);
                }
                nm1 = fmpz_fdiv_ui(fmpz_mod_ctx_modulus(ctx), m);
                if (kummer[lev] && nm1 == 1
                        && (_kummer_descent(c, W[lev], m, r, den, sqrtDn, f, state, ctx)
                            || (fmpz_mod_neg(sqrtDn, sqrtDn, ctx),
                                _kummer_descent(c, W[lev], m, r, den, sqrtDn, f, state, ctx))))
                {
                    fmpz_swap(r, c);
                    if (ecpp_verbose)
                        flint_printf("ecpp: Kummer descent in degree %wd\n", m);
                }
                else if (kummer[lev] && nm1 == (ulong) m - 1
                        && (_kummer_descent_f2(c, W[lev], m, r, den, sqrtDn, f, state, ctx)
                            || (fmpz_mod_neg(sqrtDn, sqrtDn, ctx),
                                _kummer_descent_f2(c, W[lev], m, r, den, sqrtDn, f, state, ctx))))
                {
                    fmpz_swap(r, c);
                    if (ecpp_verbose)
                        flint_printf("ecpp: Kummer descent in degree %wd through F_(n^2)\n", m);
                }
                else
                    success = ecpp_poly_root(r, f, state, ctx);
            }
            fmpz_mod_poly_clear(Wm, ctx);
        }
        if (success)
        {
            if (weber)
                _weber_to_j(j, r, D, ctx);
            else
                fmpz_set(j, r);
        }
        fmpz_mod_poly_clear(f, ctx);
        fmpz_clear(r); fmpz_clear(den); fmpz_clear(c);
        fmpz_clear(sqrtDn);
    }
#ifdef ECPP_PROFILE
    ecpp_tower_prof[3] += _tnow() - td0;
    }
#endif

    _acb_vec_clear(roots, h);
    arb_clear(sqrtD);
    acb_clear(z);
    for (lev = 0; lev < levels; lev++)
    {
        slong nd = (lev == 0) ? 1 : _level_size(degs[lev]);
        for (i = 0; i < nd; i++)
            fmpz_poly_clear(W[lev] + i);
        flint_free(W[lev]);
    }
    flint_free(W);
    flint_free(kummer);
cleanup_forms:
    flint_free(perm);
    flint_free(degs);
    qfb_array_clear(&forms, h);
    return success;
}
