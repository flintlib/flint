/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

/*
    Rational roots by p-adic lifting.

    Let g be squarefree with g(0) != 0 and let lc be its leading coefficient.
    If r = a/b is a root in lowest terms, then b | lc and a | g(0),
    so y = lc r is an integer with |y| <= |lc| R where R = min(|g(0)|, root
    bound). We choose a prime p not dividing lc such that the roots of g
    mod p are simple, find these roots, lift them to precision p^N with
    p^N > 2 |lc| R using Newton iteration, and recover the candidates
    y as symmetric remainders of lc r mod p^N. The candidates are finally
    verified exactly.

    This is much cheaper than factoring: the roots mod p can be found
    without factoring g mod p completely, and the precision needed is
    only about log(|lc| R) bits instead of the Mignotte bound.
*/

/* Starting point for the primes: large word-size primes make the
   lifting cheap (few Newton steps). */
#define ROOTS_FMPQ_PRIME_START (UWORD(1) << (FLINT_BITS - 4))

/* Number of primes to try (choosing the one with fewest roots). */
#define ROOTS_FMPQ_NUM_PRIMES 2

/* Use the rational root test with explicit enumeration of candidates when
   |g(0)| and |lc| have at most this many bits, and there are at most
   ROOTS_FMPQ_ENUM_MAX candidates. This must be at most
   SMALL_FMPZ_BITCOUNT_MAX, so that the coefficients are small fmpzs,
   and less than FLINT_BITS - 1, so that the enumeration prime (which is
   larger than 2^(FLINT_BITS - 1)) cannot divide a candidate. */
#define ROOTS_FMPQ_ENUM_BITS FLINT_MIN(40, SMALL_FMPZ_BITCOUNT_MAX)
#define ROOTS_FMPQ_ENUM_MAX 64

static ulong
fmpz_get_ui_abs_small(const fmpz_t x)
{
    FLINT_ASSERT(!COEFF_IS_MPZ(*x));
    return FLINT_ABS(*x);
}

/* Returns the number of roots of g mod p if p is suitable (p does not divide
   the leading coefficient and the roots of g mod p are simple), setting
   h to the product of (x - r) over these roots; returns -1 if p is not
   suitable. */
static slong
_roots_mod_p(nmod_poly_t h, const fmpz * g, slong len, ulong p)
{
    nmod_poly_t gp, t, u;
    slong res;

    nmod_poly_init(gp, p);
    nmod_poly_init(t, p);
    nmod_poly_init(u, p);

    nmod_poly_fit_length(gp, len);
    _fmpz_vec_get_nmod_vec(gp->coeffs, g, len, gp->mod);
    _nmod_poly_set_length(gp, len);
    _nmod_poly_normalise(gp);

    if (gp->length != len)
    {
        res = -1;
    }
    else
    {
        /* t = x^p mod gp */
        nmod_poly_reverse(u, gp, len);
        nmod_poly_inv_series(u, u, len);
        nmod_poly_powmod_x_ui_preinv(t, p, gp, u);

        /* h = gcd(gp, x^p - x) */
        nmod_poly_zero(u);
        nmod_poly_set_coeff_ui(u, 1, 1);
        nmod_poly_sub(t, t, u);
        nmod_poly_gcd(h, gp, t);

        res = nmod_poly_degree(h);

        if (res > 0)
        {
            /* require the roots to be simple roots of gp */
            nmod_poly_derivative(t, gp);
            nmod_poly_gcd(u, h, t);
            if (nmod_poly_degree(u) != 0)
                res = -1;
        }
    }

    nmod_poly_clear(gp);
    nmod_poly_clear(t);
    nmod_poly_clear(u);

    return res;
}

/* Data for early termination of the lifting: a candidate root found at
   low precision is accepted if it passes cheap necessary conditions and
   vanishes modulo two other primes. (All candidates are verified exactly
   at the end.) */
typedef struct
{
    const fmpz * g;
    slong len;
    const fmpz * R;     /* bound for |roots| */
    int integer_only;
    int pos;
    int neg;
    nmod_t q[2];
    nn_ptr gq[2];       /* g mod q[i] */
}
early_struct;

static void
early_init(early_struct * E, const fmpz * g, slong len, const fmpz_t R,
    int integer_only, int pos, int neg)
{
    slong i;
    ulong q = ROOTS_FMPQ_PRIME_START + 1000000;

    E->g = g;
    E->len = len;
    E->R = R;
    E->integer_only = integer_only;
    E->pos = pos;
    E->neg = neg;

    for (i = 0; i < 2; i++)
    {
        q = n_nextprime(q, 1);
        nmod_init(E->q + i, q);
        E->gq[i] = flint_malloc(len * sizeof(ulong));
        _fmpz_vec_get_nmod_vec(E->gq[i], g, len, E->q[i]);
    }
}

static void
early_clear(early_struct * E)
{
    flint_free(E->gq[0]);
    flint_free(E->gq[1]);
}

/* Checks whether the p-adic root r mod P reconstructs to a plausible
   rational root, which is then written to c. */
static int
early_check(fmpq_t c, const fmpz_t r, const fmpz_t P, const early_struct * E)
{
    const fmpz * g = E->g;
    slong len = E->len, i, k;
    ulong x, v, d;
    int s;

    if (E->integer_only)
    {
        fmpz_smod(fmpq_numref(c), r, P);
        fmpz_one(fmpq_denref(c));
    }
    else
    {
        if (!fmpq_reconstruct_fmpz(c, r, P))
            return 0;
        if (!fmpz_divisible(g + len - 1, fmpq_denref(c)))
            return 0;
    }

    s = fmpz_sgn(fmpq_numref(c));
    if (s == 0 || (s > 0 && !E->pos) || (s < 0 && !E->neg))
        return 0;

    /* |a / b| <= R */
    if (E->integer_only)
    {
        if (fmpz_cmpabs(fmpq_numref(c), E->R) > 0)
            return 0;
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mul(t, E->R, fmpq_denref(c));
        s = (fmpz_cmpabs(fmpq_numref(c), t) > 0);
        fmpz_clear(t);
        if (s)
            return 0;
        s = fmpz_sgn(fmpq_numref(c));
    }

    if (!fmpz_divisible(g + 0, fmpq_numref(c)))
        return 0;

    for (i = 0; i < 2; i++)
    {
        d = fmpz_get_nmod(fmpq_denref(c), E->q[i]);
        if (d == 0)
            return 0;
        x = nmod_mul(fmpz_get_nmod(fmpq_numref(c), E->q[i]),
            n_invmod(d, E->q[i].n), E->q[i]);

        v = E->gq[i][len - 1];
        for (k = len - 2; k >= 0; k--)
            v = nmod_add(nmod_mul(v, x, E->q[i]), E->gq[i][k], E->q[i]);
        if (v != 0)
            return 0;
    }

    return 1;
}

/* Lifts the n simple roots r of g mod p to roots mod p^N (in place).
   If E is not NULL, roots are checked for early termination at each
   precision; the accepted roots are written to res (their number is
   returned) and removed from r, and *n is updated. */
static slong
_lift_roots(fmpz * r, slong * n, fmpq * res, const fmpz * g, slong len,
    ulong p, slong N, const early_struct * E)
{
    slong * e;
    slong i, j, m, na, nres;
    fmpz_t P;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t G, dG;
    fmpz_poly_t gpoly;
    fmpz * y, * dy;

    na = *n;
    nres = 0;

    if (na == 0)
        return 0;

    /* precision sequence N, ceil(N/2), ..., 1 */
    e = flint_malloc(sizeof(slong) * (FLINT_BIT_COUNT(N) + 2));
    for (m = 0, e[0] = N; e[m] > 1; m++)
        e[m + 1] = (e[m] + 1) / 2;

    fmpz_init(P);
    y = _fmpz_vec_init(na);
    dy = _fmpz_vec_init(na);

    gpoly->coeffs = (fmpz *) g;
    gpoly->length = gpoly->alloc = len;

    fmpz_mod_ctx_init_ui(ctx, p);
    fmpz_mod_poly_init(G, ctx);
    fmpz_mod_poly_init(dG, ctx);

    fmpz_set_ui(P, p);

    for (j = m; ; j--)
    {
        /* here the roots are known mod P = p^e[j] */
        if (E != NULL)
        {
            for (i = 0; i < na; )
            {
                if (early_check(res + nres, r + i, P, E))
                {
                    nres++;
                    na--;
                    fmpz_swap(r + i, r + na);
                }
                else
                {
                    i++;
                }
            }
        }

        if (j == 0 || na == 0)
            break;

        /* lift from p^e[j] to p^e[j-1] */
        fmpz_set_ui(P, p);
        fmpz_pow_ui(P, P, e[j - 1]);
        fmpz_mod_ctx_set_modulus(ctx, P);

        fmpz_mod_poly_set_fmpz_poly(G, gpoly, ctx);
        fmpz_mod_poly_derivative(dG, G, ctx);

        fmpz_mod_poly_evaluate_fmpz_vec(y, G, r, na, ctx);
        fmpz_mod_poly_evaluate_fmpz_vec(dy, dG, r, na, ctx);

        for (i = 0; i < na; i++)
        {
            /* r <- r - g(r) / g'(r) */
            fmpz_mod_inv(dy + i, dy + i, ctx);
            fmpz_mod_mul(y + i, y + i, dy + i, ctx);
            fmpz_mod_sub(r + i, r + i, y + i, ctx);
        }
    }

    fmpz_mod_poly_clear(G, ctx);
    fmpz_mod_poly_clear(dG, ctx);
    fmpz_mod_ctx_clear(ctx);
    _fmpz_vec_clear(y, *n);
    _fmpz_vec_clear(dy, *n);
    fmpz_clear(P);
    flint_free(e);

    *n = na;
    return nres;
}

/* Sets d to the divisors of n and returns their number, or returns -1 if
   there are more than max. */
static slong
_small_divisors(ulong * d, ulong n, slong max)
{
    n_factor_t fac;
    slong i, j, k, num, cur;
    ulong pk;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    num = 1;
    for (i = 0; i < fac.num; i++)
    {
        num *= (fac.exp[i] + 1);
        if (num > max)
            return -1;
    }

    d[0] = 1;
    cur = 1;
    for (i = 0; i < fac.num; i++)
    {
        slong prev = cur;
        pk = 1;
        for (k = 1; k <= fac.exp[i]; k++)
        {
            pk *= fac.p[i];
            for (j = 0; j < prev; j++)
                d[cur++] = d[j] * pk;
        }
    }

    return num;
}

/* Is num/den a root of (g, len)? */
static int
_is_root_fmpq(const fmpz * g, slong len, const fmpz_t num, const fmpz_t den)
{
    fmpz_t rn, rd;
    int res;

    fmpz_init(rn);
    fmpz_init(rd);
    _fmpz_poly_evaluate_fmpq(rn, rd, g, len, num, den);
    res = fmpz_is_zero(rn);
    fmpz_clear(rn);
    fmpz_clear(rd);
    return res;
}

/* Test the candidates +/- a / b with a in d0, b in d1 (a, b coprime),
   first modulo a prime and then exactly. */
static slong
_roots_by_enumeration(fmpq * res, const fmpz * g, slong len,
    const ulong * d0, slong n0, const ulong * d1, slong n1, int pos, int neg)
{
    ulong p, x, v, * gp;
    nmod_t mod;
    slong i, j, k, num;
    int sign;
    fmpz_t a, b;

    /* p > 2^(FLINT_BITS - 1) > 2^ROOTS_FMPQ_ENUM_BITS, so p does not
       divide any candidate numerator or denominator */
    p = n_nextprime(UWORD(1) << (FLINT_BITS - 1), 1);
    nmod_init(&mod, p);
    gp = flint_malloc(len * sizeof(ulong));
    _fmpz_vec_get_nmod_vec(gp, g, len, mod);
    fmpz_init(a);
    fmpz_init(b);

    num = 0;
    for (i = 0; i < n0; i++)
    {
        for (j = 0; j < n1; j++)
        {
            if (n_gcd(d0[i], d1[j]) != 1)
                continue;

            for (sign = -1; sign <= 1; sign += 2)
            {
                if ((sign == 1 && !pos) || (sign == -1 && !neg))
                    continue;

                /* x = +/- a / b mod p (p does not divide b since b < p) */
                FLINT_ASSERT(d1[j] < p && d0[i] < p);
                x = nmod_mul(d0[i] % p, n_invmod(d1[j] % p, p), mod);
                if (sign == -1)
                    x = nmod_neg(x, mod);

                v = gp[len - 1];
                for (k = len - 2; k >= 0; k--)
                    v = nmod_add(nmod_mul(v, x, mod), gp[k], mod);

                if (v != 0)
                    continue;

                fmpz_set_ui(a, d0[i]);
                if (sign == -1)
                    fmpz_neg(a, a);
                fmpz_set_ui(b, d1[j]);

                if (_is_root_fmpq(g, len, a, b))
                {
                    fmpz_set(fmpq_numref(res + num), a);
                    fmpz_set(fmpq_denref(res + num), b);
                    num++;
                }
            }
        }
    }

    flint_free(gp);
    fmpz_clear(a);
    fmpz_clear(b);
    return num;
}

static slong
_roots_fmpq_squarefree(fmpq * res, const fmpz * g, slong len, int integer_only, int early)
{
    const fmpz * lc = g + len - 1;
    slong i, j, n, r_alloc, best_n, num, trial;
    ulong p, best_p;
    nmod_poly_t h, best_h;
    nmod_poly_factor_t fac;
    fmpz_t R, M, t, u, pN;
    fmpz * r;
    slong N;
    int s, sp, sn, pos, neg;

    if (len <= 1)
        return 0;

    if (len == 2)
    {
        if (integer_only)
        {
            if (!fmpz_divisible(g + 0, g + 1))
                return 0;
            fmpz_divexact(fmpq_numref(res), g + 0, g + 1);
            fmpz_neg(fmpq_numref(res), fmpq_numref(res));
            fmpz_one(fmpq_denref(res));
        }
        else
        {
            fmpz_neg(fmpq_numref(res), g + 0);
            fmpz_set(fmpq_denref(res), g + 1);
            fmpq_canonicalise(res);
        }
        return 1;
    }

    /* Descartes' rule of signs: can there be positive / negative roots?
       (sign changes in the coefficients of g(x) and g(-x)) */
    sp = sn = fmpz_sgn(g);
    pos = neg = 0;
    for (i = 1; i < len && !(pos && neg); i++)
    {
        s = fmpz_sgn(g + i);
        if (s == 0)
            continue;
        if (s != sp)
            pos = 1;
        if (((i % 2 == 1) ? -s : s) != sn)
            neg = 1;
    }

    if (!pos && !neg)
        return 0;

    /* If the constant and leading coefficients are small with few divisors,
       simply test the candidates a/b with a | g(0), b | lc. */
    if (fmpz_bits(g + 0) <= ROOTS_FMPQ_ENUM_BITS && fmpz_bits(lc) <= ROOTS_FMPQ_ENUM_BITS)
    {
        ulong d0[ROOTS_FMPQ_ENUM_MAX], d1[ROOTS_FMPQ_ENUM_MAX];
        slong n0, n1;

        n0 = _small_divisors(d0, fmpz_get_ui_abs_small(g + 0), ROOTS_FMPQ_ENUM_MAX);
        n1 = integer_only ? 1 : _small_divisors(d1, fmpz_get_ui_abs_small(lc), ROOTS_FMPQ_ENUM_MAX);
        if (integer_only)
            d1[0] = 1;

        if (n0 > 0 && n1 > 0 && n0 * n1 * (pos + neg) <= ROOTS_FMPQ_ENUM_MAX)
            return _roots_by_enumeration(res, g, len, d0, n0, d1, n1, pos, neg);
    }

    /* Find a suitable prime with few roots. */
    best_n = -1;
    best_p = 0;
    p = ROOTS_FMPQ_PRIME_START;

    for (trial = 0; trial < ROOTS_FMPQ_NUM_PRIMES; )
    {
        p = n_nextprime(p, 1);
        nmod_poly_init(h, p);
        n = _roots_mod_p(h, g, len, p);

        if (n >= 0 && (best_n < 0 || n < best_n))
        {
            if (best_n >= 0)
                nmod_poly_clear(best_h);
            best_n = n;
            best_p = p;
            *best_h = *h;   /* move */
        }
        else
        {
            nmod_poly_clear(h);
        }

        if (n >= 0)
            trial++;

        if (best_n == 0)
            break;
    }

    if (best_n == 0)
    {
        nmod_poly_clear(best_h);
        return 0;
    }

    p = best_p;

    /* Bound R for |roots| */
    fmpz_init(R);
    fmpz_init(M);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(pN);

    {
        fmpz_poly_t gpoly;
        gpoly->coeffs = (fmpz *) g;
        gpoly->length = gpoly->alloc = len;
        fmpz_poly_bound_roots(R, gpoly);
        fmpz_add_ui(R, R, 1);
    }
    fmpz_abs(t, g + 0);
    if (fmpz_cmp(t, R) < 0)
        fmpz_set(R, t);

    /* candidates y = lc r satisfy |y| <= M = |lc| R */
    if (integer_only)
    {
        fmpz_set(M, R);
    }
    else
    {
        fmpz_mul(M, R, lc);
        fmpz_abs(M, M);
    }

    /* p^N > 2 M */
    fmpz_mul_2exp(t, M, 1);
    N = 0;
    fmpz_one(pN);
    while (fmpz_cmp(pN, t) <= 0)
    {
        fmpz_mul_ui(pN, pN, p);
        N++;
    }

    /* roots mod p */
    nmod_poly_factor_init(fac);
    nmod_poly_roots(fac, best_h, 0);
    n = fac->num;
    /* _lift_roots moves accepted roots to the end of r and decreases n,
       so remember the allocated length */
    r_alloc = n;
    r = _fmpz_vec_init(r_alloc);
    for (i = 0; i < n; i++)
        fmpz_set_ui(r + i, nmod_neg(fac->p[i].coeffs[0], best_h->mod));
    nmod_poly_factor_clear(fac);
    nmod_poly_clear(best_h);

    if (early)
    {
        early_struct E;
        early_init(&E, g, len, R, integer_only, pos, neg);
        num = _lift_roots(r, &n, res, g, len, p, N, &E);
        early_clear(&E);
    }
    else
    {
        num = _lift_roots(r, &n, NULL, g, len, p, N, NULL);
    }

    /* reconstruct and filter the remaining candidates */
    for (i = 0; i < n; i++)
    {
        fmpz * y = r + i;

        if (!integer_only)
        {
            fmpz_mul(y, y, lc);
            fmpz_mod(y, y, pN);
        }
        fmpz_smod(y, y, pN);

        if (fmpz_cmpabs(y, M) > 0)
            continue;

        if (fmpz_sgn(y) > 0 && !pos)
            continue;
        if (fmpz_sgn(y) < 0 && !neg)
            continue;

        if (integer_only)
        {
            /* the root must divide g(0) */
            if (fmpz_is_zero(y) || !fmpz_divisible(g + 0, y))
                continue;
            fmpz_set(fmpq_numref(res + num), y);
            fmpz_one(fmpq_denref(res + num));
        }
        else
        {
            fmpz_set(fmpq_numref(res + num), y);
            fmpz_set(fmpq_denref(res + num), lc);
            fmpq_canonicalise(res + num);
            /* the numerator must divide g(0) */
            if (fmpz_is_zero(fmpq_numref(res + num)) || !fmpz_divisible(g + 0, fmpq_numref(res + num)))
                continue;
        }

        num++;
    }

    _fmpz_vec_clear(r, r_alloc);

    /* Verify the candidates: since they are distinct, all of them are roots
       iff the product of the linear factors divides g. */
    if (num != 0)
    {
        fmpz_poly_t P, Q, gpoly;
        int ok;

        fmpz_poly_init(P);
        fmpz_poly_init(Q);
        gpoly->coeffs = (fmpz *) g;
        gpoly->length = gpoly->alloc = len;

        fmpz_poly_product_roots_fmpq_vec(P, res, num);
        ok = fmpz_poly_divides(Q, gpoly, P);

        if (!ok)
        {
            for (i = j = 0; i < num; i++)
            {
                if (_is_root_fmpq(g, len, fmpq_numref(res + i), fmpq_denref(res + i)))
                {
                    fmpq_swap(res + j, res + i);
                    j++;
                }
            }

            /* A candidate accepted by early termination could have been
               wrong (extremely unlikely); then redo the computation without
               early termination, so that no root can be missed. */
            if (early && j != num)
                num = -1;
            else
                num = j;
        }

        fmpz_poly_clear(P);
        fmpz_poly_clear(Q);
    }

    fmpz_clear(R);
    fmpz_clear(M);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(pN);

    if (num == -1)
        num = _roots_fmpq_squarefree(res, g, len, integer_only, 0);

    return num;
}

slong
_fmpz_poly_roots_fmpq_squarefree(fmpq * res, const fmpz * poly, slong len)
{
    return _roots_fmpq_squarefree(res, poly, len, 0, 1);
}

slong
_fmpz_poly_roots_fmpz_squarefree(fmpz * res, const fmpz * poly, slong len)
{
    slong i, num, alloc = FLINT_MAX(len - 1, 1);
    fmpq * tmp;

    tmp = _fmpq_vec_init(alloc);
    num = _roots_fmpq_squarefree(tmp, poly, len, 1, 1);
    for (i = 0; i < num; i++)
        fmpz_swap(res + i, fmpq_numref(tmp + i));
    _fmpq_vec_clear(tmp, alloc);

    return num;
}

/* Sort roots (with multiplicities) in increasing order */
static void
_sort_roots(fmpq * res, slong * exp, slong num)
{
    slong i, j;

    /* insertion sort; num is typically small */
    for (i = 1; i < num; i++)
    {
        for (j = i; j > 0 && fmpq_cmp(res + j - 1, res + j) > 0; j--)
        {
            fmpq_swap(res + j - 1, res + j);
            if (exp != NULL)
                FLINT_SWAP(slong, exp[j - 1], exp[j]);
        }
    }
}

static slong
_fmpz_poly_roots_fmpq_generic(fmpq * res, slong * exp, const fmpz_poly_t poly, int integer_only)
{
    slong len = poly->length;
    slong v, i, j, num, n;
    fmpz_poly_t f;
    int squarefree;

    if (len == 0)
        flint_throw(FLINT_ERROR, "(%s): zero polynomial\n", __func__);

    /* roots at zero */
    for (v = 0; fmpz_is_zero(poly->coeffs + v); v++)
        ;

    num = 0;
    if (v > 0)
    {
        fmpq_zero(res);
        if (exp != NULL)
            exp[0] = v;
        num = 1;
    }

    if (len - v <= 1)
        return num;

    fmpz_poly_init(f);
    fmpz_poly_shift_right(f, poly, v);
    fmpz_poly_primitive_part(f, f);

    /* Quick squarefreeness test modulo a prime */
    {
        ulong p = n_nextprime(ROOTS_FMPQ_PRIME_START + 12345, 1);
        nmod_poly_t fp, df, g;

        nmod_poly_init(fp, p);
        nmod_poly_init(df, p);
        nmod_poly_init(g, p);
        fmpz_poly_get_nmod_poly(fp, f);
        squarefree = 0;
        if (fp->length == f->length)
        {
            nmod_poly_derivative(df, fp);
            nmod_poly_gcd(g, fp, df);
            squarefree = (nmod_poly_degree(g) == 0);
        }
        nmod_poly_clear(fp);
        nmod_poly_clear(df);
        nmod_poly_clear(g);
    }

    if (squarefree)
    {
        n = _roots_fmpq_squarefree(res + num, f->coeffs, f->length, integer_only, 1);
        if (exp != NULL)
            for (i = 0; i < n; i++)
                exp[num + i] = 1;
        num += n;
    }
    else
    {
        fmpz_poly_factor_t fac;

        fmpz_poly_factor_init(fac);
        fmpz_poly_factor_squarefree(fac, f);

        for (i = 0; i < fac->num; i++)
        {
            n = _roots_fmpq_squarefree(res + num,
                fac->p[i].coeffs, fac->p[i].length, integer_only, 1);
            if (exp != NULL)
                for (j = 0; j < n; j++)
                    exp[num + j] = fac->exp[i];
            num += n;
        }

        fmpz_poly_factor_clear(fac);
    }

    fmpz_poly_clear(f);

    _sort_roots(res, exp, num);

    return num;
}

slong
fmpz_poly_roots_fmpq(fmpq * res, slong * exp, const fmpz_poly_t poly)
{
    return _fmpz_poly_roots_fmpq_generic(res, exp, poly, 0);
}

slong
fmpz_poly_roots_fmpz(fmpz * res, slong * exp, const fmpz_poly_t poly)
{
    slong i, num, len = poly->length;
    fmpq * tmp;

    if (len == 0)
        flint_throw(FLINT_ERROR, "(%s): zero polynomial\n", __func__);

    tmp = _fmpq_vec_init(FLINT_MAX(len - 1, 1));
    num = _fmpz_poly_roots_fmpq_generic(tmp, exp, poly, 1);
    for (i = 0; i < num; i++)
        fmpz_swap(res + i, fmpq_numref(tmp + i));
    _fmpq_vec_clear(tmp, FLINT_MAX(len - 1, 1));

    return num;
}
