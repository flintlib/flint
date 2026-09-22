/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "fmpz_factor.h"
#include "qfb.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_ec.h"
#include "impl.h"

/*
    Counting the points of E(F_q).

    Three algorithms, all returning #E(F_q) = q + 1 - t with |t| <= 2 sqrt(q):

      naive   O(q) field operations; walks over every x and counts the
              roots in y. A reference implementation, and the only one of
              the three that is completely independent of the group law.

      bsgs    O(q^(1/4)) group operations, Shanks and Mestre: find a
              multiple of the order of a random point inside the Hasse
              interval by baby-step giant-step, take the exact order from
              it, and repeat with further points until only one multiple
              of the accumulated lcm is left in the interval.

      schoof  polynomial in log q; t mod l from the action of Frobenius on
              the l-torsion, for enough small primes l to pin t down by the
              Chinese remainder theorem.
*/

/* ------------------------------------------------------------------ */
/* walking over the elements of a finite field                        */
/* ------------------------------------------------------------------ */

/*
    F_q = F_p[g] with g the generator, so the elements are exactly the
    sums c_0 + c_1 g + ... + c_{d-1} g^{d-1} with c_i in [0, p). We walk
    over them with an odometer on the digits.
*/
typedef struct
{
    gr_ctx_struct * R;
    gr_ptr pows;            /* 1, g, ..., g^(d-1) */
    ulong * digits;
    slong d;
    ulong p;
}
field_iter_struct;

static void
field_iter_clear(field_iter_struct * it)
{
    if (it->pows != NULL)
        gr_heap_clear_vec(it->pows, it->d, it->R);
    flint_free(it->digits);
}

static int
field_iter_init(field_iter_struct * it, gr_ctx_t R)
{
    fmpz_t p;
    slong i, sz = R->sizeof_elem;
    int status = GR_SUCCESS;

    it->R = R;
    it->pows = NULL;
    it->digits = NULL;

    fmpz_init(p);

    if (gr_ctx_fq_prime(p, R) != GR_SUCCESS || !fmpz_abs_fits_ui(p)
            || gr_ctx_fq_degree(&it->d, R) != GR_SUCCESS || it->d < 1)
    {
        fmpz_clear(p);
        return GR_UNABLE;
    }

    it->p = fmpz_get_ui(p);
    fmpz_clear(p);

    it->digits = flint_calloc(it->d, sizeof(ulong));
    it->pows = gr_heap_init_vec(it->d, R);

    status |= gr_one(it->pows, R);

    if (it->d > 1)
    {
        gr_ptr g;
        GR_TMP_INIT(g, R);
        status |= gr_gen(g, R);

        for (i = 1; i < it->d; i++)
            status |= gr_mul(GR_ENTRY(it->pows, i, sz),
                        GR_ENTRY(it->pows, i - 1, sz), g, R);

        GR_TMP_CLEAR(g, R);
    }

    if (status != GR_SUCCESS)
        field_iter_clear(it);

    return status;
}

/* the element for the current digits */
static int
field_iter_get(gr_ptr res, const field_iter_struct * it)
{
    gr_ctx_struct * R = it->R;
    slong i, sz = R->sizeof_elem;
    int status = GR_SUCCESS;

    if (it->d == 1)
        return gr_set_ui(res, it->digits[0], R);

    status |= gr_zero(res, R);

    for (i = 0; i < it->d; i++)
        if (it->digits[i] != 0)
        {
            gr_ptr t;
            GR_TMP_INIT(t, R);
            status |= gr_mul_ui(t, GR_ENTRY(it->pows, i, sz), it->digits[i], R);
            status |= gr_add(res, res, t, R);
            GR_TMP_CLEAR(t, R);
        }

    return status;
}

/* returns 0 once it has wrapped around, that is, once every element is done */
static int
field_iter_next(field_iter_struct * it)
{
    slong i;

    for (i = 0; i < it->d; i++)
    {
        it->digits[i]++;

        if (it->digits[i] < it->p)
            return 1;

        it->digits[i] = 0;
    }

    return 0;
}

/* ------------------------------------------------------------------ */
/* naive counting                                                     */
/* ------------------------------------------------------------------ */

/*
    Refuse to walk over a field bigger than this; the caller asked for the
    O(q) algorithm, but not for it to never return.
*/
#define GR_EC_NAIVE_MAX_Q WORD(100000000)

int
gr_ec_ctx_cardinality_naive(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    field_iter_struct it;
    fmpz_t q, p, count;
    gr_ptr x, b, c, t;
    slong d;
    int char_two, status = GR_SUCCESS;

    if (gr_ctx_is_field(R) != T_TRUE)
        return GR_DOMAIN;

    fmpz_init(q);
    fmpz_init(p);
    fmpz_init(count);

    if (gr_ctx_cardinality_fmpz(q, R) != GR_SUCCESS
            || gr_ctx_fq_prime(p, R) != GR_SUCCESS
            || !fmpz_fits_si(q) || fmpz_cmp_si(q, GR_EC_NAIVE_MAX_Q) > 0)
    {
        fmpz_clear(q); fmpz_clear(p); fmpz_clear(count);
        return GR_UNABLE;
    }

    char_two = (fmpz_cmp_ui(p, 2) == 0);

    if (gr_ctx_fq_degree(&d, R) != GR_SUCCESS)
    {
        fmpz_clear(q); fmpz_clear(p); fmpz_clear(count);
        return GR_UNABLE;
    }

    status = field_iter_init(&it, R);

    if (status != GR_SUCCESS)
    {
        fmpz_clear(q); fmpz_clear(p); fmpz_clear(count);
        return status;
    }

    GR_TMP_INIT4(x, b, c, t, R);

    fmpz_one(count);            /* the point at infinity */

    do
    {
        status |= field_iter_get(x, &it);

        /* b = a1 x + a3 */
        status |= gr_mul(b, GR_EC_A1(ctx), x, R);
        status |= gr_add(b, b, GR_EC_A3(ctx), R);

        /* c = x^3 + a2 x^2 + a4 x + a6, by Horner */
        status |= gr_add(c, x, GR_EC_A2(ctx), R);
        status |= gr_mul(c, c, x, R);
        status |= gr_add(c, c, GR_EC_A4(ctx), R);
        status |= gr_mul(c, c, x, R);
        status |= gr_add(c, c, GR_EC_A6(ctx), R);

        if (status != GR_SUCCESS)
            break;

        /* number of y with y^2 + b y - c = 0 */
        if (!char_two)
        {
            /* the discriminant of the quadratic, b^2 + 4c */
            status |= gr_sqr(t, b, R);
            status |= gr_mul_ui(c, c, 4, R);
            status |= gr_add(t, t, c, R);

            if (status != GR_SUCCESS)
                break;

            if (gr_is_zero(t, R) == T_TRUE)
                fmpz_add_ui(count, count, 1);
            else if (gr_is_square(t, R) == T_TRUE)
                fmpz_add_ui(count, count, 2);
        }
        else if (gr_is_zero(b, R) == T_TRUE)
        {
            /* y^2 = c, and squaring is a bijection in characteristic 2 */
            fmpz_add_ui(count, count, 1);
        }
        else
        {
            /*
                y = b z turns y^2 + b y = c into z^2 + z = c / b^2, which
                is soluble exactly when the absolute trace vanishes, and
                then has two solutions.
            */
            fmpz_t tr;

            status |= gr_sqr(t, b, R);
            status |= gr_div(t, c, t, R);

            fmpz_init(tr);

            if (d == 1)
                status |= gr_get_fmpz(tr, t, R);
            else
                status |= gr_fq_trace(tr, t, R);

            if (status == GR_SUCCESS && fmpz_is_even(tr))
                fmpz_add_ui(count, count, 2);

            fmpz_clear(tr);

            if (status != GR_SUCCESS)
                break;
        }
    }
    while (field_iter_next(&it));

    if (status == GR_SUCCESS)
        fmpz_set(res, count);

    GR_TMP_CLEAR4(x, b, c, t, R);
    field_iter_clear(&it);
    fmpz_clear(q);
    fmpz_clear(p);
    fmpz_clear(count);

    return status;
}

/* ------------------------------------------------------------------ */
/* baby-step giant-step (Shanks, Mestre)                              */
/* ------------------------------------------------------------------ */

/*
    Points have no canonical ordering in a general finite field, so the
    baby-step table is keyed on a hash of the decimal (or polynomial)
    string of the x-coordinate. A hit is only a candidate; it is always
    confirmed by testing the resulting order against the point itself, so
    a collision costs a little time and never an answer.
*/
typedef struct
{
    ulong hash;
    slong j;
}
baby_struct;

static ulong
_hash_elem(gr_srcptr x, gr_ctx_t R, int * status)
{
    ulong h = UWORD(14695981039346656037);
    char * s;
    slong i;

    if (gr_get_str(&s, x, R) != GR_SUCCESS)
    {
        *status |= GR_UNABLE;
        return 0;
    }

    for (i = 0; s[i] != '\0'; i++)
    {
        h ^= (ulong) (unsigned char) s[i];
        h *= UWORD(1099511628211);
    }

    flint_free(s);

    return h;
}

static int
_baby_cmp(const void * a, const void * b)
{
    ulong x = ((const baby_struct *) a)->hash;
    ulong y = ((const baby_struct *) b)->hash;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* first index with this hash, or -1 */
static slong
_baby_find(const baby_struct * tab, slong len, ulong h)
{
    slong lo = 0, hi = len - 1, best = -1;

    while (lo <= hi)
    {
        slong mid = lo + (hi - lo) / 2;

        if (tab[mid].hash == h)
        {
            best = mid;
            hi = mid - 1;
        }
        else if (tab[mid].hash < h)
            lo = mid + 1;
        else
            hi = mid - 1;
    }

    return best;
}

/*
    The exact order of P, given any multiple n of it: strip one prime at a
    time for as long as the result still kills P.
*/
static int
_point_order(fmpz_t ord, const gr_ec_point_t P, const fmpz_t n, gr_ec_ctx_t ctx)
{
    fmpz_factor_t fac;
    gr_ec_point_t T;
    fmpz_t red;
    slong i;
    int status = GR_SUCCESS;

    fmpz_factor_init(fac);
    fmpz_factor(fac, n);
    fmpz_init(red);
    gr_ec_point_init(T, ctx);

    fmpz_set(ord, n);

    for (i = 0; i < fac->num && status == GR_SUCCESS; i++)
    {
        slong e;

        for (e = 0; e < fac->exp[i]; e++)
        {
            fmpz_divexact(red, ord, fac->p + i);

            status |= gr_ec_point_mul_fmpz(T, P, red, ctx);

            if (status != GR_SUCCESS)
                break;

            if (gr_ec_point_is_inf(T, ctx) == T_TRUE)
                fmpz_set(ord, red);
            else
                break;
        }
    }

    gr_ec_point_clear(T, ctx);
    fmpz_clear(red);
    fmpz_factor_clear(fac);

    return status;
}

/* how many attempts with fresh points before giving up */
#define GR_EC_BSGS_MAX_POINTS 40

int
gr_ec_ctx_cardinality_bsgs(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ec_point_t P, Q, mP, S, Sp, Sm, T;
    gr_ptr xc, yc;
    baby_struct * tab = NULL;
    flint_rand_t state;
    fmpz_t q, M, lo, hi, L, n, ord, t, cand, kl, kh;
    slong m, i, j, attempt, ngiant;
    int status = GR_SUCCESS, resolved = 0;

    if (gr_ctx_is_field(R) != T_TRUE)
        return GR_DOMAIN;

    fmpz_init(q); fmpz_init(M); fmpz_init(lo); fmpz_init(hi);
    fmpz_init(L); fmpz_init(n); fmpz_init(ord); fmpz_init(t);
    fmpz_init(cand); fmpz_init(kl); fmpz_init(kh);

    if (gr_ctx_cardinality_fmpz(q, R) != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup_ints;
    }

    /* the Hasse interval [q + 1 - 2 sqrt(q), q + 1 + 2 sqrt(q)] */
    fmpz_sqrt(M, q);
    fmpz_mul_ui(M, M, 2);
    fmpz_add_ui(M, M, 2);              /* round up, a little slack is free */

    fmpz_add_ui(lo, q, 1);
    fmpz_sub(hi, lo, M);
    fmpz_swap(lo, hi);                 /* lo = q + 1 - M */
    fmpz_add(hi, q, M);
    fmpz_add_ui(hi, hi, 1);            /* hi = q + 1 + M */

    /* m = ceil(sqrt(2M + 1)), so that i and j both range over about m */
    {
        fmpz_t w;
        fmpz_init(w);
        fmpz_mul_ui(w, M, 2);
        fmpz_add_ui(w, w, 1);
        fmpz_sqrt(w, w);
        fmpz_add_ui(w, w, 1);

        if (!fmpz_fits_si(w))
        {
            fmpz_clear(w);
            status = GR_UNABLE;         /* q far too large for BSGS */
            goto cleanup_ints;
        }

        m = fmpz_get_si(w);
        fmpz_clear(w);
    }

    flint_rand_init(state);
    flint_rand_set_seed(state, UWORD(0x9e3779b97f4a7c15), UWORD(0xbf58476d1ce4e5b9));

    gr_ec_point_init(P, ctx);
    gr_ec_point_init(Q, ctx);
    gr_ec_point_init(mP, ctx);
    gr_ec_point_init(S, ctx);
    gr_ec_point_init(Sp, ctx);
    gr_ec_point_init(Sm, ctx);
    gr_ec_point_init(T, ctx);
    GR_TMP_INIT2(xc, yc, R);

    tab = flint_malloc(m * sizeof(baby_struct));

    fmpz_one(L);

    for (attempt = 0; attempt < GR_EC_BSGS_MAX_POINTS && !resolved; attempt++)
    {
        slong ntab = 0;
        int found = 0;

        if (gr_ec_point_randtest(P, state, ctx) != GR_SUCCESS)
            continue;

        if (gr_ec_point_is_inf(P, ctx) != T_FALSE)
            continue;

        /* baby steps: the x-coordinates of jP for j = 1, ..., m - 1 */
        status |= gr_ec_point_set(S, P, ctx);

        for (j = 1; j < m && status == GR_SUCCESS; j++)
        {
            if (gr_ec_point_is_inf(S, ctx) == T_TRUE)
            {
                /* jP = O already: the order of P divides j */
                fmpz_set_si(n, j);
                found = 2;
                break;
            }

            if (gr_ec_point_get_affine(xc, yc, S, ctx) == GR_SUCCESS)
            {
                tab[ntab].hash = _hash_elem(xc, R, &status);
                tab[ntab].j = j;
                ntab++;
            }

            status |= gr_ec_point_add(S, S, P, ctx);
        }

        if (status != GR_SUCCESS)
            break;

        if (found != 2)
        {
            qsort(tab, ntab, sizeof(baby_struct), _baby_cmp);

            /* Q = (q + 1) P, and mP for the giant strides */
            fmpz_add_ui(t, q, 1);
            status |= gr_ec_point_mul_fmpz(Q, P, t, ctx);
            fmpz_set_si(t, m);
            status |= gr_ec_point_mul_fmpz(mP, P, t, ctx);

            /* Sp = Q - i (mP) and Sm = Q + i (mP), both built up a stride
               at a time so that a giant step costs one addition */
            status |= gr_ec_point_set(Sp, Q, ctx);
            status |= gr_ec_point_set(Sm, Q, ctx);

            if (status != GR_SUCCESS)
                break;

            ngiant = (slong) (fmpz_get_d(M) / m) + 2;

            for (i = 0; i <= ngiant && !found && status == GR_SUCCESS; i++)
            {
                int which;

                for (which = 0; which < 2 && !found; which++)
                {
                    gr_ec_point_struct * S_i = which ? Sm : Sp;
                    slong im = (which ? -1 : 1) * i * m;
                    slong idx;
                    ulong h;

                    if (which && i == 0)
                        continue;               /* Sp and Sm coincide */

                    if (gr_ec_point_is_inf(S_i, ctx) == T_TRUE)
                    {
                        /* Q = (i m) P with the relevant sign, so t = im */
                        fmpz_set_si(t, im);
                        fmpz_add_ui(n, q, 1);
                        fmpz_sub(n, n, t);
                        found = 1;
                        break;
                    }

                    if (gr_ec_point_get_affine(xc, yc, S_i, ctx) != GR_SUCCESS)
                        continue;

                    h = _hash_elem(xc, R, &status);
                    idx = _baby_find(tab, ntab, h);

                    /* S_i = +- jP, so t = im +- j; confirm which, if either */
                    for ( ; idx >= 0 && idx < ntab && tab[idx].hash == h && !found; idx++)
                    {
                        int sj;

                        for (sj = -1; sj <= 1 && !found; sj += 2)
                        {
                            fmpz_set_si(t, im + sj * tab[idx].j);
                            fmpz_add_ui(cand, q, 1);
                            fmpz_sub(cand, cand, t);

                            if (fmpz_cmp(cand, lo) < 0 || fmpz_cmp(cand, hi) > 0
                                    || fmpz_sgn(cand) <= 0)
                                continue;

                            status |= gr_ec_point_mul_fmpz(T, P, cand, ctx);

                            if (status == GR_SUCCESS
                                    && gr_ec_point_is_inf(T, ctx) == T_TRUE)
                            {
                                fmpz_set(n, cand);
                                found = 1;
                            }
                        }
                    }
                }

                status |= gr_ec_point_sub(Sp, Sp, mP, ctx);
                status |= gr_ec_point_add(Sm, Sm, mP, ctx);
            }
        }

        if (status != GR_SUCCESS)
            break;

        if (!found)
            continue;

        /* the exact order of P, then the lcm over all points so far */
        status |= _point_order(ord, P, n, ctx);

        if (status != GR_SUCCESS)
            break;

        fmpz_lcm(L, L, ord);

        /* is there exactly one multiple of L left in the Hasse interval? */
        fmpz_cdiv_q(kl, lo, L);
        fmpz_fdiv_q(kh, hi, L);

        if (fmpz_equal(kl, kh))
        {
            fmpz_mul(res, L, kl);
            resolved = 1;
        }
    }

    if (status == GR_SUCCESS && !resolved)
        status = GR_UNABLE;

    flint_free(tab);
    GR_TMP_CLEAR2(xc, yc, R);
    gr_ec_point_clear(T, ctx);
    gr_ec_point_clear(Sm, ctx);
    gr_ec_point_clear(Sp, ctx);
    gr_ec_point_clear(S, ctx);
    gr_ec_point_clear(mP, ctx);
    gr_ec_point_clear(Q, ctx);
    gr_ec_point_clear(P, ctx);
    flint_rand_clear(state);

cleanup_ints:
    fmpz_clear(q); fmpz_clear(M); fmpz_clear(lo); fmpz_clear(hi);
    fmpz_clear(L); fmpz_clear(n); fmpz_clear(ord); fmpz_clear(t);
    fmpz_clear(cand); fmpz_clear(kl); fmpz_clear(kh);

    return status;
}

/* ------------------------------------------------------------------ */
/* Schoof's algorithm                                                 */
/* ------------------------------------------------------------------ */

/*
    A deliberately simple Schoof: short Weierstrass y^2 = f(x) = x^3 + a4 x
    + a6 over F_q with q coprime to 6.

    Everything happens in R = F_q[x]/(psi_l), where psi_l is the l-division
    polynomial. A point of the l-torsion is written (u, v y) with u, v in R,
    which is closed under the group law because y^2 = f reduces any even
    power of y. Frobenius is phi(u, v y) = (u o xq, (v o xq) yq y), where
    xq = x^q and yq = f^((q-1)/2), both taken mod psi_l.

    The trace satisfies phi^2 - t phi + q = 0 on the l-torsion, so we form
    phi^2(P) + (q mod l) P and match it against t' phi(P) for t' = 0, 1,
    ..., (l-1)/2, taking the sign from the y-coordinate. Enough l that
    their product exceeds 4 sqrt(q) determines t by the Chinese remainder
    theorem, and #E = q + 1 - t.

    R is not a field when psi_l is reducible, so an inverse can fail. A
    real implementation would take the factor this hands it and carry on
    modulo that (this is where Elkies primes come from); this one reports
    GR_UNABLE and lets the caller fall back.
*/

/* a point of the l-torsion as (u, v y), or the identity */
typedef struct
{
    gr_poly_t u;
    gr_poly_t v;
    int is_inf;
}
tors_struct;

typedef struct
{
    gr_ctx_struct * R;
    gr_poly_preinv_t P;         /* preconditioned psi_l */
    gr_poly_t psi;              /* psi_l itself, for the xgcd in tors_inv */
    gr_poly_t f;                /* x^3 + a4 x + a6 mod psi_l */
    gr_poly_t a4;
    gr_poly_t factor;           /* a proper factor of psi, when one shows up */
    int have_factor;
}
tors_ctx_struct;

static void
tors_init(tors_struct * T, tors_ctx_struct * C)
{
    gr_poly_init(T->u, C->R);
    gr_poly_init(T->v, C->R);
    T->is_inf = 1;
}

static void
tors_clear(tors_struct * T, tors_ctx_struct * C)
{
    gr_poly_clear(T->u, C->R);
    gr_poly_clear(T->v, C->R);
}

static int
tors_set(tors_struct * D, const tors_struct * S, tors_ctx_struct * C)
{
    int status = GR_SUCCESS;
    status |= gr_poly_set(D->u, S->u, C->R);
    status |= gr_poly_set(D->v, S->v, C->R);
    D->is_inf = S->is_inf;
    return status;
}

static int
tors_mulmod(gr_poly_t res, const gr_poly_t a, const gr_poly_t b, tors_ctx_struct * C)
{
    return gr_poly_preinv_mulmod(res, a, b, C->P, C->R);
}

/* 1/a in R, or GR_UNABLE when psi_l turns out to be reducible here */
static int
tors_inv(gr_poly_t res, const gr_poly_t a, tors_ctx_struct * C)
{
    gr_poly_t g, s, t;
    int status = GR_SUCCESS;

    gr_poly_init(g, C->R);
    gr_poly_init(s, C->R);
    gr_poly_init(t, C->R);

    status |= gr_poly_xgcd(g, s, t, a, C->psi, C->R);

    if (status == GR_SUCCESS)
    {
        if (gr_poly_length(g, C->R) != 1)
        {
            /*
                The xgcd has handed us a nontrivial factor of the modulus.
                phi^2 - t phi + q vanishes on the whole l-torsion, so it
                still vanishes modulo this factor, and the caller can start
                again with the smaller modulus. (Recognising that factor as
                a kernel polynomial is exactly what turns Schoof into SEA.)
            */
            if (gr_poly_length(g, C->R) < gr_poly_length(C->psi, C->R)
                    && !C->have_factor)
            {
                if (gr_poly_set(C->factor, g, C->R) == GR_SUCCESS)
                    C->have_factor = 1;
            }

            status = GR_UNABLE;
        }
        else
            status = gr_poly_preinv_rem(res, s, C->P, C->R);
    }

    gr_poly_clear(g, C->R);
    gr_poly_clear(s, C->R);
    gr_poly_clear(t, C->R);

    return status;
}

static truth_t
tors_equal(const tors_struct * A, const tors_struct * B, tors_ctx_struct * C)
{
    if (A->is_inf || B->is_inf)
        return (A->is_inf && B->is_inf) ? T_TRUE : T_FALSE;

    if (gr_poly_equal(A->u, B->u, C->R) != T_TRUE)
        return T_FALSE;

    return gr_poly_equal(A->v, B->v, C->R);
}

static int
tors_neg(tors_struct * D, const tors_struct * S, tors_ctx_struct * C)
{
    int status = GR_SUCCESS;
    status |= gr_poly_set(D->u, S->u, C->R);
    status |= gr_poly_neg(D->v, S->v, C->R);
    D->is_inf = S->is_inf;
    return status;
}

/*
    The chord and tangent law written for (u, v y) with y^2 = f:

      distinct: w = (v2 - v1)/(u2 - u1),  u3 = w^2 f - u1 - u2,
                v3 = w (u1 - u3) - v1
      doubling: w = (3 u^2 + a4) / (2 v f), u3 = w^2 f - 2 u,
                v3 = w (u - u3) - v
*/
static int
tors_add(tors_struct * D, const tors_struct * A, const tors_struct * B,
        tors_ctx_struct * C)
{
    gr_poly_t w, num, den, u3, v3, tmp;
    int status = GR_SUCCESS;

    if (A->is_inf)
        return tors_set(D, B, C);

    if (B->is_inf)
        return tors_set(D, A, C);

    gr_poly_init(w, C->R);
    gr_poly_init(num, C->R);
    gr_poly_init(den, C->R);
    gr_poly_init(u3, C->R);
    gr_poly_init(v3, C->R);
    gr_poly_init(tmp, C->R);

    if (gr_poly_equal(A->u, B->u, C->R) == T_TRUE)
    {
        if (gr_poly_equal(A->v, B->v, C->R) != T_TRUE)
        {
            /* B = -A */
            D->is_inf = 1;
            goto cleanup;
        }

        /* doubling */
        status |= tors_mulmod(tmp, A->u, A->u, C);
        status |= gr_poly_mul_ui(num, tmp, 3, C->R);
        status |= gr_poly_add(num, num, C->a4, C->R);

        status |= tors_mulmod(den, A->v, C->f, C);
        status |= gr_poly_mul_ui(den, den, 2, C->R);
    }
    else
    {
        status |= gr_poly_sub(num, B->v, A->v, C->R);
        status |= gr_poly_sub(den, B->u, A->u, C->R);
    }

    status |= tors_inv(tmp, den, C);

    if (status != GR_SUCCESS)
        goto cleanup;

    status |= tors_mulmod(w, num, tmp, C);

    /* u3 = w^2 f - u1 - u2 */
    status |= tors_mulmod(u3, w, w, C);
    status |= tors_mulmod(u3, u3, C->f, C);
    status |= gr_poly_sub(u3, u3, A->u, C->R);
    status |= gr_poly_sub(u3, u3, B->u, C->R);

    /* v3 = w (u1 - u3) - v1 */
    status |= gr_poly_sub(tmp, A->u, u3, C->R);
    status |= tors_mulmod(v3, w, tmp, C);
    status |= gr_poly_sub(v3, v3, A->v, C->R);

    if (status == GR_SUCCESS)
    {
        status |= gr_poly_set(D->u, u3, C->R);
        status |= gr_poly_set(D->v, v3, C->R);
        D->is_inf = 0;
    }

cleanup:
    gr_poly_clear(w, C->R);
    gr_poly_clear(num, C->R);
    gr_poly_clear(den, C->R);
    gr_poly_clear(u3, C->R);
    gr_poly_clear(v3, C->R);
    gr_poly_clear(tmp, C->R);

    return status;
}

/* k P by a left-to-right binary ladder; k is small (below l) */
static int
tors_mul_ui(tors_struct * D, const tors_struct * S, ulong k, tors_ctx_struct * C)
{
    tors_struct acc;
    int status = GR_SUCCESS;
    slong i;

    tors_init(&acc, C);

    if (k != 0)
    {
        for (i = FLINT_BIT_COUNT(k) - 1; i >= 0 && status == GR_SUCCESS; i--)
        {
            status |= tors_add(&acc, &acc, &acc, C);

            if ((k >> i) & 1)
                status |= tors_add(&acc, &acc, S, C);
        }
    }

    if (status == GR_SUCCESS)
        status = tors_set(D, &acc, C);

    tors_clear(&acc, C);

    return status;
}

/*
    t mod l for an odd prime l, via phi^2 - t phi + q = 0, working modulo
    "modulus" (which is psi_l, or a factor of it that a previous attempt
    ran into). On GR_UNABLE, *factor_out is set when a usable proper factor
    turned up.
*/
static int
_schoof_trace_mod_l_mod(ulong * tl, ulong l, const fmpz_t q, gr_ec_ctx_t ctx,
        const gr_poly_t modulus, gr_poly_t factor_out, int * got_factor)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    tors_ctx_struct C;
    tors_struct P, phiP, phi2P, qP, lhs, rhs;
    gr_poly_t xq, yq, tmp, xpoly;
    fmpz_t e;
    ulong tbar;
    int status = GR_SUCCESS, done = 0;

    C.R = R;
    C.have_factor = 0;
    gr_poly_init(C.psi, R);
    gr_poly_init(C.factor, R);
    gr_poly_init(xq, R);
    gr_poly_init(yq, R);
    gr_poly_init(tmp, R);
    gr_poly_init(xpoly, R);
    gr_poly_init(C.f, R);
    gr_poly_init(C.a4, R);
    gr_poly_preinv_init(C.P, R);
    fmpz_init(e);

    status |= gr_poly_set(C.psi, modulus, R);

    if (status != GR_SUCCESS || gr_poly_length(C.psi, R) < 2)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    status |= gr_poly_preinv_set(C.P, C.psi, R);

    /* f = x^3 + a4 x + a6, and the constant polynomial a4 */
    status |= gr_poly_zero(C.f, R);
    status |= gr_poly_set_coeff_si(C.f, 3, 1, R);
    status |= gr_poly_set_coeff_scalar(C.f, 1, GR_EC_A4(ctx), R);
    status |= gr_poly_set_coeff_scalar(C.f, 0, GR_EC_A6(ctx), R);
    status |= gr_poly_preinv_rem(C.f, C.f, C.P, R);

    status |= gr_poly_zero(C.a4, R);
    status |= gr_poly_set_coeff_scalar(C.a4, 0, GR_EC_A4(ctx), R);

    /* x as a polynomial, the generic point */
    status |= gr_poly_zero(xpoly, R);
    status |= gr_poly_set_coeff_si(xpoly, 1, 1, R);

    /* xq = x^q, yq = f^((q-1)/2), both mod psi_l */
    status |= gr_poly_preinv_powmod_x_fmpz(xq, q, C.P, R);

    fmpz_sub_ui(e, q, 1);
    fmpz_fdiv_q_ui(e, e, 2);
    status |= gr_poly_preinv_powmod_fmpz_binexp(yq, C.f, e, C.P, R);

    if (status != GR_SUCCESS)
        goto cleanup;

    tors_init(&P, &C);
    tors_init(&phiP, &C);
    tors_init(&phi2P, &C);
    tors_init(&qP, &C);
    tors_init(&lhs, &C);
    tors_init(&rhs, &C);

    /*
        P = (x, y). The reduction matters: once the descent has taken the
        modulus down to a linear factor, x is no longer reduced, and two
        points with the same x-coordinate would otherwise compare unequal
        and send the group law down the "distinct" branch, where it would
        invert zero.
    */
    status |= gr_poly_set(P.u, xpoly, R);
    status |= gr_poly_preinv_rem(P.u, P.u, C.P, R);
    status |= gr_poly_one(P.v, R);
    status |= gr_poly_preinv_rem(P.v, P.v, C.P, R);
    P.is_inf = 0;

    /* phi(P) = (xq, yq y) */
    status |= gr_poly_set(phiP.u, xq, R);
    status |= gr_poly_set(phiP.v, yq, R);
    phiP.is_inf = 0;

    /* phi^2(P) = (xq o xq, (yq o xq) yq y) */
    status |= gr_poly_preinv_compose_mod(phi2P.u, xq, xq, C.P, R);
    status |= gr_poly_preinv_compose_mod(tmp, yq, xq, C.P, R);
    status |= tors_mulmod(phi2P.v, tmp, yq, &C);
    phi2P.is_inf = 0;

    /* (q mod l) P */
    status |= tors_mul_ui(&qP, &P, fmpz_fdiv_ui(q, l), &C);

    /* lhs = phi^2(P) + (q mod l) P */
    status |= tors_add(&lhs, &phi2P, &qP, &C);

    if (status != GR_SUCCESS)
        goto cleanup_points;

    if (lhs.is_inf)
    {
        *tl = 0;
        done = 1;
    }

    for (tbar = 1; tbar <= (l - 1) / 2 && !done && status == GR_SUCCESS; tbar++)
    {
        status |= tors_mul_ui(&rhs, &phiP, tbar, &C);

        if (status != GR_SUCCESS)
            break;

        if (tors_equal(&lhs, &rhs, &C) == T_TRUE)
        {
            *tl = tbar;
            done = 1;
        }
        else
        {
            status |= tors_neg(&rhs, &rhs, &C);

            if (status == GR_SUCCESS && tors_equal(&lhs, &rhs, &C) == T_TRUE)
            {
                *tl = l - tbar;
                done = 1;
            }
        }
    }

    if (status == GR_SUCCESS && !done)
        status = GR_UNABLE;

cleanup_points:
    tors_clear(&rhs, &C);
    tors_clear(&lhs, &C);
    tors_clear(&qP, &C);
    tors_clear(&phi2P, &C);
    tors_clear(&phiP, &C);
    tors_clear(&P, &C);

cleanup:
    *got_factor = 0;

    if (status != GR_SUCCESS && C.have_factor)
        if (gr_poly_set(factor_out, C.factor, R) == GR_SUCCESS)
            *got_factor = 1;

    fmpz_clear(e);
    gr_poly_clear(C.factor, R);
    gr_poly_preinv_clear(C.P, R);
    gr_poly_clear(C.a4, R);
    gr_poly_clear(C.f, R);
    gr_poly_clear(xpoly, R);
    gr_poly_clear(tmp, R);
    gr_poly_clear(yq, R);
    gr_poly_clear(xq, R);
    gr_poly_clear(C.psi, R);

    return status;
}

/* start from psi_l, and descend into a factor whenever one is offered */
static int
_schoof_trace_mod_l(ulong * tl, ulong l, const fmpz_t q, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_poly_t mod, factor;
    slong depth;
    int status, got_factor;

    gr_poly_init(mod, R);
    gr_poly_init(factor, R);

    status = gr_ec_ctx_division_poly(mod, l, ctx);

    for (depth = 0; status == GR_SUCCESS && depth < 16; depth++)
    {
        status = _schoof_trace_mod_l_mod(tl, l, q, ctx, mod, factor, &got_factor);

        if (status == GR_SUCCESS || !got_factor)
            break;

        status = gr_poly_set(mod, factor, R);
    }

    gr_poly_clear(mod, R);
    gr_poly_clear(factor, R);

    return status;
}

/* t mod 2: E has a rational point of order 2 iff gcd(x^q - x, f) != 1 */
static int
_schoof_trace_mod_2(ulong * t2, const fmpz_t q, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_poly_preinv_t P;
    gr_poly_t f, xq, xpoly, g;
    int status = GR_SUCCESS;

    gr_poly_init(f, R);
    gr_poly_init(xq, R);
    gr_poly_init(xpoly, R);
    gr_poly_init(g, R);
    gr_poly_preinv_init(P, R);

    status |= gr_poly_zero(f, R);
    status |= gr_poly_set_coeff_si(f, 3, 1, R);
    status |= gr_poly_set_coeff_scalar(f, 1, GR_EC_A4(ctx), R);
    status |= gr_poly_set_coeff_scalar(f, 0, GR_EC_A6(ctx), R);

    status |= gr_poly_preinv_set(P, f, R);
    status |= gr_poly_preinv_powmod_x_fmpz(xq, q, P, R);

    status |= gr_poly_zero(xpoly, R);
    status |= gr_poly_set_coeff_si(xpoly, 1, 1, R);
    status |= gr_poly_sub(xq, xq, xpoly, R);
    status |= gr_poly_gcd(g, xq, f, R);

    if (status == GR_SUCCESS)
        *t2 = (gr_poly_length(g, R) > 1) ? 0 : 1;

    gr_poly_preinv_clear(P, R);
    gr_poly_clear(g, R);
    gr_poly_clear(xpoly, R);
    gr_poly_clear(xq, R);
    gr_poly_clear(f, R);

    return status;
}

/* ------------------------------------------------------------------ */

int
gr_ec_ctx_cardinality_schoof(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_t q, p, bound, M, t, sq;
    ulong l, tl;
    int status = GR_SUCCESS;

    if (gr_ctx_is_field(R) != T_TRUE)
        return GR_DOMAIN;

    /* the simple version only knows the short model away from 2 and 3 */
    if (gr_ec_ctx_model(ctx) != GR_EC_SHORT_WEIERSTRASS)
        return GR_DOMAIN;

    fmpz_init(q); fmpz_init(p); fmpz_init(bound);
    fmpz_init(M); fmpz_init(t); fmpz_init(sq);

    if (gr_ctx_cardinality_fmpz(q, R) != GR_SUCCESS
            || gr_ctx_fq_prime(p, R) != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    if (fmpz_cmp_ui(p, 3) <= 0)
    {
        status = GR_DOMAIN;
        goto cleanup;
    }

    /* the primes must multiply past 4 sqrt(q) to pin t down */
    fmpz_sqrt(sq, q);
    fmpz_mul_ui(bound, sq, 4);
    fmpz_add_ui(bound, bound, 4);

    fmpz_zero(t);
    fmpz_one(M);

    status |= _schoof_trace_mod_2(&tl, q, ctx);

    if (status != GR_SUCCESS)
        goto cleanup;

    fmpz_CRT_ui(t, t, M, tl, 2, 1);
    fmpz_mul_ui(M, M, 2);

    for (l = 3; fmpz_cmp(M, bound) <= 0 && status == GR_SUCCESS; l = n_nextprime(l, 1))
    {
        if (fmpz_cmp_ui(p, l) == 0)
            continue;

        status = _schoof_trace_mod_l(&tl, l, q, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        fmpz_CRT_ui(t, t, M, tl, l, 1);
        fmpz_mul_ui(M, M, l);
    }

    if (status == GR_SUCCESS)
    {
        /* #E = q + 1 - t */
        fmpz_add_ui(res, q, 1);
        fmpz_sub(res, res, t);
    }

cleanup:
    fmpz_clear(q); fmpz_clear(p); fmpz_clear(bound);
    fmpz_clear(M); fmpz_clear(t); fmpz_clear(sq);

    return status;
}

/* ------------------------------------------------------------------ */
/* complex multiplication and supersingular curves                    */
/* ------------------------------------------------------------------ */

/*
    When E/F_p has complex multiplication by an order of small discriminant
    the trace is not something to search for. For the thirteen discriminants
    D of class number one the Hilbert class polynomial is linear, so
    j(E) = j_D is exactly the statement that E has that complex
    multiplication, and then 4p = t^2 + |D| v^2 has a solution, which
    qfb_cornacchia finds. Supersingular curves are the same story with
    t = 0, which is what happens when D is not a square modulo p.

    Cornacchia only gives |t|, and which of the twists with this
    j-invariant the curve actually is has to be decided separately. None of
    the thirteen cases needs a scalar multiplication for that.

    For j = 0 and j = 1728 the curve has sextic and quartic twists, and the
    sextic (respectively quartic) residue character of the coefficient
    picks one out, at the cost of one exponentiation. Those characters are
    classical evaluations of Jacobi sums -- see Ireland and Rosen, A
    Classical Introduction to Modern Number Theory, chapter 18.

    For the other eleven only the quadratic twist exists, and then two
    independent bits decide the sign, neither of which costs more than a
    Jacobi symbol:

      * a sign convention for |t| itself. For odd m = |D| the Jacobi symbol
        (t | m) does it, because -1 is a non-residue modulo every such m,
        so exactly one of +-t satisfies (t | m) = 1. Where m is a power of
        two (D = -8 and D = -16) a congruence on t and v takes its place.

      * which quadratic twist the curve is, which is the Legendre symbol
        (c_D a6 | p) for a constant c_D attached to the discriminant.

    That this is the shape of the answer is classical; see Rubin and
    Silverberg, Choosing the correct elliptic curve in the CM method,
    Math. Comp. 79 (2010). The constants c_D and the overall signs in
    cm_discs below were determined here by calibration against
    gr_ec_ctx_cardinality_bsgs over three thousand curves per discriminant,
    and each resulting rule was then confirmed on twenty thousand further
    curves.
*/

/* the representative of x modulo p in (-p/2, p/2] */
static void
_center(fmpz_t x, const fmpz_t p)
{
    fmpz_t h;
    fmpz_init(h);
    fmpz_tdiv_q_2exp(h, p, 1);

    if (fmpz_cmp(x, h) > 0)
        fmpz_sub(x, x, p);

    fmpz_clear(h);
}

/*
    The trace of y^2 = x^3 + a6 over F_p, which has j = 0.

    For p = 2 mod 3 the curve is supersingular. Otherwise 4p = A^2 + 27 B^2
    has a solution, normalised by A = 2 mod 3, and the sextic character of
    -108 a6 turns |A| into the trace of this particular sextic twist.
*/
static int
_cm_trace_j0(fmpz_t t, const fmpz_t a6, const fmpz_t p)
{
    fmpz_t D, s, A, B, d, e;
    int ok = 0;

    if (fmpz_fdiv_ui(p, 3) != 1)
        return fmpz_zero(t), 1;                 /* supersingular */

    fmpz_init(D); fmpz_init(s); fmpz_init(A);
    fmpz_init(B); fmpz_init(d); fmpz_init(e);

    fmpz_set_si(D, -27);
    fmpz_mod(D, D, p);

    if (fmpz_jacobi(D, p) == 1 && fmpz_sqrtmod(s, D, p)
            && qfb_cornacchia(A, B, p, -27, s))
    {
        if (fmpz_fdiv_ui(A, 3) == 1)
            fmpz_neg(A, A);

        fmpz_mul_si(d, a6, -108);
        fmpz_mod(d, d, p);

        fmpz_sub_ui(e, p, 1);
        fmpz_divexact_ui(e, e, 6);
        fmpz_powm(d, d, e, p);

        fmpz_mul(t, A, d);
        fmpz_mod(t, t, p);
        _center(t, p);
        ok = 1;
    }

    fmpz_clear(D); fmpz_clear(s); fmpz_clear(A);
    fmpz_clear(B); fmpz_clear(d); fmpz_clear(e);

    return ok;
}

/*
    The trace of y^2 = x^3 + a4 x over F_p, which has j = 1728.

    For p = 3 mod 4 the curve is supersingular. Otherwise 4p = A^2 + 4 B^2,
    normalised to the even member of the pair with A = 2 mod 8, and the
    quartic character of a4 selects the quartic twist.
*/
static int
_cm_trace_j1728(fmpz_t t, const fmpz_t a4, const fmpz_t p)
{
    fmpz_t D, s, A, B, e;
    int ok = 0;

    if (fmpz_fdiv_ui(p, 4) != 1)
        return fmpz_zero(t), 1;                 /* supersingular */

    fmpz_init(D); fmpz_init(s); fmpz_init(A); fmpz_init(B); fmpz_init(e);

    fmpz_set_si(D, -4);
    fmpz_mod(D, D, p);

    if (fmpz_jacobi(D, p) == 1 && fmpz_sqrtmod(s, D, p)
            && qfb_cornacchia(A, B, p, -4, s))
    {
        if (fmpz_fdiv_ui(A, 4) == 0)
            fmpz_set(A, B);

        if (fmpz_is_odd(A))
            fmpz_mul_2exp(A, A, 1);

        if (fmpz_fdiv_ui(A, 8) == 6)
            fmpz_neg(A, A);

        fmpz_sub_ui(e, p, 1);
        fmpz_tdiv_q_2exp(e, e, 2);
        fmpz_powm(e, a4, e, p);

        fmpz_mul(t, A, e);
        fmpz_mod(t, t, p);
        _center(t, p);
        ok = 1;
    }

    fmpz_clear(D); fmpz_clear(s); fmpz_clear(A); fmpz_clear(B); fmpz_clear(e);

    return ok;
}

/* how the sign convention for |t| is pinned down */
#define CM_SIGN_JACOBI 0        /* by (t | m), m odd */
#define CM_SIGN_D8     1        /* D = -8:  by t mod 16 and v mod 4 */
#define CM_SIGN_D16    2        /* D = -16: by t mod 8 */

/*
    The orders of discriminant D with class number one other than -3 and
    -4, the j-invariant of the corresponding curve, and the data of the
    rule above:

      Dc   the discriminant Cornacchia is run on. It is D itself except
           for -16 and -28, whose own forms do not represent every split
           p, so the maximal order is used and the normalisation picks the
           right member of the pair out.
      m    modulus of the Jacobi symbol on t, 1 when the kind says
           otherwise. For -27 the symbol modulo 27 is the symbol modulo 3.
      c    the constant in the Legendre symbol (c a6 | p).
      eps  the overall sign.
*/
typedef struct
{
    slong D;
    slong j;
    slong Dc;
    slong m;
    slong c;
    slong eps;
    int kind;
}
cm_disc_struct;

static const cm_disc_struct cm_discs[] = {
    {  -7, WORD(-3375),                  -7,   7, WORD(-2),      1, CM_SIGN_JACOBI },
    {  -8, WORD(8000),                   -8,   1, WORD(21),      1, CM_SIGN_D8     },
    { -11, WORD(-32768),                -11,  11, WORD(21),     -1, CM_SIGN_JACOBI },
    { -12, WORD(54000),                 -12,   3, WORD(22),     -1, CM_SIGN_JACOBI },
    { -16, WORD(287496),                 -4,   1, WORD(7),       1, CM_SIGN_D16    },
    { -19, WORD(-884736),               -19,  19, WORD(1),      -1, CM_SIGN_JACOBI },
    { -27, WORD(-12288000),             -27,   3, WORD(253),    -1, CM_SIGN_JACOBI },
    { -28, WORD(16581375),               -7,   7, WORD(-114),    1, CM_SIGN_JACOBI },
    { -43, WORD(-884736000),            -43,  43, WORD(21),     -1, CM_SIGN_JACOBI },
    { -67, WORD(-147197952000),         -67,  67, WORD(217),    -1, CM_SIGN_JACOBI },
    {-163, WORD(-262537412640768000),  -163, 163, WORD(185801), -1, CM_SIGN_JACOBI },
};

#define CM_NUM_DISCS (sizeof(cm_discs) / sizeof(cm_disc_struct))

/*
    The trace of a curve with j = j_D for one of the entries above, given
    its a6. Returns 0 if the trace could not be determined, which for a
    curve that really does have this j-invariant should not happen.
*/
static int
_cm_trace_disc(fmpz_t t, const cm_disc_struct * E, const fmpz_t a6,
        const fmpz_t p)
{
    fmpz_t D, s, a, b, u;
    slong sgn = E->eps;
    int ok = 0;

    fmpz_init(D); fmpz_init(s); fmpz_init(a); fmpz_init(b); fmpz_init(u);

    fmpz_set_si(D, E->Dc);
    fmpz_mod(D, D, p);

    if (fmpz_jacobi(D, p) != 1)
    {
        fmpz_zero(t);           /* inert: the reduction is supersingular */
        ok = 1;
        goto cleanup;
    }

    if (!fmpz_sqrtmod(s, D, p) || !qfb_cornacchia(a, b, p, E->Dc, s))
        goto cleanup;

    if (E->kind == CM_SIGN_JACOBI)
    {
        int e = n_jacobi((slong) fmpz_fdiv_ui(a, (ulong) E->m), (ulong) E->m);

        if (e == 0)
            goto cleanup;

        if (e < 0)
            sgn = -sgn;
    }
    else if (E->kind == CM_SIGN_D8)
    {
        /* here 4p = a^2 + 8 b^2 with a = 2 mod 4 */
        if (fmpz_fdiv_ui(a, 16) >= 8)
            sgn = -sgn;

        if (fmpz_fdiv_ui(b, 4) == 2)
            sgn = -sgn;
    }
    else
    {
        /* here 4p = a^2 + 4 b^2 and exactly one of a/2, b is odd; the
           trace is twice that one, normalised to 2 mod 4 */
        if (fmpz_fdiv_ui(a, 4) == 0)
            fmpz_swap(a, b);

        if (fmpz_is_odd(a))
            fmpz_mul_2exp(a, a, 1);

        if (fmpz_fdiv_ui(a, 8) == 6)
            sgn = -sgn;
    }

    fmpz_mul_si(u, a6, E->c);
    fmpz_mod(u, u, p);

    if (fmpz_is_zero(u))
        goto cleanup;

    if (fmpz_jacobi(u, p) < 0)
        sgn = -sgn;

    if (sgn > 0)
        fmpz_set(t, a);
    else
        fmpz_neg(t, a);

    ok = 1;

cleanup:
    fmpz_clear(D); fmpz_clear(s); fmpz_clear(a);
    fmpz_clear(b); fmpz_clear(u);

    return ok;
}

/*
    Is j(E) = J, without dividing? For y^2 = x^3 + a4 x + a6 the
    j-invariant is 6912 a4^3 / (4 a4^3 + 27 a6^2), so the question is
    whether 6912 a4^3 = J (4 a4^3 + 27 a6^2) in F_p. The numerator is the
    same for every J in the table, so the caller passes it reduced.
*/
static int
_j_equals(slong J, const fmpz_t num, const fmpz_t den, const fmpz_t p)
{
    fmpz_t v;
    int eq;

    fmpz_init(v);

    fmpz_mul_si(v, den, J);
    fmpz_mod(v, v, p);

    eq = fmpz_equal(num, v);

    fmpz_clear(v);

    return eq;
}

int
gr_ec_ctx_cardinality_cm(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_t p, a4, a6, a4c, den, t;
    slong deg, i;
    int status = GR_SUCCESS, found = 0;

    if (gr_ctx_is_field(R) != T_TRUE)
        return GR_DOMAIN;

    /* the class number one theory used here is over a prime field, and the
       twist characters below are written for the short model */
    if (gr_ctx_fq_degree(&deg, R) != GR_SUCCESS || deg != 1
            || gr_ec_ctx_model(ctx) != GR_EC_SHORT_WEIERSTRASS)
        return GR_UNABLE;

    fmpz_init(p); fmpz_init(a4); fmpz_init(a6);
    fmpz_init(a4c); fmpz_init(den); fmpz_init(t);

    if (gr_ctx_fq_prime(p, R) != GR_SUCCESS || fmpz_cmp_ui(p, 3) <= 0
            || gr_get_fmpz(a4, GR_EC_A4(ctx), R) != GR_SUCCESS
            || gr_get_fmpz(a6, GR_EC_A6(ctx), R) != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    if (fmpz_is_zero(a4))
        found = _cm_trace_j0(t, a6, p);
    else if (fmpz_is_zero(a6))
        found = _cm_trace_j1728(t, a4, p);
    else
    {
        /* the numerator 6912 a4^3 and the denominator 4 a4^3 + 27 a6^2 of
           the j-invariant, both reduced once for the whole table */
        {
            fmpz_t w;
            fmpz_init(w);

            fmpz_powm_ui(a4c, a4, 3, p);
            fmpz_mul_ui(den, a4c, 4);
            fmpz_mul(w, a6, a6);
            fmpz_mul_ui(w, w, 27);
            fmpz_add(den, den, w);
            fmpz_mod(den, den, p);

            fmpz_mul_ui(a4c, a4c, 6912);
            fmpz_mod(a4c, a4c, p);

            fmpz_clear(w);
        }

        for (i = 0; i < (slong) CM_NUM_DISCS && !found; i++)
            if (_j_equals(cm_discs[i].j, a4c, den, p))
                found = _cm_trace_disc(t, cm_discs + i, a6, p);
    }

    if (found)
    {
        fmpz_add_ui(res, p, 1);
        fmpz_sub(res, res, t);
    }
    else
        status = GR_UNABLE;

cleanup:
    fmpz_clear(p); fmpz_clear(a4); fmpz_clear(a6);
    fmpz_clear(a4c); fmpz_clear(den); fmpz_clear(t);

    return status;
}

/* ------------------------------------------------------------------ */
/* dispatch                                                           */
/* ------------------------------------------------------------------ */

/*
    Below this, walking the whole field is both quicker than setting up a
    baby-step table and free of any randomness.
*/
#define GR_EC_CARDINALITY_NAIVE_CUTOFF WORD(10000)

/*
    BSGS costs O(q^(1/4)) group operations and Schoof is polynomial in
    log q, so Schoof eventually wins; measured against each other (see
    profile/p-cardinality.c) the crossover sits near 2^80, with BSGS ahead
    below it. BSGS also needs square roots in the base ring, to find a
    point at all, so where gr_sqrt is unavailable -- mpn_mod, at the time
    of writing -- Schoof is the only one of the two that can run, and the
    order below makes that fall out on its own.
*/
#define GR_EC_CARDINALITY_SCHOOF_BITS 80

int
gr_ec_ctx_cardinality(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_t q;
    int status, schoof_first;

    if (gr_ctx_is_field(R) != T_TRUE)
        return GR_DOMAIN;

    fmpz_init(q);

    if (gr_ctx_cardinality_fmpz(q, R) != GR_SUCCESS)
    {
        fmpz_clear(q);
        return GR_UNABLE;
    }

    if (fmpz_cmp_si(q, GR_EC_CARDINALITY_NAIVE_CUTOFF) <= 0)
    {
        status = gr_ec_ctx_cardinality_naive(res, ctx);

        if (status == GR_SUCCESS)
        {
            fmpz_clear(q);
            return status;
        }
    }

    /*
        Complex multiplication first: it is a j-invariant comparison and a
        handful of scalar multiplications, and it gives up at once when the
        curve is not one of the special ones, so it costs almost nothing to
        try and saves everything when it applies.
    */
    if (gr_ec_ctx_cardinality_cm(res, ctx) == GR_SUCCESS)
    {
        fmpz_clear(q);
        return GR_SUCCESS;
    }

    schoof_first = (fmpz_bits(q) >= GR_EC_CARDINALITY_SCHOOF_BITS);

    status = schoof_first ? gr_ec_ctx_cardinality_schoof(res, ctx)
                          : gr_ec_ctx_cardinality_bsgs(res, ctx);

    if (status != GR_SUCCESS)
        status = schoof_first ? gr_ec_ctx_cardinality_bsgs(res, ctx)
                              : gr_ec_ctx_cardinality_schoof(res, ctx);

    /* a curve whose group has small exponent can leave BSGS with more than
       one candidate; fall back while walking the field is still affordable */
    if (status != GR_SUCCESS && fmpz_cmp_si(q, GR_EC_NAIVE_MAX_Q) <= 0)
        status = gr_ec_ctx_cardinality_naive(res, ctx);

    fmpz_clear(q);

    return status;
}
