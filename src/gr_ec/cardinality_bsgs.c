/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpz.h"
#include "fmpz_factor.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"
#include "fq_zech.h"
#include <math.h>
#include "ulong_extras.h"
#include "gr.h"
#include "gr_ec.h"
#include "impl.h"

/* ------------------------------------------------------------------ */
/* baby-step giant-step (Shanks, Mestre)                              */
/* ------------------------------------------------------------------ */

/*
    Points have no canonical ordering in a general finite field, so the
    baby-step table is keyed on a hash of the decimal (or polynomial)
    string of the x-coordinate (or of its integer value in a prime
    field). A hit is only a candidate; it is always
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
    fmpz_t z;

    /* the finite field types, read directly */
    if (R->which_ring == GR_CTX_FQ_ZECH)
    {
        h = ((const fq_zech_struct *) x)->value * UWORD(0x9e3779b97f4a7c15);
        return h ^ (h >> 29);
    }

    if (R->which_ring == GR_CTX_FQ_NMOD)
    {
        const nmod_poly_struct * f = x;

        for (i = 0; i < f->length; i++)
            h = (h ^ f->coeffs[i]) * UWORD(0x100000001b3);

        return h ^ (h >> 29);
    }

    if (R->which_ring == GR_CTX_FQ)
    {
        const fmpz_poly_struct * f = x;

        for (i = 0; i < f->length; i++)
            h = (h ^ fmpz_fdiv_ui(f->coeffs + i, UWORD_MAX)) * UWORD(0x100000001b3);

        return h ^ (h >> 29);
    }

    /* an element of the prime field has an integer to hash directly */
    fmpz_init(z);

    if (gr_get_fmpz(z, x, R) == GR_SUCCESS)
    {
        h = fmpz_fdiv_ui(z, UWORD_MAX) * UWORD(0x9e3779b97f4a7c15);
        fmpz_clear(z);
        return h ^ (h >> 29);
    }

    fmpz_clear(z);

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


/* ------------------------------------------------------------------ */
/* walks in affine coordinates, many at once                          */
/* ------------------------------------------------------------------ */

/*
    A baby-step giant-step search spends nearly all its time stepping
    along P, P + D, P + 2D, ... and hashing each x-coordinate. One step in
    projective coordinates plus the inversion that the affine x needs
    costs some sixty multiplications in the base ring; running many walks
    side by side in affine coordinates, their chord slopes share a single
    inversion (Montgomery's trick), and a step costs about six.
*/

/* P[i] += D for i < n, with one inversion; t has room for 2n elements */
static int
_aff_add_batch(gr_ec_aff_point_struct * P, slong n, const gr_ec_aff_point_t D,
        gr_ptr t, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem, i;
    gr_ptr u = t, w = GR_ENTRY(t, n, sz);
    gr_ptr inv, lam, nu, x3;
    char * reg;
    int is_short = (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS);
    int status = GR_SUCCESS;

    if (n == 0 || gr_ec_aff_point_is_inf(D, ctx) == T_TRUE)
        return GR_SUCCESS;

    reg = flint_malloc(n);
    GR_TMP_INIT4(inv, lam, nu, x3, R);

    /* u_i = x_i - x_D where the chord formula applies, 1 elsewhere */
    for (i = 0; i < n; i++)
    {
        reg[i] = 0;

        if (P[i].is_infinity == T_FALSE)
        {
            status |= gr_sub(GR_ENTRY(u, i, sz), GR_EC_AFF_POINT_X(P + i, ctx),
                    GR_EC_AFF_POINT_X(D, ctx), R);
            reg[i] = (gr_is_zero(GR_ENTRY(u, i, sz), R) == T_FALSE);
        }

        if (!reg[i])
            status |= gr_one(GR_ENTRY(u, i, sz), R);

        if (i == 0)
            status |= gr_set(w, u, R);
        else
            status |= gr_mul(GR_ENTRY(w, i, sz), GR_ENTRY(w, i - 1, sz),
                    GR_ENTRY(u, i, sz), R);
    }

    status |= gr_inv(inv, GR_ENTRY(w, n - 1, sz), R);

    for (i = n - 1; i >= 0 && status == GR_SUCCESS; i--)
    {
        gr_ec_aff_point_struct * Pi = P + i;

        /* lam = 1/u_i, then inv = 1/(u_0 ... u_{i-1}) */
        if (i > 0)
            status |= gr_mul(lam, inv, GR_ENTRY(w, i - 1, sz), R);
        else
            status |= gr_set(lam, inv, R);

        status |= gr_mul(inv, inv, GR_ENTRY(u, i, sz), R);

        if (!reg[i])
        {
            status |= gr_ec_aff_point_add(Pi, Pi, D, ctx);
            continue;
        }

        /* lam = (y_i - y_D) / (x_i - x_D), nu = y_i - lam x_i */
        status |= gr_sub(nu, GR_EC_AFF_POINT_Y(Pi, ctx), GR_EC_AFF_POINT_Y(D, ctx), R);
        status |= gr_mul(lam, lam, nu, R);
        status |= gr_mul(nu, lam, GR_EC_AFF_POINT_X(Pi, ctx), R);
        status |= gr_sub(nu, GR_EC_AFF_POINT_Y(Pi, ctx), nu, R);

        /* x3 = lam^2 + a1 lam - a2 - x_i - x_D */
        status |= gr_sqr(x3, lam, R);
        if (!is_short)
        {
            status |= gr_addmul(x3, GR_EC_A1(ctx), lam, R);
            status |= gr_sub(x3, x3, GR_EC_A2(ctx), R);
        }
        status |= gr_sub(x3, x3, GR_EC_AFF_POINT_X(Pi, ctx), R);
        status |= gr_sub(x3, x3, GR_EC_AFF_POINT_X(D, ctx), R);

        /* y3 = -(lam + a1) x3 - nu - a3 = -(lam + a1) x3 - (nu + a3) */
        if (!is_short)
        {
            status |= gr_add(lam, lam, GR_EC_A1(ctx), R);
            status |= gr_add(nu, nu, GR_EC_A3(ctx), R);
        }
        status |= gr_mul(lam, lam, x3, R);
        status |= gr_add(lam, lam, nu, R);
        status |= gr_neg(GR_EC_AFF_POINT_Y(Pi, ctx), lam, R);
        gr_swap(GR_EC_AFF_POINT_X(Pi, ctx), x3, R);
    }

    GR_TMP_CLEAR4(inv, lam, nu, x3, R);
    flint_free(reg);

    return status;
}

/* affine copies of projective points, with one inversion */
static int
_aff_set_batch(gr_ec_aff_point_struct * A, const gr_ec_point_struct * P,
        slong n, gr_ptr t, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem, i;
    gr_ptr u = t, w = GR_ENTRY(t, n, sz);
    gr_ptr inv, zi;
    int status = GR_SUCCESS;

    if (n == 0)
        return GR_SUCCESS;

    GR_TMP_INIT2(inv, zi, R);

    for (i = 0; i < n; i++)
    {
        truth_t inf = gr_ec_point_is_inf(P + i, ctx);

        if (inf == T_UNKNOWN)
            status |= GR_UNABLE;

        A[i].is_infinity = (inf == T_TRUE) ? T_TRUE : T_FALSE;

        if (inf == T_FALSE)
            status |= gr_set(GR_ENTRY(u, i, sz), GR_EC_POINT_Z(P + i, ctx), R);
        else
            status |= gr_one(GR_ENTRY(u, i, sz), R);

        if (i == 0)
            status |= gr_set(w, u, R);
        else
            status |= gr_mul(GR_ENTRY(w, i, sz), GR_ENTRY(w, i - 1, sz),
                    GR_ENTRY(u, i, sz), R);
    }

    status |= gr_inv(inv, GR_ENTRY(w, n - 1, sz), R);

    for (i = n - 1; i >= 0 && status == GR_SUCCESS; i--)
    {
        if (i > 0)
            status |= gr_mul(zi, inv, GR_ENTRY(w, i - 1, sz), R);
        else
            status |= gr_set(zi, inv, R);

        status |= gr_mul(inv, inv, GR_ENTRY(u, i, sz), R);

        if (A[i].is_infinity == T_FALSE)
        {
            status |= gr_mul(GR_EC_AFF_POINT_X(A + i, ctx), GR_EC_POINT_X(P + i, ctx), zi, R);
            status |= gr_mul(GR_EC_AFF_POINT_Y(A + i, ctx), GR_EC_POINT_Y(P + i, ctx), zi, R);
        }
    }

    GR_TMP_CLEAR2(inv, zi, R);

    return status;
}

/* called on every point of a walk; returning nonzero stops the walk */
typedef int (*_walk_visit_t)(void * data, slong k, slong z, ulong h,
        const gr_ec_aff_point_struct * pt, int * status);

#define GR_EC_WALK_LANES 128

/*
    Visits, for k < nu and z < ztot, the points

        P0 + sgn U[k] Bu + z Dz

    (U sorted increasingly; U = NULL means nu = 1 and U[0] = 0), handing
    visit the hash of their x-coordinate. The walks along z are cut into
    lanes so that at least GR_EC_WALK_LANES of them step together, and
    the starting points are reached from each other by short scalar
    multiplications, the U being sorted.
*/
static int
_walk(const gr_ec_point_t P0, const fmpz * U, slong nu, int sgn,
        const gr_ec_point_t Bu, const gr_ec_point_t Dz, slong ztot,
        _walk_visit_t visit, void * data, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong K = GR_EC_WALK_LANES, S, zc, nl, c0, i, z, kcur = -1;
    gr_ec_point_struct * E, * start;
    gr_ec_aff_point_struct * lane;
    gr_ec_aff_point_t Da;
    gr_ec_point_t Acur, T;
    gr_ptr t;
    fmpz_t d;
    int status = GR_SUCCESS, stop = 0;

    if (nu <= 0 || ztot <= 0)
        return GR_SUCCESS;

    /* short walks: fewer lanes, each still sharing its inversion widely */
    K = FLINT_MAX(1, FLINT_MIN(K, (nu * ztot) / 16));

    /* S lanes per k, each zc steps long */
    S = (nu >= K) ? 1 : FLINT_MIN(ztot, (K + nu - 1) / nu);
    zc = (ztot + S - 1) / S;
    nl = nu * S;

    E = flint_malloc(S * sizeof(gr_ec_point_struct));
    for (i = 0; i < S; i++)
        gr_ec_point_init(E + i, ctx);

    start = flint_malloc(K * sizeof(gr_ec_point_struct));
    lane = flint_malloc(K * sizeof(gr_ec_aff_point_struct));
    for (i = 0; i < K; i++)
    {
        gr_ec_point_init(start + i, ctx);
        gr_ec_aff_point_init(lane + i, ctx);
    }

    gr_ec_aff_point_init(Da, ctx);
    gr_ec_point_init(Acur, ctx);
    gr_ec_point_init(T, ctx);
    t = gr_heap_init_vec(2 * K, R);
    fmpz_init(d);

    /* E[s] = s zc Dz */
    status |= gr_ec_point_zero(E + 0, ctx);
    if (S > 1)
        status |= gr_ec_point_mul_si(E + 1, Dz, zc, ctx);
    for (i = 2; i < S; i++)
        status |= gr_ec_point_add(E + i, E + i - 1, E + 1, ctx);

    status |= gr_ec_aff_point_set_point(Da, Dz, ctx);
    status |= gr_ec_point_set(Acur, P0, ctx);

    for (c0 = 0; c0 < nl && status == GR_SUCCESS && !stop; c0 += K)
    {
        slong n = FLINT_MIN(K, nl - c0);

        /* the starting points of lanes c0, ..., c0 + n - 1 */
        for (i = 0; i < n && status == GR_SUCCESS; i++)
        {
            slong k = (c0 + i) / S, s = (c0 + i) % S;

            if (k != kcur)
            {
                if (U != NULL)
                {
                    if (kcur < 0)
                        fmpz_set(d, U + k);
                    else
                        fmpz_sub(d, U + k, U + kcur);

                    if (sgn < 0)
                        fmpz_neg(d, d);

                    status |= gr_ec_point_mul_fmpz(T, Bu, d, ctx);
                    status |= gr_ec_point_add(Acur, Acur, T, ctx);
                }

                kcur = k;
            }

            status |= gr_ec_point_add(start + i, Acur, E + s, ctx);
        }

        status |= _aff_set_batch(lane, start, n, t, ctx);

        for (z = 0; z < zc && status == GR_SUCCESS && !stop; z++)
        {
            for (i = 0; i < n && !stop; i++)
            {
                slong k = (c0 + i) / S, zz = ((c0 + i) % S) * zc + z;
                ulong h;

                if (zz >= ztot)
                    continue;

                if (lane[i].is_infinity == T_TRUE)
                    h = UWORD(0x5bd1e9955bd1e995);
                else
                    h = _hash_elem(GR_EC_AFF_POINT_X(lane + i, ctx), R, &status);

                stop = visit(data, k, zz, h, lane + i, &status);
            }

            if (z + 1 < zc && !stop)
                status |= _aff_add_batch(lane, n, Da, t, ctx);
        }
    }

    fmpz_clear(d);
    gr_heap_clear_vec(t, 2 * K, R);
    gr_ec_point_clear(T, ctx);
    gr_ec_point_clear(Acur, ctx);
    gr_ec_aff_point_clear(Da, ctx);

    for (i = 0; i < K; i++)
    {
        gr_ec_point_clear(start + i, ctx);
        gr_ec_aff_point_clear(lane + i, ctx);
    }
    flint_free(lane);
    flint_free(start);

    for (i = 0; i < S; i++)
        gr_ec_point_clear(E + i, ctx);
    flint_free(E);

    return status;
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

/*
    Is there exactly one N in [lo, hi] with N = r mod M and N = 0 mod L?
    The true order is always one of them, so if there is exactly one it is
    the order.
*/
static int
_unique_in_interval(fmpz_t res, const fmpz_t r, const fmpz_t M,
        const fmpz_t L, const fmpz_t lo, const fmpz_t hi)
{
    fmpz_t g, Mg, Lg, k, N, mod;
    int unique = 0;

    fmpz_init(g); fmpz_init(Mg); fmpz_init(Lg);
    fmpz_init(k); fmpz_init(N); fmpz_init(mod);

    /* N = r + M k with M k = -r mod L, i.e. (M/g) k = -r/g mod L/g */
    fmpz_gcd(g, M, L);
    fmpz_fdiv_r(k, r, g);

    if (fmpz_is_zero(k))
    {
        fmpz_divexact(Mg, M, g);
        fmpz_divexact(Lg, L, g);
        fmpz_divexact(k, r, g);
        fmpz_neg(k, k);

        if (fmpz_is_one(Lg))
            fmpz_zero(k);
        else
        {
            fmpz_invmod(N, Mg, Lg);
            fmpz_mul(k, k, N);
            fmpz_mod(k, k, Lg);
        }

        fmpz_mul(N, M, k);
        fmpz_add(N, N, r);
        fmpz_mul(mod, M, Lg);

        /* the first one at or above lo, then is the next one past hi? */
        fmpz_sub(k, lo, N);
        fmpz_cdiv_q(k, k, mod);
        fmpz_addmul(N, k, mod);

        if (fmpz_cmp(N, hi) <= 0)
        {
            fmpz_add(k, N, mod);

            if (fmpz_cmp(k, hi) > 0)
            {
                fmpz_set(res, N);
                unique = 1;
            }
        }
    }

    fmpz_clear(g); fmpz_clear(Mg); fmpz_clear(Lg);
    fmpz_clear(k); fmpz_clear(N); fmpz_clear(mod);

    return unique;
}

typedef struct
{
    baby_struct * tab;
    slong ntab;
    int found;
    gr_ec_point_struct * P;
    gr_ec_point_struct * T;
    const fmpz * c;
    const fmpz * M;
    fmpz * n;
    fmpz * cand;
    slong klo, khi, m;
    gr_ec_ctx_struct * ctx;
}
_prog_struct;

/* the baby step j B = (z + 1) B */
static int
_prog_baby(void * data, slong FLINT_UNUSED(k), slong z, ulong h,
        const gr_ec_aff_point_struct * pt, int * FLINT_UNUSED(status))
{
    _prog_struct * D = data;

    if (pt->is_infinity == T_TRUE)
    {
        /* j B = O already: the order of P divides j M */
        fmpz_mul_si(D->n, D->M, z + 1);
        D->found = 1;
        return 1;
    }

    D->tab[D->ntab].hash = h;
    D->tab[D->ntab].j = z + 1;
    D->ntab++;

    return 0;
}

/* a giant step equal to Q - im B: is it +- j B for some baby j? */
static int
_prog_giant(_prog_struct * D, slong im, ulong h,
        const gr_ec_aff_point_struct * pt, int * status)
{
    gr_ec_ctx_struct * ctx = D->ctx;
    slong idx, nj, jj;
    fmpz_t t;
    int sj;

    if (pt->is_infinity == T_TRUE)
    {
        idx = -1;
        nj = 1;
    }
    else
    {
        idx = _baby_find(D->tab, D->ntab, h);

        if (idx < 0)
            return 0;

        nj = 0;
        while (idx + nj < D->ntab && D->tab[idx + nj].hash == h)
            nj++;
    }

    fmpz_init(t);

    /* Q - im B = +- j B, so k = im +- j; confirm N_k P = O */
    for (jj = 0; jj < nj && !D->found; jj++)
    {
        slong jv = (idx < 0) ? 0 : D->tab[idx + jj].j;

        for (sj = -1; sj <= 1 && !D->found; sj += 2)
        {
            slong k = im + sj * jv;

            if (k < D->klo || k > D->khi)
                continue;

            fmpz_set_si(t, k);
            fmpz_set(D->cand, D->c);
            fmpz_submul(D->cand, t, D->M);

            if (fmpz_sgn(D->cand) <= 0)
                continue;

            *status |= gr_ec_point_mul_fmpz(D->T, D->P, D->cand, ctx);

            if (*status == GR_SUCCESS && gr_ec_point_is_inf(D->T, ctx) == T_TRUE)
            {
                fmpz_set(D->n, D->cand);
                D->found = 1;
            }
        }
    }

    fmpz_clear(t);

    return D->found;
}

/* Q - z m B */
static int
_prog_giant_minus(void * data, slong FLINT_UNUSED(k), slong z, ulong h,
        const gr_ec_aff_point_struct * pt, int * status)
{
    _prog_struct * D = data;
    return _prog_giant(D, z * D->m, h, pt, status);
}

/* Q + (z + 1) m B */
static int
_prog_giant_plus(void * data, slong FLINT_UNUSED(k), slong z, ulong h,
        const gr_ec_aff_point_struct * pt, int * status)
{
    _prog_struct * D = data;
    return _prog_giant(D, -(z + 1) * D->m, h, pt, status);
}

/*
    Baby-step giant-step for #E once the trace is known to be t0 modulo
    M, which makes the candidates N = q + 1 - t the terms of a progression
    of difference M in the Hasse interval. M = 1 is the plain algorithm of
    Shanks and Mestre; after Schoof's or Elkies' steps M is large and the
    search takes the square root of what is left rather than walking it.

    With tc the representative of t0 nearest to zero, N_k = q + 1 - tc - k M
    kills P exactly when (q + 1 - tc) P = k (M P), which is a discrete
    logarithm in base B = M P with k in a known range; the baby steps are
    j B and the giant steps Q -+ i m B, matched on the x-coordinate so that
    each covers both signs of j.

    Any hit is confirmed by N P = O and only then used, as a multiple of
    the order of P: the exact orders of the points drawn go into their
    lcm L, and the answer is taken once exactly one N of the progression
    in the Hasse interval is divisible by L. The true order is always
    such an N, so the answer is proved, whatever hash collisions happen.

    GR_UNABLE, without doing anything, when more than max_baby baby steps
    would be needed (max_baby <= 0: no limit but a slong); and when no
    point resolves the question, which happens for groups whose exponent
    is small against the interval.
*/
int
_gr_ec_cardinality_bsgs_progression(fmpz_t res, const fmpz_t t0,
        const fmpz_t M, slong max_baby, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ec_point_t P, Q, B, mB, S, Sp, Sm, T;
    gr_ptr xc, yc;
    baby_struct * tab = NULL;
    flint_rand_t state;
    fmpz_t q, W, lo, hi, L, n, ord, t, tc, c, r, cand, kmin, kmax;
    slong m, attempt, ngiant, klo, khi;
    _prog_struct prog;
    int status = GR_SUCCESS, resolved = 0;

    if (gr_ctx_is_field(R) != T_TRUE)
        return GR_DOMAIN;

    if (fmpz_sgn(M) <= 0)
        return GR_DOMAIN;

    fmpz_init(q); fmpz_init(W); fmpz_init(lo); fmpz_init(hi);
    fmpz_init(L); fmpz_init(n); fmpz_init(ord); fmpz_init(t);
    fmpz_init(tc); fmpz_init(c); fmpz_init(r); fmpz_init(cand);
    fmpz_init(kmin); fmpz_init(kmax);

    if (gr_ctx_cardinality_fmpz(q, R) != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup_ints;
    }

    /* |t| <= W: 2 sqrt(q) rounded up, a little slack is free */
    fmpz_sqrt(W, q);
    fmpz_mul_ui(W, W, 2);
    fmpz_add_ui(W, W, 2);

    fmpz_add_ui(lo, q, 1);
    fmpz_sub(lo, lo, W);
    fmpz_add_ui(hi, q, 1);
    fmpz_add(hi, hi, W);

    /* tc = t0 mod M, balanced; c = q + 1 - tc; r = c mod M */
    fmpz_smod(tc, t0, M);
    fmpz_add_ui(c, q, 1);
    fmpz_sub(c, c, tc);
    fmpz_mod(r, c, M);

    /* t = tc + k M with |t| <= W */
    fmpz_neg(kmin, W);
    fmpz_sub(kmin, kmin, tc);
    fmpz_cdiv_q(kmin, kmin, M);
    fmpz_sub(kmax, W, tc);
    fmpz_fdiv_q(kmax, kmax, M);

    if (!fmpz_fits_si(kmin) || !fmpz_fits_si(kmax))
    {
        status = GR_UNABLE;
        goto cleanup_ints;
    }

    klo = fmpz_get_si(kmin);
    khi = fmpz_get_si(kmax);

    if (klo > khi)
    {
        status = GR_UNABLE;             /* t0 mod M contradicts Hasse */
        goto cleanup_ints;
    }

    /* baby steps j < m, giant steps |i| <= K/m + 1 with K = max |k| */
    {
        ulong K = FLINT_MAX(FLINT_ABS(klo), FLINT_ABS(khi));

        m = n_sqrt(K) + 1;

        if (max_baby > 0 && m > max_baby)
        {
            status = GR_UNABLE;
            goto cleanup_ints;
        }

        ngiant = K / m + 1;
    }

    flint_rand_init(state);
    flint_rand_set_seed(state, UWORD(0x9e3779b97f4a7c15), UWORD(0xbf58476d1ce4e5b9));

    gr_ec_point_init(P, ctx);
    gr_ec_point_init(Q, ctx);
    gr_ec_point_init(B, ctx);
    gr_ec_point_init(mB, ctx);
    gr_ec_point_init(S, ctx);
    gr_ec_point_init(Sp, ctx);
    gr_ec_point_init(Sm, ctx);
    gr_ec_point_init(T, ctx);
    GR_TMP_INIT2(xc, yc, R);

    tab = flint_malloc(m * sizeof(baby_struct));

    prog.P = P;
    prog.T = T;
    prog.c = c;
    prog.M = (fmpz *) M;
    prog.n = n;
    prog.cand = cand;
    prog.klo = klo;
    prog.khi = khi;
    prog.m = m;
    prog.ctx = ctx;

    fmpz_one(L);

    for (attempt = 0; attempt < GR_EC_BSGS_MAX_POINTS && !resolved; attempt++)
    {
        int found = 0;

        if (gr_ec_point_randtest(P, state, ctx) != GR_SUCCESS)
            continue;

        if (gr_ec_point_is_inf(P, ctx) != T_FALSE)
            continue;

        status |= gr_ec_point_mul_fmpz(B, P, M, ctx);

        if (status != GR_SUCCESS)
            break;

        /* baby steps: the x-coordinates of j B for j = 1, ..., m - 1 */
        prog.tab = tab;
        prog.ntab = 0;
        prog.found = 0;
        status |= _walk(B, NULL, 1, 1, B, B, m - 1, _prog_baby, &prog, ctx);

        if (status != GR_SUCCESS)
            break;

        if (!prog.found)
        {
            qsort(tab, prog.ntab, sizeof(baby_struct), _baby_cmp);

            fmpz_set_si(t, m);
            status |= gr_ec_point_mul_fmpz(Q, P, c, ctx);
            status |= gr_ec_point_mul_fmpz(mB, B, t, ctx);

            /* giant steps Q - i (m B), i >= 0, and Q + i (m B), i >= 1 */
            status |= gr_ec_point_neg(Sm, mB, ctx);
            status |= _walk(Q, NULL, 1, 1, mB, Sm, ngiant + 1, _prog_giant_minus, &prog, ctx);

            if (!prog.found && status == GR_SUCCESS)
            {
                status |= gr_ec_point_add(Sp, Q, mB, ctx);
                status |= _walk(Sp, NULL, 1, 1, mB, mB, ngiant, _prog_giant_plus, &prog, ctx);
            }
        }

        if (status != GR_SUCCESS)
            break;

        found = prog.found;

        if (!found)
            continue;

        /* the exact order of P, then the lcm over all points so far */
        status |= _point_order(ord, P, n, ctx);

        if (status != GR_SUCCESS)
            break;

        fmpz_lcm(L, L, ord);

        resolved = _unique_in_interval(res, r, M, L, lo, hi);
    }

    if (status == GR_SUCCESS && !resolved)
        status = GR_UNABLE;

    flint_free(tab);
    GR_TMP_CLEAR2(xc, yc, R);
    gr_ec_point_clear(T, ctx);
    gr_ec_point_clear(Sm, ctx);
    gr_ec_point_clear(Sp, ctx);
    gr_ec_point_clear(S, ctx);
    gr_ec_point_clear(mB, ctx);
    gr_ec_point_clear(B, ctx);
    gr_ec_point_clear(Q, ctx);
    gr_ec_point_clear(P, ctx);
    flint_rand_clear(state);

cleanup_ints:
    fmpz_clear(q); fmpz_clear(W); fmpz_clear(lo); fmpz_clear(hi);
    fmpz_clear(L); fmpz_clear(n); fmpz_clear(ord); fmpz_clear(t);
    fmpz_clear(tc); fmpz_clear(c); fmpz_clear(r); fmpz_clear(cand);
    fmpz_clear(kmin); fmpz_clear(kmax);

    return status;
}

int
gr_ec_ctx_cardinality_bsgs(fmpz_t res, gr_ec_ctx_t ctx)
{
    fmpz_t t0, M;
    int status;

    fmpz_init(t0);
    fmpz_init_set_ui(M, 1);

    status = _gr_ec_cardinality_bsgs_progression(res, t0, M, 0, ctx);

    fmpz_clear(t0);
    fmpz_clear(M);

    return status;
}

/* ------------------------------------------------------------------ */
/* match and sort: BSGS over Atkin's candidate sets                   */
/* ------------------------------------------------------------------ */

/*
    Once t is known modulo m3 (Elkies and Schoof primes) and restricted to
    the sets T_i modulo Atkin primes l_i, the candidates are no longer a
    progression, but they still factor. Split the Atkin primes used into
    a baby side (product m1) and a giant side (product m2), M = m3 m1 m2;
    by the Chinese remainder theorem every candidate is

        t = t3 + m3 m2 u1 + m3 m1 u2 + M z,

    with u1 mod m1 determined by the residues on the baby side, u2 mod m2
    by those on the giant side, and z ranging over what the Hasse bound
    leaves. Splitting z = z1 + Z1 z2 as well, N P = O reads

        (q + 1 - t3 - M zmin) P - u1 (m3 m2 P) - z1 (M P)
            = u2 (m3 m1 P) + z2 (M Z1 P),

    and a hash table of the left sides, looked up with every right side,
    finds all the candidates that kill P (Atkin; Mueller and Lercier made
    it practical). The number of candidates is the product of the set
    sizes with the number of z, and the work is about twice its square
    root when the two sides are balanced.

    This is exact: every match is confirmed by N P = O, all of them are
    collected, and the answer is taken only when a single candidate kills
    every point drawn. The candidate sets themselves come from theorems
    (Atkin's, and Hasse's bound) applied to data checked on the way, see
    _atkin_order in elkies.c.
*/

typedef struct
{
    ulong l;
    double rho;
    slong i;
}
_atkin_sel_struct;

static int
_atkin_sel_cmp(const void * a, const void * b)
{
    double x = ((const _atkin_sel_struct *) a)->rho;
    double y = ((const _atkin_sel_struct *) b)->rho;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/*
    Which Atkin primes to use and on which side (side[i] = 0, 1 baby,
    2 giant), and how many z go to the baby side (*Z1). Returns the
    estimated number of group operations; *use_atkin says whether any
    Atkin set was worth taking. Without any, the search is a plain
    baby-step giant-step on the progression, which unlike
    _gr_ec_cardinality_bsgs_progression never needs to factor a candidate
    to prove it unique.
*/
double
_gr_ec_atkin_plan(int * side, slong * Z1, int * use_atkin,
        const gr_ec_atkin_struct * A, slong nA, const fmpz_t m3, const fmpz_t q)
{
    _atkin_sel_struct * S;
    fmpz_t w;
    double W2, M, P, C, cost, Pb, Pg, nz, root, z1;
    slong i, ns, nsel;

    fmpz_init(w);
    fmpz_sqrt(w, q);
    W2 = 4.0 * fmpz_get_d(w) + 8;      /* the width of the range of t */
    fmpz_clear(w);

    M = fmpz_get_d(m3);

    for (i = 0; i < nA; i++)
        side[i] = 0;

    *Z1 = 1;
    *use_atkin = 0;

    /* take the most informative sets first, as long as they pay */
    S = flint_malloc(FLINT_MAX(nA, 1) * sizeof(_atkin_sel_struct));

    for (i = ns = 0; i < nA; i++)
    {
        /* no information, or none that m3 does not already carry */
        if (A[i].n <= 0 || (ulong) A[i].n >= A[i].l
                || fmpz_fdiv_ui(m3, A[i].l) == 0)
            continue;

        S[ns].l = A[i].l;
        S[ns].rho = (double) A[i].n / A[i].l;
        S[ns].i = i;
        ns++;
    }

    qsort(S, ns, sizeof(_atkin_sel_struct), _atkin_sel_cmp);

    P = 1;
    C = W2 / M + 2;

    for (nsel = 0; nsel < ns; nsel++)
    {
        double M2 = M * S[nsel].l, P2 = P * A[S[nsel].i].n;
        double C2 = (W2 / M2 + 2) * P2;

        if (C2 >= 0.95 * C)
            break;

        M = M2;
        P = P2;
        C = C2;
    }

    /* the larger sets to the baby side, up to about sqrt(C) */
    root = sqrt(C);
    Pb = 1;

    for (i = nsel - 1; i >= 0; i--)
    {
        slong n = A[S[i].i].n;

        if (Pb * n <= root)
        {
            Pb *= n;
            side[S[i].i] = 1;
        }
        else
            side[S[i].i] = 2;
    }

    Pg = P / Pb;
    nz = W2 / M + 2;

    /* the z range, split to balance Pb Z1 against Pg nz / Z1 */
    z1 = sqrt(Pg * nz / Pb);
    z1 = FLINT_MAX(1.0, FLINT_MIN(z1, nz));
    *Z1 = (slong) z1;

    /* each u costs a short scalar multiplication on top */
    cost = Pb * (*Z1) + Pg * ceil(nz / (*Z1)) + 30.0 * (Pb + Pg);

    flint_free(S);

    *use_atkin = (nsel > 0);
    return cost;
}

/*
    All values u mod m (m = prod of the l_i on one side) with u = (r_i - t3)
    / c mod l_i for r_i in T_i, sorted. c is m3 times the other side's
    modulus.
*/
static fmpz *
_atkin_side_values(slong * len, fmpz_t m, const gr_ec_atkin_struct * A,
        slong nA, const int * side, int which, const fmpz_t t3,
        const fmpz_t c)
{
    fmpz * U, * V;
    slong n = 1, i, k, a;

    U = _fmpz_vec_init(1);
    fmpz_one(m);

    for (i = 0; i < nA; i++)
    {
        ulong l, inv, t3l;

        if (side[i] != which)
            continue;

        l = A[i].l;
        inv = n_invmod(fmpz_fdiv_ui(c, l), l);
        t3l = fmpz_fdiv_ui(t3, l);

        V = _fmpz_vec_init(n * A[i].n);

        for (k = 0; k < n; k++)
            for (a = 0; a < A[i].n; a++)
            {
                ulong ui = n_mulmod2(n_submod(A[i].T[a], t3l, l), inv, l);

                fmpz_CRT_ui(V + k * A[i].n + a, U + k, m, ui, l, 0);
            }

        _fmpz_vec_clear(U, n);
        U = V;
        n *= A[i].n;
        fmpz_mul_ui(m, m, l);
    }

    qsort(U, n, sizeof(fmpz), (int (*)(const void *, const void *)) fmpz_cmp);

    *len = n;
    return U;
}

typedef struct
{
    baby_struct * tab;
    slong ntab, Z1;
    const fmpz * tc;
    const fmpz * M;
    const fmpz * zmin;
    const fmpz * c1;            /* m3 m2 */
    const fmpz * c2;            /* m3 m1 */
    const fmpz * U1;
    const fmpz * U2;
    const fmpz * q;
    gr_ec_point_struct * P;
    gr_ec_point_struct * T;
    fmpz ** cand;
    slong * ncand;
    slong * alloc;
    gr_ec_ctx_struct * ctx;
}
_ms_struct;

static int
_ms_baby(void * data, slong k, slong z, ulong h,
        const gr_ec_aff_point_struct * FLINT_UNUSED(pt), int * FLINT_UNUSED(status))
{
    _ms_struct * D = data;

    D->tab[D->ntab].hash = h;
    D->tab[D->ntab].j = k * D->Z1 + z;
    D->ntab++;

    return 0;
}

/* the giant step u2 B2 + z D: every baby with its hash is a candidate */
static int
_ms_giant(void * data, slong k, slong z, ulong h,
        const gr_ec_aff_point_struct * FLINT_UNUSED(pt), int * status)
{
    _ms_struct * D = data;
    gr_ec_ctx_struct * ctx = D->ctx;
    slong idx = _baby_find(D->tab, D->ntab, h);
    fmpz_t t, e;

    if (idx < 0)
        return 0;

    fmpz_init(t);
    fmpz_init(e);

    for ( ; idx < D->ntab && D->tab[idx].hash == h; idx++)
    {
        slong k1 = D->tab[idx].j / D->Z1, z1 = D->tab[idx].j % D->Z1;

        /* t = tc + M (zmin + z1 + Z1 z) + m3 m2 u1 + m3 m1 u2 */
        fmpz_set(t, D->tc);
        fmpz_add_si(e, D->zmin, z1 + D->Z1 * z);
        fmpz_addmul(t, D->M, e);
        fmpz_addmul(t, D->c1, D->U1 + k1);
        fmpz_addmul(t, D->c2, D->U2 + k);

        /* inside the Hasse interval: t^2 <= 4 q */
        fmpz_mul(e, t, t);
        fmpz_submul_ui(e, D->q, 4);

        if (fmpz_sgn(e) > 0)
            continue;

        fmpz_add_ui(e, D->q, 1);
        fmpz_sub(e, e, t);
        *status |= gr_ec_point_mul_fmpz(D->T, D->P, e, ctx);

        if (*status == GR_SUCCESS && gr_ec_point_is_inf(D->T, ctx) == T_TRUE)
        {
            slong i;
            int dup = 0;

            for (i = 0; i < *D->ncand && !dup; i++)
                dup = fmpz_equal(*D->cand + i, e);

            if (!dup)
            {
                if (*D->ncand == *D->alloc)
                {
                    slong old = *D->alloc;
                    *D->alloc = FLINT_MAX(8, 2 * old);
                    *D->cand = flint_realloc(*D->cand, *D->alloc * sizeof(fmpz));
                    for (i = old; i < *D->alloc; i++)
                        fmpz_init(*D->cand + i);
                }

                fmpz_set(*D->cand + *D->ncand, e);
                (*D->ncand)++;
            }
        }
    }

    fmpz_clear(t);
    fmpz_clear(e);

    return 0;
}

int
_gr_ec_cardinality_match_sort(fmpz_t res, const fmpz_t t3, const fmpz_t m3,
        const gr_ec_atkin_struct * A, slong nA, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int * side;
    slong Z1, n1, n2, nz, nz2, i, ncand = 0, alloc = 0, attempt;
    int use_atkin, status = GR_SUCCESS, resolved = 0;
    fmpz * U1 = NULL, * U2 = NULL;
    fmpz * cand = NULL;
    fmpz_t q, W, tc, m1, m2, M, c, zmin, zmax, e, t, N, c1, c2;
    gr_ec_point_t P, Q, B1, B2, Cz, D, Ab, S, T;
    _ms_struct ms;
    gr_ptr xc, yc;
    baby_struct * tab = NULL;
    flint_rand_t state;

    if (gr_ctx_is_field(R) != T_TRUE)
        return GR_DOMAIN;

    fmpz_init(q); fmpz_init(W); fmpz_init(tc); fmpz_init(m1); fmpz_init(m2);
    fmpz_init(M); fmpz_init(c); fmpz_init(zmin); fmpz_init(zmax);
    fmpz_init(e); fmpz_init(t); fmpz_init(N); fmpz_init(c1); fmpz_init(c2);

    side = flint_malloc(FLINT_MAX(nA, 1) * sizeof(int));

    if (gr_ctx_cardinality_fmpz(q, R) != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup_ints;
    }

    _gr_ec_atkin_plan(side, &Z1, &use_atkin, A, nA, m3, q);

    /*
        Without Atkin sets, and while the candidates are small enough to
        factor at no cost, the progression search is cheaper: matching up
        to sign lets each giant step cover twice as many candidates.
    */
    if (!use_atkin && fmpz_bits(q) <= 64)
    {
        status = _gr_ec_cardinality_bsgs_progression(res, t3, m3, 0, ctx);
        goto cleanup_ints;
    }

    fmpz_sqrt(W, q);
    fmpz_mul_ui(W, W, 2);
    fmpz_add_ui(W, W, 2);
    fmpz_smod(tc, t3, m3);

    /* the u values of both sides; c = m3 m2 for the baby side, m3 m1 for
       the giant side, so each needs the other's modulus first */
    for (i = 0, fmpz_one(m1), fmpz_one(m2); i < nA; i++)
    {
        if (side[i] == 1) fmpz_mul_ui(m1, m1, A[i].l);
        if (side[i] == 2) fmpz_mul_ui(m2, m2, A[i].l);
    }

    fmpz_mul(c, m3, m2);
    U1 = _atkin_side_values(&n1, m1, A, nA, side, 1, tc, c);
    fmpz_mul(c, m3, m1);
    U2 = _atkin_side_values(&n2, m2, A, nA, side, 2, tc, c);
    fmpz_mul(M, m3, m1);
    fmpz_mul(M, M, m2);

    /* tc + m3 m2 u1 + m3 m1 u2 lies in [tc, tc + 2M); |t| <= W */
    fmpz_neg(zmin, W);
    fmpz_sub(zmin, zmin, tc);
    fmpz_fdiv_q(zmin, zmin, M);
    fmpz_sub_ui(zmin, zmin, 2);
    fmpz_sub(zmax, W, tc);
    fmpz_cdiv_q(zmax, zmax, M);

    fmpz_sub(e, zmax, zmin);

    if (!fmpz_fits_si(e) || fmpz_get_si(e) > WORD(1) << 40)
    {
        status = GR_UNABLE;
        goto cleanup_ints;
    }

    nz = fmpz_get_si(e) + 1;
    Z1 = FLINT_MIN(Z1, nz);
    nz2 = (nz + Z1 - 1) / Z1;

    flint_rand_init(state);
    flint_rand_set_seed(state, UWORD(0x2545f4914f6cdd1d), UWORD(0x9e3779b97f4a7c15));

    gr_ec_point_init(P, ctx);
    gr_ec_point_init(Q, ctx);
    gr_ec_point_init(B1, ctx);
    gr_ec_point_init(B2, ctx);
    gr_ec_point_init(Cz, ctx);
    gr_ec_point_init(D, ctx);
    gr_ec_point_init(Ab, ctx);
    gr_ec_point_init(S, ctx);
    gr_ec_point_init(T, ctx);
    GR_TMP_INIT2(xc, yc, R);

    tab = flint_malloc(n1 * Z1 * sizeof(baby_struct));

    ms.Z1 = Z1;
    ms.tc = tc;
    ms.M = M;
    ms.zmin = zmin;
    ms.c1 = c1;
    ms.c2 = c2;
    ms.U1 = U1;
    ms.U2 = U2;
    ms.q = q;
    ms.cand = &cand;
    ms.ncand = &ncand;
    ms.alloc = &alloc;
    ms.ctx = ctx;
    fmpz_mul(c1, m3, m2);
    fmpz_mul(c2, m3, m1);

    /* the first point finds every candidate that kills it */
    for (attempt = 0; attempt < GR_EC_BSGS_MAX_POINTS && ncand == 0 && status == GR_SUCCESS; attempt++)
    {

        if (gr_ec_point_randtest(P, state, ctx) != GR_SUCCESS
                || gr_ec_point_is_inf(P, ctx) != T_FALSE)
            continue;

        /* Q = (q + 1 - tc - M zmin) P, B1 = m3 m2 P, Cz = M P,
           B2 = m3 m1 P, D = M Z1 P */
        fmpz_add_ui(e, q, 1);
        fmpz_sub(e, e, tc);
        fmpz_submul(e, M, zmin);
        status |= gr_ec_point_mul_fmpz(Q, P, e, ctx);
        fmpz_mul(e, m3, m2);
        status |= gr_ec_point_mul_fmpz(B1, P, e, ctx);
        status |= gr_ec_point_mul_fmpz(Cz, P, M, ctx);
        fmpz_mul(e, m3, m1);
        status |= gr_ec_point_mul_fmpz(B2, P, e, ctx);
        status |= gr_ec_point_mul_si(D, Cz, Z1, ctx);

        ms.tab = tab;
        ms.ntab = 0;
        ms.P = P;
        ms.T = T;

        /* baby side: Q - u1 B1 - z1 Cz */
        status |= gr_ec_point_neg(S, Cz, ctx);
        status |= _walk(Q, U1, n1, -1, B1, S, Z1, _ms_baby, &ms, ctx);

        if (status != GR_SUCCESS)
            break;

        qsort(tab, ms.ntab, sizeof(baby_struct), _baby_cmp);

        /* giant side: u2 B2 + z2 D, every match confirmed by N P = O */
        status |= gr_ec_point_zero(Ab, ctx);
        status |= _walk(Ab, U2, n2, 1, B2, D, nz2, _ms_giant, &ms, ctx);
    }

    /* more points, on the few candidates left, until one survives */
    for ( ; attempt <= GR_EC_BSGS_MAX_POINTS && status == GR_SUCCESS && ncand > 1; attempt++)
    {
        slong keep = 0;

        if (gr_ec_point_randtest(P, state, ctx) != GR_SUCCESS)
            continue;

        for (i = 0; i < ncand && status == GR_SUCCESS; i++)
        {
            status |= gr_ec_point_mul_fmpz(T, P, cand + i, ctx);

            if (status == GR_SUCCESS && gr_ec_point_is_inf(T, ctx) == T_TRUE)
                fmpz_swap(cand + keep++, cand + i);
        }

        ncand = keep;
    }

    if (status == GR_SUCCESS && ncand == 1)
    {
        fmpz_set(res, cand);
        resolved = 1;
    }

    for (i = 0; i < alloc; i++)
        fmpz_clear(cand + i);
    flint_free(cand);
    flint_free(tab);
    GR_TMP_CLEAR2(xc, yc, R);
    gr_ec_point_clear(T, ctx);
    gr_ec_point_clear(S, ctx);
    gr_ec_point_clear(Ab, ctx);
    gr_ec_point_clear(D, ctx);
    gr_ec_point_clear(Cz, ctx);
    gr_ec_point_clear(B2, ctx);
    gr_ec_point_clear(B1, ctx);
    gr_ec_point_clear(Q, ctx);
    gr_ec_point_clear(P, ctx);
    flint_rand_clear(state);

    if (status == GR_SUCCESS && !resolved)
        status = GR_UNABLE;

cleanup_ints:
    if (U1 != NULL) _fmpz_vec_clear(U1, n1);
    if (U2 != NULL) _fmpz_vec_clear(U2, n2);
    flint_free(side);
    fmpz_clear(q); fmpz_clear(W); fmpz_clear(tc); fmpz_clear(m1); fmpz_clear(m2);
    fmpz_clear(M); fmpz_clear(c); fmpz_clear(zmin); fmpz_clear(zmax);
    fmpz_clear(e); fmpz_clear(t); fmpz_clear(N); fmpz_clear(c1); fmpz_clear(c2);

    return status;
}
