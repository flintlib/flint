/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdlib.h>
#include "fmpz.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_ec.h"
#include "impl.h"

/*
    Schoof's algorithm, for the short model y^2 = f(x) = x^3 + a4 x + a6
    over F_q with q coprime to 6.

    For each small prime l, t mod l comes from the action of Frobenius on
    the l-torsion: phi^2 - t phi + q = 0 there, so phi^2(P) + (q mod l) P
    is matched against t phi(P) for the generic point P. The arithmetic is
    in F_q[x]/(psi_l) and lives in torsion.c; Frobenius acts as
    phi(u, v y) = (u o xq, (v o xq) yq y) with xq = x^q, yq = f^((q-1)/2).

    Three things separate this from the textbook version.

      * The search for t mod l is baby-step giant-step rather than a scan
        over the l candidates. Each candidate costs a torsion addition, a
        torsion addition costs an inversion in F_q[x]/(psi_l), and that
        extended gcd is some thirty times a multiplication there and
        dominates everything else.

      * When psi_l is reducible an inversion can fail. The failed extended
        gcd hands over a proper factor, phi^2 - t phi + q still vanishes
        modulo it, and the computation restarts there, in a smaller ring.

      * The loop over l stops before the product M of the primes covers
        the Hasse interval: the 4 sqrt(q)/M candidates left form a
        progression that baby-step giant-step searches in about their
        square root of group operations. The primes at the top of the
        range cost far more than that, so this is a large saving, and it
        stays a proof (see _gr_ec_cardinality_bsgs_progression).

    SEA is the same loop with Elkies' step (elkies.c) offered first: at
    an Elkies prime, where E has a rational l-isogeny, the torsion ring is
    F_q[x]/(h) for its kernel polynomial h of degree (l - 1)/2 instead of
    (l^2 - 1)/2. An Atkin prime has no such h; Schoof's step still runs
    there for small l, and above that l is skipped, because the next
    Elkies prime and the baby-step giant-step tail are both cheaper.
*/

/*
    t mod l for an odd prime l from phi^2 - t phi + q = 0 on the l-torsion:
    the t with phi^2(P) + (q mod l) P = t phi(P), for the generic point P
    of the subgroup C's modulus describes. _gr_ec_tors_solve runs this
    modulo psi_l and takes care of descending into a factor of it.
*/
static int
_schoof_step(ulong * tl, ulong l, const fmpz_t q, gr_ec_tors_ctx_struct * C,
        gr_ec_ctx_t FLINT_UNUSED(ctx))
{
    gr_ctx_struct * R = C->R;
    gr_ec_tors_struct P, phiP, phi2P, qP, lhs;
    gr_poly_t tmp;
    int status = GR_SUCCESS, found = 0;

    _gr_ec_tors_init(&P, C);
    _gr_ec_tors_init(&phiP, C);
    _gr_ec_tors_init(&phi2P, C);
    _gr_ec_tors_init(&qP, C);
    _gr_ec_tors_init(&lhs, C);
    gr_poly_init(tmp, R);

    status |= _gr_ec_tors_frobenius(&P, &phiP, q, C);

    /* phi^2(P) = (xq o xq, (yq o xq) yq y) */
    status |= gr_poly_preinv_compose_mod(phi2P.u, phiP.u, phiP.u, C->P, R);
    status |= gr_poly_preinv_compose_mod(tmp, phiP.v, phiP.u, C->P, R);
    status |= _gr_ec_tors_mulmod(phi2P.v, tmp, phiP.v, C);
    phi2P.is_inf = 0;

    /* lhs = phi^2(P) + (q mod l) P */
    status |= _gr_ec_tors_mul_ui(&qP, &P, fmpz_fdiv_ui(q, l), C);
    status |= _gr_ec_tors_add(&lhs, &phi2P, &qP, C);

    if (status == GR_SUCCESS)
    {
        if (lhs.is_inf)
        {
            *tl = 0;
            found = 1;
        }
        else
            status |= _gr_ec_tors_discrete_log(tl, &found, &lhs, &phiP, l, C);
    }

    if (status == GR_SUCCESS && !found)
        status = GR_UNABLE;

    gr_poly_clear(tmp, R);
    _gr_ec_tors_clear(&lhs, C);
    _gr_ec_tors_clear(&qP, C);
    _gr_ec_tors_clear(&phi2P, C);
    _gr_ec_tors_clear(&phiP, C);
    _gr_ec_tors_clear(&P, C);

    return status;
}

static int
_schoof_trace_mod_l(ulong * tl, ulong l, const fmpz_t q, gr_ec_ctx_t ctx)
{
    gr_poly_t psi;
    int status;

    gr_poly_init(psi, GR_EC_ELEM_CTX(ctx));

    status = gr_ec_ctx_division_poly(psi, l, ctx);

    if (status == GR_SUCCESS)
        status = _gr_ec_tors_solve(tl, psi, l, q, _schoof_step, ctx);

    gr_poly_clear(psi, GR_EC_ELEM_CTX(ctx));

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

/*
    Once t is known modulo M, what is left of the Hasse interval is a
    progression of about 4 sqrt(q)/M candidates, which baby-step giant-step
    settles in about 2 sqrt(4 sqrt(q)/M) group operations -- see
    _gr_ec_cardinality_bsgs_progression, which also keeps the answer
    proved. The primes at the top of the range dominate the cost of the
    polynomial steps, so the loop stops paying for them as soon as that
    search is cheap enough: at most GR_EC_TAIL_BABY(q) baby steps.
*/
#ifndef GR_EC_TAIL_BABY
#define GR_EC_TAIL_BABY(bits) (WORD(1) << FLINT_MIN(20, (bits) / 16 + 7))
#endif

static slong
_schoof_tail_limit(const fmpz_t q)
{
    return GR_EC_TAIL_BABY(fmpz_bits(q));
}

/* is the search for what is left, with the Atkin sets so far, cheap? */
static int
_schoof_tail_cheap(const fmpz_t M, const gr_ec_atkin_struct * A, slong nA,
        const fmpz_t q)
{
    int * side = flint_malloc(FLINT_MAX(nA, 1) * sizeof(int));
    slong Z1;
    int use_atkin;
    double cost;

    cost = _gr_ec_atkin_plan(side, &Z1, &use_atkin, A, nA, M, q);
    flint_free(side);

    return cost <= 3.0 * _schoof_tail_limit(q);
}

/*
    The order in which SEA takes the odd primes. What an Elkies prime
    costs is dominated by the canonical modular polynomial, whose degree
    in j is v = (l - 1)/gcd(12, l - 1), and measured over the whole step
    it is close to l^2 (v + 4); what it brings is log l bits of t. The
    primes are taken by increasing cost per bit, which puts, say, 73
    (v = 6) before 71 (v = 35). Schoof's step has no such structure and
    takes them in increasing order.

    The list goes far enough for its primes to cover twice the bits of
    the Hasse interval, since about half of them turn out to be Atkin; if
    it still runs out, the caller carries on past its largest prime.
*/
typedef struct
{
    ulong l;
    double key;
}
_sea_prime_struct;

static int
_sea_prime_cmp(const void * a, const void * b)
{
    double x = ((const _sea_prime_struct *) a)->key;
    double y = ((const _sea_prime_struct *) b)->key;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static ulong *
_schoof_prime_order(slong * n, ulong * lmax, const fmpz_t bound, int use_elkies)
{
    _sea_prime_struct * P = NULL;
    ulong * ls;
    double bits = 0, need = 2.0 * fmpz_bits(bound) + 16;
    slong len = 0, alloc = 0, i;
    ulong l;

    for (l = 3; bits < need; l = n_nextprime(l, 1))
    {
        double v = (double) (l - 1) / n_gcd(12, l - 1);

        if (len == alloc)
        {
            alloc = FLINT_MAX(64, 2 * alloc);
            P = flint_realloc(P, alloc * sizeof(_sea_prime_struct));
        }

        P[len].l = l;
        P[len].key = use_elkies ? (double) l * l * (v + 4) / log(l) : (double) l;
        len++;
        bits += log(l) / log(2);
    }

    qsort(P, len, sizeof(_sea_prime_struct), _sea_prime_cmp);

    ls = flint_malloc(len * sizeof(ulong));
    *lmax = 0;

    for (i = 0; i < len; i++)
    {
        ls[i] = P[i].l;
        *lmax = FLINT_MAX(*lmax, P[i].l);
    }

    flint_free(P);
    *n = len;

    return ls;
}

#ifndef GR_EC_SEA_USE_ATKIN
#define GR_EC_SEA_USE_ATKIN 1
#endif

#ifndef GR_EC_SEA_SCHOOF_MAX_L
#define GR_EC_SEA_SCHOOF_MAX_L 11
#endif

static int
_schoof_driver(fmpz_t res, gr_ec_ctx_t ctx, int use_elkies)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_t q, p, bound, M, t, sq;
    ulong l, tl, lmax;
    ulong * ls = NULL;
    slong nls, idx, nA = 0, i;
    gr_ec_atkin_struct * A = NULL;
    int status = GR_SUCCESS, tail_usable = 1;

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

    ls = _schoof_prime_order(&nls, &lmax, bound, use_elkies);

    for (idx = 0, l = lmax; fmpz_cmp(M, bound) <= 0 && status == GR_SUCCESS; idx++)
    {
        if (idx < nls)
            l = ls[idx];
        else
            l = n_nextprime((idx == nls) ? lmax : l, 1);

        /*
            Before paying for another prime, see whether what is left of
            the trace is cheap to settle with the group law instead. The
            primes at the top of the range dominate the whole computation,
            so this is worth asking every time round.
        */
        if (tail_usable && _schoof_tail_cheap(M, A, nA, q))
        {
            if (_gr_ec_cardinality_match_sort(res, t, M, A, nA, ctx)
                    == GR_SUCCESS)
                goto cleanup;

            /* no point resolves it here -- a group of small exponent,
               most likely -- so go back to paying for primes */
            tail_usable = 0;
        }

        if (fmpz_cmp_ui(p, l) == 0)
            continue;

        if (use_elkies)
        {
            int elkies = 0, st;
            ulong r = 0;

            st = _gr_ec_elkies_trace(&tl, &elkies, GR_EC_SEA_USE_ATKIN ? &r : NULL,
                                     l, q, ctx);

            /*
                Schoof's step at an Atkin prime costs degree (l^2 - 1)/2
                arithmetic, far more than the next Elkies prime; above a
                small l it is cheaper to skip l, keeping only Atkin's
                restriction of t mod l for the final search. Only a
                genuine Atkin prime is skipped: when Elkies' step cannot
                run at all, skipping could go on forever.
            */
            if (st == GR_SUCCESS && !elkies && l > GR_EC_SEA_SCHOOF_MAX_L)
            {
                if (r != 0)
                {
                    A = flint_realloc(A, (nA + 1) * sizeof(gr_ec_atkin_struct));
                    A[nA].l = l;
                    A[nA].T = flint_malloc((l + 1) * sizeof(ulong));
                    A[nA].n = _gr_ec_atkin_candidates(A[nA].T, l, r,
                                                      fmpz_fdiv_ui(q, l));
                    nA++;
                }

                continue;
            }

            if (st != GR_SUCCESS || !elkies)
                status = _schoof_trace_mod_l(&tl, l, q, ctx);
        }
        else
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
    for (i = 0; i < nA; i++)
        flint_free(A[i].T);
    flint_free(A);
    flint_free(ls);
    fmpz_clear(q); fmpz_clear(p); fmpz_clear(bound);
    fmpz_clear(M); fmpz_clear(t); fmpz_clear(sq);

    return status;
}

int
gr_ec_ctx_cardinality_schoof(fmpz_t res, gr_ec_ctx_t ctx)
{
    return _schoof_driver(res, ctx, 0);
}

int
gr_ec_ctx_cardinality_sea(fmpz_t res, gr_ec_ctx_t ctx)
{
    return _schoof_driver(res, ctx, 1);
}
