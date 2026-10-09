/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Relations between gamma generators, found by the zero test (called
    from the relation search in certify.c before the LLL round): the
    distribution relations (reflection and Gauss multiplication) are
    theorems, so they are computed exactly rather than proposed
    numerically.
*/

#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpq_mat.h"
#include "fmpz_mpoly_q.h"
#include "arb.h"
#include "acb.h"
#include "arb_fmpz_poly.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* -------------------------------------------------------------------- */
/* gamma values at rational arguments                                   */
/* -------------------------------------------------------------------- */

/*
    Gamma generators at rational arguments p/q in (0, 1) satisfy the
    distribution relations (reflection and Gauss multiplication), which
    involve pi, integers and the numbers sin(pi k/N). These relations are
    theorems, and by the Rohrlich-Lang conjecture they are all the
    algebraic relations between such values: rather than proposing
    relations numerically (LLL), the relations between the gamma
    generators present are computed exactly at the level N = lcm(q) by
    linear algebra over Q (a reduced row echelon form with the columns of
    the other values first, so that the rows with a pivot among the
    generators are relations between generators and constants only).
    A relation with the latest generator g on top,

        g^M = prod_j Gamma_j^(a_j) pi^b prod_n n^(e_n) prod_k sin(pi k/N)^(f_k),

    (integers, Gamma_j earlier generators) makes g algebraic: a root of
    X^M - u over the tower below it, where u is formed from pi and the
    roots of unity (sin(pi k/N) = (z^k - z^(-k)) / (2i), z = exp(pi i/N)),
    which are adjoined and moved to the front when not present. The
    relation is checked numerically as a safeguard against inconsistencies.
*/


/* the rational argument of a gamma generator in (0, 1) */
static int
_gamma_gen_rational(slong * p, slong * q, const gr_tower_gen_struct * g, slong level_limit)
{
    fmpq_t r;
    int ok;

    if (g->kind != GR_TOWER_GAMMA || g->arg.mctx == NULL)
        return 0;

    fmpq_init(r);
    ok = fmpz_mpoly_q_get_fmpq(r, &g->arg.data, g->arg.mctx);
    if (ok)
    {
        fmpq_canonicalise(r);
        ok = fmpz_sgn(fmpq_numref(r)) > 0 && fmpz_cmp(fmpq_numref(r), fmpq_denref(r)) < 0 &&
             fmpz_cmp_ui(fmpq_denref(r), level_limit) <= 0;
    }
    if (ok)
    {
        *p = fmpz_get_si(fmpq_numref(r));
        *q = fmpz_get_si(fmpq_denref(r));
    }
    fmpq_clear(r);
    return ok;
}

/* Whether the modulus of a_k has rational coefficients (the generator
   depends on nothing and may be moved to the front). */
static int
_step_over_q(gr_tower_t T, slong k)
{
    const gr_poly_struct * m = gr_tower_step_minpoly(T, k);
    gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
    fmpq_t c;
    slong i;
    int ok = 1;

    fmpq_init(c);
    for (i = 0; i < m->length && ok; i++)
        ok = (gr_get_fmpq(c, gr_poly_coeff_srcptr(m, i, below), below) == GR_SUCCESS);
    fmpq_clear(c);
    return ok;
}

/*
    A generator giving a primitive (p^e)-th root of unity as a power:
    sets *d and *pw (zeta_{p^e} = gen_d^pw) and returns 1, or returns 0.
*/
static int
_root_of_unity_gen(slong * d, ulong * pw, gr_tower_t T, ulong pe)
{
    slong k;

    for (k = 1; k <= T->length; k++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_STEP(T, k - 1);
        if (g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param > 0 && (ulong) g->def_param % pe == 0)
        {
            *d = g->def_order;
            *pw = (ulong) g->def_param / pe;
            return 1;
        }
    }

    if (pe == 4)
    {
        slong di = _gr_tower_certify_find_i(T);
        if (di >= 0)
        {
            const gr_tower_gen_struct * g = T->gens + di;
            *d = di;
            *pw = (g->def_kind == GR_TOWER_ROOT_OF_UNITY) ? (ulong) g->def_param / 4 : 1;
            return 1;
        }
    }

    return 0;
}

/*
    Makes sure that the generators needed for zeta_m (prime power roots
    of unity) and, if need_pi, pi are present with definition order
    < before (the gamma generator to be eliminated). Returns 1 if the
    tower was changed (adjoining, moving to the front), 0 if everything
    is in place, -1 if this is not possible.
*/
static int
_gamma_prepare_constants(gr_tower_t T, ulong m, int need_pi, slong before)
{
    n_factor_t fac;
    slong i, d;
    ulong pw;

    if (need_pi)
    {
        d = _gr_tower_certify_find_pi(T, T->num_gens);
        if (d < 0)
        {
            if (gr_tower_adjoin_pi(T, NULL) != GR_SUCCESS)
                return -1;
            _gr_tower_move_to_front(T, T->num_gens - 1);
            return 1;
        }
        if (d >= before)
        {
            _gr_tower_move_to_front(T, d);
            return 1;
        }
    }

    if (m > 2)
    {
        n_factor_init(&fac);
        n_factor(&fac, m, 1);
        for (i = 0; i < fac.num; i++)
        {
            ulong pe = n_pow(fac.p[i], fac.exp[i]);

            if (pe == 2)
                continue;

            if (!_root_of_unity_gen(&d, &pw, T, pe))
            {
                if (gr_tower_adjoin_root_of_unity(T, pe, NULL) != GR_SUCCESS)
                    return -1;
                _gr_tower_move_to_front(T, T->num_gens - 1);
                return 1;
            }

            if (d >= before)
            {
                const gr_tower_gen_struct * g = T->gens + d;
                if (g->kind != GR_TOWER_ALGEBRAIC || !_step_over_q(T, g->index))
                    return -1;
                _gr_tower_move_to_front(T, d);
                return 1;
            }
        }
    }

    return 0;
}

/* res = zeta_m^k as a flat element (the generators must be present):
   a product of powers of the roots of prime power orders, reduced once
   (rather than a power of zeta_m, a sum when some root has an exponent
   beyond the degree of its modulus) */
static int
_flat_zeta_pow(fmpz_mpoly_q_t res, ulong m, ulong k, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    n_factor_t fac;
    fmpz_mpoly_q_t t;
    slong i, d;
    ulong pw;

    fmpz_mpoly_q_one(res, F->mctx);
    if (m <= 1)
        return 1;

    fmpz_mpoly_q_init(t, F->mctx);
    n_factor_init(&fac);
    n_factor(&fac, m, 1);

    for (i = 0; i < fac.num; i++)
    {
        ulong pe = n_pow(fac.p[i], fac.exp[i]);
        /* zeta_m = prod zeta_{p^e}^{a_p} with a_p = (m / p^e)^(-1) mod p^e */
        ulong a = n_mulmod2_preinv(n_invmod((m / pe) % pe, pe), k % pe, pe, n_preinvert_limb(pe));

        if (pe == 2)
        {
            if (a % 2 == 1)
                fmpz_mpoly_q_neg(res, res, F->mctx);
            continue;
        }

        if (!_root_of_unity_gen(&d, &pw, T, pe))
        {
            fmpz_mpoly_q_clear(t, F->mctx);
            return 0;
        }

        fmpz_mpoly_q_gen(t, GR_TOWER_FLAT_VAR_D(F, d), F->mctx);
        fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(t), fmpz_mpoly_q_numref(t), pw * a, F->mctx);
        fmpz_mpoly_q_mul(res, res, t, F->mctx);
    }

    fmpz_mpoly_q_clear(t, F->mctx);
    return (gr_tower_flat_reduce(res, F) == GR_SUCCESS);
}

static int
_flat_zeta(fmpz_mpoly_q_t res, ulong m, gr_tower_flat_t F)
{
    return _flat_zeta_pow(res, m, 1, F);
}

/* res *= x^e (e >= 0) by repeated squaring, reduced after each product
   (the powers of a sum of roots of unity, unreduced, have exponents
   and terms growing with e) */
static int
_mul_pow_reduced(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, ulong e, gr_tower_flat_t F)
{
    fmpz_mpoly_q_t b;
    int ok = 1;

    fmpz_mpoly_q_init(b, F->mctx);
    fmpz_mpoly_q_set(b, x, F->mctx);
    while (e != 0 && ok)
    {
        if (e & 1)
        {
            fmpz_mpoly_q_mul(res, res, b, F->mctx);
            ok = (gr_tower_flat_reduce(res, F) == GR_SUCCESS);
        }
        e >>= 1;
        if (e != 0 && ok)
        {
            fmpz_mpoly_q_mul(b, b, b, F->mctx);
            ok = (gr_tower_flat_reduce(b, F) == GR_SUCCESS);
        }
    }
    fmpz_mpoly_q_clear(b, F->mctx);
    return ok;
}

/* res = prod_k sin(pi k/N)^(f_k) (1 <= k <= nsin, k < N/2) as a flat
   element: sin(pi k/N) = (z^k - z^(2N - k)) / (2i), z = zeta_m^(m/(2N))
   (zeta_m and pi must be present: m = 2N, or 4N when the sum of the
   exponents is odd and N is odd) */
static int
_sine_product_flat(fmpz_mpoly_q_t res, slong N, const slong * f, slong nsin, ulong m, gr_tower_flat_t F)
{
    fmpz_mpoly_q_t x, y, z, zi;
    fmpq_t t;
    slong j, fsum = 0, ee;
    ulong s2;
    int ok = 1;

    for (j = 1; j <= nsin; j++)
        fsum += f[j];

    fmpz_mpoly_q_init(x, F->mctx);
    fmpz_mpoly_q_init(y, F->mctx);
    fmpz_mpoly_q_init(z, F->mctx);
    fmpz_mpoly_q_init(zi, F->mctx);
    fmpq_init(t);

    /* 2^(-fsum) */
    fmpq_set_si(t, 2, 1);
    fmpq_pow_si(t, t, -fsum);
    fmpz_mpoly_q_set_fmpq(res, t, F->mctx);

    /* (z^e = zeta_m^(e m / (2N)) directly) */
    ok = _flat_zeta(zi, m, F);
    s2 = m / (2 * N);

    for (j = 1; j <= nsin && ok; j++)
    {
        if (f[j] == 0)
            continue;
        if (f[j] > 0)
        {
            /* z^j - z^(-j) */
            ok = _flat_zeta_pow(x, m, j * s2, F) && _flat_zeta_pow(y, m, (2 * N - j) * s2, F);
            fmpz_mpoly_q_sub(x, x, y, F->mctx);
        }
        else
        {
            /* 1 / (z^j - z^(-j)) = z^j / (w - 1) with w = z^(2j) of
               order r = N / gcd(N, j), and 1 / (w - 1) =
               (1/r) sum_{i<r} i w^i (no algebraic denominator) */
            slong r = N / n_gcd(N, j), ii;
            fmpz_mpoly_q_t w, wp;
            fmpz_mpoly_q_init(w, F->mctx);
            fmpz_mpoly_q_init(wp, F->mctx);
            ok = _flat_zeta_pow(w, m, 2 * j * s2, F);
            fmpz_mpoly_q_zero(x, F->mctx);
            fmpz_mpoly_q_one(wp, F->mctx);
            for (ii = 1; ii < r; ii++)
            {
                fmpz_mpoly_q_mul(wp, wp, w, F->mctx);
                GR_MUST_SUCCEED(gr_tower_flat_reduce(wp, F));
                fmpz_mpoly_q_mul_si(y, wp, ii, F->mctx);
                fmpz_mpoly_q_add(x, x, y, F->mctx);
            }
            fmpz_mpoly_q_div_si(x, x, r, F->mctx);
            ok = ok && _flat_zeta_pow(y, m, j * s2, F);
            fmpz_mpoly_q_mul(x, x, y, F->mctx);
            fmpz_mpoly_q_clear(w, F->mctx);
            fmpz_mpoly_q_clear(wp, F->mctx);
        }
        if (ok)
        {
            GR_MUST_SUCCEED(gr_tower_flat_reduce(x, F));
            ok = _mul_pow_reduced(res, x, FLINT_ABS(f[j]), F);
        }
    }

    /* i^(-fsum) */
    ee = ((-fsum) % 4 + 4) % 4;
    if (ok && ee == 2)
        fmpz_mpoly_q_neg(res, res, F->mctx);
    else if (ok && (ee == 1 || ee == 3))
    {
        ok = _flat_zeta_pow(x, m, m / 4, F);         /* i */
        if (ee == 3)
            fmpz_mpoly_q_neg(x, x, F->mctx);
        fmpz_mpoly_q_mul(res, res, x, F->mctx);
    }

    if (ok && gr_tower_flat_reduce(res, F) != GR_SUCCESS)
        ok = 0;

    fmpz_mpoly_q_clear(x, F->mctx);
    fmpz_mpoly_q_clear(y, F->mctx);
    fmpz_mpoly_q_clear(z, F->mctx);
    fmpz_mpoly_q_clear(zi, F->mctx);
    fmpq_clear(t);
    return ok;
}

/*
    A generator root_{n'}(c) with delta | n' (c a positive flat element,
    the principal root) with definition order < before: sets *d and
    *mult = n' / delta (root_delta(c) = gen_d^mult) and returns 0; or
    adjoins root_delta(c), or moves the existing one, right after the
    generators c depends on, and returns 1; returns -1 if this is not
    possible.
*/
static int
_root_gen(slong * d_out, slong * mult, gr_tower_flat_t F, const fmpz_mpoly_q_t c, ulong delta, slong before)
{
    gr_tower_struct * T = F->T;
    slong k, dep;

    for (k = 1; k <= T->length; k++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_STEP(T, k - 1);
        fmpz_mpoly_q_t diff;
        int eq;

        if (g->def_kind != GR_TOWER_ROOT || g->arg.mctx == NULL || g->def_param <= 0 ||
            (ulong) g->def_param % delta != 0)
            continue;

        fmpz_mpoly_q_init(diff, F->mctx);
        gr_tower_flat_convert(diff, &g->arg.data, g->arg.mctx, F);
        fmpz_mpoly_q_sub(diff, diff, c, F->mctx);
        eq = (gr_tower_flat_reduce(diff, F) == GR_SUCCESS && fmpz_mpoly_q_is_zero(diff, F->mctx));
        fmpz_mpoly_q_clear(diff, F->mctx);
        if (!eq)
            continue;

        if (g->def_order < before)
        {
            *d_out = g->def_order;
            *mult = g->def_param / delta;
            return 0;
        }

        /* move it (only an unsplit modulus X^n - c depends on nothing
           but the generators of c) */
        dep = _gr_tower_flat_max_dep(c, F);
        if (dep + 1 >= before || gr_tower_step_minpoly(T, k)->length != g->def_param + 1)
            return -1;
        _gr_tower_move_gen(T, g->def_order, dep + 1);
        return 1;
    }

    dep = _gr_tower_flat_max_dep(c, F);
    if (dep + 1 >= before)
        return -1;

    if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(c), F->mctx) && fmpz_mpoly_is_one(fmpz_mpoly_q_denref(c), F->mctx))
    {
        fmpz_t p;
        int st;
        fmpz_init(p);
        fmpz_mpoly_get_fmpz(p, fmpz_mpoly_q_numref(c), F->mctx);
        st = gr_tower_adjoin_root_fmpz(T, p, delta, NULL);
        fmpz_clear(p);
        if (st != GR_SUCCESS)
            return -1;
    }
    else
    {
        gr_ctx_struct * top = gr_tower_field(T);
        gr_ptr x;
        int st;
        GR_TMP_INIT(x, top);
        st = gr_tower_flat_get_nested_at(x, c, T->length, F);
        if (st == GR_SUCCESS)
            st = gr_tower_adjoin_root_ui(T, x, delta, NULL);
        GR_TMP_CLEAR(x, top);
        if (st != GR_SUCCESS)
            return -1;
    }

    _gr_tower_move_gen(T, T->num_gens - 1, dep + 1);
    return 1;
}

/* whether the generator is root_n(c) for the flat element c (unsplit
   modulus X^n - c), with n >= 2 */
static int
_is_root_of(const gr_tower_gen_struct * g, gr_tower_flat_t F, const fmpz_mpoly_q_t c)
{
    fmpz_mpoly_q_t diff;
    int eq;

    if (g->kind != GR_TOWER_ALGEBRAIC || g->def_kind != GR_TOWER_ROOT || g->arg.mctx == NULL || g->def_param < 2 ||
        gr_tower_step_minpoly(F->T, g->index)->length != g->def_param + 1)
        return 0;

    fmpz_mpoly_q_init(diff, F->mctx);
    gr_tower_flat_convert(diff, &g->arg.data, g->arg.mctx, F);
    fmpz_mpoly_q_sub(diff, diff, c, F->mctx);
    eq = (gr_tower_flat_reduce(diff, F) == GR_SUCCESS && fmpz_mpoly_q_is_zero(diff, F->mctx));
    fmpz_mpoly_q_clear(diff, F->mctx);
    return eq;
}

/*
    res *= c^(num/den) (c > 0 free of algebraic denominators, den > 0)
    through the root generators root_{n_i}(c) with definition order
    < before: when den divides L = lcm(n_i), c^(num/den) =
    prod_i root_{n_i}(c)^(k_i) with sum k_i/n_i = num/den (an extended
    gcd of the L/n_i), the exponents reduced modulo n_i through powers
    of c (3^(25/74) = sqrt(3)^a root_37(3)^b). Otherwise root_{l^e}(c)
    is adjoined for a prime power l^e | den not dividing L (one at a
    time, never a composite order) or an existing one moved: returns 1
    (the tower changed), -1 if impossible, 0 when done.
*/
static int
_mul_root_power(fmpz_mpoly_q_t res, gr_tower_flat_t F, const fmpz_mpoly_q_t c, slong num, slong den, slong before)
{
    gr_tower_struct * T = F->T;
    slong g, i, nr = 0, L = 1;
    slong * dd, * nn;
    int st = 0;

    if (num == 0)
        return 0;

    g = n_gcd(FLINT_ABS(num), den);
    num /= g;
    den /= g;

    if (den == 1)
    {
        _gr_tower_certify_mul_pow_si(res, c, num, F);
        return 0;
    }

    dd = flint_malloc(sizeof(slong) * FLINT_MAX(T->num_gens, 1));
    nn = flint_malloc(sizeof(slong) * FLINT_MAX(T->num_gens, 1));

    for (i = 0; i < before && i < T->num_gens; i++)
    {
        if (_is_root_of(T->gens + i, F, c))
        {
            dd[nr] = i;
            nn[nr] = T->gens[i].def_param;
            L = (L / n_gcd(L, nn[nr])) * nn[nr];
            nr++;
        }
    }

    if (L % den != 0)
    {
        /* adjoin (or move) the missing roots one prime power at a time:
           root_74(3) next to sqrt(3) and root_37(3) in the lazy layer's
           towers would make a tower which is not a field */
        n_factor_t fac;
        slong d, mult;
        n_factor_init(&fac);
        n_factor(&fac, den, 1);
        st = -1;
        for (i = 0; i < fac.num; i++)
        {
            slong pe = n_pow(fac.p[i], fac.exp[i]);
            if (L % pe != 0)
            {
                /* sqrt of an integer as a Gauss sum in the roots of unity
                   below: c^(num/den) = sqrt(c) c^(b/m), den = 2 m */
                if (pe == 2 && fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(c), F->mctx) &&
                    fmpz_mpoly_is_one(fmpz_mpoly_q_denref(c), F->mctx))
                {
                    fmpz_t cz;
                    fmpz_mpoly_q_t s;
                    int found;
                    fmpz_init(cz);
                    fmpz_mpoly_q_init(s, F->mctx);
                    fmpz_mpoly_get_fmpz(cz, fmpz_mpoly_q_numref(c), F->mctx);
                    found = _gr_tower_cyclotomic_sqrt(s, cz, F, before);
                    if (found)
                    {
                        slong m = den / 2;
                        fmpz_mpoly_q_mul(res, res, s, F->mctx);
                        st = _mul_root_power(res, F, c, (num - m) / 2, m, before);
                    }
                    fmpz_mpoly_q_clear(s, F->mctx);
                    fmpz_clear(cz);
                    if (found)
                        break;
                }

                st = _root_gen(&d, &mult, F, c, pe, before);
                if (st == 0)
                    st = -1;    /* (found below `before`: it was collected above) */
                break;
            }
        }
        flint_free(dd);
        flint_free(nn);
        return st;
    }

    {
        /* 1 = sum u_i (L / n_i) */
        fmpz * u = _fmpz_vec_init(nr);
        fmpz_t gg, s, t, a, b, k, q, r;
        fmpz_mpoly_q_t x;

        fmpz_init(gg); fmpz_init(s); fmpz_init(t); fmpz_init(a); fmpz_init(b);
        fmpz_init(k); fmpz_init(q); fmpz_init(r);
        fmpz_mpoly_q_init(x, F->mctx);

        fmpz_set_si(gg, L / nn[0]);
        fmpz_one(u + 0);
        for (i = 1; i < nr; i++)
        {
            slong j;
            fmpz_set_si(b, L / nn[i]);
            fmpz_xgcd(a, s, t, gg, b);     /* a = s gg + t b */
            for (j = 0; j < i; j++)
                fmpz_mul(u + j, u + j, s);
            fmpz_set(u + i, t);
            fmpz_swap(gg, a);
        }

        /* (gg = 1 since gcd(L / n_i) = 1) */
        fmpz_zero(s);
        for (i = 0; i < nr; i++)
        {
            /* k_i = num (L / den) u_i, reduced modulo n_i: x_i^(k_i) =
               x_i^(r_i) c^(q_i) with x_i^(n_i) = c */
            fmpz_mul_si(k, u + i, num * (L / den));
            fmpz_set_si(t, nn[i]);
            fmpz_fdiv_qr(q, r, k, t);
            fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, dd[i]), F->mctx);
            _gr_tower_certify_mul_pow_si(res, x, fmpz_get_si(r), F);
            fmpz_add(s, s, q);
        }
        /* the powers of c at once: the q_i (from the Bezout coefficients
           u_i) may be large, their sum is num / den - sum r_i / n_i,
           small */
        if (!fmpz_is_zero(s))
            _gr_tower_certify_mul_pow_si(res, c, fmpz_get_si(s), F);

        _fmpz_vec_clear(u, nr);
        fmpz_clear(gg); fmpz_clear(s); fmpz_clear(t); fmpz_clear(a); fmpz_clear(b);
        fmpz_clear(k); fmpz_clear(q); fmpz_clear(r);
        fmpz_mpoly_q_clear(x, F->mctx);
        st = 0;
    }
    flint_free(dd);
    flint_free(nn);
    return st;
}

/*
    Elimination of the gamma generator g (definition order dtop) by
    g = prod_j Gamma_j^(a_j/M) pi^(b/M) prod_n n^(e_n/M) prod sin^(f_k/M)
    as a degree-one modulus over the tower below g, with the radicals
    as root generators (existing ones are reused: sqrt(pi), root_4(2),
    ...). Returns 1 if the tower was changed (the elimination, or a root
    adjoined or moved: the round restarts), 0 if not applicable.
*/
static int
_gamma_eliminate_linear(gr_tower_flat_t F, slong dtop, slong ns, const slong * sd, const slong * a, slong top,
    slong M, slong b, slong N, const slong * e, slong nsin, const slong * f, ulong m)
{
    gr_tower_struct * T = F->T;
    fmpz_mpoly_q_t expr, c, x;
    slong j, ds = 1, st = 0;
    slong * F2;
    int ok = 1;

    for (j = 0; j < ns; j++)
        if (j != top && a[j] % M != 0)
            return 0;

    fmpz_mpoly_q_init(expr, F->mctx);
    fmpz_mpoly_q_init(c, F->mctx);
    fmpz_mpoly_q_init(x, F->mctx);
    F2 = flint_calloc(nsin + 1, sizeof(slong));

    fmpz_mpoly_q_one(expr, F->mctx);

    for (j = 0; j < ns; j++)
    {
        if (j != top && a[j] != 0)
        {
            fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, sd[j]), F->mctx);
            _gr_tower_certify_mul_pow_si(expr, x, a[j] / M, F);
        }
    }

    /* pi */
    if (b != 0)
    {
        fmpz_mpoly_q_gen(c, GR_TOWER_FLAT_VAR_D(F, _gr_tower_certify_find_pi(T, T->num_gens)), F->mctx);
        st = _mul_root_power(expr, F, c, b, M, dtop);
    }

    /* integers: by primes */
    if (st == 0)
    {
        slong p;
        for (p = 2; p <= N && st == 0; p++)
        {
            slong E = 0, n;
            if (!n_is_prime(p))
                continue;
            for (n = 2; n <= N; n++)
            {
                if (e[n] != 0)
                {
                    slong v = 0, t = n;
                    while (t % p == 0)
                    {
                        t /= p;
                        v++;
                    }
                    E += v * e[n];
                }
            }
            if (E != 0)
            {
                fmpz_mpoly_q_set_si(c, p, F->mctx);
                st = _mul_root_power(expr, F, c, E, M, dtop);
            }
        }
    }

    /* sines: one radical root_ds(prod sin^(f_k ds / M)) */
    if (st == 0 && m > 1)
    {
        for (j = 1; j <= nsin; j++)
            if (f[j] != 0)
                ds = (ds / n_gcd(ds, M / n_gcd(FLINT_ABS(f[j]), M))) * (M / n_gcd(FLINT_ABS(f[j]), M));
        for (j = 1; j <= nsin; j++)
            F2[j] = (f[j] / (M / ds));     /* f_k ds / M, an integer */
        {
            ulong m2 = 2 * N;
            slong fs = 0;
            for (j = 1; j <= nsin; j++)
                fs += F2[j];
            if (fs % 2 != 0 && m2 % 4 != 0)
                m2 *= 2;
            if (m2 > m)
                st = -1;   /* (zeta_m2 was not prepared) */
            else
                m = m2;
        }
        if (st == 0)
            ok = _sine_product_flat(c, N, F2, nsin, m, F);
        if (st == 0 && ok)
            st = _mul_root_power(expr, F, c, 1, ds, dtop);
    }

    if (st == 0 && ok)
    {
        /* modulus X - expr */
        if (gr_tower_flat_reduce(expr, F) == GR_SUCCESS && _gr_tower_flat_eliminate_gen(dtop, expr, F))
            st = 1;
    }

    fmpz_mpoly_q_clear(expr, F->mctx);
    fmpz_mpoly_q_clear(c, F->mctx);
    fmpz_mpoly_q_clear(x, F->mctx);
    flint_free(F2);

    if (st == 1)
        return 1;
    return 0;
}

/*
    The relation of row `row` of the reduced row echelon form B of the
    relation matrix A (nrows x ncols, unknowns in the first nvar columns)
    rewritten as an integral combination of the rows of A (up to a
    common factor): the part on the unknowns is the same up to scaling,
    but the constant part is then an integral combination of the original
    relations, in which the sines have integer exponents (the reduced row
    expresses the sines in a basis of the sine values, whose coordinates
    may have denominators: fourth roots of sines where the relation is
    a duplication formula). The row is replaced, scaled so that its pivot
    is 1. Returns 1 on success, and then sets P (when not NULL) to the
    absolute value of the pivot coefficient of the integral combination
    before the scaling: the relations holding modulo 2 pi i, the scaled
    row holds modulo 2 pi i / P.
*/
/*
    The Hermite normal form of (A_u | I) for the unknown part A_u of the
    relation matrix at level N with the unknowns in their natural order
    (x_k in column k - 1): it depends only on N, and the generators of a
    tower select a column order (col_of_k), which only permutes the
    unknown columns (the rows of _gr_tower_gamma_relations do not depend
    on it). The forms are cached per thread: in a long session the
    towers come back to the same levels (a merge of two computations with
    gamma values at one level repeats the eliminations at that level),
    and the form is the expensive part of a relation (seconds at level
    210).
*/
#define GAMMA_HNF_CACHE_SIZE 4

static FLINT_TLS_PREFIX fmpz_mat_struct _gamma_hnf_cache[GAMMA_HNF_CACHE_SIZE];
static FLINT_TLS_PREFIX slong _gamma_hnf_cache_N[GAMMA_HNF_CACHE_SIZE];
static FLINT_TLS_PREFIX slong _gamma_hnf_cache_next = 0;
static FLINT_TLS_PREFIX int _gamma_hnf_cache_registered = 0;

static void
_gamma_hnf_cache_cleanup(void)
{
    slong i;
    for (i = 0; i < GAMMA_HNF_CACHE_SIZE; i++)
    {
        if (_gamma_hnf_cache_N[i] != 0)
            fmpz_mat_clear(_gamma_hnf_cache + i);
        _gamma_hnf_cache_N[i] = 0;
    }
    _gamma_hnf_cache_next = 0;
    _gamma_hnf_cache_registered = 0;
}

/* (valid until the next call) */
static const fmpz_mat_struct *
_gamma_level_hnf(slong N)
{
    slong i, j, nvar = N - 1, nrows;
    slong * col_of_k;
    fmpq_mat_t A;
    fmpz_mat_t Z;
    fmpz_mat_struct * H;

    for (i = 0; i < GAMMA_HNF_CACHE_SIZE; i++)
        if (_gamma_hnf_cache_N[i] == N)
            return _gamma_hnf_cache + i;

    if (!_gamma_hnf_cache_registered)
    {
        flint_register_cleanup_function(_gamma_hnf_cache_cleanup);
        _gamma_hnf_cache_registered = 1;
    }

    col_of_k = flint_malloc(sizeof(slong) * N);
    for (i = 1; i < N; i++)
        col_of_k[i] = i - 1;
    fmpq_mat_init(A, _gr_tower_gamma_relations_rows(N), GR_TOWER_GAMMA_NCOLS(N));
    nrows = _gr_tower_gamma_relations(A, col_of_k, N);

    fmpz_mat_init(Z, nrows, nvar + nrows);
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < nvar; j++)
            fmpz_set(fmpz_mat_entry(Z, i, j), fmpq_numref(fmpq_mat_entry(A, i, j)));
        fmpz_one(fmpz_mat_entry(Z, i, nvar + i));
    }

    i = _gamma_hnf_cache_next;
    _gamma_hnf_cache_next = (i + 1) % GAMMA_HNF_CACHE_SIZE;
    H = _gamma_hnf_cache + i;
    if (_gamma_hnf_cache_N[i] != 0)
        fmpz_mat_clear(H);
    fmpz_mat_init(H, nrows, nvar + nrows);
    /* (the entries stay small: the classical algorithm is several times
       faster on these sparse matrices than the default, which is
       tuned for dense ones; 0.3 s instead of 1.5 s at level 210) */
    fmpz_mat_hnf_classical(H, Z);
    _gamma_hnf_cache_N[i] = N;

    fmpz_mat_clear(Z);
    fmpq_mat_clear(A);
    flint_free(col_of_k);
    return H;
}

/*
    With Hc not NULL, the Hermite normal form is that of
    _gamma_level_hnf, and the unknown in column perm[j] of A is the one in
    column j of Hc (the entries of A_u being integers).
*/
static int
_gamma_integral_relation(fmpq_mat_t B, slong row, const fmpq_mat_t A, slong nrows, slong nvar, slong ncols, slong pivot, fmpz_t P, const fmpz_mat_struct * Hc, const slong * perm)
{
    fmpz_mat_t Z, H0;
    const fmpz_mat_struct * H;
    fmpq * res, * c;
    fmpz * r;
    fmpz_t den;
    slong i, j, k;
    int ok = 1;

    res = _fmpq_vec_init(nvar);
    c = _fmpq_vec_init(nrows);
    r = _fmpz_vec_init(nvar);
    fmpz_init(den);

    if (Hc != NULL)
    {
        H = Hc;
    }
    else
    {
        fmpz_mat_init(Z, nrows, nvar + nrows);
        fmpz_mat_init(H0, nrows, nvar + nrows);
        H = H0;

        for (i = 0; i < nrows && ok; i++)
        {
            for (j = 0; j < nvar && ok; j++)
            {
                const fmpq * x = fmpq_mat_entry(A, i, j);
                if (!fmpz_is_one(fmpq_denref(x)))
                    ok = 0;
                else
                    fmpz_set(fmpz_mat_entry(Z, i, j), fmpq_numref(x));
            }
            fmpz_one(fmpz_mat_entry(Z, i, nvar + i));
        }
    }

    /* the target: the unknown part of the row, as a primitive integer
       vector */
    fmpz_one(den);
    for (j = 0; j < nvar; j++)
        fmpz_lcm(den, den, fmpq_denref(fmpq_mat_entry(B, row, j)));
    for (j = 0; j < nvar; j++)
    {
        fmpz_mul(r + j, fmpq_numref(fmpq_mat_entry(B, row, j)), den);
        fmpz_divexact(r + j, r + j, fmpq_denref(fmpq_mat_entry(B, row, j)));
    }
    for (j = 0; j < nvar; j++)
        fmpq_set_fmpz(res + j, r + (perm != NULL ? perm[j] : j));

    if (ok)
    {
        if (Hc == NULL)
            fmpz_mat_hnf(H0, Z);

        /* r = sum_i c_i H_i (unknown parts), H in echelon form */
        for (i = 0; i < nrows && ok; i++)
        {
            for (k = 0; k < nvar; k++)
                if (!fmpz_is_zero(fmpz_mat_entry(H, i, k)))
                    break;
            if (k == nvar)
                break;
            fmpq_div_fmpz(c + i, res + k, fmpz_mat_entry(H, i, k));
            if (!fmpq_is_zero(c + i))
            {
                for (j = k; j < nvar; j++)
                {
                    fmpq_t t;
                    fmpq_init(t);
                    fmpq_mul_fmpz(t, c + i, fmpz_mat_entry(H, i, j));
                    fmpq_sub(res + j, res + j, t);
                    fmpq_clear(t);
                }
            }
        }
        for (j = 0; j < nvar && ok; j++)
            if (!fmpq_is_zero(res + j))
                ok = 0;
    }

    if (ok)
    {
        /* lambda = D sum_i c_i U_i, the row = lambda A / (its pivot) */
        fmpq * lam = _fmpq_vec_init(nrows);
        fmpq * v = _fmpq_vec_init(ncols);
        fmpq_t t;
        fmpq_init(t);

        for (i = 0; i < nrows; i++)
            for (k = 0; k < nrows; k++)
                if (!fmpq_is_zero(c + i) && !fmpz_is_zero(fmpz_mat_entry(H, i, nvar + k)))
                {
                    fmpq_mul_fmpz(t, c + i, fmpz_mat_entry(H, i, nvar + k));
                    fmpq_add(lam + k, lam + k, t);
                }

        for (k = 0; k < nrows; k++)
            if (!fmpq_is_zero(lam + k))
                for (j = 0; j < ncols; j++)
                    if (!fmpq_is_zero(fmpq_mat_entry(A, k, j)))
                    {
                        fmpq_mul(t, lam + k, fmpq_mat_entry(A, k, j));
                        fmpq_add(v + j, v + j, t);
                    }

        if (fmpq_is_zero(v + pivot))
            ok = 0;
        else
        {
            if (P != NULL)
            {
                /* the combination made integral: L lambda with L the
                   lcm of the denominators; its pivot coefficient L v_pivot */
                fmpz_t L;
                fmpz_init(L);
                fmpz_one(L);
                for (k = 0; k < nrows; k++)
                    fmpz_lcm(L, L, fmpq_denref(lam + k));
                fmpq_mul_fmpz(t, v + pivot, L);
                fmpz_abs(P, fmpq_numref(t));
                fmpz_clear(L);
            }
            fmpq_inv(t, v + pivot);
            for (j = 0; j < ncols; j++)
                fmpq_mul(fmpq_mat_entry(B, row, j), v + j, t);
        }

        fmpq_clear(t);
        _fmpq_vec_clear(lam, nrows);
        _fmpq_vec_clear(v, ncols);
    }

    if (Hc == NULL)
    {
        fmpz_mat_clear(Z);
        fmpz_mat_clear(H0);
    }
    _fmpq_vec_clear(res, nvar);
    _fmpq_vec_clear(c, nrows);
    _fmpz_vec_clear(r, nvar);
    fmpz_clear(den);
    return ok;
}

/* the relations at the level N between the selected gamma generators
   (definition orders sel_d, arguments sel_p / sel_q, latest first) */
static int
_gamma_round_level(gr_tower_flat_t F, slong N, const slong * sel_d, const slong * sel_p, const slong * sel_q, slong nsel)
{
    gr_tower_struct * T = F->T;
    slong * sd, * sk, * sq, ns = nsel, i, j;
    slong nvar, ncols, nrows, rank, r, top = -1, row = -1, ntopcol = 0;
    slong * col_of_k, * s_of_col;
    fmpq_mat_t A, B;
    int changed = 0;

    sd = flint_malloc(sizeof(slong) * FLINT_MAX(ns, 1));
    sk = flint_malloc(sizeof(slong) * FLINT_MAX(ns, 1));
    sq = flint_malloc(sizeof(slong) * FLINT_MAX(ns, 1));
    for (i = 0; i < ns; i++)
    {
        sd[i] = sel_d[i];
        sk[i] = sel_p[i];
        sq[i] = sel_q[i];
    }

    for (i = 0; i < ns; i++)
        sk[i] = sk[i] * (N / sq[i]);
    flint_free(sq);

    nvar = N - 1;
    ncols = GR_TOWER_GAMMA_NCOLS(N);
    col_of_k = flint_malloc(sizeof(slong) * N);
    s_of_col = flint_malloc(sizeof(slong) * nvar);

    /* column order: the other values, then the generators (latest
       first); a value present twice is a relation by itself */
    for (i = 0; i < N; i++)
        col_of_k[i] = -1;
    for (i = 0; i < nvar; i++)
        s_of_col[i] = -1;
    {
        slong c = 0;
        for (i = 1; i < N; i++)
        {
            int present = 0;
            for (j = 0; j < ns; j++)
                if (sk[j] == i)
                    present = 1;
            if (!present)
                col_of_k[i] = c++;
        }
        ntopcol = c;
        for (j = 0; j < ns; j++)
        {
            if (col_of_k[sk[j]] < 0)
            {
                col_of_k[sk[j]] = c;
                s_of_col[c] = j;
                c++;
            }
            else
            {
                /* (a duplicate: the same value as a later generator,
                   which is eliminated first) */
            }
        }
    }

    /* duplicates: the earlier of two generators with the same argument
       remains; the later one equals it */
    for (i = 0; i < ns && top < 0; i++)
        for (j = i + 1; j < ns && top < 0; j++)
            if (sk[i] == sk[j])
                top = i;    /* sd[i] > sd[j] */

    if (top >= 0)
    {
        /* g = Gamma_j: modulus X - g_j */
        slong jj;
        fmpz_mpoly_q_t g;
        for (jj = top + 1; jj < ns; jj++)
            if (sk[jj] == sk[top])
                break;
        fmpz_mpoly_q_init(g, F->mctx);
        fmpz_mpoly_q_gen(g, GR_TOWER_FLAT_VAR_D(F, sd[jj]), F->mctx);
        _gr_tower_flat_eliminate_gen(sd[top], g, F);
        fmpz_mpoly_q_clear(g, F->mctx);
        flint_free(col_of_k);
        flint_free(s_of_col);
        flint_free(sd);
        flint_free(sk);
        return 1;
    }

    fmpq_mat_init(A, _gr_tower_gamma_relations_rows(N), ncols);
    nrows = _gr_tower_gamma_relations(A, col_of_k, N);
    {
        fmpq_mat_t W;
        fmpq_mat_window_init(W, A, 0, 0, nrows, ncols);
        fmpq_mat_init(B, nrows, ncols);
        rank = fmpq_mat_rref(B, W);
        fmpq_mat_window_clear(W);
    }

    /* the first row with its pivot among the generators */
    for (r = 0; r < rank && row < 0; r++)
    {
        for (j = 0; j < ncols; j++)
            if (!fmpq_is_zero(fmpq_mat_entry(B, r, j)))
                break;
        if (j >= ntopcol && j < nvar)
        {
            row = r;
            top = s_of_col[j];
        }
        else if (j >= nvar)
            break;
    }

    /* (the constant part as an integral combination of the relations;
       if this fails, the reduced row is used as it is) */
    if (row >= 0)
    {
        slong pivot;
        for (pivot = 0; pivot < ncols; pivot++)
            if (!fmpq_is_zero(fmpq_mat_entry(B, row, pivot)))
                break;
        slong * perm = flint_malloc(sizeof(slong) * FLINT_MAX(nvar, 1));
        for (j = 0; j < nvar; j++)
            perm[j] = col_of_k[j + 1];
        (void) _gamma_integral_relation(B, row, A, nrows, nvar, ncols, pivot, NULL, _gamma_level_hnf(N), perm);
        flint_free(perm);
    }

    if (row >= 0)
    {
        /* g^M = prod_j Gamma_j^(a_j) pi^b prod n^(e_n) prod sin^(f_k) */
        slong cpi = GR_TOWER_GAMMA_CPI(N), cint = GR_TOWER_GAMMA_CINT(N), csin = GR_TOWER_GAMMA_CSIN(N), nsin = (N - 1) / 2;
        slong M, b, fsum = 0, dtop = sd[top];
        slong * a = flint_calloc(ns, sizeof(slong));
        slong * e = flint_calloc(N + 1, sizeof(slong));
        slong * f = flint_calloc(nsin + 1, sizeof(slong));
        fmpz_t den;
        fmpq_t t;
        int ok = 1, need_pi, prep;
        ulong m = 1;

        fmpz_init(den);
        fmpq_init(t);
        fmpz_one(den);
        for (j = 0; j < ncols; j++)
            fmpz_lcm(den, den, fmpq_denref(fmpq_mat_entry(B, row, j)));
        ok = fmpz_fits_si(den) && fmpz_cmp_ui(den, 1000) <= 0;
        M = ok ? fmpz_get_si(den) : 0;

        for (j = ntopcol; j < nvar && ok; j++)
        {
            if (s_of_col[j] != top && !fmpq_is_zero(fmpq_mat_entry(B, row, j)))
            {
                fmpq_mul_si(t, fmpq_mat_entry(B, row, j), -M);
                a[s_of_col[j]] = fmpz_get_si(fmpq_numref(t));
            }
        }
        fmpq_mul_si(t, fmpq_mat_entry(B, row, cpi), -M);
        b = fmpz_get_si(fmpq_numref(t));
        for (j = 2; j <= N; j++)
        {
            fmpq_mul_si(t, fmpq_mat_entry(B, row, cint + j - 2), -M);
            e[j] = fmpz_get_si(fmpq_numref(t));
        }
        for (j = 1; j <= nsin; j++)
        {
            fmpq_mul_si(t, fmpq_mat_entry(B, row, csin + j - 1), -M);
            f[j] = fmpz_get_si(fmpq_numref(t));
            fsum += f[j];
        }

        /* safeguard: the relation holds numerically (logarithms of
           positive reals) */
        if (ok)
        {
            arb_t s, v, w;
            slong prec = 128;
            arb_init(s);
            arb_init(v);
            arb_init(w);
            for (j = 0; j < ns; j++)
            {
                slong c = (j == top) ? M : -a[j];
                if (c == 0)
                    continue;
                arb_set_si(v, sk[j]);
                arb_div_si(v, v, N, prec);
                arb_lgamma(v, v, prec);
                arb_mul_si(v, v, c, prec);
                arb_add(s, s, v, prec);
            }
            arb_const_pi(v, prec);
            arb_log(v, v, prec);
            arb_mul_si(v, v, b, prec);
            arb_sub(s, s, v, prec);
            for (j = 2; j <= N; j++)
            {
                if (e[j] == 0)
                    continue;
                arb_set_si(v, j);
                arb_log(v, v, prec);
                arb_mul_si(v, v, e[j], prec);
                arb_sub(s, s, v, prec);
            }
            for (j = 1; j <= nsin; j++)
            {
                if (f[j] == 0)
                    continue;
                arb_set_si(w, j);
                arb_div_si(w, w, N, prec);
                arb_sin_pi(v, w, prec);
                arb_log(v, v, prec);
                arb_mul_si(v, v, f[j], prec);
                arb_sub(s, s, v, prec);
            }
            ok = arb_contains_zero(s);
            arb_clear(s);
            arb_clear(v);
            arb_clear(w);
        }

        /* the constants: pi, and zeta_m with z = exp(pi i/N) a power of
           it (and i, when the power of 2i is odd) */
        need_pi = (b != 0);
        if (nsin > 0)
            for (j = 1; j <= nsin; j++)
                if (f[j] != 0)
                    m = 2 * N;
        /* (with i, for odd powers of 2i in any of the products formed) */
        if (m > 1 && m % 4 != 0)
            m = 2 * m;

        if (ok)
        {
            prep = _gamma_prepare_constants(T, m, need_pi, dtop);
            if (prep != 0)
            {
                ok = 0;
                changed = (prep > 0);
            }
        }

        if (ok && _gamma_eliminate_linear(F, dtop, ns, sd, a, top, M, b, N, e, nsin, f, m))
        {
            ok = 0;
            changed = 1;
        }

        if (ok)
        {
            fmpz_mpoly_q_t u, x;
            fmpz_mpoly_q_struct * mod;
            fmpq_t R;
            slong len, cc = M;

            gr_tower_flat_ensure(F);
            fmpz_mpoly_q_init(u, F->mctx);
            fmpz_mpoly_q_init(x, F->mctx);
            fmpq_init(R);

            /* rational part */
            fmpq_one(R);
            for (j = 2; j <= N; j++)
            {
                if (e[j] != 0)
                {
                    fmpq_set_si(t, j, 1);
                    fmpq_pow_si(t, t, e[j]);
                    fmpq_mul(R, R, t);
                }
            }
            fmpz_mpoly_q_set_fmpq(u, R, F->mctx);

            for (j = 0; j < ns; j++)
            {
                if (j != top && a[j] != 0)
                {
                    fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, sd[j]), F->mctx);
                    _gr_tower_certify_mul_pow_si(u, x, a[j], F);
                }
            }

            if (b != 0)
            {
                fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, _gr_tower_certify_find_pi(T, T->num_gens)), F->mctx);
                _gr_tower_certify_mul_pow_si(u, x, b, F);
            }

            if (m > 1)
            {
                ok = _sine_product_flat(x, N, f, nsin, m, F);
                if (ok)
                    fmpz_mpoly_q_mul(u, u, x, F->mctx);
            }

            if (ok && gr_tower_flat_reduce(u, F) != GR_SUCCESS)
                ok = 0;
            if (ok && gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(u), F))
                ok = (gr_tower_flat_rationalize(u, F) == GR_SUCCESS);

            if (ok)
            {
                /* modulus X^M - u (an irreducible factor when u is
                   rational, or a power times a rational) */
                mod = flint_malloc(sizeof(fmpz_mpoly_q_struct) * (cc + 1));
                for (i = 0; i <= cc; i++)
                    fmpz_mpoly_q_init(mod + i, F->mctx);
                fmpz_mpoly_q_neg(mod + 0, u, F->mctx);
                fmpz_mpoly_q_one(mod + cc, F->mctx);
                len = cc + 1;
                if (cc > 1 && fmpz_mpoly_q_is_fmpq(u, F->mctx))
                    len = _gr_tower_certify_rational_radical_factor(mod, cc, u, T, T->gens[dtop].index, F);
                else if (cc > 1)
                    len = _gr_tower_certify_radical_factor(mod, cc, u, T, T->gens[dtop].index, F);

                if (_gr_tower_make_algebraic(T, dtop, mod, len, F->mctx,
                        (len == 2) ? GR_TOWER_STATUS_PROVEN : GR_TOWER_STATUS_DYNAMIC))
                    changed = 1;

                for (i = 0; i <= cc; i++)
                    fmpz_mpoly_q_clear(mod + i, F->mctx);
                flint_free(mod);
            }

            fmpz_mpoly_q_clear(u, F->mctx);
            fmpz_mpoly_q_clear(x, F->mctx);
            fmpq_clear(R);
        }

        flint_free(a);
        flint_free(e);
        flint_free(f);
        fmpz_clear(den);
        fmpq_clear(t);
    }

    fmpq_mat_clear(A);
    fmpq_mat_clear(B);
    flint_free(col_of_k);
    flint_free(s_of_col);
    flint_free(sd);
    flint_free(sk);
    return changed;
}

/* -------------------------------------------------------------------- */
/* gamma values on rational affine lines                                 */
/* -------------------------------------------------------------------- */

/*
    Gamma generators whose arguments are z_j = a_j w + b_j for a common
    w (irrational: sqrt(2), pi, ...) with a_j, b_j rational satisfy the
    shift, reflection and multiplication relations along the line:

        Gamma(z + 1) = z Gamma(z),
        Gamma(z) Gamma(1 - z) = pi / sin(pi z),
        Gamma(n z) = (2 pi)^((1-n)/2) n^(n z - 1/2) prod_{j<n} Gamma(z + j/n).

    The lines are found numerically (an integer relation between z_j,
    z_g and 1) and verified by an exact zero test. With w normalized so
    that the a_j are coprime integers, the unknowns are
    x(a, k) = log Gamma(a w + k/N) for a in S (the divisors of the a_j,
    with both signs) and 0 <= k < N; the constants are log sin(pi(a w + c)),
    log(a w + c) (field elements, from shifts), log p, w log p (from
    n^(n a w)) and log pi. As for rational arguments, a reduced row
    echelon form with the unknowns first, then the generators (latest
    first), gives the relations between generators and constants, made
    integral by an HNF. Since the additive relations hold only modulo
    2 pi i, the elimination uses their multiplicative form and the
    result is checked numerically against the eliminated generator.
    The constants enter the tower through exp(pi i w) (the sines),
    exp(w log p / d), log p, pi, roots of unity and root generators.
*/

#define GL_SIN 0
#define GL_LIN 1
#define GL_PRIME 2
#define GL_WLOG 3
#define GL_PI 4
#define GL_HLIN 5       /* polygamma lines: (a w + c)^(-s) */
#define GL_HCOT 6       /* zeta(s, a w + c) + (-1)^s zeta(s, 1 - a w - c) */


static void
_qvec_set(fmpq * r, const fmpq * x, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        fmpq_set(r + i, x + i);
}

static void
_qvec_zero(fmpq * r, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        fmpq_zero(r + i);
}
#define GAMMA_LINE_UNKNOWNS_LIMIT 600

/*
    By the multiplication formula (and the reflection formula for a < 0),
    Gamma(a w + b) is a constant times the product of the Gamma(w + r) for
    r in R = {(b + i)/a mod 1 : 0 <= i < |a|} (likewise the Hurwitz zeta
    function, with a sum), and the systems of the lines relate the values
    at w + r for distinct r mod 1 only through such constants. A generator
    with a residue in no other R_j thus takes no part in a relation; this
    returns 1 when every generator has one (there is then no relation:
    two values on a line from unrelated computations of a long session,
    say Gamma(4w) and Gamma(3w - 19/8), would otherwise cost a large
    reduced row echelon form). (A relation derived from the rows of a
    system maps to a cancellation of the formal sums of the Gamma(w + r),
    with the sign of a_j, since every row does; the residues need not be
    multiples of 1/N, the unknowns of the system at other residues being
    then unrelated to the others.)
*/
static int
_line_all_private(slong ns, const slong * a, const fmpq * b)
{
    slong * cnt, * res, * off, i, j, k, tot = 0, N = 1;
    int all_private = 1;
    fmpq_t t;
    fmpz_t u;

    /* the residues in (1/N) Z / Z */
    for (j = 0; j < ns; j++)
    {
        slong q;
        if (!fmpz_fits_si(fmpq_denref(b + j)))
            return 0;
        q = fmpz_get_si(fmpq_denref(b + j)) * FLINT_ABS(a[j]);
        N = (N / n_gcd(N, q)) * q;
        if (N > 1000000)
            return 0;
        tot += FLINT_ABS(a[j]);
    }

    cnt = flint_calloc(N, sizeof(slong));
    res = flint_malloc(sizeof(slong) * FLINT_MAX(tot, 1));
    off = flint_malloc(sizeof(slong) * (ns + 1));
    fmpq_init(t);
    fmpz_init(u);

    off[0] = 0;
    for (j = 0; j < ns && all_private; j++)
    {
        slong aa = FLINT_ABS(a[j]);
        off[j + 1] = off[j] + aa;
        for (i = 0; i < aa; i++)
        {
            fmpq_add_si(t, b + j, i);
            fmpq_mul_si(t, t, N);
            fmpz_mul_si(fmpq_denref(t), fmpq_denref(t), a[j]);
            fmpq_canonicalise(t);
            if (!fmpz_is_one(fmpq_denref(t)))
            {
                all_private = 0;    /* (not reached) */
                break;
            }
            fmpz_mod_ui(u, fmpq_numref(t), N);
            res[off[j] + i] = fmpz_get_si(u);
        }
        /* each residue once per generator */
        for (i = 0; i < aa && all_private; i++)
        {
            for (k = 0; k < i; k++)
                if (res[off[j] + k] == res[off[j] + i])
                    break;
            if (k == i)
                cnt[res[off[j] + i]]++;
        }
    }

    for (j = 0; j < ns && all_private; j++)
    {
        int private_res = 0;
        for (i = off[j]; i < off[j + 1] && !private_res; i++)
            if (cnt[res[i]] == 1)
                private_res = 1;
        if (!private_res)
            all_private = 0;
    }

    fmpq_clear(t);
    fmpz_clear(u);
    flint_free(cnt);
    flint_free(res);
    flint_free(off);
    return all_private;
}

/*
    The common part of the systems of a line (gamma and Hurwitz zeta):
    the multipliers Sa (the divisors of the |a_j|, with both signs) and
    the level N (the denominators of the b_j and the ratios of the
    multipliers). Returns 0 (Sa not allocated) when the system would be
    too large, or has no relation (_line_all_private).
*/
static int
_line_setup(slong ** Sa_out, slong * nSa_out, slong * N_out, slong ns, const slong * a, const fmpq * b, gr_tower_flat_t F)
{
    slong limit = GR_TOWER_OPTION(F->T, GR_TOWER_OPT_GAMMA_LINE_LIMIT);
    slong * Sa, nSa = 0, N = 1, i, j, n;

    for (j = 0; j < ns; j++)
        if (FLINT_ABS(a[j]) > limit || a[j] == 0 ||
            !fmpz_fits_si(fmpq_denref(b + j)) || fmpz_cmp_ui(fmpq_denref(b + j), 1000) > 0)
            return 0;

    Sa = flint_malloc(sizeof(slong) * 4 * limit);
    for (j = 0; j < ns; j++)
    {
        slong aa = FLINT_ABS(a[j]), dd, q;
        for (dd = 1; dd <= aa; dd++)
        {
            if (aa % dd == 0)
            {
                for (i = 0; i < nSa; i++)
                    if (Sa[i] == dd)
                        break;
                if (i == nSa)
                {
                    Sa[nSa++] = dd;
                    Sa[nSa++] = -dd;
                }
            }
        }
        q = fmpz_get_si(fmpq_denref(b + j));
        N = (N / n_gcd(N, q)) * q;
    }
    for (i = 0; i < nSa; i++)
        for (j = 0; j < nSa; j++)
            if (Sa[i] > 0 && Sa[j] > Sa[i] && Sa[j] % Sa[i] == 0)
            {
                n = Sa[j] / Sa[i];
                N = (N / n_gcd(N, n)) * n;
            }

    if (nSa * N > GAMMA_LINE_UNKNOWNS_LIMIT || _line_all_private(ns, a, b))
    {
        flint_free(Sa);
        return 0;
    }

    *Sa_out = Sa;
    *nSa_out = nSa;
    *N_out = N;
    return 1;
}

/* b = mm + kk/N with 0 <= kk < N, and the index ia of the multiplier aj */
static void
_line_position(slong * mm, slong * kk, slong * ia, const fmpq_t b, slong aj, slong N, const slong * Sa, slong nSa)
{
    fmpz_t fl;
    fmpq_t beta;
    fmpz_init(fl);
    fmpq_init(beta);
    fmpz_fdiv_q(fl, fmpq_numref(b), fmpq_denref(b));
    fmpq_sub_fmpz(beta, b, fl);
    *mm = fmpz_get_si(fl);
    *kk = fmpz_get_si(fmpq_numref(beta)) * (N / fmpz_get_si(fmpq_denref(beta)));
    for (*ia = 0; *ia < nSa; (*ia)++)
        if (Sa[*ia] == aj)
            break;
    fmpz_clear(fl);
    fmpq_clear(beta);
}

typedef struct
{
    int type;
    slong a;        /* SIN, LIN */
    fmpq_t c;       /* SIN, LIN: a w + c */
    slong p;        /* PRIME, WLOG */
}
gl_const_struct;

typedef struct
{
    slong ncols0;               /* unknowns + generators */
    slong nconst;
    slong alloc_const;
    gl_const_struct * consts;
    slong nrows;
    slong alloc_rows;
    fmpq ** rows0;              /* dense parts over the unknowns and generators */
    fmpq ** rowsc;              /* constant parts (length alloc_const at creation; grown) */
    slong * rowsc_len;
    fmpq * cur0;                /* the row being built */
    fmpq * curc;
    slong curc_alloc;
}
gl_system_struct;

static void
_gl_init(gl_system_struct * S, slong ncols0)
{
    S->ncols0 = ncols0;
    S->nconst = 0;
    S->alloc_const = 0;
    S->consts = NULL;
    S->nrows = 0;
    S->alloc_rows = 0;
    S->rows0 = NULL;
    S->rowsc = NULL;
    S->rowsc_len = NULL;
    S->cur0 = _fmpq_vec_init(ncols0);
    S->curc_alloc = 64;
    S->curc = _fmpq_vec_init(S->curc_alloc);
}

static void
_gl_clear(gl_system_struct * S)
{
    slong i;
    for (i = 0; i < S->nconst; i++)
        fmpq_clear(S->consts[i].c);
    flint_free(S->consts);
    for (i = 0; i < S->nrows; i++)
    {
        _fmpq_vec_clear(S->rows0[i], S->ncols0);
        _fmpq_vec_clear(S->rowsc[i], S->rowsc_len[i]);
    }
    flint_free(S->rows0);
    flint_free(S->rowsc);
    flint_free(S->rowsc_len);
    _fmpq_vec_clear(S->cur0, S->ncols0);
    _fmpq_vec_clear(S->curc, S->curc_alloc);
}

static slong
_gl_const(gl_system_struct * S, int type, slong a, const fmpq_t c, slong p)
{
    slong i;

    for (i = 0; i < S->nconst; i++)
    {
        gl_const_struct * k = S->consts + i;
        if (k->type != type)
            continue;
        if ((type == GL_SIN || type == GL_LIN || type == GL_HLIN || type == GL_HCOT) && k->a == a && fmpq_equal(k->c, c))
            return i;
        if ((type == GL_PRIME || type == GL_WLOG) && k->p == p)
            return i;
        if (type == GL_PI)
            return i;
    }

    if (S->nconst == S->alloc_const)
    {
        S->alloc_const = FLINT_MAX(16, 2 * S->alloc_const);
        S->consts = flint_realloc(S->consts, sizeof(gl_const_struct) * S->alloc_const);
    }
    if (S->nconst >= S->curc_alloc)
    {
        fmpq * t = _fmpq_vec_init(2 * S->curc_alloc);
        _qvec_set(t, S->curc, S->curc_alloc);
        _fmpq_vec_clear(S->curc, S->curc_alloc);
        S->curc = t;
        S->curc_alloc *= 2;
    }

    S->consts[S->nconst].type = type;
    S->consts[S->nconst].a = a;
    fmpq_init(S->consts[S->nconst].c);
    if (c != NULL)
        fmpq_set(S->consts[S->nconst].c, c);
    S->consts[S->nconst].p = p;
    return S->nconst++;
}

static void
_gl_add0(gl_system_struct * S, slong col, slong num, slong den)
{
    fmpq_t t;
    fmpq_init(t);
    fmpq_set_si(t, num, den);
    fmpq_add(S->cur0 + col, S->cur0 + col, t);
    fmpq_clear(t);
}

static void
_gl_addc(gl_system_struct * S, int type, slong a, const fmpq_t c, slong p, const fmpq_t v)
{
    slong i = _gl_const(S, type, a, c, p);
    fmpq_add(S->curc + i, S->curc + i, v);
}

static void
_gl_addc_si(gl_system_struct * S, int type, slong a, const fmpq_t c, slong p, slong num, slong den)
{
    fmpq_t v;
    fmpq_init(v);
    fmpq_set_si(v, num, den);
    _gl_addc(S, type, a, c, p, v);
    fmpq_clear(v);
}

/* log(a w + c) for c = k/N + t */
static void
_gl_add_lin(gl_system_struct * S, slong a, slong k, slong N, slong t, slong sign)
{
    fmpq_t c;
    fmpq_init(c);
    fmpq_set_si(c, k + t * N, N);
    _gl_addc_si(S, GL_LIN, a, c, 0, sign, 1);
    fmpq_clear(c);
}

static void
_gl_finish_row(gl_system_struct * S)
{
    slong i;
    int nonzero = 0;

    for (i = 0; i < S->ncols0 && !nonzero; i++)
        nonzero = !fmpq_is_zero(S->cur0 + i);
    for (i = 0; i < S->nconst && !nonzero; i++)
        nonzero = !fmpq_is_zero(S->curc + i);

    if (nonzero)
    {
        if (S->nrows == S->alloc_rows)
        {
            S->alloc_rows = FLINT_MAX(16, 2 * S->alloc_rows);
            S->rows0 = flint_realloc(S->rows0, sizeof(fmpq *) * S->alloc_rows);
            S->rowsc = flint_realloc(S->rowsc, sizeof(fmpq *) * S->alloc_rows);
            S->rowsc_len = flint_realloc(S->rowsc_len, sizeof(slong) * S->alloc_rows);
        }
        S->rows0[S->nrows] = _fmpq_vec_init(S->ncols0);
        _qvec_set(S->rows0[S->nrows], S->cur0, S->ncols0);
        S->rowsc_len[S->nrows] = FLINT_MAX(S->nconst, 1);
        S->rowsc[S->nrows] = _fmpq_vec_init(S->rowsc_len[S->nrows]);
        _qvec_set(S->rowsc[S->nrows], S->curc, S->nconst);
        S->nrows++;
    }

    _qvec_zero(S->cur0, S->ncols0);
    _qvec_zero(S->curc, S->curc_alloc);
}

/* the matrix of the system (A initialised here), the constants in the
   columns cperm (NULL: in their order) after the others */
static void
_gl_matrix(fmpq_mat_t A, const gl_system_struct * S, const slong * cperm)
{
    slong r, j, n0 = S->ncols0;
    fmpq_mat_init(A, S->nrows, n0 + S->nconst);
    for (r = 0; r < S->nrows; r++)
    {
        for (j = 0; j < n0; j++)
            fmpq_set(fmpq_mat_entry(A, r, j), S->rows0[r] + j);
        for (j = 0; j < S->nconst && j < S->rowsc_len[r]; j++)
            fmpq_set(fmpq_mat_entry(A, r, n0 + (cperm != NULL ? cperm[j] : j)), S->rowsc[r] + j);
    }
}

/* the flat element i (the generator must be present) */
static int
_gl_flat_i(fmpz_mpoly_q_t res, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong di = _gr_tower_certify_find_i(T);
    const gr_tower_gen_struct * g;

    if (di < 0)
        return 0;
    g = T->gens + di;
    fmpz_mpoly_q_gen(res, GR_TOWER_FLAT_VAR_D(F, di), F->mctx);
    if (g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param != 4)
    {
        fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), g->def_param / 4, F->mctx);
        GR_MUST_SUCCEED(gr_tower_flat_reduce(res, F));
    }
    return 1;
}

/*
    A generator of kind (EXP or LOG, or one which became algebraic with
    that definition) with argument u and definition order < before:
    returns 0 and sets *d; or adjoins it (moved right after the
    generators of u) and returns 1; -1 if impossible.
*/
static int
_gl_trans_gen(slong * d_out, gr_tower_flat_t F, int kind, const fmpz_mpoly_q_t u, slong before)
{
    gr_tower_struct * T = F->T;
    slong d, dep;
    int st;

    for (d = 0; d < T->num_gens; d++)
    {
        const gr_tower_gen_struct * g = T->gens + d;
        fmpz_mpoly_q_t diff;
        int eq;

        if (g->arg.mctx == NULL || !(g->kind == kind || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == kind)))
            continue;

        fmpz_mpoly_q_init(diff, F->mctx);
        gr_tower_flat_convert(diff, &g->arg.data, g->arg.mctx, F);
        fmpz_mpoly_q_sub(diff, diff, u, F->mctx);
        eq = (gr_tower_flat_reduce(diff, F) == GR_SUCCESS && fmpz_mpoly_q_is_zero(diff, F->mctx));
        fmpz_mpoly_q_clear(diff, F->mctx);

        if (eq)
        {
            if (d < before)
            {
                *d_out = d;
                return 0;
            }
            /* an existing transcendental generator after `before` is moved
               right after the generators of its argument (the generators
               between do not matter: only those depending on it must come
               after it, and they do) */
            dep = _gr_tower_flat_max_dep(u, F);
            if (g->kind != kind || dep + 1 >= before)
                return -1;
            _gr_tower_move_gen(T, d, dep + 1);
            return 1;
        }
    }

    dep = _gr_tower_flat_max_dep(u, F);
    if (dep + 1 >= before)
        return -1;

    if (kind == GR_TOWER_EXP)
        st = gr_tower_adjoin_exp_flat(T, u, F->mctx, NULL);
    else
        st = gr_tower_adjoin_log_flat(T, u, F->mctx, NULL);
    if (st != GR_SUCCESS)
        return -1;

    _gr_tower_move_gen(T, T->num_gens - 1, dep + 1);
    return 1;
}

/* the elimination data of one relation */
typedef struct
{
    slong dtop;
    slong ns;
    const slong * sd;           /* the generators of the line */
    slong * e;                  /* exponents of the other generators (times scale) */
    fmpq * kappa;               /* exponents of the constants (times scale) */
    gl_system_struct * S;
    fmpz_mpoly_q_t w;           /* the base of the line (flat) */
    slong N;
}
gl_elim_struct;

/*
    expr = prod Gamma_j^(e_j) prod C_t^(kappa_t), with the generators
    needed adjoined first (returns 1 then: restart), -1 if impossible or
    an exponent which must be an integer is not, 0 when expr is set.
*/
static int
_gl_build(fmpz_mpoly_q_t expr, gr_tower_flat_t F, gl_elim_struct * E)
{
    gr_tower_struct * T = F->T;
    gl_system_struct * S = E->S;
    fmpz_mpoly_q_t x, y, c, piw;
    slong j, t, st = 0, dpi = -1;
    int need_sin = 0, need_pi = 0;
    ulong m = 1;

    for (t = 0; t < S->nconst; t++)
    {
        if (fmpq_is_zero(E->kappa + t))
            continue;
        if ((S->consts[t].type == GL_SIN || S->consts[t].type == GL_LIN) &&
            !fmpz_is_one(fmpq_denref(E->kappa + t)))
            return -1;
        if (!fmpz_fits_si(fmpq_numref(E->kappa + t)) || !fmpz_fits_si(fmpq_denref(E->kappa + t)) ||
            fmpz_bits(fmpq_numref(E->kappa + t)) > 10)
            return -1;
        if (S->consts[t].type == GL_SIN)
            need_sin = 1;
        if (S->consts[t].type == GL_PI)
            need_pi = 1;
    }

    /* pi, i and the roots of unity zeta_{2N} (with i) */
    if (need_sin)
    {
        need_pi = 1;
        m = 2 * E->N;
        if (m % 4 != 0)
            m *= 2;
    }
    st = _gamma_prepare_constants(T, need_sin ? m : 4, need_pi, E->dtop);
    if (st != 0)
        return st;

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(x, F->mctx);
    fmpz_mpoly_q_init(y, F->mctx);
    fmpz_mpoly_q_init(c, F->mctx);
    fmpz_mpoly_q_init(piw, F->mctx);

    if (need_pi)
        dpi = _gr_tower_certify_find_pi(T, T->num_gens);

    /* the transcendental constants: log p, exp(w log p / d), exp(pi i w) */
    for (t = 0; t < S->nconst && st == 0; t++)
    {
        gl_const_struct * k = S->consts + t;
        slong dl, de;

        if (fmpq_is_zero(E->kappa + t) || k->type != GL_WLOG)
            continue;

        fmpz_mpoly_q_set_si(c, k->p, F->mctx);
        /* (log p depends on nothing: adjoined at the front) */
        st = _gl_trans_gen(&dl, F, GR_TOWER_LOG, c, T->num_gens);
        if (st != 0)
            break;

        /* exp(w log p / d) */
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, dl), F->mctx);
        gr_tower_flat_convert(c, E->w, F->mctx, F);
        fmpz_mpoly_q_mul(x, x, c, F->mctx);
        fmpz_mpoly_q_div_fmpz(x, x, fmpq_denref(E->kappa + t), F->mctx);
        GR_MUST_SUCCEED(gr_tower_flat_reduce(x, F));
        st = _gl_trans_gen(&de, F, GR_TOWER_EXP, x, E->dtop);
        if (st != 0)
            break;
        fmpz_mpoly_q_gen(y, GR_TOWER_FLAT_VAR_D(F, de), F->mctx);
        _gr_tower_certify_mul_pow_si(expr, y, fmpz_get_si(fmpq_numref(E->kappa + t)), F);
    }

    if (st == 0 && need_sin)
    {
        /* E = exp(pi i w) */
        slong de;
        gr_tower_flat_ensure(F);
        if (!_gl_flat_i(x, F))
            st = -1;
        else
        {
            fmpz_mpoly_q_gen(y, GR_TOWER_FLAT_VAR_D(F, dpi), F->mctx);
            fmpz_mpoly_q_mul(x, x, y, F->mctx);
            gr_tower_flat_convert(c, E->w, F->mctx, F);
            fmpz_mpoly_q_mul(piw, x, c, F->mctx);
            GR_MUST_SUCCEED(gr_tower_flat_reduce(piw, F));
            st = _gl_trans_gen(&de, F, GR_TOWER_EXP, piw, E->dtop);
        }

        if (st == 0)
        {
            fmpz_mpoly_q_t Eg, zi, z;
            fmpz_mpoly_q_init(Eg, F->mctx);
            fmpz_mpoly_q_init(zi, F->mctx);
            fmpz_mpoly_q_init(z, F->mctx);
            fmpz_mpoly_q_gen(Eg, GR_TOWER_FLAT_VAR_D(F, de), F->mctx);

            if (!_flat_zeta(zi, m, F))
                st = -1;

            for (t = 0; t < S->nconst && st == 0; t++)
            {
                gl_const_struct * k = S->consts + t;
                slong kk, e2;
                fmpz_t num;

                if (fmpq_is_zero(E->kappa + t) || k->type != GL_SIN)
                    continue;

                /* sin(pi(a w + c)) = (E^a z - E^(-a) z^(-1)) / (2i), z = e^(pi i c),
                   c = kk/N: z = zeta_m^(kk m / (2N)) */
                fmpz_init(num);
                fmpz_mul_si(num, fmpq_numref(k->c), (slong) (m / 2));
                fmpz_divexact(num, num, fmpq_denref(k->c));
                kk = fmpz_get_si(num);    /* c m / 2 */
                fmpz_clear(num);
                kk = ((kk % (slong) m) + (slong) m) % (slong) m;

                fmpz_mpoly_q_one(x, F->mctx);
                _gr_tower_certify_mul_pow_si(x, Eg, k->a, F);
                fmpz_mpoly_q_one(z, F->mctx);
                _gr_tower_certify_mul_pow_si(z, zi, kk, F);
                fmpz_mpoly_q_mul(x, x, z, F->mctx);

                fmpz_mpoly_q_one(y, F->mctx);
                _gr_tower_certify_mul_pow_si(y, Eg, -k->a, F);
                fmpz_mpoly_q_one(z, F->mctx);
                _gr_tower_certify_mul_pow_si(z, zi, ((slong) m - kk) % (slong) m, F);
                fmpz_mpoly_q_mul(y, y, z, F->mctx);

                fmpz_mpoly_q_sub(x, x, y, F->mctx);

                /* / (2i) = * (-i/2) */
                fmpz_mpoly_q_one(z, F->mctx);
                _gr_tower_certify_mul_pow_si(z, zi, 3 * (m / 4), F);
                fmpz_mpoly_q_mul(x, x, z, F->mctx);
                fmpz_mpoly_q_div_si(x, x, 2, F->mctx);
                GR_MUST_SUCCEED(gr_tower_flat_reduce(x, F));

                e2 = fmpz_get_si(fmpq_numref(E->kappa + t));
                _gr_tower_certify_mul_pow_si(expr, x, e2, F);
                GR_MUST_SUCCEED(gr_tower_flat_reduce(expr, F));
            }

            fmpz_mpoly_q_clear(Eg, F->mctx);
            fmpz_mpoly_q_clear(zi, F->mctx);
            fmpz_mpoly_q_clear(z, F->mctx);
        }
    }

    /* field elements a w + c */
    for (t = 0; t < S->nconst && st == 0; t++)
    {
        gl_const_struct * k = S->consts + t;
        if (fmpq_is_zero(E->kappa + t) || k->type != GL_LIN)
            continue;
        gr_tower_flat_convert(x, E->w, F->mctx, F);
        fmpz_mpoly_q_mul_si(x, x, k->a, F->mctx);
        fmpz_mpoly_q_add_fmpq(x, x, k->c, F->mctx);
        _gr_tower_certify_mul_pow_si(expr, x, fmpz_get_si(fmpq_numref(E->kappa + t)), F);
    }

    /* pi and primes, through root generators */
    for (t = 0; t < S->nconst && st == 0; t++)
    {
        gl_const_struct * k = S->consts + t;
        if (fmpq_is_zero(E->kappa + t) || (k->type != GL_PI && k->type != GL_PRIME))
            continue;
        if (k->type == GL_PI)
            fmpz_mpoly_q_gen(c, GR_TOWER_FLAT_VAR_D(F, dpi), F->mctx);
        else
            fmpz_mpoly_q_set_si(c, k->p, F->mctx);
        st = _mul_root_power(expr, F, c, fmpz_get_si(fmpq_numref(E->kappa + t)),
                fmpz_get_si(fmpq_denref(E->kappa + t)), E->dtop);
    }

    /* the other gamma generators */
    for (j = 0; j < E->ns && st == 0; j++)
    {
        if (E->sd[j] == E->dtop || E->e[j] == 0)
            continue;
        fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, E->sd[j]), F->mctx);
        _gr_tower_certify_mul_pow_si(expr, x, E->e[j], F);
    }

    if (st == 0 && gr_tower_flat_reduce(expr, F) != GR_SUCCESS)
        st = -1;

    fmpz_mpoly_q_clear(x, F->mctx);
    fmpz_mpoly_q_clear(y, F->mctx);
    fmpz_mpoly_q_clear(c, F->mctx);
    fmpz_mpoly_q_clear(piw, F->mctx);
    return st;
}

/*
    The relations on the line of the generators sd[0..ns) (latest first),
    with z_j = a_j w + b_j: eliminates one generator (returns 1, also
    when the tower changed otherwise), or returns 0.
*/
static int
_gamma_line_level(gr_tower_flat_t F, slong ns, const slong * sd, const slong * a, const fmpq * b,
    const fmpz_mpoly_q_t w, slong depth)
{
    slong * Sa, nSa = 0, N = 1, i, j, k, n, ia, nX, nG, ncols, row = -1, top = -1, rank, r;
    gl_system_struct S;
    fmpq_mat_t A, B;
    int changed = 0;

    (void) depth;

    if (!_line_setup(&Sa, &nSa, &N, ns, a, b, F))
        return 0;

    nX = nSa * N;
    nG = ns;
    _gl_init(&S, nX + nG);

#define XCOL(ia, kk) ((ia) * N + (kk))

    /* generators: G_j = x(a_j, beta_j) + shifts */
    for (j = 0; j < ns; j++)
    {
        slong mm, kk;
        _line_position(&mm, &kk, &ia, b + j, a[j], N, Sa, nSa);
        _gl_add0(&S, nX + j, 1, 1);
        _gl_add0(&S, XCOL(ia, kk), -1, 1);
        if (mm > 0)
            for (i = 0; i < mm; i++)
                _gl_add_lin(&S, a[j], kk, N, i, -1);
        else
            for (i = 0; i < -mm; i++)
                _gl_add_lin(&S, a[j], kk, N, i - (-mm), 1);
        _gl_finish_row(&S);
    }

    /* reflection: x(a, k) + x(-a, N - k) (+ log(-a w) for k = 0)
       - log pi + log sin(pi(a w + k/N)) = 0 */
    for (ia = 0; ia < nSa; ia++)
    {
        slong ib;
        if (Sa[ia] < 0)
            continue;
        for (ib = 0; ib < nSa; ib++)
            if (Sa[ib] == -Sa[ia])
                break;
        for (k = 0; k < N; k++)
        {
            fmpq_t c;
            fmpq_init(c);
            fmpq_set_si(c, k, N);
            _gl_add0(&S, XCOL(ia, k), 1, 1);
            if (k == 0)
            {
                _gl_add0(&S, XCOL(ib, 0), 1, 1);
                _gl_add_lin(&S, -Sa[ia], 0, N, 0, 1);
            }
            else
                _gl_add0(&S, XCOL(ib, N - k), 1, 1);
            _gl_addc_si(&S, GL_PI, 0, NULL, 0, -1, 1);
            _gl_addc_si(&S, GL_SIN, Sa[ia], c, 0, 1, 1);
            _gl_finish_row(&S);
            fmpq_clear(c);
        }
    }

    /* multiplication: x(na, r) + shifts - sum_j (x(a, k_j) + shifts)
       - (1-n)/2 (log 2 + log pi) - (n k/N - 1/2) log n - n a w log n = 0 */
    for (ia = 0; ia < nSa; ia++)
    {
        for (j = 0; j < nSa; j++)
        {
            slong ib = j, aa = Sa[ia], qq, rr, jj, p, mfac;
            if (Sa[ib] == aa || Sa[ib] % aa != 0 || Sa[ib] / aa < 2)
                continue;
            n = Sa[ib] / aa;
            if (N % n != 0)
                continue;
            for (k = 0; k < N; k++)
            {
                qq = (n * k) / N;
                rr = (n * k) % N;
                _gl_add0(&S, XCOL(ib, rr), 1, 1);
                for (i = 0; i < qq; i++)
                    _gl_add_lin(&S, Sa[ib], rr, N, i, 1);
                for (jj = 0; jj < n; jj++)
                {
                    slong kj = k + jj * (N / n);
                    if (kj >= N)
                    {
                        kj -= N;
                        _gl_add0(&S, XCOL(ia, kj), -1, 1);
                        _gl_add_lin(&S, aa, kj, N, 0, -1);
                    }
                    else
                        _gl_add0(&S, XCOL(ia, kj), -1, 1);
                }
                _gl_addc_si(&S, GL_PRIME, 0, NULL, 2, -(1 - n), 2);
                _gl_addc_si(&S, GL_PI, 0, NULL, 0, -(1 - n), 2);
                mfac = n;
                for (p = 2; mfac > 1; p++)
                {
                    slong v = 0;
                    while (mfac % p == 0)
                    {
                        mfac /= p;
                        v++;
                    }
                    if (v != 0)
                    {
                        /* -(n k/N - 1/2) v log p = (N - 2 n k) v / (2N) log p */
                        _gl_addc_si(&S, GL_PRIME, 0, NULL, p, (N - 2 * n * k) * v, 2 * N);
                        _gl_addc_si(&S, GL_WLOG, 0, NULL, p, -n * aa * v, 1);
                    }
                }
                _gl_finish_row(&S);
            }
        }
    }

#undef XCOL

    /* the matrix, constants ordered: sines, linear factors, primes,
       w log p, pi */
    {
        slong * cperm = flint_malloc(sizeof(slong) * FLINT_MAX(S.nconst, 1));
        slong c = 0, ty;
        for (ty = GL_SIN; ty <= GL_PI; ty++)
            for (i = 0; i < S.nconst; i++)
                if (S.consts[i].type == ty)
                    cperm[i] = c++;

        ncols = nX + nG + S.nconst;
        _gl_matrix(A, &S, cperm);

        fmpq_mat_init(B, S.nrows, ncols);
        rank = fmpq_mat_rref(B, A);

        for (r = 0; r < rank && row < 0; r++)
        {
            for (j = 0; j < ncols; j++)
                if (!fmpq_is_zero(fmpq_mat_entry(B, r, j)))
                    break;
            if (j >= nX && j < nX + nG)
            {
                row = r;
                top = j - nX;
            }
            else if (j >= nX + nG)
                break;
        }

        if (row >= 0)
        {
            gl_elim_struct E;
            fmpq_t t;
            fmpz_t P;
            int ok;

            /* (the relations hold modulo 2 pi i: the elimination needs
               the integral combination, whose pivot P bounds the
               ambiguity of the exponentiated relation to a P-th root of
               unity, which the numerical check below excludes) */
            fmpz_init(P);
            ok = _gamma_integral_relation(B, row, A, S.nrows, nX + nG, ncols, nX + top, P, NULL, NULL);

            /* log G_top = - sum r_j log G_j - sum kappa_t C_t */
            fmpq_init(t);
            E.dtop = sd[top];
            E.ns = ns;
            E.sd = sd;
            E.e = flint_calloc(ns, sizeof(slong));
            E.kappa = _fmpq_vec_init(FLINT_MAX(S.nconst, 1));
            E.S = &S;
            E.N = N;
            fmpz_mpoly_q_init(E.w, F->mctx);
            fmpz_mpoly_q_set(E.w, w, F->mctx);

            for (j = 0; j < ns && ok; j++)
            {
                if (j == top)
                    continue;
                fmpq_neg(t, fmpq_mat_entry(B, row, nX + j));
                if (!fmpz_is_one(fmpq_denref(t)) || !fmpz_fits_si(fmpq_numref(t)))
                    ok = 0;
                else
                    E.e[j] = fmpz_get_si(fmpq_numref(t));
            }
            for (j = 0; j < S.nconst; j++)
                fmpq_neg(E.kappa + j, fmpq_mat_entry(B, row, nX + nG + cperm[j]));

            if (ok)
            {
                fmpz_mpoly_q_t expr;
                gr_tower_flat_ensure(F);
                fmpz_mpoly_q_init(expr, F->mctx);
                fmpz_mpoly_q_one(expr, F->mctx);
                {
                    int st = _gl_build(expr, F, &E);
                    if (st == 1)
                        changed = 1;
                    else if (st == 0)
                    {
                        /* the value of expr is that of the generator
                           up to a P-th root of unity: the enclosures
                           decide, with more accurate bits than P has */
                        if (_gr_tower_flat_matches_gen_bits(expr, E.dtop, (slong) fmpz_bits(P) + 3, F) &&
                            _gr_tower_flat_eliminate_gen(E.dtop, expr, F))
                            changed = 1;
                    }
                }
                fmpz_mpoly_q_clear(expr, F->mctx);
            }

            flint_free(E.e);
            _fmpq_vec_clear(E.kappa, FLINT_MAX(S.nconst, 1));
            fmpz_mpoly_q_clear(E.w, F->mctx);
            fmpq_clear(t);
            fmpz_clear(P);
        }

        flint_free(cperm);
        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
    }

    _gl_clear(&S);
    flint_free(Sa);
    return changed;
}

/* -------------------------------------------------------------------- */
/* Hurwitz zeta (polygamma) values on rational affine lines             */
/* -------------------------------------------------------------------- */

/*
    Polygamma generators psi^(m)(z_j) = (-1)^s (s-1)! zeta(s, z_j)
    (s = m + 1 >= 2) with z_j = a_j w + b_j on a common line satisfy,
    with y(a, k) = zeta(s, a w + k/N),

        zeta(s, z + 1) = zeta(s, z) - z^(-s),
        zeta(s, z) + (-1)^s zeta(s, 1 - z) = (-1)^(s-1) pi^s P_{s-1}(cot(pi z)) / (s-1)!,
        sum_{j<n} zeta(s, z + j/n) = n^s zeta(s, n z).

    The digamma function (m = 0) is included with s = 1 and y = -psi
    (the same shift and reflection relations, with P_0(x) = x), the
    multiplication theorem then having the constant term n log n
    (sum_{j<n} psi(z + j/n) = n psi(n z) - n log n), through the
    generators log p.

    These are exact (no 2 pi i ambiguity, no exponential factors), so the
    system is the additive analogue of the gamma lines: unknowns y(a, k),
    the generators, and the constants (a w + c)^(-s) (field elements) and
    the reflection terms E(a, c) (through exp(pi i w), pi, i and roots of
    unity: cot(pi u) = i (Y + 1) (1 + Y + ... + Y^(r-1)) / (Y^r - 1) with
    Y = exp(2 pi i u) = exp(pi i w)^(2a) v, v a root of unity of order r,
    so that Y^r = exp(pi i w)^(2 a r) and no algebraic denominator
    appears). Every relation with a generator pivot is eliminated; the
    result is checked numerically against the generator.
*/

/* the reflection term E(a, c), c = k/N, as a flat element (E_w = exp(pi i w),
   zeta_m with N | m, 4 | m, and pi present) */
static int
_hl_cot_term(fmpz_mpoly_q_t res, gr_tower_flat_t F, slong s, slong a, slong k, slong N, ulong m,
    const fmpz_mpoly_q_t Ew, slong dpi, const fmpz_poly_t P)
{
    fmpz_mpoly_q_t zi, v, Y, Yp, S, x, cot, t;
    slong r = N / n_gcd(N, k), i;
    fmpz_t f;
    int ok;

    fmpz_mpoly_q_init(zi, F->mctx);
    fmpz_mpoly_q_init(v, F->mctx);
    fmpz_mpoly_q_init(Y, F->mctx);
    fmpz_mpoly_q_init(Yp, F->mctx);
    fmpz_mpoly_q_init(S, F->mctx);
    fmpz_mpoly_q_init(x, F->mctx);
    fmpz_mpoly_q_init(cot, F->mctx);
    fmpz_mpoly_q_init(t, F->mctx);
    fmpz_init(f);

    ok = _flat_zeta(zi, m, F);
    if (ok)
    {
        /* v = exp(2 pi i k/N) = zeta_m^((m/N) k), Y = E_w^(2a) v */
        fmpz_mpoly_q_one(v, F->mctx);
        _gr_tower_certify_mul_pow_si(v, zi, (m / N) * k, F);
        GR_MUST_SUCCEED(gr_tower_flat_reduce(v, F));
        fmpz_mpoly_q_set(Y, v, F->mctx);
        _gr_tower_certify_mul_pow_si(Y, Ew, 2 * a, F);
        GR_MUST_SUCCEED(gr_tower_flat_reduce(Y, F));

        /* S = 1 + Y + ... + Y^(r-1) */
        fmpz_mpoly_q_one(S, F->mctx);
        fmpz_mpoly_q_one(Yp, F->mctx);
        for (i = 1; i < r; i++)
        {
            fmpz_mpoly_q_mul(Yp, Yp, Y, F->mctx);
            GR_MUST_SUCCEED(gr_tower_flat_reduce(Yp, F));
            fmpz_mpoly_q_add(S, S, Yp, F->mctx);
        }

        /* cot = i (Y + 1) S / (E_w^(2 a r) - 1) */
        fmpz_mpoly_q_add_si(cot, Y, 1, F->mctx);
        fmpz_mpoly_q_mul(cot, cot, S, F->mctx);
        fmpz_mpoly_q_one(t, F->mctx);
        _gr_tower_certify_mul_pow_si(t, zi, m / 4, F);
        fmpz_mpoly_q_mul(cot, cot, t, F->mctx);
        fmpz_mpoly_q_one(t, F->mctx);
        _gr_tower_certify_mul_pow_si(t, Ew, 2 * a * r, F);
        fmpz_mpoly_q_sub_si(t, t, 1, F->mctx);
        fmpz_mpoly_q_div(cot, cot, t, F->mctx);
        ok = (gr_tower_flat_reduce(cot, F) == GR_SUCCESS);
    }

    /* P(cot) by Horner */
    fmpz_mpoly_q_zero(x, F->mctx);
    for (i = fmpz_poly_degree(P); i >= 0 && ok; i--)
    {
        fmpz_mpoly_q_mul(x, x, cot, F->mctx);
        fmpz_mpoly_q_add_fmpz(x, x, P->coeffs + i, F->mctx);
        ok = (gr_tower_flat_reduce(x, F) == GR_SUCCESS);
    }

    if (ok)
    {
        /* times (-1)^(s-1) pi^s / (s-1)! */
        fmpz_mpoly_q_gen(t, GR_TOWER_FLAT_VAR_D(F, dpi), F->mctx);
        _gr_tower_certify_mul_pow_si(x, t, s, F);
        fmpz_fac_ui(f, s - 1);
        if (s % 2 == 0)
            fmpz_neg(f, f);
        fmpz_mpoly_q_div_fmpz(res, x, f, F->mctx);
        ok = (gr_tower_flat_reduce(res, F) == GR_SUCCESS);
    }

    fmpz_mpoly_q_clear(zi, F->mctx);
    fmpz_mpoly_q_clear(v, F->mctx);
    fmpz_mpoly_q_clear(Y, F->mctx);
    fmpz_mpoly_q_clear(Yp, F->mctx);
    fmpz_mpoly_q_clear(S, F->mctx);
    fmpz_mpoly_q_clear(x, F->mctx);
    fmpz_mpoly_q_clear(cot, F->mctx);
    fmpz_mpoly_q_clear(t, F->mctx);
    fmpz_clear(f);
    return ok;
}

/* the shift terms: + sign * sum (a w + c + i)^(-s) */
static void
_hl_add_lin(gl_system_struct * S, slong a, slong k, slong N, slong i, slong sign)
{
    fmpq_t c;
    fmpq_init(c);
    fmpq_set_si(c, k + i * N, N);
    _gl_addc_si(S, GL_HLIN, a, c, 0, sign, 1);
    fmpq_clear(c);
}

/*
    The relations on the line of the polygamma generators of order
    s - 1 sd[0..ns) (latest first), z_j = a_j w + b_j: eliminates the
    generators with pivots (returns 1, also when the tower changed
    otherwise), or returns 0.
*/
static int
_hurwitz_line_level(gr_tower_flat_t F, slong s, slong ns, const slong * sd, const slong * a, const fmpq * b,
    const fmpz_mpoly_q_t w)
{
    gr_tower_struct * T = F->T;
    slong * Sa, nSa = 0, N = 1, i, j, k, n, ia, nX, nG, ncols, rank, r, nrel = 0;
    slong * rows, * tops;
    gl_system_struct S;
    fmpq_mat_t A, B;
    fmpz_t ns_pow, fac;
    fmpq_t alpha;
    int changed = 0;

    if (!_line_setup(&Sa, &nSa, &N, ns, a, b, F))
        return 0;

    nX = nSa * N;
    nG = ns;
    _gl_init(&S, nX + nG);
    fmpz_init(ns_pow);
    fmpz_init(fac);
    fmpq_init(alpha);

    /* zeta(s, z) = alpha psi^(s-1)(z), alpha = 1 / ((-1)^s (s-1)!) */
    fmpz_fac_ui(fac, s - 1);
    if (s % 2 == 1)
        fmpz_neg(fac, fac);
    fmpq_set_fmpz(alpha, fac);
    fmpq_inv(alpha, alpha);

#define XCOL(ia, kk) ((ia) * N + (kk))

    /* generators: alpha G_j - y(a_j, kk) + shifts = 0 */
    for (j = 0; j < ns; j++)
    {
        slong mm, kk;
        _line_position(&mm, &kk, &ia, b + j, a[j], N, Sa, nSa);
        fmpq_add(S.cur0 + nX + j, S.cur0 + nX + j, alpha);
        _gl_add0(&S, XCOL(ia, kk), -1, 1);
        /* zeta(s, u + mm) = zeta(s, u) - sum_{i<mm} (u + i)^(-s);
           zeta(s, u - n) = zeta(s, u) + sum_{i=1}^{n} (u - i)^(-s) */
        if (mm > 0)
            for (i = 0; i < mm; i++)
                _hl_add_lin(&S, a[j], kk, N, i, 1);
        else
            for (i = 1; i <= -mm; i++)
                _hl_add_lin(&S, a[j], kk, N, -i, -1);
        _gl_finish_row(&S);
    }

    /* reflection: y(a, k) + (-1)^s zeta(s, 1 - a w - k/N) - E(a, k/N) = 0 */
    for (ia = 0; ia < nSa; ia++)
    {
        slong ib, sg = (s % 2 == 0) ? 1 : -1;
        if (Sa[ia] < 0)
            continue;
        for (ib = 0; ib < nSa; ib++)
            if (Sa[ib] == -Sa[ia])
                break;
        for (k = 0; k < N; k++)
        {
            fmpq_t c;
            fmpq_init(c);
            fmpq_set_si(c, k, N);
            _gl_add0(&S, XCOL(ia, k), 1, 1);
            if (k == 0)
            {
                /* zeta(s, -a w + 1) = y(-a, 0) - (-a w)^(-s) */
                _gl_add0(&S, XCOL(ib, 0), sg, 1);
                _hl_add_lin(&S, -Sa[ia], 0, N, 0, -sg);
            }
            else
                _gl_add0(&S, XCOL(ib, N - k), sg, 1);
            _gl_addc_si(&S, GL_HCOT, Sa[ia], c, 0, -1, 1);
            _gl_finish_row(&S);
            fmpq_clear(c);
        }
    }

    /* multiplication: sum_{j<n} zeta(s, a w + (k + j N/n)/N) - n^s zeta(s, n a w + n k/N) = 0 */
    for (ia = 0; ia < nSa; ia++)
    {
        for (j = 0; j < nSa; j++)
        {
            slong ib = j, aa = Sa[ia], qq, rr, jj;
            if (Sa[ib] == aa || Sa[ib] % aa != 0 || Sa[ib] / aa < 2)
                continue;
            n = Sa[ib] / aa;
            if (N % n != 0)
                continue;
            fmpz_set_ui(ns_pow, n);
            fmpz_pow_ui(ns_pow, ns_pow, s);
            for (k = 0; k < N; k++)
            {
                fmpq_t v;
                fmpq_init(v);
                for (jj = 0; jj < n; jj++)
                {
                    slong kj = k + jj * (N / n);
                    if (kj >= N)
                    {
                        /* zeta(s, u + 1) = zeta(s, u) - u^(-s) */
                        kj -= N;
                        _gl_add0(&S, XCOL(ia, kj), 1, 1);
                        _hl_add_lin(&S, aa, kj, N, 0, -1);
                    }
                    else
                        _gl_add0(&S, XCOL(ia, kj), 1, 1);
                }
                qq = (n * k) / N;
                rr = (n * k) % N;
                /* (s = 1, with y = -psi: - n log n) */
                if (s == 1)
                {
                    n_factor_t fac;
                    slong f;
                    n_factor_init(&fac);
                    n_factor(&fac, n, 1);
                    for (f = 0; f < fac.num; f++)
                        _gl_addc_si(&S, GL_PRIME, 0, NULL, fac.p[f], -n * (slong) fac.exp[f], 1);
                }
                /* - n^s (y(na, rr) - sum_{i<qq} (n a w + rr/N + i)^(-s)) */
                fmpq_set_fmpz(v, ns_pow);
                fmpq_sub(S.cur0 + XCOL(ib, rr), S.cur0 + XCOL(ib, rr), v);
                for (i = 0; i < qq; i++)
                {
                    fmpq_t c;
                    fmpq_init(c);
                    fmpq_set_si(c, rr + i * N, N);
                    _gl_addc(&S, GL_HLIN, Sa[ib], c, 0, v);
                    fmpq_clear(c);
                }
                _gl_finish_row(&S);
                fmpq_clear(v);
            }
        }
    }

#undef XCOL

    ncols = nX + nG + S.nconst;
    _gl_matrix(A, &S, NULL);
    fmpq_mat_init(B, S.nrows, ncols);
    rank = fmpq_mat_rref(B, A);

    rows = flint_malloc(sizeof(slong) * FLINT_MAX(rank, 1));
    tops = flint_malloc(sizeof(slong) * FLINT_MAX(rank, 1));
    for (r = 0; r < rank; r++)
    {
        for (j = 0; j < ncols; j++)
            if (!fmpq_is_zero(fmpq_mat_entry(B, r, j)))
                break;
        if (j >= nX && j < nX + nG)
        {
            rows[nrel] = r;
            tops[nrel] = j - nX;
            nrel++;
        }
        else if (j >= nX + nG)
            break;
    }

    /* the constants: pi, zeta_m (m = lcm(N, 4)) and exp(pi i w), placed
       before the earliest generator to be eliminated */
    if (nrel > 0)
    {
        slong dmin = WORD_MAX, de = -1, dpi;
        int any_cot = 0, st;
        ulong m = (N / n_gcd(N, 4)) * 4;
        fmpz_poly_t P;
        fmpz_mpoly_q_t Ew, x;

        for (i = 0; i < nrel; i++)
        {
            dmin = FLINT_MIN(dmin, sd[tops[i]]);
            for (j = 0; j < S.nconst; j++)
                if (S.consts[j].type == GL_HCOT && !fmpq_is_zero(fmpq_mat_entry(B, rows[i], nX + nG + j)))
                    any_cot = 1;
        }

        st = 0;
        if (any_cot)
        {
            st = _gamma_prepare_constants(T, m, 1, dmin);
            if (st == 0)
            {
                /* exp(pi i w) */
                gr_tower_flat_ensure(F);
                fmpz_mpoly_q_init(x, F->mctx);
                if (!_gl_flat_i(x, F))
                    st = -1;
                else
                {
                    fmpz_mpoly_q_t y;
                    fmpz_mpoly_q_init(y, F->mctx);
                    fmpz_mpoly_q_gen(y, GR_TOWER_FLAT_VAR_D(F, _gr_tower_certify_find_pi(T, T->num_gens)), F->mctx);
                    fmpz_mpoly_q_mul(x, x, y, F->mctx);
                    gr_tower_flat_convert(y, w, F->mctx, F);
                    fmpz_mpoly_q_mul(x, x, y, F->mctx);
                    GR_MUST_SUCCEED(gr_tower_flat_reduce(x, F));
                    st = _gl_trans_gen(&de, F, GR_TOWER_EXP, x, dmin);
                    fmpz_mpoly_q_clear(y, F->mctx);
                }
                fmpz_mpoly_q_clear(x, F->mctx);
            }
        }
        /* (s = 1) the logarithms of the primes, likewise */
        for (j = 0; j < S.nconst && st == 0; j++)
        {
            int used = 0;
            if (S.consts[j].type != GL_PRIME)
                continue;
            for (i = 0; i < nrel; i++)
                if (!fmpq_is_zero(fmpq_mat_entry(B, rows[i], nX + nG + j)))
                    used = 1;
            if (used)
            {
                slong dl;
                gr_tower_flat_ensure(F);
                fmpz_mpoly_q_init(x, F->mctx);
                fmpz_mpoly_q_set_si(x, S.consts[j].p, F->mctx);
                st = _gl_trans_gen(&dl, F, GR_TOWER_LOG, x, dmin);
                fmpz_mpoly_q_clear(x, F->mctx);
            }
        }

        if (st != 0)
        {
            changed = (st > 0);
            nrel = 0;
        }

        fmpz_poly_init(P);
        _gr_tower_cot_derivative_poly(P, s - 1);

        for (i = 0; i < nrel; i++)
        {
            slong row = rows[i], top = tops[i], dtop = sd[top];
            fmpq_t bt, cj;
            fmpz_mpoly_q_t expr, y;
            int ok = 1;

            gr_tower_flat_ensure(F);
            fmpq_init(bt);
            fmpq_init(cj);
            fmpz_mpoly_q_init(expr, F->mctx);
            fmpz_mpoly_q_init(y, F->mctx);
            fmpz_mpoly_q_init(x, F->mctx);
            fmpz_mpoly_q_init(Ew, F->mctx);
            dpi = any_cot ? _gr_tower_certify_find_pi(T, T->num_gens) : -1;
            if (any_cot)
                fmpz_mpoly_q_gen(Ew, GR_TOWER_FLAT_VAR_D(F, de), F->mctx);

            /* G_top = -(1/B_top) (sum_{j != top} B_j G_j + sum_t B_t C_t) */
            fmpq_inv(bt, fmpq_mat_entry(B, row, nX + top));
            fmpq_neg(bt, bt);

            fmpz_mpoly_q_zero(expr, F->mctx);
            for (j = 0; j < ns && ok; j++)
            {
                if (j == top || fmpq_is_zero(fmpq_mat_entry(B, row, nX + j)))
                    continue;
                if (sd[j] >= dtop)
                {
                    ok = 0;
                    break;
                }
                fmpq_mul(cj, fmpq_mat_entry(B, row, nX + j), bt);
                fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, sd[j]), F->mctx);
                fmpz_mpoly_q_mul_fmpq(x, x, cj, F->mctx);
                fmpz_mpoly_q_add(expr, expr, x, F->mctx);
            }

            for (j = 0; j < S.nconst && ok; j++)
            {
                gl_const_struct * kc = S.consts + j;
                if (fmpq_is_zero(fmpq_mat_entry(B, row, nX + nG + j)))
                    continue;
                fmpq_mul(cj, fmpq_mat_entry(B, row, nX + nG + j), bt);
                if (kc->type == GL_PRIME)
                {
                    /* log p (present before dmin) */
                    slong dl;
                    fmpz_mpoly_q_set_si(y, kc->p, F->mctx);
                    if (_gl_trans_gen(&dl, F, GR_TOWER_LOG, y, dtop) != 0)
                        ok = 0;
                    else
                        fmpz_mpoly_q_gen(x, GR_TOWER_FLAT_VAR_D(F, dl), F->mctx);
                }
                else if (kc->type == GL_HLIN)
                {
                    /* (a w + c)^(-s) */
                    gr_tower_flat_convert(y, w, F->mctx, F);
                    fmpz_mpoly_q_mul_si(y, y, kc->a, F->mctx);
                    fmpz_mpoly_q_add_fmpq(y, y, kc->c, F->mctx);
                    fmpz_mpoly_q_one(x, F->mctx);
                    _gr_tower_certify_mul_pow_si(x, y, -s, F);
                }
                else
                {
                    /* E(a, c), c = k/N */
                    fmpz_t kk;
                    fmpz_init(kk);
                    fmpz_mul_si(kk, fmpq_numref(kc->c), N);
                    fmpz_divexact(kk, kk, fmpq_denref(kc->c));
                    ok = _hl_cot_term(x, F, s, kc->a, fmpz_get_si(kk), N, m, Ew, dpi, P);
                    fmpz_clear(kk);
                }
                fmpz_mpoly_q_mul_fmpq(x, x, cj, F->mctx);
                fmpz_mpoly_q_add(expr, expr, x, F->mctx);
            }

            if (ok)
                ok = (gr_tower_flat_reduce(expr, F) == GR_SUCCESS);

            /* safeguard: the value of expr is that of the generator */
            if (ok)
                ok = _gr_tower_flat_matches_gen(expr, dtop, F);

            if (ok && _gr_tower_flat_eliminate_gen(dtop, expr, F))
                changed = 1;

            fmpq_clear(bt);
            fmpq_clear(cj);
            fmpz_mpoly_q_clear(expr, F->mctx);
            fmpz_mpoly_q_clear(y, F->mctx);
            fmpz_mpoly_q_clear(x, F->mctx);
            fmpz_mpoly_q_clear(Ew, F->mctx);
        }

        fmpz_poly_clear(P);
    }

    flint_free(rows);
    flint_free(tops);
    fmpq_mat_clear(A);
    fmpq_mat_clear(B);
    _gl_clear(&S);
    fmpz_clear(ns_pow);
    fmpz_clear(fac);
    fmpq_clear(alpha);
    flint_free(Sa);
    return changed;
}

/* the argument of a generator as a flat element (0 if none) */
static int
_gen_arg(fmpz_mpoly_q_t z, const gr_tower_gen_struct * g, gr_tower_flat_t F)
{
    if (g->arg.mctx == NULL)
        return 0;
    gr_tower_flat_convert(z, &g->arg.data, g->arg.mctx, F);
    return 1;
}

/*
    The lines: for each generator g of the given kind (gamma, or
    polygamma of order param) at an irrational argument (latest first),
    the other such generators whose arguments are c z_g + d (c, d
    rational: an integer relation found numerically, verified exactly).
*/
static int
_line_round(gr_tower_flat_t F, slong limit, slong depth, int kind, slong param)
{
    gr_tower_struct * T = F->T;
    slong * cand, nc = 0, t, j;
    int changed = 0;
    ulong version = T->version, structure_version = T->structure_version;

    cand = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));
    for (j = limit - 1; j >= 0; j--)
    {
        const gr_tower_gen_struct * g = T->gens + j;
        if (g->kind == kind && g->arg.mctx != NULL && (kind != GR_TOWER_POLYGAMMA || g->def_param == param))
        {
            /* (rational arguments are not on a line) */
            if (!(fmpz_mpoly_q_is_fmpq(&g->arg.data, g->arg.mctx)))
                cand[nc++] = j;
        }
    }

    for (t = 0; t < nc && !changed; t++)
    {
        slong * sd = flint_malloc(sizeof(slong) * nc);
        fmpq * cc = _fmpq_vec_init(nc);
        fmpq * dd = _fmpq_vec_init(nc);
        slong ns = 0;
        fmpz_mpoly_q_t zg, zj, x;
        acb_ptr vec = _acb_vec_init(3);
        slong prec = 128;

        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(zg, F->mctx);
        fmpz_mpoly_q_init(zj, F->mctx);
        fmpz_mpoly_q_init(x, F->mctx);

        _gen_arg(zg, T->gens + cand[t], F);
        sd[ns] = cand[t];
        fmpq_one(cc + ns);
        fmpq_zero(dd + ns);
        ns++;

        if (gr_tower_flat_get_acb(vec + 1, zg, prec, F) == GR_SUCCESS)
        {
            acb_one(vec + 2);
            for (j = 0; j < nc && !changed; j++)
            {
                fmpz_mat_t rel;
                slong nr, i;

                if (j == t)
                    continue;
                _gen_arg(zj, T->gens + cand[j], F);
                if (gr_tower_flat_get_acb(vec + 0, zj, prec, F) != GR_SUCCESS)
                    continue;

                nr = _gr_tower_certify_lindep_all(rel, vec, 3, prec);
                for (i = 0; i < nr; i++)
                {
                    const fmpz * rr = fmpz_mat_row(rel, i);
                    truth_t z;

                    if (fmpz_is_zero(rr + 0) || fmpz_is_zero(rr + 1) ||
                        fmpz_bits(rr + 0) > 12 || fmpz_bits(rr + 1) > 12 || fmpz_bits(rr + 2) > 12)
                        continue;

                    /* r0 z_j + r1 z_g + r2 = 0, exactly */
                    fmpz_mpoly_q_mul_fmpz(x, zj, rr + 0, F->mctx);
                    {
                        fmpz_mpoly_q_t y;
                        fmpz_mpoly_q_init(y, F->mctx);
                        fmpz_mpoly_q_mul_fmpz(y, zg, rr + 1, F->mctx);
                        fmpz_mpoly_q_add(x, x, y, F->mctx);
                        fmpz_mpoly_q_add_fmpz(x, x, rr + 2, F->mctx);
                        fmpz_mpoly_q_clear(y, F->mctx);
                    }
                    z = _gr_tower_decide_zero_flat(x, F, FLINT_MAX(cand[t], cand[j]), depth + 1);
                    if (T->version != version || T->structure_version != structure_version)
                    {
                        changed = 1;
                        break;
                    }
                    if (z == T_TRUE)
                    {
                        sd[ns] = cand[j];
                        fmpq_set_fmpz_frac(cc + ns, rr + 1, rr + 0);
                        fmpq_neg(cc + ns, cc + ns);
                        fmpq_set_fmpz_frac(dd + ns, rr + 2, rr + 0);
                        fmpq_neg(dd + ns, dd + ns);
                        ns++;
                        break;
                    }
                }
                fmpz_mat_clear(rel);
            }
        }

        if (!changed && ns >= 2)
        {
            /* w = z_g G / D with integer coprime a_j = c_j D / G */
            fmpz_t D, G;
            slong * a = flint_malloc(sizeof(slong) * ns);
            int ok = 1;
            fmpz_init(D);
            fmpz_init(G);
            fmpz_one(D);
            for (j = 0; j < ns; j++)
                fmpz_lcm(D, D, fmpq_denref(cc + j));
            fmpz_zero(G);
            for (j = 0; j < ns; j++)
            {
                fmpz_t u;
                fmpz_init(u);
                fmpz_mul(u, fmpq_numref(cc + j), D);
                fmpz_divexact(u, u, fmpq_denref(cc + j));
                fmpz_gcd(G, G, u);
                fmpz_clear(u);
            }
            for (j = 0; j < ns && ok; j++)
            {
                fmpz_t u;
                fmpz_init(u);
                fmpz_mul(u, fmpq_numref(cc + j), D);
                fmpz_divexact(u, u, fmpq_denref(cc + j));
                fmpz_divexact(u, u, G);
                if (!fmpz_fits_si(u) || FLINT_ABS(fmpz_get_si(u)) > GR_TOWER_OPTION(F->T, GR_TOWER_OPT_GAMMA_LINE_LIMIT))
                    ok = 0;
                else
                    a[j] = fmpz_get_si(u);
                fmpz_clear(u);
            }
            if (ok)
            {
                fmpz_mpoly_q_mul_fmpz(x, zg, G, F->mctx);
                fmpz_mpoly_q_div_fmpz(x, x, D, F->mctx);
                if (kind == GR_TOWER_GAMMA)
                    changed = _gamma_line_level(F, ns, sd, a, dd, x, depth);
                else
                    changed = _hurwitz_line_level(F, param + 1, ns, sd, a, dd, x);
            }
            flint_free(a);
            fmpz_clear(D);
            fmpz_clear(G);
        }

        fmpz_mpoly_q_clear(zg, F->mctx);
        fmpz_mpoly_q_clear(zj, F->mctx);
        fmpz_mpoly_q_clear(x, F->mctx);
        _acb_vec_clear(vec, 3);
        flint_free(sd);
        _fmpq_vec_clear(cc, nc);
        _fmpq_vec_clear(dd, nc);
    }

    flint_free(cand);
    return changed;
}

/*
    The gamma generators at rational arguments with definition order
    < limit are grouped by levels N <= the option
    GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT: for
    each generator (latest first) as a candidate for elimination, the
    level is the lcm of its denominator and those of the other
    generators which keep it within the limit (a tower may contain
    values at many unrelated levels).
*/
int
_gr_tower_gamma_round(gr_tower_flat_t F, slong limit, slong depth)
{
    slong level_limit = GR_TOWER_OPTION(F->T, GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT);
    gr_tower_struct * T = F->T;
    slong * ad, * ap, * aq, * bd, * bp, * bq, * tried;
    slong na = 0, ntried = 0, t, j;
    int changed = 0;

    (void) depth;

    ad = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));
    ap = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));
    aq = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));
    bd = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));
    bp = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));
    bq = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));
    tried = flint_malloc(sizeof(slong) * FLINT_MAX(limit, 1));

    for (j = limit - 1; j >= 0; j--)
    {
        slong p, q;
        if (_gamma_gen_rational(&p, &q, T->gens + j, level_limit))
        {
            ad[na] = j;
            ap[na] = p;
            aq[na] = q;
            na++;
        }
    }

    for (t = 0; t < na && !changed; t++)
    {
        slong N = aq[t], nb = 0, i;
        int seen = 0;

        for (j = 0; j < na; j++)
        {
            slong M = (N / n_gcd(N, aq[j])) * aq[j];
            if (j != t && M <= level_limit)
                N = M;
        }

        for (i = 0; i < ntried; i++)
            if (tried[i] == N)
                seen = 1;
        if (seen)
            continue;
        tried[ntried++] = N;

        for (j = 0; j < na; j++)
        {
            if (N % aq[j] == 0)
            {
                bd[nb] = ad[j];
                bp[nb] = ap[j];
                bq[nb] = aq[j];
                nb++;
            }
        }

        changed = _gamma_round_level(F, N, bd, bp, bq, nb);
    }

    flint_free(ad); flint_free(ap); flint_free(aq);
    flint_free(bd); flint_free(bp); flint_free(bq);
    flint_free(tried);

    if (!changed)
        changed = _line_round(F, limit, depth, GR_TOWER_GAMMA, 0);

    return changed;
}


/* helpers shared with hurwitz_relations.c */
int
_gr_tower_special_prepare_constants(gr_tower_t T, ulong m, int need_pi, slong before)
{
    return _gamma_prepare_constants(T, m, need_pi, before);
}

int
_gr_tower_special_flat_zeta(fmpz_mpoly_q_t res, ulong m, gr_tower_flat_t F)
{
    return _flat_zeta(res, m, F);
}

/*
    res = E(zeta_m) as a flat element (the prime power roots of unity must
    be present; returns 0 otherwise): zeta_m^i is a monomial in their
    generators, with exponents reduced modulo the orders, so the sum is
    formed term by term and reduced once.
*/
int
_gr_tower_special_flat_cyclotomic_eval(fmpz_mpoly_q_t res, const fmpq_poly_t E, ulong m, gr_tower_flat_t F)
{
    gr_tower_struct * T = F->T;
    slong nvars = F->mctx->minfo->nvars, nf = 0, i, j;
    slong var[FLINT_BITS];
    ulong ex[FLINT_BITS], ord[FLINT_BITS];
    int neg = 0, ok = 1;
    n_factor_t fac;
    ulong * exp;
    fmpz_t c;

    n_factor_init(&fac);
    if (m > 1)
        n_factor(&fac, m, 1);

    for (i = 0; i < fac.num && ok; i++)
    {
        ulong pe = n_pow(fac.p[i], fac.exp[i]);
        ulong a = n_invmod((m / pe) % pe, pe), pw, o;
        slong d;

        if (pe == 2)
        {
            neg = (a % 2 == 1);
            continue;
        }
        if (!_root_of_unity_gen(&d, &pw, T, pe))
        {
            ok = 0;
            break;
        }
        o = (T->gens[d].def_kind == GR_TOWER_ROOT_OF_UNITY) ? (ulong) T->gens[d].def_param : 4;
        for (j = 0; j < nf; j++)
            if (var[j] == GR_TOWER_FLAT_VAR_D(F, d))
                break;
        if (j == nf)
        {
            var[nf] = GR_TOWER_FLAT_VAR_D(F, d);
            ex[nf] = 0;
            ord[nf] = o;
            nf++;
        }
        ex[j] = (ex[j] + n_mulmod2(pw % o, a % o, o)) % o;
    }

    if (!ok)
        return 0;

    exp = flint_calloc(nvars, sizeof(ulong));
    fmpz_init(c);
    fmpz_mpoly_zero(fmpz_mpoly_q_numref(res), F->mctx);
    for (i = 0; i < fmpq_poly_length(E); i++)
    {
        if (fmpz_is_zero(E->coeffs + i))
            continue;
        for (j = 0; j < nf; j++)
            exp[var[j]] = n_mulmod2(ex[j], (ulong) i % ord[j], ord[j]);
        if (neg && (i % 2 == 1))
            fmpz_neg(c, E->coeffs + i);
        else
            fmpz_set(c, E->coeffs + i);
        fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_q_numref(res), c, exp, F->mctx);
    }
    fmpz_mpoly_sort_terms(fmpz_mpoly_q_numref(res), F->mctx);
    fmpz_mpoly_combine_like_terms(fmpz_mpoly_q_numref(res), F->mctx);
    fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(res), E->den, F->mctx);
    fmpz_mpoly_q_canonicalise(res, F->mctx);
    flint_free(exp);
    fmpz_clear(c);

    return gr_tower_flat_reduce(res, F) == GR_SUCCESS;
}

/* polygamma generators on lines, for each order present */
int
_gr_tower_hurwitz_line_round(gr_tower_flat_t F, slong limit, slong depth)
{
    gr_tower_struct * T = F->T;
    slong j, i;
    int changed = 0;

    for (j = limit - 1; j >= 0 && !changed; j--)
    {
        const gr_tower_gen_struct * g = T->gens + j;
        int seen = 0;
        if (g->kind != GR_TOWER_POLYGAMMA || g->def_param < 0 || g->def_param + 1 > 64)
            continue;
        /* (each order once: the latest generator of that order) */
        for (i = j + 1; i < limit && !seen; i++)
            if (T->gens[i].kind == GR_TOWER_POLYGAMMA && T->gens[i].def_param == g->def_param)
                seen = 1;
        if (!seen)
            changed = _line_round(F, limit, depth, GR_TOWER_POLYGAMMA, g->def_param);
    }

    return changed;
}

int
_gr_tower_special_root_gen(slong * d, slong * mult, gr_tower_flat_t F, const fmpz_mpoly_q_t c, ulong delta, slong before)
{
    return _root_gen(d, mult, F, c, delta, before);
}

/* (for the dilogarithm relations: a generator may also be moved or
   adjoined right before position `before`, the generator to eliminate;
   the caller restarts its round when 1 is returned, so the shifted
   positions are recomputed) */
int
_gr_tower_special_trans_gen(slong * d, gr_tower_flat_t F, int kind, const fmpz_mpoly_q_t u, slong before)
{
    return _gl_trans_gen(d, F, kind, u, before + 1);
}

POP_OPTIONS
