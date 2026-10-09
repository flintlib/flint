/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Lazy fields: radicals (square roots and n-th roots of rational numbers). */

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "gr_tower/lazy_impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* trial division by the first GR_TOWER_OPT_SMOOTH_LIMIT primes (the
   option is a number of primes); 1 when complete, otherwise the last
   factor is the unfactored cofactor */
static int
_gr_tower_lazy_factor_trial(fmpz_factor_t fac, const fmpz_t m, gr_ctx_t ctx)
{
    return _gr_tower_fmpz_factor_trial(fac, NULL, m, LAZY(ctx)->options[GR_TOWER_OPT_SMOOTH_LIMIT]);
}

/*
    The square root of the negative rational number c = -A B^2 / D^2 with
    A squarefree (as far as trial division and perfect square detection
    tell): B sqrt(-A) / D with sqrt(-A) a single generator, a root of
    x^2 + A (hash-consed with the algebraic numbers), except for A = 1
    (i) and A = 3 (2 zeta_3 + 1), which are roots of unity. Returns
    GR_UNABLE to leave the split form i sqrt(A) to the caller.
*/
static int
_gr_tower_lazy_sqrt_neg_fmpq(gr_tower_lazy_elem_t res, const fmpq_t c, gr_ctx_t ctx)
{
    fmpz_t m, A, B, t;
    fmpz_factor_t fac;
    fmpq_t s;
    slong i;
    int complete, status = GR_SUCCESS;

    fmpz_init(m);
    fmpz_init(A);
    fmpz_init(B);
    fmpz_init(t);
    fmpq_init(s);

    /* sqrt(-a/b) = sqrt(-a b) / b */
    fmpz_mul(m, fmpq_numref(c), fmpq_denref(c));
    fmpz_neg(m, m);

    fmpz_factor_init(fac);
    complete = _gr_tower_lazy_factor_trial(fac, m, ctx);

    fmpz_one(A);
    fmpz_one(B);
    for (i = 0; i < fac->num; i++)
    {
        fmpz_set(t, fac->p + i);
        if (!complete && i == fac->num - 1 && fmpz_is_square(t))
        {
            fmpz_sqrt(t, t);
            fmpz_pow_ui(t, t, fac->exp[i]);
            fmpz_mul(B, B, t);
            continue;
        }
        fmpz_pow_ui(t, fac->p + i, fac->exp[i] / 2);
        fmpz_mul(B, B, t);
        if (fac->exp[i] % 2)
            fmpz_mul(A, A, fac->p + i);
    }
    fmpz_factor_clear(fac);

    /* B / b */
    fmpz_set(fmpq_numref(s), B);
    fmpz_set(fmpq_denref(s), fmpq_denref(c));
    fmpq_canonicalise(s);

    if (fmpz_is_one(A))
    {
        status = _gr_tower_lazy_root_of_unity(res, 1, 4, ctx);
    }
    else if (fmpz_equal_ui(A, 3))
    {
        /* sqrt(-3) = 2 zeta_3 + 1 */
        status = _gr_tower_lazy_root_of_unity(res, 1, 3, ctx);
        status |= _gr_tower_lazy_add(res, res, res, ctx);
        {
            gr_tower_lazy_elem_struct one;
            _gr_tower_lazy_init(&one, ctx);
            status |= _gr_tower_lazy_one(&one, ctx);
            status |= _gr_tower_lazy_add(res, res, &one, ctx);
            _gr_tower_lazy_clear(&one, ctx);
        }
    }
    else
    {
        gr_tower_flat_struct * F;
        slong gid;
        qqbar_t z;

        /* (the root of x^2 + A with the enclosure i sqrt(A), whose real
           part is exactly zero) */
        qqbar_init(z);
        fmpz_neg(t, A);
        {
            fmpz_poly_t p;
            acb_t w;
            fmpz_poly_init(p);
            acb_init(w);
            fmpz_poly_set_coeff_fmpz(p, 0, A);
            fmpz_poly_set_coeff_ui(p, 2, 1);
            acb_set_fmpz(w, t);
            acb_sqrt(w, w, QQBAR_DEFAULT_PREC);
            status = qqbar_set_fmpz_poly_root(z, p, w, 2 * QQBAR_DEFAULT_PREC) ? GR_SUCCESS : GR_UNABLE;
            fmpz_poly_clear(p);
            acb_clear(w);
        }
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_qqbar_tower(&F, &gid, z, GR_TOWER_ROOT, 2, t, ctx);
        if (status == GR_SUCCESS)
            _gr_tower_lazy_set_gen_d(res, F, gr_tower_gid_order(F->T, gid), ctx);
        qqbar_clear(z);
    }

    if (status == GR_SUCCESS && !fmpq_is_one(s))
        status = gr_mul_fmpq(res, res, s, ctx);

    fmpz_clear(m);
    fmpz_clear(A);
    fmpz_clear(B);
    fmpz_clear(t);
    fmpq_clear(s);
    return status;
}

int
_gr_tower_lazy_root_fmpq(gr_tower_lazy_elem_t res, const fmpq_t c, ulong n, gr_ctx_t ctx)
{
    fmpz_t m, t;
    fmpz_factor_t fac;
    gr_tower_lazy_elem_struct u;
    slong i;
    int complete;
    int status = GR_SUCCESS;

    if (fmpq_is_zero(c))
        return _gr_tower_lazy_zero(res, ctx);

    if (fmpq_sgn(c) < 0 && n == 2 && !LAZY(ctx)->options[GR_TOWER_OPT_SPLIT_IMAGINARY])
    {
        status = _gr_tower_lazy_sqrt_neg_fmpq(res, c, ctx);
        if (status != GR_UNABLE)
            return status;
    }

    if (fmpq_sgn(c) < 0)
    {
        fmpq_t d;
        fmpq_init(d);
        fmpq_neg(d, c);
        status = _gr_tower_lazy_root_fmpq(res, d, n, ctx);
        if (status == GR_SUCCESS)
        {
            _gr_tower_lazy_init(&u, ctx);
            status = _gr_tower_lazy_root_of_unity(&u, 1, 2 * n, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_mul(res, res, &u, ctx);
            _gr_tower_lazy_clear(&u, ctx);
        }
        fmpq_clear(d);
        return status;
    }

    fmpz_init(m);
    fmpz_init(t);
    fmpz_pow_ui(m, fmpq_denref(c), n - 1);
    fmpz_mul(m, m, fmpq_numref(c));

    fmpz_factor_init(fac);
    complete = _gr_tower_lazy_factor_trial(fac, m, ctx);

    if (n == 2 && complete && LAZY(ctx)->options[GR_TOWER_OPT_COMPOSITE_RADICALS])
    {
        /* sqrt(m) / den with sqrt(m) = B sqrt(A), A squarefree, sqrt(A)
           one generator */
        fmpz_t A, B;
        fmpq_t r;
        fmpz_init(A);
        fmpz_init(B);
        fmpq_init(r);
        fmpz_one(A);
        fmpz_one(B);
        for (i = 0; i < fac->num; i++)
        {
            fmpz_pow_ui(t, fac->p + i, fac->exp[i] / 2);
            fmpz_mul(B, B, t);
            if (fac->exp[i] % 2 == 1)
                fmpz_mul(A, A, fac->p + i);
        }
        fmpz_set(fmpq_numref(r), B);
        fmpz_set(fmpq_denref(r), fmpq_denref(c));
        fmpq_canonicalise(r);
        if (fmpz_is_one(A))
            status = _gr_tower_lazy_set_fmpq(res, r, ctx);
        else
        {
            status = _gr_tower_lazy_prime_power_root(res, A, 2, 1, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_mul_fmpq(res, res, r, ctx);
        }
        fmpz_clear(A);
        fmpz_clear(B);
        fmpq_clear(r);
        fmpz_factor_clear(fac);
        fmpz_clear(m);
        fmpz_clear(t);
        return status;
    }

    {
        fmpq_t inv;
        fmpq_init(inv);
        fmpq_one(inv);
        fmpz_set(fmpq_denref(inv), fmpq_denref(c));
        status = _gr_tower_lazy_set_fmpq(res, inv, ctx);
        fmpq_clear(inv);
    }

    _gr_tower_lazy_init(&u, ctx);

    for (i = 0; i < fac->num && status == GR_SUCCESS; i++)
    {
        fmpz_t p;
        ulong e = fac->exp[i], q, r, k;

        fmpz_init_set(p, fac->p + i);

        if (!complete && i == fac->num - 1)
        {
            while ((k = fmpz_is_perfect_power(t, p)) > 1)
            {
                fmpz_set(p, t);
                e *= k;
            }
        }

        q = e / n;
        r = e % n;

        if (q > 0)
        {
            fmpz_pow_ui(t, p, q);
            status |= _gr_tower_lazy_set_fmpz(&u, t, ctx);
            status |= _gr_tower_lazy_mul(res, res, &u, ctx);
        }

        if (r > 0 && status == GR_SUCCESS)
        {
            status = _gr_tower_lazy_const_root(&u, p, n, ctx);
            if (status == GR_SUCCESS && r > 1)
                status = gr_pow_ui(&u, &u, r, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_mul(res, res, &u, ctx);
        }

        fmpz_clear(p);
    }

    _gr_tower_lazy_clear(&u, ctx);
    fmpz_factor_clear(fac);
    fmpz_clear(m);
    fmpz_clear(t);
    return status;
}

/*
    The principal p-th root of xt (the top field of T) when it lies in the
    field: the roots of X^p - xt by factoring (all the steps of T proven),
    the principal one identified by the isolating enclosure of the
    principal root (zx: the enclosure of xt, as for the adjunction).
    Returns 1 if found (res set).
*/
static int
_principal_root_in_field(gr_ptr res, gr_srcptr xt, ulong p, const acb_t zx, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_poly_t q;
    gr_vec_t roots;
    fmpz_vec_t mult;
    acb_t zp, w;
    slong i, prec, found = -1;
    int ok = 0;

    gr_poly_init(q, top);
    gr_vec_init(roots, 0, top);
    fmpz_vec_init(mult, 0);
    acb_init(zp);
    acb_init(w);

    if (gr_poly_set_coeff_si(q, p, 1, top) == GR_SUCCESS &&
        gr_neg(gr_poly_coeff_ptr(q, 0, top), xt, top) == GR_SUCCESS &&
        _gr_tower_principal_root_enclosure(zp, xt, p, zx, T) == GR_SUCCESS &&
        gr_tower_poly_roots(roots, mult, q, T->length, T) == GR_SUCCESS)
    {
        /* (zp contains exactly one root of X^p - xt) */
        for (prec = GR_TOWER_DEFAULT_PREC; prec <= GR_TOWER_OPTION(T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT) && found < 0 && roots->length > 0; prec *= 2)
            for (i = 0; i < roots->length && found < 0; i++)
                if (gr_tower_get_acb(w, gr_vec_entry_srcptr(roots, i, top), prec, T) == GR_SUCCESS && acb_contains(zp, w))
                    found = i;
        if (found >= 0)
            ok = (gr_set(res, gr_vec_entry_srcptr(roots, found, top), top) == GR_SUCCESS);
    }

    acb_clear(zp);
    acb_clear(w);
    fmpz_vec_clear(mult);
    gr_vec_clear(roots, top);
    gr_poly_clear(q, top);
    return ok;
}

/* the bound on p [K : F_0] for the check for powers: the option, 3/8 of
   it with transcendental generators (over which the norms of the
   factorization are multivariate), halved for each further one */
static slong
_power_check_limit(gr_tower_t T)
{
    slong lim = GR_TOWER_OPTION(T, GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT);
    if (T->num_trans > 0)
        lim = (T->num_trans > FLINT_BITS / 2) ? 0 : ((lim / 8) * 3) >> (T->num_trans - 1);
    return lim;
}

int
_gr_tower_lazy_root_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, ulong n, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * x = (gr_tower_lazy_elem_struct *) x_in;
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    gr_ctx_struct * top;
    gr_ptr xt, r;
    gr_poly_t q;
    acb_t z, zx;
    int status, real_cut = 0;
    ulong outer = 1;     /* (> 1: the result is the principal outer-th root of the one found) */
    int power_check;
    slong power_limit;

    if (n == 0)
        return GR_DOMAIN;

    if (n == 1)
        return _gr_tower_lazy_set(res, x, ctx);

    x = _gr_tower_lazy_flat_view(x);

    if (x->F != LAZY(ctx)->trivial &&
        fmpz_mpoly_q_is_fmpq(&x->elem.flat.data, x->elem.flat.mctx))
    {
        /* a rational number in an unrelated tower: move it to the trivial
           tower so that the root gets a tower of its own */
        gr_tower_lazy_elem_struct t;
        fmpq_t c;
        fmpq_init(c);
        (void) fmpz_mpoly_q_get_fmpq(c, &x->elem.flat.data, x->elem.flat.mctx);
        _gr_tower_lazy_init(&t, ctx);
        fmpq_set(&t.elem.q, c);
        status = _gr_tower_lazy_root_ui(res, &t, n, ctx);
        _gr_tower_lazy_clear(&t, ctx);
        fmpq_clear(c);
        return status;
    }

    if (x->F != LAZY(ctx)->trivial && x->level > 0 &&
        x->F->T->num_gens - x->level >= LAZY_ROOT_FORK_SUFFIX)
    {
        /* x in the prefix of a tower continuing with many generators
           (those of other computations, in a long session): the root
           goes into a copy of the prefix, rather than after them (the
           root would carry them in its prefix, and the relation
           searches with its later uses would visit them) */
        gr_tower_lazy_elem_struct t;
        _gr_tower_lazy_init(&t, ctx);
        _gr_tower_lazy_prefix_copy(&t, x, ctx);
        status = _gr_tower_lazy_root_ui(res, &t, n, ctx);
        _gr_tower_lazy_clear(&t, ctx);
        return status;
    }

    F = x->F;
    if (F == LAZY(ctx)->trivial)
    {
        fmpq_t c;

        fmpq_init(c);
        (void) fmpz_mpoly_q_get_fmpq(c, &x->elem.flat.data, x->elem.flat.mctx);
        status = _gr_tower_lazy_root_fmpq(res, c, n, ctx);
        fmpq_clear(c);
        return status;
    }

    T = F->T;

    /* c times a monomial in roots of unity and roots of integers: the
       principal root is root_n(c) times a root of unity (a root of a
       root of unity is a root of unity, in the registry:
       root_3(zeta_5^3) = zeta_15^(-2)) times roots of the integers of
       higher orders (sqrt(4 sqrt 2) = 2 root_4(2)) */
    {
        const fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(&x->elem.flat.data);
        if (num->length == 1 && fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(&x->elem.flat.data), x->elem.flat.mctx))
        {
            ulong * exp = flint_malloc(sizeof(ulong) * F->cap);
            fmpq_t r, c;
            slong d, num_rad = 0;
            int ok = 1;
            /* the radical factors: p_i^(e_i / (q_i n)) */
            fmpz * rad_p = _fmpz_vec_init(FLINT_MAX(T->num_gens, 1));
            ulong * rad_q = flint_malloc(sizeof(ulong) * FLINT_MAX(T->num_gens, 1));
            ulong * rad_e = flint_malloc(sizeof(ulong) * FLINT_MAX(T->num_gens, 1));

            fmpq_init(r);
            fmpq_init(c);
            fmpz_mpoly_get_term_exp_ui(exp, num, 0, x->elem.flat.mctx);
            for (d = 0; d < T->num_gens && ok; d++)
            {
                slong v = GR_TOWER_FLAT_VAR_D(F, d);
                const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
                fmpz_t p;
                ulong q;
                int kind;
                if (exp[v] == 0)
                    continue;
                fmpz_init(p);
                kind = _gr_tower_gen_const_root(g, p, &q);
                if (kind == 1)
                {
                    fmpq_t t;
                    fmpq_init(t);
                    fmpq_set_si(t, exp[v], g->def_param);
                    fmpq_add(r, r, t);
                    fmpq_clear(t);
                }
                else if (kind == 2 && g->status == GR_TOWER_STATUS_PROVEN)
                {
                    fmpz_set(rad_p + num_rad, p);
                    rad_q[num_rad] = q;
                    rad_e[num_rad] = exp[v];
                    num_rad++;
                }
                else
                    ok = 0;
                fmpz_clear(p);
            }
            if (ok && (!fmpq_is_zero(r) || num_rad > 0))
            {
                gr_tower_lazy_elem_struct u;
                fmpz_t tt;

                fmpz_mpoly_get_term_coeff_fmpz(fmpq_numref(c), num, 0, x->elem.flat.mctx);
                fmpz_mpoly_get_fmpz(fmpq_denref(c), fmpz_mpoly_q_denref(&x->elem.flat.data), x->elem.flat.mctx);
                fmpq_canonicalise(c);
                if (fmpq_sgn(c) < 0)
                {
                    /* -1 = exp(2 pi i / 2) */
                    fmpq_t h;
                    fmpq_neg(c, c);
                    fmpq_init(h);
                    fmpq_set_si(h, 1, 2);
                    fmpq_add(r, r, h);
                    fmpq_clear(h);
                }
                /* the principal argument: r in (-1/2, 1/2] (as a fraction of a turn) */
                fmpz_init(tt);
                fmpz_fdiv_q(tt, fmpq_numref(r), fmpq_denref(r));
                fmpz_submul(fmpq_numref(r), tt, fmpq_denref(r));   /* r in [0, 1) */
                {
                    fmpq_t h;
                    fmpq_init(h);
                    fmpq_set_si(h, 1, 2);
                    if (fmpq_cmp(r, h) > 0)
                        fmpq_sub_si(r, r, 1);
                    fmpq_clear(h);
                }
                fmpz_clear(tt);
                /* the root: exp(2 pi i (r / n)) root_n(c) */
                fmpz_mul_ui(fmpq_denref(r), fmpq_denref(r), n);
                fmpq_canonicalise(r);

                _gr_tower_lazy_init(&u, ctx);
                if (fmpz_fits_si(fmpq_numref(r)) && fmpz_abs_fits_ui(fmpq_denref(r)))
                    status = _gr_tower_lazy_root_of_unity(&u, fmpz_get_si(fmpq_numref(r)), fmpz_get_ui(fmpq_denref(r)), ctx);
                else
                    status = GR_UNABLE;
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_root_fmpq(res, c, n, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_mul(res, res, &u, ctx);
                /* the radicals: p^(e/(q n)) = root_{q n}(p)^e */
                for (d = 0; d < num_rad && status == GR_SUCCESS; d++)
                {
                    fmpq_t pq;
                    ulong qn = rad_q[d] * n;
                    fmpq_init(pq);
                    fmpz_set(fmpq_numref(pq), rad_p + d);
                    fmpz_one(fmpq_denref(pq));
                    status = _gr_tower_lazy_root_fmpq(&u, pq, qn, ctx);
                    if (status == GR_SUCCESS && rad_e[d] > 1)
                        status = gr_pow_ui(&u, &u, rad_e[d], ctx);
                    if (status == GR_SUCCESS)
                        status = _gr_tower_lazy_mul(res, res, &u, ctx);
                    fmpq_clear(pq);
                }
                _gr_tower_lazy_clear(&u, ctx);
                flint_free(exp);
                fmpq_clear(r);
                fmpq_clear(c);
                _fmpz_vec_clear(rad_p, FLINT_MAX(T->num_gens, 1));
                flint_free(rad_q);
                flint_free(rad_e);
                if (status == GR_SUCCESS)
                    return status;
                /* (otherwise: the general path) */
                status = GR_SUCCESS;
                x = _gr_tower_lazy_flat_view(x);
                F = x->F;
                T = F->T;
            }
            else
            {
                flint_free(exp);
                fmpq_clear(r);
                fmpq_clear(c);
                _fmpz_vec_clear(rad_p, FLINT_MAX(T->num_gens, 1));
                flint_free(rad_q);
                flint_free(rad_e);
            }
        }
    }

    /* a positive rational content c is taken out: root_n(c y) =
       root_n(c) root_n(y) for c > 0 (the principal branch), and root_n(c)
       is structured (sqrt(4 pi) = 2 sqrt(pi)) */
    {
        fmpz_t cn, cd;
        fmpz_init(cn);
        fmpz_init(cd);
        _fmpz_vec_content(cn, fmpz_mpoly_q_numref(&x->elem.flat.data)->coeffs, fmpz_mpoly_q_numref(&x->elem.flat.data)->length);
        _fmpz_vec_content(cd, fmpz_mpoly_q_denref(&x->elem.flat.data)->coeffs, fmpz_mpoly_q_denref(&x->elem.flat.data)->length);
        if (!fmpz_is_zero(cn) && !(fmpz_is_one(cn) && fmpz_is_one(cd)))
        {
            /* (a content in the denominator is possible only when the
               fraction is not canonical; the numerator's is the usual case) */
            fmpq_t c;
            gr_tower_lazy_elem_struct y, u;
            fmpq_init(c);
            fmpz_set(fmpq_numref(c), cn);
            fmpz_set(fmpq_denref(c), cd);
            fmpq_canonicalise(c);
            if (!fmpq_is_one(c))
            {
                _gr_tower_lazy_init(&y, ctx);
                _gr_tower_lazy_init(&u, ctx);
                status = _gr_tower_lazy_set(&y, x, ctx);
                if (status == GR_SUCCESS)
                {
                    fmpz_mpoly_q_div_fmpq(&y.elem.flat.data, &y.elem.flat.data, c, y.elem.flat.mctx);
                    status = _gr_tower_lazy_root_ui(&u, &y, n, ctx);
                }
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_root_fmpq(res, c, n, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_mul(res, res, &u, ctx);
                _gr_tower_lazy_clear(&y, ctx);
                _gr_tower_lazy_clear(&u, ctx);
                fmpq_clear(c);
                fmpz_clear(cn);
                fmpz_clear(cd);
                return status;
            }
            fmpq_clear(c);
        }
        fmpz_clear(cn);
        fmpz_clear(cd);
    }

    /* a square root of a flat element whose numerator and denominator
       are squares as polynomials (a square written out, (a + b)^2 =
       a^2 + 2 a b + b^2): the root up to sign, without searching the
       tower (the lattice search is costly when the tower has steps of
       degree other than 2) */
    if (n == 2)
    {
        fmpz_mpoly_q_t sq;
        fmpz_mpoly_t mn;
        int ok, imag = 0;
        x = _gr_tower_lazy_flat_view(x);
        fmpz_mpoly_q_init(sq, x->elem.flat.mctx);
        fmpz_mpoly_init(mn, x->elem.flat.mctx);
        ok = fmpz_mpoly_sqrt(fmpz_mpoly_q_denref(sq), fmpz_mpoly_q_denref(&x->elem.flat.data), x->elem.flat.mctx) &&
             !fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(&x->elem.flat.data), x->elem.flat.mctx);
        if (ok && !fmpz_mpoly_sqrt(fmpz_mpoly_q_numref(sq), fmpz_mpoly_q_numref(&x->elem.flat.data), x->elem.flat.mctx))
        {
            /* -(a square): i times its root */
            fmpz_mpoly_neg(mn, fmpz_mpoly_q_numref(&x->elem.flat.data), x->elem.flat.mctx);
            ok = fmpz_mpoly_sqrt(fmpz_mpoly_q_numref(sq), mn, x->elem.flat.mctx);
            imag = 1;
        }
        fmpz_mpoly_clear(mn, x->elem.flat.mctx);
        if (ok)
        {
            acb_t zx, zs;
            slong prec;
            int sign = 0;
            fmpz_mpoly_q_canonicalise(sq, x->elem.flat.mctx);
            acb_init(zx);
            acb_init(zs);
            /* the principal root: the sign agreeing with sqrt(x) */
            for (prec = GR_TOWER_DEFAULT_PREC; prec <= LAZY(ctx)->options[GR_TOWER_OPT_PREC_LIMIT] && sign == 0; prec *= 2)
            {
                if (gr_tower_lazy_get_acb(zx, x, prec, ctx) != GR_SUCCESS ||
                    gr_tower_flat_get_acb(zs, sq, prec, x->F) != GR_SUCCESS)
                    break;
                if (imag)
                    acb_mul_onei(zs, zs);
                acb_sqrt(zx, zx, prec);
                if (acb_overlaps(zx, zs))
                {
                    acb_neg(zs, zs);
                    if (!acb_overlaps(zx, zs))
                        sign = 1;
                }
                else
                {
                    acb_neg(zs, zs);
                    if (acb_overlaps(zx, zs))
                        sign = -1;
                    else
                        break;
                }
            }
            acb_clear(zx);
            acb_clear(zs);
            if (sign != 0)
            {
                gr_tower_flat_struct * G = x->F;
                int st = GR_SUCCESS;
                if (sign < 0)
                    fmpz_mpoly_q_neg(sq, sq, x->elem.flat.mctx);
                _gr_tower_lazy_install(res, G, sq, ctx);
                fmpz_mpoly_q_clear(sq, G->mctx);
                if (imag)
                {
                    gr_tower_lazy_elem_struct ii;
                    _gr_tower_lazy_init(&ii, ctx);
                    st = _gr_tower_lazy_i(&ii, ctx);
                    if (st == GR_SUCCESS)
                        st = _gr_tower_lazy_mul(res, res, &ii, ctx);
                    _gr_tower_lazy_clear(&ii, ctx);
                }
                return st;
            }
        }
        fmpz_mpoly_q_clear(sq, x->elem.flat.mctx);
    }

    /* the same root already adjoined to this tower? (or a root of a
       multiple order: root_n(y) = root_(k n)(y)^k for principal roots,
       sqrt(y) = root_4(y)^2) */
    {
        slong d;
        for (d = 0; d < T->num_gens; d++)
        {
            gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
            fmpz_mpoly_q_t diff;
            truth_t t;
            slong gid;

            if (g->def_kind != GR_TOWER_ROOT || g->def_param < (slong) n ||
                g->def_param % (slong) n != 0 || g->arg.mctx == NULL)
                continue;

            gid = g->gid;
            x = _gr_tower_lazy_flat_view(x);
            fmpz_mpoly_q_init(diff, F->mctx);
            gr_tower_flat_convert(diff, &g->arg.data, g->arg.mctx, F);
            fmpz_mpoly_q_sub(diff, diff, &x->elem.flat.data, F->mctx);
            t = gr_tower_flat_num_is_zero(diff, F);
            fmpz_mpoly_q_clear(diff, F->mctx);

            if (t == T_TRUE)
            {
                ulong k = g->def_param / n;
                _gr_tower_lazy_set_gen_d(res, F, gr_tower_gid_order(T, gid), ctx);
                return (k == 1) ? GR_SUCCESS : gr_pow_ui(res, res, k, ctx);
            }
        }
    }

    /* a negative real number whose enclosure does not have an exactly
       zero imaginary part (a real number expressed through complex
       generators): the principal root is decided from the exact
       realness, rather than from an enclosure straddling the cut (the
       test comes first: it may change the tower) */
    if (LAZY(ctx)->conj_depth == 0)
    {
        acb_t w;
        acb_init(w);
        if (_gr_tower_lazy_get_acb_impl(w, x, GR_TOWER_DEFAULT_PREC, ctx) == GR_SUCCESS &&
            arb_is_negative(acb_realref(w)) && arb_contains_zero(acb_imagref(w)) && !arb_is_zero(acb_imagref(w)))
            real_cut = (gr_tower_lazy_is_real(x, ctx) == T_TRUE);
        acb_clear(w);
    }

    /* (roots are adjoined to x's tower even when it is long, up to
       LAZY_ROOT_FORK_SUFFIX: the lattice search sees all of it, and
       finds roots among its later generators, sqrt(pi (5 - 2 sqrt 6)) =
       sqrt(3 pi) - sqrt(2 pi) say, which a copy of the prefix would
       have to rediscover) */
    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;

    /* the check for powers below (in fields of degree at most the limit
       over p, p the smallest prime factor of n: _power_check_limit) needs
       the steps proven: the proofs are attempted once per version of the
       moduli, before the nested form of x is taken (they may refine the
       tower) */
    {
        ulong p = 2;
        while (n % p != 0)
            p++;
        power_limit = _power_check_limit(T);
        power_check = (gr_tower_degree(T) <= power_limit / (slong) p);
    }
    if (power_check)
    {
        slong k;
        for (k = 1; k <= T->length; k++)
        {
            if (GR_TOWER_STEP(T, k - 1)->status != GR_TOWER_STATUS_PROVEN)
            {
                gr_tower_prove_modular(T, GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES));
                break;
            }
        }
    }

    top = gr_tower_field(T);

    /* nested representation of x in the top field */
    GR_TMP_INIT(xt, top);
    {
        slong k = _gr_tower_lazy_alg_level(x);
        gr_ctx_struct * Fk = gr_tower_field_at(T, k);
        gr_ptr u;
        GR_TMP_INIT(u, Fk);
        status = gr_tower_flat_get_nested_at(u, &x->elem.flat.data, k, x->F);
        if (status == GR_SUCCESS)
            status = gr_tower_promote(xt, u, k, T->length, T);
        GR_TMP_CLEAR(u, Fk);
    }
    if (status != GR_SUCCESS)
    {
        GR_TMP_CLEAR(xt, top);
        return status;
    }

    gr_poly_init(q, top);
    acb_init(z);
    acb_init(zx);
    r = gr_heap_init(top);

    status = gr_poly_set_coeff_si(q, n, 1, top);
    status |= gr_neg(gr_poly_coeff_ptr(q, 0, top), xt, top);
    _gr_poly_normalise(q, top);
    if (status == GR_SUCCESS)
        status = _gr_tower_get_acb_accurate(zx, xt, GR_TOWER_DEFAULT_PREC, T);

    if (status == GR_SUCCESS && real_cut && arb_contains_zero(acb_imagref(zx)))
        arb_zero(acb_imagref(zx));

    /* an enclosure of the principal root isolating it from the other
       roots (GR_UNABLE when x is not known to be real and its enclosure
       straddles the negative real axis) */
    if (status == GR_SUCCESS)
        status = _gr_tower_principal_root_enclosure(z, xt, n, zx, T);

    if (status == GR_SUCCESS)
    {
        int found = 0;
        int st;

        /* exact search through quadratic steps first, then the lattice
           search for towers of moderate degree */
        st = (n == 2) ? gr_tower_sqrt(r, xt, T) : GR_UNABLE;
        if (st == GR_SUCCESS)
            found = 1;
        else if (st == GR_UNABLE && gr_tower_degree(T) <= GR_TOWER_OPTION(T, GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT) &&
                 !gr_tower_poly_no_roots_modular(q, T, GR_TOWER_OPTION(T, GR_TOWER_OPT_NO_ROOTS_TRIES)))
        {
            found = (gr_tower_express_limit(r, q, z, GR_TOWER_MERGE_EXPRESS_PREC(T, gr_tower_degree(T)), T) == GR_SUCCESS);
        }

        if (found)
        {
            fmpz_mpoly_q_t t;
            gr_tower_flat_ensure(F);
            fmpz_mpoly_q_init(t, F->mctx);
            status = gr_tower_flat_set_nested_at(t, r, T->length, F);
            if (status == GR_SUCCESS)
            {
                _gr_tower_lazy_install(res, F, t, ctx);
                res->reduced_version = 0;
            }
            fmpz_mpoly_q_clear(t, F->mctx);
        }
        else
        {
            int known = -1;

            /*
                The radical a likely p-th power for a prime p | n (by the
                places of the proof of irreducibility of X^n - x, which
                fails then): the principal p-th root y found in the field
                (not left to the zero tests, as a generator with a
                reducible modulus), and root_n(x) = root_(n/p)(y) for
                principal roots (the argument of y is in (-pi/p, pi/p]).
            */
            if (power_check)
            {
                fmpz_mpoly_q_t arg;
                const fmpz_mpoly_ctx_struct * actx;
                ulong p = 0;

                gr_tower_flat_ensure(&T->flat);
                actx = T->flat.mctx;
                fmpz_mpoly_q_init(arg, actx);
                if (gr_tower_flat_set_nested_at(arg, xt, T->length, &T->flat) == GR_SUCCESS)
                    known = _gr_tower_binomial_modular_evidence(T, arg, actx, n, GR_TOWER_OPTION(T, GR_TOWER_OPT_MODULAR_TRIES), &p)
                            ? GR_TOWER_STATUS_PROVEN : GR_TOWER_STATUS_DYNAMIC;
                fmpz_mpoly_q_clear(arg, actx);

                if (p != 0 && (slong) p * gr_tower_degree(T) <= power_limit &&
                    _principal_root_in_field(r, xt, p, zx, T))
                {
                    fmpz_mpoly_q_t t;
                    gr_tower_flat_ensure(F);
                    fmpz_mpoly_q_init(t, F->mctx);
                    status = gr_tower_flat_set_nested_at(t, r, T->length, F);
                    if (status == GR_SUCCESS)
                    {
                        _gr_tower_lazy_install(res, F, t, ctx);
                        res->reduced_version = 0;
                        outer = n / p;
                    }
                    fmpz_mpoly_q_clear(t, F->mctx);
                    known = -2;
                }
            }

            if (known != -2)
            {
                status = _gr_tower_adjoin_root_ui_enclosure(T, xt, n, zx, known, NULL);
                if (status == GR_SUCCESS)
                {
                    _gr_tower_lazy_new_def(GR_TOWER_STEP(T, T->length - 1), T, ctx);
                    _gr_tower_lazy_set_gen(res, F, T->length, ctx);
                }
            }
        }
    }

    gr_heap_clear(r, top);
    gr_poly_clear(q, top);
    acb_clear(z);
    acb_clear(zx);
    GR_TMP_CLEAR(xt, top);

    if (status == GR_SUCCESS && outer > 1)
        status = _gr_tower_lazy_root_ui(res, res, outer, ctx);

    return status;
}

POP_OPTIONS
