/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Lazy fields: exponentials, logarithms and transcendental generators. */

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "gr_tower/lazy_impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* -------------------------------------------------------------------- */
/* transcendental functions                                              */
/* -------------------------------------------------------------------- */

/* Sets res to the element x + y i of the tower Q(i) of i (hash-consed). */
static int
_gr_tower_lazy_set_gaussian(gr_tower_lazy_elem_t res, const fmpq_t x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct t;
    int status;

    if (fmpq_is_zero(y))
        return _gr_tower_lazy_set_fmpq(res, x, ctx);

    _gr_tower_lazy_init(&t, ctx);
    status = _gr_tower_lazy_i(&t, ctx);
    if (status == GR_SUCCESS)
        status = gr_mul_fmpq(&t, &t, y, ctx);
    if (status == GR_SUCCESS)
        status = gr_add_fmpq(res, &t, x, ctx);
    _gr_tower_lazy_clear(&t, ctx);
    return status;
}

/* Looks up or creates the hash-consed generator for a constant argument
   c = x + y i (y = 0 except for logarithms of Gaussian rationals), with
   the parameter param for special functions. Returns the status of the
   adjunction (special function values may fail: poles). */
static int
_gr_tower_lazy_const_trans2_param(gr_tower_lazy_elem_t res, int kind, slong param, const fmpq_t x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_flat_struct * F;
    slong i, d;
    int status;

    for (i = 0; i < L->num_const_trans; i++)
    {
        if (L->const_trans[i].kind == kind && L->const_trans[i].param == param && L->const_trans[i].nxs == 0 &&
            fmpq_equal(&L->const_trans[i].x, x) && fmpq_equal(&L->const_trans[i].y, y))
        {
            _gr_tower_lazy_set_gen_d(res, L->const_trans[i].F, gr_tower_gid_order(L->const_trans[i].F->T, L->const_trans[i].gid), ctx);
            return GR_SUCCESS;
        }
    }

    if (GR_TOWER_KIND_IS_SPECIAL(kind))
    {
        fmpz_mpoly_q_t u;
        F = _gr_tower_lazy_new_tower(ctx);
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(u, F->mctx);
        fmpz_mpoly_q_set_fmpq(u, x, F->mctx);
        status = gr_tower_adjoin_special_flat(F->T, kind, param, u, F->mctx, NULL);
        fmpz_mpoly_q_clear(u, F->mctx);
        if (status != GR_SUCCESS)
        {
            /* (the new tower is empty: it is collected) */
            return status;
        }
    }
    else if (kind == GR_TOWER_PI || fmpq_is_zero(y))
    {
        F = _gr_tower_lazy_new_tower(ctx);
        if (kind == GR_TOWER_PI)
            status = gr_tower_adjoin_pi(F->T, NULL);
        else
        {
            fmpz_mpoly_q_t u;
            gr_tower_flat_ensure(F);
            fmpz_mpoly_q_init(u, F->mctx);
            fmpz_mpoly_q_set_fmpq(u, x, F->mctx);
            if (kind == GR_TOWER_EXP)
                status = gr_tower_adjoin_exp_flat(F->T, u, F->mctx, NULL);
            else if (kind == GR_TOWER_TAN)
                status = gr_tower_adjoin_tan_flat(F->T, u, F->mctx, NULL);
            else if (kind == GR_TOWER_ATAN)
                status = gr_tower_adjoin_atan_flat(F->T, u, F->mctx, NULL);
            else
                status = gr_tower_adjoin_log_flat(F->T, u, F->mctx, NULL);
            fmpz_mpoly_q_clear(u, F->mctx);
        }
    }
    else
    {
        /* a copy of the tower of i, extended by the generator */
        gr_tower_lazy_elem_struct g;
        gr_tower_flat_struct * Fi;

        _gr_tower_lazy_init(&g, ctx);
        status = _gr_tower_lazy_set_gaussian(&g, x, y, ctx);
        if (status != GR_SUCCESS)
        {
            _gr_tower_lazy_clear(&g, ctx);
            return status;
        }
        Fi = g.F;
        F = _gr_tower_lazy_new_tower(ctx);
        gr_tower_set(F->T, Fi->T);
        gr_tower_flat_ensure(F);
        {
            fmpz_mpoly_q_t u;
            fmpz_mpoly_q_init(u, F->mctx);
            _gr_tower_lazy_update(&g);
            _gr_tower_flat_transport(u, &g.elem.flat.data, Fi, F);
            if (kind == GR_TOWER_EXP)
                status = gr_tower_adjoin_exp_flat(F->T, u, F->mctx, NULL);
            else
                status = gr_tower_adjoin_log_flat(F->T, u, F->mctx, NULL);
            fmpz_mpoly_q_clear(u, F->mctx);
        }
        _gr_tower_lazy_clear(&g, ctx);
    }
    if (status != GR_SUCCESS)
        return status;   /* (the new tower is empty: it is collected) */

    d = F->T->num_gens - 1;
    _gr_tower_lazy_new_def(GR_TOWER_GEN(F->T, d), F->T, ctx);

    if (L->num_const_trans == L->alloc_const_trans)
    {
        L->alloc_const_trans = FLINT_MAX(4, 2 * L->alloc_const_trans);
        L->const_trans = flint_realloc(L->const_trans, L->alloc_const_trans * sizeof(gr_tower_lazy_const_trans_entry_struct));
    }
    L->const_trans[L->num_const_trans].kind = kind;
    L->const_trans[L->num_const_trans].param = param;
    L->const_trans[L->num_const_trans].nxs = 0;
    L->const_trans[L->num_const_trans].xs = NULL;
    fmpq_init(&L->const_trans[L->num_const_trans].x);
    fmpq_init(&L->const_trans[L->num_const_trans].y);
    fmpq_set(&L->const_trans[L->num_const_trans].x, x);
    fmpq_set(&L->const_trans[L->num_const_trans].y, y);
    L->const_trans[L->num_const_trans].F = F;
    F->gc |= GR_TOWER_GC_PINNED;
    L->const_trans[L->num_const_trans].gid = GR_TOWER_GEN(F->T, d)->gid;
    L->const_trans[L->num_const_trans].def_id = GR_TOWER_GEN(F->T, d)->def_id;
    L->num_const_trans++;

    _gr_tower_lazy_set_gen_d(res, F, d, ctx);
    return GR_SUCCESS;
}

static int
_gr_tower_lazy_const_trans2(gr_tower_lazy_elem_t res, int kind, const fmpq_t x, const fmpq_t y, gr_ctx_t ctx)
{
    return _gr_tower_lazy_const_trans2_param(res, kind, 0, x, y, ctx);
}

static int
_gr_tower_lazy_const_trans(gr_tower_lazy_elem_t res, int kind, const fmpq_t c, gr_ctx_t ctx)
{
    fmpq_t y;
    int status;
    fmpq_init(y);
    status = _gr_tower_lazy_const_trans2(res, kind, c, y, ctx);
    fmpq_clear(y);
    return status;
}

/* -------------------------------------------------------------------- */
/* logarithms of Gaussian rationals                                      */
/* -------------------------------------------------------------------- */

/*
    The logarithms of rational numbers span a vector space with the
    canonical basis {log p : p prime}, and those of Gaussian rationals
    add one generator log(a + b i) per prime p = a^2 + b^2 (a > b > 0),
    plus pi i. A logarithm log(A + B i) of a Gaussian integer with a
    smooth norm is decomposed over this basis (trial division by the
    first GR_TOWER_OPT_SMOOTH_LIMIT primes; a rough cofactor is kept as
    an opaque generator), the
    multiple of pi i / 2 coming from units and branches being determined
    numerically (it is an exact integer). This makes relations such as
    Machin-type formulas hold in the representation itself.
*/

/* a > b > 0 with a^2 + b^2 = p (p = 1 mod 4) */
static void
_gaussian_prime(ulong * a, ulong * b, ulong p)
{
    /* Cornacchia: from a square root x of -1 mod p, the first remainder
       below sqrt(p) in the Euclidean algorithm on (p, x) is a */
    ulong x = n_sqrtmod(p - 1, p), r0 = p, r1 = x, r;

    while (r1 * r1 > p)
    {
        r = r0 % r1;
        r0 = r1;
        r1 = r;
    }

    *a = r1;
    *b = n_sqrt(p - r1 * r1);
    if (*a < *b)
        FLINT_SWAP(ulong, *a, *b);
}

/* Divides A + B i by a + b i if the quotient is a Gaussian integer;
   returns 1 on success. */
static int
_gaussian_divexact(fmpz_t A, fmpz_t B, const fmpz_t a, const fmpz_t b)
{
    fmpzi_t x, y, q, r;
    int ok;

    fmpzi_init(x); fmpzi_init(y); fmpzi_init(q); fmpzi_init(r);
    fmpz_set(fmpzi_realref(x), A); fmpz_set(fmpzi_imagref(x), B);
    fmpz_set(fmpzi_realref(y), a); fmpz_set(fmpzi_imagref(y), b);
    fmpzi_divrem(q, r, x, y);
    ok = fmpzi_is_zero(r);
    if (ok)
    {
        fmpz_swap(A, fmpzi_realref(q));
        fmpz_swap(B, fmpzi_imagref(q));
    }
    fmpzi_clear(x); fmpzi_clear(y); fmpzi_clear(q); fmpzi_clear(r);
    return ok;
}

/* res += e * log(x + y i), through the hash-consed generator */
static int
_add_const_log(gr_tower_lazy_elem_t res, slong e, const fmpq_t x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct t;
    int status;

    if (e == 0)
        return GR_SUCCESS;

    _gr_tower_lazy_init(&t, ctx);
    status = _gr_tower_lazy_const_trans2(&t, GR_TOWER_LOG, x, y, ctx);
    if (status == GR_SUCCESS)
        status = gr_mul_si(&t, &t, e, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_add(res, res, &t, ctx);
    _gr_tower_lazy_clear(&t, ctx);
    return status;
}

/* Sets res to log(u + v i) for rationals u, v (not both zero). */
static int
_gr_tower_lazy_log_gaussian(gr_tower_lazy_elem_t res, const fmpq_t u, const fmpq_t v, gr_ctx_t ctx)
{
    fmpz_t A, B, D, N, a, b;
    fmpq_t qx, qy;
    gr_tower_lazy_elem_struct S;
    ulong p, extra[2 * FLINT_MAX_FACTORS_IN_LIMB];
    slong ip, nextra = 0, iextra = 0;
    n_primes_t iter;
    int status = GR_SUCCESS;

    fmpz_init(A); fmpz_init(B); fmpz_init(D); fmpz_init(N); fmpz_init(a); fmpz_init(b);
    fmpq_init(qx); fmpq_init(qy);
    _gr_tower_lazy_init(&S, ctx);

    /* u + v i = (A + B i) / D */
    fmpz_lcm(D, fmpq_denref(u), fmpq_denref(v));
    fmpz_divexact(A, D, fmpq_denref(u)); fmpz_mul(A, A, fmpq_numref(u));
    fmpz_divexact(B, D, fmpq_denref(v)); fmpz_mul(B, B, fmpq_numref(v));

    /* factor the norm of A + B i, and D, by trial division */
    fmpz_mul(N, A, A);
    fmpz_addmul(N, B, B);

    /* the primes: the first GR_TOWER_OPT_SMOOTH_LIMIT ones by trial
       division, then the prime factors of the cofactors of N and D when
       these fit in a word (n_factor); larger cofactors stay opaque */
    n_primes_init(iter);
    for (ip = 0; status == GR_SUCCESS; ip++)
    {
        slong e_D = 0, e_pi = 0, e_pibar = 0, e_p = 0;

        if (fmpz_is_one(N) && fmpz_is_one(D))
            break;

        if (ip < LAZY(ctx)->options[GR_TOWER_OPT_SMOOTH_LIMIT])
            p = n_primes_next(iter);
        else
        {
            if (ip == LAZY(ctx)->options[GR_TOWER_OPT_SMOOTH_LIMIT])
            {
                n_factor_t fN, fD;
                slong i2;
                n_factor_init(&fN);
                n_factor_init(&fD);
                /* (for a real argument, N = A^2 with the cofactor of A) */
                if (fmpz_is_zero(B))
                {
                    fmpz_t r;
                    fmpz_init(r);
                    fmpz_sqrt(r, N);
                    if (fmpz_abs_fits_ui(r) && !fmpz_is_one(r))
                        n_factor(&fN, fmpz_get_ui(r), 1);
                    fmpz_clear(r);
                }
                else if (fmpz_abs_fits_ui(N) && !fmpz_is_one(N))
                    n_factor(&fN, fmpz_get_ui(N), 1);
                if (fmpz_abs_fits_ui(D) && !fmpz_is_one(D))
                    n_factor(&fD, fmpz_get_ui(D), 1);
                nextra = 0;
                for (i2 = 0; i2 < fN.num; i2++)
                    extra[nextra++] = fN.p[i2];
                for (i2 = 0; i2 < fD.num; i2++)
                {
                    slong i3;
                    for (i3 = 0; i3 < nextra; i3++)
                        if (extra[i3] == fD.p[i2])
                            break;
                    if (i3 == nextra)
                        extra[nextra++] = fD.p[i2];
                }
                iextra = 0;
            }
            if (iextra >= nextra)
                break;
            p = extra[iextra++];
            if (p > (UWORD(1) << (FLINT_BITS - 2)))
                break;   /* (the slong arithmetic below) */
        }

        /* D: real prime logs */
        while (fmpz_divisible_si(D, p))
        {
            fmpz_divexact_si(D, D, p);
            e_D++;
        }

        if (fmpz_divisible_si(N, p))
        {
            if (p == 2)
            {
                /* (1 + i): log(1 + i) = log(2)/2 + pi i/4 (the pi i part
                   is absorbed in the numerical branch term below) */
                fmpz_one(a); fmpz_one(b);
                while (_gaussian_divexact(A, B, a, b))
                {
                    e_pi++;
                    fmpz_divexact_si(N, N, 2);
                }
                fmpq_set_si(qx, 2, 1); fmpq_zero(qy);
                /* e_pi * log(2)/2, as e_pi/2 * log 2: use the coefficient on the element */
                {
                    gr_tower_lazy_elem_struct t;
                    fmpq_t half;
                    _gr_tower_lazy_init(&t, ctx);
                    fmpq_init(half);
                    fmpq_set_si(half, e_pi, 2);
                    status |= _gr_tower_lazy_const_trans2(&t, GR_TOWER_LOG, qx, qy, ctx);
                    status |= gr_mul_fmpq(&t, &t, half, ctx);
                    status |= _gr_tower_lazy_add(&S, &S, &t, ctx);
                    fmpq_clear(half);
                    _gr_tower_lazy_clear(&t, ctx);
                }
                e_pi = 0;
            }
            else if (p % 4 == 3)
            {
                fmpz_set_ui(a, p); fmpz_zero(b);
                while (_gaussian_divexact(A, B, a, b))
                {
                    e_p++;
                    fmpz_divexact_si(N, N, p);
                    fmpz_divexact_si(N, N, p);
                }
            }
            else
            {
                ulong ga, gb;
                _gaussian_prime(&ga, &gb, p);
                fmpz_set_ui(a, ga); fmpz_set_ui(b, gb);
                while (_gaussian_divexact(A, B, a, b))
                {
                    e_pi++;
                    fmpz_divexact_si(N, N, p);
                }
                fmpz_neg(b, b);
                while (_gaussian_divexact(A, B, a, b))
                {
                    e_pibar++;
                    fmpz_divexact_si(N, N, p);
                }
                /* log(a - b i) = log p - log(a + b i) */
                fmpq_set_ui(qx, ga, 1); fmpq_set_ui(qy, gb, 1);
                status |= _add_const_log(&S, e_pi - e_pibar, qx, qy, ctx);
                e_p += e_pibar;
            }
        }

        fmpq_set_ui(qx, p, 1); fmpq_zero(qy);
        status |= _add_const_log(&S, e_p - e_D, qx, qy, ctx);
    }
    n_primes_clear(iter);

    /* rough cofactors: opaque generators */
    if (!fmpz_is_one(D))
    {
        fmpz_set(fmpq_numref(qx), D); fmpz_one(fmpq_denref(qx)); fmpq_zero(qy);
        status |= _add_const_log(&S, -1, qx, qy, ctx);
    }
    if (!fmpz_is_one(N))
    {
        /* A + B i is now a Gaussian integer of rough norm times a unit;
           take it as it is (the unit is accounted for numerically) */
        fmpz_set(fmpq_numref(qx), A); fmpz_one(fmpq_denref(qx));
        fmpz_set(fmpq_numref(qy), B); fmpz_one(fmpq_denref(qy));
        /* normalize to the first quadrant by multiplying by a unit */
        while (fmpz_sgn(fmpq_numref(qx)) < 0 || fmpz_sgn(fmpq_numref(qy)) < 0)
        {
            fmpz_t t;
            fmpz_init(t);
            /* multiply by i: (x + y i) i = -y + x i */
            fmpz_set(t, fmpq_numref(qx));
            fmpz_neg(fmpq_numref(qx), fmpq_numref(qy));
            fmpz_set(fmpq_numref(qy), t);
            fmpz_clear(t);
        }
        status |= _add_const_log(&S, 1, qx, qy, ctx);
    }

    /* the difference log(u + v i) - S is an integer multiple of pi i / 4
       (units, branches, and the factor 1 + i whose argument is pi/4) */
    if (status == GR_SUCCESS)
    {
        acb_t z, w;
        arb_t t;
        fmpz_t m;
        slong prec;
        int ok = 0;

        acb_init(z); acb_init(w); arb_init(t); fmpz_init(m);
        for (prec = 64; prec <= 4096 && !ok; prec *= 2)
        {
            arb_set_fmpq(acb_realref(z), u, prec);
            arb_set_fmpq(acb_imagref(z), v, prec);
            acb_log(z, z, prec);
            if (gr_tower_lazy_get_acb(w, &S, prec, ctx) != GR_SUCCESS)
                break;
            acb_sub(z, z, w, prec);
            /* z = m pi i / 4 */
            arb_const_pi(t, prec);
            arb_div(t, acb_imagref(z), t, prec);
            arb_mul_2exp_si(t, t, 2);
            if (arb_contains_zero(acb_realref(z)) && arb_get_unique_fmpz(m, t))
                ok = 1;
        }
        if (ok && !fmpz_is_zero(m))
        {
            gr_tower_lazy_elem_struct pii;
            fmpq_t c;
            _gr_tower_lazy_init(&pii, ctx);
            fmpq_init(c);
            status |= _gr_tower_lazy_pi(&pii, ctx);
            {
                gr_tower_lazy_elem_struct ii;
                _gr_tower_lazy_init(&ii, ctx);
                status |= _gr_tower_lazy_i(&ii, ctx);
                status |= _gr_tower_lazy_mul(&pii, &pii, &ii, ctx);
                _gr_tower_lazy_clear(&ii, ctx);
            }
            fmpz_set(fmpq_numref(c), m); fmpz_set_ui(fmpq_denref(c), 4); fmpq_canonicalise(c);
            status |= gr_mul_fmpq(&pii, &pii, c, ctx);
            status |= _gr_tower_lazy_add(&S, &S, &pii, ctx);
            fmpq_clear(c);
            _gr_tower_lazy_clear(&pii, ctx);
        }
        if (!ok)
            status = GR_UNABLE;
        acb_clear(z); acb_clear(w); arb_clear(t); fmpz_clear(m);
    }

    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_set(res, &S, ctx);

    _gr_tower_lazy_clear(&S, ctx);
    fmpz_clear(A); fmpz_clear(B); fmpz_clear(D); fmpz_clear(N); fmpz_clear(a); fmpz_clear(b);
    fmpq_clear(qx); fmpq_clear(qy);
    return status;
}

int
_gr_tower_lazy_pi(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    fmpq_t c;
    int status;
    fmpq_init(c);
    status = _gr_tower_lazy_const_trans(res, GR_TOWER_PI, c, ctx);
    fmpq_clear(c);
    return status;
}

/*
    Whether x = u + v i with rational u, v, for an element of a tower
    without transcendental generators: u = (x + conj(x)) / 2 and
    v = (x - conj(x)) / (2 i) must be rational numbers (which is decided
    exactly in the tower).
*/
static int
_gr_tower_lazy_get_gaussian(fmpq_t u, fmpq_t v, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct xc, a, b;
    int status;

    _gr_tower_lazy_init(&xc, ctx);
    _gr_tower_lazy_init(&a, ctx);
    _gr_tower_lazy_init(&b, ctx);

    status = _gr_tower_lazy_conj(&xc, x, ctx);
    if (status == GR_SUCCESS)
    {
        status = _gr_tower_lazy_add(&a, x, &xc, ctx);
        status |= gr_div_ui(&a, &a, 2, ctx);
        status |= _gr_tower_lazy_get_fmpq(u, &a, ctx);
    }
    if (status == GR_SUCCESS)
    {
        status = _gr_tower_lazy_sub(&b, x, &xc, ctx);
        status |= _gr_tower_lazy_i(&a, ctx);
        status |= gr_mul_ui(&a, &a, 2, ctx);
        status |= _gr_tower_lazy_div(&b, &b, &a, ctx);
        status |= _gr_tower_lazy_get_fmpq(v, &b, ctx);
    }

    _gr_tower_lazy_clear(&xc, ctx);
    _gr_tower_lazy_clear(&a, ctx);
    _gr_tower_lazy_clear(&b, ctx);
    return (status == GR_SUCCESS) ? GR_SUCCESS : GR_UNABLE;
}

/*
    Whether x = r pi i for a rational number r (set), decided exactly
    from the flat representation: x / pi must be a constant multiple of
    i, which is a generator or a power of a root of unity generator.
*/
static int
_gr_tower_lazy_pi_i_multiple(fmpq_t r, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    fmpz_mpoly_q_t y;
    fmpz_mpoly_t g;
    slong d, dpi = -1, di = -1, vpi;
    int ok = 0;

    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;

    for (d = 0; d < T->num_gens; d++)
    {
        const gr_tower_gen_struct * gen = GR_TOWER_GEN(T, d);
        if (gen->kind == GR_TOWER_PI)
            dpi = d;
        else if (gen->kind == GR_TOWER_ALGEBRAIC &&
                 ((gen->def_kind == GR_TOWER_ROOT_OF_UNITY && gen->def_param % 4 == 0) || _gr_tower_lazy_gen_is_i(T, d)))
            di = d;
    }

    if (dpi < 0 || di < 0)
        return 0;

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(y, F->mctx);
    fmpz_mpoly_init(g, F->mctx);
    fmpz_mpoly_q_set(y, &x->elem.flat.data, F->mctx);
    GR_MUST_SUCCEED(gr_tower_flat_reduce(y, F));

    vpi = GR_TOWER_FLAT_VAR_D(F, dpi);

    if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(y), F->mctx))
    {
        /* y = x / pi must not involve pi */
        fmpz_mpoly_gen(g, vpi, F->mctx);
        if (fmpz_mpoly_divides(fmpz_mpoly_q_numref(y), fmpz_mpoly_q_numref(y), g, F->mctx) &&
            fmpz_mpoly_degree_si(fmpz_mpoly_q_numref(y), vpi, F->mctx) <= 0)
        {
            /* y (-i) must be a constant */
            const gr_tower_gen_struct * gen = GR_TOWER_GEN(T, di);
            fmpz_mpoly_gen(g, GR_TOWER_FLAT_VAR_D(F, di), F->mctx);
            if (gen->def_kind == GR_TOWER_ROOT_OF_UNITY && gen->def_param != 4)
                fmpz_mpoly_pow_ui(g, g, 3 * (gen->def_param / 4), F->mctx);   /* -i = i^3 */
            else
                fmpz_mpoly_neg(g, g, F->mctx);
            fmpz_mpoly_mul(fmpz_mpoly_q_numref(y), fmpz_mpoly_q_numref(y), g, F->mctx);
            GR_MUST_SUCCEED(gr_tower_flat_reduce(y, F));
            if (fmpz_mpoly_q_get_fmpq(r, y, F->mctx))
            {
                fmpq_canonicalise(r);
                ok = 1;
            }
        }
    }

    fmpz_mpoly_clear(g, F->mctx);
    fmpz_mpoly_q_clear(y, F->mctx);
    return ok;
}

/* whether the n-th roots of unity are in the field of T (n <= 2, or a
   root of unity generator of order divisible by n, or i for n = 4) */
static int
_gr_tower_lazy_has_root_of_unity(gr_tower_t T, ulong n)
{
    slong d;
    if (n <= 2)
        return 1;
    for (d = 0; d < T->num_gens; d++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
        if (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param % n == 0)
            return 1;
        if (n == 4 && g->kind == GR_TOWER_ALGEBRAIC && _gr_tower_lazy_gen_is_i(T, d))
            return 1;
    }
    return 0;
}

/*
    For x (flat) affine in pi, pi not in the denominator: sets v to the
    variable of pi, initialises a to the coefficient of pi in the
    numerator, and returns 1; otherwise returns 0 (a not initialised).
*/
static int
_pi_coeff(slong * v, fmpz_mpoly_t a, const gr_tower_lazy_elem_struct * x)
{
    gr_tower_struct * T = x->F->T;
    const fmpz_mpoly_struct * num = fmpz_mpoly_q_numref(&x->elem.flat.data);
    const fmpz_mpoly_struct * den = fmpz_mpoly_q_denref(&x->elem.flat.data);
    slong d, dpi = -1;
    ulong one = 1;

    for (d = 0; d < T->num_gens; d++)
        if (T->gens[d].kind == GR_TOWER_PI)
            dpi = d;
    if (dpi < 0)
        return 0;

    *v = GR_TOWER_FLAT_VAR_D(x->F, dpi);
    if (fmpz_mpoly_degree_si(den, *v, x->elem.flat.mctx) > 0 || fmpz_mpoly_degree_si(num, *v, x->elem.flat.mctx) != 1)
        return 0;

    fmpz_mpoly_init(a, x->elem.flat.mctx);
    fmpz_mpoly_get_coeff_vars_ui(a, num, v, &one, 1, x->elem.flat.mctx);
    return 1;
}

/*
    x = r pi i + y with r rational nonzero, the part of x affine in pi
    whose coefficient is free of transcendental generators being r pi i
    (y is free of pi except for terms pi t with t involving other
    transcendental generators: x = pi (i log(3) - 9 i) gives r = -9 and
    y = pi i log(3)): sets r and y and returns 1, otherwise 0. Used so
    that exponentials are not adjoined with pi i in their arguments:
    exp(r pi i + y) is a root of unity times exp(y).
*/
static int
_gr_tower_lazy_pi_i_split(fmpq_t r, gr_tower_lazy_elem_t y, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    slong v;
    fmpz_mpoly_t a;
    const fmpz_mpoly_struct * den;
    int ok = 0;

    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;
    den = fmpz_mpoly_q_denref(&x->elem.flat.data);

    if (!_pi_coeff(&v, a, x))
        return 0;

    /* only the part of the coefficient of pi free of transcendental
       generators can be a rational multiple of i (in pi (i log(3) - 9 i),
       the part -9 pi i): keep those terms, provided the denominator is
       free of transcendental generators */
    if (T->num_trans > 1)
    {
        slong j, k, len = fmpz_mpoly_length(a, x->elem.flat.mctx);
        ulong * e = flint_malloc(sizeof(ulong) * x->elem.flat.mctx->minfo->nvars);
        fmpz_mpoly_t b;
        int den_ok = 1;

        for (j = 1; j <= T->num_trans; j++)
        {
            slong tv = GR_TOWER_FLAT_TVAR(F, j);
            if (tv != v && fmpz_mpoly_degree_si(den, tv, x->elem.flat.mctx) > 0)
                den_ok = 0;
        }

        fmpz_mpoly_init(b, x->elem.flat.mctx);
        for (k = 0; k < len && den_ok; k++)
        {
            int keep = 1;
            fmpz_mpoly_get_term_exp_ui(e, a, k, x->elem.flat.mctx);
            for (j = 1; j <= T->num_trans && keep; j++)
                if (e[GR_TOWER_FLAT_TVAR(F, j)] != 0)
                    keep = 0;
            if (keep)
            {
                fmpz_mpoly_t term;
                fmpz_mpoly_init(term, x->elem.flat.mctx);
                fmpz_mpoly_get_term(term, a, k, x->elem.flat.mctx);
                fmpz_mpoly_add(b, b, term, x->elem.flat.mctx);
                fmpz_mpoly_clear(term, x->elem.flat.mctx);
            }
        }
        if (den_ok)
            fmpz_mpoly_swap(a, b, x->elem.flat.mctx);
        else
            fmpz_mpoly_zero(a, x->elem.flat.mctx);
        fmpz_mpoly_clear(b, x->elem.flat.mctx);
        flint_free(e);

        if (fmpz_mpoly_is_zero(a, x->elem.flat.mctx))
        {
            fmpz_mpoly_clear(a, x->elem.flat.mctx);
            return 0;
        }
    }

    {
        gr_tower_lazy_elem_struct w;
        fmpz_mpoly_q_t t;

        /* w = a pi / den */
        _gr_tower_lazy_init(&w, ctx);
        fmpz_mpoly_q_init(t, x->elem.flat.mctx);
        fmpz_mpoly_set(fmpz_mpoly_q_numref(t), a, x->elem.flat.mctx);
        fmpz_mpoly_set(fmpz_mpoly_q_denref(t), den, x->elem.flat.mctx);
        fmpz_mpoly_q_canonicalise(t, x->elem.flat.mctx);
        {
            fmpz_mpoly_q_t pi;
            fmpz_mpoly_q_init(pi, x->elem.flat.mctx);
            fmpz_mpoly_q_gen(pi, v, x->elem.flat.mctx);
            fmpz_mpoly_q_mul(t, t, pi, x->elem.flat.mctx);
            fmpz_mpoly_q_clear(pi, x->elem.flat.mctx);
        }
        _gr_tower_lazy_fresh(&w, F, T->num_gens, ctx);
        fmpz_mpoly_q_swap(&w.elem.flat.data, t, w.elem.flat.mctx);
        _gr_tower_lazy_shrink(&w);

        if (_gr_tower_lazy_pi_i_multiple(r, &w, ctx) && !fmpq_is_zero(r))
        {
            /* y = x - w */
            fmpz_mpoly_q_sub(t, &x->elem.flat.data, &w.elem.flat.data, x->elem.flat.mctx);
            _gr_tower_lazy_fresh(y, F, T->num_gens, ctx);
            fmpz_mpoly_q_swap(&y->elem.flat.data, t, y->elem.flat.mctx);
            _gr_tower_lazy_shrink(y);
            ok = 1;
        }

        fmpz_mpoly_q_clear(t, x->elem.flat.mctx);
        _gr_tower_lazy_clear(&w, ctx);
    }

    fmpz_mpoly_clear(a, x->elem.flat.mctx);
    return ok;
}

/*
    If x = c exp(u) for a positive rational c and a generator exp(u)
    (possibly one which became algebraic) with Im(u) in (-pi, pi]
    (verified numerically: an enclosure of Im(u) which is exactly zero
    or lies strictly inside the interval), returns the definition order
    of the generator and sets c; otherwise returns -1.
*/
static slong
_gr_tower_lazy_exp_gen_multiple(fmpq_t c, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    const fmpz_mpoly_struct * num;
    ulong * exp;
    slong d, v, found = -1, i;

    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;
    num = fmpz_mpoly_q_numref(&x->elem.flat.data);

    if (fmpz_mpoly_length(num, F->mctx) != 1 || !fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(&x->elem.flat.data), F->mctx))
        return -1;

    exp = flint_malloc(sizeof(ulong) * F->cap);
    fmpz_mpoly_get_term_exp_ui(exp, num, 0, F->mctx);
    v = -1;
    for (i = 0; i < F->cap; i++)
    {
        if (exp[i] == 1 && v == -1)
            v = i;
        else if (exp[i] != 0)
        {
            v = -2;
            break;
        }
    }
    flint_free(exp);
    if (v < 0)
        return -1;

    for (d = 0; d < T->num_gens; d++)
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
        int is_exp = (g->kind == GR_TOWER_EXP) || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_EXP);
        if (GR_TOWER_FLAT_VAR_D(F, d) == v && is_exp && g->arg.mctx != NULL)
            found = d;
    }
    if (found < 0)
        return -1;

    fmpz_mpoly_get_term_coeff_fmpz(fmpq_numref(c), num, 0, F->mctx);
    fmpz_mpoly_get_fmpz(fmpq_denref(c), fmpz_mpoly_q_denref(&x->elem.flat.data), F->mctx);
    fmpq_canonicalise(c);
    if (fmpq_sgn(c) <= 0)
        return -1;

    /* Im(u) in (-pi, pi] */
    {
        const gr_tower_gen_struct * g = GR_TOWER_GEN(T, found);
        fmpz_mpoly_q_t u;
        acb_t z;
        arb_t pi;
        int ok = 0;

        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(u, F->mctx);
        gr_tower_flat_convert(u, &g->arg.data, g->arg.mctx, F);
        acb_init(z);
        arb_init(pi);
        if (gr_tower_flat_get_acb(z, u, GR_TOWER_DEFAULT_PREC, F) == GR_SUCCESS)
        {
            if (arb_is_zero(acb_imagref(z)))
                ok = 1;
            else
            {
                arb_const_pi(pi, GR_TOWER_DEFAULT_PREC);
                if (arb_lt(acb_imagref(z), pi))
                {
                    arb_neg(pi, pi);
                    if (arb_gt(acb_imagref(z), pi))
                        ok = 1;
                }
            }
        }
        acb_clear(z);
        arb_clear(pi);
        fmpz_mpoly_q_clear(u, F->mctx);
        if (!ok)
            return -1;
    }

    return found;
}

/*
    exp(x) for x a rational linear combination of logarithms of positive
    rational numbers (generators with constant arguments) and of pi i:
    exp(sum c_j log(p_j) + r pi i) = prod p_j^{c_j} exp(r pi i) is a product
    of structured radicals and a root of unity. Logarithms of algebraic
    numbers (log(2+i), say) are allowed too: p^c = exp(c log p) is the
    principal power, a power of the principal root. Returns 1 (with the
    status in *status) if x has that form, 0 otherwise.
*/
static int
_gr_tower_lazy_exp_of_logs(gr_tower_lazy_elem_t res, gr_tower_lazy_elem_t x, int * status, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    fmpz_mpoly_q_t y;
    fmpz_mpoly_t d, g;
    fmpz_t D;
    fmpq_t c, r, p;
    gr_tower_lazy_elem_struct t, u, w, acc;
    slong j, count = 0;
    int ok = 1, have_pi_i = 0, algebraic;

    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;
    if (T->num_trans == 0)
        return 0;

    gr_tower_flat_ensure(F);
    fmpz_mpoly_q_init(y, F->mctx);
    fmpz_mpoly_q_set(y, &x->elem.flat.data, F->mctx);
    GR_MUST_SUCCEED(gr_tower_flat_reduce(y, F));
    if (!fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(y), F->mctx))
    {
        fmpz_mpoly_q_clear(y, F->mctx);
        return 0;
    }

    fmpz_mpoly_init(d, F->mctx);
    fmpz_mpoly_init(g, F->mctx);
    fmpz_init(D);
    fmpq_init(c);
    fmpq_init(r);
    fmpq_init(p);
    _gr_tower_lazy_init(&t, ctx);
    _gr_tower_lazy_init(&u, ctx);
    _gr_tower_lazy_init(&w, ctx);
    _gr_tower_lazy_init(&acc, ctx);
    fmpz_mpoly_get_fmpz(D, fmpz_mpoly_q_denref(y), F->mctx);
    *status = _gr_tower_lazy_one(&acc, ctx);   /* res may alias x */

    for (j = 0; j < T->num_gens && ok && *status == GR_SUCCESS; j++)
    {
        const gr_tower_gen_struct * gen = GR_TOWER_GEN(T, j);
        slong v;

        if (gen->kind != GR_TOWER_LOG || gen->arg.mctx == NULL)
            continue;
        if (fmpz_mpoly_q_get_fmpq(p, &gen->arg.data, gen->arg.mctx))
        {
            if (fmpq_sgn(p) <= 0)
                continue;
            algebraic = 0;
        }
        else
        {
            /* the logarithm of an algebraic number */
            _gr_tower_lazy_set_flat(&w, F, &gen->arg.data, gen->arg.mctx, ctx);
            if (_gr_tower_lazy_is_algebraic_repr(&w) != T_TRUE)
                continue;
            algebraic = 1;
        }

        v = GR_TOWER_FLAT_VAR_D(F, j);
        if (fmpz_mpoly_degree_si(fmpz_mpoly_q_numref(y), v, F->mctx) <= 0)
            continue;
        if (fmpz_mpoly_degree_si(fmpz_mpoly_q_numref(y), v, F->mctx) != 1)
        {
            ok = 0;
            break;
        }
        fmpz_mpoly_derivative(d, fmpz_mpoly_q_numref(y), v, F->mctx);
        if (!fmpz_mpoly_is_fmpz(d, F->mctx))
        {
            ok = 0;
            break;
        }
        fmpz_mpoly_get_fmpz(fmpq_numref(c), d, F->mctx);
        fmpz_set(fmpq_denref(c), D);
        fmpq_canonicalise(c);
        if (!fmpz_fits_si(fmpq_numref(c)) || !fmpz_abs_fits_ui(fmpq_denref(c)))
        {
            ok = 0;
            break;
        }

        /* y -= c log(p) */
        fmpz_mpoly_gen(g, v, F->mctx);
        fmpz_mpoly_mul(g, g, d, F->mctx);
        fmpz_mpoly_sub(fmpz_mpoly_q_numref(y), fmpz_mpoly_q_numref(y), g, F->mctx);

        /* res *= p^c */
        if (algebraic)
        {
            /* (in a copy of the tower: the root must not be adjoined to
               x's tower, in whose context y lives) */
            _gr_tower_lazy_prefix_copy(&t, &w, ctx);
            *status = _gr_tower_lazy_root_ui(&t, &t, fmpz_get_ui(fmpq_denref(c)), ctx);
        }
        else
            *status = _gr_tower_lazy_root_fmpq(&t, p, fmpz_get_ui(fmpq_denref(c)), ctx);
        if (*status == GR_SUCCESS)
            *status = gr_pow_si(&t, &t, fmpz_get_si(fmpq_numref(c)), ctx);
        if (*status == GR_SUCCESS)
            *status = _gr_tower_lazy_mul(&acc, &acc, &t, ctx);
        count++;
    }

    if (ok && *status == GR_SUCCESS)
    {
        /* the remainder must be zero or a rational multiple of pi i */
        if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(y), F->mctx))
        {
            ok = (count > 0);
        }
        else
        {
            _gr_tower_lazy_set_flat(&u, F, y, F->mctx, ctx);
            if (_gr_tower_lazy_pi_i_multiple(r, &u, ctx))
            {
                fmpq_div_2exp(r, r, 1);
                _gr_tower_lazy_fmpq_frac_part(r);
                if (fmpz_fits_si(fmpq_numref(r)) && fmpz_abs_fits_ui(fmpq_denref(r)))
                {
                    *status = _gr_tower_lazy_root_of_unity(&t, fmpz_get_si(fmpq_numref(r)), fmpz_get_ui(fmpq_denref(r)), ctx);
                    if (*status == GR_SUCCESS)
                        *status = _gr_tower_lazy_mul(&acc, &acc, &t, ctx);
                    have_pi_i = 1;
                }
                else
                    ok = 0;
            }
            else
                ok = 0;
        }
    }

    ok = ok && (count > 0 || have_pi_i);
    if (ok && *status == GR_SUCCESS)
        _gr_tower_lazy_swap(res, &acc, ctx);

    fmpz_mpoly_clear(d, F->mctx);
    fmpz_mpoly_clear(g, F->mctx);
    fmpz_mpoly_q_clear(y, F->mctx);
    fmpz_clear(D);
    fmpq_clear(c);
    fmpq_clear(r);
    fmpq_clear(p);
    _gr_tower_lazy_clear(&t, ctx);
    _gr_tower_lazy_clear(&u, ctx);
    _gr_tower_lazy_clear(&w, ctx);
    _gr_tower_lazy_clear(&acc, ctx);
    return ok;
}

/* If x = exp(2 pi i r) with -1/2 < r <= 1/2 of denominator at most
   the option GR_TOWER_OPT_ROOT_OF_UNITY_ORDER_LIMIT (the angle is
   recognized numerically as the simplest fraction of a turn in the
   enclosure, and the identity verified exactly), sets r and returns 1. */
static int
_gr_tower_lazy_root_of_unity_angle(fmpq_t r, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    acb_t z;
    arb_t t, a;
    slong prec = 128;
    int found = 0;

    acb_init(z);
    arb_init(t);
    arb_init(a);

    if (gr_tower_lazy_get_acb(z, x, prec, ctx) == GR_SUCCESS && acb_rel_accuracy_bits(z) > 60)
    {
        /* |x| = 1 and the angle a rational number with a small denominator */
        acb_abs(t, z, prec);
        arb_sub_ui(t, t, 1, prec);
        if (arb_contains_zero(t) && mag_cmp_2exp_si(arb_radref(t), -50) < 0)
        {
            arb_t pi;
            fmpq_t lo, hi;
            fmpz_t e;

            arb_init(pi);
            fmpq_init(lo);
            fmpq_init(hi);
            fmpz_init(e);
            acb_arg(a, z, prec);
            arb_const_pi(pi, prec);
            arb_div(a, a, pi, prec);
            arb_mul_2exp_si(a, a, -1);       /* a = arg / (2 pi), in (-1/2, 1/2] */

            /* the simplest fraction in the enclosure of a */
            arb_get_interval_fmpz_2exp(fmpq_numref(lo), fmpq_numref(hi), e, a);
            fmpz_one(fmpq_denref(lo));
            fmpz_one(fmpq_denref(hi));
            if (fmpz_sgn(e) >= 0)
            {
                fmpq_mul_2exp(lo, lo, fmpz_get_si(e));
                fmpq_mul_2exp(hi, hi, fmpz_get_si(e));
            }
            else
            {
                fmpq_div_2exp(lo, lo, -fmpz_get_si(e));
                fmpq_div_2exp(hi, hi, -fmpz_get_si(e));
            }
            fmpq_simplest_between(r, lo, hi);
            found = (fmpz_cmp_ui(fmpq_denref(r), LAZY(ctx)->options[GR_TOWER_OPT_ROOT_OF_UNITY_ORDER_LIMIT]) <= 0);

            arb_clear(pi);
            fmpq_clear(lo);
            fmpq_clear(hi);
            fmpz_clear(e);
        }
    }

    if (found)
    {
        /* verify x = exp(2 pi i r) exactly */
        gr_tower_lazy_elem_struct u;
        int st;

        _gr_tower_lazy_init(&u, ctx);
        st = _gr_tower_lazy_root_of_unity(&u, fmpz_get_si(fmpq_numref(r)), fmpz_get_ui(fmpq_denref(r)), ctx);
        if (st != GR_SUCCESS || _gr_tower_lazy_equal(&u, x, ctx) != T_TRUE)
            found = 0;
        _gr_tower_lazy_clear(&u, ctx);
    }

    acb_clear(z);
    arb_clear(t);
    arb_clear(a);
    return found;
}

int
_gr_tower_lazy_root_of_unity_angle_locked(fmpq_t r, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct t;
    int found;
    _gr_tower_lazy_lock(ctx);
    _gr_tower_lazy_init(&t, ctx);
    found = (_gr_tower_lazy_set(&t, x, ctx) == GR_SUCCESS) && _gr_tower_lazy_root_of_unity_angle(r, &t, ctx);
    _gr_tower_lazy_clear(&t, ctx);
    _gr_tower_lazy_unlock(ctx);
    return found;
}

/*
    res = E(zeta_m) for E in Q[x], formed directly in the flat
    representation of the tower of zeta_m (the powers of zeta_m are
    monomials in the prime power roots; no field operations)
*/
int
_gr_tower_lazy_cyclotomic_eval(gr_tower_lazy_elem_t res, const fmpq_poly_t E, ulong m, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct z, r;
    gr_tower_flat_struct * F;
    fmpz_mpoly_q_t acc;
    int status;

    _gr_tower_lazy_lock(ctx);
    _gr_tower_lazy_init(&z, ctx);
    _gr_tower_lazy_init(&r, ctx);

    status = _gr_tower_lazy_root_of_unity(&z, 1, m, ctx);

    if (status == GR_SUCCESS)
    {
        _gr_tower_lazy_update(&z);
        F = z.F;
        fmpz_mpoly_q_init(acc, F->mctx);
        /* (the prime power roots are now present) */
        if (!_gr_tower_special_flat_cyclotomic_eval(acc, E, m, F))
            status = GR_UNABLE;

        if (status == GR_SUCCESS)
            _gr_tower_lazy_set_flat(&r, F, acc, F->mctx, ctx);

        fmpz_mpoly_q_clear(acc, F->mctx);
    }

    if (status == GR_SUCCESS)
        _gr_tower_lazy_swap((gr_tower_lazy_elem_struct *) res, &r, ctx);

    _gr_tower_lazy_clear(&z, ctx);
    _gr_tower_lazy_clear(&r, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
_gr_tower_lazy_log_root_of_unity(gr_tower_lazy_elem_t res, gr_tower_lazy_elem_t x, int * status, gr_ctx_t ctx)
{
    fmpq_t r;
    int found;

    fmpq_init(r);
    found = _gr_tower_lazy_root_of_unity_angle(r, x, ctx);

    if (found)
    {
        /* res = 2 r pi i */
        gr_tower_lazy_elem_struct u;
        int st;

        _gr_tower_lazy_init(&u, ctx);
        st = _gr_tower_lazy_pi(res, ctx);
        if (st == GR_SUCCESS)
            st = _gr_tower_lazy_i(&u, ctx);
        if (st == GR_SUCCESS)
            st = _gr_tower_lazy_mul(res, res, &u, ctx);
        if (st == GR_SUCCESS)
        {
            fmpq_mul_2exp(r, r, 1);
            fmpz_mpoly_q_mul_fmpq(&res->elem.flat.data, &res->elem.flat.data, r, res->elem.flat.mctx);
        }
        *status = st;
        _gr_tower_lazy_clear(&u, ctx);
    }

    fmpq_clear(r);
    return found;
}

/*
    The transcendental generator kind(x) (exp, log, tan, atan) for x:
    the generator of x's tower with that definition if there is one,
    otherwise a new one (followed by a cheap relation search, which may
    eliminate it at once: exp(log(x)) = x, tan(2x) over tan(x)).

    When x lives in a prefix of its tower which the tower continues with
    at least LAZY_FORK_SUFFIX more generators, the new generator goes
    into a fresh copy of that prefix instead (keeping the definition id
    and name of an existing generator with that definition), so that
    computations not involving the rest of the tower do not pay for its
    generators and relations (a long session accumulates them in shared
    towers); merges recognize the shared generators by their definitions.
*/
static int
_gr_tower_lazy_trans_gen_param(gr_tower_lazy_elem_t res, gr_tower_lazy_elem_t x, int kind, slong param, gr_ctx_t ctx);

int
_gr_tower_lazy_trans_gen(gr_tower_lazy_elem_t res, gr_tower_lazy_elem_t x, int kind, gr_ctx_t ctx)
{
    return _gr_tower_lazy_trans_gen_param(res, x, kind, 0, ctx);
}

/* the generator of the tower of x (by its gid) kind(param) at x, or -1
   (exact comparisons, which may restructure the tower) */
static slong
_trans_gen_find(gr_tower_lazy_elem_t x, int kind, slong param)
{
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    slong j, found_gid = -1;
    truth_t t;

    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;

    for (j = 0; j < T->num_gens && found_gid < 0; j++)
    {
        gr_tower_gen_struct * g = GR_TOWER_GEN(T, j);
        fmpz_mpoly_q_t d;
        slong gid = g->gid;

        if (!((g->kind == kind) || (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == kind)) || g->arg.mctx == NULL)
            continue;
        if (GR_TOWER_KIND_IS_SPECIAL(kind) && g->def_param != param)
            continue;

        /* zero tests may restructure the tower: keep x current */
        x = _gr_tower_lazy_flat_view(x);
        fmpz_mpoly_q_init(d, F->mctx);
        gr_tower_flat_convert(d, &g->arg.data, g->arg.mctx, F);
        fmpz_mpoly_q_sub(d, d, &x->elem.flat.data, F->mctx);
        t = gr_tower_flat_num_is_zero(d, F);
        fmpz_mpoly_q_clear(d, F->mctx);

        if (t == T_TRUE)
            found_gid = gid;
    }

    /* (a relation found by a later zero test may have eliminated it) */
    if (found_gid >= 0 && gr_tower_gid_order(T, found_gid) < 0)
        found_gid = -1;

    return found_gid;
}

static int
_gr_tower_lazy_trans_gen_param(gr_tower_lazy_elem_t res, gr_tower_lazy_elem_t x, int kind, slong param, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    slong found_gid = -1;
    ulong def_id = 0;
    char * name = NULL;
    int status;

    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;
    found_gid = _trans_gen_find(x, kind, param);

    /* the same definition in another tower */
    if (found_gid < 0)
    {
        def_id = _gr_tower_lazy_find_def(&name, kind, param, x, 1, T, ctx);
        x = _gr_tower_lazy_flat_view(x);
        /* (the comparisons may have merged x into another tower: a
           generator there with this definition) */
        if (x->F != F)
        {
            F = x->F;
            T = F->T;
            found_gid = _trans_gen_find(x, kind, param);
            if (found_gid >= 0 && name != NULL)
            {
                flint_free(name);
                name = NULL;
                def_id = 0;
            }
        }
    }

    x = _gr_tower_lazy_flat_view(x);

    if (T->num_gens - x->level < LAZY_FORK_SUFFIX || x->level == 0)
    {
        if (found_gid >= 0)
        {
            /* (the zero tests may have restructured the tower: the
               generator is identified by its gid) */
            _gr_tower_lazy_set_gen_d(res, F, gr_tower_gid_order(T, found_gid), ctx);
            return GR_SUCCESS;
        }
    }
    else
    {
        /* a copy of the prefix */
        gr_tower_flat_struct * G;
        fmpz_mpoly_q_t u;

        if (found_gid >= 0)
        {
            const gr_tower_gen_struct * g = GR_TOWER_GEN(T, gr_tower_gid_order(T, found_gid));
            def_id = g->def_id;
            name = flint_malloc(strlen(g->name) + 1);
            strcpy(name, g->name);
        }

        {
            gr_tower_lazy_elem_struct y;
            _gr_tower_lazy_init(&y, ctx);
            _gr_tower_lazy_prefix_copy(&y, x, ctx);
            G = y.F;
            _gr_tower_lazy_update(&y);
            fmpz_mpoly_q_init(u, G->mctx);
            fmpz_mpoly_q_set(u, &y.elem.flat.data, G->mctx);
            _gr_tower_lazy_clear(&y, ctx);
        }

        F = G;
        T = G->T;
        if (GR_TOWER_KIND_IS_SPECIAL(kind))
            status = gr_tower_adjoin_special_flat(T, kind, param, u, G->mctx, NULL);
        else if (kind == GR_TOWER_EXP)
            status = gr_tower_adjoin_exp_flat(T, u, G->mctx, NULL);
        else if (kind == GR_TOWER_LOG)
            status = gr_tower_adjoin_log_flat(T, u, G->mctx, NULL);
        else if (kind == GR_TOWER_TAN)
            status = gr_tower_adjoin_tan_flat(T, u, G->mctx, NULL);
        else
            status = gr_tower_adjoin_atan_flat(T, u, G->mctx, NULL);
        fmpz_mpoly_q_clear(u, G->mctx);   /* (the context of G at the time) */

        if (status == GR_SUCCESS)
        {
            gr_tower_gen_struct * ng = GR_TOWER_GEN(T, T->num_gens - 1);
            slong gid = ng->gid;
            _gr_tower_lazy_set_def(ng, T, def_id, name, ctx);
            _gr_tower_search_relations(F, GR_TOWER_DEFAULT_PREC);
            _gr_tower_lazy_set_gen_d(res, F, gr_tower_gid_order(T, gid), ctx);
        }
        flint_free(name);
        return status;
    }

    if (GR_TOWER_KIND_IS_SPECIAL(kind))
        status = gr_tower_adjoin_special_flat(T, kind, param, &x->elem.flat.data, x->elem.flat.mctx, NULL);
    else if (kind == GR_TOWER_EXP)
        status = gr_tower_adjoin_exp_flat(T, &x->elem.flat.data, x->elem.flat.mctx, NULL);
    else if (kind == GR_TOWER_LOG)
        status = gr_tower_adjoin_log_flat(T, &x->elem.flat.data, x->elem.flat.mctx, NULL);
    else if (kind == GR_TOWER_TAN)
        status = gr_tower_adjoin_tan_flat(T, &x->elem.flat.data, x->elem.flat.mctx, NULL);
    else
        status = gr_tower_adjoin_atan_flat(T, &x->elem.flat.data, x->elem.flat.mctx, NULL);

    if (status != GR_SUCCESS)
    {
        flint_free(name);
        return status;
    }

    {
        slong gid = GR_TOWER_GEN(T, T->num_gens - 1)->gid;

        _gr_tower_lazy_set_def(GR_TOWER_GEN(T, T->num_gens - 1), T, def_id, name, ctx);
        flint_free(name);

        /* a cheap relation search at adjunction time keeps the towers
           small (exp(log(x)) becomes x at once); it may reorder the
           generators */
        _gr_tower_search_relations(F, GR_TOWER_DEFAULT_PREC);

        _gr_tower_lazy_set_gen_d(res, F, gr_tower_gid_order(T, gid), ctx);
    }
    return GR_SUCCESS;
}

int
_gr_tower_lazy_exp_log(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, int kind, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_lazy_elem_struct * x = (gr_tower_lazy_elem_struct *) x_in;
    truth_t t;
    int status, alg;

    x = _gr_tower_lazy_flat_view(x);

    /* a rational number represented in some tower (pi/pi, say) is
       treated as a rational number */
    if (x->F != L->trivial &&
        fmpz_mpoly_q_is_fmpq(&x->elem.flat.data, x->elem.flat.mctx))
    {
        gr_tower_lazy_elem_struct t;
        fmpq_t c;
        fmpq_init(c);
        (void) fmpz_mpoly_q_get_fmpq(c, &x->elem.flat.data, x->elem.flat.mctx);
        _gr_tower_lazy_init(&t, ctx);
        fmpq_set(&t.elem.q, c);
        status = _gr_tower_lazy_exp_log(res, &t, kind, ctx);
        _gr_tower_lazy_clear(&t, ctx);
        fmpq_clear(c);
        return status;
    }

    /* log(c exp(u)) = log(c) + u for a positive rational c and
       Im(u) in (-pi, pi]: decided from the definition, since exp(u) may
       be numerically indistinguishable from 1 (or 0) */
    if (kind == GR_TOWER_LOG && x->F != L->trivial)
    {
        fmpq_t c;
        slong d;
        int st;

        fmpq_init(c);
        d = _gr_tower_lazy_exp_gen_multiple(c, x, ctx);
        if (d >= 0)
        {
            gr_tower_lazy_elem_struct u;
            const gr_tower_gen_struct * g = GR_TOWER_GEN(x->F->T, d);

            _gr_tower_lazy_init(&u, ctx);
            _gr_tower_lazy_set_flat(&u, x->F, &g->arg.data, g->arg.mctx, ctx);
            if (fmpq_is_one(c))
                st = _gr_tower_lazy_set(res, &u, ctx);
            else
            {
                fmpq_t zero;
                fmpq_init(zero);
                st = _gr_tower_lazy_log_gaussian(res, c, zero, ctx);
                if (st == GR_SUCCESS)
                    st = _gr_tower_lazy_add(res, res, &u, ctx);
                fmpq_clear(zero);
            }
            _gr_tower_lazy_clear(&u, ctx);
            fmpq_clear(c);
            return st;
        }
        fmpq_clear(c);
    }

    /* trivial values */
    t = _gr_tower_lazy_is_zero(x, ctx);
    if (t == T_UNKNOWN)
        return GR_UNABLE;
    if (t == T_TRUE)
        return (kind == GR_TOWER_EXP) ? _gr_tower_lazy_one(res, ctx) : GR_DOMAIN;

    if (kind == GR_TOWER_LOG)
    {
        t = _gr_tower_lazy_is_one(x, ctx);
        if (t == T_UNKNOWN)
            return GR_UNABLE;
        if (t == T_TRUE)
            return _gr_tower_lazy_zero(res, ctx);
    }

    /* constants are hash-consed; logarithms of rational numbers are
       decomposed over the logarithms of primes */
    if (x->F == L->trivial)
    {
        fmpq_t c;
        fmpq_init(c);
        (void) fmpz_mpoly_q_get_fmpq(c, &x->elem.flat.data, x->elem.flat.mctx);
        if (kind == GR_TOWER_LOG)
        {
            fmpq_t zero;
            fmpq_init(zero);
            status = _gr_tower_lazy_log_gaussian(res, c, zero, ctx);
            fmpq_clear(zero);
        }
        else
            status = _gr_tower_lazy_const_trans(res, kind, c, ctx);
        fmpq_clear(c);
        return status;
    }

    /* x may live in a tower with transcendental generators without
       using them (i parsed in a context whose generators include a
       logarithm, say) */
    alg = (kind == GR_TOWER_LOG) && (x->F->T->num_trans == 0 ||
                                     _gr_tower_lazy_is_algebraic_repr(x) == T_TRUE);

    /* logarithms of Gaussian rationals likewise */
    if (kind == GR_TOWER_LOG && alg)
    {
        fmpq_t u, v;
        int st;

        fmpq_init(u); fmpq_init(v);
        st = _gr_tower_lazy_get_gaussian(u, v, x, ctx);
        if (st == GR_SUCCESS)
            st = _gr_tower_lazy_log_gaussian(res, u, v, ctx);
        fmpq_clear(u); fmpq_clear(v);
        if (st == GR_SUCCESS)
            return st;
    }

    /* the logarithm of a negative real number whose enclosure does not
       have an exactly zero imaginary part (a real number expressed
       through complex generators, cos(5/3) say) would be evaluated
       numerically across the branch cut: log(x) = log(-x) + pi i */
    if (kind == GR_TOWER_LOG && L->conj_depth == 0)
    {
        acb_t z;
        int cut = 0;
        acb_init(z);
        if (gr_tower_lazy_get_acb(z, x, GR_TOWER_DEFAULT_PREC, ctx) == GR_SUCCESS &&
            arb_is_negative(acb_realref(z)) && arb_contains_zero(acb_imagref(z)) && !arb_is_zero(acb_imagref(z)))
            cut = (gr_tower_lazy_is_real(x, ctx) == T_TRUE);
        acb_clear(z);

        if (cut)
        {
            gr_tower_lazy_elem_struct t;
            _gr_tower_lazy_init(&t, ctx);
            status = _gr_tower_lazy_neg(&t, x, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_exp_log(res, &t, GR_TOWER_LOG, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_pi(&t, ctx);
            if (status == GR_SUCCESS)
            {
                gr_tower_lazy_elem_struct ii;
                _gr_tower_lazy_init(&ii, ctx);
                status = _gr_tower_lazy_i(&ii, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_mul(&t, &t, &ii, ctx);
                _gr_tower_lazy_clear(&ii, ctx);
            }
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_add(res, res, &t, ctx);
            _gr_tower_lazy_clear(&t, ctx);
            return status;
        }
    }

    /* the logarithm of a root of unity (recognized numerically, then
       verified exactly: (sqrt(3) + i)/2, say) is a rational multiple
       of pi i */
    if (kind == GR_TOWER_LOG && alg)
    {
        int st;
        if (_gr_tower_lazy_log_root_of_unity(res, x, &st, ctx))
            return st;
    }

    /* exponentials of rational combinations of logarithms of positive
       rationals and of pi i are products of radicals and roots of unity */
    if (kind == GR_TOWER_EXP)
    {
        int st;
        if (_gr_tower_lazy_exp_of_logs(res, x, &st, ctx))
            return st;
    }

    /* exp(r pi i + y) = exp(r pi i) exp(y), a root of unity times an
       exponential whose argument does not involve pi (otherwise the
       relations between exp(y) and exp(r pi i + y) are radicals of high
       degree), when that root of unity is already in the tower (a new
       one would increase its degree) */
    if (kind == GR_TOWER_EXP)
    {
        gr_tower_lazy_elem_struct y, z;
        fmpq_t r;
        int st = 0;

        fmpq_init(r);
        _gr_tower_lazy_init(&y, ctx);
        _gr_tower_lazy_init(&z, ctx);
        if (_gr_tower_lazy_pi_i_split(r, &y, x, ctx) && _gr_tower_lazy_is_zero(&y, ctx) == T_FALSE &&
            (fmpq_div_2exp(r, r, 1), _gr_tower_lazy_fmpq_frac_part(r), fmpz_fits_si(fmpq_numref(r)) && fmpz_abs_fits_ui(fmpq_denref(r))) &&
            _gr_tower_lazy_has_root_of_unity(x->F->T, fmpz_get_ui(fmpq_denref(r))))
        {
            status = _gr_tower_lazy_root_of_unity(&z, fmpz_get_si(fmpq_numref(r)), fmpz_get_ui(fmpq_denref(r)), ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_exp_log(res, &y, GR_TOWER_EXP, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_mul(res, res, &z, ctx);
            st = 1;
        }
        fmpq_clear(r);
        _gr_tower_lazy_clear(&y, ctx);
        _gr_tower_lazy_clear(&z, ctx);
        if (st)
            return status;
    }

    /* exponentials of rational multiples of pi i are roots of unity */
    if (kind == GR_TOWER_EXP)
    {
        fmpq_t r;
        int st = 0;

        fmpq_init(r);
        if (_gr_tower_lazy_pi_i_multiple(r, x, ctx))
        {
            /* exp(r pi i) = exp(2 pi i (r/2)) */
            fmpq_div_2exp(r, r, 1);
            _gr_tower_lazy_fmpq_frac_part(r);
            if (fmpz_fits_si(fmpq_numref(r)) && fmpz_abs_fits_ui(fmpq_denref(r)))
            {
                status = _gr_tower_lazy_root_of_unity(res, fmpz_get_si(fmpq_numref(r)), fmpz_get_ui(fmpq_denref(r)), ctx);
                st = 1;
            }
        }
        fmpq_clear(r);
        if (st)
            return status;
    }

    return _gr_tower_lazy_trans_gen(res, x, kind, ctx);
}

/*
    tan(x) or atan(x) (kind GR_TOWER_TAN or GR_TOWER_ATAN) as a real
    generator, for a real x (which the caller ensures), x not an odd
    multiple of pi/2 for tan: hash-consed for rational x, otherwise the
    generator of x's tower with that argument, or a new one (followed by
    a cheap relation search, which finds for instance tan(2x) in terms
    of tan(x) when both are present).
*/
static int
_gr_tower_lazy_trig(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, int kind, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * x = (gr_tower_lazy_elem_struct *) x_in;
    truth_t t;
    int status;

    x = _gr_tower_lazy_flat_view(x);

    t = _gr_tower_lazy_is_zero(x, ctx);
    if (t == T_UNKNOWN)
        return GR_UNABLE;
    if (t == T_TRUE)
        return _gr_tower_lazy_zero(res, ctx);

    /* odd functions: the generator of a positive argument */
    {
        acb_t z;
        int neg;
        acb_init(z);
        neg = (_gr_tower_lazy_get_acb_impl(z, x, GR_TOWER_DEFAULT_PREC, ctx) == GR_SUCCESS) && arb_is_negative(acb_realref(z));
        acb_clear(z);
        if (neg)
        {
            gr_tower_lazy_elem_struct y;
            _gr_tower_lazy_init(&y, ctx);
            status = _gr_tower_lazy_neg(&y, x, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_trig(res, &y, kind, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_neg(res, res, ctx);
            _gr_tower_lazy_clear(&y, ctx);
            return status;
        }
    }

    /* rational arguments (possibly represented in some tower) */
    if (fmpz_mpoly_q_is_fmpq(&x->elem.flat.data, x->elem.flat.mctx))
    {
        fmpq_t c;
        fmpq_init(c);
        (void) fmpz_mpoly_q_get_fmpq(c, &x->elem.flat.data, x->elem.flat.mctx);
        status = _gr_tower_lazy_const_trans(res, kind, c, ctx);
        fmpq_clear(c);
        return status;
    }

    return _gr_tower_lazy_trans_gen(res, x, kind, ctx);
}

/*
    x = r pi + y with r rational and y free of pi (x affine in pi with a
    rational coefficient, pi not in the denominator): sets r and y and
    returns 1; returns 0 otherwise (r = 0 when pi does not occur).
*/
static int
_gr_tower_lazy_pi_part(fmpq_t r, gr_tower_lazy_elem_t y, const gr_tower_lazy_elem_t x_in, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * x = (gr_tower_lazy_elem_struct *) x_in;
    gr_tower_flat_struct * F;
    gr_tower_struct * T;
    slong v;
    fmpz_mpoly_t a;
    const fmpz_mpoly_struct * den;
    int ok = 0;

    fmpq_zero(r);
    x = _gr_tower_lazy_flat_view(x);
    F = x->F;
    T = F->T;
    den = fmpz_mpoly_q_denref(&x->elem.flat.data);

    if (!_pi_coeff(&v, a, x))
        return 0;

    /* the coefficient of pi: a / den must be rational */
    if (fmpz_mpoly_is_fmpz(a, x->elem.flat.mctx) && fmpz_mpoly_is_fmpz(den, x->elem.flat.mctx))
    {
        fmpz_mpoly_q_t t;
        fmpz_mpoly_get_fmpz(fmpq_numref(r), a, x->elem.flat.mctx);
        fmpz_mpoly_get_fmpz(fmpq_denref(r), den, x->elem.flat.mctx);
        fmpq_canonicalise(r);

        /* y = x - r pi */
        fmpz_mpoly_q_init(t, x->elem.flat.mctx);
        fmpz_mpoly_q_gen(t, v, x->elem.flat.mctx);
        fmpz_mpoly_q_mul_fmpq(t, t, r, x->elem.flat.mctx);
        fmpz_mpoly_q_sub(t, &x->elem.flat.data, t, x->elem.flat.mctx);
        _gr_tower_lazy_fresh(y, F, T->num_gens, ctx);
        fmpz_mpoly_q_swap(&y->elem.flat.data, t, y->elem.flat.mctx);
        _gr_tower_lazy_shrink(y);
        fmpz_mpoly_q_clear(t, x->elem.flat.mctx);
        ok = 1;
    }
    fmpz_mpoly_clear(a, x->elem.flat.mctx);
    return ok;
}

int
_gr_tower_lazy_pi_part_locked(fmpq_t r, gr_tower_lazy_elem_t y, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int res;
    _gr_tower_lazy_lock(ctx);
    res = _gr_tower_lazy_pi_part(r, y, x, ctx);
    _gr_tower_lazy_unlock(ctx);
    return res;
}

int
_gr_tower_lazy_trig_locked(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, int kind, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    status = _gr_tower_lazy_trig(res, x, kind, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* Whether x is represented by a rational constant (then set c). */
static int
_gr_tower_lazy_rational_repr(fmpq_t c, gr_tower_lazy_elem_t x)
{
    x = _gr_tower_lazy_flat_view(x);
    return fmpz_mpoly_q_get_fmpq(c, &x->elem.flat.data, x->elem.flat.mctx);
}

int
_gr_tower_lazy_rational_repr_locked(fmpq_t c, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int res;
    _gr_tower_lazy_lock(ctx);
    if (((const gr_tower_lazy_elem_struct *) x)->shallow)
    {
        gr_tower_lazy_elem_struct t;
        _gr_tower_lazy_init(&t, ctx);
        res = (_gr_tower_lazy_set(&t, x, ctx) == GR_SUCCESS) && _gr_tower_lazy_rational_repr(c, &t);
        _gr_tower_lazy_clear(&t, ctx);
    }
    else
        res = _gr_tower_lazy_rational_repr(c, (gr_tower_lazy_elem_struct *) x);
    _gr_tower_lazy_unlock(ctx);
    return res;
}

/*
    The special function generator kind(x) (with parameter param), without
    any canonicalization of the argument (see lazy_special.c): hash-consed
    for a rational x, otherwise the generator of x's tower with that
    definition, or a new one. x is ignored for named constants.
*/
static int
_gr_tower_lazy_special_gen(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, int kind, slong param, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * x = (gr_tower_lazy_elem_struct *) x_in;
    fmpq_t c, zero;
    int status;

    fmpq_init(c);
    fmpq_init(zero);

    if (kind == GR_TOWER_CONSTANT)
        status = _gr_tower_lazy_const_trans2_param(res, kind, param, zero, zero, ctx);
    else if (_gr_tower_lazy_rational_repr(c, x))
        status = _gr_tower_lazy_const_trans2_param(res, kind, param, c, zero, ctx);
    else
        status = _gr_tower_lazy_trans_gen_param(res, x, kind, param, ctx);

    fmpq_clear(c);
    fmpq_clear(zero);
    return status;
}

int
_gr_tower_lazy_special_gen_locked(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, int kind, slong param, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    if (kind != GR_TOWER_CONSTANT && ((const gr_tower_lazy_elem_struct *) x)->shallow)
    {
        gr_tower_lazy_elem_struct t;
        _gr_tower_lazy_init(&t, ctx);
        status = _gr_tower_lazy_set(&t, x, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_special_gen(res, &t, kind, param, ctx);
        _gr_tower_lazy_clear(&t, ctx);
    }
    else
        status = _gr_tower_lazy_special_gen(res, x, kind, param, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
_gr_tower_lazy_exp(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    return _gr_tower_lazy_exp_log(res, x, GR_TOWER_EXP, ctx);
}

int
_gr_tower_lazy_log(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    return _gr_tower_lazy_exp_log(res, x, GR_TOWER_LOG, ctx);
}

POP_OPTIONS
