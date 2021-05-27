/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


/**************** product of roots *******************************************/

void fq_nmod_poly_product_roots(
    fq_nmod_poly_t master,
    const fq_nmod_struct * monomials,
    slong mlength,
    const fq_nmod_ctx_t ctx);

void n_fq_poly_product_roots_n_fq(
    n_poly_t master,
    const mp_limb_t * monomials,
    slong mlength,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    fq_nmod_poly_t p;
    fq_nmod_struct * m = FLINT_ARRAY_ALLOC(mlength, fq_nmod_struct);

    fq_nmod_poly_init(p, ctx);
    for (i = 0; i < mlength; i++)
    {
        fq_nmod_init(m + i, ctx);
        n_fq_get_fq_nmod(m + i, monomials + d*i, ctx);
    }

    fq_nmod_poly_product_roots(p, m, mlength, ctx);

    n_fq_poly_set_fq_nmod_poly(master, p, ctx);

    fq_nmod_poly_clear(p, ctx);
    for (i = 0; i < mlength; i++)
        fq_nmod_clear(m + i, ctx);
    flint_free(m);
}

/* return the largest degree */
slong n_polyun_product_roots(
    n_polyun_t M,
    const n_polyun_t H,
    nmod_t ctx)
{
    slong i, max_length = 0;

    n_polyun_fit_length(M, H->length);
    M->length = H->length;
    for (i = 0; i < H->length; i++)
    {
        slong len = H->coeffs[i].length;
        M->exps[i] = H->exps[i];
        max_length = FLINT_MAX(max_length, len);
        n_poly_mod_product_roots_nmod_vec(M->coeffs + i, H->coeffs[i].coeffs, len, ctx);
    }

    return max_length;
}

slong n_fq_polyun_product_roots(
    n_fq_polyun_t M,
    const n_fq_polyun_t H,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong i, max_length = 0;

    n_polyun_fit_length(M, H->length);
    M->length = H->length;
    for (i = 0; i < H->length; i++)
    {
        slong len = H->coeffs[i].length;
        M->exps[i] = H->exps[i];
        max_length = FLINT_MAX(max_length, len);
        n_fq_poly_product_roots_n_fq(M->coeffs + i, H->coeffs[i].coeffs, len, ctx, St);
    }

    return max_length;
}


/******************** evaluation *********************************************/
/*
    return dot(cur, coeffs)
    and multiply cur pointwise by inc
*/

mp_limb_t _nmod_zip_eval_step(
    mp_limb_t * cur,            /* in Fp */
    const mp_limb_t * inc,      /* in Fp */
    const mp_limb_t * coeffs,   /* in Fp */
    slong length,
    nmod_t ctx)
{
    slong i;
    ulong t0, t1, t2, p0, p1;
    t2 = t1 = t0 = 0;
    for (i = 0; i < length; i++)
    {
        umul_ppmm(p1, p0, cur[i], coeffs[i]);
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
        cur[i] = nmod_mul(cur[i], inc[i], ctx);
    }
    NMOD_RED3(t0, t2, t1, t0, ctx);
    return t0;
}


void _n_fq_zip_eval_step(
    mp_limb_t * res,            /* in Fq: size d */
    mp_limb_t * cur,            /* in Fq: size d*length */
    const mp_limb_t * inc,      /* in Fq: size d*length */
    const mp_limb_t * coeffs,   /* in Fq: size d*length */
    slong length,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    mp_limb_t * tmp, * sum;
    TMP_INIT;

    if (length < 1)
    {
        _n_fq_zero(res, d);
        return;
    }

    TMP_START;
    tmp = (mp_limb_t *) TMP_ALLOC(8*d*sizeof(mp_limb_t));
    sum = tmp + 4*d;

    i = 0;
    _n_fq_mul2(sum, cur + d*i, coeffs + d*i, ctx);
    _n_fq_mul(cur + d*i, cur + d*i, inc + d*i, ctx, tmp);
    for (i = 1; i < length; i++)
    {
        _n_fq_madd2(sum, cur + d*i, coeffs + d*i, ctx, tmp);
        _n_fq_mul(cur + d*i, cur + d*i, inc + d*i, ctx, tmp);
    }
    _n_fq_reduce2(res, sum, ctx, tmp);

    TMP_END;
}


void _n_fqp_zip_eval_step(
    mp_limb_t * res,            /* in Fq: size d */
    mp_limb_t * cur,            /* in Fp: size length */
    const mp_limb_t * inc,      /* in Fp: size length */
    const mp_limb_t * coeffs,   /* in Fq: size d*length */
    slong length,
    slong d,
    nmod_t mod)
{
    slong i, j;
    mp_limb_t p0, p1;
    mp_limb_t * tmp;
    TMP_INIT;

    if (length < 1)
    {
        _n_fq_zero(res, d);
        return;
    }

    TMP_START;
    tmp = (mp_limb_t *) TMP_ALLOC(3*d*sizeof(mp_limb_t));

    i = 0;

    for (j = 0; j < d; j++)
    {
        umul_ppmm(tmp[3*j+1], tmp[3*j+0], cur[i], (coeffs + d*i)[j]);
        tmp[3*j+2] = 0;
    }
    cur[i] = nmod_mul(cur[i], inc[i], mod);

    for (i = 1; i < length; i++)
    {
        for (j = 0; j < d; j++)
        {
            umul_ppmm(p1, p0, cur[i], (coeffs + d*i)[j]);
            add_sssaaaaaa(tmp[3*j+2], tmp[3*j+1], tmp[3*j+0],
                          tmp[3*j+2], tmp[3*j+1], tmp[3*j+0], 0, p1, p0);
        }
        cur[i] = nmod_mul(cur[i], inc[i], mod);
    }

    for (j = 0; j < d; j++)
        NMOD_RED3(res[j], tmp[3*j+2], tmp[3*j+1], tmp[3*j+0], mod);

    TMP_END;
}


/************** vandermonde solving ******************************************/

/*
    return 
        -1: singular
        0:  inconsistent
        1:  success
*/
int _nmod_zip_vand_solve(
    mp_limb_t * coeffs,             /* in Fp: size mlength */
    const mp_limb_t * monomials,    /* in Fp: size mlength */
    slong mlength,
    const mp_limb_t * evals,        /* in Fp: size elength */
    slong elength,
    const mp_limb_t * master,       /* in Fp: size mlength + 1 */
    mp_limb_t * scratch,            /* in Fp: size mlength */
    nmod_t ctx)
{
    slong i, j;
    mp_limb_t V, V0, V1, V2, T, S, r, p0, p1;

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = monomials[i];
        for (j = mlength; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, ctx), master[j], ctx);
            S = nmod_add(nmod_mul(r, S, ctx), T, ctx);
            umul_ppmm(p1, p0, evals[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, 0, p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, ctx), master[0], ctx) == 0);
        NMOD_RED3(V, V2, V1, V0, ctx);
        S = nmod_mul(S, r, ctx); /* shift is one */
        if (S == 0)
            return -1;
        coeffs[i] = nmod_mul(V, nmod_inv(S, ctx), ctx);
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        scratch[j] = nmod_pow_ui(monomials[j], mlength, ctx);

    for (i = mlength; i < elength; i++)
    {
        V0 = V1 = V2 = S = 0;
        for (j = 0; j < mlength; j++)
        {
            scratch[j] = nmod_mul(scratch[j], monomials[j], ctx);
            umul_ppmm(p1, p0, coeffs[j], scratch[j]);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, 0, p1, p0);
        }
        NMOD_RED3(V, V2, V1, V0, ctx);
        if (V != evals[i])
            return 0;
    }
    return 1;
}

int _n_fq_zip_vand_solve(
    mp_limb_t * coeffs,             /* in Fq: size d*mlength */
    const mp_limb_t * monomials,    /* in Fq: size d*mlength */
    slong mlength,
    const mp_limb_t * evals,        /* in Fq: size d*elength */
    slong elength,
    const mp_limb_t * master,       /* in Fq: size d*(mlength + 1) */
    mp_limb_t * scratch,            /* in Fq: size d*mlength */
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    int success;
    slong i, j;
    mp_limb_t * tmp = FLINT_ARRAY_ALLOC(12*d, mp_limb_t);
    mp_limb_t * V = tmp + 6*d;
    mp_limb_t * V0 = V + d;
    mp_limb_t * T = V0 + d;
    mp_limb_t * S = T + d;
    mp_limb_t * r = S + d;
    mp_limb_t * p0 = r + d;

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        _n_fq_zero(V0, d);
        _n_fq_zero(T, d);
        _n_fq_zero(S, d);
        _n_fq_set(r, monomials + d*i, d);
        for (j = mlength; j > 0; j--)
        {
            _n_fq_mul(T, r, T, ctx, tmp);
            _n_fq_add(T, T, master + d*j, d, mod);

            _n_fq_mul(S, r, S, ctx, tmp);
            _n_fq_add(S, S, T, d, mod);

            _n_fq_mul(p0, evals + d*(j - 1), T, ctx, tmp);
            _n_fq_add(V0, V0, p0, d, mod);
        }
        /* roots[i] should be a root of master */
#if FLINT_WANT_ASSERT
        _n_fq_mul(p0, r, T, ctx, tmp);
        _n_fq_add(p0, p0, master + d*0, d, mod);
        FLINT_ASSERT(_n_fq_is_zero(p0, d));
#endif
        _n_fq_set(V, V0, d);
        _n_fq_mul(S, S, r, ctx, tmp);
        if (_n_fq_is_zero(S, d))
        {
            success = -1;
            goto cleanup;
        }

        _n_fq_inv(p0, S, ctx, tmp);
        _n_fq_mul(p0, V, p0, ctx, tmp);
        _n_fq_set(coeffs + d*i, p0, d);
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
    {
        _n_fq_set(p0, monomials + d*j, d);
        _n_fq_pow_ui(scratch + d*j, p0, mlength, ctx);
    }

    for (i = mlength; i < elength; i++)
    {
        _n_fq_zero(V0, d);
        _n_fq_zero(S, d);
        for (j = 0; j < mlength; j++)
        {
            _n_fq_set(p0, monomials + d*j, d);
            _n_fq_mul(scratch + d*j, scratch + d*j, p0, ctx, tmp);
            _n_fq_set(p0, coeffs + d*j, d);
            _n_fq_mul(p0, p0, scratch + d*j, ctx, tmp);
            _n_fq_add(V0, V0, p0, d, mod);
        }
        _n_fq_set(V, V0, d);
        if (!_n_fq_equal(V, evals + d*i, d))
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    flint_free(tmp);

    return success;
}


int _n_fqp_zip_vand_solve(
    mp_limb_t * coeffs,             /* in Fq: size d*mlength */
    const mp_limb_t * monomials,    /* in Fp: size mlength */
    slong mlength,
    const mp_limb_t * evals,        /* in Fq: size d*elength */
    slong elength,
    const mp_limb_t * master,       /* in Fp: size (mlength + 1) */
    mp_limb_t * scratch,            /* in Fp: size mlength */
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    int success;
    slong i, j, k;
    mp_limb_t * tmp = FLINT_ARRAY_ALLOC(d*20, mp_limb_t);
    mp_limb_t * V = tmp + 6*d;
    mp_limb_t * V0 = V + d;
    mp_limb_t * T = V0 + d;
    mp_limb_t * S = T + d;
    mp_limb_t * r = S + d;
    mp_limb_t * p0 = r + d;
    mp_limb_t * V_p = p0 + d;
    mp_limb_t r_p, T_p, S_p;

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        _n_fq_zero(V0, d);
        _n_fq_zero(T, d);
        _n_fq_zero(S, d);

        _nmod_vec_zero(V_p, 3*d);
        T_p = S_p = 0;
        r_p = monomials[i];
        for (j = mlength; j > 0; j--)
        {
            T_p = nmod_add(nmod_mul(r_p, T_p, mod), master[j], mod);
            S_p = nmod_add(nmod_mul(r_p, S_p, mod), T_p, mod);
            for (k = 0; k < d; k++)
            {
                mp_limb_t hi, lo;
                umul_ppmm(hi, lo, T_p, (evals + d*(j - 1))[k]);
                add_sssaaaaaa(V_p[3*k+2], V_p[3*k+1], V_p[3*k+0],
                              V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], 0, hi, lo);

            }
        }

        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r_p, T_p, mod), master[0], mod) == 0);

        S_p = nmod_mul(S_p, r_p, mod); /* shift is one */
        if (S == 0)
            return -1;

        S_p = nmod_inv(S_p, mod);

        for (k = 0; k < d; k++)
        {
            mp_limb_t vk;
            NMOD_RED3(vk, V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], mod);
            (coeffs + d*i)[k] = nmod_mul(vk, S_p, mod);
        }
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        scratch[j] = nmod_pow_ui(monomials[j], mlength, mod);

    for (i = mlength; i < elength; i++)
    {
        _nmod_vec_zero(V_p, 3*d);
        S_p = 0;
        for (j = 0; j < mlength; j++)
        {
            scratch[j] = nmod_mul(scratch[j], monomials[j], mod);
            for (k = 0; k < d; k++)
            {
                mp_limb_t hi, lo;
                umul_ppmm(hi, lo, scratch[j], (coeffs + d*j)[k]);
                add_sssaaaaaa(V_p[3*k+2], V_p[3*k+1], V_p[3*k+0],
                              V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], 0, hi, lo);
            }
        }

        for (k = 0; k < d; k++)
        {
            mp_limb_t vk;
            NMOD_RED3(vk, V_p[3*k+2], V_p[3*k+1], V_p[3*k+0], mod);
            if (vk != (evals + d*i)[k])
            {
                success = 0;
                goto cleanup;
            }
        }
    }

    success = 1;

cleanup:

    flint_free(tmp);

    return success;
}

