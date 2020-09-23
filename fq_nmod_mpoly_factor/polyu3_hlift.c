/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


#if FLINT_WANT_ASSERT
static void fq_nmod_polyu_get_fq_nmod_polyun(
    n_polyu_t A,
    const n_polyun_t B,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, j;
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = B->terms[i].coeff->length - 1; j >= 0; j--)
        {
            if (_n_fq_is_zero(B->terms[i].coeff->coeffs + d*j, d))
                continue;
            n_polyu_fit_length(A, d*(A->length + 1));
            _n_fq_set(A->coeffs + d*A->length, B->terms[i].coeff->coeffs + d*j, d);
            A->exps[A->length] = B->terms[i].exp + j;
            A->length++;
        }
    }
}

static void fq_nmod_polyu_sort_terms(
    n_polyu_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, j;
    for (i = 1; i < A->length; i++)
    for (j = i; j > 0 && A->exps[j - 1] < A->exps[j]; j--)
    {
        ULONG_SWAP(A->exps[j - 1], A->exps[j]);
        _n_fq_swap(A->coeffs + d*(j - 1), A->coeffs + d*(j), d);
    }
    return;
}

static void fq_nmod_polyu_combine_like_terms(
    n_polyu_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong in, out;

    out = -1;

    for (in = 0; in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 && A->exps[out] == A->exps[in])
        {
            _n_fq_add(A->coeffs + d*out, A->coeffs + d*out, A->coeffs + d*in, ctx);
        }
        else
        {
            if (out < 0 || !_n_fq_is_zero(A->coeffs + d*out, d))
                out++;

            if (out != in)
            {
                A->exps[out] = A->exps[in];
                _n_fq_set(A->coeffs + d*out, A->coeffs + d*in, d);
            }
        }
    }

    if (out < 0 || !_n_fq_is_zero(A->coeffs + d*out, d))
        out++;

    A->length = out;
}

static int fq_nmod_polyu_equal(
    n_polyu_t A,
    n_polyu_t B,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    if (A->length != B->length)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!_n_fq_equal(A->coeffs + d*i, B->coeffs + d*i, d))
            return 0;
        if (A->exps[i] != B->exps[i])
            return 0;
    }
    return 1;
}


static void fq_nmod_polyu_mul(
    n_polyu_t A,
    const n_polyu_t B,
    const n_polyu_t C,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Ai, Bi, Ci;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_polyu_fit_length(A, d*B->length*C->length);
    Ai = 0;
    for (Bi = 0; Bi < B->length; Bi++)
    for (Ci = 0; Ci < C->length; Ci++)
    {
        A->exps[Ai] = B->exps[Bi] + C->exps[Ci];
        n_fq_mul(A->coeffs + d*Ai, B->coeffs + d*Bi, C->coeffs + d*Ci, ctx);
        Ai++;
    }
    A->length = Ai;
    fq_nmod_polyu_sort_terms(A, ctx);
    fq_nmod_polyu_combine_like_terms(A, ctx);
}

#endif



void fq_nmod_polyu3_interp_reduce_bpoly(
    n_bpoly_t Ap,
    const n_polyu_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    slong cur0, cur1, e0, e1, e2;
    fq_nmod_t p, t, q;

    fq_nmod_init(p, ctx);
    fq_nmod_init(t, ctx);
    fq_nmod_init(q, ctx);

    n_bpoly_zero(Ap);

    FLINT_ASSERT(A->length > 0);

    i = 0;

    cur0 = extract_exp(A->exps[i], 2, 3);
    cur1 = extract_exp(A->exps[i], 1, 3);
    e2   = extract_exp(A->exps[i], 0, 3);

    fq_nmod_pow_ui(t, alpha, e2, ctx);
    n_fq_get_fq_nmod(q, A->coeffs + d*i, ctx);
    fq_nmod_mul(t, t, q, ctx);

    for (i = 1; i < A->length; i++)
    {
        e0 = extract_exp(A->exps[i], 2, 3);
        e1 = extract_exp(A->exps[i], 1, 3);
        e2 = extract_exp(A->exps[i], 0, 3);

        FLINT_ASSERT(e0 <= cur0);
        if (e0 < cur0 || e1 < cur1)
        {
            n_bpoly_fq_set_coeff_fq_nmod(Ap, cur0, cur1, t, ctx);
            fq_nmod_zero(t, ctx);
        }
        else
        {
            FLINT_ASSERT(e0 == cur0);
            FLINT_ASSERT(e1 == cur1);
        }

        cur0 = e0;
        cur1 = e1;

        fq_nmod_pow_ui(p, alpha, e2, ctx);
        n_fq_get_fq_nmod(q, A->coeffs + d*i, ctx);
        fq_nmod_mul(p, p, q, ctx);
        fq_nmod_add(t, t, p, ctx);
    }

    n_bpoly_fq_set_coeff_fq_nmod(Ap, cur0, cur1, t, ctx);

    fq_nmod_clear(p, ctx);
    fq_nmod_clear(t, ctx);
    fq_nmod_clear(q, ctx);
}


/*
    T(x0, x1, x2) is in F[x2][x0, x1]
    A(x0, x1) are B(x0, x1) are in F[x0, x1]
    set T so that
        T(x0, x1, x2) == A(x0, x1) mod (x2 - alpha)
        T(x0, x1, x2) == B(x0, x1) mod (x2 + alpha)
*/

void fq_nmod_polyu3n_interp_lift_sm_bpoly(
    slong * lastdeg,
    n_polyun_t T,
    const n_bpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong lastlength = 0;
    slong Ti;
    slong Ai, j;

    Ti = 0;

    for (Ai = A->length - 1; Ai >= 0; Ai--)
    {
        n_poly_struct * Ac = A->coeffs + Ai;
        for (j = Ac->length - 1; j >= 0; j--)
        {
            if (_n_fq_is_zero(Ac->coeffs + d*j, d))
                continue;
            n_polyun_fit_length(T, Ti + 1);
            T->terms[Ti].exp = pack_exp3(Ai, j, 0);
            n_poly_fq_set_n_fq(T->terms[Ti].coeff, Ac->coeffs + d*j, ctx);
            Ti++;
            lastlength = 1;
        }
    }

    T->length = Ti;

    *lastdeg = lastlength - 1;

    return;
}

/*
    F is in Fq[x2][x0, x1]
    A is in Fq[x0, x1]
*/
int fq_nmod_polyu3n_interp_crt_sm_bpoly(
    slong * lastdeg,
    n_polyun_t F,
    n_polyun_t T,
    const n_bpoly_t A,
    const n_poly_t modulus,
    const mp_limb_t * alpha,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    int changed = 0;
    slong lastlength = 0;
    n_polyun_term_struct * Tterms;
    slong Ti;
    n_polyun_term_struct * Fterms = F->terms;
    slong Flen = F->length;
    slong Fi;
    const n_poly_struct * Acoeffs = A->coeffs;
    slong Ai, ai;
    n_poly_t tp;
    mp_limb_t * v = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_poly_init(tp);

    FLINT_ASSERT(n_bpoly_fq_is_canonical(A, ctx));
    FLINT_ASSERT(n_polyun_fq_is_canonical(F, ctx));

    n_polyun_fit_length(T, FLINT_MAX(Flen, A->length));
    Tterms = T->terms;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(A->coeffs + Ai);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            n_polyun_fit_length(T, Ti + extra + 1);
            Tterms = T->terms;
        }

        FLINT_ASSERT(Fi >= Flen || Fterms[Fi].coeff->length > 0);
        FLINT_ASSERT(Ai < 0 || !_n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));

        if (Fi < Flen && Ai >= 0 && Fterms[Fi].exp == pack_exp3(Ai, ai, 0))
        {
            /* F term ok, A term ok */
            n_poly_fq_evaluate_n_fq(v, Fterms[Fi].coeff, alpha, ctx);
            _n_fq_sub(v, Acoeffs[Ai].coeffs + d*ai, v, ctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
                n_poly_fq_scalar_mul_n_fq(tp, modulus, v, ctx);
                n_poly_fq_add(Tterms[Ti].coeff, Fterms[Fi].coeff, tp, ctx);
            }
            else
            {
                n_poly_fq_set(Tterms[Ti].coeff, Fterms[Fi].coeff, ctx);
            }
            Tterms[Ti].exp = Fterms[Fi].exp;

            Fi++;

            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Fi < Flen && (Ai < 0 || Fterms[Fi].exp > pack_exp3(Ai, ai, 0)))
        {
            /* F term ok, A term missing */
            n_poly_fq_evaluate_n_fq(v, Fterms[Fi].coeff, alpha, ctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
                n_poly_fq_scalar_mul_n_fq(tp, modulus, v, ctx);
                n_poly_fq_sub(Tterms[Ti].coeff, Fterms[Fi].coeff, tp, ctx);
            }
            else
            {
                n_poly_fq_set(Tterms[Ti].coeff, Fterms[Fi].coeff, ctx);
            }

            Tterms[Ti].exp = Fterms[Fi].exp;

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai >= 0 && (Fi >= Flen || Fterms[Fi].exp < pack_exp3(Ai, ai, 0)));
            /* F term missing, Aterm ok */
            changed = 1;
            n_poly_fq_scalar_mul_n_fq(Tterms[Ti].coeff, modulus, Acoeffs[Ai].coeffs + d*ai, ctx);
            Tterms[Ti].exp = pack_exp3(Ai, ai, 0);

            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }

        FLINT_ASSERT(!n_poly_is_zero(Tterms[Ti].coeff));
        lastlength = FLINT_MAX(lastlength, Tterms[Ti].coeff->length);
        Ti++;
    }
    T->length = Ti;

    if (changed)
        n_polyun_swap(T, F);

    FLINT_ASSERT(n_polyun_fq_is_canonical(F, ctx));

    n_poly_clear(tp);

    flint_free(v);

    *lastdeg = lastlength - 1;
    return changed;
}



int n_polyu3_fq_hlift(
    slong r,
    n_polyun_struct * BB,
    n_polyu_t A,
    n_polyu_struct * B,
    const fq_nmod_t beta,
    slong degree_inner, /* required degree in X (var 1) */
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    int success;
    slong i, j;
    n_polyun_t T;
    n_bpoly_struct * Bp;
    n_bpoly_t Ap;
    n_poly_t modulus;
    fq_nmod_t alpha;
    slong * BBdegZ;
    slong AdegY, AdegX, AdegZ;
    slong bad_primes_left;
    mp_limb_t * alpha_ = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    mp_limb_t * c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    fq_nmod_init(alpha, ctx);

    FLINT_ASSERT(n_polyu_fq_is_canonical(A, ctx));
    for (i = 0; i < r; i++)
        FLINT_ASSERT(n_polyu_fq_is_canonical(B + i, ctx));

    BBdegZ = FLINT_ARRAY_ALLOC(r, slong);
    Bp = FLINT_ARRAY_ALLOC(r, n_bpoly_struct);
    for (i = 0; i < r; i++)
        n_bpoly_init(Bp + i);

    n_polyun_init(T);
    n_poly_init(modulus);
    n_bpoly_init(Ap);

    n_polyu3_degrees(&AdegY, &AdegX, &AdegZ, A);
    if (AdegX != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    _n_poly_fq_one(modulus, d);

    fq_nmod_zero(alpha, ctx);

    bad_primes_left = FLINT_MAX(5, AdegZ);

choose_prime:

    if (fq_nmod_next(alpha, ctx) == 0)
    {
        success = -1;
        goto cleanup;
    }

    n_fq_set_fq_nmod(alpha_, alpha, ctx);

    fq_nmod_polyu3_interp_reduce_bpoly(Ap, A, alpha, ctx);
    for (i = 0; i < r; i++)
        fq_nmod_polyu3_interp_reduce_bpoly(Bp + i, B + i, alpha, ctx);

    if (r < 3)
        success = n_bpoly_fq_hlift2(Ap, Bp + 0, Bp + 1, beta, degree_inner, ctx);
    else
        success = n_bpoly_fq_hlift(r, Ap, Bp, beta, degree_inner, ctx);

    if (success < 1)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    if (n_poly_degree(modulus) > 0)
    {
        n_poly_fq_evaluate_n_fq(c, modulus, alpha_, ctx);
        n_fq_inv(c, c, ctx);
        n_poly_fq_scalar_mul_n_fq(modulus, modulus, c, ctx);

        for (i = 0; i < r; i++)
        {
            fq_nmod_polyu3n_interp_crt_sm_bpoly(BBdegZ + i, BB + i, T,
                                                 Bp + i, modulus, alpha_, ctx);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            fq_nmod_polyu3n_interp_lift_sm_bpoly(BBdegZ + i, BB + i, Bp + i, ctx);
        }
    }

    n_poly_fq_shift_left_scalar_submul(modulus, 1, alpha_, ctx);

    j = BBdegZ[0];
    for (i = 1; i < r; i++)
        j += BBdegZ[i];

    if (j > AdegZ)
    {
        success = 0;
        goto cleanup;
    }

    if (n_poly_degree(modulus) <= AdegZ)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success == 1)
    {
        n_polyu_t T1, T2, T3;
        n_polyu_init(T1);
        n_polyu_init(T2);
        n_polyu_init(T3);
        fq_nmod_polyu_get_fq_nmod_polyun(T2, BB + 0, ctx);
        fq_nmod_polyu_get_fq_nmod_polyun(T3, BB + 1, ctx);
        fq_nmod_polyu_mul(T1, T2, T3, ctx);
        for (i = 2; i < r; i++)
        {
            fq_nmod_polyu_get_fq_nmod_polyun(T3, BB + i, ctx);
            fq_nmod_polyu_mul(T2, T1, T3, ctx);
            n_polyu_swap(T2, T1);
        }
        FLINT_ASSERT(fq_nmod_polyu_equal(A, T1, ctx));
        n_polyu_clear(T1);
        n_polyu_clear(T2);
        n_polyu_clear(T3);
    }
#endif

    n_polyun_clear(T);
    n_bpoly_clear(Ap);

    for (i = 0; i < r; i++)
        n_bpoly_clear(Bp + i);
    flint_free(BBdegZ);
    flint_free(Bp);

    n_poly_clear(modulus);
    fq_nmod_clear(alpha, ctx);
    flint_free(alpha_);
    flint_free(c);

    return success;
}
