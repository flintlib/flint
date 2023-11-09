/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"


#ifdef FLINT_WANT_ASSERT
static void fq_zech_polyu_set_fq_zech_polyun(
    fq_zech_polyu_t A,
    const fq_zech_polyun_t B,
    const fq_zech_ctx_t ctx)
{
    slong i, j;
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = B->coeffs[i].length - 1; j >= 0; j--)
        {
            if (fq_zech_is_zero(B->coeffs[i].coeffs + j, ctx))
                continue;
            fq_zech_polyu_fit_length(A, A->length + 1, ctx);
            fq_zech_set(A->coeffs + A->length, B->coeffs[i].coeffs + j, ctx);
            A->exps[A->length] = B->exps[i] + j;
            A->length++;
        }
    }
}

static void fq_zech_polyu_sort_terms(
    fq_zech_polyu_t A,
    const fq_zech_ctx_t ctx)
{
    slong i, j;
    for (i = 1; i < A->length; i++)
    for (j = i; j > 0 && A->exps[j - 1] < A->exps[j]; j--)
    {
        FLINT_SWAP(ulong, A->exps[j - 1], A->exps[j]);
        fq_zech_swap(A->coeffs + j - 1, A->coeffs + j, ctx);
    }
    return;
}

static void fq_zech_polyu_combine_like_terms(
    fq_zech_polyu_t A,
    const fq_zech_ctx_t ctx)
{
    slong in, out;

    out = -1;

    for (in = 0; in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 && A->exps[out] == A->exps[in])
        {
            fq_zech_add(A->coeffs + out, A->coeffs + out, A->coeffs + in, ctx);
        }
        else
        {
            if (out < 0 || !fq_zech_is_zero(A->coeffs + out, ctx))
                out++;

            if (out != in)
            {
                A->exps[out] = A->exps[in];
                fq_zech_set(A->coeffs + out, A->coeffs + in, ctx);
            }
        }
    }

    if (out < 0 || !fq_zech_is_zero(A->coeffs + out, ctx))
        out++;

    A->length = out;
}

static int fq_zech_polyu_equal(
    fq_zech_polyu_t A,
    fq_zech_polyu_t B,
    const fq_zech_ctx_t ctx)
{
    slong i;
    if (A->length != B->length)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!fq_zech_equal(A->coeffs + i, B->coeffs + i, ctx))
            return 0;
        if (A->exps[i] != B->exps[i])
            return 0;
    }
    return 1;
}


static void fq_zech_polyu_mul(
    fq_zech_polyu_t A,
    const fq_zech_polyu_t B,
    const fq_zech_polyu_t C,
    const fq_zech_ctx_t ctx)
{
    slong Ai, Bi, Ci;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fq_zech_polyu_fit_length(A, B->length*C->length, ctx);
    Ai = 0;
    for (Bi = 0; Bi < B->length; Bi++)
    for (Ci = 0; Ci < C->length; Ci++)
    {
        A->exps[Ai] = B->exps[Bi] + C->exps[Ci];
        fq_zech_mul(A->coeffs + Ai, B->coeffs + Bi, C->coeffs + Ci, ctx);
        Ai++;
    }
    A->length = Ai;
    fq_zech_polyu_sort_terms(A, ctx);
    fq_zech_polyu_combine_like_terms(A, ctx);
}

#endif



void fq_zech_polyu3_interp_reduce_bpoly(
    fq_zech_bpoly_t Ap,
    const fq_zech_polyu_t A,
    const fq_zech_t alpha,
    const fq_zech_ctx_t ctx)
{
    slong i;
    slong cur0, cur1, e0, e1, e2;
    fq_zech_t p, t, q;

    fq_zech_init(p, ctx);
    fq_zech_init(t, ctx);
    fq_zech_init(q, ctx);

    fq_zech_bpoly_zero(Ap, ctx);

    FLINT_ASSERT(A->length > 0);

    i = 0;

    cur0 = extract_exp(A->exps[i], 2, 3);
    cur1 = extract_exp(A->exps[i], 1, 3);
    e2   = extract_exp(A->exps[i], 0, 3);

    fq_zech_pow_ui(t, alpha, e2, ctx);
    fq_zech_set(q, A->coeffs + i, ctx);
    fq_zech_mul(t, t, q, ctx);

    for (i = 1; i < A->length; i++)
    {
        e0 = extract_exp(A->exps[i], 2, 3);
        e1 = extract_exp(A->exps[i], 1, 3);
        e2 = extract_exp(A->exps[i], 0, 3);

        FLINT_ASSERT(e0 <= cur0);
        if (e0 < cur0 || e1 < cur1)
        {
            fq_zech_bpoly_set_coeff_fq_zech(Ap, cur0, cur1, t, ctx);
            fq_zech_zero(t, ctx);
        }
        else
        {
            FLINT_ASSERT(e0 == cur0);
            FLINT_ASSERT(e1 == cur1);
        }

        cur0 = e0;
        cur1 = e1;

        fq_zech_pow_ui(p, alpha, e2, ctx);
        fq_zech_set(q, A->coeffs + i, ctx);
        fq_zech_mul(p, p, q, ctx);
        fq_zech_add(t, t, p, ctx);
    }

    fq_zech_bpoly_set_coeff_fq_zech(Ap, cur0, cur1, t, ctx);

    fq_zech_clear(p, ctx);
    fq_zech_clear(t, ctx);
    fq_zech_clear(q, ctx);
}


/*
    T(x0, x1, x2) is in F[x2][x0, x1]
    A(x0, x1) are B(x0, x1) are in F[x0, x1]
    set T so that
        T(x0, x1, x2) == A(x0, x1) mod (x2 - alpha)
        T(x0, x1, x2) == B(x0, x1) mod (x2 + alpha)
*/

void fq_zech_polyu3n_interp_lift_sm_bpoly(
    slong * lastdeg,
    fq_zech_polyun_t T,
    const fq_zech_bpoly_t A,
    const fq_zech_ctx_t ctx)
{
    slong lastlength = 0;
    slong Ti;
    slong Ai, j;

    Ti = 0;

    for (Ai = A->length - 1; Ai >= 0; Ai--)
    {
        fq_zech_poly_struct * Ac = A->coeffs + Ai;
        for (j = Ac->length - 1; j >= 0; j--)
        {
            if (fq_zech_is_zero(Ac->coeffs + j, ctx))
                continue;
            fq_zech_polyun_fit_length(T, Ti + 1, ctx);
            T->exps[Ti] = pack_exp3(Ai, j, 0);
            fq_zech_poly_set_fq_zech(T->coeffs + Ti, Ac->coeffs + j, ctx);
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
int fq_zech_polyu3n_interp_crt_sm_bpoly(
    slong * lastdeg,
    fq_zech_polyun_t F,
    fq_zech_polyun_t T,
    const fq_zech_bpoly_t A,
    const fq_zech_poly_t modulus,
    const fq_zech_t alpha,
    const fq_zech_ctx_t ctx)
{
    int changed = 0;
    slong lastlength = 0;
    fq_zech_poly_struct * Tcoeffs;
    ulong * Texps;
    slong Ti;
    fq_zech_poly_struct * Fcoeffs = F->coeffs;
    ulong * Fexps = F->exps;
    slong Flen = F->length;
    slong Fi;
    const fq_zech_poly_struct * Acoeffs = A->coeffs;
    slong Ai, ai;
    fq_zech_poly_t tp;
    fq_zech_t v;

    fq_zech_init(v, ctx);
    fq_zech_poly_init(tp, ctx);

    FLINT_ASSERT(fq_zech_bpoly_is_canonical(A, ctx));
    FLINT_ASSERT(fq_zech_polyun_is_canonical(F, ctx));

    fq_zech_polyun_fit_length(T, FLINT_MAX(Flen, A->length), ctx);
    Tcoeffs = T->coeffs;
    Texps = T->exps;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : fq_zech_poly_degree(A->coeffs + Ai, ctx);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            fq_zech_polyun_fit_length(T, Ti + extra + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        FLINT_ASSERT(Fi >= Flen || Fcoeffs[Fi].length > 0);
        FLINT_ASSERT(Ai < 0 || !fq_zech_is_zero(Acoeffs[Ai].coeffs + ai, ctx));

        if (Fi < Flen && Ai >= 0 && Fexps[Fi] == pack_exp3(Ai, ai, 0))
        {
            /* F term ok, A term ok */
            fq_zech_poly_evaluate_fq_zech(v, Fcoeffs + Fi, alpha, ctx);
            fq_zech_sub(v, Acoeffs[Ai].coeffs + ai, v, ctx);
            if (!fq_zech_is_zero(v, ctx))
            {
                changed = 1;
                fq_zech_poly_scalar_mul_fq_zech(tp, modulus, v, ctx);
                fq_zech_poly_add(Tcoeffs + Ti, Fcoeffs + Fi, tp, ctx);
            }
            else
            {
                fq_zech_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx);
            }
            Texps[Ti] = Fexps[Fi];

            Fi++;

            do {
                ai--;
            } while (ai >= 0 && fq_zech_is_zero(Acoeffs[Ai].coeffs + ai, ctx));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = fq_zech_poly_degree(Acoeffs + Ai, ctx);
            }
        }
        else if (Fi < Flen && (Ai < 0 || Fexps[Fi] > pack_exp3(Ai, ai, 0)))
        {
            /* F term ok, A term missing */
            fq_zech_poly_evaluate_fq_zech(v, Fcoeffs + Fi, alpha, ctx);
            if (!fq_zech_is_zero(v, ctx))
            {
                changed = 1;
                fq_zech_poly_scalar_mul_fq_zech(tp, modulus, v, ctx);
                fq_zech_poly_sub(Tcoeffs + Ti, Fcoeffs + Fi, tp, ctx);
            }
            else
            {
                fq_zech_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx);
            }

            Texps[Ti] = Fexps[Fi];

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai >= 0 && (Fi >= Flen || Fexps[Fi] < pack_exp3(Ai, ai, 0)));
            /* F term missing, Aterm ok */
            changed = 1;
            fq_zech_poly_scalar_mul_fq_zech(Tcoeffs + Ti, modulus, Acoeffs[Ai].coeffs + ai, ctx);
            Texps[Ti] = pack_exp3(Ai, ai, 0);

            do {
                ai--;
            } while (ai >= 0 && fq_zech_is_zero(Acoeffs[Ai].coeffs + ai, ctx));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = fq_zech_poly_degree(Acoeffs + Ai, ctx);
            }
        }

        FLINT_ASSERT(!fq_zech_poly_is_zero(Tcoeffs + Ti, ctx));
        lastlength = FLINT_MAX(lastlength, Tcoeffs[Ti].length);
        Ti++;
    }
    T->length = Ti;

    if (changed)
        fq_zech_polyun_swap(T, F, ctx);

    FLINT_ASSERT(fq_zech_polyun_is_canonical(F, ctx));

    fq_zech_poly_clear(tp, ctx);
    fq_zech_clear(v, ctx);

    *lastdeg = lastlength - 1;
    return changed;
}

int fq_zech_next(
    fq_zech_t a,
    const fq_zech_ctx_t ctx);

void fq_zech_poly_shift_left_scalar_submul(
    fq_zech_poly_t modulus,
    slong k,
    const fq_zech_t alpha,
    const fq_zech_ctx_t ctx)
{
    fq_zech_poly_t t;
    fq_zech_poly_init(t, ctx);
    fq_zech_poly_scalar_mul_fq_zech(t, modulus, alpha, ctx);
    fq_zech_poly_shift_left(modulus, modulus, k, ctx);
    fq_zech_poly_sub(modulus, modulus, t, ctx);
    fq_zech_poly_clear(t, ctx);
}


int fq_zech_polyu3_hlift(
    slong r,
    fq_zech_polyun_struct * BB,
    fq_zech_polyu_t A,
    fq_zech_polyu_struct * B,
    const fq_zech_t beta,
    slong degree_inner, /* required degree in X (var 1) */
    const fq_zech_ctx_t ctx)
{
    int success;
    slong i, j;
    fq_zech_polyun_t T;
    fq_zech_bpoly_struct * Bp;
    fq_zech_bpoly_t Ap;
    fq_zech_poly_t modulus;
    fq_zech_t alpha;
    slong * BBdegZ;
    slong AdegY, AdegX, AdegZ;
    slong bad_primes_left;
    fq_zech_t c;

    fq_zech_init(c, ctx);
    fq_zech_init(alpha, ctx);

    FLINT_ASSERT(fq_zech_polyu_is_canonical(A, ctx));
    for (i = 0; i < r; i++)
        FLINT_ASSERT(fq_zech_polyu_is_canonical(B + i, ctx));

    BBdegZ = FLINT_ARRAY_ALLOC(r, slong);
    Bp = FLINT_ARRAY_ALLOC(r, fq_zech_bpoly_struct);
    for (i = 0; i < r; i++)
        fq_zech_bpoly_init(Bp + i, ctx);

    fq_zech_polyun_init(T, ctx);
    fq_zech_poly_init(modulus, ctx);
    fq_zech_bpoly_init(Ap, ctx);

    fq_zech_polyu3_degrees(&AdegY, &AdegX, &AdegZ, A);
    if (AdegX != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    fq_zech_poly_one(modulus, ctx);

    fq_zech_zero(alpha, ctx);

    bad_primes_left = FLINT_MAX(5, AdegZ);

choose_prime:

    if (fq_zech_next(alpha, ctx) == 0)
    {
        success = -1;
        goto cleanup;
    }

    fq_zech_polyu3_interp_reduce_bpoly(Ap, A, alpha, ctx);
    for (i = 0; i < r; i++)
        fq_zech_polyu3_interp_reduce_bpoly(Bp + i, B + i, alpha, ctx);

    if (r < 3)
        success = fq_zech_bpoly_hlift2(Ap, Bp + 0, Bp + 1, beta, degree_inner, ctx);
    else
        success = fq_zech_bpoly_hlift(r, Ap, Bp, beta, degree_inner, ctx);

    if (success < 1)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    if (fq_zech_poly_degree(modulus, ctx) > 0)
    {
        fq_zech_poly_evaluate_fq_zech(c, modulus, alpha, ctx);
        fq_zech_inv(c, c, ctx);
        fq_zech_poly_scalar_mul_fq_zech(modulus, modulus, c, ctx);

        for (i = 0; i < r; i++)
        {
            fq_zech_polyu3n_interp_crt_sm_bpoly(BBdegZ + i, BB + i, T,
                                                  Bp + i, modulus, alpha, ctx);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            fq_zech_polyu3n_interp_lift_sm_bpoly(BBdegZ + i, BB + i, Bp + i, ctx);
        }
    }

    fq_zech_poly_shift_left_scalar_submul(modulus, 1, alpha, ctx);

    j = BBdegZ[0];
    for (i = 1; i < r; i++)
        j += BBdegZ[i];

    if (j > AdegZ)
    {
        success = 0;
        goto cleanup;
    }

    if (fq_zech_poly_degree(modulus, ctx) <= AdegZ)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

#ifdef FLINT_WANT_ASSERT
    if (success == 1)
    {
        fq_zech_polyu_t T1, T2, T3;
        fq_zech_polyu_init(T1, ctx);
        fq_zech_polyu_init(T2, ctx);
        fq_zech_polyu_init(T3, ctx);
        fq_zech_polyu_set_fq_zech_polyun(T2, BB + 0, ctx);
        fq_zech_polyu_set_fq_zech_polyun(T3, BB + 1, ctx);
        fq_zech_polyu_mul(T1, T2, T3, ctx);
        for (i = 2; i < r; i++)
        {
            fq_zech_polyu_set_fq_zech_polyun(T3, BB + i, ctx);
            fq_zech_polyu_mul(T2, T1, T3, ctx);
            fq_zech_polyu_swap(T2, T1, ctx);
        }
        FLINT_ASSERT(fq_zech_polyu_equal(A, T1, ctx));
        fq_zech_polyu_clear(T1, ctx);
        fq_zech_polyu_clear(T2, ctx);
        fq_zech_polyu_clear(T3, ctx);
    }
#endif

    fq_zech_polyun_clear(T, ctx);
    fq_zech_bpoly_clear(Ap, ctx);

    for (i = 0; i < r; i++)
        fq_zech_bpoly_clear(Bp + i, ctx);
    flint_free(BBdegZ);
    flint_free(Bp);

    fq_zech_poly_clear(modulus, ctx);
    fq_zech_clear(alpha, ctx);
    fq_zech_clear(c, ctx);

    return success;
}
