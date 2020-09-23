/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"


int fmpz_mpolyu_is_canonical(const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (A->length > A->alloc)
    {
        return 0;
    }

    for (i = 0; i < A->length; i++)
    {
        if (   fmpz_mpoly_is_zero(A->coeffs + i, ctx)
            || !fmpz_mpoly_is_canonical(A->coeffs + i, ctx))
        {
            return 0;
        }

        if (i > 0 && A->exps[i - 1] <= A->exps[i])
        {
            return 0;
        }
    }

    return 1;
}

void fmpz_mpolyu_init(fmpz_mpolyu_t A, flint_bitcnt_t bits,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void fmpz_mpolyu_clear(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void fmpz_mpolyu_swap(fmpz_mpolyu_t A, fmpz_mpolyu_t B,
                                                   const fmpz_mpoly_ctx_t uctx)
{
   fmpz_mpolyu_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mpolyu_zero(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fmpz_mpoly_clear(A->coeffs + i, uctx);
        fmpz_mpoly_init(A->coeffs + i, uctx);
    }
    A->length = 0;
}



void fmpz_mpolyu_print_pretty(
    const fmpz_mpolyu_t poly,
    const char ** x,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fmpz_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")*X^%wd", poly->exps[i]);
    }
}

void fmpz_mpolyuu_print_pretty(
    const fmpz_mpolyu_t poly,
    const char ** x,
    slong nmainvars,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/nmainvars);

    if (poly->length == 0)
        flint_printf("0");

    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fmpz_mpoly_print_pretty(poly->coeffs + i, x, ctx);
        flint_printf(")");
        for (j = nmainvars - 1; j >= 0; j--)
        {
            flint_printf("*X%wd^%wd", nmainvars - 1 - j,
                    mask & (poly->exps[i] >> (FLINT_BITS/nmainvars*j)));
        }
    }
}


void fmpz_mpolyu_fit_length(fmpz_mpolyu_t A, slong length,
                                                   const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mpoly_struct *) flint_malloc(
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        } else
        {
            A->exps = (ulong *) flint_realloc(A->exps,
                                                      new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(fmpz_mpoly_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpoly_init(A->coeffs + i, uctx);
            fmpz_mpoly_fit_bits(A->coeffs + i, A->bits, uctx);
            (A->coeffs + i)->bits = A->bits;
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mpolyu_one(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx)
{
    fmpz_mpolyu_fit_length(A, WORD(1), uctx);
    A->exps[0] = UWORD(0);
    fmpz_mpoly_one(A->coeffs + 0, uctx);
    A->length = WORD(1);
}


int fmpz_equal_upto_unit(
    const fmpz_t a,
    const fmpz_t b)
{
    if (fmpz_equal(a, b))
        return 1;

    if (fmpz_cmpabs(a, b) == 0)
        return -1;

    return 0;
}

int fmpz_mpoly_equal_upto_unit(
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    int res;
    slong i, n = A->length;

    if (A->length != B->length)
        return 0;

    if (A->length < 1)
        return 1;

    if (mpoly_monomials_cmp(A->exps, A->bits, B->exps, B->bits, n, ctx->minfo) != 0)
        return 0;

    i = 0;

    res = fmpz_equal_upto_unit(A->coeffs + i, B->coeffs + i);
    if (res == 0)
        return 0;

    for (i++; i < n; i++)
    {
        int res2 = fmpz_equal_upto_unit(A->coeffs + i, B->coeffs + i);
        if (res2 == 0 || res != res2)
            return 0;
    }

    return res;
}


int fmpz_mpolyu_equal_upto_unit(const fmpz_mpolyu_t A, const fmpz_mpolyu_t B,
                                                   const fmpz_mpoly_ctx_t ctx)
{
    int res;
    slong i;

    if (A->length != B->length)
        return 0;

    if (A->length < 1)
        return 1;

    for (i = 0; i < A->length; i++)
    {
        if (A->exps[i] != B->exps[i])
            return 0;
    }

    i = 0;

    res = fmpz_mpoly_equal_upto_unit(A->coeffs + i, B->coeffs + i, ctx);
    if (res == 0)
        return 0;

    for (i++; i < A->length; i++)
    {
        int res2 = fmpz_mpoly_equal_upto_unit(A->coeffs + i, B->coeffs + i, ctx);
        if (res2 == 0 || res != res2)
            return 0;
    }

    return res;
}

void fmpz_mpolyu_inner_degrees_si(
    slong * degs,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, * t;
    TMP_INIT;

    if (A->length < 1)
    {
        for (j = 0; j < ctx->minfo->nvars; j++)
            degs[j] = -1;
        return;
    }

    TMP_START;

    t = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    
    fmpz_mpoly_degrees_si(degs, A->coeffs + 0, ctx);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_degrees_si(t, A->coeffs + i, ctx);
        for (j = 0; j < ctx->minfo->nvars; j++)
            degs[j] = FLINT_MAX(degs[j], t[j]);
    }

    TMP_END;
}

void fmpz_mpolyu_set(fmpz_mpolyu_t A, const fmpz_mpolyu_t B,
                                                   const fmpz_mpoly_ctx_t uctx)
{
    slong i;
    fmpz_mpoly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Alen, Blen;

    Alen = 0;
    Blen = B->length;
    fmpz_mpolyu_fit_length(A, Blen, uctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        fmpz_mpoly_set(Acoeff + Alen, Bcoeff + i, uctx);
        Aexp[Alen++] = Bexp[i];
    }
    Alen = Blen;

    /* demote remaining coefficients */
    for (i = Alen; i < A->length; i++)
    {
        fmpz_mpoly_clear(Acoeff + i, uctx);
        fmpz_mpoly_init(Acoeff + i, uctx);
    }
    A->length = Alen;
}


/* if the coefficient doesn't exist, a new one is created (and set to zero) */
fmpz_mpoly_struct * _fmpz_mpolyu_get_coeff(
    fmpz_mpolyu_t A,
    ulong pow,
    const fmpz_mpoly_ctx_t uctx)
{
    slong i, j, a, b;
    fmpz_mpoly_struct * xk;

    a = 0;
    b = A->length;

    if (b == 0 || pow > A->exps[0])
    {
        i = 0;
        goto create_new;
    }

    if (pow == A->exps[b - 1])
    {
        return A->coeffs + b - 1;
    }

try_again:

    if (b - a < 8)
    {
        for (i = a; i < b && (A->exps[i] >= pow); i++)
        {
            if (A->exps[i] == pow)
            {
                return A->coeffs + i;
            }
        }
        goto create_new;
    }

    i = a + (b - a)/2;
    if (A->exps[i] == pow)
    {
        return A->coeffs + i;
    }
    else if (A->exps[i] > pow)
    {
        a = i;
    }
    else
    {
        b = i;
    }
    goto try_again;

create_new: /* new at position i */

    fmpz_mpolyu_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        fmpz_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
    }

    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);

    return xk;
}


/*
    Convert B to A using the variable permutation perm.
    The uctx (m vars) should be the context of A.
    The ctx (n vars) should be the context of B.

    operation on each term:

    for 0 <= k < m
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]
*/
void fmpz_mpoly_to_mpoly_perm_deflate_threaded_pool(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t uctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);

    TMP_START;

    uexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    NA = mpoly_words_per_exp(A->bits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    fmpz_mpoly_fit_length(A, B->length, uctx);
    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            if (stride[l] > 1)
            {
                FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == UWORD(0));
                uexps[k] = (Bexps[l] - shift[l]) / stride[l];
            }
            else
            {
                uexps[k] = Bexps[l] - shift[l];
            }
        }

        fmpz_set(A->coeffs + j, B->coeffs + j);
        mpoly_set_monomial_ui(A->exps + NA*j, uexps, A->bits, uctx->minfo);
    }
    A->length = B->length;

    fmpz_mpoly_sort_terms(A, uctx);

    TMP_END;
}






typedef struct
{
    volatile slong index;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong length;
    fmpz_mpoly_struct * coeffs;
    const fmpz_mpoly_ctx_struct * ctx;
}
_sort_arg_struct;

typedef _sort_arg_struct _sort_arg_t[1];


static void _worker_sort(void * varg)
{
    _sort_arg_struct * arg = (_sort_arg_struct *) varg;
    slong i;

get_next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&arg->mutex);
#endif
    i = arg->index;
    arg->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&arg->mutex);
#endif

    if (i >= arg->length)
        goto cleanup;

    fmpz_mpoly_sort_terms(arg->coeffs + i, arg->ctx);

    goto get_next_index;

cleanup:

    return;
}



typedef struct
{
    fmpz_mpoly_t poly;
    slong threadidx;
}
_arrayconvertu_base_elem_struct;


typedef struct
{
    const fmpz_mpoly_ctx_struct * ctx, * uctx;
    slong degbx;
    const slong * perm;
    const ulong * shift, * stride;
    flint_bitcnt_t Abits;
    const fmpz_mpoly_struct * B;
    _arrayconvertu_base_elem_struct * array;
    slong nthreads;
}
_arrayconvertu_base_struct;

typedef _arrayconvertu_base_struct _arrayconvertu_base_t[1];


typedef struct
{
    slong idx;
    _arrayconvertu_base_struct * base;
}
_arrayconvertu_worker_arg_struct;


void _arrayconvertu_worker(void * varg)
{
    _arrayconvertu_worker_arg_struct * arg = (_arrayconvertu_worker_arg_struct *) varg;
    _arrayconvertu_base_struct * base = arg->base;
    const fmpz_mpoly_ctx_struct * uctx = base->uctx;
    const fmpz_mpoly_ctx_struct * ctx = base->ctx;
    const slong * perm = base->perm;
    const ulong * shift = base->shift;
    const ulong * stride = base->stride;
    const ulong shiftx  = shift[perm[0]];
    const ulong stridex = stride[perm[0]];
    const fmpz_mpoly_struct * B = base->B;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    slong xoffset, xshift;
    slong j, k, l, arrayidx;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fmpz_mpoly_struct * Ac;
    TMP_INIT;

    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(base->Abits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    mpoly_gen_offset_shift_sp(&xoffset, &xshift, perm[0], B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        ulong xexp = ((B->exps[NB*j + xoffset] >> xshift) & mask) - shiftx;
        if (stridex != 1)
        {
            FLINT_ASSERT((xexp % stridex) == 0);
            xexp /= stridex;
        }
        arrayidx = xexp;

        FLINT_ASSERT(arrayidx < base->degbx);
        if (base->array[arrayidx].threadidx == arg->idx)
        {
            mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
            for (k = 0; k < m + 1; k++)
            {
                l = perm[k];
                if (stride[l] == 1)
                {
                    uexps[k] = (Bexps[l] - shift[l]);
                }
                else
                {
                    FLINT_ASSERT(stride[l] != 0);
                    FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == 0);
                    uexps[k] = (Bexps[l] - shift[l]) / stride[l];
                }
            }
            Ac = base->array[arrayidx].poly;
            fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
            fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
            mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 1, base->Abits, uctx->minfo);
            Ac->length++;
        }
    }

    /* sort all of ours */
    for (arrayidx = base->degbx - 1; arrayidx >= 0; arrayidx--)
    {
        if (base->array[arrayidx].threadidx == arg->idx)
        {
            fmpz_mpoly_sort_terms(base->array[arrayidx].poly, uctx);
        }
    }

    TMP_END;
}



/*
    Convert B to A using the variable permutation perm.
    The uctx (m vars) should be the context of the coefficients of A.
    The ctx (n vars) should be the context of B.

    operation on each term:

    for 0 <= k < m + 1
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]

    the most significant main variable uses k = 0
    the coefficients of A use variables k = 1 ... m
    maxexps if it exists is supposed to be a degree bound on B
*/
void fmpz_mpoly_to_mpolyu_perm_deflate_threaded_pool(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t uctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const ulong * maxexps, /* nullptr is ok */
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong limit = 1000;
    slong degbx;
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fmpz_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);

    fmpz_mpolyu_zero(A, uctx);

    /* strict degree bounds on the result */
    degbx = limit + 1;
    if (maxexps != NULL)
    {
        degbx = (maxexps[perm[0]] - shift[perm[0]])/stride[perm[0]] + 1;
    }

    if (degbx <= limit)
    {
        _arrayconvertu_base_t base;
        _arrayconvertu_worker_arg_struct * args;

        base->ctx = ctx;
        base->uctx = uctx;
        base->degbx = degbx;
        base->perm = perm;
        base->shift = shift;
        base->stride = stride;
        base->Abits = A->bits;
        base->B = B;
        base->nthreads = num_handles + 1;
        base->array = (_arrayconvertu_base_elem_struct *) flint_malloc(
                                degbx*sizeof(_arrayconvertu_base_elem_struct));
        for (i = degbx - 1; i >= 0; i--)
        {
            base->array[i].threadidx = i % base->nthreads;
            fmpz_mpoly_init3(base->array[i].poly, 0, A->bits, uctx);
        }
        args = (_arrayconvertu_worker_arg_struct *) flint_malloc(
                     base->nthreads*sizeof(_arrayconvertu_worker_arg_struct));

        for (i = 0; i < num_handles; i++)
        {
            args[i].idx = i;
            args[i].base = base;
            thread_pool_wake(global_thread_pool, handles[i], 0,
                                             _arrayconvertu_worker, &args[i]);
        }
        i = num_handles;
        args[i].idx = i;
        args[i].base = base;
        _arrayconvertu_worker(&args[i]);
        for (i = 0; i < num_handles; i++)
        {
            thread_pool_wait(global_thread_pool, handles[i]);
        }

        A->length = 0;
        for (i = degbx - 1; i >= 0; i--)
        {
            if (base->array[i].poly->length > 0)
            {
                fmpz_mpolyu_fit_length(A, A->length + 1, uctx);
                A->exps[A->length] = i;
                fmpz_mpoly_swap(A->coeffs + A->length, base->array[i].poly, uctx);
                A->length++;
            }
            fmpz_mpoly_clear(base->array[i].poly, uctx);
        }

        flint_free(base->array);
        flint_free(args);
    }
    else
    {
        TMP_START;

        uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
        Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
        NA = mpoly_words_per_exp(A->bits, uctx->minfo);
        NB = mpoly_words_per_exp(B->bits, ctx->minfo);

        for (j = 0; j < B->length; j++)
        {
            mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
            for (k = 0; k < m + 1; k++)
            {
                l = perm[k];
                FLINT_ASSERT(stride[l] != UWORD(0));
                if (stride[l] > 1)
                {
                    FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == UWORD(0));
                    uexps[k] = (Bexps[l] - shift[l]) / stride[l];
                }
                else
                {
                    uexps[k] = Bexps[l] - shift[l];
                }
            }
            Ac = _fmpz_mpolyu_get_coeff(A, uexps[0], uctx);
            FLINT_ASSERT(Ac->bits == A->bits);

            fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
            fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
            mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 1, A->bits, uctx->minfo);
            Ac->length++;
        }

        if (num_handles > 0)
        {
            _sort_arg_t arg;

#if FLINT_USES_PTHREAD
            pthread_mutex_init(&arg->mutex, NULL);
#endif
	    arg->index = 0;
            arg->coeffs = A->coeffs;
            arg->length = A->length;
            arg->ctx = uctx;

            for (i = 0; i < num_handles; i++)
            {
                thread_pool_wake(global_thread_pool, handles[i], 0, _worker_sort, arg);
            }
            _worker_sort(arg);
            for (i = 0; i < num_handles; i++)
            {
                thread_pool_wait(global_thread_pool, handles[i]);
            }

#if FLINT_USES_PTHREAD
            pthread_mutex_destroy(&arg->mutex);
#endif
	}
        else
        {
            for (i = 0; i < A->length; i++)
            {
                fmpz_mpoly_sort_terms(A->coeffs + i, uctx);
            }
        }

        TMP_END;
    }
}



/*
    Convert B to A using the variable permutation vector perm.
    This function inverts fmpz_mpoly_to_mpoly_perm_deflate.
    A will be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/

void fmpz_mpoly_from_mpoly_perm_inflate(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t uctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong j, k, l;
    slong NA, NB;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * uexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp_sp(B->bits, uctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, B->length, NA);
    for (j = 0; j < B->length; j++)
    {
        fmpz_set(Acoeff + j, B->coeffs + j);
        mpoly_get_monomial_ui(uexps, B->exps + NB*j, B->bits, uctx->minfo);
        for (l = 0; l < n; l++)
        {
            Aexps[l] = shift[l];
        }
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            Aexps[l] += stride[l]*uexps[k];
        }
        mpoly_set_monomial_ui(Aexp + NA*j, Aexps, Abits, ctx->minfo);
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, B->length, ctx);

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}



/*
    Convert B to A using the variable permutation vector perm.
    This function inverts fmpz_mpoly_to_mpolyu_perm_deflate.
    A will be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m + 1
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fmpz_mpoly_from_mpolyu_perm_inflate(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_ctx_t uctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, k, l;
    slong NA, NB;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * uexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA);
        FLINT_ASSERT(Bc->bits == B->bits);

        for (j = 0; j < Bc->length; j++)
        {
            fmpz_set(Acoeff + Alen + j, Bc->coeffs + j);
            mpoly_get_monomial_ui(uexps + 1, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[0] = B->exps[i];
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k < m + 1; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }
            mpoly_set_monomial_ui(Aexp + NA*(Alen + j), Aexps, Abits, ctx->minfo);
        }
        Alen += Bc->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}


typedef struct
{
    fmpz_mpoly_t poly;
    slong threadidx;
}
_arrayconvertuu_base_elem_struct;


typedef struct
{
    const fmpz_mpoly_ctx_struct * ctx, * uctx;
    slong degbx, degby;
    const slong * perm;
    const ulong * shift, * stride;
    flint_bitcnt_t Abits;
    const fmpz_mpoly_struct * B;
    _arrayconvertuu_base_elem_struct * array;
    slong nthreads;
}
_arrayconvertuu_base_struct;

typedef _arrayconvertuu_base_struct _arrayconvertuu_base_t[1];


typedef struct
{
    slong idx;
    _arrayconvertuu_base_struct * base;
}
_arrayconvertuu_worker_arg_struct;


void _arrayconvertuu_worker(void * varg)
{
    _arrayconvertuu_worker_arg_struct * arg = (_arrayconvertuu_worker_arg_struct *) varg;
    _arrayconvertuu_base_struct * base = arg->base;
    const fmpz_mpoly_ctx_struct * uctx = base->uctx;
    const fmpz_mpoly_ctx_struct * ctx = base->ctx;
    const slong * perm = base->perm;
    const ulong * shift = base->shift;
    const ulong * stride = base->stride;
    const ulong shiftx  = shift[perm[0]];
    const ulong stridex = stride[perm[0]];
    const ulong shifty  = shift[perm[1]];
    const ulong stridey = stride[perm[1]];
    const fmpz_mpoly_struct * B = base->B;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);
    slong xoffset, xshift;
    slong yoffset, yshift;
    slong j, k, l, arrayidx;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fmpz_mpoly_struct * Ac;
    TMP_INIT;

    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(base->Abits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    mpoly_gen_offset_shift_sp(&xoffset, &xshift, perm[0], B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&yoffset, &yshift, perm[1], B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        ulong xexp = ((B->exps[NB*j + xoffset] >> xshift) & mask) - shiftx;
        ulong yexp = ((B->exps[NB*j + yoffset] >> yshift) & mask) - shifty;
        if (stridex != 1 || stridey != 1)
        {
            FLINT_ASSERT((xexp % stridex) == 0);
            FLINT_ASSERT((yexp % stridey) == 0);
            xexp /= stridex;
            yexp /= stridey;
        }
        arrayidx = xexp*base->degby + yexp;

        FLINT_ASSERT(arrayidx < base->degbx*base->degby);
        if (base->array[arrayidx].threadidx == arg->idx)
        {
            mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
            for (k = 0; k < m + 2; k++)
            {
                l = perm[k];
                if (stride[l] == 1)
                {
                    uexps[k] = (Bexps[l] - shift[l]);
                }
                else
                {
                    FLINT_ASSERT(stride[l] != 0);
                    FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == 0);
                    uexps[k] = (Bexps[l] - shift[l]) / stride[l];
                }
            }
            Ac = base->array[arrayidx].poly;
            fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
            fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
            mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 2, base->Abits, uctx->minfo);
            Ac->length++;
        }
    }

    /* sort all of ours */
    for (arrayidx = base->degbx*base->degby - 1; arrayidx >= 0; arrayidx--)
    {
        if (base->array[arrayidx].threadidx == arg->idx)
        {
            fmpz_mpoly_sort_terms(base->array[arrayidx].poly, uctx);
        }
    }

    TMP_END;
}

/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.

    operation on each term:

    for 0 <= k < m + 2
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]

    the most significant main variable uses Aexp[0]
    the least significant main variable uses Aexp[1]
    the coefficients of A use variables Aexp[2], ..., Aexp[m + 1]
    maxexps if it exists is supposed to be a degree bound on B
*/
void fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t uctx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const ulong * maxexps, /* nullptr is ok */
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong limit = 1000; /* limit*limit should not overflow a slong */
    slong degbx, degby;
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fmpz_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(FLINT_BIT_COUNT(limit) < FLINT_BITS/2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);

    fmpz_mpolyu_zero(A, uctx);

    /* strict degree bounds on the result */
    degbx = limit + 1;
    degby = limit + 1;
    if (maxexps != NULL)
    {
        degbx = (maxexps[perm[0]] - shift[perm[0]])/stride[perm[0]] + 1;
        degby = (maxexps[perm[1]] - shift[perm[1]])/stride[perm[1]] + 1;
    }

    if (degbx <= limit && degby <= limit && degbx*degby <= limit)
    {
        _arrayconvertuu_base_t base;
        _arrayconvertuu_worker_arg_struct * args;

        base->ctx = ctx;
        base->uctx = uctx;
        base->degbx = degbx;
        base->degby = degby;
        base->perm = perm;
        base->shift = shift;
        base->stride = stride;
        base->Abits = A->bits;
        base->B = B;
        base->nthreads = num_handles + 1;
        base->array = (_arrayconvertuu_base_elem_struct *) flint_malloc(
                         degbx*degby*sizeof(_arrayconvertuu_base_elem_struct));
        for (i = degbx*degby - 1; i >= 0; i--)
        {
            base->array[i].threadidx = i % base->nthreads;
            fmpz_mpoly_init3(base->array[i].poly, 0, A->bits, uctx);
        }
        args = (_arrayconvertuu_worker_arg_struct *) flint_malloc(
                     base->nthreads*sizeof(_arrayconvertuu_worker_arg_struct));

        for (i = 0; i < num_handles; i++)
        {
            args[i].idx = i;
            args[i].base = base;
            thread_pool_wake(global_thread_pool, handles[i], 0,
                                             _arrayconvertuu_worker, &args[i]);
        }
        i = num_handles;
        args[i].idx = i;
        args[i].base = base;
        _arrayconvertuu_worker(&args[i]);
        for (i = 0; i < num_handles; i++)
        {
            thread_pool_wait(global_thread_pool, handles[i]);
        }

        A->length = 0;
        for (i = degbx - 1; i >= 0; i--)
        for (j = degby - 1; j >= 0; j--)
        {
            slong off = i*degby + j;
            if (base->array[off].poly->length > 0)
            {
                fmpz_mpolyu_fit_length(A, A->length + 1, uctx);
                A->exps[A->length] = (i << (FLINT_BITS/2)) + j;
                fmpz_mpoly_swap(A->coeffs + A->length, base->array[off].poly, uctx);
                A->length++;
            }
            fmpz_mpoly_clear(base->array[off].poly, uctx);
        }

        flint_free(base->array);
        flint_free(args);
    }
    else
    {
        TMP_START;

        uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
        Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
        NA = mpoly_words_per_exp(A->bits, uctx->minfo);
        NB = mpoly_words_per_exp(B->bits, ctx->minfo);

        for (j = 0; j < B->length; j++)
        {
            mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
            for (k = 0; k < m + 2; k++)
            {
                l = perm[k];
                if (stride[l] == 1)
                {
                    uexps[k] = (Bexps[l] - shift[l]);
                }
                else
                {
                    FLINT_ASSERT(stride[l] != 0);
                    FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == 0);
                    uexps[k] = (Bexps[l] - shift[l]) / stride[l];
                }
            }
            FLINT_ASSERT(FLINT_BIT_COUNT(uexps[0]) < FLINT_BITS/2);
            FLINT_ASSERT(FLINT_BIT_COUNT(uexps[1]) < FLINT_BITS/2);
            Ac = _fmpz_mpolyu_get_coeff(A, (uexps[0] << (FLINT_BITS/2)) + uexps[1], uctx);
            FLINT_ASSERT(Ac->bits == A->bits);

            fmpz_mpoly_fit_length(Ac, Ac->length + 1, uctx);
            fmpz_set(Ac->coeffs + Ac->length, B->coeffs + j);
            mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 2, A->bits, uctx->minfo);
            Ac->length++;
        }

        if (num_handles > 0)
        {
            _sort_arg_t arg;

#if FLINT_USES_PTHREAD
            pthread_mutex_init(&arg->mutex, NULL);
#endif
            arg->index = 0;
            arg->coeffs = A->coeffs;
            arg->length = A->length;
            arg->ctx = uctx;

            for (i = 0; i < num_handles; i++)
            {
                thread_pool_wake(global_thread_pool, handles[i], 0, _worker_sort, arg);
            }
            _worker_sort(arg);
            for (i = 0; i < num_handles; i++)
            {
                thread_pool_wait(global_thread_pool, handles[i]);
            }

#if FLINT_USES_PTHREAD
            pthread_mutex_destroy(&arg->mutex);
#endif
	}
        else
        {
            for (i = 0; i < A->length; i++)
            {
                fmpz_mpoly_sort_terms(A->coeffs + i, uctx);
            }
        }

        TMP_END;
    }
}


/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m + 2
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fmpz_mpoly_from_mpolyuu_perm_inflate( /* only for 2 main vars */
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_ctx_t uctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, k, l;
    slong NA, NB;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * uexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    fmpz_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * Bc = B->coeffs + i;
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA);
        FLINT_ASSERT(Bc->bits == B->bits);

        for (j = 0; j < Bc->length; j++)
        {
            fmpz_set(Acoeff + Alen + j, Bc->coeffs + j);
            mpoly_get_monomial_ui(uexps + 2, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[0] = B->exps[i] >> (FLINT_BITS/2);
            uexps[1] = B->exps[i] & ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2));
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k < m + 2; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }
            mpoly_set_monomial_ui(Aexp + NA*(Alen + j), Aexps, Abits, ctx->minfo);
        }
        Alen += Bc->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_mpoly_sort_terms(A, ctx);
    TMP_END;
}


void fmpz_mpolyu_fmpz_content(fmpz_t c, fmpz_mpolyu_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    fmpz_zero(c);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            fmpz_gcd(c, c, (A->coeffs + i)->coeffs + j);
            if (fmpz_is_one(c))
                break;
        }
    }
}


void fmpz_mpolyu_mul_fmpz(
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    fmpz_t c,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(!fmpz_is_zero(c));
    FLINT_ASSERT(A->bits == B->bits);
    fmpz_mpolyu_fit_length(A, B->length, ctx);

    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fmpz_mpoly_scalar_mul_fmpz(A->coeffs + i, B->coeffs + i, c, ctx);
        FLINT_ASSERT((A->coeffs + i)->bits == B->bits);
    }
    A->length = B->length;
}


void fmpz_mpolyu_divexact_fmpz(
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    fmpz_t c,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(!fmpz_is_zero(c));
    FLINT_ASSERT(A->bits == B->bits);
    fmpz_mpolyu_fit_length(A, B->length, ctx);

    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fmpz_mpoly_scalar_divexact_fmpz(A->coeffs + i, B->coeffs + i, c, ctx);
        FLINT_ASSERT((A->coeffs + i)->bits == B->bits);
    }
    A->length = B->length;
}


/*
    The bit counts of A, B and c must all agree for this division
    If saveB is zero, B may be clobbered by this operation.
*/
void fmpz_mpolyu_divexact_mpoly(
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B, int saveB,
    fmpz_mpoly_t c,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    flint_bitcnt_t exp_bits = B->bits;
    ulong * cmpmask;
    TMP_INIT;

    FLINT_ASSERT(exp_bits == A->bits);
    FLINT_ASSERT(exp_bits == B->bits);
    FLINT_ASSERT(exp_bits == c->bits);

    if (saveB == 0 && fmpz_mpoly_is_one(c, ctx))
    {
        fmpz_mpolyu_swap(A, B, ctx);
        return;
    }

    TMP_START;

    fmpz_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * poly1 = A->coeffs + i;
        fmpz_mpoly_struct * poly2 = B->coeffs + i;
        fmpz_mpoly_struct * poly3 = c;
        A->exps[i] = B->exps[i];

        fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        FLINT_ASSERT(poly2->length > 0);

        len = _fmpz_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask);
        FLINT_ASSERT(len > 0);
        _fmpz_mpoly_set_length(poly1, len, ctx);
    }
    A->length = B->length;

    TMP_END;
}

void fmpz_mpolyu_divexact_mpoly_inplace(fmpz_mpolyu_t A, fmpz_mpoly_t c,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, len;
    flint_bitcnt_t bits;
    ulong * cmpmask;
    fmpz_mpoly_t t;
    TMP_INIT;

    FLINT_ASSERT(c->length > 0);

    if (fmpz_mpoly_is_fmpz(c, ctx))
    {
        if (fmpz_is_one(c->coeffs + 0))
            return;
        for (i = 0; i < A->length; i++)
            _fmpz_vec_scalar_divexact_fmpz(A->coeffs[i].coeffs, A->coeffs[i].coeffs,
                                           A->coeffs[i].length, c->coeffs + 0);
        return;
    }

    bits = A->bits;
    FLINT_ASSERT(bits == c->bits);

    fmpz_mpoly_init3(t, 0, bits, ctx);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (i = A->length - 1; i >= 0; i--)
    {
        fmpz_mpoly_struct * poly1 = t;
        fmpz_mpoly_struct * poly2 = A->coeffs + i;
        fmpz_mpoly_struct * poly3 = c;

        FLINT_ASSERT(poly2->bits == bits);

        len = _fmpz_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, bits, N,
                                                  cmpmask);
        FLINT_ASSERT(len > 0);
        poly1->length = len;
        fmpz_mpoly_swap(poly2, poly1, ctx);
    }

    TMP_END;

    fmpz_mpoly_clear(t, ctx);
}



/*
    The bit counts of A, B and c must all agree for this multiplication
*/
void fmpz_mpolyu_mul_mpoly(fmpz_mpolyu_t A, fmpz_mpolyu_t B,
                                    fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    flint_bitcnt_t exp_bits;
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    exp_bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == c->bits);

    fmpz_mpolyu_fit_length(A, B->length, ctx);

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_struct * poly1 = A->coeffs + i;
        fmpz_mpoly_struct * poly2 = B->coeffs + i;
        fmpz_mpoly_struct * poly3 = c;
        A->exps[i] = B->exps[i];

        fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
        fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
        poly1->bits = exp_bits;

        len = _fmpz_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, exp_bits, N,
                                                  cmpmask);

        _fmpz_mpoly_set_length(poly1, len, ctx);

    }
    A->length = B->length;

    TMP_END;
}

void fmpz_mpolyu_mul_mpoly_inplace(fmpz_mpolyu_t A, fmpz_mpoly_t c,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong len;
    slong N;
    flint_bitcnt_t bits;
    ulong * cmpmask;
    fmpz_mpoly_t t;
    TMP_INIT;

    FLINT_ASSERT(c->length > 0);

    if (fmpz_mpoly_is_fmpz(c, ctx))
    {
        if (fmpz_is_one(c->coeffs + 0))
            return;
        for (i = 0; i < A->length; i++)
            _fmpz_vec_scalar_mul_fmpz(A->coeffs[i].coeffs, A->coeffs[i].coeffs,
                                           A->coeffs[i].length, c->coeffs + 0);
        return;
    }

    bits = A->bits;
    FLINT_ASSERT(bits == c->bits);

    fmpz_mpoly_init3(t, 0, bits, ctx);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (i = A->length - 1; i >= 0; i--)
    {
        fmpz_mpoly_struct * poly1 = t;
        fmpz_mpoly_struct * poly2 = A->coeffs + i;
        fmpz_mpoly_struct * poly3 = c;

        FLINT_ASSERT(poly2->bits == bits);

        len = _fmpz_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, poly2->exps, poly2->length,
                              poly3->coeffs, poly3->exps, poly3->length, bits, N,
                                                  cmpmask);
        FLINT_ASSERT(len > 0);
        poly1->length = len;
        fmpz_mpoly_swap(poly2, poly1, ctx);
    }

    TMP_END;

    fmpz_mpoly_clear(t, ctx);
}


void fmpz_mpolyu_shift_right(fmpz_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}


void fmpz_mpolyu_shift_left(fmpz_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((slong)(A->exps[i] + s) >= 0);
        A->exps[i] += s;
    }
}


void fmpz_mpolyu_content_fmpz(
    fmpz_t g,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;

    fmpz_zero(g);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_struct * Ac = A->coeffs + i;
        for (j = 0; j < Ac->length; j++)
        {
            fmpz_gcd(g, g, Ac->coeffs + j);
            if (fmpz_is_one(g))
                return;
        }
    }
}


int fmpz_mpolyu_content_mpoly_threaded_pool(
    fmpz_mpoly_t g,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, j;
    int success;
    flint_bitcnt_t bits = A->bits;

    FLINT_ASSERT(g->bits == bits);

    if (A->length < 2)
    {
        if (A->length == 0)
        {
            fmpz_mpoly_zero(g, ctx);
        }
        else
        {
            fmpz_mpoly_set(g, A->coeffs + 0, ctx);
        }

        FLINT_ASSERT(g->bits == bits);
        return 1;
    }

    j = 0;
    for (i = 1; i < A->length; i++)
    {
        if ((A->coeffs + i)->length < (A->coeffs + j)->length)
        {
            j = i;
        }
    }

    if (j == 0)
    {
        j = 1;
    }
    success = _fmpz_mpoly_gcd_threaded_pool(g, bits, A->coeffs + 0,
                                     A->coeffs + j, ctx, handles, num_handles);
    if (!success)
    {
        return 0;
    }

    for (i = 1; i < A->length; i++)
    {
        if (i == j)
        {
            continue;
        }
        success = _fmpz_mpoly_gcd_threaded_pool(g, bits, g,
                                     A->coeffs + i, ctx, handles, num_handles);
        FLINT_ASSERT(g->bits == bits);
        if (!success)
        {
            return 0;
        }
    }

    return 1;
}

