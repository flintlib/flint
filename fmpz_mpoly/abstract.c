/*
    Copyright (C) 2021 Brent Baccala

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

/* This subroutine adds a set of polynomials by pulling from an input
 * function into a heap structure and writing an output polynomial.
 *
 * output_function is called every time we output a term
 *
 * input_function
 *
 *    copies an exponent vector and a coefficient into the provided pointer arguments
 *
 */

void _default_output_function(void * poly, slong index, const flint_bitcnt_t bits,
                              ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    fmpz_mpoly_struct * A = (fmpz_mpoly_struct *) poly;

    fmpz_mpoly_fit_length(A, index, ctx);
    mpoly_monomial_set(A->exps + index*N, exp, N);
    fmpz_swap(A->coeffs + index, coeff);
    _fmpz_mpoly_set_length(A, index, ctx);
}

void _default_input_function(void * poly, slong index, const flint_bitcnt_t bits,
                             ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    fmpz_mpoly_struct * B = (fmpz_mpoly_struct *) poly;

    mpoly_monomial_set(exp, B->exps + index*N, N);
    fmpz_set(coeff, B->coeffs + index);
}

void fmpz_mpoly_abstract_add(
    void * A,
    void ** Blist,
    const slong Blen,
    const flint_bitcnt_t bits,
    const fmpz_mpoly_ctx_t ctx,
    const void (* input_function)(void * poly, slong index, const flint_bitcnt_t bits,
                                  ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx),
    const void (* output_function)(void * poly, slong index, const flint_bitcnt_t bits,
                                   ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx))
{
    slong i,j,k;
    slong N;
    ulong * cmpmask;
    slong heaplen;
    ulong * exps;
    fmpz * coeffs;
    ulong * Blength;
    int first;
    TMP_INIT;

    TMP_START;

    if (input_function == NULL)
        input_function = _default_input_function;
    if (output_function == NULL)
        output_function = _default_output_function;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (k = 0; ; k++) /* round up numterms to a power of 2 */
    {
        if (Blen <= (1<<k))
            break;
    }

    /* (1<<k) is the number of input slots we need at the end of the heap */
    /* (2<<k) is the required size of the heap */

    heaplen = 2<<k;

    coeffs = (fmpz *) TMP_ALLOC(heaplen * sizeof(fmpz));
    exps = (ulong *) TMP_ALLOC(heaplen * N * sizeof(ulong));

    for (i=0; i<heaplen; i++)
    {
        fmpz_init(coeffs + i);
    }

    Blength = (ulong *) TMP_ALLOC(Blen * sizeof(ulong));

    /* start by filling in the back half of the heap with data */

    for (i=0; i < Blen; i++)
    {
        Blength[i] = 0;
        input_function(Blist[i], Blength[i] ++, bits, exps + N*(heaplen/2 + i), coeffs + (heaplen/2 + i), ctx);
    }

    /* now run the heap algorithm on the front half to fill it in */

    for (k = heaplen/2 - 1; k > 0; k--)
    {
        i = k;
        while ((j = HEAP_LEFT(i)) < heaplen)
        {
            if (fmpz_is_zero(coeffs + j)
                || (!fmpz_is_zero(coeffs + j + 1)
                    && !mpoly_monomial_gt(exps + N*j, exps + N*(j+1), N, cmpmask)))
                j ++;
            fmpz_set(coeffs + i, coeffs + j);
            mpoly_monomial_set(exps + N*i, exps + N*j, N);
            i = j;

            if ((j >= heaplen/2) && (j - heaplen/2 < Blen) && !fmpz_is_zero(coeffs + j))
                input_function(Blist[j - heaplen/2], Blength[j - heaplen/2] ++, bits, exps + N*j, coeffs + j, ctx);
        }
    }

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

    while (! fmpz_is_zero(coeffs + 1))
    {
        k ++;

        mpoly_monomial_set(exps + N*0, exps + N*1, N);

        first = 1;

        while (!fmpz_is_zero(coeffs + 1) && mpoly_monomial_equal(exps + N*0, exps + N*1, N))
        {
            if (first)
                fmpz_swap(coeffs + 0, coeffs + 1);
            else
                fmpz_add(coeffs + 0, coeffs + 0, coeffs + 1);

            first = 0;

            /* pop lead item from heap and propagate items through the heap */
            i = 1;

            while ((j = HEAP_LEFT(i)) < heaplen)
            {
                if (fmpz_is_zero(coeffs + j)
                    || (!fmpz_is_zero(coeffs + j + 1)
                        && !mpoly_monomial_gt(exps + N*j, exps + N*(j+1), N, cmpmask)))
                    j ++;
                fmpz_set(coeffs + i, coeffs + j);
                mpoly_monomial_set(exps + N*i, exps + N*j, N);
                i = j;

                if ((j >= heaplen/2) && (j - heaplen/2 < Blen) && !fmpz_is_zero(coeffs + j))
                    input_function(Blist[j - heaplen/2], Blength[j - heaplen/2], bits, exps + N*j, coeffs + j, ctx);
            }
        }

        if (fmpz_is_zero(coeffs + 0))
            k --;
        else
            output_function(A, k, bits, exps + 0, coeffs + 0, ctx);
    }

    for (i=0; i<heaplen; i++)
    {
        fmpz_clear(coeffs + i);
    }

    TMP_END;
}

/*
 * output_function is called every time we're ready to output a term
 *
 *    what does it do if the output buffer is full?
 *
 * input_function
 *
 *    sets pointers to a block of exponents and coefficients and returns the length
 *
 *    what does it do if the input buffer is empty?
 *        block?
 *        return 0?
 *
 *    how do we know when a block of exponents and coefficients can be freed?  who frees it?
 *
 * simplified input_function
 *
 *    copies an exponent vector and a coefficient into the provided pointer arguments
 *
 * chain of abstract operations
 *
 *  [mul] -\
 *          \
 *  [mul] ---+- [add] -- [buffer] -+- [mul] -- [write]
 *          /                     /
 *  [mul] -/             [poly] -/
 *
 *
 * highest priority is the final write; always run it if we're able
 * then run final mul if it's got valid input
 * then run add if it's got valid input
 * then run first mul functions until they've buffered enough to run the add,
 *    prioritizing the ones with the least data in their output buffers
 *
 * array of structures with block processing functions and pointers to their input and output buffers
 *
 * another possibility it to use threading at the Python level.  Once the add's output buffer is sufficiently
 * full, start another thread to execute the final mul
 */
