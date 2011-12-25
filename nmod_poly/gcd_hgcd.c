/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

#define NMOD_POLY_HGCD_CUTOFF 60        /* HGCD: Basecase -> Recursion */
#define NMOD_POLY_GCD_CUTOFF 184        /* GCD:  Euclidean -> HGCD */
#define NMOD_POLY_SMALL_GCD_CUTOFF 174  /* GCD (small n): Euclidean -> HGCD */


static __inline__ 
void _nmod_poly_attach(nmod_poly_t output, nmod_poly_t input)
{
   output->length = input->length;
   output->coeffs = input->coeffs;
   output->mod = input->mod;
}

static __inline__ 
void nmod_poly_attach(nmod_poly_t output, nmod_poly_t input)
{
   _nmod_poly_attach(output, input);
}

/*
   Attach input shifted right by n to output
*/

static __inline__ 
void _nmod_poly_attach_shift(nmod_poly_t output, 
                   nmod_poly_t input, long n)
{
   if (input->length >= n) output->length = input->length - n;
   else output->length = 0;
   output->coeffs = input->coeffs + n;
   output->mod = input->mod;
}

static __inline__ 
void nmod_poly_attach_shift(nmod_poly_t output, 
                   nmod_poly_t input, long n)
{
   _nmod_poly_attach_shift(output, input, n);
}

/*
   Attach input to first n coefficients of input
*/

static __inline__ 
void _nmod_poly_attach_truncate(nmod_poly_t output, 
                     nmod_poly_t input, long n)
{
   if (input->length < n) output->length = input->length;
   else output->length = n;
   output->coeffs = input->coeffs;
   output->mod = input->mod;
}

static __inline__ 
void nmod_poly_attach_truncate(nmod_poly_t output, 
                     nmod_poly_t input, long n)
{
   _nmod_poly_attach_truncate(output, input, n);
}

long nmod_poly_half_gcd_iter(nmod_poly_mat_t res, 
    nmod_poly_t a2, nmod_poly_t b2, const nmod_poly_t a, const nmod_poly_t b)
{
    const long m = a->length / 2;

    nmod_poly_mat_one(res);
    nmod_poly_set(a2, a);
    nmod_poly_set(b2, b);

    if (b->length < m + 1) 
    {
        return 1L;
    }
    else
    {
        nmod_poly_t Q, temp;
        long sign = 1L;

        nmod_poly_init_preinv(Q, a->mod.n, a->mod.ninv);
        nmod_poly_init_preinv(temp, a->mod.n, a->mod.ninv);

        while (b2->length >= m + 1)
        {
            nmod_poly_divrem(Q, a2, a2, b2);
            nmod_poly_swap(a2, b2);

            nmod_poly_mul(temp, Q, nmod_poly_mat_entry(res,1,0));
            nmod_poly_add(temp, nmod_poly_mat_entry(res,1,1), temp);
            nmod_poly_swap(nmod_poly_mat_entry(res,1,1), nmod_poly_mat_entry(res,1,0));
            nmod_poly_swap(nmod_poly_mat_entry(res,1,0), temp);

            nmod_poly_mul(temp, Q, nmod_poly_mat_entry(res,0,0));
            nmod_poly_add(temp, nmod_poly_mat_entry(res,0,1), temp);
            nmod_poly_swap(nmod_poly_mat_entry(res,0,1), nmod_poly_mat_entry(res,0,0));
            nmod_poly_swap(nmod_poly_mat_entry(res,0,0), temp);
            sign = -sign;
        }

        nmod_poly_clear(temp);
        nmod_poly_clear(Q);

        return sign;
    }
}

long nmod_poly_half_gcd(nmod_poly_mat_t R, nmod_poly_t A, nmod_poly_t B, 
    const nmod_poly_t a, const nmod_poly_t b, int flag)
{
    const long m = a->length/2;

    if (b->length < m + 1)
    {
        nmod_poly_mat_one(R);
        nmod_poly_set(A, a);
        nmod_poly_set(B, b);
        return 1L;
    }
    else
    {
        long signR;
        nmod_poly_t a0, b0;
        nmod_poly_t temp, a2, b2; 
        nmod_poly_t a3, b3, a4, b4, s, t;

        _nmod_poly_attach_shift(a0, a, m);
        _nmod_poly_attach_shift(b0, b, m);

        nmod_poly_init_preinv(temp, a->mod.n, a->mod.ninv);
        nmod_poly_init_preinv(a2, a->mod.n, a->mod.ninv);
        nmod_poly_init_preinv(b2, a->mod.n, a->mod.ninv);

        nmod_poly_init_preinv(a3, a->mod.n, a->mod.ninv);
        nmod_poly_init_preinv(b3, a->mod.n, a->mod.ninv);

        if (a0->length < NMOD_POLY_HGCD_CUTOFF) 
            signR = nmod_poly_half_gcd_iter(R, a3, b3, a0, b0);
        else 
            signR = nmod_poly_half_gcd(R, a3, b3, a0, b0, 1);

        nmod_poly_attach_truncate(s, a, m);
        nmod_poly_attach_truncate(t, b, m);

        nmod_poly_mul(b2, nmod_poly_mat_entry(R,1,0), s);
        nmod_poly_mul(temp, nmod_poly_mat_entry(R,0,0), t);

        if (signR < 0L) 
            nmod_poly_sub(b2, b2, temp);
        else 
            nmod_poly_sub(b2, temp, b2);

        nmod_poly_fit_length(b2, m + b3->length);

        long i;
        for (i = b2->length; i < m + b3->length; i++) 
            b2->coeffs[i] = 0L;

        nmod_poly_attach_shift(b4, b2, m);
        b4->alloc = FLINT_MAX(b4->length, b3->length);
        nmod_poly_add(b4, b4, b3);
        b2->length = FLINT_MAX(m + b3->length, b2->length);
        _nmod_poly_normalise(b2);

        nmod_poly_mul(a2, nmod_poly_mat_entry(R,1,1), s);
        nmod_poly_mul(temp, nmod_poly_mat_entry(R,0,1), t);

        if (signR < 0L) 
            nmod_poly_sub(a2, temp, a2);
        else 
            nmod_poly_sub(a2, a2, temp);

        nmod_poly_fit_length(a2, m + a3->length);
        for (i = a2->length; i < m + a3->length; i++) 
            a2->coeffs[i] = 0L;
        nmod_poly_attach_shift(a4, a2, m);
        a4->alloc = FLINT_MAX(a4->length, a3->length);
        nmod_poly_add(a4, a4, a3);
        a2->length = FLINT_MAX(m + a3->length, a2->length);
        _nmod_poly_normalise(a2);

        if (b2->length < m + 1)
        {
            nmod_poly_set(A, a2);
            nmod_poly_set(B, b2);
            nmod_poly_clear(temp);
            nmod_poly_clear(a2);
            nmod_poly_clear(b2);
            nmod_poly_clear(a3);
            nmod_poly_clear(b3);
            return signR;
        }

        nmod_poly_t q, d;
        nmod_poly_init(q, a->mod.n);
        nmod_poly_init(d, a->mod.n);

        nmod_poly_divrem(q, d, a2, b2);

        long k = 2*m - b2->length + 1;

        nmod_poly_t c0, d0;
        _nmod_poly_attach_shift(c0, b2, k);
        _nmod_poly_attach_shift(d0, d, k);
        nmod_poly_mat_t S;
        nmod_poly_mat_init(S, 2, 2, a->mod.n);

        long signS;

        if (c0->length < NMOD_POLY_HGCD_CUTOFF) 
            signS = nmod_poly_half_gcd_iter(S, a3, b3, c0, d0);
        else 
            signS = nmod_poly_half_gcd(S, a3, b3, c0, d0, 1);

        nmod_poly_attach_truncate(s, b2, k);
        nmod_poly_attach_truncate(t, d, k);

        nmod_poly_mul(B, nmod_poly_mat_entry(S,1,0), s);
        nmod_poly_mul(temp, nmod_poly_mat_entry(S,0,0), t);

        if (signS < 0L) 
            nmod_poly_sub(B, B, temp);
        else 
            nmod_poly_sub(B, temp, B);

        nmod_poly_fit_length(B, k + b3->length);
        for (i = B->length; i < k + b3->length; i++) 
            B->coeffs[i] = 0L;
        nmod_poly_attach_shift(b4, B, k);
        b4->alloc = FLINT_MAX(b4->length, b3->length);
        nmod_poly_add(b4, b4, b3);
        B->length = FLINT_MAX(k + b3->length, B->length);
        _nmod_poly_normalise(B);

        nmod_poly_mul(A, nmod_poly_mat_entry(S,1,1), s);
        nmod_poly_mul(temp, nmod_poly_mat_entry(S,0,1), t);

        if (signS < 0L) 
            nmod_poly_sub(A, temp, A);
        else 
            nmod_poly_sub(A, A, temp);

        nmod_poly_fit_length(A, k + a3->length);
        for (i = A->length; i < k + a3->length; i++) 
            A->coeffs[i] = 0L;
        nmod_poly_attach_shift(a4, A, k);
        a4->alloc = FLINT_MAX(a4->length, a3->length);
        nmod_poly_add(a4, a4, a3);
        A->length = FLINT_MAX(k + a3->length, A->length);
        _nmod_poly_normalise(A);

if (flag) {

        nmod_poly_swap(nmod_poly_mat_entry(S,0,0), nmod_poly_mat_entry(S,1,0));
        nmod_poly_swap(nmod_poly_mat_entry(S,0,1), nmod_poly_mat_entry(S,1,1));
        nmod_poly_mul(temp, nmod_poly_mat_entry(S,1,0), q);
        nmod_poly_add(nmod_poly_mat_entry(S,0,0), nmod_poly_mat_entry(S,0,0), temp);
        nmod_poly_mul(temp, nmod_poly_mat_entry(S,1,1), q);
        nmod_poly_add(nmod_poly_mat_entry(S,0,1), nmod_poly_mat_entry(S,0,1), temp);

        nmod_poly_mat_mul(R, R, S);
}

        nmod_poly_mat_clear(S);
        nmod_poly_clear(temp);
        nmod_poly_clear(a3); 
        nmod_poly_clear(b3);
        nmod_poly_clear(q);
        nmod_poly_clear(d);
        nmod_poly_clear(a2);
        nmod_poly_clear(b2);

        return - (signR * signS);
    }
}

void nmod_poly_half_gcd_no_matrix(nmod_poly_t a_out, nmod_poly_t b_out, nmod_poly_t a, nmod_poly_t b)
{
    nmod_poly_mat_t mat;

    nmod_poly_mat_init(mat, 2, 2, a->mod.n);
    nmod_poly_half_gcd(mat, a_out, b_out, a, b, 0);
    
    nmod_poly_mat_clear(mat);
}

void nmod_poly_gcd_hgcd(nmod_poly_t res, nmod_poly_t f, const nmod_poly_t g)
{
    const nmod_t mod  = f->mod;
    const long CUTOFF = FLINT_BIT_COUNT(mod.n) <= 8 ? 
    NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;
    nmod_poly_t h, j, r;

    if (f->length == 0)
    {
        if (g->length == 0) nmod_poly_zero(res);
        else nmod_poly_make_monic(res, g);
        return;
    }
    if (g->length == 0)
    {
        nmod_poly_make_monic(res, f);
        return;
    }

    nmod_poly_init_preinv(r, mod.n, mod.ninv);

    nmod_poly_rem(r, f, g);
    if (r->length == 0)
    {
        nmod_poly_make_monic(res, g);
        nmod_poly_clear(r);
        return;
    }

    nmod_poly_init_preinv(j, mod.n, mod.ninv);
    nmod_poly_init_preinv(h, mod.n, mod.ninv);

    nmod_poly_half_gcd_no_matrix(h, j, g, r);

    while (j->length != 0)
    {
        nmod_poly_rem(r, h, j);

        if (r->length == 0)
        {
            nmod_poly_make_monic(res, j);
            goto gcd_exit;
        }

        if (j->length < CUTOFF)
        {
            nmod_poly_gcd_euclidean(res, j, r);
            goto gcd_exit;
        }

        nmod_poly_half_gcd_no_matrix(h, j, j, r);
    }

    nmod_poly_make_monic(res, h);

  gcd_exit: 

    nmod_poly_clear(j);
    nmod_poly_clear(h);
    nmod_poly_clear(r);
}
