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

    Copyright (C) 2008, Martin Albrecht
    Copyright (C) 2008, 2009 William Hart.
    Copyright (C) 2010, Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


void
nmod_mat_mul_strassen(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    long a, b, c;
    long anr, anc, bnr, bnc;

    nmod_mat_t A11, A12, A21, A22;
    nmod_mat_t B11, B12, B21, B22;
    nmod_mat_t C11, C12, C21, C22;
    nmod_mat_t X1, X2;

    a = A->r;
    b = A->c;
    c = B->c;

    if (a <= 4 || b <= 4 || c <= 4)
    {
        nmod_mat_mul(C, A, B);
        return;
    }

    anr = a / 2;
    anc = b / 2;
    bnr = anc;
    bnc = c / 2;

    nmod_mat_window_init(A11, A, 0, 0, anr, anc);
    nmod_mat_window_init(A12, A, 0, anc, anr, 2*anc);
    nmod_mat_window_init(A21, A, anr, 0, 2*anr, anc);
    nmod_mat_window_init(A22, A, anr, anc, 2*anr, 2*anc);

    nmod_mat_window_init(B11, B, 0, 0, bnr, bnc);
    nmod_mat_window_init(B12, B, 0, bnc, bnr, 2*bnc);
    nmod_mat_window_init(B21, B, bnr, 0, 2*bnr, bnc);
    nmod_mat_window_init(B22, B, bnr, bnc, 2*bnr, 2*bnc);

    nmod_mat_window_init(C11, C, 0, 0, anr, bnc);
    nmod_mat_window_init(C12, C, 0, bnc, anr, 2*bnc);
    nmod_mat_window_init(C21, C, anr, 0, 2*anr, bnc);
    nmod_mat_window_init(C22, C, anr, bnc, 2*anr, 2*bnc);

    nmod_mat_init(X1, anr, FLINT_MAX(bnc, anc), A->mod.n);
    nmod_mat_init(X2, anc, bnc, A->mod.n);

    X1->c = anc;

    /*
        See Jean-Guillaume Dumas, Clement Pernet, Wei Zhou; "Memory
        efficient scheduling of Strassen-Winograd's matrix multiplication
        algorithm"; http://arxiv.org/pdf/0707.2347v3 for reference on the
        used operation scheduling.
    */

    nmod_mat_sub(X1, A11, A21);
    nmod_mat_sub(X2, B22, B12);
    nmod_mat_mul(C21, X1, X2);

    nmod_mat_add(X1, A21, A22);
    nmod_mat_sub(X2, B12, B11);
    nmod_mat_mul(C22, X1, X2);

    nmod_mat_sub(X1, X1, A11);
    nmod_mat_sub(X2, B22, X2);
    nmod_mat_mul(C12, X1, X2);

    nmod_mat_sub(X1, A12, X1);
    nmod_mat_mul(C11, X1, B22);

    X1->c = bnc;
    nmod_mat_mul(X1, A11, B11);

    nmod_mat_add(C12, X1, C12);
    nmod_mat_add(C21, C12, C21);
    nmod_mat_add(C12, C12, C22);
    nmod_mat_add(C22, C21, C22);
    nmod_mat_add(C12, C12, C11);
    nmod_mat_sub(X2, X2, B21);
    nmod_mat_mul(C11, A22, X2);

    nmod_mat_clear(X2);

    nmod_mat_sub(C21, C21, C11);
    nmod_mat_mul(C11, A12, B21);

    nmod_mat_add(C11, X1, C11);

    nmod_mat_clear(X1);

    nmod_mat_window_clear(A11);
    nmod_mat_window_clear(A12);
    nmod_mat_window_clear(A21);
    nmod_mat_window_clear(A22);

    nmod_mat_window_clear(B11);
    nmod_mat_window_clear(B12);
    nmod_mat_window_clear(B21);
    nmod_mat_window_clear(B22);

    nmod_mat_window_clear(C11);
    nmod_mat_window_clear(C12);
    nmod_mat_window_clear(C21);
    nmod_mat_window_clear(C22);

    if (c > 2*bnc) /* A by last col of B -> last col of C */
    {
        nmod_mat_t Bc, Cc;
        nmod_mat_window_init(Bc, B, 0, 2*bnc, b, c);
        nmod_mat_window_init(Cc, C, 0, 2*bnc, a, c);
        nmod_mat_mul(Cc, A, Bc);
        nmod_mat_window_clear(Bc);
        nmod_mat_window_clear(Cc);
    }

    if (a > 2*anr) /* last row of A by B -> last row of C */
    {
        nmod_mat_t Ar, Cr;
        nmod_mat_window_init(Ar, A, 2*anr, 0, a, b);
        nmod_mat_window_init(Cr, C, 2*anr, 0, a, c);
        nmod_mat_mul(Cr, Ar, B);
        nmod_mat_window_clear(Ar);
        nmod_mat_window_clear(Cr);
    }

    if (b > 2*anc) /* last col of A by last row of B -> C */
    {
        nmod_mat_t Ac, Br, Cb;
        nmod_mat_window_init(Ac, A, 0, 2*anc, 2*anr, b);
        nmod_mat_window_init(Br, B, 2*bnr, 0, b, 2*bnc);
        nmod_mat_window_init(Cb, C, 0, 0, 2*anr, 2*bnc);
        nmod_mat_addmul(Cb, Cb, Ac, Br);
        nmod_mat_window_clear(Ac);
        nmod_mat_window_clear(Br);
        nmod_mat_window_clear(Cb);
    }
}
