/*
    Copyright (C) 2016 Aaditya Thakkar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "flint.h"

void fmpz_mat_mul_strassen(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong a, b, c;
    slong anr, anc, bnr, bnc;

    fmpz_mat_t A11, A12, A21, A22;
    fmpz_mat_t B11, B12, B21, B22;
    fmpz_mat_t C11, C12, C21, C22;
    fmpz_mat_t X1, X2;

    a = A->r;
    b = A->c;
    c = B->c;

    if (a <= 4 || b <= 4 || c <= 4)
    {
        fmpz_mat_mul(C, A, B);
        return;
    }

    anr = a / 2;
    anc = b / 2;
    bnr = anc;
    bnc = c / 2;

    fmpz_mat_window_init(A11, A, 0, 0, anr, anc);
    fmpz_mat_window_init(A12, A, 0, anc, anr, 2*anc);
    fmpz_mat_window_init(A21, A, anr, 0, 2*anr, anc);
    fmpz_mat_window_init(A22, A, anr, anc, 2*anr, 2*anc);

    fmpz_mat_window_init(B11, B, 0, 0, bnr, bnc);
    fmpz_mat_window_init(B12, B, 0, bnc, bnr, 2*bnc);
    fmpz_mat_window_init(B21, B, bnr, 0, 2*bnr, bnc);
    fmpz_mat_window_init(B22, B, bnr, bnc, 2*bnr, 2*bnc);

    fmpz_mat_window_init(C11, C, 0, 0, anr, bnc);
    fmpz_mat_window_init(C12, C, 0, bnc, anr, 2*bnc);
    fmpz_mat_window_init(C21, C, anr, 0, 2*anr, bnc);
    fmpz_mat_window_init(C22, C, anr, bnc, 2*anr, 2*bnc);

    fmpz_mat_init(X1, anr, FLINT_MAX(bnc, anc));
    fmpz_mat_init(X2, anc, bnc);

    X1->c = anc;

    fmpz_mat_sub(X1, A11, A21);
    fmpz_mat_sub(X2, B22, B12);
    fmpz_mat_mul(C21, X1, X2);

    fmpz_mat_add(X1, A21, A22);
    fmpz_mat_sub(X2, B12, B11);
    fmpz_mat_mul(C22, X1, X2);

    fmpz_mat_sub(X1, X1, A11);
    fmpz_mat_sub(X2, B22, X2);
    fmpz_mat_mul(C12, X1, X2);

    fmpz_mat_sub(X1, A12, X1);
    fmpz_mat_mul(C11, X1, B22);

    X1->c = bnc;
    fmpz_mat_mul(X1, A11, B11);
    fmpz_mat_add(C12, X1, C12);
    fmpz_mat_add(C21, C12, C21);
    fmpz_mat_add(C12, C12, C22);
    fmpz_mat_add(C22, C21, C22);
    fmpz_mat_add(C12, C12, C11);
    fmpz_mat_sub(X2, X2, B21);
    fmpz_mat_mul(C11, A22, X2);

    fmpz_mat_clear(X2);

    fmpz_mat_sub(C21, C21, C11);
    fmpz_mat_mul(C11, A12, B21);

    fmpz_mat_add(C11, X1, C11);

    X1->c = FLINT_MAX(bnc, anc);
    fmpz_mat_clear(X1);

    fmpz_mat_window_clear(A11);
    fmpz_mat_window_clear(A12);
    fmpz_mat_window_clear(A21);
    fmpz_mat_window_clear(A22);

    fmpz_mat_window_clear(B11);
    fmpz_mat_window_clear(B12);
    fmpz_mat_window_clear(B21);
    fmpz_mat_window_clear(B22);

    fmpz_mat_window_clear(C11);
    fmpz_mat_window_clear(C12);
    fmpz_mat_window_clear(C21);
    fmpz_mat_window_clear(C22);

    if (c > 2*bnc) 
    {
        fmpz_mat_t Bc, Cc;
        fmpz_mat_window_init(Bc, B, 0, 2*bnc, b, c);
        fmpz_mat_window_init(Cc, C, 0, 2*bnc, a, c);
        fmpz_mat_mul(Cc, A, Bc);
        fmpz_mat_window_clear(Bc);
        fmpz_mat_window_clear(Cc);
    }

    if (a > 2*anr) 
    {
        fmpz_mat_t Ar, Cr;
        fmpz_mat_window_init(Ar, A, 2*anr, 0, a, b);
        fmpz_mat_window_init(Cr, C, 2*anr, 0, a, c);
        fmpz_mat_mul(Cr, Ar, B);
        fmpz_mat_window_clear(Ar);
        fmpz_mat_window_clear(Cr);
    }

    if (b > 2*anc) 
    {
        fmpz_mat_t Ac, Br, Cb, tmp;
        slong mt, nt;

        fmpz_mat_window_init(Ac, A, 0, 2*anc, 2*anr, b);
        fmpz_mat_window_init(Br, B, 2*bnr, 0, b, 2*bnc);
        fmpz_mat_window_init(Cb, C, 0, 0, 2*anr, 2*bnc);
        
        mt = Ac->r;
        nt = Br->c;
        
        fmpz_mat_init(tmp, mt, nt);
        fmpz_mat_mul(tmp, Ac, Br);
        fmpz_mat_add(Cb, Cb, tmp); 
        fmpz_mat_clear(tmp);
        fmpz_mat_window_clear(Ac);
        fmpz_mat_window_clear(Br);
        fmpz_mat_window_clear(Cb);
    }
}
