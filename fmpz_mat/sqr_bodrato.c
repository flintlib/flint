/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2015 Anubhav Srivastava

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

#define E fmpz_mat_entry

void
fmpz_mat_sqr_bodrato(fmpz_mat_t B, const fmpz_mat_t A)
{
    slong n = A->r;
    
    if (n == 0)
    {
        return;
    }
    else if (n == 1)
    {
        fmpz_mul(E(B, 0, 0), E(A, 0, 0), E(A, 0, 0));
    }
    else if (n == 2)
    {
        fmpz_add(E(B, 0, 0), E(A, 0, 0), E(A, 1, 1));

        fmpz_mul(E(B, 0, 1), E(A, 0, 1), E(B, 0, 0));
        fmpz_mul(E(B, 1, 0), E(A, 1, 0), E(B, 0, 0));

        fmpz_mul(E(B, 0, 0), E(A, 0, 0), E(A, 0, 0));
        fmpz_mul(E(B, 1, 1), E(A, 0, 1), E(A, 1, 0));
        fmpz_add(E(B, 0, 0), E(B, 0, 0), E(B, 1, 1));

        fmpz_addmul(E(B, 1, 1), E(A, 1, 1), E(A, 1, 1));
    }
    else if (n == 3)
    {
        fmpz_t temp23;
        
        fmpz_init(temp23);
       
        fmpz_mul(E(B, 2, 2), E(A, 0, 2), E(A, 2, 0));
        fmpz_mul(E(B, 1, 1), E(A, 0, 1), E(A, 1, 0));
        fmpz_mul(temp23, E(A, 1, 2), E(A, 2, 1));

        fmpz_add(E(B, 0, 0), E(B, 2, 2), E(B, 1, 1));
        fmpz_addmul(E(B, 0, 0), E(A, 0, 0), E(A, 0, 0));
       
        fmpz_add(E(B, 1, 1), E(B, 1, 1), temp23);
        fmpz_addmul(E(B, 1, 1), E(A, 1, 1), E(A, 1, 1));
        
        fmpz_add(E(B, 2, 2), E(B, 2, 2), temp23);
        fmpz_addmul(E(B, 2, 2), E(A, 2, 2), E(A, 2, 2));
      
        
        fmpz_add(E(B, 1, 2), E(A, 0, 0), E(A, 1, 1));
        fmpz_add(E(B, 2, 1), E(A, 0, 0), E(A, 2, 2));
        fmpz_add(temp23, E(A, 1, 1), E(A, 2, 2));
        
        fmpz_mul(E(B, 0, 1), E(B, 1, 2), E(A, 0, 1));
        fmpz_addmul(E(B, 0, 1), E(A, 0, 2), E(A, 2, 1));

        fmpz_mul(E(B, 0, 2), E(B, 2, 1), E(A, 0, 2));
        fmpz_addmul(E(B, 0, 2), E(A, 0, 1), E(A, 1, 2));     
 
        fmpz_mul(E(B, 1, 0), E(B, 1, 2), E(A, 1, 0));
        fmpz_addmul(E(B, 1, 0), E(A, 2, 0), E(A, 1, 2));

        fmpz_mul(E(B, 1, 2), temp23, E(A, 1, 2));
        fmpz_addmul(E(B, 1, 2), E(A, 1, 0), E(A, 0, 2));
 
        fmpz_mul(E(B, 2, 0), E(B, 2, 1), E(A, 2, 0));
        fmpz_addmul(E(B, 2, 0), E(A, 2, 1), E(A, 1, 0));
 
        fmpz_mul(E(B, 2, 1), temp23, E(A, 2, 1));
        fmpz_addmul(E(B, 2, 1), E(A, 0, 1), E(A, 2, 0));

        fmpz_clear(temp23);
    }
    else
    {
    slong a;
    slong anr;

    fmpz_mat_t A11, A12, A21, A22;
    fmpz_mat_t C11, C12, C21, C22;
    fmpz_mat_t X1, X2;

    a = A->r;

    anr = a / 2;

    fmpz_mat_window_init(A11, A, 0, 0, anr, anr);
    fmpz_mat_window_init(A12, A, 0, anr, anr, 2*anr);
    fmpz_mat_window_init(A21, A, anr, 0, 2*anr, anr);
    fmpz_mat_window_init(A22, A, anr, anr, 2*anr, 2*anr);

    fmpz_mat_window_init(C11, B, 0, 0, anr, anr);
    fmpz_mat_window_init(C12, B, 0, anr, anr, 2*anr);
    fmpz_mat_window_init(C21, B, anr, 0, 2*anr, anr);
    fmpz_mat_window_init(C22, B, anr, anr, 2*anr, 2*anr);

    fmpz_mat_init(X1, anr, anr);
    fmpz_mat_init(X2, anr, anr);

    fmpz_mat_add(X1, A22, A12);
    fmpz_mat_sqr(C21, X1);

    fmpz_mat_sub(X1, A22, A21);
    fmpz_mat_sqr(C22, X1);

    fmpz_mat_add(X1, X1, A12);
    fmpz_mat_sqr(C11, X1);

    fmpz_mat_sub(X1, X1, A11);
    fmpz_mat_mul(C12, X1, A12);
    fmpz_mat_add(C12, C12, C22);

    fmpz_mat_mul(X2, A12, A21);
    fmpz_mat_add(C11, C11, X2);
    fmpz_mat_sub(C12, C11, C12);
    fmpz_mat_sub(C11, C21, C11);
    fmpz_mat_mul(C21, A21, X1);

    fmpz_mat_clear(X1);

    fmpz_mat_sub(C21, C11, C21);
    fmpz_mat_add(C22, C22, C11);
    fmpz_mat_sqr(C11, A11);
    fmpz_mat_add(C11, C11, X2); /* C11 = P4 + P5 */

    fmpz_mat_clear(X2);

    fmpz_mat_window_clear(A11);
    fmpz_mat_window_clear(A12);
    fmpz_mat_window_clear(A21);
    fmpz_mat_window_clear(A22);

    fmpz_mat_window_clear(C11);
    fmpz_mat_window_clear(C12);
    fmpz_mat_window_clear(C21);
    fmpz_mat_window_clear(C22);

    if (a > 2*anr)
    {
        fmpz_mat_t Bc, Cc;
        fmpz_mat_window_init(Bc, A, 0, 2*anr, a, a);
        fmpz_mat_window_init(Cc, B, 0, 2*anr, a, a);
        fmpz_mat_mul(Cc, A, Bc);
        fmpz_mat_window_clear(Bc);
        fmpz_mat_window_clear(Cc);
    }

    if (a > 2*anr)
    {
        fmpz_mat_t Ar, Cr;
        fmpz_mat_window_init(Ar, A, 2*anr, 0, a, a-1);
        fmpz_mat_window_init(Cr, B, 2*anr, 0, a, a-1);
        fmpz_mat_mul(Cr, Ar, A);
        fmpz_mat_window_clear(Ar);
        fmpz_mat_window_clear(Cr);
    }

    if (a > 2*anr)
    {
        fmpz_mat_t Ac, Br, Cb, tmp;

        fmpz_mat_window_init(Ac, A, 0, 2*anr, 2*anr, a);
        fmpz_mat_window_init(Br, A, 2*anr, 0, a, 2*anr);
        fmpz_mat_window_init(Cb, B, 0, 0, 2*anr, 2*anr);

        fmpz_mat_init(tmp, 2*anr, 2*anr);
        fmpz_mat_mul(tmp, Ac, Br);
        fmpz_mat_add(Cb, Cb, tmp);
        fmpz_mat_clear(tmp);
        fmpz_mat_window_clear(Ac);
        fmpz_mat_window_clear(Br);
        fmpz_mat_window_clear(Cb);
    }

    }
}
