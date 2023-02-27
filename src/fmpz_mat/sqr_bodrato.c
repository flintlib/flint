/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2015 Anubhav Srivastava

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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
        fmpz_t t, u;

        fmpz_init(t);
        fmpz_init(u);

        fmpz_add(t, E(A, 0, 0), E(A, 1, 1));
        fmpz_mul(u, E(A, 0, 1), E(A, 1, 0));

        fmpz_mul(E(B, 0, 0), E(A, 0, 0), E(A, 0, 0));
        fmpz_add(E(B, 0, 0), E(B, 0, 0), u);

        fmpz_mul(E(B, 1, 1), E(A, 1, 1), E(A, 1, 1));
        fmpz_add(E(B, 1, 1), E(B, 1, 1), u);

        fmpz_mul(E(B, 0, 1), E(A, 0, 1), t);
        fmpz_mul(E(B, 1, 0), E(A, 1, 0), t);

        fmpz_clear(t);
        fmpz_clear(u);
    }
    else if (n == 3)
    {
        fmpz_t temp13, temp12, temp23;
        
        fmpz_init(temp13);
        fmpz_init(temp12);
        fmpz_init(temp23);
       
        fmpz_mul(temp13, E(A, 0, 2), E(A, 2, 0));
        fmpz_mul(temp12, E(A, 0, 1), E(A, 1, 0));
        fmpz_mul(temp23, E(A, 1, 2), E(A, 2, 1));

        fmpz_add(E(B, 0, 0), temp13, temp12);
        fmpz_addmul(E(B, 0, 0), E(A, 0, 0), E(A, 0, 0));
       
        fmpz_add(E(B, 1, 1), temp23, temp12);
        fmpz_addmul(E(B, 1, 1), E(A, 1, 1), E(A, 1, 1));
        
        fmpz_add(E(B, 2, 2), temp13, temp23);
        fmpz_addmul(E(B, 2, 2), E(A, 2, 2), E(A, 2, 2));
      
        
        fmpz_add(temp12, E(A, 0, 0), E(A, 1, 1));
        fmpz_add(temp13, E(A, 0, 0), E(A, 2, 2));
        fmpz_add(temp23, E(A, 1, 1), E(A, 2, 2));
        
        fmpz_mul(E(B, 0, 1), temp12, E(A, 0, 1));
        fmpz_addmul(E(B, 0, 1), E(A, 0, 2), E(A, 2, 1));

        fmpz_mul(E(B, 0, 2), temp13, E(A, 0, 2));
        fmpz_addmul(E(B, 0, 2), E(A, 0, 1), E(A, 1, 2));     
 
        fmpz_mul(E(B, 1, 0), temp12, E(A, 1, 0));
        fmpz_addmul(E(B, 1, 0), E(A, 2, 0), E(A, 1, 2));

        fmpz_mul(E(B, 1, 2), temp23, E(A, 1, 2));
        fmpz_addmul(E(B, 1, 2), E(A, 1, 0), E(A, 0, 2));
 
        fmpz_mul(E(B, 2, 0), temp13, E(A, 2, 0));
        fmpz_addmul(E(B, 2, 0), E(A, 2, 1), E(A, 1, 0));
 
        fmpz_mul(E(B, 2, 1), temp23, E(A, 2, 1));
        fmpz_addmul(E(B, 2, 1), E(A, 0, 1), E(A, 2, 0));

        fmpz_clear(temp13);
        fmpz_clear(temp23);
        fmpz_clear(temp12);
    }
    else
    {

        slong i,j;

        fmpz_mat_t window11, window12, window21, window22;
        fmpz_mat_t s1, s2, s3;
        fmpz_mat_t p1, p2, p3, p5, p6;

        slong m = n, x, iseven = 1; 

        if (n % 2 == 1)
        {
            m = n - 1;
            iseven = 0;
        }

        fmpz_mat_init(s1, m/2, m/2);
        fmpz_mat_init(s2, m/2, m/2);
        fmpz_mat_init(s3, m/2, m/2);
        fmpz_mat_init(p1, m/2, m/2);
        fmpz_mat_init(p2, m/2, m/2);
        fmpz_mat_init(p3, m/2, m/2);
        fmpz_mat_init(p5, m/2, m/2);
        fmpz_mat_init(p6, m/2, m/2);

        fmpz_mat_window_init(window11, A, 0, 0, m/2, m/2);
        fmpz_mat_window_init(window12, A, 0, m/2, m/2, m);
        fmpz_mat_window_init(window21, A, m/2, 0, m, m/2);
        fmpz_mat_window_init(window22, A, m/2, m/2, m, m);

        fmpz_mat_add(s1, window22, window12);
        fmpz_mat_sqr(p1, s1);

        fmpz_mat_sub(s2, window22, window21);
        fmpz_mat_sqr(p2, s2);

        fmpz_mat_add(s3, s2, window12);
        fmpz_mat_sqr(p3, s3);    

        fmpz_mat_sub(s1, s3, window11);
        fmpz_mat_mul(p6, s1, window12);
        fmpz_mat_mul(s3, window21, s1);

        fmpz_mat_mul(p5, window12, window21);
        fmpz_mat_add(s1, p3, p5);
        fmpz_mat_sub(s2, p1, s1);

        if (iseven == 1)
        {
            for (i = n/2; i < n; i++)
            {
                for (j = 0; j < n/2; j++)
                {
                    fmpz_sub(E(B, i, j), E(s2, i - n/2, j), E(s3, i - n/2, j));
                }
            }

            fmpz_mat_sub(s3, s1, p2);
            fmpz_mat_sqr(s1, window11);

            for (i = 0; i < n/2; i++)
            {
                for (j = 0; j < n/2; j++)
                {
                    fmpz_add(E(B, i, j), E(s1, i, j), E(p5, i, j));
                }
            }

            for (i = n/2; i < n; i++)
            {
                for (j = n/2; j < n; j++)
                {
                    fmpz_add(E(B, i, j), E(p2, i - n/2, j - n/2), E(s2, i - n/2, j - n/2));
                }
            }

            for (i = 0; i < n/2; i++)
            {
                for (j = n/2; j < n; j++)
                {
                    fmpz_sub(E(B, i, j), E(s3, i, j - n/2), E(p6, i, j - n/2) );
                }
            }
        }
        else
        {
            for (i = 0; i < n; i++)
            {
                fmpz_mul(E(B, n - 1, i), E(A, n - 1, 0), E(A, 0, i));
                for (x = 1; x < n; x++)
                {
                    fmpz_addmul(E(B, n - 1, i), E(A, n - 1, x), E(A, x, i));
                }
            }

            for (i = 0; i < n; i++)
            {
                fmpz_mul(E(B, i, n - 1), E(A, 0, n - 1), E(A, i, 0));
                for (x = 1; x < n; x++)
                {
                    fmpz_addmul(E(B, i, n - 1), E(A, x, n - 1), E(A, i, x));
                }
            }


            for (i = m/2; i < m; i++)
            {
                for (j = 0; j < m/2; j++)
                {
                    fmpz_sub(E(B, i, j), E(s2, i - m/2, j), E(s3, i - m/2, j)); 
                    fmpz_addmul(E(B, i, j), E(A, i, n - 1), E(A, n - 1, j));
                }
            }

            fmpz_mat_sub(s3, s1, p2);
            fmpz_mat_sqr(s1, window11);


            for (i = 0; i < m/2; i++)
            {
                for (j = 0; j < m/2; j++)
                {
                    fmpz_add(E(B, i, j), E(s1, i, j), E(p5, i, j)); 
                    fmpz_addmul(E(B, i, j), E(A, i, n - 1), E(A, n - 1, j));
                }
            }

            for (i = m/2; i < m; i++)
            {
                for (j = m/2; j < m; j++)
                {
                    fmpz_add(E(B, i, j), E(p2, i - m/2, j - m/2), E(s2, i - m/2, j - m/2)); 
                    fmpz_addmul(E(B, i, j), E(A, i, n - 1), E(A, n - 1, j));
                }
            }

            for (i = 0; i < m/2; i++)
            {
                for (j = m/2; j < m; j++)
                {
                    fmpz_sub(E(B, i, j), E(s3, i, j - m/2), E(p6, i, j - m/2)); 
                    fmpz_addmul(E(B, i, j), E(A, i, n - 1), E(A, n - 1, j));
                }
            }

        }

        fmpz_mat_window_clear(window11);
        fmpz_mat_window_clear(window12);
        fmpz_mat_window_clear(window21);
        fmpz_mat_window_clear(window22);
        fmpz_mat_clear(s1);
        fmpz_mat_clear(s2);
        fmpz_mat_clear(s3);
        fmpz_mat_clear(p1);
        fmpz_mat_clear(p2);
        fmpz_mat_clear(p3);
        fmpz_mat_clear(p5);
        fmpz_mat_clear(p6);
    }
}
