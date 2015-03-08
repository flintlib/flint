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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpz_mat.h"

#define E fmpz_mat_entry

void
fmpz_mat_sqr_bodrato(fmpz_mat_t B, const fmpz_mat_t A)
{
    slong n = A->r;
    slong i,j;
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
    else
    { 
        if (n%2 == 0)
        {
            fmpz_mat_t window11, window12, window21, window22;
            fmpz_mat_t s1, s2, s3, s4;
            fmpz_mat_t p1, p2, p3, p4, p5, p6, p7;

            fmpz_mat_init(s1, n/2, n/2);
            fmpz_mat_init(s2, n/2, n/2);
            fmpz_mat_init(s3, n/2, n/2);
            fmpz_mat_init(s4, n/2, n/2);
            fmpz_mat_init(p1, n/2, n/2);
            fmpz_mat_init(p2, n/2, n/2);
            fmpz_mat_init(p3, n/2, n/2);
            fmpz_mat_init(p4, n/2, n/2);
            fmpz_mat_init(p5, n/2, n/2);
            fmpz_mat_init(p6, n/2, n/2);
            fmpz_mat_init(p7, n/2, n/2);

            fmpz_mat_window_init(window11, A, 0, 0, n/2, n/2);
            fmpz_mat_window_init(window12, A, 0, n/2, n/2, n);
            fmpz_mat_window_init(window21, A, n/2, 0, n, n/2);
            fmpz_mat_window_init(window22, A, n/2, n/2, n, n);
          
            
            fmpz_mat_add(s1, window22, window12);
            fmpz_mat_sub(s2, window22, window21);
            fmpz_mat_add(s3, s2, window12);
            fmpz_mat_sub(s4, s3, window11);
       

            fmpz_mat_sqr_bodrato(p1, s1);
            fmpz_mat_sqr_bodrato(p2, s2);
            fmpz_mat_sqr_bodrato(p3, s3);    
            fmpz_mat_sqr_bodrato(p4, window11);
            fmpz_mat_mul(p5, window12, window21);
            fmpz_mat_mul(p6, s4, window12);
            fmpz_mat_mul(p7, window21, s4);
            

            fmpz_mat_zero(s1);
            fmpz_mat_zero(s2);
            fmpz_mat_zero(s3);
           
            fmpz_mat_add(s1, p3, p5);
            fmpz_mat_sub(s2, p1, s1);
            fmpz_mat_sub(s3, s1, p2);
            
            fmpz_mat_zero(s1);
            fmpz_mat_zero(p1);
            fmpz_mat_zero(p3);
            fmpz_mat_zero(s4);
            
            fmpz_mat_add(p1, p4, p5);
            fmpz_mat_sub(s1, s3, p6);
            fmpz_mat_sub(p3, s2, p7);
            fmpz_mat_add(s4, p2, s2);
            
            for (i = 0; i < n/2; ++i)
            {
                for (j = 0; j < n/2; ++j)
                {
                    fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(p1, i, j));
                }
            }

            for (i = n/2; i < n; ++i)
            {
                for (j = 0; j < n/2; ++j)
                {
                    fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(p3, i-n/2, j));
                }
            }

            for (i = 0; i < n/2; ++i)
            {
                for (j = n/2; j < n; ++j)
                {
                    fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(s1, i, j-n/2));
                }
            }

            for (i = n/2; i < n; ++i)
            {
                for (j = n/2; j < n; ++j)
                {
                    fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(s4, i-n/2, j-n/2));
                }
            }
   
            fmpz_mat_window_clear(window11);
            fmpz_mat_window_clear(window12);
            fmpz_mat_window_clear(window21);
            fmpz_mat_window_clear(window22);
 
            fmpz_mat_clear(s1);
            fmpz_mat_clear(s2);
            fmpz_mat_clear(s3);
            fmpz_mat_clear(s4);
            fmpz_mat_clear(p1);
            fmpz_mat_clear(p2);
            fmpz_mat_clear(p3);
            fmpz_mat_clear(p4);
            fmpz_mat_clear(p5);
            fmpz_mat_clear(p6);
            fmpz_mat_clear(p7);

        }
        else
        {
            fmpz_mat_t window_A;
            fmpz_mat_window_init(window_A, A, 0, 0, n-1, n-1);
            
            slong m = window_A->r, x;
    
            fmpz_t sum, val, temp;
 
            fmpz_mat_t window11, window12, window21, window22;
            fmpz_mat_t s1, s2, s3, s4;
            fmpz_mat_t p1, p2, p3, p4, p5, p6, p7;

            fmpz_mat_init(s1, m/2, m/2);
            fmpz_mat_init(s2, m/2, m/2);
            fmpz_mat_init(s3, m/2, m/2);
            fmpz_mat_init(s4, m/2, m/2);
            fmpz_mat_init(p1, m/2, m/2);
            fmpz_mat_init(p2, m/2, m/2);
            fmpz_mat_init(p3, m/2, m/2);
            fmpz_mat_init(p4, m/2, m/2);
            fmpz_mat_init(p5, m/2, m/2);
            fmpz_mat_init(p6, m/2, m/2);
            fmpz_mat_init(p7, n/2, n/2);

            fmpz_mat_window_init(window11, window_A, 0, 0, m/2, m/2);
            fmpz_mat_window_init(window12, window_A, 0, m/2, m/2, m);
            fmpz_mat_window_init(window21, window_A, m/2, 0, m, m/2);
            fmpz_mat_window_init(window22, window_A, m/2, m/2, m, m);
          
            
            fmpz_mat_add(s1, window22, window12);
            fmpz_mat_sub(s2, window22, window21);
            fmpz_mat_add(s3, s2, window12);
            fmpz_mat_sub(s4, s3, window11);
       

            fmpz_mat_sqr_bodrato(p1, s1);
            fmpz_mat_sqr_bodrato(p2, s2);
            fmpz_mat_sqr_bodrato(p3, s3);    
            fmpz_mat_sqr_bodrato(p4, window11);
            fmpz_mat_mul(p5, window12, window21);
            fmpz_mat_mul(p6, s4, window12);
            fmpz_mat_mul(p7, window21, s4);
            

            fmpz_mat_zero(s1);
            fmpz_mat_zero(s2);
            fmpz_mat_zero(s3);
           
            fmpz_mat_add(s1, p3, p5);
            fmpz_mat_sub(s2, p1, s1);
            fmpz_mat_sub(s3, s1, p2);
            
            fmpz_mat_zero(s1);
            fmpz_mat_zero(p1);
            fmpz_mat_zero(p3);
            fmpz_mat_zero(s4);
            
            fmpz_mat_add(p1, p4, p5);
            fmpz_mat_sub(s1, s3, p6);
            fmpz_mat_sub(p3, s2, p7);
            fmpz_mat_add(s4, p2, s2);
            
            for (i = 0; i < m/2; ++i)
            {
                for (j = 0; j < m/2; ++j)
                {
                    fmpz_init(temp);
                    fmpz_init(val);

                    fmpz_mul(val, fmpz_mat_entry(A, i, n-1), fmpz_mat_entry(A, n-1, j));
                    fmpz_add(temp, fmpz_mat_entry(p1, i, j), val);
                    fmpz_set(fmpz_mat_entry(B, i, j), temp);
                    
                    fmpz_clear(temp);
                    fmpz_clear(val);
                }
            }

            for (i = m/2; i < m; ++i)
            {
                for (j = 0; j < m/2; ++j)
                {
                    fmpz_init(temp);
                    fmpz_init(val);

                    fmpz_mul(val, fmpz_mat_entry(A, i, n-1), fmpz_mat_entry(A, n-1, j));
                    fmpz_add(temp, fmpz_mat_entry(p3, i-m/2, j), val);
                    fmpz_set(fmpz_mat_entry(B, i, j), temp);
                    
                    fmpz_clear(temp);
                    fmpz_clear(val);

                }
            }

            for (i = 0; i < m/2; ++i)
            {
                for (j = m/2; j < m; ++j)
                {
                    fmpz_init(temp);
                    fmpz_init(val);

                    fmpz_mul(val, fmpz_mat_entry(A, i, n-1), fmpz_mat_entry(A, n-1, j));
                    fmpz_add(temp, fmpz_mat_entry(s1, i, j-m/2), val);
                    fmpz_set(fmpz_mat_entry(B, i, j), temp);
                    
                    fmpz_clear(temp);
                    fmpz_clear(val);

                }
            }

            for (i = m/2; i < m; ++i)
            {
                for (j = m/2; j < m; ++j)
                {
                    fmpz_init(temp);
                    fmpz_init(val);

                    fmpz_mul(val, fmpz_mat_entry(A, i, n-1), fmpz_mat_entry(A, n-1, j));
                    fmpz_add(temp, fmpz_mat_entry(s4, i-m/2, j-m/2), val);
                    fmpz_set(fmpz_mat_entry(B, i, j), temp);
                    
                    fmpz_clear(temp);
                    fmpz_clear(val);
                }
            }
   
            fmpz_mat_window_clear(window11);
            fmpz_mat_window_clear(window12);
            fmpz_mat_window_clear(window21);
            fmpz_mat_window_clear(window22);
 
            fmpz_mat_clear(s1);
            fmpz_mat_clear(s2);
            fmpz_mat_clear(s3);
            fmpz_mat_clear(s4);
            fmpz_mat_clear(p1);
            fmpz_mat_clear(p2);
            fmpz_mat_clear(p3);
            fmpz_mat_clear(p4);
            fmpz_mat_clear(p5);
            fmpz_mat_clear(p6);
            fmpz_mat_clear(p7);


            /* Matrix Peeling */
            for (i = 0; i < n; ++i)
            {
                fmpz_init(sum);
                for (x = 0; x < n; ++x)
                {
                    fmpz_init(val);
                    fmpz_mul(val, fmpz_mat_entry(A, n-1, x), fmpz_mat_entry(A, x, i));    
                    fmpz_add(sum, sum, val); 
                    fmpz_clear(val);
                }
                fmpz_set(fmpz_mat_entry(B, n-1, i), sum);
                fmpz_clear(sum);
            }

            for (i = 0; i < n; ++i)
            {
                fmpz_init(sum);
                for (x = 0; x < n; ++x)
                {
                    fmpz_init(val);
                    fmpz_mul(val, fmpz_mat_entry(A, x, n-1), fmpz_mat_entry(A, i, x));    
                    fmpz_add(sum, sum, val);
                    fmpz_clear(val);
                }
                fmpz_set(fmpz_mat_entry(B, i, n-1), sum);
                fmpz_clear(sum);
            }

            fmpz_mat_window_clear(window_A);
           
           
            /*
            fmpz_mat_t D;
            fmpz_mat_init(D, n, n);
            fmpz_mat_mul(D, A, A);

            flint_printf("\nA:\n");
            fmpz_mat_print_pretty(A);
            flint_printf("\nB:\n");
            fmpz_mat_print_pretty(B);
            flint_printf("\nC:\n");
            fmpz_mat_print_pretty(D);
            flint_printf("\n");
            fmpz_mat_clear(D);*/

        }
    }
}
