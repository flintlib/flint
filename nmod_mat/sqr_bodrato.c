#include "nmod_mat.h"


void
nmod_mat_sqr_bodrato(nmod_mat_t B, const nmod_mat_t A)
{
    slong n = A->r;
    slong i,j;
    
    nmod_t mod;

    mod = A->mod;

    if (n == 0)
    {
        return;
    }
    else if (n == 1)
    {
        nmod_mat_entry(B, 0, 0) = nmod_mul(nmod_mat_entry(A, 0, 0), nmod_mat_entry(A, 0, 0), mod);
    }
    else if (n == 2)
    {
        mp_limb_t t, u;

        t = nmod_add(nmod_mat_entry(A, 0, 0), nmod_mat_entry(A, 1, 1), mod);
        u = nmod_mul(nmod_mat_entry(A, 0, 1), nmod_mat_entry(A, 1, 0), mod);

        nmod_mat_entry(B, 0, 0) = nmod_mul(nmod_mat_entry(A, 0, 0), nmod_mat_entry(A, 0, 0), mod);
        nmod_mat_entry(B, 0, 0) = nmod_add(nmod_mat_entry(B, 0, 0), u, mod);

        nmod_mat_entry(B, 1, 1) = nmod_mul(nmod_mat_entry(A, 1, 1), nmod_mat_entry(A, 1, 1), mod);
        nmod_mat_entry(B, 1, 1) = nmod_add(nmod_mat_entry(B, 1, 1), u, mod);

        nmod_mat_entry(B, 0, 1) = nmod_mul(nmod_mat_entry(A, 0, 1), t, mod);
        nmod_mat_entry(B, 1, 0) = nmod_mul(nmod_mat_entry(A, 1, 0), t, mod);

    }
    else
    {
        nmod_mat_mul(B, A, A);
    
        nmod_mat_t window_A, window11, window12, window21, window22;
        nmod_mat_t s1, s2, s3, s4;
        nmod_mat_t p1, p2, p3, p4, p5, p6, p7;
        
        mp_limb_t sum, val;

        slong m = n, x, iseven = 1; 
        
        if (n % 2 == 1)
        {
            m = n - 1;
            iseven = 0;
        }

        nmod_mat_init(s1, m/2, m/2, A->mod.n);
        nmod_mat_init(s2, m/2, m/2, A->mod.n);
        nmod_mat_init(s3, m/2, m/2, A->mod.n);
        nmod_mat_init(s4, m/2, m/2, A->mod.n);
        nmod_mat_init(p1, m/2, m/2, A->mod.n);
        nmod_mat_init(p2, m/2, m/2, A->mod.n);
        nmod_mat_init(p3, m/2, m/2, A->mod.n);
        nmod_mat_init(p4, m/2, m/2, A->mod.n);
        nmod_mat_init(p5, m/2, m/2, A->mod.n);
        nmod_mat_init(p6, m/2, m/2, A->mod.n);
        nmod_mat_init(p7, m/2, m/2, A->mod.n);

        nmod_mat_window_init(window_A, A, 0, 0, m, m);
        nmod_mat_window_init(window11, window_A, 0, 0, m/2, m/2);
        nmod_mat_window_init(window12, window_A, 0, m/2, m/2, m);
        nmod_mat_window_init(window21, window_A, m/2, 0, m, m/2);
        nmod_mat_window_init(window22, window_A, m/2, m/2, m, m);

        nmod_mat_add(s1, window22, window12);
        nmod_mat_sub(s2, window22, window21);
        nmod_mat_add(s3, s2, window12);
        nmod_mat_sub(s4, s3, window11);

        nmod_mat_sqr(p1, s1);
        nmod_mat_sqr(p2, s2);
        nmod_mat_sqr(p3, s3);    
        nmod_mat_sqr(p4, window11);
        nmod_mat_mul(p5, window12, window21);
        nmod_mat_mul(p6, s4, window12);
        nmod_mat_mul(p7, window21, s4);

        nmod_mat_clear(s1);
        nmod_mat_clear(s2);
        nmod_mat_clear(s3);

        nmod_mat_init(s1, m/2, m/2, A->mod.n);
        nmod_mat_init(s2, m/2, m/2, A->mod.n);
        nmod_mat_init(s3, m/2, m/2, A->mod.n);
    

        nmod_mat_add(s1, p3, p5);
        nmod_mat_sub(s2, p1, s1);
        nmod_mat_sub(s3, s1, p2);

        nmod_mat_clear(s1);
        nmod_mat_clear(p1);
        nmod_mat_clear(p3);
        nmod_mat_clear(s4);

        nmod_mat_init(s1, m/2, m/2, A->mod.n);
        nmod_mat_init(p1, m/2, m/2, A->mod.n);
        nmod_mat_init(p3, m/2, m/2, A->mod.n);
        nmod_mat_init(s4, m/2, m/2, A->mod.n);
    

        nmod_mat_add(p1, p4, p5);
        nmod_mat_sub(s1, s3, p6);
        nmod_mat_sub(p3, s2, p7);
        nmod_mat_add(s4, p2, s2);

        if (iseven == 1)
        {
            for (i = 0; i < n/2; ++i)
            {
                for (j = 0; j < n/2; ++j)
                {
                    nmod_mat_entry(B, i, j) = nmod_mat_entry(p1, i, j);
                }
            }

            for (i = n/2; i < n; ++i)
            {
                for (j = 0; j < n/2; ++j)
                {
                    nmod_mat_entry(B, i, j) = nmod_mat_entry(p3, i - n/2, j);
                }
            }

            for (i = 0; i < n/2; ++i)
            {
                for (j = n/2; j < n; ++j)
                {
                    nmod_mat_entry(B, i, j) = nmod_mat_entry(s1, i, j - n/2);
                }
            }

            for (i = n/2; i < n; ++i)
            {
                for (j = n/2; j < n; ++j)
                {
                    nmod_mat_entry(B, i, j) = nmod_mat_entry(s4, i - n/2, j - n/2);
                }
            }
        }
        else
        {
            nmod_mat_t temp_A, cache_A;
            
            nmod_mat_init(temp_A, n, n, A->mod.n);
            nmod_mat_init(cache_A, n, n, A->mod.n);

            
            nmod_mat_set(temp_A, A);


            for (i = 0; i < n; ++i)
            {
                for (j = 0; j < n; ++j)
                {
                    nmod_mat_entry(cache_A, i, j) = nmod_mul(nmod_mat_entry(A, i, n - 1), nmod_mat_entry(A, n - 1, j), mod); 
                }
            }

            for (i = 0; i < n; ++i)
            {
                sum = 0;
                for (x = 0; x < n; ++x)
                {
                    val = nmod_mul(nmod_mat_entry(temp_A, n - 1, x), nmod_mat_entry(temp_A, x, i), mod);
                    sum = nmod_add(sum, val, mod);
                }
                nmod_mat_entry(B, n - 1, i) = sum;
            }

            for (i = 0; i < n; ++i)
            {
                sum = 0;
                for (x = 0; x < n; ++x)
                {
                    val = nmod_mul(nmod_mat_entry(temp_A, x, n - 1), nmod_mat_entry(temp_A, i, x), mod);
                    sum = nmod_add(sum, val, mod);
                }
                nmod_mat_entry(B, i, n - 1) = sum;
            }

            for (i = 0; i < m/2; ++i)
            {
                for (j = 0; j < m/2; ++j)
                {
                    nmod_mat_entry(B, i, j) = nmod_add(nmod_mat_entry(p1, i, j), nmod_mat_entry(cache_A, i, j), mod);
                }
            }
            for (i = m/2; i < m; ++i)
            {
                for (j = 0; j < m/2; ++j)
                {
                    nmod_mat_entry(B, i, j) = nmod_add(nmod_mat_entry(p3, i - m/2, j), nmod_mat_entry(cache_A, i, j), mod);
                }
            }
            for (i = 0; i < m/2; ++i)
            {
                for (j = m/2; j < m; ++j)
                {
                    nmod_mat_entry(B, i, j) = nmod_add(nmod_mat_entry(s1, i, j - m/2), nmod_mat_entry(cache_A, i, j), mod);
                }
            }
            for (i = m/2; i < m; ++i)
            {
                for (j = m/2; j < m; ++j)
                {
                    nmod_mat_entry(B, i, j) = nmod_add(nmod_mat_entry(s4, i - m/2, j - m/2), nmod_mat_entry(cache_A, i, j), mod);
                }
            }

            nmod_mat_clear(temp_A);
            nmod_mat_clear(cache_A);

        }

        nmod_mat_window_clear(window_A);
        nmod_mat_window_clear(window11);
        nmod_mat_window_clear(window12);
        nmod_mat_window_clear(window21);
        nmod_mat_window_clear(window22);
        nmod_mat_clear(s1);
        nmod_mat_clear(s2);
        nmod_mat_clear(s3);
        nmod_mat_clear(s4);
        nmod_mat_clear(p1);
        nmod_mat_clear(p2);
        nmod_mat_clear(p3);
        nmod_mat_clear(p4);
        nmod_mat_clear(p5);
        nmod_mat_clear(p6);
        nmod_mat_clear(p7);
        
    }
}
