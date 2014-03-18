#include "fmpz_mat.h"

void fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong i,i2,j,k,l;
    fmpz_mat_init_set(H,A);

    for (i = A->r - 1, k = A->c - 1, l = (i - k)*(i > k); i + 1 != l; i--, k--)
    {
        int row_finished = 1;
        for (j = 0; (j < k) && row_finished; j++)
        {
            row_finished = fmpz_is_zero(fmpz_mat_entry(H,i,j));
        }
        if (row_finished)
        {
            if (fmpz_sgn(fmpz_mat_entry(H,i,k)) < 0)
            {
                for (i2 = 0; i2 < A->r; i2++)
                {
                    fmpz_neg(fmpz_mat_entry(H,i2,k),fmpz_mat_entry(H,i2,k));
                }
            }
            if (fmpz_is_zero(fmpz_mat_entry(H,i,k)))
            {
                k++;
            }
            else
            {
                for (j = k + 1; j < A->c; j++)
                {
                    slong q;
                    fmpz_fdiv_q(&q,fmpz_mat_entry(H,i,j),fmpz_mat_entry(H,i,k));
                    for (i2 = 0; i2 < A->r; i2++)
                    {
                        fmpz_submul(fmpz_mat_entry(H,i2,j), &q,
                                fmpz_mat_entry(H,i2,k));
                    }
                }
            }
        }
        else
        {
            slong j0,min;
            min = 0;
            for (j = 0; j <= k; j++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(H,i,j)))
                    continue;
                if (fmpz_is_zero(&min) || 
                        fmpz_cmpabs(fmpz_mat_entry(H,i,j), &min) < 0)
                {
                    j0 = j;
                    fmpz_abs(&min,fmpz_mat_entry(H,i,j));
                }
            }
            if (fmpz_cmp(&j0, &k) < 0)
            {
                for (i2 = 0; i2 < A->r; i2++)
                {
                    fmpz_swap(fmpz_mat_entry(H,i2,j0), fmpz_mat_entry(H,i2,k));
                }
            }
            if (fmpz_sgn(fmpz_mat_entry(H,i,k)) < 0)
            {
                for (i2 = 0; i2 < A->r; i2++)
                {
                    fmpz_neg(fmpz_mat_entry(H,i2,k), fmpz_mat_entry(H,i2,k));
                }
            }
            for (j = 0; j <= k - 1; j++)
            {
                slong q;
                fmpz_fdiv_q(&q,fmpz_mat_entry(H,i,j), fmpz_mat_entry(H,i,k));
                for (i2 = 0; i2 < A->r; i2++)
                {
                    fmpz_submul(fmpz_mat_entry(H,i2,j), &q,
                            fmpz_mat_entry(H,i2,k));
                }
            }
            i++;
            k++;
        }
    }
    /* move output to left */
    /* for (j = 0; j < A->c - k - 1; j++)
    {
        for (i = 0; i < A->r; i++)
        {
            fmpz_set(fmpz_mat_entry(H,i,j), fmpz_mat_entry(H,i,j + k + 1));
        }
    }
    for (; j < A->c; j++)
    {
        for (i = 0; i < A->r; i++)
        {
            fmpz_zero(fmpz_mat_entry(H,i,j));
        }
    } */
}
