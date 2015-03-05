#include "fmpz_mat.h"

void
fmpz_mat_concat_vertical(fmpz_mat_t mat1, fmpz_mat_t mat2)
{
    slong i,j;
    slong r1 = mat1->r;
    slong c1 = mat1->c;
    slong r2 = mat2->r;
    slong c2 = mat2->c;
    
    if(c1!=c2)
        return;

    for(i=0;i<r2;++i)
    {
        for(j=0;j<c2;++j)
        {
            fmpz_set(fmpz_mat_entry(mat1, (i+r1), j), fmpz_mat_entry(mat2, i, j));
        }
    }
}

void
fmpz_mat_concat_horizontal(fmpz_mat_t mat1, fmpz_mat_t mat2)
{
    slong i,j;
    slong r1 = mat1->r;
    slong c1 = mat1->c;
    slong r2 = mat2->r;
    slong c2 = mat2->c;
    
    if(r1!=r2)
        return;

    for(i=0;i<r2;++i)
    {
        for(j=0;j<c2;++j)
        {
            fmpz_set(fmpz_mat_entry(mat1, i, (j+c1)), fmpz_mat_entry(mat2, i, j));
        }
    }
}
