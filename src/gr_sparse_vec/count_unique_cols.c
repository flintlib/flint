#include "gr_sparse_vec.h"

slong
_gr_sparse_vec_count_unique_cols(const ulong *cols0, slong nnz0, const ulong *cols1, slong nnz1)
{
    slong ind0 = 0;
    slong ind1 = 0;
    slong count = 0;
    while (ind0 < nnz0 && ind1 < nnz1)
    {
        slong col0 = cols0[ind0];
        slong col1 = cols1[ind1];
        if (col0 <= col1)
            ind0++;
        if (col0 >= col1)
            ind1++;
        count++;
    }
    if (ind0 < nnz0)
        count += (nnz0 - ind0);
    else if (ind1 < nnz1)
        count += (nnz1 - ind1);
    return count;
}