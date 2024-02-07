#include "gr_sparse_vec.h"

slong _gr_sparse_vec_count_unique_inds(const ulong *inds0, slong nnz0, const ulong *inds1, slong nnz1)
{
    slong ptr0 = 0;
    slong ptr1 = 0;
    slong count = 0;
    while (ptr0 < nnz0 && ptr1 < nnz1)
    {
        slong ind0 = inds0[ptr0];
        slong ind1 = inds1[ptr1];
        if (ind0 <= ind1)
            ptr0++;
        if (ind0 >= ind1)
            ptr1++;
        count++;
    }
    if (ptr0 < nnz0)
        count += (nnz0 - ptr0);
    else if (ptr1 < nnz1)
        count += (nnz1 - ptr1);
    return count;
}