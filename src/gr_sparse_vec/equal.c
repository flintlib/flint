#include "gr_sparse_vec.h"

truth_t
gr_sparse_vec_equal(const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx)
{
    slong i1, i2, sz;
    truth_t cur_test;
    truth_t ret = T_TRUE;
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);

    if (vec1->length != vec2->length)
        return T_FALSE;
    sz = ctx->sizeof_elem;

    for (i1 = 0, i2 = 0; i1 < vec1->nnz && i2 < vec2->nnz; )
    {
        if (vec1->inds[i1] < vec2->inds[i2])
        {
            // In vector => either known or maybe nonzero
            if (is_zero(GR_ENTRY(vec1->nzs, i1, sz), ctx) == T_FALSE)
                return T_FALSE;
            else
                ret = T_UNKNOWN; // Have maybe zero vs known zero
            i1++;
        }
        else if (vec1->inds[i1] > vec2->inds[i2])
        {
            if (is_zero(GR_ENTRY(vec2->nzs, i2, sz), ctx) == T_FALSE)
                return T_FALSE;
            else
                ret = T_UNKNOWN; // Have maybe zero vs known zero
            i2++;
        }
        else
        {
            cur_test = gr_equal(GR_ENTRY(vec1->nzs, i1, sz), GR_ENTRY(vec2->nzs, i2, sz), ctx);
            if (cur_test == T_FALSE)
                return T_FALSE;
            else if (cur_test == T_UNKNOWN)
                ret = T_UNKNOWN;
            i1++; i2++;
        }
    }
    
    return ret;
}
