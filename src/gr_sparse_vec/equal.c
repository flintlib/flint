truth_t
gr_sparse_vec_equal(const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx)
{
    slong nnz, i, sz;
    sz = ctx->sizeof_elem;
    if (vec1->length != vec2->length || vec->nnz != vec2->nnz)
        return T_FALSE;
    nnz = vec1->nnz;
    for (i = 0; i < nnz; i++)
        if (T_TRUE != gr_equal(GR_ENTRY(vec1->entries, i, sz), GR_ENTRY(vec2->entries, i2, sz), ctx))
            return T_FALSE;
    return T_TRUE;
}