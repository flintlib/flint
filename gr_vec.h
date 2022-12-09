#ifndef GR_VEC_H
#define GR_VEC_H

#ifdef GR_VEC_INLINES_C
#define GR_VEC_INLINE FLINT_DLL
#else
#define GR_VEC_INLINE static __inline__
#endif

#include "gr.h"

#ifdef __cplusplus
 extern "C" {
#endif

void gr_vec_init(gr_vec_t vec, slong len, gr_ctx_t ctx);
void gr_vec_clear(gr_vec_t vec, gr_ctx_t ctx);

GR_VEC_INLINE gr_ptr
gr_vec_entry_ptr(gr_vec_t vec, slong i, gr_ctx_t ctx)
{
    return GR_ENTRY(vec->entries, i, ctx->sizeof_elem);
}

GR_VEC_INLINE slong gr_vec_length(const gr_vec_t vec, gr_ctx_t ctx)
{
    return vec->length;
}

void gr_vec_fit_length(gr_vec_t vec, slong len, gr_ctx_t ctx);
void gr_vec_set_length(gr_vec_t vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_set(gr_vec_t res, const gr_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_append(gr_vec_t vec, gr_srcptr f, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
