#ifndef GR_MAT_H
#define GR_MAT_H

#ifdef GR_MAT_INLINES_C
#define GR_MAT_INLINE FLINT_DLL
#else
#define GR_MAT_INLINE static __inline__
#endif

#include "flint/fmpz_mat.h"
#include "gr.h"
#include "gr_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    gr_ptr entries;
    slong r;
    slong c;
    gr_ptr * rows;
}
gr_mat_struct;

typedef gr_mat_struct gr_mat_t[1];

#define GR_MAT_ENTRY(mat,i,j,sz) GR_ENTRY((mat)->rows[i], j, sz)
#define gr_mat_nrows(mat, ctx) ((mat)->r)
#define gr_mat_ncols(mat, ctx) ((mat)->c)

int gr_mat_init(gr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx);
int gr_mat_clear(gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_swap(gr_mat_t mat1, gr_mat_t mat2, gr_ctx_t ctx);
int gr_mat_randtest(gr_mat_t mat, flint_rand_t state, void * options, gr_ctx_t ctx);
truth_t gr_mat_equal(const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
int gr_mat_zero(gr_mat_t res, gr_ctx_t ctx);
truth_t gr_mat_is_zero(const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_one(const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_set_si(gr_mat_t res, slong v, gr_ctx_t ctx);
int gr_mat_one(gr_mat_t res, gr_ctx_t ctx);
int gr_mat_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_set_ui(gr_mat_t res, ulong v, gr_ctx_t ctx);
int gr_mat_set_si(gr_mat_t res, slong v, gr_ctx_t ctx);
int gr_mat_set_fmpz(gr_mat_t res, const fmpz_t v, gr_ctx_t ctx);
int gr_mat_set_fmpq(gr_mat_t res, const fmpq_t v, gr_ctx_t ctx);
int gr_mat_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_swap_entrywise(gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
int gr_mat_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
int gr_mat_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
int gr_mat_print(const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_mul_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
