#ifndef GR_MAT_H
#define GR_MAT_H

#ifdef GR_MAT_INLINES_C
#define GR_MAT_INLINE FLINT_DLL
#else
#define GR_MAT_INLINE static __inline__
#endif

#include "flint/fmpz_mat.h"
#include "flint/fmpq_mat.h"
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

void gr_mat_init(gr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx);
void gr_mat_clear(gr_mat_t mat, gr_ctx_t ctx);

GR_MAT_INLINE void
gr_mat_swap(gr_mat_t mat1, gr_mat_t mat2, gr_ctx_t ctx)
{
    if (mat1 != mat2)
    {
        gr_mat_t tmp;
        *tmp = *mat1;
        *mat1 = *mat2;
        *mat2 = *tmp;
    }
}

void gr_mat_window_init(gr_mat_t window, const gr_mat_t mat, slong r1, slong c1, slong r2, slong c2, gr_ctx_t ctx);

GR_MAT_INLINE void
gr_mat_window_clear(gr_mat_t window, gr_ctx_t ctx)
{
    flint_free(window->rows);
}

int gr_mat_write(gr_stream_t out, const gr_mat_t mat, gr_ctx_t ctx);

GR_INLINE int
gr_mat_print(const gr_mat_t mat, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_mat_write(out, mat, ctx);
}

int gr_mat_randtest(gr_mat_t mat, flint_rand_t state, gr_ctx_t ctx);
int gr_mat_randops(gr_mat_t mat, flint_rand_t state, slong count, gr_ctx_t ctx);

GR_INLINE truth_t
gr_mat_is_empty(const gr_mat_t mat, gr_ctx_t ctx)
{
    return ((mat->r == 0) || (mat->c == 0)) ? T_TRUE : T_FALSE;
}

GR_INLINE truth_t
gr_mat_is_square(const gr_mat_t mat, gr_ctx_t ctx)
{
    return (mat->r == mat->c) ? T_TRUE : T_FALSE;
}

truth_t gr_mat_equal(const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
truth_t gr_mat_is_zero(const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_one(const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx);

int gr_mat_zero(gr_mat_t res, gr_ctx_t ctx);
int gr_mat_one(gr_mat_t res, gr_ctx_t ctx);
int gr_mat_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_set_scalar(gr_mat_t res, gr_srcptr c, gr_ctx_t ctx);
int gr_mat_set_ui(gr_mat_t res, ulong v, gr_ctx_t ctx);
int gr_mat_set_si(gr_mat_t res, slong v, gr_ctx_t ctx);
int gr_mat_set_fmpz(gr_mat_t res, const fmpz_t v, gr_ctx_t ctx);
int gr_mat_set_fmpq(gr_mat_t res, const fmpq_t v, gr_ctx_t ctx);

int gr_mat_set_fmpz_mat(gr_mat_t res, const fmpz_mat_t mat, gr_ctx_t ctx);
int gr_mat_set_fmpq_mat(gr_mat_t res, const fmpq_mat_t mat, gr_ctx_t ctx);

int gr_mat_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_swap_entrywise(gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
int gr_mat_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);
int gr_mat_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx);

/* todo: test, wrap; div; more conversions */
int gr_mat_add_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
int gr_mat_sub_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
int gr_mat_mul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
int gr_mat_addmul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
int gr_mat_submul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);
int gr_mat_div_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx);

int gr_mat_mul_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int gr_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

/* todo */
GR_INLINE int
gr_mat_sqr(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_mul(res, mat, mat, ctx);
}

int gr_mat_find_nonzero_pivot(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx);

int gr_mat_lu_recursive(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);
int gr_mat_lu_classical(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);
int gr_mat_lu(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);

int gr_mat_fflu(slong * res_rank, slong * P, gr_mat_t LU, gr_ptr den, const gr_mat_t A, int rank_check, gr_ctx_t ctx);

int gr_mat_solve_fflu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int gr_mat_solve_lu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int gr_mat_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

int gr_mat_solve_fflu_precomp(gr_mat_t X, const slong * perm, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int gr_mat_solve_lu_precomp(gr_mat_t X, const slong * perm, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

int gr_mat_solve_den_fflu(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int gr_mat_solve_den(gr_mat_t X, gr_ptr den, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

int gr_mat_det_bareiss(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
int gr_mat_det_berkowitz(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
int gr_mat_det_lu(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
int gr_mat_det_cofactor(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);
int gr_mat_det(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx);

int gr_mat_rank_lu(slong * rank, const gr_mat_t A, gr_ctx_t ctx);
int gr_mat_rank_fflu(slong * rank, const gr_mat_t A, gr_ctx_t ctx);
int gr_mat_rank(slong * rank, const gr_mat_t A, gr_ctx_t ctx);

int gr_mat_ones(gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_pascal(gr_mat_t mat, int triangular, gr_ctx_t ctx);
int gr_mat_stirling(gr_mat_t mat, int kind, gr_ctx_t ctx);
int gr_mat_hilbert(gr_mat_t mat, gr_ctx_t ctx);
/* todo: hadamard, dft, dct */

int gr_mat_transpose(gr_mat_t B, const gr_mat_t A, gr_ctx_t ctx);

int gr_mat_solve_tril_classical(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);
int gr_mat_solve_tril_recursive(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);
int gr_mat_solve_tril(gr_mat_t X, const gr_mat_t L, const gr_mat_t B, int unit, gr_ctx_t ctx);

int gr_mat_solve_triu_classical(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx);
int gr_mat_solve_triu_recursive(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx);
int gr_mat_solve_triu(gr_mat_t X, const gr_mat_t U, const gr_mat_t B, int unit, gr_ctx_t ctx);

int gr_mat_trace(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);

int _gr_mat_charpoly_berkowitz(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_charpoly_berkowitz(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx);

int _gr_mat_charpoly_danilevsky_inplace(gr_ptr res, gr_mat_t mat, gr_ctx_t ctx);
int _gr_mat_charpoly_danilevsky(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_charpoly_danilevsky(gr_poly_t res, const gr_mat_t mat, gr_ctx_t ctx);

int _gr_mat_charpoly_faddeev(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_charpoly_faddeev(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx);

int _gr_mat_charpoly_faddeev_bsgs(gr_ptr res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_charpoly_faddeev_bsgs(gr_poly_t res, gr_mat_t adj, const gr_mat_t mat, gr_ctx_t ctx);

int _gr_mat_charpoly_hessenberg(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_charpoly_hessenberg(gr_poly_t cp, const gr_mat_t mat, gr_ctx_t ctx);

int gr_mat_hessenberg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_hessenberg_gauss(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
int gr_mat_hessenberg_householder(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx);
truth_t gr_mat_is_hessenberg(const gr_mat_t mat, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
