#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <string.h>

typedef struct {
    slong base_nvars;
    slong fibre_nvars;
    slong max_diff_order;
    slong nvars;
    slong nvars_with_derivatives;
    const char** base_vars;
    const char** vars;
    slong* vars_derivatives;
} jet_ctx_struct;

typedef jet_ctx_struct jet_ctx_t[1];

const char *jet_variable_from_multi_index(slong *multi_index, slong num_derivatives, slong base_nvars, const char **base_vars, const char *fibre_var) {
    /* TODO: Handle longer base and fibre variable names? */
    char* var = flint_malloc(sizeof(char) * (2 + num_derivatives + 1));
    var[0] = fibre_var[0];
    var[1] = '_';
    size_t idx = 2;
    for (slong i = 0; i < base_nvars; i++) {
        for (slong j = 0; j < multi_index[i]; j++) {
            var[idx++] = base_vars[i][0];
        }
    }
    var[idx] = '\0';
    return var;
}

void jet_ctx_init(jet_ctx_t jctx, const char** base_vars, slong base_nvars, const char** fibre_vars, slong fibre_nvars, slong max_diff_order)
{
    fmpz_t nvars, tmp, nvars_with_derivatives, num_multi_indices, num_multi_indices_prev, num_vars_prev;
    size_t pos, idx, local_offset;
    slong num_derivatives, fibre_idx, i, rem;
    slong *cur, *cur_integral, *jet_multi_indices, *jet_multi_indices_prev = NULL, *jet_multi_indices_prev_prev = NULL;

    jctx->base_nvars = base_nvars;
    jctx->fibre_nvars = fibre_nvars;
    jctx->max_diff_order = max_diff_order;
    jctx->base_vars = base_vars; /* NOTE: Not copying */

    /* Count the number of jet variables: u, v, u_x, u_y, u_z, u_xx, etc. */
    fmpz_init(nvars);
    fmpz_init(tmp);
    fmpz_zero(nvars);
    /* TODO: Binomial sum closed form? */
    for (slong i = 0; i <= max_diff_order; i++) {
        fmpz_bin_uiui(tmp, (ulong)(base_nvars + i - 1), (ulong)i);
        fmpz_add(nvars, nvars, tmp);
    }
    fmpz_mul_si(nvars, nvars, fibre_nvars);
    fmpz_set(nvars_with_derivatives, nvars);
    fmpz_submul_si(nvars_with_derivatives, tmp, fibre_nvars);
    fmpz_clear(tmp);

    /* Generate jet variable names. */

    jctx->nvars = fmpz_get_ui(nvars);
    jctx->nvars_with_derivatives = fmpz_get_si(nvars_with_derivatives);
    jctx->vars = flint_malloc(sizeof(char*) * (size_t) fmpz_get_ui(nvars));
    jctx->vars_derivatives = flint_malloc(sizeof(slong) * (size_t) (base_nvars * jctx->nvars_with_derivatives));

    /*
     * Generating integer vectors, based on SageMath's integer_vectors_nk_fast_iter
     * https://github.com/sagemath/sage/blob/931cc5e8dc03580c9fecbd0535544a725e9baeb7/src/sage/combinat/integer_vector.py#L1742
     */
    idx = 0;
    cur = flint_malloc(sizeof(slong) * base_nvars);
    cur_integral = flint_malloc(sizeof(slong) * base_nvars);
    jet_multi_indices = flint_calloc(base_nvars, sizeof(slong)); /* NOTE: Multi-indices with all zero entries. */
    fmpz_init(num_multi_indices);
    fmpz_init(num_multi_indices_prev);
    fmpz_init(num_vars_prev);
    fmpz_one(num_multi_indices);
    fmpz_zero(num_vars_prev);
    for (i = 0; i < fibre_nvars; i++) {
        /* TODO: Support longer fibre variable lengths? */
        jctx->vars[idx] = flint_malloc(sizeof(char));
        strncpy((char*)jctx->vars[idx], fibre_vars[i], 1);
        idx++;
    }
    for (num_derivatives = 1; num_derivatives <= max_diff_order; num_derivatives++) {
        jet_multi_indices_prev_prev = jet_multi_indices_prev;
        if (jet_multi_indices_prev_prev != NULL) flint_free(jet_multi_indices_prev_prev);
        jet_multi_indices_prev = jet_multi_indices;
        fmpz_addmul_si(num_vars_prev, num_multi_indices_prev, fibre_nvars);
        fmpz_set(num_multi_indices_prev, num_multi_indices);
        fmpz_bin_uiui(num_multi_indices, (ulong)(base_nvars + num_derivatives - 1), (ulong)num_derivatives);
        jet_multi_indices = flint_malloc(sizeof(slong) * fibre_nvars * (size_t)fmpz_get_ui(num_multi_indices));
        local_offset = 0;
        for (fibre_idx = 0; fibre_idx < fibre_nvars; fibre_idx++) {
            /* Integer vectors of length nvars that sum to num_derivatives */
            cur[0] = num_derivatives + 1; /* NOTE: Artificial + 1 to make the loop work. */
            for (i = 1; i < base_nvars; i++) {
                cur[i] = 0;
            }

            pos = 0;
            rem = -1; /* NOTE: Artificial -1 to make the loop work. */
            while (1) {
                if (cur[pos] == 0) {
                    if (pos == 0) {
                        break;
                    } else {
                        pos--;
                        continue;
                    }
                }
                cur[pos]--;
                rem++;

                int rollover = (pos + 2 == (size_t) base_nvars);
                if (rollover) {
                    cur[pos + 1] = rem;
                } else {
                    pos++;
                    cur[pos] = rem; /* Guaranteed to be at least 1. */
                }

                /* Store the multi-index (for simplicity, only do it for the total derivatives of the first fibre variable) */
                if (fibre_idx == 0) {
                    memcpy(jet_multi_indices + local_offset, cur, sizeof(slong) * base_nvars);
                    local_offset += base_nvars;
                }

                jctx->vars[idx] = jet_variable_from_multi_index(cur, num_derivatives, base_nvars, base_vars, fibre_vars[fibre_idx]);

                /* Find all the "integrals" of the jet variable, if any */
                memcpy(cur_integral, cur, sizeof(slong) * base_nvars);
                for (i = 0; i < base_nvars; i++) {
                    if (cur[i] == 0)
                        continue;
                    int found = 0;
                    cur_integral[i]--;
                    for (size_t j = 0; j < (size_t)fmpz_get_ui(num_multi_indices_prev); j++) {
                        if (memcmp(jet_multi_indices_prev + j*base_nvars, cur_integral, sizeof(slong) * base_nvars) == 0) {
                            /* Store the derivative */
                            size_t j_global = fmpz_get_ui(num_vars_prev) + fibre_idx * fmpz_get_ui(num_multi_indices_prev) + j;
                            jctx->vars_derivatives[i * jctx->nvars_with_derivatives + j_global] = idx;
                            found = 1;
                            break;
                        }
                    }
                    (void)found; /* Silence compiler warning when asserts are disabled. */
                    FLINT_ASSERT(found);
                    cur_integral[i]++;
                }

                idx++;

                if (rollover) {
                    cur[pos + 1] = 0;
                } else {
                    rem = 0;
                }
            }
        }
    }
    flint_free(cur);
    flint_free(cur_integral);
    flint_free(jet_multi_indices);
    flint_free(jet_multi_indices_prev);
    fmpz_clear(num_multi_indices);
    fmpz_clear(num_multi_indices_prev);

    fmpz_clear(nvars);
}

void jet_ctx_clear(jet_ctx_t jctx)
{
    for (slong i = 0; i < jctx->nvars; i++) {
        flint_free((void*)jctx->vars[i]);
    }
    flint_free(jctx->vars);
    flint_free(jctx->vars_derivatives);
}

void _fmpq_mpoly_jet_total_derivative(fmpq_mpoly_t res, const fmpq_mpoly_t f, slong base_var, fmpq_mpoly_ctx_t ctx, jet_ctx_t jctx)
{
    fmpq_mpoly_zero(res, ctx);
    fmpq_mpoly_t tmp;
    fmpq_mpoly_init(tmp, ctx);
    fmpq_mpoly_t v;
    fmpq_mpoly_init(v, ctx);
    for (slong i = 0; i < jctx->nvars_with_derivatives; i++)
    {
        fmpq_mpoly_derivative(tmp, f, i, ctx);
        if (!fmpq_mpoly_is_zero(tmp, ctx))
        {
            fmpq_mpoly_gen(v, jctx->vars_derivatives[base_var * jctx->nvars_with_derivatives + i], ctx);
            fmpq_mpoly_mul(tmp, tmp, v, ctx);
            fmpq_mpoly_add(res, res, tmp, ctx);
        }
    }
    fmpq_mpoly_clear(tmp, ctx);
    fmpq_mpoly_clear(v, ctx);
}

void fmpq_mpoly_jet_total_derivative(fmpq_mpoly_t res, const fmpq_mpoly_t f, slong base_var, fmpq_mpoly_ctx_t ctx, jet_ctx_t jctx)
{
    if (res == f)
    {
        fmpq_mpoly_t tmp;
        fmpq_mpoly_init(tmp, ctx);
        _fmpq_mpoly_jet_total_derivative(tmp, f, base_var, ctx, jctx);
        fmpq_mpoly_swap(tmp, res, ctx);
        fmpq_mpoly_clear(tmp, ctx);
    }
    else
    {
        _fmpq_mpoly_jet_total_derivative(res, f, base_var, ctx, jctx);
    }
}

void print_total_derivative(const fmpq_mpoly_t f, slong base_var, fmpq_mpoly_ctx_t ctx, jet_ctx_t jctx)
{
    fmpq_mpoly_t df;
    fmpq_mpoly_init(df, ctx);
    flint_printf("d/d%s(", jctx->base_vars[base_var]);
    fmpq_mpoly_print_pretty(f, jctx->vars, ctx);
    flint_printf(") = ");
    fmpq_mpoly_jet_total_derivative(df, f, base_var, ctx, jctx);
    fmpq_mpoly_print_pretty(df, jctx->vars, ctx);
    fmpq_mpoly_clear(df, ctx);
}

int main(int argc, char* argv[])
{
    slong base_nvars = 3;
    slong fibre_nvars = 2;
    const char* base_vars[] = {"x", "y", "z"};
    const char* fibre_vars[] = {"u", "v"};
    slong max_diff_order = 3;

    ordering_t ord = ORD_DEGREVLEX;
    fmpq_mpoly_ctx_t ctx;
    jet_ctx_t jctx;
    fmpq_mpoly_t v, f;

    jet_ctx_init(jctx, base_vars, base_nvars, fibre_vars, fibre_nvars, max_diff_order);

    fmpq_mpoly_ctx_init(ctx, jctx->nvars, ord);
    fmpq_mpoly_init(v, ctx);
    fmpq_mpoly_init(f, ctx);

    for (slong base_var = 0; base_var < jctx->base_nvars; base_var++)
    {
        for (slong jet_var = 0; jet_var < jctx->nvars_with_derivatives; jet_var++) {
            fmpq_mpoly_gen(v, jet_var, ctx);
            print_total_derivative(v, base_var, ctx, jctx); flint_printf("; ");
        }
        flint_printf("\n");
    }
    flint_printf("\n");

    fmpq_mpoly_set_str_pretty(f, "u^2", jctx->vars, ctx);
    print_total_derivative(f, 0, ctx, jctx); flint_printf("\n");
    fmpq_mpoly_set_str_pretty(f, "u^3", jctx->vars, ctx);
    print_total_derivative(f, 1, ctx, jctx); flint_printf("\n");
    fmpq_mpoly_set_str_pretty(f, "u*v", jctx->vars, ctx);
    print_total_derivative(f, 2, ctx, jctx); flint_printf("\n");
    fmpq_mpoly_set_str_pretty(f, "u_x*v_y", jctx->vars, ctx);
    print_total_derivative(f, 2, ctx, jctx); flint_printf("\n");

    fmpq_mpoly_clear(v, ctx);
    fmpq_mpoly_clear(f, ctx);
    fmpq_mpoly_ctx_clear(ctx);

    jet_ctx_clear(jctx);

    flint_cleanup_master();

    return 0;
}
