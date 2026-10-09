/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>
#include "acb.h"
#include "fmpz_factor.h"
#include "ulong_extras.h"
#include "fmpz_poly.h"
#include "fmpz_mpoly_q.h"
#include "gr_poly.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* the options: name, default, range (lo <= value <= hi), with
   GR_TOWER_OPT_FLAGS_MAX standing for "any bit combination" */
typedef struct
{
    const char * name;
    slong value;
    slong lo, hi;
}
option_info_struct;

#define OPT_ANY WORD_MAX

static const option_info_struct option_info[GR_TOWER_OPT_NUM_OPTIONS] =
{
    {"verbose", 0, 0, 2},
    {"print_flags", GR_TOWER_PRINT_SYMBOLIC | GR_TOWER_PRINT_DEFS, 0, GR_TOWER_PRINT_NUMERIC | GR_TOWER_PRINT_SYMBOLIC | GR_TOWER_PRINT_DEFS},
    {"print_digits", GR_TOWER_PRINT_DIGITS_DEFAULT, 1, OPT_ANY},
    {"prec_limit", 256, 64, OPT_ANY},
    {"certify_prec_limit", 4096, 64, OPT_ANY},
    {"numeric_prec_limit", 65536, 64, OPT_ANY},
    {"smooth_limit", 3512, 1, OPT_ANY},
    {"express_degree_limit", 64, 0, OPT_ANY},
    {"express_prec", 256, 64, OPT_ANY},
    {"trager_degree_limit", 48, 0, OPT_ANY},
    {"factor_degree_limit", 128, 0, OPT_ANY},
    {"roots_factor_degree_limit", 64, 0, OPT_ANY},
    {"modular_tries", 6, 0, OPT_ANY},
    {"modular_terms_limit", 50000, 0, OPT_ANY},
    {"no_roots_tries", 24, 0, OPT_ANY},
    {"rationalize_limit", 2000, 0, OPT_ANY},
    {"relation_cost_limit", 200000, 0, OPT_ANY},
    {"sqrt_budget", 2000, 0, OPT_ANY},
    {"gauss_sum_limit", 10000, 0, 1000000},
    {"dense_limit", 100000, 0, OPT_ANY},
    {"inv_dense_degree_limit", 512, 0, OPT_ANY},
    {"inv_linear_degree_limit", 64, 0, OPT_ANY},
    {"deferred_ideal_size", 12, 0, OPT_ANY},
    {"grow_threshold", 6, 0, OPT_ANY},
    {"root_of_unity_order_limit", 360, 1, OPT_ANY},
    {"cyclotomic_order_limit", 100000, 1, OPT_ANY},
    {"realify_order_limit", 480, 1, OPT_ANY},
    {"trig_pi_limit", 240, 1, OPT_ANY},
    {"trig_algebraic_limit", 1000, 1, OPT_ANY},
    {"trig_form", GR_TOWER_TRIG_EXPONENTIAL, 0, 1},
    {"split_imaginary", 0, 0, 1},
    {"composite_radicals", 0, 0, 1},
    {"cyclotomic_degree_limit", 0, 0, OPT_ANY},
    {"same_origin_degree_limit", 120, 0, OPT_ANY},
    {"annihilating_flat_limit", 1024, 0, OPT_ANY},
    {"gamma_lattice_limit", 36, 0, OPT_ANY},
    {"gamma_line_limit", 12, 0, OPT_ANY},
    {"hurwitz_lattice_limit", 240, 0, OPT_ANY},
    {"hurwitz_weight_limit", 64, 0, OPT_ANY},
    {"special_relation_level_limit", 240, 0, OPT_ANY},
    {"gauss_digamma_limit", 30, 0, OPT_ANY},
    {"dense_form_degree_limit", 1048576, 0, OPT_ANY},
    {"dense_form_sparsity", 8, 0, OPT_ANY},
    {"primitive_degree_limit", 0, 0, OPT_ANY},
    {"split_degree_limit", 4096, 0, OPT_ANY},
    {"minpoly_degree_limit", 512, 0, OPT_ANY},
    {"inv_dense_alg", 0, 0, 2},
    {"power_check_degree_limit", 64, 0, OPT_ANY},
};

const slong gr_tower_default_options[GR_TOWER_OPT_NUM_OPTIONS] =
{
    0, GR_TOWER_PRINT_SYMBOLIC | GR_TOWER_PRINT_DEFS, GR_TOWER_PRINT_DIGITS_DEFAULT,
    256, 4096, 65536, 3512, 64, 256, 48, 128, 64, 6, 50000, 24, 2000, 200000, 2000, 10000,
    100000, 512, 64, 12, 6, 360, 100000, 480, 240, 1000, GR_TOWER_TRIG_EXPONENTIAL, 0, 0, 0,
    120, 1024, 36, 12, 240, 64, 240, 30, 1048576, 8, 0, 4096, 512, 0, 64,
};

const char *
gr_tower_option_name(slong option)
{
    if (option < 0 || option >= GR_TOWER_OPT_NUM_OPTIONS)
        return NULL;
    return option_info[option].name;
}

slong
gr_tower_option_find(const char * name)
{
    slong k;
    for (k = 0; k < GR_TOWER_OPT_NUM_OPTIONS; k++)
        if (strcmp(option_info[k].name, name) == 0)
            return k;
    return -1;
}

slong
gr_tower_option_default(slong option)
{
    return option_info[option].value;
}

int
gr_tower_option_valid(slong option, slong value)
{
    if (option < 0 || option >= GR_TOWER_OPT_NUM_OPTIONS)
        return 0;
    return value >= option_info[option].lo && value <= option_info[option].hi;
}

void
gr_tower_options_init(slong * options)
{
    slong k;
    for (k = 0; k < GR_TOWER_OPT_NUM_OPTIONS; k++)
        options[k] = option_info[k].value;
}

int
gr_tower_options_set(slong * options, slong option, slong value)
{
    if (!gr_tower_option_valid(option, value))
        return GR_DOMAIN;
    options[option] = value;
    return GR_SUCCESS;
}

void
gr_tower_set_options(gr_tower_t T, const slong * options)
{
    T->options = (options != NULL) ? options : gr_tower_default_options;
}

void
gr_tower_init(gr_tower_t T, gr_ctx_t base)
{
    /* (the capabilities of the constant field: only QQ is implemented) */
    if (base->which_ring != GR_CTX_FMPQ)
        flint_throw(FLINT_ERROR, "(%s): towers over this field are not implemented (the constant field must be QQ)\n", __func__);
    T->caps = GR_TOWER_CAP_RATIONAL | GR_TOWER_CAP_FACTOR | GR_TOWER_CAP_PLACES |
              GR_TOWER_CAP_EMBEDDING | GR_TOWER_CAP_INVOLUTION | GR_TOWER_CAP_ORDERED;
    T->options = gr_tower_default_options;
    T->consts = base;
    T->base = base;
    T->num_gens = 0;
    T->alloc_gens = 0;
    T->gens = NULL;
    T->length = 0;
    T->alg = NULL;
    T->num_trans = 0;
    T->trans = NULL;
    T->next_gid = 0;
    T->version = 0;
    T->moduli_version = 0;
    T->structure_version = 0;
    T->keep_retired = 1;
    T->retired = NULL;
    T->num_retired = 0;
    gr_tower_flat_init(&T->flat, T, 8);
}

slong
gr_tower_gid_order(const gr_tower_t T, slong gid)
{
    slong d;

    /* usually the identity */
    if (gid >= 0 && gid < T->num_gens && T->gens[gid].gid == gid)
        return gid;

    for (d = 0; d < T->num_gens; d++)
        if (T->gens[d].gid == gid)
            return d;

    return -1;
}

static void
_gr_tower_gen_clear(gr_tower_gen_struct * g)
{
    if (g->ctx != NULL)
    {
        gr_ctx_clear(g->ctx);
        flint_free(g->ctx);
    }
    if (g->arg.mctx != NULL)
        fmpz_mpoly_q_clear(&g->arg.data, g->arg.mctx);
    _gr_tower_gen_xargs_clear(g);
    if (g->origin != NULL)
    {
        fmpz_poly_clear(g->origin);
        flint_free(g->origin);
    }
    acb_clear(&g->enclosure);
    flint_free(g->name);
}

/* Renames a generator (updating the variable names of the nested
   contexts which print it). */
void
_gr_tower_gen_rename(gr_tower_t T, gr_tower_gen_struct * g, const char * name)
{
    flint_free(g->name);
    g->name = flint_malloc(strlen(name) + 1);
    strcpy(g->name, name);

    if (g->ctx != NULL)
        GR_MUST_SUCCEED(gr_ctx_set_gen_name(g->ctx, g->name));

    if (g->kind != GR_TOWER_ALGEBRAIC && T->num_trans > 0 && T->base != T->consts)
    {
        char ** names = flint_malloc(sizeof(char *) * T->num_trans);
        slong i;
        for (i = 0; i < T->num_trans; i++)
            names[i] = GR_TOWER_TRANS(T, i)->name;
        GR_MUST_SUCCEED(gr_ctx_set_gen_names(T->base, (const char **) names));
        flint_free(names);
    }
}

/* Records the integer polynomial p as one the generator is a root of. */
void
_gr_tower_gen_set_origin(gr_tower_gen_struct * g, const fmpz_poly_t p)
{
    if (g->origin == NULL)
    {
        g->origin = flint_malloc(sizeof(fmpz_poly_struct));
        fmpz_poly_init(g->origin);
    }
    fmpz_poly_set(g->origin, p);
}

static void
_gr_tower_drop_ctx(gr_tower_t T, gr_ctx_struct * ctx, int retire)
{
    if (retire)
    {
        T->retired = flint_realloc(T->retired, sizeof(gr_ctx_struct *) * (T->num_retired + 1));
        T->retired[T->num_retired++] = ctx;
    }
    else
    {
        gr_ctx_clear(ctx);
        flint_free(ctx);
    }
}

/* Detaches the nested contexts (the step contexts and the base context,
   if owned) from the tower, keeping the generator records. The contexts
   are cleared, or with retire set kept until the tower is cleared, so
   that elements of them held by the caller (including the argument of
   the adjunction that caused the rebuild) can still be used and
   cleared. They are detached from the top down, which is the order in
   which they can be cleared (each step context refers to the ones
   below it). */
static void
_gr_tower_clear_nested(gr_tower_t T, int retire)
{
    slong k;

    for (k = T->length - 1; k >= 0; k--)
    {
        gr_tower_gen_struct * g = GR_TOWER_STEP(T, k);
        if (g->ctx != NULL)
        {
            _gr_tower_drop_ctx(T, g->ctx, retire);
            g->ctx = NULL;
        }
    }

    if (T->base != T->consts)
    {
        _gr_tower_drop_ctx(T, T->base, retire);
        T->base = T->consts;
    }
}

void
gr_tower_clear(gr_tower_t T)
{
    slong d;

    /* nested contexts must be cleared from the top down */
    _gr_tower_clear_nested(T, 0);

    /* each batch of retired contexts was retired from the top down, and
       the batches are independent of each other */
    for (d = 0; d < T->num_retired; d++)
    {
        gr_ctx_clear(T->retired[d]);
        flint_free(T->retired[d]);
    }
    flint_free(T->retired);
    T->retired = NULL;
    T->num_retired = 0;

    for (d = 0; d < T->num_gens; d++)
        _gr_tower_gen_clear(T->gens + d);
    flint_free(T->gens);
    flint_free(T->alg);
    flint_free(T->trans);

    gr_tower_flat_clear(&T->flat);
}

slong
gr_tower_prefix_length(const gr_tower_t T, slong p)
{
    slong d, n = 0;
    if (p < 0 || p > T->num_gens)
        p = T->num_gens;
    for (d = 0; d < p; d++)
        n += (T->gens[d].kind == GR_TOWER_ALGEBRAIC);
    return n;
}

slong
gr_tower_prefix_num_trans(const gr_tower_t T, slong p)
{
    slong d, n = 0;
    if (p < 0 || p > T->num_gens)
        p = T->num_gens;
    for (d = 0; d < p; d++)
        n += (T->gens[d].kind != GR_TOWER_ALGEBRAIC);
    return n;
}

static char *
_copy_name(const char * name, const char * prefix, slong index)
{
    char * s;

    if (name == NULL)
    {
        char buf[32];
        flint_sprintf(buf, "%s%wd", prefix, index);
        s = flint_malloc(strlen(buf) + 1);
        strcpy(s, buf);
    }
    else
    {
        s = flint_malloc(strlen(name) + 1);
        strcpy(s, name);
    }

    return s;
}

/* Appends a blank generator record (def_order = position). */
static gr_tower_gen_struct *
_gr_tower_push_gen(gr_tower_t T, int kind, int status, const char * name, const char * default_prefix, slong default_index)
{
    gr_tower_gen_struct * g;

    if (T->num_gens == T->alloc_gens)
    {
        T->alloc_gens = FLINT_MAX(4, 2 * T->alloc_gens);
        T->gens = flint_realloc(T->gens, T->alloc_gens * sizeof(gr_tower_gen_struct));
        T->alg = flint_realloc(T->alg, T->alloc_gens * sizeof(slong));
        T->trans = flint_realloc(T->trans, T->alloc_gens * sizeof(slong));
    }

    /* names are unique within a tower: a name already in use (e.g.
       copied from another tower) is replaced by a default one */
    if (name != NULL)
    {
        slong i;
        for (i = 0; i < T->num_gens && name != NULL; i++)
            if (strcmp(T->gens[i].name, name) == 0)
                name = NULL;
    }

    g = T->gens + T->num_gens;
    g->kind = kind;
    g->def_kind = kind;
    g->def_param = 0;
    g->status = status;
    g->index = 0;
    g->ctx = NULL;
    g->arg.mctx = NULL;
    g->num_xargs = 0;
    g->xargs = NULL;
    g->origin = NULL;
    acb_init(&g->enclosure);
    g->enclosure_prec = 0;
    /* (a default name is made unique too) */
    if (name == NULL)
    {
        char buf[32];
        slong i;
        for (;; default_index++)
        {
            flint_sprintf(buf, "%s%wd", default_prefix, default_index);
            for (i = 0; i < T->num_gens; i++)
                if (strcmp(T->gens[i].name, buf) == 0)
                    break;
            if (i == T->num_gens)
                break;
        }
    }
    g->name = _copy_name(name, default_prefix, default_index);
    g->def_id = 0;
    g->def_order = T->num_gens;
    g->gid = T->next_gid++;
    g->proof_version = 0;
    g->real = 0;
    T->num_gens++;

    return g;
}

/* Registers the generator in the index arrays (appending: its index
   becomes the last one of its kind). */
static void
_gr_tower_register(gr_tower_t T, gr_tower_gen_struct * g)
{
    if (g->kind == GR_TOWER_ALGEBRAIC)
    {
        T->alg[T->length] = g->def_order;
        g->index = ++T->length;
    }
    else
    {
        T->trans[T->num_trans] = g->def_order;
        g->index = ++T->num_trans;
    }
}

/* Sets the status of a generator, keeping the field flag of its step
   context in sync (a step is a field when its modulus is proven
   irreducible; otherwise it only pretends to be one). */
void
_gr_tower_gen_set_status(gr_tower_gen_struct * g, int status)
{
    g->status = status;
    if (g->kind == GR_TOWER_ALGEBRAIC && g->ctx != NULL && g->ctx->which_ring == GR_CTX_GR_POLY_QUOTIENT)
        GR_MUST_SUCCEED(gr_ctx_set_is_field(g->ctx, (status == GR_TOWER_STATUS_PROVEN) ? T_TRUE : T_UNKNOWN));
}

static void
_gr_tower_set_step_ctx(gr_tower_t T, gr_tower_gen_struct * g, const gr_poly_t m)
{
    gr_ctx_struct * below = gr_tower_field_at(T, g->index - 1);

    g->ctx = flint_malloc(sizeof(gr_ctx_struct));
    gr_ctx_init_gr_poly_quotient(g->ctx, below, m);
    /* the step pretends to be a field (dynamic evaluation); it is one
       when the modulus is proven irreducible */
    GR_MUST_SUCCEED(gr_ctx_set_is_pretend_field(g->ctx, T_TRUE));
    GR_MUST_SUCCEED(gr_ctx_set_is_field(g->ctx, (g->status == GR_TOWER_STATUS_PROVEN) ? T_TRUE : T_UNKNOWN));
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(g->ctx, g->name));
}

/*
    Appends a step F_k = F_{k-1}[x]/(m) with the given enclosure, without
    any verification. Returns the new generator record.
*/
gr_tower_gen_struct *
_gr_tower_push_step(gr_tower_t T, const gr_poly_t m, const acb_t z, slong prec, int status, const char * name)
{
    gr_tower_gen_struct * g;

    g = _gr_tower_push_gen(T, GR_TOWER_ALGEBRAIC, status, name, "a", T->length + 1);
    _gr_tower_register(T, g);
    _gr_tower_set_step_ctx(T, g, m);
    acb_set(&g->enclosure, z);
    g->enclosure_prec = prec;
    T->version++;

    return g;
}

/*
    Rebuilds the base field context and the nested chain from flat moduli:
    moduli[d] (of length lens[d]) for each algebraic generator with
    definition order d, as flat elements in the current context of the
    flat machinery. The nested contexts must have been cleared. The index
    arrays must be up to date.
*/
static void
_gr_tower_build_nested(gr_tower_t T, fmpz_mpoly_q_struct ** moduli, const slong * lens)
{
    gr_tower_flat_struct * F = &T->flat;
    slong k, i;
    slong saved_length = T->length;

    if (T->num_trans > 0)
    {
        char ** names = flint_malloc(sizeof(char *) * T->num_trans);
        T->base = flint_malloc(sizeof(gr_ctx_struct));
        gr_ctx_init_fmpz_mpoly_q(T->base, T->num_trans, ORD_LEX);
        for (i = 0; i < T->num_trans; i++)
            names[i] = GR_TOWER_TRANS(T, i)->name;
        GR_MUST_SUCCEED(gr_ctx_set_gen_names(T->base, (const char **) names));
        flint_free(names);
    }

    T->version++;
    T->moduli_version++;

    /* the ideal is rebuilt step by step: the conversion of the modulus of
       step k needs the reductions of steps < k */
    T->length = 0;
    gr_tower_flat_ensure(F);

    for (k = 1; k <= saved_length; k++)
    {
        gr_tower_gen_struct * g = T->gens + T->alg[k - 1];
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        gr_poly_t m;

        gr_poly_init(m, below);
        gr_poly_fit_length(m, lens[g->def_order], below);
        for (i = 0; i < lens[g->def_order]; i++)
        {
            int st = gr_tower_flat_get_nested_at(gr_poly_coeff_ptr(m, i, below), moduli[g->def_order] + i, k - 1, F);
            if (st != GR_SUCCESS)
            {
                flint_throw(FLINT_ERROR, "(%s): cannot convert a modulus to the nested representation (status %d)\n", __func__, st);
            }
        }
        _gr_poly_set_length(m, lens[g->def_order], below);

        _gr_tower_set_step_ctx(T, g, m);
        T->length = k;
        gr_tower_flat_ensure(F);

        gr_poly_clear(m, below);
    }
}

/* Converts the nested moduli of all algebraic generators to flat form
   (indexed by definition order; NULL for other generators). */
static void
_gr_tower_collect_moduli(gr_tower_t T, fmpz_mpoly_q_struct *** moduli_out, slong ** lens_out)
{
    gr_tower_flat_struct * F = &T->flat;
    fmpz_mpoly_q_struct ** moduli;
    slong * lens;
    slong k, i;

    gr_tower_flat_ensure(F);
    moduli = flint_calloc(FLINT_MAX(T->num_gens, 1), sizeof(fmpz_mpoly_q_struct *));
    lens = flint_calloc(FLINT_MAX(T->num_gens, 1), sizeof(slong));

    for (k = 1; k <= T->length; k++)
    {
        gr_tower_gen_struct * g = GR_TOWER_STEP(T, k - 1);
        const gr_poly_struct * m = gr_poly_quotient_ctx_modulus(g->ctx);
        slong d = g->def_order;

        lens[d] = m->length;
        moduli[d] = flint_malloc(sizeof(fmpz_mpoly_q_struct) * m->length);
        for (i = 0; i < m->length; i++)
        {
            fmpz_mpoly_q_init(moduli[d] + i, F->mctx);
            GR_MUST_SUCCEED(gr_tower_flat_set_nested_at(moduli[d] + i,
                gr_poly_coeff_srcptr(m, i, gr_tower_field_at(T, k - 1)), k - 1, F));
        }
    }

    *moduli_out = moduli;
    *lens_out = lens;
}

static void
_gr_tower_free_moduli(gr_tower_t T, fmpz_mpoly_q_struct ** moduli, slong * lens)
{
    slong d, i;
    for (d = 0; d < T->num_gens; d++)
    {
        if (moduli[d] != NULL)
        {
            for (i = 0; i < lens[d]; i++)
                fmpz_mpoly_q_clear(moduli[d] + i, T->flat.mctx);
            flint_free(moduli[d]);
        }
    }
    flint_free(moduli);
    flint_free(lens);
}

/* Rebuilds the nested chain after a structural change which does not
   change the moduli of the algebraic generators (e.g. a new
   transcendental generator). */
static void
_gr_tower_rebase(gr_tower_t T)
{
    fmpz_mpoly_q_struct ** moduli;
    slong * lens;

    _gr_tower_collect_moduli(T, &moduli, &lens);
    _gr_tower_clear_nested(T, T->keep_retired);
    _gr_tower_build_nested(T, moduli, lens);
    _gr_tower_free_moduli(T, moduli, lens);
}

/* Rebuilds the index arrays from the kinds of the generators. */
static void
_gr_tower_reindex(gr_tower_t T)
{
    slong d;

    T->length = 0;
    T->num_trans = 0;
    for (d = 0; d < T->num_gens; d++)
    {
        gr_tower_gen_struct * g = T->gens + d;
        if (g->kind == GR_TOWER_ALGEBRAIC)
        {
            T->alg[T->length] = d;
            g->index = ++T->length;
        }
        else
        {
            T->trans[T->num_trans] = d;
            g->index = ++T->num_trans;
        }
    }
}

/*
    Turns the generator with definition order d into an algebraic
    generator with the monic modulus given by the flat elements
    m[0], ..., m[len - 1] (in mctx), which may involve only generators
    with lower definition order (for an algebraic generator, this replaces
    its modulus, e.g. by a linear one expressing it in terms of other
    generators). The enclosure is kept. Elements in the nested
    representation held by the caller become invalid; flat elements
    remain valid.
*/
int
_gr_tower_make_algebraic(gr_tower_t T, slong d, const fmpz_mpoly_q_struct * m, slong len, const fmpz_mpoly_ctx_t mctx, int status)
{
    gr_tower_gen_struct * g = T->gens + d;
    fmpz_mpoly_q_struct ** moduli;
    fmpz_mpoly_q_struct * mc;
    slong * lens, i;

    /* the coefficients of the modulus must convert to the nested
       representation once the chain is rebuilt, which is guaranteed when
       their denominators are free of algebraic generators (they are then
       units of the base field); others are rationalized first, on a
       copy in the tower's own flat context, and the relation is not
       applied when that fails */
    gr_tower_flat_ensure(&T->flat);
    mc = flint_malloc(sizeof(fmpz_mpoly_q_struct) * len);
    for (i = 0; i < len; i++)
    {
        fmpz_mpoly_q_init(mc + i, T->flat.mctx);
        gr_tower_flat_convert(mc + i, m + i, mctx, &T->flat);
        if (gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(mc + i), &T->flat) &&
            gr_tower_flat_rationalize(mc + i, &T->flat) != GR_SUCCESS)
        {
            slong j;
            for (j = 0; j <= i; j++)
                fmpz_mpoly_q_clear(mc + j, T->flat.mctx);
            flint_free(mc);
            return 0;
        }
    }
    m = mc;
    mctx = T->flat.mctx;

    if (len == 2 && GR_TOWER_OPTION(T, GR_TOWER_OPT_VERBOSE))
    {
        /* (the relation found: g = -m[0] / m[1]) */
        fmpz_mpoly_q_t e;
        char * str;
        fmpz_mpoly_q_init(e, mctx);
        fmpz_mpoly_q_div(e, m + 0, m + 1, mctx);
        fmpz_mpoly_q_neg(e, e, mctx);
        str = _gr_tower_flat_get_str(e, mctx, T);
        flint_printf("gr_tower: %s = %s\n", g->name, str);
        flint_free(str);
        fmpz_mpoly_q_clear(e, mctx);
    }

    /* the enclosure is needed once the generator is algebraic (it is
       what selects the root of the modulus) */
    if (g->kind != GR_TOWER_ALGEBRAIC && g->enclosure_prec == 0)
    {
        acb_t z;
        int st;
        acb_init(z);
        st = gr_tower_trans_get_acb(z, T, g->index, GR_TOWER_DEFAULT_PREC);
        acb_clear(z);
        if (st != GR_SUCCESS)
        {
            for (i = 0; i < len; i++)
                fmpz_mpoly_q_clear(mc + i, T->flat.mctx);
            flint_free(mc);
            return 0;
        }
    }

    _gr_tower_collect_moduli(T, &moduli, &lens);

    if (moduli[d] != NULL)
    {
        for (i = 0; i < lens[d]; i++)
            fmpz_mpoly_q_clear(moduli[d] + i, T->flat.mctx);
        flint_free(moduli[d]);
    }

    /* a generator which was transcendental may have carried the proof of
       irreducibility of later moduli over the base field */
    if (g->kind != GR_TOWER_ALGEBRAIC)
    {
        for (i = d + 1; i < T->num_gens; i++)
            if (T->gens[i].kind == GR_TOWER_ALGEBRAIC && T->gens[i].status == GR_TOWER_STATUS_PROVEN && lens[i] > 2)
                _gr_tower_gen_set_status(T->gens + i, GR_TOWER_STATUS_DYNAMIC);
    }

    lens[d] = len;
    moduli[d] = flint_malloc(sizeof(fmpz_mpoly_q_struct) * len);
    for (i = 0; i < len; i++)
    {
        fmpz_mpoly_q_init(moduli[d] + i, T->flat.mctx);
        gr_tower_flat_convert(moduli[d] + i, m + i, mctx, &T->flat);
    }

    _gr_tower_clear_nested(T, T->keep_retired);

    g->kind = GR_TOWER_ALGEBRAIC;
    g->status = status;
    _gr_tower_reindex(T);

    _gr_tower_build_nested(T, moduli, lens);
    _gr_tower_free_moduli(T, moduli, lens);

    /* (the copies were made in the flat context before the rebuild,
       which may have replaced it, keeping the old one alive: they are
       cleared in that context) */
    for (i = 0; i < len; i++)
        fmpz_mpoly_q_clear(mc + i, (fmpz_mpoly_ctx_struct *) mctx);
    flint_free(mc);
    return 1;
}

/*
    Moves the generator with definition order d to position p <= d of
    the definition order. The generator must not depend on any generator
    with definition order in [p, d). The flat contexts are replaced
    (elements in the old ones remain convertible) and the nested chain
    is rebuilt.
*/
void
_gr_tower_move_gen(gr_tower_t T, slong d, slong p)
{
    fmpz_mpoly_q_struct ** moduli;
    slong * lens;
    gr_tower_gen_struct g;
    slong i;
    fmpz_mpoly_ctx_struct * old_mctx;

    if (d == p)
        return;

    _gr_tower_collect_moduli(T, &moduli, &lens);
    old_mctx = T->flat.mctx;

    _gr_tower_clear_nested(T, T->keep_retired);

    /* reorder the records */
    g = T->gens[d];
    for (i = d; i > p; i--)
        T->gens[i] = T->gens[i - 1];
    T->gens[p] = g;
    for (i = 0; i < T->num_gens; i++)
        T->gens[i].def_order = i;
    {
        fmpz_mpoly_q_struct * mg = moduli[d];
        slong lg = lens[d];
        for (i = d; i > p; i--)
        {
            moduli[i] = moduli[i - 1];
            lens[i] = lens[i - 1];
        }
        moduli[p] = mg;
        lens[p] = lg;
    }

    T->structure_version++;
    _gr_tower_reindex(T);

    /* the generators which now come after the moved one are defined
       over a larger field: proofs of irreducibility of their moduli no
       longer apply, unless a structural rule covers the new prefix */
    if (g.kind == GR_TOWER_ALGEBRAIC)   /* a transcendental extension preserves irreducibility */
    {
        for (i = p + 1; i <= d; i++)
        {
            gr_tower_gen_struct * h = T->gens + i;
            if (h->kind == GR_TOWER_ALGEBRAIC && h->status == GR_TOWER_STATUS_PROVEN && lens[i] > 2 &&
                !_gr_tower_structural_proven(h, T, i))
                _gr_tower_gen_set_status(h, GR_TOWER_STATUS_DYNAMIC);
        }
    }

    /* new flat layout; the moduli are converted from the old context */
    _gr_tower_flat_grow(&T->flat);
    for (i = 0; i < T->num_gens; i++)
    {
        slong j;
        for (j = 0; j < lens[i]; j++)
        {
            fmpz_mpoly_q_t t;
            fmpz_mpoly_q_init(t, T->flat.mctx);
            gr_tower_flat_convert(t, moduli[i] + j, old_mctx, &T->flat);
            fmpz_mpoly_q_clear(moduli[i] + j, old_mctx);
            fmpz_mpoly_q_init(moduli[i] + j, T->flat.mctx);
            fmpz_mpoly_q_swap(moduli[i] + j, t, T->flat.mctx);
            fmpz_mpoly_q_clear(t, T->flat.mctx);
        }
    }

    _gr_tower_build_nested(T, moduli, lens);
    _gr_tower_free_moduli(T, moduli, lens);
}

/*
    Sets the generators with definition orders d[0], ..., d[n - 1] to the
    values v[0], ..., v[n - 1] (in mctx): they become algebraic with the
    linear moduli X - v[i] (the generators of the values, and those of
    their moduli and arguments, may come later in the definition order),
    and the generators are reordered by the stable topological order of
    the dependencies (moduli, arguments), with one rebuild of the nested
    chain. Returns 0, leaving the tower unchanged, if the dependencies
    are cyclic (a value depending on its own generator), a value cannot
    be rationalized or a generator cannot be evaluated. Flat elements
    remain valid (the flat contexts are replaced).
*/
int
_gr_tower_set_linear_gens(gr_tower_t T, slong n, const slong * dd, const fmpz_mpoly_q_struct * v, const fmpz_mpoly_ctx_t mctx)
{
    gr_tower_flat_struct * F = &T->flat;
    fmpz_mpoly_q_struct ** moduli, * mc;
    fmpz_mpoly_ctx_struct * old_mctx;
    slong * lens, * order, * newpos, N = T->num_gens, i, j, e, f, nvars;
    char * dep, * placed, * target, * was_alg;
    int * used;
    gr_tower_gen_struct * gens;
    int ok = 1;

    if (n == 0)
        return 1;

    /* the enclosures of the transcendental targets (they select the
       root of the modulus once algebraic) */
    for (i = 0; i < n; i++)
    {
        gr_tower_gen_struct * g = T->gens + dd[i];
        if (g->kind != GR_TOWER_ALGEBRAIC && g->enclosure_prec == 0)
        {
            acb_t z;
            acb_init(z);
            ok = (gr_tower_trans_get_acb(z, T, g->index, GR_TOWER_DEFAULT_PREC) == GR_SUCCESS);
            acb_clear(z);
            if (!ok)
                return 0;
        }
    }

    _gr_tower_collect_moduli(T, &moduli, &lens);
    old_mctx = F->mctx;
    nvars = old_mctx->minfo->nvars;

    /* the values, rationalized, in the current context */
    mc = flint_malloc(sizeof(fmpz_mpoly_q_struct) * n);
    for (i = 0; i < n; i++)
        fmpz_mpoly_q_init(mc + i, old_mctx);
    for (i = 0; i < n && ok; i++)
    {
        gr_tower_flat_convert(mc + i, v + i, mctx, F);
        if (gr_tower_flat_has_alg_var(fmpz_mpoly_q_denref(mc + i), F) &&
            gr_tower_flat_rationalize(mc + i, F) != GR_SUCCESS)
            ok = 0;
    }

    /* the dependencies with the new moduli: dep[e N + f], e depends on f */
    dep = flint_calloc(N * N + 1, 1);
    target = flint_calloc(N + 1, 1);
    used = flint_malloc(sizeof(int) * FLINT_MAX(nvars, 1));
    for (i = 0; i < n; i++)
        target[dd[i]] = 1;
    for (e = 0; e < N && ok; e++)
    {
        const gr_tower_gen_struct * g = T->gens + e;
        slong na = _gr_tower_gen_num_args(g), m;
        slong nterms = target[e] ? 1 : lens[e] + na;

        for (m = 0; m < nterms; m++)
        {
            fmpz_mpoly_q_t a;
            slong w, q;
            fmpz_mpoly_q_init(a, old_mctx);
            if (target[e])
            {
                for (q = 0; q < n; q++)
                    if (dd[q] == e)
                        fmpz_mpoly_q_set(a, mc + q, old_mctx);
            }
            else if (m < lens[e])
                fmpz_mpoly_q_set(a, moduli[e] + m, old_mctx);
            else
            {
                const gr_tower_flat_elem_struct * ar = _gr_tower_gen_arg_ptr(g, m - lens[e]);
                if (ar->mctx != NULL)
                    gr_tower_flat_convert(a, &ar->data, ar->mctx, F);
            }
            for (q = 0; q < 2; q++)
            {
                fmpz_mpoly_used_vars(used, q ? fmpz_mpoly_q_denref(a) : fmpz_mpoly_q_numref(a), old_mctx);
                for (w = 0; w < nvars; w++)
                {
                    f = F->cap - 1 - w;
                    if (used[w] && f >= 0 && f < N && f != e)
                        dep[e * N + f] = 1;
                    else if (used[w] && f == e)
                        ok = 0;    /* (a value involving its own generator) */
                }
            }
            fmpz_mpoly_q_clear(a, old_mctx);
        }
    }
    flint_free(used);

    /* the stable topological order */
    order = flint_malloc(sizeof(slong) * (N + 1));
    newpos = flint_malloc(sizeof(slong) * (N + 1));
    placed = flint_calloc(N + 1, 1);
    for (i = 0; i < N && ok; i++)
    {
        for (e = 0; e < N; e++)
        {
            if (placed[e])
                continue;
            for (f = 0; f < N; f++)
                if (dep[e * N + f] && !placed[f])
                    break;
            if (f == N)
                break;
        }
        if (e == N)
            ok = 0;
        else
        {
            order[i] = e;
            newpos[e] = i;
            placed[e] = 1;
        }
    }

    if (!ok)
    {
        for (i = 0; i < n; i++)
            fmpz_mpoly_q_clear(mc + i, old_mctx);
        flint_free(mc);
        flint_free(dep);
        flint_free(target);
        flint_free(order);
        flint_free(newpos);
        flint_free(placed);
        _gr_tower_free_moduli(T, moduli, lens);
        return 0;
    }

    /* the new moduli */
    for (i = 0; i < n; i++)
    {
        e = dd[i];
        if (moduli[e] != NULL)
        {
            for (j = 0; j < lens[e]; j++)
                fmpz_mpoly_q_clear(moduli[e] + j, old_mctx);
            flint_free(moduli[e]);
        }
        lens[e] = 2;
        moduli[e] = flint_malloc(sizeof(fmpz_mpoly_q_struct) * 2);
        fmpz_mpoly_q_init(moduli[e] + 0, old_mctx);
        fmpz_mpoly_q_init(moduli[e] + 1, old_mctx);
        fmpz_mpoly_q_neg(moduli[e] + 0, mc + i, old_mctx);
        fmpz_mpoly_q_one(moduli[e] + 1, old_mctx);
    }
    for (i = 0; i < n; i++)
        fmpz_mpoly_q_clear(mc + i, old_mctx);
    flint_free(mc);

    _gr_tower_clear_nested(T, T->keep_retired);

    /* the irreducibility of a modulus (of degree > 2) proven over the
       field of the generators before it holds over a field with
       transcendental generators added only: otherwise checked again,
       unless a structural rule covers the new prefix */
    was_alg = flint_malloc(N + 1);
    for (e = 0; e < N; e++)
        was_alg[e] = (T->gens[e].kind == GR_TOWER_ALGEBRAIC);
    for (e = 0; e < N; e++)
    {
        int changed = 0;
        for (f = 0; f < N && !changed; f++)
        {
            int bo = (f < e), bn = (newpos[f] < newpos[e]);
            changed = (bo && !bn) || (bn && target[f]) || (bn && !bo && was_alg[f]);
        }
        placed[e] = changed && was_alg[e] && !target[e] &&
            T->gens[e].status == GR_TOWER_STATUS_PROVEN && lens[e] > 2;
    }

    /* reorder the records and the moduli */
    gens = flint_malloc(sizeof(gr_tower_gen_struct) * (N + 1));
    {
        fmpz_mpoly_q_struct ** m2 = flint_calloc(N + 1, sizeof(fmpz_mpoly_q_struct *));
        slong * l2 = flint_calloc(N + 1, sizeof(slong));
        char * ch = flint_malloc(N + 1);
        for (i = 0; i < N; i++)
        {
            gens[i] = T->gens[order[i]];
            m2[i] = moduli[order[i]];
            l2[i] = lens[order[i]];
            ch[i] = placed[order[i]];   /* (the flags of the checks) */
            if (target[order[i]])
            {
                gens[i].kind = GR_TOWER_ALGEBRAIC;
                gens[i].status = GR_TOWER_STATUS_PROVEN;
            }
        }
        for (i = 0; i < N; i++)
        {
            T->gens[i] = gens[i];
            T->gens[i].def_order = i;
            moduli[i] = m2[i];
            lens[i] = l2[i];
        }
        T->structure_version++;
        _gr_tower_reindex(T);
        for (i = 0; i < N; i++)
            if (ch[i] && !_gr_tower_structural_proven(T->gens + i, T, i))
                _gr_tower_gen_set_status(T->gens + i, GR_TOWER_STATUS_DYNAMIC);
        flint_free(m2);
        flint_free(l2);
        flint_free(ch);
    }
    flint_free(gens);
    flint_free(was_alg);
    flint_free(dep);
    flint_free(target);
    flint_free(order);
    flint_free(newpos);
    flint_free(placed);

    /* new flat layout; the moduli are converted from the old context */
    _gr_tower_flat_grow(F);
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < lens[i]; j++)
        {
            fmpz_mpoly_q_t t;
            fmpz_mpoly_q_init(t, F->mctx);
            gr_tower_flat_convert(t, moduli[i] + j, old_mctx, F);
            fmpz_mpoly_q_clear(moduli[i] + j, old_mctx);
            fmpz_mpoly_q_init(moduli[i] + j, F->mctx);
            fmpz_mpoly_q_swap(moduli[i] + j, t, F->mctx);
            fmpz_mpoly_q_clear(t, F->mctx);
        }
    }

    _gr_tower_build_nested(T, moduli, lens);
    _gr_tower_free_moduli(T, moduli, lens);
    return 1;
}

static gr_tower_gen_struct *
_gr_tower_push_trans_multi(gr_tower_t T, int kind, const fmpz_mpoly_q_struct * u, slong nargs, const fmpz_mpoly_ctx_struct * u_mctx, int status, const char * name)
{
    gr_tower_gen_struct * g;
    gr_tower_flat_struct * F = &T->flat;
    slong i;

    gr_tower_flat_ensure(F);

    g = _gr_tower_push_gen(T, kind, status, name, "t", T->num_trans + 1);
    _gr_tower_register(T, g);

    /* the arguments are stored after growing the flat context if needed */
    _gr_tower_flat_grow(F);
    if (u != NULL && nargs >= 1)
    {
        g->arg.mctx = F->mctx;
        fmpz_mpoly_q_init(&g->arg.data, g->arg.mctx);
        gr_tower_flat_convert(&g->arg.data, u, u_mctx, F);
    }
    if (u != NULL && nargs > 1)
    {
        g->num_xargs = nargs - 1;
        g->xargs = flint_malloc(sizeof(gr_tower_flat_elem_struct) * (nargs - 1));
        for (i = 1; i < nargs; i++)
        {
            g->xargs[i - 1].mctx = F->mctx;
            fmpz_mpoly_q_init(&g->xargs[i - 1].data, F->mctx);
            gr_tower_flat_convert(&g->xargs[i - 1].data, u + i, u_mctx, F);
        }
    }

    _gr_tower_rebase(T);

    return g;
}

static gr_tower_gen_struct *
_gr_tower_push_trans(gr_tower_t T, int kind, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_struct * u_mctx, int status, const char * name)
{
    return _gr_tower_push_trans_multi(T, kind, u, 1, u_mctx, status, name);
}

int
gr_tower_adjoin_pi(gr_tower_t T, const char * name)
{
    slong d;
    /* (one pi per tower: two proven transcendental generators of the same
       value would make the zero test of their difference wrong) */
    for (d = 0; d < T->num_gens; d++)
        if (T->gens[d].def_kind == GR_TOWER_PI)
            return GR_DOMAIN;
    /* (pi is algebraically independent of the algebraic generators
       before it; a second transcendental generator could not get this
       status by transcendence alone) */
    _gr_tower_push_trans(T, GR_TOWER_PI, NULL, NULL, (T->num_trans == 0) ? GR_TOWER_STATUS_INDEPENDENT : GR_TOWER_STATUS_SCHANUEL, name == NULL ? "pi" : name);
    return GR_SUCCESS;
}

int
gr_tower_adjoin_free(gr_tower_t T, const char * name)
{
    /* a formal variable: independent of everything by definition, with
       no numerical value (the enclosure stays indeterminate) and no
       relation search */
    gr_tower_gen_struct * g = _gr_tower_push_trans(T, GR_TOWER_FREE, NULL, NULL, GR_TOWER_STATUS_INDEPENDENT, name);
    acb_indeterminate(&g->enclosure);
    return GR_SUCCESS;
}

/* Zero test of u - c for a flat element u (c = 0 or 1), using the
   complete zero test of the nested representation. */
static truth_t
_flat_is_const(const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, slong c, gr_tower_t T)
{
    gr_ctx_struct * top;
    gr_ptr v;
    fmpz_mpoly_q_t w;
    truth_t res;

    gr_tower_flat_ensure(&T->flat);
    fmpz_mpoly_q_init(w, T->flat.mctx);
    gr_tower_flat_convert(w, u, mctx, &T->flat);
    if (c != 0)
        fmpz_mpoly_q_sub_si(w, w, c, T->flat.mctx);

    top = gr_tower_field(T);
    GR_TMP_INIT(v, top);
    if (gr_tower_flat_get_nested_at(v, w, T->length, &T->flat) != GR_SUCCESS)
        res = T_UNKNOWN;
    else
        res = gr_tower_is_zero(v, T);
    GR_TMP_CLEAR(v, top);
    fmpz_mpoly_q_clear(w, T->flat.mctx);

    return res;
}

/*
    The trivial cases exp(0) and log(1) are rejected: these would be
    generators equal to a constant, contradicting their status.
*/
int
gr_tower_adjoin_exp_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
{
    truth_t zero = _flat_is_const(u, mctx, 0, T);

    if (zero == T_UNKNOWN)
        return GR_UNABLE;
    if (zero == T_TRUE)
        return GR_DOMAIN;

    {
        /* pushing may replace the flat context; old contexts stay alive */
        fmpz_mpoly_q_t w;
        fmpz_mpoly_ctx_struct * wctx;

        gr_tower_flat_ensure(&T->flat);
        wctx = T->flat.mctx;
        fmpz_mpoly_q_init(w, wctx);
        gr_tower_flat_convert(w, u, mctx, &T->flat);
        _gr_tower_push_trans(T, GR_TOWER_EXP, w, wctx, GR_TOWER_STATUS_SCHANUEL, name);
        fmpz_mpoly_q_clear(w, wctx);
    }
    return GR_SUCCESS;
}

int
gr_tower_adjoin_log_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
{
    truth_t zero, one;

    zero = _flat_is_const(u, mctx, 0, T);
    if (zero == T_UNKNOWN)
        return GR_UNABLE;
    if (zero == T_TRUE)
        return GR_DOMAIN;

    one = _flat_is_const(u, mctx, 1, T);
    if (one == T_UNKNOWN)
        return GR_UNABLE;
    if (one == T_TRUE)
        return GR_DOMAIN;

    {
        fmpz_mpoly_q_t w;
        fmpz_mpoly_ctx_struct * wctx;

        gr_tower_flat_ensure(&T->flat);
        wctx = T->flat.mctx;
        fmpz_mpoly_q_init(w, wctx);
        gr_tower_flat_convert(w, u, mctx, &T->flat);
        _gr_tower_push_trans(T, GR_TOWER_LOG, w, wctx, GR_TOWER_STATUS_SCHANUEL, name);
        fmpz_mpoly_q_clear(w, wctx);
    }
    return GR_SUCCESS;
}

/* tan(u) or atan(u) for a real nonzero u */
static int
_adjoin_trig_flat(gr_tower_t T, int kind, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
{
    truth_t zero = _flat_is_const(u, mctx, 0, T);
    fmpz_mpoly_q_t w;
    fmpz_mpoly_ctx_struct * wctx;

    if (zero == T_UNKNOWN)
        return GR_UNABLE;
    if (zero == T_TRUE)
        return GR_DOMAIN;

    gr_tower_flat_ensure(&T->flat);
    wctx = T->flat.mctx;
    fmpz_mpoly_q_init(w, wctx);
    gr_tower_flat_convert(w, u, mctx, &T->flat);

    /* a numerically non-real argument is rejected */
    {
        acb_t z;
        int nonreal;
        acb_init(z);
        nonreal = (gr_tower_flat_get_acb(z, w, GR_TOWER_DEFAULT_PREC, &T->flat) == GR_SUCCESS) &&
                  !arb_contains_zero(acb_imagref(z));
        acb_clear(z);
        if (nonreal)
        {
            fmpz_mpoly_q_clear(w, wctx);
            return GR_DOMAIN;
        }
    }

    /* (the realness of u is the caller's promise, recorded on the
       generator: its enclosures are made exactly real) */
    _gr_tower_push_trans(T, kind, w, wctx, GR_TOWER_STATUS_SCHANUEL, name)->real = 1;
    fmpz_mpoly_q_clear(w, wctx);
    return GR_SUCCESS;
}

int
gr_tower_adjoin_tan_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
{
    return _adjoin_trig_flat(T, GR_TOWER_TAN, u, mctx, name);
}

int
gr_tower_adjoin_atan_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
{
    return _adjoin_trig_flat(T, GR_TOWER_ATAN, u, mctx, name);
}

/* the nested-argument versions of the adjunctions: the argument u (an
   element of the top field; NULL for a constant) is converted to the
   flat representation */
static int
_adjoin_with_arg(gr_tower_t T, int kind, slong param, gr_srcptr u, const char * name)
{
    fmpz_mpoly_q_t w;
    fmpz_mpoly_ctx_struct * wctx;
    int status = GR_SUCCESS;

    gr_tower_flat_ensure(&T->flat);
    wctx = T->flat.mctx;
    fmpz_mpoly_q_init(w, wctx);
    if (u != NULL)
        status = gr_tower_flat_set_nested_at(w, u, T->length, &T->flat);
    if (status == GR_SUCCESS)
    {
        if (kind == GR_TOWER_EXP)
            status = gr_tower_adjoin_exp_flat(T, w, wctx, name);
        else if (kind == GR_TOWER_LOG)
            status = gr_tower_adjoin_log_flat(T, w, wctx, name);
        else if (kind == GR_TOWER_TAN || kind == GR_TOWER_ATAN)
            status = _adjoin_trig_flat(T, kind, w, wctx, name);
        else
            status = gr_tower_adjoin_special_flat(T, kind, param, w, wctx, name);
    }
    fmpz_mpoly_q_clear(w, wctx);
    return status;
}

int gr_tower_adjoin_tan(gr_tower_t T, gr_srcptr u, const char * name) { return _adjoin_with_arg(T, GR_TOWER_TAN, 0, u, name); }
int gr_tower_adjoin_atan(gr_tower_t T, gr_srcptr u, const char * name) { return _adjoin_with_arg(T, GR_TOWER_ATAN, 0, u, name); }

/* special function values (check = 0: the value is known to be valid,
   as when a generator of another tower is absorbed) */
static int
_adjoin_special_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name, int check);

int
gr_tower_adjoin_special_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
{
    return _adjoin_special_flat(T, kind, param, u, mctx, name, 1);
}

int
_gr_tower_adjoin_special_flat_nocheck(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
{
    return _adjoin_special_flat(T, kind, param, u, mctx, name, 0);
}

static int
_adjoin_special_multi_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_struct * u, slong nargs,
    const fmpz_mpoly_ctx_t mctx, const char * name, int check)
{
    fmpz_mpoly_q_struct * w;
    fmpz_mpoly_ctx_struct * wctx;
    gr_tower_gen_struct * g;
    int status = GR_SUCCESS;
    slong i;

    if (!GR_TOWER_KIND_IS_SPECIAL(kind))
        return GR_DOMAIN;
    if (kind == GR_TOWER_CONSTANT)
        nargs = 0;
    else if (nargs != _gr_tower_special_num_args(kind, param))
        return GR_DOMAIN;

    gr_tower_flat_ensure(&T->flat);
    wctx = T->flat.mctx;
    w = flint_malloc(sizeof(fmpz_mpoly_q_struct) * FLINT_MAX(nargs, 1));
    for (i = 0; i < nargs; i++)
    {
        fmpz_mpoly_q_init(w + i, wctx);
        gr_tower_flat_convert(w + i, u + i, mctx, &T->flat);
    }

    if (nargs > 0)
    {
        /* u = 0: a trivial value (erf, W_0, Li_s, pFq) or a pole */
        if (check && (kind == GR_TOWER_ERF || kind == GR_TOWER_LAMBERTW || kind == GR_TOWER_POLYLOG ||
            kind == GR_TOWER_GAMMA || kind == GR_TOWER_POLYGAMMA || kind == GR_TOWER_HYPGEOM))
        {
            truth_t zero = _flat_is_const(w, wctx, 0, T);
            if (zero == T_UNKNOWN)
                status = GR_UNABLE;
            else if (zero == T_TRUE)
                status = GR_DOMAIN;
        }

        /* the value must be computable */
        if (status == GR_SUCCESS)
        {
            acb_ptr z;
            acb_t v;
            z = _acb_vec_init(nargs);
            acb_init(v);
            for (i = 0; i < nargs && status == GR_SUCCESS; i++)
                if (gr_tower_flat_get_acb(z + i, w + i, GR_TOWER_DEFAULT_PREC, &T->flat) != GR_SUCCESS)
                    status = GR_UNABLE;
            if (status == GR_SUCCESS &&
                _gr_tower_special_eval_multi_flags(v, kind, param, z, nargs,
                    (kind == GR_TOWER_HYPGEOM && GR_TOWER_HYPGEOM_P(param) == 2 && GR_TOWER_HYPGEOM_Q(param) == 1 && nargs == 4) ?
                        _gr_tower_hypgeom_2f1_flags(w + 1, w + 2, w + 3, wctx) : 0,
                    GR_TOWER_DEFAULT_PREC) != GR_SUCCESS)
                status = GR_UNABLE;
            _acb_vec_clear(z, nargs);
            acb_clear(v);
        }
    }

    if (status == GR_SUCCESS)
    {
        g = _gr_tower_push_trans_multi(T, kind, (nargs == 0) ? NULL : w, nargs, wctx,
            (kind == GR_TOWER_LAMBERTW) ? GR_TOWER_STATUS_SCHANUEL : GR_TOWER_STATUS_CONJECTURAL, name);
        g->def_param = param;
    }

    for (i = 0; i < nargs; i++)
        fmpz_mpoly_q_clear(w + i, wctx);
    flint_free(w);
    return status;
}

static int
_adjoin_special_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name, int check)
{
    return _adjoin_special_multi_flat(T, kind, param, u, (kind == GR_TOWER_CONSTANT) ? 0 : 1, mctx, name, check);
}

int
gr_tower_adjoin_special_multi_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_struct * u, slong nargs, const fmpz_mpoly_ctx_t mctx, const char * name)
{
    return _adjoin_special_multi_flat(T, kind, param, u, nargs, mctx, name, 1);
}

int
_gr_tower_adjoin_special_multi_flat_nocheck(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_struct * u, slong nargs, const fmpz_mpoly_ctx_t mctx, const char * name)
{
    return _adjoin_special_multi_flat(T, kind, param, u, nargs, mctx, name, 0);
}

int
gr_tower_adjoin_special(gr_tower_t T, int kind, slong param, gr_srcptr u, const char * name)
{
    return _adjoin_with_arg(T, kind, param, (kind == GR_TOWER_CONSTANT) ? NULL : u, name);
}

int gr_tower_adjoin_exp(gr_tower_t T, gr_srcptr u, const char * name) { return _adjoin_with_arg(T, GR_TOWER_EXP, 0, u, name); }
int gr_tower_adjoin_log(gr_tower_t T, gr_srcptr u, const char * name) { return _adjoin_with_arg(T, GR_TOWER_LOG, 0, u, name); }

/*
    Copies the first p generators of T. All data is transported through
    the flat representation; the nested chain of res is built at the end.
*/
/*
    Sets res to the tower of the generators of T marked in mark[] (by
    definition order), which must be closed under the dependencies of the
    definitions (the generators occurring in the moduli and arguments of
    marked generators are marked). The statuses carry over: a modulus
    irreducible over the field of all the preceding generators is
    irreducible over the smaller field of the marked ones.
*/
void
gr_tower_set_subset(gr_tower_t res, const gr_tower_t T_, const int * mark)
{
    gr_tower_struct * T = (gr_tower_struct *) T_;
    fmpz_mpoly_q_struct ** moduli, ** rmoduli;
    slong * lens, * rlens, * order_map;
    slong d, i, p, nd;

    if (res == T)
        flint_throw(FLINT_ERROR, "(%s): aliasing not supported\n", __func__);

    gr_tower_clear(res);
    gr_tower_init(res, T->consts);
    res->options = T->options;
    res->keep_retired = T->keep_retired;

    _gr_tower_collect_moduli(T, &moduli, &lens);

    order_map = flint_malloc(sizeof(slong) * FLINT_MAX(T->num_gens, 1));
    for (d = 0, p = 0; d < T->num_gens; d++)
        order_map[d] = mark[d] ? p++ : -1;

    /* generator records */
    for (d = 0; d < T->num_gens; d++)
    {
        const gr_tower_gen_struct * g = T->gens + d;
        gr_tower_gen_struct * ng;

        if (!mark[d])
            continue;

        ng = _gr_tower_push_gen(res, g->kind, g->status, g->name, "x", order_map[d]);
        ng->def_kind = g->def_kind;
        ng->def_param = g->def_param;
        _gr_tower_register(res, ng);
        acb_set(&ng->enclosure, &g->enclosure);
        ng->enclosure_prec = g->enclosure_prec;
        ng->def_id = g->def_id;
        ng->real = g->real;
        if (g->origin != NULL)
            _gr_tower_gen_set_origin(ng, g->origin);
    }

    /* flat data, once the capacity of res is known */
    _gr_tower_flat_grow(&res->flat);

    rmoduli = flint_calloc(FLINT_MAX(p, 1), sizeof(fmpz_mpoly_q_struct *));
    rlens = flint_calloc(FLINT_MAX(p, 1), sizeof(slong));

    for (d = 0; d < T->num_gens; d++)
    {
        const gr_tower_gen_struct * g = T->gens + d;
        gr_tower_gen_struct * ng;

        if (!mark[d])
            continue;

        nd = order_map[d];
        ng = res->gens + nd;

        if (g->arg.mctx != NULL)
        {
            fmpz_mpoly_q_t u;
            fmpz_mpoly_q_init(u, T->flat.mctx);
            gr_tower_flat_convert(u, &g->arg.data, g->arg.mctx, &T->flat);
            ng->arg.mctx = res->flat.mctx;
            fmpz_mpoly_q_init(&ng->arg.data, ng->arg.mctx);
            _gr_tower_flat_transport_map(&ng->arg.data, u, &T->flat, &res->flat, order_map);
            fmpz_mpoly_q_clear(u, T->flat.mctx);
        }

        if (g->num_xargs > 0)
        {
            fmpz_mpoly_q_t u;
            ng->num_xargs = g->num_xargs;
            ng->xargs = flint_malloc(sizeof(gr_tower_flat_elem_struct) * g->num_xargs);
            fmpz_mpoly_q_init(u, T->flat.mctx);
            for (i = 0; i < g->num_xargs; i++)
            {
                gr_tower_flat_convert(u, &g->xargs[i].data, g->xargs[i].mctx, &T->flat);
                ng->xargs[i].mctx = res->flat.mctx;
                fmpz_mpoly_q_init(&ng->xargs[i].data, ng->xargs[i].mctx);
                _gr_tower_flat_transport_map(&ng->xargs[i].data, u, &T->flat, &res->flat, order_map);
            }
            fmpz_mpoly_q_clear(u, T->flat.mctx);
        }

        if (moduli[d] != NULL)
        {
            rlens[nd] = lens[d];
            rmoduli[nd] = flint_malloc(sizeof(fmpz_mpoly_q_struct) * lens[d]);
            for (i = 0; i < lens[d]; i++)
            {
                fmpz_mpoly_q_init(rmoduli[nd] + i, res->flat.mctx);
                _gr_tower_flat_transport_map(rmoduli[nd] + i, moduli[d] + i, &T->flat, &res->flat, order_map);
            }
        }
    }

    _gr_tower_build_nested(res, rmoduli, rlens);

    _gr_tower_free_moduli(res, rmoduli, rlens);
    _gr_tower_free_moduli(T, moduli, lens);
    flint_free(order_map);
}

void
gr_tower_set_prefix(gr_tower_t res, const gr_tower_t T, slong p)
{
    int * mark;
    slong d;

    if (res == T)
        flint_throw(FLINT_ERROR, "(%s): aliasing not supported\n", __func__);

    if (p < 0 || p > T->num_gens)
        p = T->num_gens;

    res->options = T->options;
    mark = flint_malloc(sizeof(int) * FLINT_MAX(T->num_gens, 1));
    for (d = 0; d < T->num_gens; d++)
        mark[d] = (d < p);
    gr_tower_set_subset(res, T, mark);
    flint_free(mark);
}

void
gr_tower_set(gr_tower_t res, const gr_tower_t T)
{
    if (res == T)
        return;
    gr_tower_set_prefix(res, T, T->num_gens);
}

void
gr_tower_gen_set_def_id(gr_tower_t T, slong d, ulong def_id)
{
    GR_TOWER_GEN(T, d)->def_id = def_id;
}

slong
gr_tower_find_def_order(const gr_tower_t T, ulong def_id)
{
    slong d;

    if (def_id == 0)
        return -1;

    for (d = 0; d < T->num_gens; d++)
        if (T->gens[d].def_id == def_id)
            return d;

    return -1;
}

slong
gr_tower_find_def(const gr_tower_t T, ulong def_id)
{
    slong d = gr_tower_find_def_order(T, def_id);
    if (d >= 0 && T->gens[d].kind == GR_TOWER_ALGEBRAIC)
        return T->gens[d].index;
    return 0;
}

slong
gr_tower_find_trans_def(const gr_tower_t T, ulong def_id)
{
    slong d = gr_tower_find_def_order(T, def_id);
    if (d >= 0 && T->gens[d].kind != GR_TOWER_ALGEBRAIC)
        return T->gens[d].index;
    return 0;
}

gr_tower_flat_struct *
gr_tower_flat(gr_tower_t T)
{
    gr_tower_flat_ensure(&T->flat);
    return &T->flat;
}

slong
gr_tower_step_degree(const gr_tower_t T, slong k)
{
    return gr_poly_quotient_ctx_degree(GR_TOWER_STEP(T, k - 1)->ctx);
}

/* Trial division of n by the first num_primes primes, completed by
   n_factor when the remaining cofactor fits in a word. Returns 1 when
   the factorization is complete. Otherwise the unfactored part is stored
   in cofactor when it is not NULL, and appended to fac (with exponent 1)
   when it is. */
int
_gr_tower_fmpz_factor_trial(fmpz_factor_t fac, fmpz_t cofactor, const fmpz_t n, slong num_primes)
{
    fmpz_t c;
    slong i;
    int complete;

    num_primes = FLINT_MAX(num_primes, 1);

    /* (fmpz_factor_trial_range does not report the cofactor) */
    complete = fmpz_factor_trial_range(fac, n, 0, num_primes);

    fmpz_init(c);
    if (!complete)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_abs(c, n);
        for (i = 0; i < fac->num; i++)
        {
            fmpz_pow_ui(t, fac->p + i, fac->exp[i]);
            fmpz_divexact(c, c, t);
        }
        fmpz_clear(t);

        if (fmpz_abs_fits_ui(c))
        {
            n_factor_t nf;
            n_factor_init(&nf);
            n_factor(&nf, fmpz_get_ui(c), 1);
            for (i = 0; i < nf.num; i++)
                _fmpz_factor_append_ui(fac, nf.p[i], nf.exp[i]);
            fmpz_one(c);
            complete = 1;
        }
    }
    else
        fmpz_one(c);

    if (cofactor != NULL)
        fmpz_swap(cofactor, c);
    else if (!fmpz_is_one(c))
        _fmpz_factor_append(fac, c, 1);

    fmpz_clear(c);
    return complete;
}

slong
gr_tower_degree_at(const gr_tower_t T, slong k)
{
    slong j, d = 1;

    /* saturating: the degree of a large tower is only compared against
       small limits */
    for (j = 1; j <= k; j++)
    {
        slong e = gr_tower_step_degree(T, j);
        if (d > WORD_MAX / e)
            return WORD_MAX;
        d *= e;
    }
    return d;
}

slong
gr_tower_degree(const gr_tower_t T)
{
    slong k, d = 1;

    /* saturating: the degree of a large tower is only compared against
       small limits */
    for (k = 1; k <= T->length; k++)
    {
        slong e = gr_tower_step_degree(T, k);
        if (d > WORD_MAX / e)
            return WORD_MAX;
        d *= e;
    }

    return d;
}

const gr_poly_struct *
gr_tower_step_minpoly(const gr_tower_t T, slong k)
{
    return gr_poly_quotient_ctx_modulus(GR_TOWER_STEP(T, k - 1)->ctx);
}

static int
_write_enclosure(gr_stream_t out, const acb_t z)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write_free(out, arb_get_str(acb_realref(z), 10, ARB_STR_NO_RADIUS));
    if (!arb_is_zero(acb_imagref(z)))
    {
        status |= gr_stream_write(out, " + ");
        status |= gr_stream_write_free(out, arb_get_str(acb_imagref(z), 10, ARB_STR_NO_RADIUS));
        status |= gr_stream_write(out, "*I");
    }
    return status;
}

int
gr_tower_write(gr_stream_t out, const gr_tower_t T)
{
    int status = GR_SUCCESS;
    slong d;

    status |= gr_stream_write(out, "Tower of degree ");
    status |= gr_stream_write_si(out, gr_tower_degree(T));
    status |= gr_stream_write(out, " over ");
    status |= gr_ctx_write(out, T->base);
    status |= gr_stream_write(out, "\n");

    for (d = 0; d < T->num_gens; d++)
    {
        const gr_tower_gen_struct * g = T->gens + d;

        status |= gr_stream_write(out, "  ");
        status |= gr_stream_write(out, g->name);
        status |= gr_stream_write(out, " = ");
        if (g->enclosure_prec > 0 || g->kind == GR_TOWER_ALGEBRAIC)
            status |= _write_enclosure(out, &g->enclosure);
        else
        {
            /* a transcendental generator not yet evaluated */
            acb_t z;
            acb_init(z);
            if (gr_tower_trans_get_acb(z, (gr_tower_struct *) T, g->index, GR_TOWER_DEFAULT_PREC) == GR_SUCCESS)
                status |= _write_enclosure(out, z);
            else
                status |= gr_stream_write(out, "?");
            acb_clear(z);
        }

        if (g->def_kind != GR_TOWER_ALGEBRAIC)
        {
            if (g->def_kind == GR_TOWER_ROOT)
            {
                status |= gr_stream_write(out, "  root_");
                status |= gr_stream_write_si(out, g->def_param);
                status |= gr_stream_write(out, "(");
            }
            else if (g->def_kind == GR_TOWER_ROOT_OF_UNITY)
            {
                status |= gr_stream_write(out, "  exp(2 pi i / ");
                status |= gr_stream_write_si(out, g->def_param);
                status |= gr_stream_write(out, ")");
            }
            else if (g->def_kind == GR_TOWER_TAN_PI)
            {
                status |= gr_stream_write(out, "  tan(pi / ");
                status |= gr_stream_write_si(out, g->def_param);
                status |= gr_stream_write(out, ")");
            }
            else if (GR_TOWER_KIND_IS_SPECIAL(g->def_kind))
            {
                char * s = _gr_tower_gen_args_str(g, (gr_tower_struct *) T);
                status |= gr_stream_write(out, "  ");
                status |= _gr_tower_special_write(out, g->def_kind, g->def_param, (s != NULL) ? s : "");
                flint_free(s);
            }
            else
                status |= gr_stream_write(out, (g->def_kind == GR_TOWER_EXP) ? "  exp(" : (g->def_kind == GR_TOWER_LOG) ? "  log(" :
                    (g->def_kind == GR_TOWER_TAN) ? "  tan(" : (g->def_kind == GR_TOWER_ATAN) ? "  atan(" : (g->def_kind == GR_TOWER_PI) ? "  pi" : "  free");
            if (g->arg.mctx != NULL && !GR_TOWER_KIND_IS_SPECIAL(g->def_kind))
            {
                char * s = _gr_tower_flat_get_str(&g->arg.data, g->arg.mctx, (gr_tower_struct *) T);
                status |= gr_stream_write_free(out, s);
                status |= gr_stream_write(out, ")");
            }
        }

        if (g->kind == GR_TOWER_ALGEBRAIC)
        {
            slong k = g->index;
            status |= gr_stream_write(out, (g->def_kind != GR_TOWER_ALGEBRAIC) ? "  =  root of  " : "  root of  ");
            status |= gr_poly_write(out, gr_poly_quotient_ctx_modulus(g->ctx), g->name, gr_tower_field_at(T, k - 1));
            status |= gr_stream_write(out, (g->status == GR_TOWER_STATUS_PROVEN) ? "  (irreducible)\n" : "  (dynamic)\n");
        }
        else
        {
            status |= gr_stream_write(out, (g->status == GR_TOWER_STATUS_PROVEN) ? "  (transcendental)\n" : "  (conjecturally transcendental)\n");
        }
    }

    return status;
}

void
gr_tower_print(const gr_tower_t T)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    GR_MUST_SUCCEED(gr_tower_write(out, T));
}

void
_gr_tower_move_to_front(gr_tower_t T, slong d)
{
    _gr_tower_move_gen(T, d, 0);
}

POP_OPTIONS
