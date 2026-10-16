/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"
#include "gr.h"

typedef struct
{
    char * s;
    size_t len;
    size_t alloc;
}
_strbuf;

static void
_strbuf_append(_strbuf * b, const char * t)
{
    size_t n = strlen(t);

    if (b->len + n + 1 > b->alloc)
    {
        b->alloc = FLINT_MAX(2 * b->alloc, b->len + n + 1);
        b->s = flint_realloc(b->s, b->alloc);
    }

    memcpy(b->s + b->len, t, n + 1);
    b->len += n;
}

static void
_strbuf_append_free(_strbuf * b, char * t)
{
    _strbuf_append(b, t);
    flint_free(t);
}

char *
fmpz_mpoly_q_factored_get_str_pretty(const fmpz_mpoly_q_factored_t x, const char ** vars, const fmpz_mpoly_ctx_t ctx)
{
    const fmpz_mpoly_factor_struct * den = fmpz_mpoly_q_factored_denref(x);
    _strbuf b;
    slong i;
    int need_star;

    b.alloc = 16;
    b.len = 0;
    b.s = flint_malloc(b.alloc);
    b.s[0] = '\0';

    if (den->num == 0 && fmpz_is_one(den->constant))
    {
        _strbuf_append_free(&b, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_factored_numref(x), vars, ctx));
        return b.s;
    }

    _strbuf_append(&b, "(");
    _strbuf_append_free(&b, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_factored_numref(x), vars, ctx));
    _strbuf_append(&b, ")/");

    if (den->num == 0)
    {
        _strbuf_append_free(&b, fmpz_get_str(NULL, 10, den->constant));
        return b.s;
    }

    _strbuf_append(&b, "(");

    need_star = 0;
    if (!fmpz_is_one(den->constant))
    {
        _strbuf_append_free(&b, fmpz_get_str(NULL, 10, den->constant));
        need_star = 1;
    }

    for (i = 0; i < den->num; i++)
    {
        if (need_star)
            _strbuf_append(&b, "*");

        _strbuf_append(&b, "(");
        _strbuf_append_free(&b, fmpz_mpoly_get_str_pretty(den->poly + i, vars, ctx));
        _strbuf_append(&b, ")");

        if (!fmpz_is_one(den->exp + i))
        {
            _strbuf_append(&b, "^");
            _strbuf_append_free(&b, fmpz_get_str(NULL, 10, den->exp + i));
        }

        need_star = 1;
    }

    _strbuf_append(&b, ")");

    return b.s;
}

void
fmpz_mpoly_q_factored_print_pretty(const fmpz_mpoly_q_factored_t x, const char ** vars, const fmpz_mpoly_ctx_t ctx)
{
    char * s = fmpz_mpoly_q_factored_get_str_pretty(x, vars, ctx);
    flint_printf("%s", s);
    flint_free(s);
}

int
fmpz_mpoly_q_factored_set_str_pretty(fmpz_mpoly_q_factored_t res, const char * s, const char ** vars, const fmpz_mpoly_ctx_t ctx)
{
    gr_ctx_t grctx;
    int ret;

    gr_ctx_init_fmpz_mpoly_q_factored(grctx, ctx->minfo->nvars, ctx->minfo->ord);
    if (vars != NULL)
        GR_MUST_SUCCEED(gr_ctx_set_gen_names(grctx, vars));
    ret = (GR_SUCCESS == gr_set_str(res, s, grctx)) ? 0 : -1;
    gr_ctx_clear(grctx);

    if (ret != 0)
        fmpz_mpoly_q_factored_zero(res, ctx);

    return ret;
}
