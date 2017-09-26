/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include <assert.h>


void _fmpz_mpoly_parse_pretty_pop(fmpz_mpoly_struct ** estack,
                        slong * ostack, slong * _ei, slong * _oi,
                                    const fmpz_mpoly_ctx_t ctx, int all) {
    slong ei = *_ei;
    slong oi = *_oi;

    while (oi > 0 && ((all && (ostack[oi - 1] == '+' || ostack[oi - 1] == '-'))
           || ostack[oi - 1] == '*' || ostack[oi - 1] == 100 + '-'
                                    || ostack[oi - 1] == 100 + '+')) {
        oi--;
        if (ostack[oi] == '+') {
            --ei;
            fmpz_mpoly_add(estack[ei - 1], estack[ei - 1], estack[ei], ctx);
            fmpz_mpoly_clear(estack[ei], ctx);
        } else if (ostack[oi] == '-') {
            --ei;
            fmpz_mpoly_sub(estack[ei - 1], estack[ei - 1], estack[ei], ctx);
            fmpz_mpoly_clear(estack[--ei], ctx);
        } else if (ostack[oi] == '*') {
            --ei;
            fmpz_mpoly_mul_johnson(estack[ei - 1], estack[ei - 1], estack[ei], ctx);
            fmpz_mpoly_clear(estack[ei], ctx);
        } else if (ostack[oi] == 100 + '-') {
            fmpz_mpoly_neg(estack[ei - 1], estack[ei - 1], ctx);
        }
    }

    *_ei = ei;
    *_oi = oi;
}

void _fmpz_mpoly_parse_pretty_fit_estack(fmpz_mpoly_struct *** estack,
                                                      slong ei, slong * ealloc)
{
    if (ei >= *ealloc) {
        slong new_ealloc = ei + 8;

        (*estack) = (fmpz_mpoly_struct **) flint_realloc(*estack,
                                        new_ealloc*sizeof(fmpz_mpoly_struct*));
        for (; ei < new_ealloc; ei++) {
            (*estack)[ei] = (fmpz_mpoly_struct*)
                                       flint_malloc(sizeof(fmpz_mpoly_struct));            
        }
        *ealloc = new_ealloc;
    }
}

void _fmpz_mpoly_parse_pretty_fit_ostack(slong ** ostack,
                                                      slong oi, slong * oalloc)
{
    if (oi >= *oalloc) {
        slong new_oalloc = oi + 8;

        (*ostack) = (slong*) flint_realloc(*ostack, new_oalloc*sizeof(slong));
        *oalloc = new_oalloc;
    }
}


int _fmpz_mpoly_parse_pretty(fmpz_mpoly_t poly, const char * s, slong sn,
                                         char ** x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_struct ** estack;
    slong * ostack;
    slong ealloc, oalloc;
    slong ei, oi;
    const char *end = s + sn;
    int expecting = 1;
    fmpz_t c;
    slong l, k;
    int deg, rev, ret;
    slong nvars;

    degrev_from_ord(deg, rev, ctx->ord);
    nvars = ctx->n - deg;

    ealloc = 8;
    estack = (fmpz_mpoly_struct **) flint_malloc(ealloc*sizeof(fmpz_mpoly_struct*));
    for (k = 0; k < ealloc; k++) {
        estack[k] = (fmpz_mpoly_struct*) flint_malloc(sizeof(fmpz_mpoly_struct));
    }
    oalloc = 8;
    ostack = (slong*) flint_malloc(ealloc*sizeof(slong));

    ei = 0;
    oi = 0;

    fmpz_init(c);

get_next_char:

    if ('0' <= *s && *s <= '9') {
        if (!(expecting & 1)) {
            goto failed;
        }
        fmpz_set_ui(c, *s++ - '0');
        while (s < end && '0' <= *s && *s <= '9') {
            fmpz_mul_ui(c, c, WORD(10));
            fmpz_add_ui(c, c, *s++ - '0');
        }
        _fmpz_mpoly_parse_pretty_fit_estack(&estack, ei, &ealloc);
        fmpz_mpoly_init(estack[ei], ctx);
        fmpz_mpoly_set_fmpz(estack[ei], c, ctx);
        ei++;
        expecting = 2;

    } else if (*s == '^') {
        s++;
        if (!(expecting & 2)) {
            goto failed;
        }
        if (s < end && '0' <= *s && *s <= '9') {
            k = *s++ - '0';
            while (s < end && '0' <= *s && *s <= '9') {
                k = 10*k + *s++ - '0';
                if (k < 0) {
                    goto failed;
                }
            }
            fmpz_mpoly_pow_fps(estack[ei - 1], estack[ei - 1], k, ctx);
            expecting = 2;
        } else {
            goto failed;
        }

    } else if ((*s == '+' || *s == '-') && (expecting & 2)) {
        /* infix */
        _fmpz_mpoly_parse_pretty_pop(estack, ostack, &ei, &oi, ctx, 1);
        _fmpz_mpoly_parse_pretty_fit_ostack(&ostack, oi, &oalloc);
        ostack[oi++] = *s++;
        expecting = 1;

    } else if ((*s == '+' || *s == '-') && !(expecting & 2)) {
        /* unary */
        _fmpz_mpoly_parse_pretty_fit_ostack(&ostack, oi, &oalloc);
        ostack[oi++] = 100 + *s++;
        expecting = 1;

    } else if (*s == '*') {
        if (!(expecting & 2)) {
            goto failed;
        }
        _fmpz_mpoly_parse_pretty_pop(estack, ostack, &ei, &oi, ctx, 0);
        _fmpz_mpoly_parse_pretty_fit_ostack(&ostack, oi, &oalloc);
        ostack[oi++] = *s++;
        expecting = 1;

    } else if (*s == ' ') {
        s++;

    } else if (*s == '(') {
        if (!(expecting & 1)) {
            goto failed;
        }
        _fmpz_mpoly_parse_pretty_fit_ostack(&ostack, oi, &oalloc);
        ostack[oi++] = *s++;
        expecting = 1;
        
    } else if (*s == ')') {
        _fmpz_mpoly_parse_pretty_pop(estack, ostack, &ei, &oi, ctx, 1);
        s++;
        if (oi>0 && ostack[oi-1] == '(') {
            oi--;
        } else {
            goto failed;
        }
        expecting = 2;

    } else {
        /* must be a variable */
        if (!(expecting & 1)) {
            goto failed;
        }
        for (k = 0; k < nvars; k++) {
            l = strlen(x[k]);
            if ((end - s >= l) && (strncmp(s, x[k], l) == 0)) {
                break;
            }
        }
        if (k >= nvars) {
            goto failed;
        }
        _fmpz_mpoly_parse_pretty_fit_estack(&estack, ei, &ealloc);
        fmpz_mpoly_init(estack[ei], ctx);
        fmpz_mpoly_gen(estack[ei], k, ctx);
        ei++;
        s += l;
        expecting = 2;
    }
    if (s < end) {
        goto get_next_char;
    }

    _fmpz_mpoly_parse_pretty_pop(estack, ostack, &ei, &oi, ctx, 1);
    if (ei != 1 || oi != 0) {
        goto failed;
    }

    fmpz_mpoly_swap(poly, estack[0], ctx);
    fmpz_mpoly_clear(estack[0], ctx);

    ret = 0;

done:
    for (k = 0; k < ealloc; k++) {
        flint_free(estack[k]);
    }

    fmpz_clear(c);
    flint_free(ostack);
    flint_free(estack);
    return ret;


failed:
    ret = 1;
    fmpz_mpoly_set_ui(poly, 0, ctx);
    goto done;
}


int fmpz_mpoly_set_str_pretty(fmpz_mpoly_t poly, const char * str,
                                 const char** x_in, const fmpz_mpoly_ctx_t ctx)
{

    int ret, deg, rev;
    slong i, nvars;
    char ** x = (char **) x_in;
    TMP_INIT;

    TMP_START;
    degrev_from_ord(deg, rev, ctx->ord);
    nvars = ctx->n - deg;
    if (x == NULL) {
        x = (char **) TMP_ALLOC(nvars*sizeof(char *));
        for (i = 0; i < nvars; i++) {
            x[i] = (char *) TMP_ALLOC(22*sizeof(char));
            flint_sprintf(x[i], "x%wd", i + 1);
        }
    }
    ret = _fmpz_mpoly_parse_pretty(poly, str, strlen(str), x, ctx);
    TMP_END;
    return ret;
}
