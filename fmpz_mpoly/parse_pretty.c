/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fmpz_mpoly.h"
#include <assert.h>


int _fmpz_mpoly_parse_pretty_pop(fmpz_mpoly_geobucket_struct ** estack,
                                 slong * ostack, slong * _ei, slong * _oi,
                                           const fmpz_mpoly_ctx_t ctx, int all)
{
    int ret = 0;
    slong ei = *_ei;
    slong oi = *_oi;

    while (oi > 0 && ((all && (ostack[oi - 1] == '+' || ostack[oi - 1] == '-'))
                      || ostack[oi - 1] == '*'
                      || ostack[oi - 1] == '/'
                      || ostack[oi - 1] == 100 + '-'
                      || ostack[oi - 1] == 100 + '+'
                     )
          )
    {
        oi--;
        if (ostack[oi] == '+')
        {
            if (ei < 2)
                goto failed;
            fmpz_mpoly_geobucket_add_inplace(estack[ei - 2], estack[ei - 1], ctx);
            fmpz_mpoly_geobucket_clear(estack[--ei], ctx);
        } else if (ostack[oi] == '-')
        {
            if (ei < 2)
                goto failed;
            fmpz_mpoly_geobucket_sub_inplace(estack[ei - 2], estack[ei - 1], ctx);
            fmpz_mpoly_geobucket_clear(estack[--ei], ctx);
        } else if (ostack[oi] == '*')
        {
            if (ei < 2)
                goto failed;
            fmpz_mpoly_geobucket_mul_inplace(estack[ei - 2], estack[ei - 1], ctx);
            fmpz_mpoly_geobucket_clear(estack[--ei], ctx);
        } else if (ostack[oi] == '/')
        {
            if (ei < 2)
                goto failed;
            if (!fmpz_mpoly_geobucket_divides_inplace(estack[ei - 2], estack[ei - 1], ctx))
                goto failed;
            fmpz_mpoly_geobucket_clear(estack[--ei], ctx);
        } else if (ostack[oi] == 100 + '-')
        {
            if (ei < 1)
                goto failed;
            fmpz_mpoly_geobucket_neg_inplace(estack[ei - 1], ctx);
        } else {
            /* must be unary + */
            if (ei < 1)
                goto failed;
        }
    }

done:
    *_ei = ei;
    *_oi = oi;
    return ret;

failed:
    ret = -1;
    goto done;
}

void _fmpz_mpoly_parse_pretty_fit_estack(fmpz_mpoly_geobucket_struct *** estack,
                                                      slong ei, slong * ealloc)
{
    if (ei >= *ealloc)
    {
        slong new_ealloc = ei + 8;

        (* estack) = (fmpz_mpoly_geobucket_struct **) flint_realloc(* estack,
                              new_ealloc*sizeof(fmpz_mpoly_geobucket_struct*));
        for (; ei < new_ealloc; ei++)
        {
            (* estack)[ei] = (fmpz_mpoly_geobucket_struct *)
                             flint_malloc(sizeof(fmpz_mpoly_geobucket_struct));
        }
        *ealloc = new_ealloc;
    }
}

void _fmpz_mpoly_parse_pretty_fit_ostack(slong ** ostack,
                                                      slong oi, slong * oalloc)
{
    if (oi >= *oalloc)
    {
        slong new_oalloc = oi + 8;

        (*ostack) = (slong *) flint_realloc(*ostack, new_oalloc*sizeof(slong));
        *oalloc = new_oalloc;
    }
}

const char * _fmpz_mpoly_parse_pretty_int(const char * s, const char * end,
                                                           fmpz_t c, int * ret)
{
    char * buffer, * v;
    TMP_INIT;
    TMP_START;
    v = buffer = (char *) TMP_ALLOC((end - s + 1)*sizeof(char));
    while (s < end && '0' <= *s && *s <= '9')
        *v++ = *s++;
    *v++ = '\0';
    *ret = fmpz_set_str(c, buffer, 10);
    TMP_END;
    return s;
}

int _fmpz_mpoly_parse_pretty(fmpz_mpoly_t poly, const char * s, slong sn,
                                         char ** x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_geobucket_struct ** estack;
    slong * ostack;
    slong ealloc = 8, oalloc = 8;
    slong ei = 0, oi = 0;
    const char *end = s + sn;
    int expecting = 1;
    fmpz_t c;
    slong l, k;
    int ret;

    fmpz_init(c);
    ostack = (slong *) flint_malloc(ealloc*sizeof(slong));
    estack = (fmpz_mpoly_geobucket_struct **) flint_malloc(ealloc
                                       * sizeof(fmpz_mpoly_geobucket_struct*));
    for (k = 0; k < ealloc; k++)
    {
        estack[k] = (fmpz_mpoly_geobucket_struct*) flint_malloc(
                                          sizeof(fmpz_mpoly_geobucket_struct));
    }

    while (s < end)
    {
        if ('0' <= *s && *s <= '9')
        {
            s = _fmpz_mpoly_parse_pretty_int(s, end, c, &ret);
            if (!(expecting & 1) || ret)
                goto failed;
            _fmpz_mpoly_parse_pretty_fit_estack(&estack, ei, &ealloc);
            fmpz_mpoly_geobucket_init(estack[ei], ctx);
            fmpz_mpoly_geobucket_set_fmpz(estack[ei++], c, ctx);
            expecting = 2;

        } else if (*s == '^')
        {
            s = _fmpz_mpoly_parse_pretty_int(++s, end, c, &ret);
            if (!(expecting & 2) || ret)
                goto failed;
            fmpz_mpoly_geobucket_pow_fmpz_inplace(estack[ei - 1], c, ctx);
            expecting = 2;

        } else if ((*s == '+' || *s == '-') && (expecting & 2))
        {
            /* infix */
            if (_fmpz_mpoly_parse_pretty_pop(estack, ostack, &ei, &oi, ctx, 1))
                goto failed;
            _fmpz_mpoly_parse_pretty_fit_ostack(&ostack, oi, &oalloc);
            ostack[oi++] = *s++;
            expecting = 1;

        } else if ((*s == '+' || *s == '-') && (expecting & 1))
        {
            /* unary */
            _fmpz_mpoly_parse_pretty_fit_ostack(&ostack, oi, &oalloc);
            ostack[oi++] = 100 + *s++;
            expecting = 1;

        } else if (*s == '*')
        {
            if (!(expecting & 2)
                 || _fmpz_mpoly_parse_pretty_pop(estack, ostack, &ei, &oi, ctx, 0))
                goto failed;
            _fmpz_mpoly_parse_pretty_fit_ostack(&ostack, oi, &oalloc);
            ostack[oi++] = *s++;
            expecting = 1;

        } else if (*s == '/')
        {
            if (!(expecting & 2)
                 || _fmpz_mpoly_parse_pretty_pop(estack, ostack, &ei, &oi, ctx, 0))
                goto failed;
            _fmpz_mpoly_parse_pretty_fit_ostack(&ostack, oi, &oalloc);
            ostack[oi++] = *s++;
            expecting = 1;

        } else if (*s == ' ')
        {
            s++;

        } else if (*s == '(')
        {
            if (!(expecting & 1))
                goto failed;
            _fmpz_mpoly_parse_pretty_fit_ostack(&ostack, oi, &oalloc);
            ostack[oi++] = *s++;
            expecting = 1;

        } else if (*s == ')')
        {
            if (_fmpz_mpoly_parse_pretty_pop(estack, ostack, &ei, &oi, ctx, 1)
                     || oi < 1 || ostack[--oi] != '(')
                goto failed;
            s++;
            expecting = 2;

        } else {
            /* must be a variable */
            slong var = -WORD(1);
            slong matched_length = -WORD(1);
            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                l = strlen(x[k]);
                if ((end - s >= l) && (strncmp(s, x[k], l) == 0)
                                   && l > matched_length)
                {
                    var = k;
                    matched_length = l;
                }
            }

            if (!(expecting & 1) || var < 0)
                goto failed;

            _fmpz_mpoly_parse_pretty_fit_estack(&estack, ei, &ealloc);
            fmpz_mpoly_geobucket_init(estack[ei], ctx);
            fmpz_mpoly_geobucket_gen(estack[ei], var, ctx);
            ei++;
            s += matched_length;
            expecting = 2;
        }
    }

    if (_fmpz_mpoly_parse_pretty_pop(estack, ostack, &ei, &oi, ctx, 1)
               || ei != 1 || oi != 0)
        goto failed;

    ret = 0;
    fmpz_mpoly_geobucket_empty(poly, estack[0], ctx);
    fmpz_mpoly_geobucket_clear(estack[0], ctx);

done:
    for (k = 0; k < ealloc; k++)
        flint_free(estack[k]);
    flint_free(ostack);
    flint_free(estack);
    fmpz_clear(c);
    return ret;

failed:
    ret = -1;
    fmpz_mpoly_set_ui(poly, 0, ctx);
    for (k = 0; k < ei; k++)
        fmpz_mpoly_geobucket_clear(estack[k], ctx);
    goto done;
}


int fmpz_mpoly_set_str_pretty(fmpz_mpoly_t poly, const char * str,
                                 const char** x_in, const fmpz_mpoly_ctx_t ctx)
{

    int ret;
    slong i, nvars = ctx->minfo->nvars;
    char ** x = (char **) x_in;
    TMP_INIT;

    TMP_START;
    if (x == NULL)
    {
        x = (char **) TMP_ALLOC(nvars*sizeof(char *));
        for (i = 0; i < nvars; i++)
        {
            x[i] = (char *) TMP_ALLOC(22*sizeof(char));
            flint_sprintf(x[i], "x%wd", i + 1);
        }
    }
    ret = _fmpz_mpoly_parse_pretty(poly, str, strlen(str), x, ctx);
    TMP_END;
    return ret;
}
