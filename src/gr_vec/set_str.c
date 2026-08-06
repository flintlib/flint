/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "gr_vec.h"

int
_gr_str_is_space(char c)
{
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

/* Scan [s, end) for the next top-level ',' (at bracket depth 0), treating
   commas nested inside (), [] or {} as part of an entry. Returns a pointer to
   the separating comma, or end if there is none. Sets *err on unbalanced
   delimiters. */
const char *
_gr_str_find_sep(const char * s, const char * end, int * err)
{
    slong depth = 0;

    for (; s < end; s++)
    {
        char c = *s;

        if (c == '(' || c == '[' || c == '{')
            depth++;
        else if (c == ')' || c == ']' || c == '}')
        {
            if (depth == 0)
            {
                *err = 1;
                return end;
            }
            depth--;
        }
        else if (c == ',' && depth == 0)
            return s;
    }

    if (depth != 0)
        *err = 1;

    return end;
}

/* Locate the content strictly between a leading 'open' bracket and its matching
   'close' bracket. Characters following the closing bracket are ignored. */
static int
_gr_str_bracket_content(const char ** content_begin, const char ** content_end,
        const char * s, char open, char close)
{
    const char * end = s + strlen(s);
    const char * p;
    slong depth;

    while (s < end && _gr_str_is_space(*s))
        s++;

    if (s >= end || *s != open)
        return GR_UNABLE;

    depth = 0;
    for (p = s; p < end; p++)
    {
        char c = *p;
        if (c == '(' || c == '[' || c == '{')
            depth++;
        else if (c == ')' || c == ']' || c == '}')
        {
            depth--;
            if (depth == 0)
                break;
        }
    }

    if (p >= end || *p != close)
        return GR_UNABLE;

    *content_begin = s + 1;
    *content_end = p;
    return GR_SUCCESS;
}

/* As above, but only whitespace may follow the closing bracket. */
int
_gr_str_strip_brackets(const char ** content_begin, const char ** content_end,
        const char * s, char open, char close)
{
    const char * q;
    int status = _gr_str_bracket_content(content_begin, content_end, s, open, close);

    if (status != GR_SUCCESS)
        return status;

    for (q = *content_end + 1; *q != '\0'; q++)
        if (!_gr_str_is_space(*q))
            return GR_UNABLE;

    return GR_SUCCESS;
}

/* True if [s, end) is empty or all whitespace. */
int
_gr_str_blank(const char * s, const char * end)
{
    for (; s < end; s++)
        if (!_gr_str_is_space(*s))
            return 0;
    return 1;
}

/* Trim leading and trailing whitespace from [*s, *end). */
void
_gr_str_trim(const char ** s, const char ** end)
{
    const char * a = *s;
    const char * b = *end;

    while (a < b && _gr_str_is_space(*a))
        a++;
    while (b > a && _gr_str_is_space(b[-1]))
        b--;

    *s = a;
    *end = b;
}

/* Count top-level comma-separated entries in [lo, hi); 0 if blank. */
static int
_gr_str_count_in_range(slong * count, const char * lo, const char * hi)
{
    slong n;
    const char * seg;

    if (_gr_str_blank(lo, hi))
    {
        *count = 0;
        return GR_SUCCESS;
    }

    n = 0;
    seg = lo;
    while (1)
    {
        int err = 0;
        const char * sep = _gr_str_find_sep(seg, hi, &err);
        if (err)
            return GR_UNABLE;
        n++;
        if (sep == hi)
            break;
        seg = sep + 1;
    }

    *count = n;
    return GR_SUCCESS;
}

int
gr_vec_str_count_entries(slong * count, const char * s, gr_ctx_t FLINT_UNUSED(ctx))
{
    const char * lo;
    const char * hi;
    int status = _gr_str_bracket_content(&lo, &hi, s, '[', ']');

    if (status != GR_SUCCESS)
        return status;

    return _gr_str_count_in_range(count, lo, hi);
}

int
_gr_vec_set_str(gr_ptr res, const char * s, slong len, gr_ctx_t ctx)
{
    const char * lo;
    const char * hi;
    const char * seg;
    char * buf;
    slong i, sz = ctx->sizeof_elem;
    int status;

    status = _gr_str_bracket_content(&lo, &hi, s, '[', ']');
    if (status != GR_SUCCESS)
        return status;

    if (_gr_str_blank(lo, hi))
        return (len == 0) ? GR_SUCCESS : GR_DOMAIN;

    buf = flint_malloc((slong) (hi - lo) + 1);

    status = GR_SUCCESS;
    seg = lo;
    i = 0;
    while (1)
    {
        int err = 0;
        const char * sep = _gr_str_find_sep(seg, hi, &err);
        const char * a = seg;
        const char * b = sep;
        slong elen;

        if (err)
        {
            status = GR_UNABLE;
            break;
        }

        if (i >= len)   /* more entries than expected */
        {
            status = GR_DOMAIN;
            break;
        }

        _gr_str_trim(&a, &b);
        elen = (slong) (b - a);
        if (elen <= 0)   /* empty entry */
        {
            status = GR_UNABLE;
            break;
        }

        memcpy(buf, a, elen);
        buf[elen] = '\0';

        status = gr_set_str(GR_ENTRY(res, i, sz), buf, ctx);
        if (status != GR_SUCCESS)
            break;

        i++;
        if (sep == hi)
            break;
        seg = sep + 1;
    }

    flint_free(buf);

    if (status == GR_SUCCESS && i != len)   /* fewer entries than expected */
        status = GR_DOMAIN;

    return status;
}

int
gr_vec_set_str(gr_vec_t vec, const char * s, int resize, gr_ctx_t ctx)
{
    const char * lo;
    const char * hi;
    slong count;
    int status;

    /* Strict validation: a single bracketed list with no trailing junk. */
    status = _gr_str_strip_brackets(&lo, &hi, s, '[', ']');
    if (status != GR_SUCCESS)
        return status;

    status = _gr_str_count_in_range(&count, lo, hi);
    if (status != GR_SUCCESS)
        return status;

    if (resize)
        gr_vec_set_length(vec, count, ctx);
    else if (gr_vec_length(vec, ctx) != count)
        return GR_DOMAIN;

    status = _gr_vec_set_str(vec->entries, s, count, ctx);

    if (status != GR_SUCCESS && resize)
        gr_vec_set_length(vec, 0, ctx);

    return status;
}
