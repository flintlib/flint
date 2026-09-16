/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include "fmpz.h"
#include "fmpz_mod.h"
#include "ecpp.h"

void
ecpp_cert_init(ecpp_cert_t cert)
{
    cert->steps = NULL;
    cert->num = 0;
    cert->alloc = 0;
}

void
ecpp_cert_clear(ecpp_cert_t cert)
{
    slong i;

    for (i = 0; i < cert->alloc; i++)
    {
        ecpp_step_struct * s = cert->steps + i;
        fmpz_clear(s->n);
        fmpz_clear(s->a);
        fmpz_clear(s->b);
        fmpz_clear(s->m);
        fmpz_clear(s->q);
        fmpz_clear(s->x);
        fmpz_clear(s->y);
    }

    flint_free(cert->steps);
}

/* append an (initialised, zero) step and return it */
ecpp_step_struct *
ecpp_cert_push(ecpp_cert_t cert)
{
    ecpp_step_struct * s;

    if (cert->num == cert->alloc)
    {
        slong i, new_alloc = FLINT_MAX(2 * cert->alloc, 8);

        cert->steps = flint_realloc(cert->steps, new_alloc * sizeof(ecpp_step_struct));
        for (i = cert->alloc; i < new_alloc; i++)
        {
            s = cert->steps + i;
            fmpz_init(s->n);
            fmpz_init(s->a);
            fmpz_init(s->b);
            fmpz_init(s->m);
            fmpz_init(s->q);
            fmpz_init(s->x);
            fmpz_init(s->y);
            s->D = 0;
        }
        cert->alloc = new_alloc;
    }

    return cert->steps + cert->num++;
}

void
ecpp_cert_pop(ecpp_cert_t cert)
{
    if (cert->num > 0)
        cert->num--;
}

/*
    Text formats. ECPP_CERT_FORMAT_FLINT:

        ecpp certificate, k steps
        [0] n = ...
            D = ...
            a = ...
            b = ...
            m = ...
            q = ...
            x = ...
            y = ...
        ...

    ECPP_CERT_FORMAT_PARI, the primecert format of PARI/GP: a vector of
    [N, t, s, a, [x, y]] with N + 1 - t = s q the order of the curve
    y^2 = x^3 + a x + b (b = y^2 - x^3 - a x mod N) and q the next N.
*/

static void
_cert_write(char ** str, slong * len, slong * alloc, const char * piece)
{
    slong l = strlen(piece);
    if (*len + l + 1 > *alloc)
    {
        *alloc = FLINT_MAX(2 * *alloc, *len + l + 1);
        *str = flint_realloc(*str, *alloc);
    }
    memcpy(*str + *len, piece, l + 1);
    *len += l;
}

static void
_cert_write_fmpz(char ** str, slong * len, slong * alloc, const fmpz_t x)
{
    char * t = fmpz_get_str(NULL, 10, x);
    _cert_write(str, len, alloc, t);
    flint_free(t);
}

char *
ecpp_cert_get_str(const ecpp_cert_t cert, int format)
{
    slong i, len = 0, alloc = 256;
    char * str = flint_malloc(alloc);
    char buf[64];
    str[0] = '\0';

    if (format == ECPP_CERT_FORMAT_PARI)
    {
        fmpz_t t, s;
        fmpz_init(t);
        fmpz_init(s);
        _cert_write(&str, &len, &alloc, "[");
        for (i = 0; i < cert->num; i++)
        {
            const ecpp_step_struct * st = cert->steps + i;
            fmpz_add_ui(t, st->n, 1);
            fmpz_sub(t, t, st->m);              /* t = n + 1 - m */
            fmpz_divexact(s, st->m, st->q);     /* s = m / q */
            _cert_write(&str, &len, &alloc, i > 0 ? ", [" : "[");
            _cert_write_fmpz(&str, &len, &alloc, st->n);
            _cert_write(&str, &len, &alloc, ", ");
            _cert_write_fmpz(&str, &len, &alloc, t);
            _cert_write(&str, &len, &alloc, ", ");
            _cert_write_fmpz(&str, &len, &alloc, s);
            _cert_write(&str, &len, &alloc, ", ");
            _cert_write_fmpz(&str, &len, &alloc, st->a);
            _cert_write(&str, &len, &alloc, ", [");
            _cert_write_fmpz(&str, &len, &alloc, st->x);
            _cert_write(&str, &len, &alloc, ", ");
            _cert_write_fmpz(&str, &len, &alloc, st->y);
            _cert_write(&str, &len, &alloc, "]]");
        }
        _cert_write(&str, &len, &alloc, "]");
        fmpz_clear(t);
        fmpz_clear(s);
    }
    else
    {
        flint_sprintf(buf, "ecpp certificate, %wd steps\n", cert->num);
        _cert_write(&str, &len, &alloc, buf);
        for (i = 0; i < cert->num; i++)
        {
            const ecpp_step_struct * st = cert->steps + i;
            flint_sprintf(buf, "[%wd] n = ", i);
            _cert_write(&str, &len, &alloc, buf);
            _cert_write_fmpz(&str, &len, &alloc, st->n);
            flint_sprintf(buf, "\n    D = %wd\n    a = ", st->D);
            _cert_write(&str, &len, &alloc, buf);
            _cert_write_fmpz(&str, &len, &alloc, st->a);
            _cert_write(&str, &len, &alloc, "\n    b = ");
            _cert_write_fmpz(&str, &len, &alloc, st->b);
            _cert_write(&str, &len, &alloc, "\n    m = ");
            _cert_write_fmpz(&str, &len, &alloc, st->m);
            _cert_write(&str, &len, &alloc, "\n    q = ");
            _cert_write_fmpz(&str, &len, &alloc, st->q);
            _cert_write(&str, &len, &alloc, "\n    x = ");
            _cert_write_fmpz(&str, &len, &alloc, st->x);
            _cert_write(&str, &len, &alloc, "\n    y = ");
            _cert_write_fmpz(&str, &len, &alloc, st->y);
            _cert_write(&str, &len, &alloc, "\n");
        }
    }
    return str;
}

void
ecpp_cert_fprint(FILE * file, const ecpp_cert_t cert, int format)
{
    char * str = ecpp_cert_get_str(cert, format);
    fputs(str, file);
    if (format == ECPP_CERT_FORMAT_PARI)
        fputc('\n', file);
    flint_free(str);
}

void
ecpp_cert_print(const ecpp_cert_t cert, int format)
{
    ecpp_cert_fprint(stdout, cert, format);
}

/* reads the decimal integer at *p (optional sign), advances *p; 0 on failure */
static int
_read_int(fmpz_t x, const char ** p)
{
    const char * q;
    char * buf;
    slong l;
    int ok;
    while (**p == ' ' || **p == '\n' || **p == '\t' || **p == '\r')
        (*p)++;
    q = *p;
    if (*q == '-' || *q == '+')
        q++;
    if (*q < '0' || *q > '9')
        return 0;
    while (*q >= '0' && *q <= '9')
        q++;
    l = q - *p;
    buf = flint_malloc(l + 1);
    memcpy(buf, *p, l);
    buf[l] = '\0';
    ok = (fmpz_set_str(x, buf, 10) == 0);
    flint_free(buf);
    *p = q;
    return ok;
}

/* skips whitespace and the expected character; 0 if absent */
static int
_expect(const char ** p, char c)
{
    while (**p == ' ' || **p == '\n' || **p == '\t' || **p == '\r')
        (*p)++;
    if (**p != c)
        return 0;
    (*p)++;
    return 1;
}

int
ecpp_cert_set_str(ecpp_cert_t cert, const char * str)
{
    const char * p = str;
    int ok = 1;

    cert->num = 0;
    while (*p == ' ' || *p == '\n' || *p == '\t' || *p == '\r')
        p++;

    if (*p == '[')
    {
        /* PARI/GP: [[N, t, s, a, [x, y]], ...] */
        fmpz_t t, s, u;
        fmpz_mod_ctx_t ctx;
        fmpz_init(t); fmpz_init(s); fmpz_init(u);
        p++;
        while (ok)
        {
            ecpp_step_struct * st;
            while (*p == ' ' || *p == '\n' || *p == '\t' || *p == '\r')
                p++;
            if (*p == ']')
            {
                p++;
                break;
            }
            if (cert->num > 0 && !_expect(&p, ','))
            {
                ok = 0;
                break;
            }
            if (!_expect(&p, '['))
            {
                ok = 0;
                break;
            }
            st = ecpp_cert_push(cert);
            ok = _read_int(st->n, &p) && _expect(&p, ',') && _read_int(t, &p) && _expect(&p, ',')
              && _read_int(s, &p) && _expect(&p, ',') && _read_int(st->a, &p) && _expect(&p, ',')
              && _expect(&p, '[') && _read_int(st->x, &p) && _expect(&p, ',') && _read_int(st->y, &p)
              && _expect(&p, ']') && _expect(&p, ']');
            if (!ok || fmpz_cmp_ui(st->n, 3) < 0 || fmpz_sgn(s) <= 0)
            {
                ok = 0;
                break;
            }
            /* m = n + 1 - t, q = m / s (must be exact), b = y^2 - x^3 - a x */
            fmpz_add_ui(st->m, st->n, 1);
            fmpz_sub(st->m, st->m, t);
            fmpz_fdiv_qr(st->q, u, st->m, s);
            if (!fmpz_is_zero(u))
            {
                ok = 0;
                break;
            }
            fmpz_mod_ctx_init(ctx, st->n);
            fmpz_mod_set_fmpz(st->a, st->a, ctx);
            fmpz_mod_set_fmpz(st->x, st->x, ctx);
            fmpz_mod_set_fmpz(st->y, st->y, ctx);
            fmpz_mod_mul(u, st->x, st->x, ctx);
            fmpz_mod_mul(u, u, st->x, ctx);         /* x^3 */
            fmpz_mod_mul(st->b, st->y, st->y, ctx);
            fmpz_mod_sub(st->b, st->b, u, ctx);
            fmpz_mod_mul(u, st->a, st->x, ctx);
            fmpz_mod_sub(st->b, st->b, u, ctx);
            fmpz_mod_ctx_clear(ctx);
            st->D = 0;                              /* not part of the format */
        }
        fmpz_clear(t); fmpz_clear(s); fmpz_clear(u);
    }
    else
    {
        /* FLINT: "ecpp certificate, k steps" then the labelled steps */
        const char * labels[8] = {"n =", "D =", "a =", "b =", "m =", "q =", "x =", "y ="};
        fmpz_t Dz;
        fmpz_init(Dz);
        p = strstr(p, "[0]");
        while (ok && p != NULL)
        {
            ecpp_step_struct * st = ecpp_cert_push(cert);
            fmpz * fields[8] = {st->n, Dz, st->a, st->b, st->m, st->q, st->x, st->y};
            slong f;
            for (f = 0; f < 8 && ok; f++)
            {
                const char * q = strstr(p, labels[f]);
                if (q == NULL)
                {
                    ok = 0;
                    break;
                }
                q += 3;
                ok = _read_int(fields[f], &q);
                p = q;
            }
            if (ok)
                st->D = fmpz_get_si(Dz);
            p = strstr(p, "\n[");
        }
        fmpz_clear(Dz);
    }

    if (!ok || cert->num == 0)
    {
        cert->num = 0;
        return 0;
    }
    return 1;
}
