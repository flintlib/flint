/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <ctype.h>
#include "arb.h"

#define RADIUS_DIGITS 3

char *
_arb_condense_digits(char * s, slong n)
{
    slong i, j, run, out;
    char * res;

    res = flint_malloc(strlen(s) + 128); /* space for some growth */
    out = 0;

    for (i = 0; s[i] != '\0'; )
    {
        if (isdigit(s[i]))
        {
            run = 0;

            for (j = 0; isdigit(s[i + j]); j++)
                run++;

            if (run > 3 * n)
            {
                for (j = 0; j < n; j++)
                {
                    res[out] = s[i + j];
                    out++;
                }

                out += flint_sprintf(res + out, "{...%wd digits...}", run - 2 * n);

                for (j = run - n; j < run; j++)
                {
                    res[out] = s[i + j];
                    out++;
                }
            }
            else
            {
                for (j = 0; j < run; j++)
                {
                    res[out] = s[i + j];
                    out++;
                }
            }

            i += run;
        }
        else
        {
            res[out] = s[i];
            i++;
            out++;
        }
    }

    res[out] = '\0';
    res = flint_realloc(res, strlen(res) + 1);

    flint_free(s);
    return res;
}

/* Format (digits=d, exponent=e) as floating-point or fixed-point.
   Reallocates the input and mutates the exponent. */
void
_arb_digits_as_float_str(char ** d, fmpz_t e, slong minfix, slong maxfix)
{
    slong i, n, alloc, dotpos;

    /* do nothing with 0 or something non-numerical */
    if (!((*d)[0] >= '1' && (*d)[0] <= '9'))
        return;

    n = strlen(*d);

    fmpz_add_ui(e, e, n - 1);

    /* fixed-point or integer format */
    /* we require e < n - 1; otherwise we would have to insert trailing zeros
      [todo: could allow e < n, if printing integers without radix point] */
    if (fmpz_cmp_si(e, minfix) >= 0 && fmpz_cmp_si(e, maxfix) <= 0 &&
        fmpz_cmp_si(e, n - 1) < 0)
    {
        slong exp = *e;

        /* 0.000xxx */
        if (exp < 0)
        {
            /* 0. + (-1-exp) zeros + digits + null terminator */
            alloc = 2 + (-1-exp) + n + 1;

            *d = flint_realloc(*d, alloc);

            /* copy in reverse order, including null terminator */
            for (i = n; i >= 0; i--)
                (*d)[2 + (-1-exp) + i] = (*d)[i];

            for (i = 0; i < 2 + (-1-exp); i++)
                (*d)[i] = (i == 1) ? '.' : '0';
        }
        else    /* xxx.yyy --- must have dotpos < n - 1 */
        {
            dotpos = exp + 1;
            alloc = n + 2;  /* space for . and null terminator */

            (*d) = flint_realloc(*d, alloc);

            /* copy fractional part in reverse order, including null */
            for (i = n; i >= dotpos; i--)
                (*d)[i + 1] = (*d)[i];

            (*d)[dotpos] = '.';
        }
    }
    else
    {
        /* format as xe+zzz or x.yyye+zzz */
        alloc = n + 1 + 2 + fmpz_sizeinbase(e, 10) + 1;
        *d = flint_realloc(*d, alloc);

        /* insert . */
        if (n > 1)
        {
            /* copy fractional part in reverse order */
            for (i = n; i >= 1; i--)
                (*d)[i + 1] = (*d)[i];

            (*d)[1] = '.';
        }

        (*d)[n + (n > 1)] = 'e';

        if (fmpz_sgn(e) >= 0)
        {
            (*d)[n + (n > 1) + 1] = '+';
        }
        else
        {
            (*d)[n + (n > 1) + 1] = '-';
            fmpz_neg(e, e);
        }

        fmpz_get_str((*d) + n + (n > 1) + 2, 10, e);  /* writes null byte */
    }
}

/* Rounds a string of decimal digits (null-terminated).
   to length at most n. The rounding mode
   can be ARF_RND_DOWN, ARF_RND_UP or ARF_RND_NEAR.
   The string is overwritten in-place, truncating it as necessary.
   The input should not have a leading sign or leading zero digits,
   but can have trailing zero digits.

   Computes shift and error such that

   int(input) = int(output) * 10^shift + error

   exactly.
*/
void
_arb_digits_round_inplace(char * s, flint_bitcnt_t * shift, fmpz_t error, slong n, arf_rnd_t rnd)
{
    slong i, m;
    int up;

    if (n < 1)
    {
        flint_throw(FLINT_ERROR, "_arb_digits_round_inplace: require n >= 1\n");
    }

    m = strlen(s);

    if (m <= n)
    {
        *shift = 0;
        fmpz_zero(error);
        return;
    }

    /* always round down */
    if (rnd == ARF_RND_DOWN)
    {
        up = 0;
    }
    else if (rnd == ARF_RND_UP) /* round up if tail is nonzero */
    {
        up = 0;

        for (i = n; i < m; i++)
        {
            if (s[i] != '0')
            {
                up = 1;
                break;
            }
        }
    }
    else /* round to nearest (up on tie -- todo: round-to-even?) */
    {
        up = (s[n] >= '5' && s[n] <= '9');
    }

    if (!up)
    {
        /* simply truncate */
        fmpz_set_str(error, s + n, 10);
        s[n] = '\0';
        *shift = m - n;
    }
    else
    {
        int digit, borrow, carry;

        /* error = 10^(m-n) - s[n:], where s[n:] is nonzero */
        /* i.e. 10s complement the truncated digits */
        borrow = 0;

        for (i = m - 1; i >= n; i--)
        {
            digit = 10 - (s[i] - '0') - borrow;

            if (digit == 10)
            {
                digit = 0;
                borrow = 0;
            }
            else
            {
                borrow = 1;
            }

            s[i] = digit + '0';
        }

        if (!borrow)
        {
            flint_throw(FLINT_ERROR, "expected borrow!\n");
        }

        fmpz_set_str(error, s + n, 10);
        fmpz_neg(error, error);

        /* add 1 ulp to the leading digits */
        carry = 1;

        for (i = n - 1; i >= 0; i--)
        {
            digit = (s[i] - '0') + carry;

            if (digit > 9)
            {
                digit = 0;
                carry = 1;
            }
            else
            {
                carry = 0;
            }

            s[i] = digit + '0';
        }

        /* carry-out -- only possible if we started with all 9s,
           so now the rest will be 0s which we don't have to shift explicitly */
        if (carry)
        {
            s[0] = '1';
            *shift = m - n + 1;
        }
        else
        {
            *shift = m - n;
        }

        s[n] = '\0'; /* truncate */
    }
}

void
arb_get_str_parts(int * negative, char **mid_digits, fmpz_t mid_exp,
                                  char **rad_digits, fmpz_t rad_exp,
                                  const arb_t x, slong n, int more)
{
    fmpz_t mid, rad, exp, err;
    slong good;
    flint_bitcnt_t shift;

    if (!arb_is_finite(x))
    {
        *negative = 0;

        fmpz_zero(mid_exp);
        *mid_digits = flint_malloc(4);
        if (arf_is_nan(arb_midref(x)))
            strcpy(*mid_digits, "nan");
        else
            strcpy(*mid_digits, "0");

        fmpz_zero(rad_exp);
        *rad_digits = flint_malloc(4);
        strcpy(*rad_digits, "inf");

        return;
    }

    fmpz_init(mid);
    fmpz_init(rad);
    fmpz_init(exp);
    fmpz_init(err);

    /* heuristic part */
    if (!more)
    {
        good = arb_rel_accuracy_bits(x) * 0.30102999566398119521 + 2;
        n = FLINT_MIN(n, good);
    }

    arb_get_fmpz_mid_rad_10exp(mid, rad, exp, x, FLINT_MAX(n, 1));
    *negative = arf_sgn(arb_midref(x)) < 0;
    fmpz_abs(mid, mid);

    *mid_digits = fmpz_get_str(NULL, 10, mid);
    *rad_digits = NULL;

    /* Truncate further so that 1 ulp error can be guaranteed (rigorous part)
       Note: mid cannot be zero here if n >= 1 and rad != 0. */
    if (n >= 1 && !(more || fmpz_is_zero(rad)))
    {
        slong lenmid, lenrad, rem;

        *rad_digits = fmpz_get_str(NULL, 10, rad);

        lenmid = strlen(*mid_digits);
        lenrad = strlen(*rad_digits);

        if (lenmid > lenrad)
        {
            /* we will truncate at n or n-1 */
            good = lenmid - lenrad;

            /* rounding to nearest can add at most 0.5 ulp */
            /* look at first omitted digit */
            rem = ((*mid_digits)[good]) - '0';
            if (rem < 5)
                rem = rem + 1;
            else
                rem = 10 - rem;

            /* and include the leading digit of the radius */
            rem = rem + ((*rad_digits)[0] - '0') + 1;

            /* if error is <= 1.0 ulp, we get to keep the extra digit */
            if (rem > 10)
                good -= 1;

            n = FLINT_MIN(n, good);
        }
        else
        {
            n = 0;
        }

        /* todo: avoid recomputing? */
        flint_free(*rad_digits);
    }

    /* no accurate digits -- output 0 +/- rad */
    if (n < 1)
    {
        fmpz_add(rad, rad, mid);
        fmpz_zero(mid);
        strcpy(*mid_digits, "0");  /* must have space already! */
    }
    else
    {
        _arb_digits_round_inplace(*mid_digits, &shift, err, n, ARF_RND_NEAR);
        fmpz_add_ui(mid_exp, exp, shift);
        fmpz_abs(err, err);
        fmpz_add(rad, rad, err);
    }

    /* write radius */
    if (fmpz_is_zero(rad))
    {
        *rad_digits = fmpz_get_str(NULL, 10, rad);
        fmpz_zero(rad_exp);
    }
    else
    {
        *rad_digits = fmpz_get_str(NULL, 10, rad);
        _arb_digits_round_inplace(*rad_digits, &shift, err, RADIUS_DIGITS, ARF_RND_UP);
        fmpz_add_ui(rad_exp, exp, shift);
    }

    fmpz_clear(mid);
    fmpz_clear(rad);
    fmpz_clear(exp);
    fmpz_clear(err);
}

char * arb_get_str(const arb_t x, slong n, ulong flags)
{
    char * res;
    char * mid_digits;
    char * rad_digits;
    int negative, more, skip_rad, skip_mid;
    fmpz_t mid_exp;
    fmpz_t rad_exp;
    slong condense;

    if (arb_is_zero(x))
    {
        res = flint_malloc(2);
        strcpy(res, "0");
        return res;
    }

    more = flags & ARB_STR_MORE;
    condense = flags / ARB_STR_CONDENSE;

    if (!arb_is_finite(x))
    {
        res = flint_malloc(10);

        if (arf_is_nan(arb_midref(x)))
            strcpy(res, "nan");
        else
            strcpy(res, "[+/- inf]");

        return res;
    }

    fmpz_init(mid_exp);
    fmpz_init(rad_exp);

    arb_get_str_parts(&negative, &mid_digits, mid_exp, &rad_digits, rad_exp, x, n, more);

    if ((flags & ARB_STR_NO_RADIUS) && mid_digits[0] == '0')
    {
        fmpz_add_ui(rad_exp, rad_exp, strlen(rad_digits));

        res = flint_malloc(fmpz_sizeinbase(rad_exp, 10) + 4);
        res[0] = '0';
        res[1] = 'e';
        if (fmpz_sgn(rad_exp) >= 0)
        {
            res[2] = '+';
            fmpz_get_str(res + 3, 10, rad_exp);
        }
        else
        {
            fmpz_get_str(res + 2, 10, rad_exp);
        }
    }
    else
    {
        skip_mid = mid_digits[0] == '0';
        skip_rad = (rad_digits[0] == '0') || ((flags & ARB_STR_NO_RADIUS) && !skip_mid);

        _arb_digits_as_float_str(&mid_digits, mid_exp, -4, FLINT_MAX(6, n - 1));
        _arb_digits_as_float_str(&rad_digits, rad_exp, -2, 2);

        if (skip_rad)
        {
            res = flint_malloc(strlen(mid_digits) + 2);

            if (negative)
                strcpy(res, "-");
            else
                strcpy(res, "");

            strcat(res, mid_digits);
        }
        else if (skip_mid)
        {
            res = flint_malloc(strlen(rad_digits) + 7);

            strcpy(res, "[+/- ");
            strcat(res, rad_digits);
            strcat(res, "]");
        }
        else
        {
            res = flint_malloc(strlen(mid_digits) + strlen(rad_digits) + 9);

            strcpy(res, "[");

            if (negative)
                strcat(res, "-");

            strcat(res, mid_digits);
            strcat(res, " +/- ");
            strcat(res, rad_digits);
            strcat(res, "]");
        }
    }

    if (condense)
        res = _arb_condense_digits(res, condense);

    flint_free(mid_digits);
    flint_free(rad_digits);

    fmpz_clear(mid_exp);
    fmpz_clear(rad_exp);

    return res;
}

