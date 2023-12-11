/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <ctype.h> /* isdigit */
#include <stdint.h> /* intmax_t */
#include <stdio.h>
#include <string.h> /* memcpy, memcmp and strchr */
#include <stdarg.h>
#include <wchar.h> /* wchar_t and wint_t */
#include "nmod_types.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpq_types.h"
#include "fmpq.h"
#include "arf_types.h"
#include "arb.h"
#include "acb.h"

/* Helper functions **********************************************************/

#define STRING_SIZE(x) (sizeof(x) - 1)
#define STRING_LENGTH(x) (STRING_SIZE(x) / sizeof(char))

/* NOTE: We have to be careful here to not alter `ip'. */
static inline void __arb_neg_readonly(arb_ptr op, arb_srcptr ip)
{
    arf_ptr op_arf;

    op->mid = ip->mid;
    op->rad = ip->rad;

    op_arf = &(op->mid);

    if (arf_is_special(op_arf))
    {
        if (arf_is_pos_inf(op_arf))
        {
            ARF_EXP(op_arf) = ARF_EXP_NEG_INF;
        }
        else if (arf_is_neg_inf(op_arf))
        {
            ARF_EXP(op_arf) = ARF_EXP_POS_INF;
        }
    }
    else
    {
        ARF_NEG(op_arf);
    }
}

static char * __mag_get_str(mag_srcptr ip, slong digits)
{
    arf_t ip_arf;
    char * str;

    arf_init(ip_arf);

    arf_set_mag(ip_arf, ip);
    str = arf_get_str(ip_arf, digits);

    arf_clear(ip_arf);

    return str;
}

static int __never_is(const void * ip)
{
    return 0;
}

/* is_pm1 */

static int __ulong_is_pm1(const void * ip)
{
    return *((const ulong *) ip) == 1;
}

static int __arb_is_pm1(const void * ip)
{
    return arb_is_one(ip) || arb_equal_si(ip, -1);
}

static int __acb_is_pm1(const void * ip)
{
    return acb_is_one(ip) || acb_equal_si(ip, -1);
}

/* is_zero */

static int __ulong_is_zero(const void * ip)
{
    return *((ulong *) ip) == 0;
}

/* is_neg */

static int __fmpz_is_neg(const void * ip)
{
    return fmpz_sgn(ip) < 0;
}

static int __fmpq_is_neg(const void * ip)
{
    return __fmpz_is_neg(fmpq_numref((const fmpq *) ip));
}

/* NOTE: The following function is not checking if a complex number is negative
 * mathematically. This simply a helper function to check whether it should
 * prepend with a minus or plus in the polynomial printing between terms. */
static int __acb_is_neg(const void * ip)
{
    /* We consider a complex number z negative in this sense if it is on one of
     * the following forms:
     *
     * - z is purely real and its real part is negative, or
     *
     * - z is purely imaginary and its imaginary part is negative.
     *
     * We could also consider it to be negative if both the real and imaginary
     * part are negative. However, I (Albin) do not think this is a good idea as
     * it makes too complex visually. */
    acb_srcptr zp = ip;
    arb_srcptr realzp = acb_realref(zp), imagzp = acb_imagref(zp);

    return (arb_is_zero(imagzp) && arb_is_negative(realzp))
        || (arb_is_zero(realzp) && arb_is_negative(imagzp));
}

/* Base types used in vectors, matrices and polynomials **********************/

typedef enum
{
    ulong_type = 0,
    fmpz_type,
    fmpq_type,
    arb_type,
    acb_type,
    slong_type,
    mag_type,
    arf_type
} flint_type_t;

static inline size_t flint_type_size_in_chars(flint_type_t type)
{
    if (type == ulong_type)
        return sizeof(ulong) / sizeof(char);
    else if (type == fmpz_type)
        return sizeof(fmpz) / sizeof(char);
    else if (type == fmpq_type)
        return sizeof(fmpq) / sizeof(char);
    else if (type == arb_type)
        return sizeof(arb_struct) / sizeof(char);
    else if (type == acb_type)
        return sizeof(acb_struct) / sizeof(char);
    else if (type == slong_type)
        return sizeof(slong) / sizeof(char);
    else if (type == mag_type)
        return sizeof(mag_struct) / sizeof(char);
    else /* if (type == arf_type) */
        return sizeof(arf_struct) / sizeof(char);
}

/* Declaring local printing functions ****************************************/

#define FLAG_NONE (0)
#define FLAG_NEG (1)
#define FLAG_PAREN (1 << 1) /* Indicates printing for polynomials */

#define FLAG_IS_NEG(flag) ((flag) & FLAG_NEG)
#define FLAG_IS_PAREN(flag) ((flag) & FLAG_PAREN)

static size_t __ulong_fprint(FILE *, const ulong *, int);
static size_t __slong_fprint(FILE *, const slong *, int);
static size_t __fmpz_fprint(FILE *, const fmpz *, int);
static size_t __fmpq_fprint(FILE *, const fmpq *, int);
static size_t __mag_fprint(FILE *, mag_srcptr);
static size_t __arf_fprint(FILE *, arf_srcptr);
static size_t __arb_fprint(FILE *, arb_srcptr, int);
static size_t __acb_fprint(FILE *, acb_srcptr, int);
static size_t __nmod_fprint(FILE *, nmod_t);
static size_t __fmpz_mod_ctx_fprint(FILE *, const fmpz_mod_ctx_struct *);
static size_t __flint_vec_fprint(FILE *, const void *, slong, flint_type_t);
static size_t __flint_mat_fprint(FILE *, const void *, flint_type_t);
static size_t __flint_poly_fprint(FILE *, const void *, flint_type_t);

/* flint_vfprintf and friends ************************************************/

/* TODO: Add options for compact/spacious printing. */

#define IS_FLINT_BASE_TYPE(ip, str) (memcmp(ip, str, sizeof(str) - sizeof(char)) == 0)
#define IS_FLINT_TYPE(ip, str) (memcmp(ip, str "}", sizeof(str)) == 0)

/* Reference used for checks: https://en.cppreference.com/w/c/io/fprintf */
#define IS_PRINTF_FLAG(chr) \
    (   (chr) == '-'        \
     || (chr) == '+'        \
     || (chr) == ' '        \
     || (chr) == '#'        \
     || (chr) == '0')

#define JUMP_FLAGS(str)     \
do                          \
{                           \
    (str)++;                \
} while (IS_PRINTF_FLAG(*(str)))

#define JUMP_MINIMAL_FIELD_WIDTH_WITH_POP(str, vlist) \
do                              \
{                               \
    if (isdigit(*(str)))        \
    {                           \
        (str)++;                \
        while (isdigit(*(str))) \
            (str)++;            \
    }                           \
    else if (*(str) == '*')     \
    {                           \
        va_arg((vlist), int);   \
        (str)++;                \
    }                           \
} while (0)

#define JUMP_PRECISION_WITH_POP(str, vlist) \
do                              \
{                               \
    if (*(str) == '.')          \
    {                           \
        (str)++;                \
        if (isdigit(*(str)))    \
        {                       \
            (str)++;            \
            while (isdigit(*(str))) \
                (str)++;        \
        }                       \
        else if (*(str) == '*') \
        {                       \
            va_arg((vlist), int); \
            (str)++;            \
        }                       \
    }                           \
} while (0)

#define _IS_PRINTF_INTEGERFMT(chr)  \
    (  (chr) == 'd' || (chr) == 'i' \
    || (chr) == 'o'                 \
    || (chr) == 'x' || (chr) == 'X' \
    || (chr) == 'u')

#define _IS_PRINTF_FLOATFMT(chr)    \
    (  (chr) == 'f' || (chr) == 'F' \
    || (chr) == 'e' || (chr) == 'E' \
    || (chr) == 'a' || (chr) == 'A' \
    || (chr) == 'g' || (chr) == 'G')

/* "Generic" ones */
#define IS_PRINTF_CHARFMT(str) \
    ((str)[0] == 'h' && (str)[1] == 'h' && _IS_PRINTF_INTEGERFMT((str)[2]))

#define IS_PRINTF_SHORTFMT(str) \
    ((str)[0] == 'h' && _IS_PRINTF_INTEGERFMT((str)[1]))

#define IS_PRINTF_INTFMT(str) \
    ((str)[0] == 'c' || _IS_PRINTF_INTEGERFMT((str)[0]))

#define IS_PRINTF_LONGFMT(str) \
    ((str)[0] == 'l' && _IS_PRINTF_INTEGERFMT((str)[1]))

#define IS_PRINTF_LONGLONGFMT(str) \
    ((str)[0] == 'l' && (str)[1] == 'l' && _IS_PRINTF_INTEGERFMT((str)[2]))

#define IS_PRINTF_INTMAXFMT(str) \
    ((str)[0] == 'j' &&  _IS_PRINTF_INTEGERFMT((str)[1]))

#define IS_PRINTF_SIZEFMT(str) \
    ((str)[0] == 'z' &&  _IS_PRINTF_INTEGERFMT((str)[1]))

#define IS_PRINTF_PTRDIFFFMT(str) \
    ((str)[0] == 't' &&  _IS_PRINTF_INTEGERFMT((str)[1]))

#define IS_FLINT_PRINTF_ULONGFMT(str) \
    ((str)[0] == 'w' &&  _IS_PRINTF_INTEGERFMT((str)[1]))

#define IS_PRINTF_DOUBLEFMT(str)    \
    (_IS_PRINTF_FLOATFMT((str)[0])  \
     || ((str)[0] == 'l' && _IS_PRINTF_FLOATFMT((str)[1])))

#define IS_PRINTF_LONGDOUBLEFMT(str) \
    ((str)[0] == 'L' && _IS_PRINTF_FLOATFMT((str)[1]))

/* "Special" ones */
#define IS_PRINTF_POINTERFMT(str) \
    ((str)[0] == 'p')

#define IS_PRINTF_CHARPFMT(str) \
    ((str)[0] == 's')

#define IS_PRINTF_WINTFMT(str) \
    ((str)[0] == 'l' && (str)[1] == 'c')

#define IS_PRINTF_WCHARPFMT(str) \
    ((str)[0] == 'l' && (str)[1] == 's')

int flint_vfprintf(FILE * fs, const char * ip, va_list vlist)
{
    size_t iplen;
    const char * ipcur;
    char * op, * opcur;
    size_t res;
    va_list vlist_cpy;
    int tmp;
    TMP_INIT;

    res = 0;
    iplen = strlen(ip);
    TMP_START;

#if defined(_LONG_LONG_LIMB)
    /*
       If mp_limb_t is long long, then

         `%(format args...)w' -> `%(format args...)ll'.

       The highest ratio between length of input and output string after this
       conversion is 4 / 3, which is obtained if `ip = "%wd...%wd"'.
     */
    op = TMP_ALLOC(sizeof(char) * (iplen + iplen / 3 + 1));
#else
    /* Same length due to length of "w" is equal to length of "l". */
    op = TMP_ALLOC(sizeof(char) * (iplen + 1));
#endif

    opcur = op;

    while (1)
    {
continue_while_from_flint_type:
        ipcur = ip;
        va_copy(vlist_cpy, vlist);

continue_while:
        ipcur = strchr(ipcur, '%');

        if (ipcur == NULL)
            break; /* Reached end of string */

        /* Check if "%%" */
        if (ipcur[1] == '%')
        {
            ipcur += 2;
            goto continue_while;
        }

        /* Check if "%{FLINT_TYPE}" */
        if (ipcur[1] == '{')
            goto print_flint_type;

        /* Check if "%(format args...)w" */
        JUMP_FLAGS(ipcur);
        /* NOTE: If the minimal field with and/or precision is specified, but
         * the format is invalid, the following pops will be invalid as well. */
        JUMP_MINIMAL_FIELD_WIDTH_WITH_POP(ipcur, vlist);
        JUMP_PRECISION_WITH_POP(ipcur, vlist);

        if (IS_FLINT_PRINTF_ULONGFMT(ipcur))
        {
            size_t cpsz;

            ipcur += 2; /* To include 'w' and following format specifier */
            cpsz = ipcur - ip;

            memcpy(opcur, ip, sizeof(char) * cpsz);
            ip = ipcur;

#if defined(_LONG_LONG_LIMB)
            opcur += cpsz + 1;
            opcur[-1] = opcur[-2];
            opcur[-2] = 'l';
            opcur[-3] = 'l';
#else
            opcur += cpsz;
            opcur[-2] = 'l';
#endif
            /* Pop entry from vlist (still present in vlist_cpy) */
            va_arg(vlist, ulong);
        }
        else if (IS_PRINTF_INTFMT(ipcur)
                || IS_PRINTF_CHARFMT(ipcur)
                || IS_PRINTF_SHORTFMT(ipcur))
            va_arg(vlist, int);
        else if (IS_PRINTF_DOUBLEFMT(ipcur))
            va_arg(vlist, double);
        else if (IS_PRINTF_SIZEFMT(ipcur))
            va_arg(vlist, size_t);
        else if (IS_PRINTF_CHARPFMT(ipcur))
            va_arg(vlist, char *);
        else if (IS_PRINTF_LONGFMT(ipcur))
            va_arg(vlist, long int);
        else if (IS_PRINTF_LONGLONGFMT(ipcur))
            va_arg(vlist, long long int);
        else if (IS_PRINTF_INTMAXFMT(ipcur))
            va_arg(vlist, intmax_t);
        else if (IS_PRINTF_PTRDIFFFMT(ipcur))
            va_arg(vlist, ptrdiff_t);
        else if (IS_PRINTF_LONGDOUBLEFMT(ipcur))
            va_arg(vlist, long double);
        else if (IS_PRINTF_POINTERFMT(ipcur))
            va_arg(vlist, void *);
        else if (IS_PRINTF_WCHARPFMT(ipcur))
            va_arg(vlist, wchar_t *);
        else if (IS_PRINTF_WINTFMT(ipcur))
        {
            /* NOTE: MinGW defines wint_t as unsigned short int */
            if (sizeof(wint_t) >= sizeof(int))
                va_arg(vlist, wint_t);
            else
                va_arg(vlist, int);
        }

        goto continue_while;
    }

    {
        size_t cpsz = strlen(ip);
        memcpy(opcur, ip, sizeof(char) * cpsz);
        opcur[cpsz] = '\0';
    }

    tmp = vfprintf(fs, op, vlist_cpy);
    if (tmp < 0)
        res = tmp;
    else
        res += tmp;

end:
    TMP_END;
    va_end(vlist_cpy);

    return (int) res;

print_flint_type:
    /* Print now to be able to print FLINT types */
    {
        size_t cpsz = ipcur - ip;
        memcpy(opcur, ip, sizeof(char) * cpsz);
        opcur[cpsz] = '\0';
    }
    tmp = vfprintf(fs, op, vlist_cpy);
    if (tmp < 0)
    {
        res = tmp;
        goto end;
    }
    else
        res += tmp;

    /* vlist_cpy is now invalid */
    va_end(vlist_cpy);

    opcur = op;
    ip = ipcur + 2; /* Now `ip' points to "FLINT_TYPE..." */

    if (IS_FLINT_BASE_TYPE(ip, "ulong"))
    {
        if (IS_FLINT_TYPE(ip, "ulong"))
        {
            res += fprintf(fs, WORD_FMT "u", va_arg(vlist, ulong));
            ip += STRING_LENGTH("ulong}");
        }
        else if (IS_FLINT_TYPE(ip, "ulong*"))
        {
            const ulong * vec = va_arg(vlist, const ulong *);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_fprint(fs, vec, len, ulong_type);
            ip += STRING_LENGTH("ulong*}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "slong"))
    {
        if (IS_FLINT_TYPE(ip, "slong"))
        {
            res += fprintf(fs, WORD_FMT "d", va_arg(vlist, slong));
            ip += STRING_LENGTH("slong}");
        }
        else if (IS_FLINT_TYPE(ip, "slong*"))
        {
            const slong * vec = va_arg(vlist, const slong *);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_fprint(fs, vec, len, slong_type);
            ip += STRING_LENGTH("slong*}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "nmod"))
    {
        if (IS_FLINT_TYPE(ip, "nmod"))
        {
            res += __nmod_fprint(fs, va_arg(vlist, nmod_t));
            ip += STRING_LENGTH("nmod}");
        }
        else if (IS_FLINT_TYPE(ip, "nmod_mat"))
        {
            res += __flint_mat_fprint(fs, va_arg(vlist, const nmod_mat_struct *), ulong_type);
            ip += STRING_LENGTH("nmod_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "nmod_poly"))
        {
            res += __flint_poly_fprint(fs, va_arg(vlist, const nmod_poly_struct *), ulong_type);
            ip += STRING_LENGTH("nmod_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "fmpz")) /* fmpz or fmpz_mod base type */
    {
        if (IS_FLINT_TYPE(ip, "fmpz"))
        {
            res += __fmpz_fprint(fs, va_arg(vlist, const fmpz *), FLAG_NONE);
            ip += STRING_LENGTH("fmpz}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz*"))
        {
            const fmpz * vec = va_arg(vlist, const fmpz *);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_fprint(fs, vec, len, fmpz_type);
            ip += STRING_LENGTH("fmpz*}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_mat"))
        {
            res += __flint_mat_fprint(fs, va_arg(vlist, const fmpz_mat_struct *), fmpz_type);
            ip += STRING_LENGTH("fmpz_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_poly"))
        {
            res += __flint_poly_fprint(fs, va_arg(vlist, const fmpz_poly_struct *), fmpz_type);
            ip += STRING_LENGTH("fmpz_poly}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_mod_ctx"))
        {
            res += __fmpz_mod_ctx_fprint(fs, va_arg(vlist, const fmpz_mod_ctx_struct *));
            ip += STRING_LENGTH("fmpz_mod_ctx}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_mod_mat"))
        {
            /* Print as if fmpz_mat */
            res += __flint_mat_fprint(fs, va_arg(vlist, const fmpz_mod_mat_struct *)->mat, fmpz_type);
            ip += STRING_LENGTH("fmpz_mod_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_mod_poly"))
        {
            /* Print as if fmpz_poly */
            res += __flint_poly_fprint(fs, va_arg(vlist, const fmpz_mod_poly_struct *), fmpz_type);
            ip += STRING_LENGTH("fmpz_mod_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "fmpq"))
    {
        if (IS_FLINT_TYPE(ip, "fmpq"))
        {
            res += __fmpq_fprint(fs, va_arg(vlist, const fmpq *), FLAG_NONE);
            ip += STRING_LENGTH("fmpq}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpq*"))
        {
            const fmpq * vec = va_arg(vlist, const fmpq *);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_fprint(fs, vec, len, fmpq_type);
            ip += STRING_LENGTH("fmpq*}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpq_mat"))
        {
            res += __flint_mat_fprint(fs, va_arg(vlist, const fmpq_mat_struct *), fmpq_type);
            ip += STRING_LENGTH("fmpq_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpq_poly"))
        {
            res += __flint_poly_fprint(fs, va_arg(vlist, const fmpq_poly_struct *), fmpq_type);
            ip += STRING_LENGTH("fmpq_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "arf"))
    {
        if (IS_FLINT_TYPE(ip, "arf"))
        {
            res += __arf_fprint(fs, va_arg(vlist, arf_srcptr));
            ip += STRING_LENGTH("arf}");
        }
        else if (IS_FLINT_TYPE(ip, "arf*"))
        {
            arf_srcptr vec = va_arg(vlist, arf_srcptr);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_fprint(fs, vec, len, arf_type);
            ip += STRING_LENGTH("arf*}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "mag"))
    {
        if (IS_FLINT_TYPE(ip, "mag"))
        {
            res += __mag_fprint(fs, va_arg(vlist, mag_srcptr));
            ip += STRING_LENGTH("mag}");
        }
        else if (IS_FLINT_TYPE(ip, "mag*"))
        {
            mag_srcptr vec = va_arg(vlist, mag_srcptr);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_fprint(fs, vec, len, mag_type);
            ip += STRING_LENGTH("mag*}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "arb"))
    {
        if (IS_FLINT_TYPE(ip, "arb"))
        {
            res += __arb_fprint(fs, va_arg(vlist, arb_srcptr), FLAG_NONE);
            ip += STRING_LENGTH("arb}");
        }
        else if (IS_FLINT_TYPE(ip, "arb*"))
        {
            arb_srcptr vec = va_arg(vlist, arb_srcptr);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_fprint(fs, vec, len, arb_type);
            ip += STRING_LENGTH("arb*}");
        }
        else if (IS_FLINT_TYPE(ip, "arb_mat"))
        {
            res += __flint_mat_fprint(fs, va_arg(vlist, const arb_mat_struct *), arb_type);
            ip += STRING_LENGTH("arb_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "arb_poly"))
        {
            res += __flint_poly_fprint(fs, va_arg(vlist, const arb_poly_struct *), arb_type);
            ip += STRING_LENGTH("arb_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "acb"))
    {
        if (IS_FLINT_TYPE(ip, "acb"))
        {
            res += __acb_fprint(fs, va_arg(vlist, acb_srcptr), FLAG_NONE);
            ip += STRING_LENGTH("acb}");
        }
        else if (IS_FLINT_TYPE(ip, "acb*"))
        {
            acb_srcptr vec = va_arg(vlist, acb_srcptr);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_fprint(fs, vec, len, acb_type);
            ip += STRING_LENGTH("acb*}");
        }
        else if (IS_FLINT_TYPE(ip, "acb_mat"))
        {
            res += __flint_mat_fprint(fs, va_arg(vlist, const acb_mat_struct *), acb_type);
            ip += STRING_LENGTH("acb_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "acb_poly"))
        {
            res += __flint_poly_fprint(fs, va_arg(vlist, const acb_poly_struct *), acb_type);
            ip += STRING_LENGTH("acb_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else
    {
printpercentcurlybracket:
        /* Invalid use of "%{FLINT_TYPE}". As we are currently pointed to
         * "FLINT_TYPE}", we let fprintf take care of printing "%{". */
#ifdef __GNUC__
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wformat"
#endif
        tmp = fprintf(fs, "%{");
#ifdef __GNUC__
# pragma GCC diagnostic pop
#endif
        if (tmp < 0)
        {
            res = tmp;
            goto end;
        }
        else
            res += tmp;
    }
    goto continue_while_from_flint_type;
}

int flint_printf(const char * str, ...)
{
   va_list vlist;
   int ret;

   va_start(vlist, str);
   ret = flint_vfprintf(stdout, str, vlist);
   va_end(vlist);

   return ret;
}

int flint_fprintf(FILE * fs, const char * str, ...)
{
   va_list vlist;
   int ret;

   va_start(vlist, str);
   ret = flint_vfprintf(fs, str, vlist);
   va_end(vlist);

   return ret;
}

int flint_vprintf(const char * str, va_list vlist)
{
    return flint_vfprintf(stdout, str, vlist);
}

/* Type specific printing ****************************************************/

/* TODO: Move these to their respective module? */

/* TODO: Add precision input to Arb type printing functions */

/* TODO: Add option to print in different basis? */

static size_t __ulong_fprint(FILE * fs, const ulong * ip, int flag)
{
    return fprintf(fs, WORD_FMT "u", *ip);
}

static size_t __slong_fprint(FILE * fs, const slong * ip, int flag)
{
    return fprintf(fs, WORD_FMT "d", *ip);
}

#define BASE 10
static size_t __fmpz_fprint(FILE * fs, const fmpz * ip, int flag)
{
    size_t res = 0;
    char * str;
    size_t skipminus = __fmpz_is_neg(ip) ? FLAG_IS_NEG(flag) : 0;

    str = fmpz_get_str(NULL, BASE, ip);
    res += fwrite(str + skipminus, sizeof(char), strlen(str + skipminus), fs);

    flint_free(str);

    return res;
}
#undef BASE

static size_t __fmpq_fprint(FILE * fs, const fmpq * ip, int flag)
{
    size_t res = 0;

    /* NOTE: We do not care about parentheses here, and __fmpz_fprint does not
     * care either. */
    res += __fmpz_fprint(fs, fmpq_numref(ip), flag);
    if (!fmpz_is_one(fmpq_denref(ip)))
    {
        res += fwrite(" / ", sizeof(char), STRING_LENGTH(" / "), fs);
        res += __fmpz_fprint(fs, fmpq_denref(ip), FLAG_NONE);
    }

    return res;
}

#define DIGITS 6
static size_t __mag_fprint(FILE * fs, mag_srcptr ip)
{
    size_t res;
    char * str;

#if DIGITS != 6
# error Change fwrite below
#endif
    if (mag_is_zero(ip))
        return fwrite("0.00000", sizeof(char), STRING_SIZE("0.00000"), fs);

    str = __mag_get_str(ip, DIGITS);
    res = fwrite(str, sizeof(char), strlen(str), fs);

    flint_free(str);

    return res;
}

static size_t __arf_fprint(FILE * fs, arf_srcptr ip)
{
    size_t res;
    char * str;

#if DIGITS != 6
# error Change fwrite below
#endif
    if (arf_is_zero(ip))
        return fwrite("0.00000", sizeof(char), STRING_SIZE("0.00000"), fs);

    str = arf_get_str(ip, DIGITS);
    res = fwrite(str, sizeof(char), strlen(str), fs);

    flint_free(str);

    return res;
}

#define MAX_INT_SIZE 64
/* NOTE: If arb is an integer, we print it as one. */
static size_t __arb_fprint(FILE * fs, arb_srcptr ip, int flag)
{
    size_t res;

    if (arb_is_zero(ip))
    {
        return (fputc('0', fs) != EOF);
    }
    else if (arb_is_int(ip) && ARF_EXP(arb_midref(ip)) <= MAX_INT_SIZE)
    {
        /* NOTE: Only print as integer if ip < 2^64. */
        fmpz_t fip;

        fmpz_init(fip);

        /* NOTE: The conversion should be exact, so we do not care about its
         * return value. */
        arf_get_fmpz(fip, arb_midref(ip), ARF_RND_DOWN);
        /* NOTE: We do not care about parentheses here, and __fmpz_fprint does not
         * care either. */
        res = __fmpz_fprint(fs, fip, flag);

        fmpz_clear(fip);
    }
    else
    {
        char * str;
        arb_struct ip2;

        if (FLAG_IS_NEG(flag))
            __arb_neg_readonly(&ip2, ip);
        else
            ip2 = *ip;

        str = arb_get_str(&ip2, DIGITS, 0);
        res = fwrite(str, sizeof(char), strlen(str), fs);

        flint_free(str);
    }

    return res;
}
#undef MAX_INT_SIZE
#undef DIGITS

static size_t __acb_fprint(FILE * fs, acb_srcptr ip, int flag)
{
    size_t res = 0;
    int realiszero, imagiszero;

    realiszero = arb_is_zero(acb_realref(ip));
    imagiszero = arb_is_zero(acb_imagref(ip));

    if (realiszero && imagiszero)
        return (fputc('0', fs) != EOF);

    /* Only print parentheses if both real and imaginary part is non-zero. */
    if (FLAG_IS_PAREN(flag) && !realiszero && !imagiszero)
        res += (fputc('(', fs) != EOF);

    /* Print real part if non-zero */
    if (!realiszero)
        res += __arb_fprint(fs, acb_realref(ip), flag);

    /* Print imaginary part if non-zero */
    if (!imagiszero)
    {
        int imagisneg = arb_is_negative(acb_imagref(ip));

        if (!realiszero)
            res += fwrite(FLAG_IS_NEG(flag) ^ imagisneg ? " - " : " + ", sizeof(char), STRING_SIZE(" - "), fs);

        /* If imaginary part of ip is \pm 1, then we only print \pm i. */
        if (!__arb_is_pm1(acb_imagref(ip)))
        {
            res += __arb_fprint(fs, acb_imagref(ip), (FLAG_IS_NEG(flag) ^ imagisneg) & (~realiszero));
            res += fwrite(" * ", sizeof(char), STRING_LENGTH(" * "), fs);
        }

        res += (fputc('i', fs) != EOF);
    }

    if (FLAG_IS_PAREN(flag) && !realiszero && !imagiszero)
        res += (fputc(')', fs) != EOF);

    return res;
}

static size_t __nmod_fprint(FILE * fs, nmod_t ip)
{
    return fprintf(fs, "mod " WORD_FMT "u", ip.n);
}

static size_t __fmpz_mod_ctx_fprint(FILE * fs, const fmpz_mod_ctx_struct * ip)
{
    size_t res = 0;

    res += fwrite("mod ", sizeof(char), STRING_LENGTH("mod "), fs);
    res += fmpz_fprint(fs, fmpz_mod_ctx_modulus(ip));

    return res;
}

/* Generic printing **********************************************************/

/* TODO: Add non-compact printing? For this it would be nice with option for
 * variable indentation and aligned columns. */

/* TODO: The square brackets used as delimiters for vectors and matrices are
 * also used for printing by arb_get_str. Think if we want to change something
 * here. */

/* TODO: Add option to specify generator for polynomials */

/* TODO: If it is possible to obtain terminal width and height, it would be nice
 * to be able to only print the beginning and ending element of vectors,
 * matrices and polynomials (such as what Julia does). This should be the
 * default. */

typedef size_t (* print_func_t)(FILE *, const void *, int);

static print_func_t print_functions[] =
{
    (print_func_t) __ulong_fprint,
    (print_func_t) __fmpz_fprint,
    (print_func_t) __fmpq_fprint,
    (print_func_t) __arb_fprint,
    (print_func_t) __acb_fprint,
    (print_func_t) __slong_fprint, /* NOTE: These print functions are only  */
    (print_func_t) __mag_fprint,   /* used for vectors, not other composite */
    (print_func_t) __arf_fprint    /* types. */
};

static size_t __flint_vec_fprint(FILE * fs, const void * ip, slong len, flint_type_t type)
{
    size_t res = 0;
    slong ix;
    const char * vec = ip;
    size_t type_size = flint_type_size_in_chars(type);
    print_func_t print = print_functions[type];

    res += (fputc('[', fs) != EOF);

    if (len > 0)
        res += print(fs, vec, FLAG_NONE);

    for (ix = 1; ix < len; ix++)
    {
        res += fwrite(", ", sizeof(char), STRING_LENGTH(", "), fs);
        res += print(fs, vec + type_size * ix, FLAG_NONE);
    }

    res += (fputc(']', fs) != EOF);

    return res;
}

/* NOTE: This function relies on the fact that the layout of
 * [fmpz/fmpq/arb/acb]_mat_struct are all on the form (pointer, slong, slong,
 * pointer). */
static size_t __flint_mat_fprint(FILE * fs, const void * ip, flint_type_t type)
{
    size_t res = 0;
    slong ix;
    slong nr, nc;
    const void ** rows;

    rows = (const void **) ((const fmpz_mat_struct *) ip)->rows;
    nr = ((const fmpz_mat_struct *) ip)->r;
    nc = ((const fmpz_mat_struct *) ip)->c;

    if (nr == 0 || nc == 0)
        return fprintf(fs, WORD_FMT "d by " WORD_FMT "d empty matrix", nr, nc);

    res += (fputc('[', fs) != EOF);
    res += __flint_vec_fprint(fs, rows[0], nc, type);

    for (ix = 1; ix < nr; ix++)
    {
        res += fwrite(", ", sizeof(char), STRING_LENGTH(", "), fs);
        res += __flint_vec_fprint(fs, rows[ix], nc, type);
    }

    res += (fputc(']', fs) != EOF);

    return res;
}

typedef int (* is_func_t)(const void *);

static is_func_t is_pm1_functions[] =
{
    __ulong_is_pm1,
    (is_func_t) fmpz_is_pm1,
    __never_is, /* is not used for fmpq */
    __arb_is_pm1,
    __acb_is_pm1
};

static is_func_t is_zero_functions[] =
{
    __ulong_is_zero,
    (is_func_t) fmpz_is_zero,
    __never_is, /* is not used for fmpq */
    (is_func_t) arb_is_zero,
    (is_func_t) acb_is_zero
};

static is_func_t is_neg_functions[] =
{
    __never_is, /* is not used for ulong */
    __fmpz_is_neg,
    __never_is, /* is not used for fmpq */
    (is_func_t) arb_is_negative,
    __acb_is_neg
};

/* NOTE: This function relies on the fact that the layout of
 * [fmpz/arb/acb]_poly_struct are all on the form (pointer, slong, slong) and
 * fmpq_poly_struct on the form (pointer, slong, slong, pointer). */
static size_t __flint_poly_fprint(FILE * fs, const void * ip, flint_type_t type)
{
    size_t res = 0;
    slong ix;
    slong len;

    len = ((const fmpz_poly_struct *) ip)->length;

    if (len == 0)
        return fputc('0', fs) != EOF;

    if (type != fmpq_type)
    {
        size_t type_size = flint_type_size_in_chars(type);
        const char * coeffs = (const char *) ((const fmpz_poly_struct *) ip)->coeffs;
        print_func_t print = print_functions[type];
        is_func_t is_pm1 = is_pm1_functions[type];
        is_func_t is_zero = is_zero_functions[type];
        is_func_t is_neg = is_neg_functions[type];

        if (len == 1)
            return print(fs, coeffs, FLAG_NONE);

        /* Leading coefficient cannot be zero */
        if (!is_pm1(coeffs + type_size * (len - 1)))
        {
            res += print(fs, coeffs + type_size * (len - 1), FLAG_PAREN);
            res += fwrite(" * ", sizeof(char), STRING_LENGTH(" * "), fs);
        }
        else if (is_neg(coeffs + type_size * (len - 1)))
            res += (fputc('-', fs) != EOF);
        if (len != 2)
            res += fprintf(fs, "x^" WORD_FMT "d", len - 1);
        else
            res += (fputc('x', fs) != EOF);

        for (ix = len - 2; ix > 0; ix--)
        {
            if (!is_zero(coeffs + type_size * ix))
            {
                res += fwrite(is_neg(coeffs + type_size * ix) ? " - " : " + ", sizeof(char), STRING_LENGTH(" - "), fs);
                if (!is_pm1(coeffs + type_size * ix))
                {
                    res += print(fs, coeffs + type_size * ix, is_neg(coeffs + type_size * ix) | FLAG_PAREN);
                    res += fwrite(" * ", sizeof(char), STRING_LENGTH(" * "), fs);
                }
                if (ix != 1)
                    res += fprintf(fs, "x^" WORD_FMT "d", ix);
                else
                    res += (fputc('x', fs) != EOF);
            }
        }

        if (!is_zero(coeffs + 0))
        {
            res += fwrite(is_neg(coeffs + 0) ? " - " : " + ", sizeof(char), STRING_LENGTH(" - "), fs);
            res += print(fs, coeffs + 0, is_neg(coeffs + 0));
        }
    }
    else
    {
        /* fmpq_poly is special as it is an fmpz_poly with a denomitator
         * strapped onto it */
        const fmpz * coeffs = ((const fmpq_poly_struct *) ip)->coeffs;
        const fmpz * den = ((const fmpq_poly_struct *) ip)->den;
        fmpq_t canonical;

        fmpq_init(canonical);
        fmpq_set_fmpz_frac(canonical, coeffs + len - 1, den);

        if (len == 1)
        {
            res += __fmpq_fprint(fs, canonical, FLAG_NONE);
            fmpq_clear(canonical);
            return res;
        }

        /* Leading coefficient cannot be zero */
        if (!fmpq_is_pm1(canonical))
        {
            res += __fmpq_fprint(fs, canonical, FLAG_NONE);
            res += fwrite(" * ", sizeof(char), STRING_LENGTH(" * "), fs);
        }
        else if (__fmpq_is_neg(canonical))
            res += (fputc('-', fs) != EOF);
        if (len != 2)
            res += fprintf(fs, "x^" WORD_FMT "d", len - 1);
        else
            res += (fputc('x', fs) != EOF);

        for (ix = len - 2; ix > 0; ix--)
        {
            if (!fmpz_is_zero(coeffs + ix))
            {
                fmpq_set_fmpz_frac(canonical, coeffs + ix, den);
                res += fwrite(__fmpq_is_neg(canonical) ? " - " : " + ", sizeof(char), STRING_LENGTH(" - "), fs);
                if (!fmpq_is_pm1(canonical))
                {
                    res += __fmpq_fprint(fs, canonical, __fmpq_is_neg(canonical));
                    res += fwrite(" * ", sizeof(char), STRING_LENGTH(" * "), fs);
                }
                if (ix != 1)
                    res += fprintf(fs, "x^" WORD_FMT "d", ix);
                else
                    res += (fputc('x', fs) != EOF);
            }
        }

        fmpq_set_fmpz_frac(canonical, coeffs + 0, den);
        if (!fmpq_is_zero(canonical))
        {
            res += fwrite(__fmpq_is_neg(canonical) ? " - " : " + ", sizeof(char), STRING_LENGTH(" - "), fs);
            res += __fmpq_fprint(fs, canonical, __fmpq_is_neg(canonical));
        }

        fmpq_clear(canonical);
    }

    return res;
}
