/*
    Copyright (C) 2023 Albin Ahlbäck
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Shared implementation for flint_vfprintf-style functions.
    This file is intentionally included multiple times with different sinks.
*/

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

static int __never_is(const void * FLINT_UNUSED(ip))
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

static size_t __ulong_print(FLINT_VPRINTF_OUT_T *, const ulong *, int);
static size_t __slong_print(FLINT_VPRINTF_OUT_T *, const slong *, int);
static size_t __fmpz_print(FLINT_VPRINTF_OUT_T *, const fmpz *, int);
static size_t __fmpq_print(FLINT_VPRINTF_OUT_T *, const fmpq *, int);
static size_t __mpz_print(FLINT_VPRINTF_OUT_T *, mpz_srcptr);
static size_t __mpq_print(FLINT_VPRINTF_OUT_T *, mpq_srcptr);
static size_t __mag_print(FLINT_VPRINTF_OUT_T *, mag_srcptr, int);
static size_t __arf_print(FLINT_VPRINTF_OUT_T *, arf_srcptr, int);
static size_t __arb_print(FLINT_VPRINTF_OUT_T *, arb_srcptr, int);
static size_t __acb_print(FLINT_VPRINTF_OUT_T *, acb_srcptr, int);
static size_t __nmod_print(FLINT_VPRINTF_OUT_T *, nmod_t);
static size_t __fmpz_mod_ctx_print(FLINT_VPRINTF_OUT_T *, const fmpz_mod_ctx_struct *);
static size_t __flint_vec_print(FLINT_VPRINTF_OUT_T *, const void *, slong, flint_type_t);
static size_t __flint_mat_print(FLINT_VPRINTF_OUT_T *, const void *, flint_type_t);
static size_t __flint_poly_print(FLINT_VPRINTF_OUT_T *, const void *, flint_type_t);

/* flint_vfprintf and friends ************************************************/

/* TODO: Add options for compact/spacious printing. */

#define IS_FLINT_BASE_TYPE(ip, str) (strncmp(ip, str, STRING_LENGTH(str)) == 0)
#define IS_FLINT_TYPE(ip, str) (strncmp(ip, str "}", STRING_LENGTH(str) + 1) == 0)

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

int FLINT_VPRINTF_FUNCTION(FLINT_VPRINTF_FUNCTION_ARGS, const char * ip, va_list vlist)
{
    size_t iplen;
    const char * ipcur;
    char * op, * opcur;
    size_t res;
    va_list vlist_cpy;
    int tmp;
    FLINT_VPRINTF_OUT_T out_state;
    FLINT_VPRINTF_OUT_T * out = &out_state;
    TMP_INIT;

    FLINT_VPRINTF_INIT(out);

    res = 0;
    iplen = strlen(ip);
    TMP_START;

#if FLINT_LONG_LONG
    /*
       If ulong is long long, then

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

#if FLINT_LONG_LONG
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

    tmp = FLINT_VPRINTF_VPRINTF(out, op, vlist_cpy);
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
    tmp = FLINT_VPRINTF_VPRINTF(out, op, vlist_cpy);
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
            res += FLINT_VPRINTF_PRINTF(out, WORD_FMT "u", va_arg(vlist, ulong));
            ip += STRING_LENGTH("ulong}");
        }
        else if (IS_FLINT_TYPE(ip, "ulong*"))
        {
            const ulong * vec = va_arg(vlist, const ulong *);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_print(out, vec, len, ulong_type);
            ip += STRING_LENGTH("ulong*}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "slong"))
    {
        if (IS_FLINT_TYPE(ip, "slong"))
        {
            res += FLINT_VPRINTF_PRINTF(out, WORD_FMT "d", va_arg(vlist, slong));
            ip += STRING_LENGTH("slong}");
        }
        else if (IS_FLINT_TYPE(ip, "slong*"))
        {
            const slong * vec = va_arg(vlist, const slong *);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_print(out, vec, len, slong_type);
            ip += STRING_LENGTH("slong*}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "nmod"))
    {
        if (IS_FLINT_TYPE(ip, "nmod"))
        {
            res += __nmod_print(out, va_arg(vlist, nmod_t));
            ip += STRING_LENGTH("nmod}");
        }
        else if (IS_FLINT_TYPE(ip, "nmod_mat"))
        {
            res += __flint_mat_print(out, va_arg(vlist, const nmod_mat_struct *), ulong_type);
            ip += STRING_LENGTH("nmod_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "nmod_poly"))
        {
            res += __flint_poly_print(out, va_arg(vlist, const nmod_poly_struct *), ulong_type);
            ip += STRING_LENGTH("nmod_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "fmpz")) /* fmpz or fmpz_mod base type */
    {
        if (IS_FLINT_TYPE(ip, "fmpz"))
        {
            res += __fmpz_print(out, va_arg(vlist, const fmpz *), FLAG_NONE);
            ip += STRING_LENGTH("fmpz}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz*"))
        {
            const fmpz * vec = va_arg(vlist, const fmpz *);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_print(out, vec, len, fmpz_type);
            ip += STRING_LENGTH("fmpz*}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_mat"))
        {
            res += __flint_mat_print(out, va_arg(vlist, const fmpz_mat_struct *), fmpz_type);
            ip += STRING_LENGTH("fmpz_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_poly"))
        {
            res += __flint_poly_print(out, va_arg(vlist, const fmpz_poly_struct *), fmpz_type);
            ip += STRING_LENGTH("fmpz_poly}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_mod_ctx"))
        {
            res += __fmpz_mod_ctx_print(out, va_arg(vlist, const fmpz_mod_ctx_struct *));
            ip += STRING_LENGTH("fmpz_mod_ctx}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_mod_mat"))
        {
            /* Print as if fmpz_mat */
            res += __flint_mat_print(out, va_arg(vlist, const fmpz_mod_mat_struct *), fmpz_type);
            ip += STRING_LENGTH("fmpz_mod_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpz_mod_poly"))
        {
            /* Print as if fmpz_poly */
            res += __flint_poly_print(out, va_arg(vlist, const fmpz_mod_poly_struct *), fmpz_type);
            ip += STRING_LENGTH("fmpz_mod_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "fmpq"))
    {
        if (IS_FLINT_TYPE(ip, "fmpq"))
        {
            res += __fmpq_print(out, va_arg(vlist, const fmpq *), FLAG_NONE);
            ip += STRING_LENGTH("fmpq}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpq*"))
        {
            const fmpq * vec = va_arg(vlist, const fmpq *);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_print(out, vec, len, fmpq_type);
            ip += STRING_LENGTH("fmpq*}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpq_mat"))
        {
            res += __flint_mat_print(out, va_arg(vlist, const fmpq_mat_struct *), fmpq_type);
            ip += STRING_LENGTH("fmpq_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "fmpq_poly"))
        {
            res += __flint_poly_print(out, va_arg(vlist, const fmpq_poly_struct *), fmpq_type);
            ip += STRING_LENGTH("fmpq_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "arf"))
    {
        if (IS_FLINT_TYPE(ip, "arf"))
        {
            res += __arf_print(out, va_arg(vlist, arf_srcptr), 0);
            ip += STRING_LENGTH("arf}");
        }
        else if (IS_FLINT_TYPE(ip, "arf*"))
        {
            arf_srcptr vec = va_arg(vlist, arf_srcptr);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_print(out, vec, len, arf_type);
            ip += STRING_LENGTH("arf*}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "mag"))
    {
        if (IS_FLINT_TYPE(ip, "mag"))
        {
            res += __mag_print(out, va_arg(vlist, mag_srcptr), 0);
            ip += STRING_LENGTH("mag}");
        }
        else if (IS_FLINT_TYPE(ip, "mag*"))
        {
            mag_srcptr vec = va_arg(vlist, mag_srcptr);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_print(out, vec, len, mag_type);
            ip += STRING_LENGTH("mag*}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "arb"))
    {
        if (IS_FLINT_TYPE(ip, "arb"))
        {
            res += __arb_print(out, va_arg(vlist, arb_srcptr), FLAG_NONE);
            ip += STRING_LENGTH("arb}");
        }
        else if (IS_FLINT_TYPE(ip, "arb*"))
        {
            arb_srcptr vec = va_arg(vlist, arb_srcptr);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_print(out, vec, len, arb_type);
            ip += STRING_LENGTH("arb*}");
        }
        else if (IS_FLINT_TYPE(ip, "arb_mat"))
        {
            res += __flint_mat_print(out, va_arg(vlist, const arb_mat_struct *), arb_type);
            ip += STRING_LENGTH("arb_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "arb_poly"))
        {
            res += __flint_poly_print(out, va_arg(vlist, const arb_poly_struct *), arb_type);
            ip += STRING_LENGTH("arb_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_BASE_TYPE(ip, "acb"))
    {
        if (IS_FLINT_TYPE(ip, "acb"))
        {
            res += __acb_print(out, va_arg(vlist, acb_srcptr), FLAG_NONE);
            ip += STRING_LENGTH("acb}");
        }
        else if (IS_FLINT_TYPE(ip, "acb*"))
        {
            acb_srcptr vec = va_arg(vlist, acb_srcptr);
            slong len = va_arg(vlist, slong);
            res += __flint_vec_print(out, vec, len, acb_type);
            ip += STRING_LENGTH("acb*}");
        }
        else if (IS_FLINT_TYPE(ip, "acb_mat"))
        {
            res += __flint_mat_print(out, va_arg(vlist, const acb_mat_struct *), acb_type);
            ip += STRING_LENGTH("acb_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "acb_poly"))
        {
            res += __flint_poly_print(out, va_arg(vlist, const acb_poly_struct *), acb_type);
            ip += STRING_LENGTH("acb_poly}");
        }
        else
            goto printpercentcurlybracket;
    }
    else if (IS_FLINT_TYPE(ip, "mpz"))
    {
        res += __mpz_print(out, va_arg(vlist, mpz_srcptr));
        ip += STRING_LENGTH("mpz}");
    }
    else if (IS_FLINT_TYPE(ip, "mpq"))
    {
        res += __mpq_print(out, va_arg(vlist, mpq_srcptr));
        ip += STRING_LENGTH("mpq}");
    }
    else if (IS_FLINT_TYPE(ip, "truth"))
    {
        truth_t t = va_arg(vlist, truth_t);
        if (t == T_TRUE) res += FLINT_VPRINTF_PRINTF(out, "T_TRUE");
        else if (t == T_FALSE) res += FLINT_VPRINTF_PRINTF(out, "T_FALSE");
        else res += FLINT_VPRINTF_PRINTF(out, "T_UNKNOWN");
        ip += STRING_LENGTH("truth}");
    }
    else if (IS_FLINT_BASE_TYPE(ip, "gr"))
    {
        gr_stream_t gr_out;
        FLINT_VPRINTF_GR_STREAM_INIT(gr_out, out);

        if (IS_FLINT_TYPE(ip, "gr"))
        {
            gr_srcptr elem = va_arg(vlist, gr_srcptr);
            gr_ctx_struct * ctx = va_arg(vlist, gr_ctx_struct *);
            GR_MUST_SUCCEED(gr_write(gr_out, elem, ctx));
            res += gr_out->len;
            ip += STRING_LENGTH("gr}");
        }
        else if (IS_FLINT_TYPE(ip, "gr*"))
        {
            gr_srcptr elem = va_arg(vlist, gr_srcptr);
            slong len = va_arg(vlist, slong);
            gr_ctx_struct * ctx = va_arg(vlist, gr_ctx_struct *);
            GR_MUST_SUCCEED(_gr_vec_write(gr_out, elem, len, ctx));
            res += gr_out->len;
            ip += STRING_LENGTH("gr*}");
        }
        else if (IS_FLINT_TYPE(ip, "gr_poly"))
        {
            const gr_poly_struct * elem = va_arg(vlist, const gr_poly_struct *);
            gr_ctx_struct * ctx = va_arg(vlist, gr_ctx_struct *);
            GR_MUST_SUCCEED(gr_poly_write(gr_out, elem, "x", ctx));
            res += gr_out->len;
            ip += STRING_LENGTH("gr_poly}");
        }
        else if (IS_FLINT_TYPE(ip, "gr_ore_poly"))
        {
            const gr_ore_poly_struct * elem = va_arg(vlist, const gr_ore_poly_struct *);
            gr_ctx_struct * ctx = va_arg(vlist, gr_ctx_struct *);
            GR_MUST_SUCCEED(gr_ore_poly_write(gr_out, elem, ctx));
            res += gr_out->len;
            ip += STRING_LENGTH("gr_ore_poly}");
        }
        else if (IS_FLINT_TYPE(ip, "gr_mat"))
        {
            const gr_mat_struct * elem = va_arg(vlist, const gr_mat_struct *);
            gr_ctx_struct * ctx = va_arg(vlist, gr_ctx_struct *);
            GR_MUST_SUCCEED(_gr_mat_write(gr_out, elem, 0, ctx));
            res += gr_out->len;
            ip += STRING_LENGTH("gr_mat}");
        }
        else if (IS_FLINT_TYPE(ip, "gr_ctx"))
        {
            gr_ctx_struct * ctx = va_arg(vlist, gr_ctx_struct *);
            GR_MUST_SUCCEED(gr_ctx_write(gr_out, ctx));
            res += gr_out->len;
            ip += STRING_LENGTH("gr_ctx}");
        }
        else
            goto printpercentcurlybracket;

        FLINT_VPRINTF_GR_STREAM_FLUSH(gr_out, out);
    }
    else
    {
printpercentcurlybracket:
        /* Invalid use of "%{FLINT_TYPE}". As we are currently pointed to
         * "FLINT_TYPE}", we let FLINT_VPRINTF_PRINTF take care of printing "%{". */
DIAGNOSTIC_PUSH
DIAGNOSTIC_IGNORE_FORMAT
        tmp = FLINT_VPRINTF_PRINTF(out, "%{");
DIAGNOSTIC_POP
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

/* Type specific printing ****************************************************/

/* TODO: Move these to their respective module? */

/* TODO: Add precision input to Arb type printing functions */

/* TODO: Add option to print in different basis? */

static size_t __ulong_print(FLINT_VPRINTF_OUT_T * out, const ulong * ip, int FLINT_UNUSED(flag))
{
    return FLINT_VPRINTF_PRINTF(out, WORD_FMT "u", *ip);
}

static size_t __slong_print(FLINT_VPRINTF_OUT_T * out, const slong * ip, int FLINT_UNUSED(flag))
{
    return FLINT_VPRINTF_PRINTF(out, WORD_FMT "d", *ip);
}

#define BASE 10
static size_t __fmpz_print(FLINT_VPRINTF_OUT_T * out, const fmpz * ip, int flag)
{
    size_t res = 0;
    char * str;
    size_t skipminus = __fmpz_is_neg(ip) ? FLAG_IS_NEG(flag) : 0;

    str = fmpz_get_str(NULL, BASE, ip);
    res += FLINT_VPRINTF_WRITE(str + skipminus, strlen(str + skipminus), out);

    flint_free(str);

    return res;
}
#undef BASE

static size_t __mpz_print(FLINT_VPRINTF_OUT_T * out, mpz_srcptr ip)
{
    fmpz_t fip;

    if ((ip->_mp_size == 1 || ip->_mp_size == -1) && ip->_mp_d[0] <= COEFF_MAX)
        *fip = ip->_mp_size > 0 ? ip->_mp_d[0] : -ip->_mp_d[0];
    else
        *fip = PTR_TO_COEFF(ip);

    return __fmpz_print(out, fip, FLAG_NONE);
}

static size_t __fmpq_print(FLINT_VPRINTF_OUT_T * out, const fmpq * ip, int flag)
{
    size_t res = 0;

    /* NOTE: We do not care about parentheses here, and __fmpz_print does not
     * care either. */
    res += __fmpz_print(out, fmpq_numref(ip), flag);
    if (!fmpz_is_one(fmpq_denref(ip)))
    {
        res += FLINT_VPRINTF_WRITE(" / ", STRING_LENGTH(" / "), out);
        res += __fmpz_print(out, fmpq_denref(ip), FLAG_NONE);
    }

    return res;
}

static size_t __mpq_print(FLINT_VPRINTF_OUT_T * out, mpq_srcptr ip)
{
    fmpq_t fip;

    if ((ip->_mp_num._mp_size == 1 || ip->_mp_num._mp_size == -1) && ip->_mp_num._mp_d[0] <= COEFF_MAX)
        *fmpq_numref(fip) = ip->_mp_num._mp_size > 0 ? ip->_mp_num._mp_d[0] : -ip->_mp_num._mp_d[0];
    else
        *fmpq_numref(fip) = PTR_TO_COEFF(&(ip->_mp_num));

    if ((ip->_mp_den._mp_size == 1 || ip->_mp_den._mp_size == -1) && ip->_mp_den._mp_d[0] <= COEFF_MAX)
        *fmpq_denref(fip) = ip->_mp_den._mp_size > 0 ? ip->_mp_den._mp_d[0] : -ip->_mp_den._mp_d[0];
    else
        *fmpq_denref(fip) = PTR_TO_COEFF(&(ip->_mp_den));

    return __fmpq_print(out, fip, FLAG_NONE);
}

#define DIGITS 6
static size_t __mag_print(FLINT_VPRINTF_OUT_T * out, mag_srcptr ip, int FLINT_UNUSED(flag))
{
    size_t res;
    char * str;

#if DIGITS != 6
# error Change FLINT_VPRINTF_WRITE below
#endif
    if (mag_is_zero(ip))
        return FLINT_VPRINTF_WRITE("0.00000", STRING_SIZE("0.00000"), out);

    str = __mag_get_str(ip, DIGITS);
    res = FLINT_VPRINTF_WRITE(str, strlen(str), out);

    flint_free(str);

    return res;
}

static size_t __arf_print(FLINT_VPRINTF_OUT_T * out, arf_srcptr ip, int FLINT_UNUSED(flag))
{
    size_t res;
    char * str;

#if DIGITS != 6
# error Change FLINT_VPRINTF_WRITE below
#endif
    if (arf_is_zero(ip))
        return FLINT_VPRINTF_WRITE("0.00000", STRING_SIZE("0.00000"), out);

    str = arf_get_str(ip, DIGITS);
    res = FLINT_VPRINTF_WRITE(str, strlen(str), out);

    flint_free(str);

    return res;
}

#define MAX_INT_SIZE 64
/* NOTE: If arb is an integer, we print it as one. */
static size_t __arb_print(FLINT_VPRINTF_OUT_T * out, arb_srcptr ip, int flag)
{
    size_t res;

    if (arb_is_zero(ip))
    {
        return (FLINT_VPRINTF_PUTC('0', out) != FLINT_VPRINTF_PUTC_ERRVAL);
    }
    else if (arb_is_int(ip) && ARF_EXP(arb_midref(ip)) <= MAX_INT_SIZE)
    {
        /* NOTE: Only print as integer if ip < 2^64. */
        fmpz_t fip;

        fmpz_init(fip);

        /* NOTE: The conversion should be exact, so we do not care about its
         * return value. */
        arf_get_fmpz(fip, arb_midref(ip), ARF_RND_DOWN);
        /* NOTE: We do not care about parentheses here, and __fmpz_print does not
         * care either. */
        res = __fmpz_print(out, fip, flag);

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
        res = FLINT_VPRINTF_WRITE(str, strlen(str), out);

        flint_free(str);
    }

    return res;
}
#undef MAX_INT_SIZE
#undef DIGITS

static size_t __acb_print(FLINT_VPRINTF_OUT_T * out, acb_srcptr ip, int flag)
{
    size_t res = 0;
    int realiszero, imagiszero;

    realiszero = arb_is_zero(acb_realref(ip));
    imagiszero = arb_is_zero(acb_imagref(ip));

    if (realiszero && imagiszero)
        return (FLINT_VPRINTF_PUTC('0', out) != FLINT_VPRINTF_PUTC_ERRVAL);

    /* Only print parentheses if both real and imaginary part is non-zero. */
    if (FLAG_IS_PAREN(flag) && !realiszero && !imagiszero)
        res += (FLINT_VPRINTF_PUTC('(', out) != FLINT_VPRINTF_PUTC_ERRVAL);

    /* Print real part if non-zero */
    if (!realiszero)
        res += __arb_print(out, acb_realref(ip), flag);

    /* Print imaginary part if non-zero */
    if (!imagiszero)
    {
        int imagisneg = arb_is_negative(acb_imagref(ip));

        if (!realiszero)
            res += FLINT_VPRINTF_WRITE(FLAG_IS_NEG(flag) ^ imagisneg ? " - " : " + ", STRING_SIZE(" - "), out);

        /* If imaginary part of ip is \pm 1, then we only print \pm i. */
        if (!__arb_is_pm1(acb_imagref(ip)))
        {
            res += __arb_print(out, acb_imagref(ip), (FLAG_IS_NEG(flag) ^ imagisneg) & (~realiszero));
            res += FLINT_VPRINTF_WRITE(" * ", STRING_LENGTH(" * "), out);
        }

        res += (FLINT_VPRINTF_PUTC('i', out) != FLINT_VPRINTF_PUTC_ERRVAL);
    }

    if (FLAG_IS_PAREN(flag) && !realiszero && !imagiszero)
        res += (FLINT_VPRINTF_PUTC(')', out) != FLINT_VPRINTF_PUTC_ERRVAL);

    return res;
}

static size_t __nmod_print(FLINT_VPRINTF_OUT_T * out, nmod_t ip)
{
    return FLINT_VPRINTF_PRINTF(out, "mod " WORD_FMT "u", ip.n);
}

static size_t __fmpz_mod_ctx_print(FLINT_VPRINTF_OUT_T * out, const fmpz_mod_ctx_struct * ip)
{
    size_t res = 0;

    res += FLINT_VPRINTF_WRITE("mod ", STRING_LENGTH("mod "), out);
    res += __fmpz_print(out, fmpz_mod_ctx_modulus(ip), FLAG_NONE);

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

typedef size_t (* print_func_t)(FLINT_VPRINTF_OUT_T *, const void *, int);

static print_func_t print_functions[] =
{
    (print_func_t) __ulong_print,
    (print_func_t) __fmpz_print,
    (print_func_t) __fmpq_print,
    (print_func_t) __arb_print,
    (print_func_t) __acb_print,
    (print_func_t) __slong_print, /* NOTE: These print functions are only  */
    (print_func_t) __mag_print,   /* used for vectors, not other composite */
    (print_func_t) __arf_print    /* types. */
};

static size_t __flint_vec_print(FLINT_VPRINTF_OUT_T * out, const void * ip, slong len, flint_type_t type)
{
    size_t res = 0;
    slong ix;
    const char * vec = ip;
    size_t type_size = flint_type_size_in_chars(type);
    print_func_t print = print_functions[type];

    res += (FLINT_VPRINTF_PUTC('[', out) != FLINT_VPRINTF_PUTC_ERRVAL);

    if (len > 0)
        res += print(out, vec, FLAG_NONE);

    for (ix = 1; ix < len; ix++)
    {
        res += FLINT_VPRINTF_WRITE(", ", STRING_LENGTH(", "), out);
        res += print(out, vec + type_size * ix, FLAG_NONE);
    }

    res += (FLINT_VPRINTF_PUTC(']', out) != FLINT_VPRINTF_PUTC_ERRVAL);

    return res;
}

/* NOTE: This function relies on the fact that the layout of
 * [fmpz/fmpq/arb/acb]_mat_struct are all on the form (pointer, slong, slong,
 * pointer). */
static size_t __flint_mat_print(FLINT_VPRINTF_OUT_T * out, const void * ip, flint_type_t type)
{
    size_t res = 0;
    slong ix;
    slong nr, nc;
    const char * entries;
    slong stride;

    entries = (const char *) ((const fmpz_mat_struct *) ip)->entries;
    nr = ((const fmpz_mat_struct *) ip)->r;
    nc = ((const fmpz_mat_struct *) ip)->c;
    stride = ((const fmpz_mat_struct *) ip)->stride * flint_type_size_in_chars(type);

    if (nr == 0 || nc == 0)
        return FLINT_VPRINTF_PRINTF(out, WORD_FMT "d by " WORD_FMT "d empty matrix", nr, nc);

    res += (FLINT_VPRINTF_PUTC('[', out) != FLINT_VPRINTF_PUTC_ERRVAL);
    res += __flint_vec_print(out, entries, nc, type);

    for (ix = 1; ix < nr; ix++)
    {
        res += FLINT_VPRINTF_WRITE(", ", STRING_LENGTH(", "), out);
        res += __flint_vec_print(out, entries + ix * stride, nc, type);
    }

    res += (FLINT_VPRINTF_PUTC(']', out) != FLINT_VPRINTF_PUTC_ERRVAL);

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
static size_t __flint_poly_print(FLINT_VPRINTF_OUT_T * out, const void * ip, flint_type_t type)
{
    size_t res = 0;
    slong ix;
    slong len;

    len = ((const fmpz_poly_struct *) ip)->length;

    if (len == 0)
        return FLINT_VPRINTF_PUTC('0', out) != FLINT_VPRINTF_PUTC_ERRVAL;

    if (type != fmpq_type)
    {
        size_t type_size = flint_type_size_in_chars(type);
        const char * coeffs = (const char *) ((const fmpz_poly_struct *) ip)->coeffs;
        print_func_t print = print_functions[type];
        is_func_t is_pm1 = is_pm1_functions[type];
        is_func_t is_zero = is_zero_functions[type];
        is_func_t is_neg = is_neg_functions[type];

        if (len == 1)
            return print(out, coeffs, FLAG_NONE);

        /* Leading coefficient cannot be zero */
        if (!is_pm1(coeffs + type_size * (len - 1)))
        {
            res += print(out, coeffs + type_size * (len - 1), FLAG_PAREN);
            res += FLINT_VPRINTF_WRITE(" * ", STRING_LENGTH(" * "), out);
        }
        else if (is_neg(coeffs + type_size * (len - 1)))
            res += (FLINT_VPRINTF_PUTC('-', out) != FLINT_VPRINTF_PUTC_ERRVAL);
        if (len != 2)
            res += FLINT_VPRINTF_PRINTF(out, "x^" WORD_FMT "d", len - 1);
        else
            res += (FLINT_VPRINTF_PUTC('x', out) != FLINT_VPRINTF_PUTC_ERRVAL);

        for (ix = len - 2; ix > 0; ix--)
        {
            if (!is_zero(coeffs + type_size * ix))
            {
                res += FLINT_VPRINTF_WRITE(is_neg(coeffs + type_size * ix) ? " - " : " + ", STRING_LENGTH(" - "), out);
                if (!is_pm1(coeffs + type_size * ix))
                {
                    res += print(out, coeffs + type_size * ix, is_neg(coeffs + type_size * ix) | FLAG_PAREN);
                    res += FLINT_VPRINTF_WRITE(" * ", STRING_LENGTH(" * "), out);
                }
                if (ix != 1)
                    res += FLINT_VPRINTF_PRINTF(out, "x^" WORD_FMT "d", ix);
                else
                    res += (FLINT_VPRINTF_PUTC('x', out) != FLINT_VPRINTF_PUTC_ERRVAL);
            }
        }

        if (!is_zero(coeffs + 0))
        {
            res += FLINT_VPRINTF_WRITE(is_neg(coeffs + 0) ? " - " : " + ", STRING_LENGTH(" - "), out);
            res += print(out, coeffs + 0, is_neg(coeffs + 0));
        }
    }
    else
    {
        /* fmpq_poly is special as it is an fmpz_poly with a denominator
         * strapped onto it */
        const fmpz * coeffs = ((const fmpq_poly_struct *) ip)->coeffs;
        const fmpz * den = ((const fmpq_poly_struct *) ip)->den;
        fmpq_t canonical;

        fmpq_init(canonical);
        fmpq_set_fmpz_frac(canonical, coeffs + len - 1, den);

        if (len == 1)
        {
            res += __fmpq_print(out, canonical, FLAG_NONE);
            fmpq_clear(canonical);
            return res;
        }

        /* Leading coefficient cannot be zero */
        if (!fmpq_is_pm1(canonical))
        {
            res += __fmpq_print(out, canonical, FLAG_NONE);
            res += FLINT_VPRINTF_WRITE(" * ", STRING_LENGTH(" * "), out);
        }
        else if (__fmpq_is_neg(canonical))
            res += (FLINT_VPRINTF_PUTC('-', out) != FLINT_VPRINTF_PUTC_ERRVAL);
        if (len != 2)
            res += FLINT_VPRINTF_PRINTF(out, "x^" WORD_FMT "d", len - 1);
        else
            res += (FLINT_VPRINTF_PUTC('x', out) != FLINT_VPRINTF_PUTC_ERRVAL);

        for (ix = len - 2; ix > 0; ix--)
        {
            if (!fmpz_is_zero(coeffs + ix))
            {
                fmpq_set_fmpz_frac(canonical, coeffs + ix, den);
                res += FLINT_VPRINTF_WRITE(__fmpq_is_neg(canonical) ? " - " : " + ", STRING_LENGTH(" - "), out);
                if (!fmpq_is_pm1(canonical))
                {
                    res += __fmpq_print(out, canonical, __fmpq_is_neg(canonical));
                    res += FLINT_VPRINTF_WRITE(" * ", STRING_LENGTH(" * "), out);
                }
                if (ix != 1)
                    res += FLINT_VPRINTF_PRINTF(out, "x^" WORD_FMT "d", ix);
                else
                    res += (FLINT_VPRINTF_PUTC('x', out) != FLINT_VPRINTF_PUTC_ERRVAL);
            }
        }

        fmpq_set_fmpz_frac(canonical, coeffs + 0, den);
        if (!fmpq_is_zero(canonical))
        {
            res += FLINT_VPRINTF_WRITE(__fmpq_is_neg(canonical) ? " - " : " + ", STRING_LENGTH(" - "), out);
            res += __fmpq_print(out, canonical, __fmpq_is_neg(canonical));
        }

        fmpq_clear(canonical);
    }

    return res;
}
