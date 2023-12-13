/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdint.h> /* intmax_t */
#include <string.h>
#include <wchar.h> /* wchar_t and wint_t */
#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpq.h"
#include "mag.h"
#include "arf.h"
#include "arb.h"
#include "acb.h"
#include "nmod_vec.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "nmod_mat.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "nmod_poly.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "arb_poly.h"
#include "acb_poly.h"

#define STR(x) TEMPLATE_STR(x)

static void check_strings(const char * str1, const char * str2)
{
    const char * str1cur, * str2cur;
    size_t lnsz1, lnsz2;
    int line;

    /* Check line by line to make debugging easier */
    str1cur = str1;
    str2cur = str2;
    line = 0;
    while (1)
    {
        /* NOTE: The newline or null character will occur at strXcur[lnszX] */
        lnsz1 = strcspn(str1cur, "\n");
        lnsz2 = strcspn(str2cur, "\n");

        if (lnsz1 != lnsz2)
            flint_throw(FLINT_TEST_FAIL,
                    "Line size at line %d differs.\n"
                    "Expected: %zu\n"
                    "Got:      %zu\n"
                    "\n"
                    "Expected string at line %d:\n"
                    "\n"
                    "%.*s\n"
                    "\n"
                    "Got string at line %d:\n"
                    "\n"
                    "%.*s\n",
                    line,
                    lnsz1,
                    lnsz2,
                    line,
                    lnsz1, str1cur,
                    line,
                    lnsz2, str2cur);

        if (memcmp(str1cur, str2cur, sizeof(char) * lnsz1))
            flint_throw(FLINT_TEST_FAIL,
                    "String at line %d differs.\n"
                    "\n"
                    "Expected:\n"
                    "\n"
                    "%.*s\n"
                    "\n"
                    "Got:\n"
                    "\n"
                    "%.*s\n",
                    line,
                    lnsz1, str1cur,
                    lnsz2, str2cur);

        if (str1cur[lnsz1] == '\0' && str2cur[lnsz2] == '\0')
            break;
        else if (str1cur[lnsz1] == '\0' || str2cur[lnsz2] == '\0')
            flint_throw(FLINT_TEST_FAIL,
                    "Reached end of string in str%d at line %d.\n"
                    "\n"
                    "Expected string at line %d:\n"
                    "\n"
                    "%.*s\n"
                    "\n"
                    "Got string at line %d:\n"
                    "\n"
                    "%.*s\n",
                    str1cur[lnsz1] == '\0' ? 1 : 2, line,
                    line,
                    lnsz1, str1cur,
                    line,
                    lnsz2, str2cur);

        str1cur += lnsz1 + 1;
        str2cur += lnsz2 + 1;
        line++;
    }
}

/* Basetypes *****************************************************************/
#define NMOD_MOD 5
#define NMOD_INIT(x) nmod_init(&(x), NMOD_MOD)
#define NMOD_CLEAR(x) do { } while (0)
#define NMOD_SET(x) do { } while (0)
#define NMOD_STRING "mod " STR(NMOD_MOD)

#define FMPZ1_ENTRY 1337
#define FMPZ1_INIT(x) fmpz_init(x)
#define FMPZ1_CLEAR(x) fmpz_clear(x)
#define FMPZ1_SET(x) fmpz_set_si(x, FMPZ1_ENTRY)
#define FMPZ1_STRING STR(FMPZ1_ENTRY)

#define FMPZ2_ENTRY -187236816283761872638123
#define FMPZ2_INIT(x) fmpz_init(x)
#define FMPZ2_CLEAR(x) fmpz_clear(x)
#define FMPZ2_SET(x) fmpz_set_str(x, STR(FMPZ2_ENTRY), 10)
#define FMPZ2_STRING STR(FMPZ2_ENTRY)

#define FMPZ_MOD_CTX_ENTRY 808
#define FMPZ_MOD_CTX_INIT(x) fmpz_mod_ctx_init_ui(x, FMPZ_MOD_CTX_ENTRY)
#define FMPZ_MOD_CTX_CLEAR(x) fmpz_mod_ctx_clear(x)
#define FMPZ_MOD_CTX_SET(x) do { } while (0)
#define FMPZ_MOD_CTX_STRING "mod " STR(FMPZ_MOD_CTX_ENTRY)

#define FMPQ1_ENTRY 1337
#define FMPQ1_INIT(x) fmpq_init(x)
#define FMPQ1_CLEAR(x) fmpq_clear(x)
#define FMPQ1_SET(x) fmpq_set_si(x, FMPQ1_ENTRY, 1)
#define FMPQ1_STRING STR(FMPQ1_ENTRY)

#define FMPQ2_ENTRY_N -1024
#define FMPQ2_ENTRY_D 3
#define FMPQ2_INIT(x) fmpq_init(x)
#define FMPQ2_CLEAR(x) fmpq_clear(x)
#define FMPQ2_SET(x) fmpq_set_si(x, FMPQ2_ENTRY_N, FMPQ2_ENTRY_D)
#define FMPQ2_STRING STR(FMPQ2_ENTRY_N) " / " STR(FMPQ2_ENTRY_D)

#define ARF1_INIT(x) arf_init(x)
#define ARF1_CLEAR(x) arf_clear(x)
#define ARF1_SET(x) do { } while (0)
#define ARF1_STRING "0.00000"

#define ARF2_ENTRY 3.00000
#define ARF2_INIT(x) arf_init(x)
#define ARF2_CLEAR(x) arf_clear(x)
#define ARF2_SET(x) arf_set_d(x, ARF2_ENTRY)
#define ARF2_STRING STR(ARF2_ENTRY)

#define ARF3_ENTRY 3.14150
#define ARF3_INIT(x) arf_init(x)
#define ARF3_CLEAR(x) arf_clear(x)
#define ARF3_SET(x) arf_set_d(x, ARF3_ENTRY)
#define ARF3_STRING STR(ARF3_ENTRY)

#define MAG1_INIT(x) mag_init(x)
#define MAG1_CLEAR(x) mag_clear(x)
#define MAG1_SET(x) do { } while (0)
#define MAG1_STRING "0.00000"

#define MAG2_ENTRY 3.00000
#define MAG2_INIT(x) mag_init(x)
#define MAG2_CLEAR(x) mag_clear(x)
#define MAG2_SET(x) mag_set_d(x, MAG2_ENTRY)
#define MAG2_STRING STR(MAG2_ENTRY)

#define MAG3_ENTRY 3.14150
#define MAG3_INIT(x) mag_init(x)
#define MAG3_CLEAR(x) mag_clear(x)
#define MAG3_SET(x) mag_set_d(x, MAG3_ENTRY)
#define MAG3_STRING STR(MAG3_ENTRY)

#define ARB1_INIT(x) arb_init(x)
#define ARB1_CLEAR(x) arb_clear(x)
#define ARB1_SET(x) do { } while (0)
#define ARB1_STRING "0"

#define ARB2_ENTRY 3
#define ARB2_INIT(x) arb_init(x)
#define ARB2_CLEAR(x) arb_clear(x)
#define ARB2_SET(x) arb_set_si(x, ARB2_ENTRY)
#define ARB2_STRING STR(ARB2_ENTRY)

#define ARB3_ENTRY "3.1415926535897 +/- 0.0000000000001"
#define ARB3_PREC 30
#define ARB3_INIT(x) arb_init(x)
#define ARB3_CLEAR(x) arb_clear(x)
#define ARB3_SET(x) arb_set_str(x, ARB3_ENTRY, 30)
#define ARB3_STRING "[3.14159 +/- 2.66e-6]"

#define ACB1_INIT(x) acb_init(x)
#define ACB1_CLEAR(x) acb_clear(x)
#define ACB1_SET(x) do { } while (0)
#define ACB1_STRING "0"

#define ACB2_ENTRY 3
#define ACB2_INIT(x) acb_init(x)
#define ACB2_CLEAR(x) acb_clear(x)
#define ACB2_SET(x) acb_set_si(x, ACB2_ENTRY)
#define ACB2_STRING STR(ACB2_ENTRY)

#define ACB3_ENTRY 3
#define ACB3_INIT(x) acb_init(x)
#define ACB3_CLEAR(x) acb_clear(x)
#define ACB3_SET(x) arb_set_si(acb_imagref(x), ACB3_ENTRY)
#define ACB3_STRING STR(ACB3_ENTRY) " * i"

#define ACB4_ENTRY -3
#define ACB4_INIT(x) acb_init(x)
#define ACB4_CLEAR(x) acb_clear(x)
#define ACB4_SET(x) arb_set_si(acb_imagref(x), ACB4_ENTRY)
#define ACB4_STRING STR(ACB4_ENTRY) " * i"

#define ACB5_ENTRY_RE 5
#define ACB5_ENTRY_IM -3
#define ACB5_INIT(x) acb_init(x)
#define ACB5_CLEAR(x) acb_clear(x)
#define ACB5_SET(x) \
do                  \
{                   \
    arb_set_si(acb_realref(x), ACB5_ENTRY_RE); \
    arb_set_si(acb_imagref(x), ACB5_ENTRY_IM); \
} while (0)
#define ACB5_STRING "5 - 3 * i"

#define ACB6_ENTRY "3.1415926535897 +/- 0.0000000000001"
#define ACB6_PREC 30
#define ACB6_INIT(x) acb_init(x)
#define ACB6_CLEAR(x) acb_clear(x)
#define ACB6_SET(x) arb_set_str(acb_realref(x), ACB6_ENTRY, ACB6_PREC)
#define ACB6_STRING "[3.14159 +/- 2.66e-6]"

#define ACB7_ENTRY "-3.1415926535897 +/- 0.0000000000001"
#define ACB7_PREC 30
#define ACB7_INIT(x) acb_init(x)
#define ACB7_CLEAR(x) acb_clear(x)
#define ACB7_SET(x) arb_set_str(acb_imagref(x), ACB7_ENTRY, ACB7_PREC)
#define ACB7_STRING "[-3.14159 +/- 2.66e-6] * i"

#define ACB8_ENTRY_RE "-1.337 +/- 0.001"
#define ACB8_ENTRY_IM "-3.1415926535897 +/- 0.0000000000001"
#define ACB8_PREC 30
#define ACB8_INIT(x) acb_init(x)
#define ACB8_CLEAR(x) acb_clear(x)
#define ACB8_SET(x) \
do                  \
{                   \
    arb_set_str(acb_realref(x), ACB8_ENTRY_RE, ACB8_PREC); \
    arb_set_str(acb_imagref(x), ACB8_ENTRY_IM, ACB8_PREC); \
} while (0)
#define ACB8_STRING "[-1.34 +/- 4.01e-3] - [3.14159 +/- 2.66e-6] * i"

/* Vectors *******************************************************************/
/* NOTE: The lengths has to be put into flint_fprintf as slongs. GCC and MSVC
 * handles it without specifying them as slongs, but at least Clang messes it
 * up. */
#define SLONG_VEC_LEN WORD(4)
#define SLONG_VEC_INIT(x) do { } while (0)
#define SLONG_VEC_CLEAR(x) do { } while (0)
#define SLONG_VEC_ENTRIES ((slong[]) {WORD(-1), WORD(0), WORD(1), WORD(2)})
#define SLONG_VEC_SET(x) memcpy((x), SLONG_VEC_ENTRIES, sizeof(slong) * SLONG_VEC_LEN)
#define SLONG_VEC_STRING "[-1, 0, 1, 2]"

#define NMOD_VEC_LEN WORD(3)
#define NMOD_VEC_INIT(x) do { (x) = _nmod_vec_init(NMOD_VEC_LEN); } while (0)
#define NMOD_VEC_CLEAR(x) _nmod_vec_clear(x)
#define NMOD_VEC_ENTRIES ((ulong[]) {UWORD(1), UWORD(2), UWORD(3)})
#define NMOD_VEC_SET(x) memcpy((x), NMOD_VEC_ENTRIES, sizeof(ulong) * NMOD_VEC_LEN)
#define NMOD_VEC_STRING "[1, 2, 3]"

#define FMPZ_VEC_LEN WORD(2)
#define FMPZ_VEC_INIT(x) do { (x) = _fmpz_vec_init(FMPZ_VEC_LEN); } while (0)
#define FMPZ_VEC_CLEAR(x) _fmpz_vec_clear(x, FMPZ_VEC_LEN)
#define FMPZ_VEC_ENTRIES ((fmpz[]) {WORD(10), WORD(-2)})
#define FMPZ_VEC_SET(x) memcpy(x, FMPZ_VEC_ENTRIES, sizeof(fmpz) * FMPZ_VEC_LEN)
#define FMPZ_VEC_STRING "[10, -2]"

/* NOTE: The entries has to be in a canonical form */
#define FMPQ_VEC_LEN WORD(3)
#define FMPQ_VEC_INIT(x) do { (x) = _fmpq_vec_init(FMPQ_VEC_LEN); } while (0)
#define FMPQ_VEC_CLEAR(x) _fmpq_vec_clear(x, FMPQ_VEC_LEN)
#define FMPQ_VEC_ENTRIES ((fmpq[]) {{WORD(10), WORD(3)}, {WORD(-2), WORD(1)}, {WORD(-2000), WORD(7)}})
#define FMPQ_VEC_SET(x) memcpy(x, FMPQ_VEC_ENTRIES, sizeof(fmpq) * FMPQ_VEC_LEN)
#define FMPQ_VEC_STRING "[10 / 3, -2, -2000 / 7]"

#define ARB_VEC_LEN WORD(2)
#define ARB_VEC_INIT(x) do { (x) = _arb_vec_init(ARB_VEC_LEN); } while (0)
#define ARB_VEC_CLEAR(x) _arb_vec_clear(x, ARB_VEC_LEN)
#define ARB_VEC_PREC 30
#define ARB_VEC_ENTRY(ix) (ix == 0) ? "12.1 +/- 1" : "3"
#define ARB_VEC_SET(x)                      \
do                                          \
{                                           \
    slong jx;                               \
    for (jx = 0; jx < ARB_VEC_LEN; jx++)    \
        arb_set_str((x) + jx, ARB_VEC_ENTRY(jx), ARB_VEC_PREC); \
} while (0)
#define ARB_VEC_STRING "[[1e+1 +/- 3.11], 3]"

#define ACB_VEC_LEN WORD(1)
#define ACB_VEC_INIT(x) do { (x) = _acb_vec_init(ACB_VEC_LEN); } while (0)
#define ACB_VEC_CLEAR(x) _acb_vec_clear(x, ACB_VEC_LEN)
#define ACB_VEC_PREC 30
#define ACB_VEC_ENTRY_REAL(ix) "12.1 +/- 1"
#define ACB_VEC_ENTRY_IMAG(ix) "3"
#define ACB_VEC_SET(x)                      \
do                                          \
{                                           \
    slong jx;                               \
    for (jx = 0; jx < ACB_VEC_LEN; jx++)    \
    {                                       \
        arb_set_str(acb_realref((x) + jx), ACB_VEC_ENTRY_REAL(jx), ACB_VEC_PREC); \
        arb_set_str(acb_imagref((x) + jx), ACB_VEC_ENTRY_IMAG(jx), ACB_VEC_PREC); \
    }                                       \
} while (0)
#define ACB_VEC_STRING "[[1e+1 +/- 3.11] + 3 * i]"

/* Matrices ******************************************************************/
#define FMPZ_MAT_EMPTY_R 8
#define FMPZ_MAT_EMPTY_C 0
#define FMPZ_MAT_EMPTY_INIT(mat) fmpz_mat_init(mat, FMPZ_MAT_EMPTY_R, FMPZ_MAT_EMPTY_C)
#define FMPZ_MAT_EMPTY_CLEAR(mat) fmpz_mat_clear(mat)
#define FMPZ_MAT_EMPTY_SET(mat) do { } while (0)
#define FMPZ_MAT_EMTPY_STRING STR(FMPZ_MAT_EMPTY_R) " by " STR(FMPZ_MAT_EMPTY_C) " empty matrix"

#define NMOD_MAT_R 3
#define NMOD_MAT_C 2
#define NMOD_MAT_MOD 8
#define NMOD_MAT_WINDOW_R1 1
#define NMOD_MAT_WINDOW_R2 3
#define NMOD_MAT_WINDOW_C1 1
#define NMOD_MAT_WINDOW_C2 2
#define NMOD_MAT_WITH_WINDOW_INIT(mat, mat_window)  \
do                                                  \
{                                                   \
    nmod_mat_init(mat, NMOD_MAT_R, NMOD_MAT_C, NMOD_MAT_MOD); \
    nmod_mat_window_init(mat_window, mat,           \
            NMOD_MAT_WINDOW_R1, NMOD_MAT_WINDOW_C1, \
            NMOD_MAT_WINDOW_R2, NMOD_MAT_WINDOW_C2);\
} while (0)
#define NMOD_MAT_WITH_WINDOW_CLEAR(mat, mat_window) \
do                                                  \
{                                                   \
    nmod_mat_clear(mat);                            \
    nmod_mat_window_clear(mat_window);              \
} while (0)
#define NMOD_MAT_WITH_WINDOW_ENTRIES \
    ((ulong[]) {UWORD(0), UWORD(1), UWORD(2), UWORD(3), UWORD(4), UWORD(5)})
#define NMOD_MAT_WITH_WINDOW_SET(mat, mat_window) \
    memcpy((mat)->entries, NMOD_MAT_WITH_WINDOW_ENTRIES, sizeof(ulong) * NMOD_MAT_R * NMOD_MAT_C)
#define NMOD_MAT_WITH_WINDOW_STRING \
    "[[3], [5]]"

#define FMPZ_MAT_R 2
#define FMPZ_MAT_C 3
#define FMPZ_MAT_INIT(mat) fmpz_mat_init(mat, FMPZ_MAT_R, FMPZ_MAT_C)
#define FMPZ_MAT_CLEAR(mat) fmpz_mat_clear(mat)
#define FMPZ_MAT_ENTRIES \
    ((fmpz[]) {WORD(0), WORD(1), WORD(2), WORD(3), WORD(4), WORD(5)})
#define FMPZ_MAT_SET(mat) \
    memcpy((mat)->entries, FMPZ_MAT_ENTRIES, sizeof(fmpz) * FMPZ_MAT_R * FMPZ_MAT_C)
#define FMPZ_MAT_STRING \
    "[[0, 1, 2], [3, 4, 5]]"

/* Polynomials ***************************************************************/
#define NMOD_POLY_ZERO_MOD 6
#define NMOD_POLY_ZERO_INIT(x) nmod_poly_init(x, NMOD_POLY_ZERO_MOD)
#define NMOD_POLY_ZERO_CLEAR(x) nmod_poly_clear(x)
#define NMOD_POLY_ZERO_SET(x) do { } while (0)
#define NMOD_POLY_ZERO_STRING "0"

#define NMOD_POLY_CONSTANT_MOD 1338
#define NMOD_POLY_CONSTANT_COEFF 1337
#define NMOD_POLY_CONSTANT_INIT(x) nmod_poly_init2(x, NMOD_POLY_CONSTANT_MOD, 1)
#define NMOD_POLY_CONSTANT_CLEAR(x) nmod_poly_clear(x)
#define NMOD_POLY_CONSTANT_SET(x) nmod_poly_set_coeff_ui(x, 0, NMOD_POLY_CONSTANT_COEFF)
#define NMOD_POLY_CONSTANT_STRING STR(NMOD_POLY_CONSTANT_COEFF)

#define NMOD_POLY_MOD 9
#define NMOD_POLY_LENGTH 4
#define NMOD_POLY_INIT(x) nmod_poly_init2(x, NMOD_POLY_MOD, NMOD_POLY_LENGTH)
#define NMOD_POLY_CLEAR(x) nmod_poly_clear(x)
#define NMOD_POLY_ENTRIES \
    ((ulong[]) {UWORD(3), UWORD(0), UWORD(1), UWORD(3)})
#define NMOD_POLY_SET(x)    \
do                          \
{                           \
    _nmod_poly_set_length(x, NMOD_POLY_LENGTH); \
    memcpy((x)->coeffs, NMOD_POLY_ENTRIES, sizeof(ulong) * NMOD_POLY_LENGTH); \
} while (0)
#define NMOD_POLY_STRING "3 * x^3 + x^2 + 3"

#define FMPZ_POLY_LENGTH 4
#define FMPZ_POLY_INIT(x) fmpz_poly_init2(x, FMPZ_POLY_LENGTH)
#define FMPZ_POLY_CLEAR(x) fmpz_poly_clear(x)
#define FMPZ_POLY_ENTRIES \
    ((fmpz[]) {WORD(0), WORD(-3), WORD(-1), WORD(-1)})
#define FMPZ_POLY_SET(x)    \
do                          \
{                           \
    _fmpz_poly_set_length(x, FMPZ_POLY_LENGTH); \
    memcpy((x)->coeffs, FMPZ_POLY_ENTRIES, sizeof(fmpz) * FMPZ_POLY_LENGTH); \
} while (0)
#define FMPZ_POLY_STRING "-x^3 - x^2 - 3 * x"

#define FMPQ_POLY_ZERO_INIT(x) fmpq_poly_init(x)
#define FMPQ_POLY_ZERO_CLEAR(x) fmpq_poly_clear(x)
#define FMPQ_POLY_ZERO_SET(x) do { } while (0)
#define FMPQ_POLY_ZERO_STRING "0"

#define FMPQ_POLY_CONSTANT_COEFF "-1337/10937"
#define FMPQ_POLY_CONSTANT_INIT(x) fmpq_poly_init2(x, 1)
#define FMPQ_POLY_CONSTANT_CLEAR(x) fmpq_poly_clear(x)
#define FMPQ_POLY_CONSTANT_SET(x) fmpq_poly_set_str(x, "1  " FMPQ_POLY_CONSTANT_COEFF)
#define FMPQ_POLY_CONSTANT_STRING "-1337 / 10937"

#define FMPQ_POLY1_LENGTH 4
#define FMPQ_POLY1_INIT(x) fmpq_poly_init2(x, FMPQ_POLY1_LENGTH)
#define FMPQ_POLY1_CLEAR(x) fmpq_poly_clear(x)
#define FMPQ_POLY1_SET(x) \
    fmpq_poly_set_str(x, STR(FMPQ_POLY1_LENGTH) "  " \
            "-3/10 -3 0 -1")
#define FMPQ_POLY1_STRING "-x^3 - 3 * x - 3 / 10"

#define FMPQ_POLY2_LENGTH 4
#define FMPQ_POLY2_INIT(x) fmpq_poly_init2(x, FMPQ_POLY2_LENGTH)
#define FMPQ_POLY2_CLEAR(x) fmpq_poly_clear(x)
#define FMPQ_POLY2_SET(x) \
    fmpq_poly_set_str(x, STR(FMPQ_POLY2_LENGTH) "  " \
            "0 10/3 -3 8")
#define FMPQ_POLY2_STRING "8 * x^3 - 3 * x^2 + 10 / 3 * x"

#define ARB_POLY_LENGTH 5
#define ARB_POLY_PREC 30
#define ARB_POLY_INIT(x) arb_poly_init2(x, ARB_POLY_LENGTH)
#define ARB_POLY_CLEAR(x) arb_poly_clear(x)
#define ARB_POLY_SET(x) \
do                      \
{                       \
    arb_poly_set_coeff_si(x, ARB_POLY_LENGTH - 1, 1); \
    arb_poly_set_coeff_si(x, 3, -1); \
    arb_set_str((x)->coeffs + 2, "12.1 +/- 1", ARB_POLY_PREC); \
    arb_set_str((x)->coeffs + 1, "-12.1 +/- 1", ARB_POLY_PREC); \
    arb_set_str((x)->coeffs + 0, "-12.1 +/- 1", ARB_POLY_PREC); \
} while (0)
#define ARB_POLY_STRING "x^4 - x^3 + [1e+1 +/- 3.11] * x^2 - [1e+1 +/- 3.11] * x - [1e+1 +/- 3.11]"

#define ACB_POLY_LENGTH 4
#define ACB_POLY_PREC 30
#define ACB_POLY_INIT(x) acb_poly_init2(x, ACB_POLY_LENGTH)
#define ACB_POLY_CLEAR(x) acb_poly_clear(x)
#define ACB_POLY_SET(x) \
do                      \
{                       \
    _acb_poly_set_length(x, ACB_POLY_LENGTH); \
    arb_set_str(acb_realref((x)->coeffs + ACB_POLY_LENGTH - 1), "-12.1 +/- 1", ACB_POLY_PREC); \
    arb_set_si(acb_imagref((x)->coeffs + ACB_POLY_LENGTH - 1), -1); \
    acb_poly_set_coeff_si(x, 2, 1); \
    acb_poly_set_coeff_si(x, 1, -1); \
    acb_poly_set_coeff_si(x, 0, 0); \
} while (0)
#define ACB_POLY_STRING "([-1e+1 +/- 3.11] - i) * x^3 + x^2 - x"

static void test_composite_string(FILE * fs)
{
    int res1, res2;
    char * str1, * str2;

    ulong xulong1 = UWORD(93112), xulong2 = UWORD(8721);
    slong xslong = WORD(-1982);
    double xdouble = 0.123456;
    char xchar = 'p';
    short xshort = 1872;
    int xint = -13214;
    size_t xsize = 41121;
    char xcharp[] = "julafton";
    wint_t xwint = L'a';
    long int xlong = -1872381273;
    long long int xlonglong = LLONG_MIN;
    intmax_t xintmax = INTMAX_MAX;
    ptrdiff_t xptrdiff = 10241;
    long double xlongdouble = 123.12398128738172381287L;
    void * xpointer = (void *) 123;
    wchar_t xwcharp[] = L"a b c";

    nmod_t xnmod;
    fmpz_t xfmpz1, xfmpz2;
    fmpz_mod_ctx_t xfmpz_mod_ctx;
    fmpq_t xfmpq1, xfmpq2;
    arf_t xarf1, xarf2, xarf3;
    mag_t xmag1, xmag2, xmag3;
    arb_t xarb1, xarb2, xarb3;
    acb_t xacb1, xacb2, xacb3, xacb4, xacb5, xacb6, xacb7, xacb8;

    slong xslong_vec[SLONG_VEC_LEN];
    mp_ptr xnmod_vec;
    fmpz * xfmpz_vec;
    fmpq * xfmpq_vec;
    arb_ptr xarb_vec;
    acb_ptr xacb_vec;

    /* Matrix printing relies on vector printing, so no need to check other
     * types */
    fmpz_mat_t empty_matrix;
    nmod_mat_t xnmod_mat; nmod_mat_t xnmod_mat_window;
    fmpz_mat_t xfmpz_mat;

    /* NOTE: We need extra checks with fmpq_poly as it is treated differntly in
     * __flint_poly_fprint. */
    nmod_poly_t xnmod_poly_zero, xnmod_poly_constant, xnmod_poly;
    fmpz_poly_t xfmpz_poly;
    fmpq_poly_t xfmpq_poly_zero, xfmpq_poly_constant, xfmpq_poly1, xfmpq_poly2;
    arb_poly_t xarb_poly;
    acb_poly_t xacb_poly;

    /* Initialize ************************************************************/
    NMOD_INIT(xnmod);
    FMPZ1_INIT(xfmpz1);
    FMPZ2_INIT(xfmpz2);
    FMPZ_MOD_CTX_INIT(xfmpz_mod_ctx);
    FMPQ1_INIT(xfmpq1);
    FMPQ2_INIT(xfmpq2);
    ARF1_INIT(xarf1);
    ARF2_INIT(xarf2);
    ARF3_INIT(xarf3);
    MAG1_INIT(xmag1);
    MAG2_INIT(xmag2);
    MAG3_INIT(xmag3);
    ARB1_INIT(xarb1);
    ARB2_INIT(xarb2);
    ARB3_INIT(xarb3);
    ACB1_INIT(xacb1);
    ACB2_INIT(xacb2);
    ACB3_INIT(xacb3);
    ACB4_INIT(xacb4);
    ACB5_INIT(xacb5);
    ACB6_INIT(xacb6);
    ACB7_INIT(xacb7);
    ACB8_INIT(xacb8);

    SLONG_VEC_INIT(xslong_vec);
    NMOD_VEC_INIT(xnmod_vec);
    FMPZ_VEC_INIT(xfmpz_vec);
    FMPQ_VEC_INIT(xfmpq_vec);
    ARB_VEC_INIT(xarb_vec);
    ACB_VEC_INIT(xacb_vec);

    NMOD_MAT_WITH_WINDOW_INIT(xnmod_mat, xnmod_mat_window);
    FMPZ_MAT_INIT(xfmpz_mat);
    FMPZ_MAT_EMPTY_INIT(empty_matrix);

    NMOD_POLY_ZERO_INIT(xnmod_poly_zero);
    NMOD_POLY_CONSTANT_INIT(xnmod_poly_constant);
    NMOD_POLY_INIT(xnmod_poly);
    FMPZ_POLY_INIT(xfmpz_poly);
    FMPQ_POLY_ZERO_INIT(xfmpq_poly_zero);
    FMPQ_POLY_CONSTANT_INIT(xfmpq_poly_constant);
    FMPQ_POLY1_INIT(xfmpq_poly1);
    FMPQ_POLY2_INIT(xfmpq_poly2);
    ARB_POLY_INIT(xarb_poly);
    ACB_POLY_INIT(xacb_poly);

    /* Set *******************************************************************/
    NMOD_SET(xnmod);
    FMPZ1_SET(xfmpz1);
    FMPZ2_SET(xfmpz2);
    FMPZ_MOD_CTX_SET(xfmpz_mod_ctx);
    FMPQ1_SET(xfmpq1);
    FMPQ2_SET(xfmpq2);
    ARF1_SET(xarf1);
    ARF2_SET(xarf2);
    ARF3_SET(xarf3);
    MAG1_SET(xmag1);
    MAG2_SET(xmag2);
    MAG3_SET(xmag3);
    ARB1_SET(xarb1);
    ARB2_SET(xarb2);
    ARB3_SET(xarb3);
    ACB1_SET(xacb1);
    ACB2_SET(xacb2);
    ACB3_SET(xacb3);
    ACB4_SET(xacb4);
    ACB5_SET(xacb5);
    ACB6_SET(xacb6);
    ACB7_SET(xacb7);
    ACB8_SET(xacb8);

    SLONG_VEC_SET(xslong_vec);
    NMOD_VEC_SET(xnmod_vec);
    FMPZ_VEC_SET(xfmpz_vec);
    FMPQ_VEC_SET(xfmpq_vec);
    ARB_VEC_SET(xarb_vec);
    ACB_VEC_SET(xacb_vec);

    NMOD_MAT_WITH_WINDOW_SET(xnmod_mat, xnmod_mat_window);
    FMPZ_MAT_SET(xfmpz_mat);
    FMPZ_MAT_EMPTY_SET(empty_matrix);

    NMOD_POLY_ZERO_SET(xnmod_poly_zero);
    NMOD_POLY_CONSTANT_SET(xnmod_poly_constant);
    NMOD_POLY_SET(xnmod_poly);
    FMPZ_POLY_SET(xfmpz_poly);
    FMPQ_POLY_ZERO_SET(xfmpq_poly_zero);
    FMPQ_POLY_CONSTANT_SET(xfmpq_poly_constant);
    FMPQ_POLY1_SET(xfmpq_poly1);
    FMPQ_POLY2_SET(xfmpq_poly2);
    ARB_POLY_SET(xarb_poly);
    ACB_POLY_SET(xacb_poly);

    /* Print *****************************************************************/
#define STR_SIZE 10000 /* 10 kB should suffice. */
    str1 = flint_calloc(STR_SIZE, sizeof(char));
    str2 = flint_calloc(STR_SIZE, sizeof(char));

#if defined(_LONG_LONG_LIMB)
# define ULONG_SLONG_STR \
            "Here we print a ulong: %020llu\n" \
            "Here we print a slong: %20lld\n"
#else
# define ULONG_SLONG_STR \
            "Here we print a ulong: %020lu\n" \
            "Here we print a slong: %20ld\n"
#endif

    res1 = snprintf(str1, STR_SIZE,
            "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do "
            "eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut "
            "enim ad minim veniam, quis nostrud exercitation ullamco laboris "
            "nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor "
            "in reprehenderit in voluptate velit esse cillum dolore eu fugiat "
            "nulla pariatur. Excepteur sint occaecat cupidatat non proident, "
            "sunt in culpa qui officia deserunt mollit anim id est laborum.\n"
            "\n"
            "Here are some %% characters: %%, %% and %%%%\n"
            "\n"
            "We will now start printing primitive types...\n"
            "\n"
            ULONG_SLONG_STR
            "Here we print a ulong in hexadecimal: " WORD_FMT "x\n"
            "Here we print a double with flags: %05.2f\n"
            "Here we print another double: %le\n"
            "\n"
            "We will now start printing FLINT types...\n"
            "\n"
            "ulong: " WORD_FMT "u\n"
            "slong: " WORD_FMT "d\n"
            "nmod: " NMOD_STRING "\n"
            "small fmpz: " FMPZ1_STRING "\n"
            "big fmpz: " FMPZ2_STRING "\n"
            "fmpz_mod_ctx: " FMPZ_MOD_CTX_STRING "\n"
            "integer fmpq: " FMPQ1_STRING "\n"
            "fmpq: " FMPQ2_STRING "\n"
            "zero arf: " ARF1_STRING "\n"
            "integer arf: " ARF2_STRING "\n"
            "arf: " ARF3_STRING "\n"
            "zero mag: " MAG1_STRING "\n"
            "integer mag: " MAG2_STRING "\n"
            "mag: " MAG3_STRING "\n"
            "zero arb: " ARB1_STRING "\n"
            "integer arb: " ARB2_STRING "\n"
            "arb: " ARB3_STRING "\n"
            "zero acb: " ACB1_STRING "\n"
            "half gaussian acb: " ACB2_STRING "\n"
            "other half gaussian acb: " ACB3_STRING "\n"
            "another half gaussian acb: " ACB4_STRING "\n"
            "gaussian integer acb: " ACB5_STRING "\n"
            "real acb: " ACB6_STRING "\n"
            "imaginary acb: " ACB7_STRING "\n"
            "acb: " ACB8_STRING "\n"
            "\n"
            "We intersect with some other primitive types...\n"
            "\n"
            "char: % 20hhd\n"
            "short: %-20hx\n"
            "int: %.*i\n"
            "size_t: %zo\n"
            "char *: %.20s\n"
            "wint_t: %lc\n"
            "long int: %.*lX\n"
            "long long int: %3llu\n"
            "intmax_t: %ji\n"
            "ptrdiff_t: %tu\n"
            "long double: %+33.12LG\n"
            "pointer: %p\n"
            "wchar_t *: %ls\n"
            "\n"
            "And now we go back to printing FLINT types...\n"
            "\n"
            "slong_vec: " SLONG_VEC_STRING "\n"
            "nmod_vec: " NMOD_VEC_STRING "\n"
            "fmpz_vec: " FMPZ_VEC_STRING "\n"
            "fmpq_vec: " FMPQ_VEC_STRING "\n"
            "arb_vec: " ARB_VEC_STRING "\n"
            "acb_vec: " ACB_VEC_STRING "\n"
            "\n"
            "empty fmpz_mat: " FMPZ_MAT_EMTPY_STRING "\n"
            "window nmod_mat: " NMOD_MAT_WITH_WINDOW_STRING "\n"
            "fmpz_mat: " FMPZ_MAT_STRING "\n"
            "\n"
            "zero nmod_poly: " NMOD_POLY_ZERO_STRING "\n"
            "constant nmod_poly: " NMOD_POLY_CONSTANT_STRING "\n"
            "nmod_poly: " NMOD_POLY_STRING "\n"
            "fmpz_poly: " FMPZ_POLY_STRING "\n"
            "zero fmpq_poly: " FMPQ_POLY_ZERO_STRING "\n"
            "constant fmpq_poly: " FMPQ_POLY_CONSTANT_STRING "\n"
            "fmpq_poly (1): " FMPQ_POLY1_STRING "\n"
            "fmpq_poly (2): " FMPQ_POLY2_STRING "\n"
            "arb_poly: " ARB_POLY_STRING "\n"
            "acb_poly: " ACB_POLY_STRING "\n",
            xulong1,
            xslong,
            xulong2,
            xdouble,
            xdouble,
            xulong1,
            xslong,
            xchar,
            xshort,
            10, xint,
            xsize,
            xcharp,
            xwint,
            4, xlong,
            xlonglong,
            xintmax,
            xptrdiff,
            xlongdouble,
            xpointer,
            xwcharp);

#undef ULONG_SLONG_STR

    if (res1 == STR_SIZE - 1)
        flint_throw(FLINT_TEST_FAIL,
                "Could not print expected string into str1.\n");

    if (res1 < 0)
        flint_throw(FLINT_TEST_FAIL,
                "Negative return value from snprintf.\n");

    res2 = flint_fprintf(fs,
            "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do "
            "eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut "
            "enim ad minim veniam, quis nostrud exercitation ullamco laboris "
            "nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor "
            "in reprehenderit in voluptate velit esse cillum dolore eu fugiat "
            "nulla pariatur. Excepteur sint occaecat cupidatat non proident, "
            "sunt in culpa qui officia deserunt mollit anim id est laborum.\n"
            "\n"
            "Here are some %% characters: %%, %% and %%%%\n"
            "\n"
            "We will now start printing primitive types...\n"
            "\n"
            "Here we print a ulong: %020wu\n"
            "Here we print a slong: %20wd\n"
            "Here we print a ulong in hexadecimal: %wx\n"
            "Here we print a double with flags: %05.2f\n"
            "Here we print another double: %le\n"
            "\n"
            "We will now start printing FLINT types...\n"
            "\n"
            "ulong: %{ulong}\n"
            "slong: %{slong}\n"
            "nmod: %{nmod}\n"
            "small fmpz: %{fmpz}\n"
            "big fmpz: %{fmpz}\n"
            "fmpz_mod_ctx: %{fmpz_mod_ctx}\n"
            "integer fmpq: %{fmpq}\n"
            "fmpq: %{fmpq}\n"
            "zero arf: %{arf}\n"
            "integer arf: %{arf}\n"
            "arf: %{arf}\n"
            "zero mag: %{mag}\n"
            "integer mag: %{mag}\n"
            "mag: %{mag}\n"
            "zero arb: %{arb}\n"
            "integer arb: %{arb}\n"
            "arb: %{arb}\n"
            "zero acb: %{acb}\n"
            "half gaussian acb: %{acb}\n"
            "other half gaussian acb: %{acb}\n"
            "another half gaussian acb: %{acb}\n"
            "gaussian integer acb: %{acb}\n"
            "real acb: %{acb}\n"
            "imaginary acb: %{acb}\n"
            "acb: %{acb}\n"
            "\n"
            "We intersect with some other primitive types...\n"
            "\n"
            "char: % 20hhd\n"
            "short: %-20hx\n"
            "int: %.*i\n"
            "size_t: %zo\n"
            "char *: %.20s\n"
            "wint_t: %lc\n"
            "long int: %.*lX\n"
            "long long int: %3llu\n"
            "intmax_t: %ji\n"
            "ptrdiff_t: %tu\n"
            "long double: %+33.12LG\n"
            "pointer: %p\n"
            "wchar_t *: %ls\n"
            "\n"
            "And now we go back to printing FLINT types...\n"
            "\n"
            "slong_vec: %{slong*}\n"
            "nmod_vec: %{ulong*}\n"
            "fmpz_vec: %{fmpz*}\n"
            "fmpq_vec: %{fmpq*}\n"
            "arb_vec: %{arb*}\n"
            "acb_vec: %{acb*}\n"
            "\n"
            "empty fmpz_mat: %{fmpz_mat}\n"
            "window nmod_mat: %{nmod_mat}\n"
            "fmpz_mat: %{fmpz_mat}\n"
            "\n"
            "zero nmod_poly: %{nmod_poly}\n"
            "constant nmod_poly: %{nmod_poly}\n"
            "nmod_poly: %{nmod_poly}\n"
            "fmpz_poly: %{fmpz_poly}\n"
            "zero fmpq_poly: %{fmpq_poly}\n"
            "constant fmpq_poly: %{fmpq_poly}\n"
            "fmpq_poly (1): %{fmpq_poly}\n"
            "fmpq_poly (2): %{fmpq_poly}\n"
            "arb_poly: %{arb_poly}\n"
            "acb_poly: %{acb_poly}\n",
            xulong1,
            xslong,
            xulong2,
            xdouble,
            xdouble,
            xulong1,
            xslong,
            xnmod,
            xfmpz1,
            xfmpz2,
            xfmpz_mod_ctx,
            xfmpq1,
            xfmpq2,
            xarf1,
            xarf2,
            xarf3,
            xmag1,
            xmag2,
            xmag3,
            xarb1,
            xarb2,
            xarb3,
            xacb1,
            xacb2,
            xacb3,
            xacb4,
            xacb5,
            xacb6,
            xacb7,
            xacb8,
            xchar,
            xshort,
            10, xint,
            xsize,
            xcharp,
            xwint,
            4, xlong,
            xlonglong,
            xintmax,
            xptrdiff,
            xlongdouble,
            xpointer,
            xwcharp,
            xslong_vec, SLONG_VEC_LEN,
            xnmod_vec, NMOD_VEC_LEN,
            xfmpz_vec, FMPZ_VEC_LEN,
            xfmpq_vec, FMPQ_VEC_LEN,
            xarb_vec, ARB_VEC_LEN,
            xacb_vec, ACB_VEC_LEN,
            empty_matrix,
            xnmod_mat_window,
            xfmpz_mat,
            xnmod_poly_zero,
            xnmod_poly_constant,
            xnmod_poly,
            xfmpz_poly,
            xfmpq_poly_zero,
            xfmpq_poly_constant,
            xfmpq_poly1,
            xfmpq_poly2,
            xarb_poly,
            xacb_poly);

    if (res2 > STR_SIZE - 1)
        flint_throw(FLINT_TEST_FAIL,
                "Printed more than expected into fs.\n");

    if (res2 < 0)
        flint_throw(FLINT_TEST_FAIL,
                "Negative return value from flint_fprintf.\n");

#undef STR_SIZE

    /* Check *****************************************************************/
    fseek(fs, 0, SEEK_SET);
    if (fread(str2, sizeof(char), res2, fs) != res2)
        flint_throw(FLINT_TEST_FAIL,
                "Could not read %d bytes from filestream.\n",
                res2);

    /* Check strings line by line */
    check_strings(str1, str2);

    /* Check return value */
    if (res1 != res2)
        flint_throw(FLINT_TEST_FAIL,
                "Result from flint_fprintf differed.\n"
                "Expected: %d\n"
                "Got:      %d\n",
                res1, res2);

    /* Clear *****************************************************************/
    flint_free(str1);
    flint_free(str2);

    NMOD_CLEAR(xnmod);
    FMPZ1_CLEAR(xfmpz1);
    FMPZ2_CLEAR(xfmpz2);
    FMPZ_MOD_CTX_CLEAR(xfmpz_mod_ctx);
    FMPQ1_CLEAR(xfmpq1);
    FMPQ2_CLEAR(xfmpq2);
    ARF1_CLEAR(xarf1);
    ARF2_CLEAR(xarf2);
    ARF3_CLEAR(xarf3);
    MAG1_CLEAR(xmag1);
    MAG2_CLEAR(xmag2);
    MAG3_CLEAR(xmag3);
    ARB1_CLEAR(xarb1);
    ARB2_CLEAR(xarb2);
    ARB3_CLEAR(xarb3);
    ACB1_CLEAR(xacb1);
    ACB2_CLEAR(xacb2);
    ACB3_CLEAR(xacb3);
    ACB4_CLEAR(xacb4);
    ACB5_CLEAR(xacb5);
    ACB6_CLEAR(xacb6);
    ACB7_CLEAR(xacb7);
    ACB8_CLEAR(xacb8);

    SLONG_VEC_CLEAR(xslong_vec);
    NMOD_VEC_CLEAR(xnmod_vec);
    FMPZ_VEC_CLEAR(xfmpz_vec);
    FMPQ_VEC_CLEAR(xfmpq_vec);
    ARB_VEC_CLEAR(xarb_vec);
    ACB_VEC_CLEAR(xacb_vec);

    NMOD_MAT_WITH_WINDOW_CLEAR(xnmod_mat, xnmod_mat_window);
    FMPZ_MAT_CLEAR(xfmpz_mat);
    FMPZ_MAT_EMPTY_CLEAR(empty_matrix);

    NMOD_POLY_ZERO_CLEAR(xnmod_poly_zero);
    NMOD_POLY_CONSTANT_CLEAR(xnmod_poly_constant);
    NMOD_POLY_CLEAR(xnmod_poly);
    FMPZ_POLY_CLEAR(xfmpz_poly);
    FMPQ_POLY_ZERO_CLEAR(xfmpq_poly_zero);
    FMPQ_POLY_CONSTANT_CLEAR(xfmpq_poly_constant);
    FMPQ_POLY1_CLEAR(xfmpq_poly1);
    FMPQ_POLY2_CLEAR(xfmpq_poly2);
    ARB_POLY_CLEAR(xarb_poly);
    ACB_POLY_CLEAR(xacb_poly);
}

#if 0
# define test_invalid_printing test_invalid_printing
static void test_invalid_printing(FILE * fs)
{
    int res1, res2;
    size_t nflag1[6], nflag2[6];
    char * str1, * str2;

#define STR_SIZE 1000
    str1 = flint_calloc(STR_SIZE, sizeof(char));
    str2 = flint_calloc(STR_SIZE, sizeof(char));

#ifdef __GNUC__
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wformat"
#endif
    res1 = snprintf(str1, STR_SIZE,
            "We will now print some invalid strings:\n"
            "\n"
            "%{tjenare\n%zn"
            "%{tjenixen}\n"
            "%{lagom}\n%zn"
            "%.12 lld\n"
            "%.12 ld\n%zn"
            "%.12 wd\n%zn"
            "%\n%zn"
            "%{fmpz_this_is_not_valid}\n%zn",
            nflag1 + 0, nflag1 + 1, nflag1 + 2,
            nflag1 + 3, nflag1 + 4, nflag1 + 5);

    if (res1 == STR_SIZE - 1)
        flint_throw(FLINT_TEST_FAIL,
                "Could not print expected string into str1.\n");

    res2 = flint_fprintf(fs,
            "We will now print some invalid strings:\n"
            "\n"
            "%{tjenare\n%zn"
            "%{tjenixen}\n"
            "%{lagom}\n%zn"
            "%.12 lld\n"
            "%.12 ld\n%zn"
            "%.12 wd\n%zn"
            "%\n%zn"
            "%{fmpz_this_is_not_valid}\n%zn",
            nflag2 + 0, nflag2 + 1, nflag2 + 2,
            nflag2 + 3, nflag2 + 4, nflag2 + 5);

    if (res2 > STR_SIZE - 1)
        flint_throw(FLINT_TEST_FAIL,
                "Printed more than expected into fs.\n");

#ifdef __GNUC__
# pragma GCC diagnostic pop
#endif

#undef STR_SIZE

    if (res1 < 0 && res2 < 0)
    {
        /* Do nothing */
    }
    else if (res1 < 0 || res2 < 0)
        flint_throw(FLINT_TEST_FAIL,
                "Return value for %s was negative while return value for %s was not.\n",
                (res1 < 0) ? "sprintf" : "flint_fprintf",
                (res2 < 0) ? "sprintf" : "flint_fprintf");
    else
    {
        slong ix;

        fseek(fs, 0, SEEK_SET);
        if (fread(str2, sizeof(char), res2, fs) != res2)
            flint_throw(FLINT_TEST_FAIL,
                    "Could not read %d bytes from filestream.\n",
                    res2);

        check_strings(str1, str2);

        if (res1 != res2)
            flint_throw(FLINT_TEST_FAIL,
                    "Result from flint_fprintf differed.\n"
                    "Expected: %d\n"
                    "Got:      %d\n",
                    res1, res2);

        /* Check that nflags work as they should */
        for (ix = 0; ix < sizeof(nflag1) / sizeof(size_t); ix++)
            if (nflag1[ix] != nflag2[ix])
                flint_throw(FLINT_TEST_FAIL,
                        "The " WORD_FMT "d'th nflag differed.\n"
                        "Expected: %d\n"
                        "Got:      %d\n",
                        ix, nflag1[ix], nflag2[ix]);
    }

    flint_free(str1);
    flint_free(str2);
}
#else
# define test_invalid_printing(fs)
#endif

#define TMP_FILENAME "tmp"

TEST_FUNCTION_START(flint_fprintf, state)
{
    FILE * fs;

    fs = fopen(TMP_FILENAME, "w+");
    if (fs == NULL)
        flint_throw(FLINT_TEST_FAIL, "Could not open temporary file \"" TMP_FILENAME "\"\n");

    test_composite_string(fs);

    fs = freopen(TMP_FILENAME, "w+", fs);
    if (fs == NULL)
        flint_throw(FLINT_TEST_FAIL, "Could not reopen temporary file \"" TMP_FILENAME "\"\n");

    test_invalid_printing(fs);

    if (fclose(fs))
        flint_throw(FLINT_TEST_FAIL, "Could not close temporary file \"" TMP_FILENAME "\"\n");

    if (remove(TMP_FILENAME))
        flint_throw(FLINT_TEST_FAIL, "Could not remove temporary file \"" TMP_FILENAME "\"\n");

    TEST_FUNCTION_END(state);
}

#if HAVE_UNISTD_H && FLINT_COVERAGE
#include <unistd.h>
TEST_FUNCTION_START(flint_printf, state)
{
    FILE * fs;
    ulong xulong = 9812;
    slong xslong = -123;
    char str[128];
    int res;
    int original_stdout_fd = dup(fileno(stdout));

    fflush(stdout);
    fs = freopen(TMP_FILENAME, "w", stdout);
    if (fs == NULL)
        flint_throw(FLINT_TEST_FAIL, "Could not redirect stdout to \"" TMP_FILENAME "\"\n");

    res = flint_printf(
            "Tjingeling!\n"
            "ulong = %wu\n"
            "slong = %wd\n",
            xulong, xslong);
    fputc('\0', stdout);
    fflush(stdout);

    if (res < 0)
        flint_throw(FLINT_TEST_FAIL, "res = %d\n", res);

    if (dup2(original_stdout_fd, fileno(stdout)) == -1)
        flint_throw(FLINT_TEST_FAIL, "Could not restore stdout\n");

    fs = fopen(TMP_FILENAME, "r");
    if (fs == NULL)
        flint_throw(FLINT_TEST_FAIL, "Could not open \"" TMP_FILENAME "\"\n");

    res = fread(str, sizeof(char), res + 1, fs);

    if (strcmp(str, "Tjingeling!\nulong = 9812\nslong = -123\n"))
        flint_throw(FLINT_TEST_FAIL, "Strings not equal.\n");

    if (fclose(fs))
        flint_throw(FLINT_TEST_FAIL, "Could not close temporary file \"" TMP_FILENAME "\"\n");

    if (remove(TMP_FILENAME))
        flint_throw(FLINT_TEST_FAIL, "Could not remove temporary file \"" TMP_FILENAME "\"\n");

    TEST_FUNCTION_END(state);
}
#else
TEST_FUNCTION_START(flint_printf, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif

#undef TMP_FILENAME
