/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_factor.h"

#define TMP_FILENAME "tmp"
#define STRLEN(x) (sizeof(x) - sizeof(char))
#define RESULT_STRING   \
    "0\n"               \
    "1\n"               \
    "-1\n"              \
    "3 * 5 * 7\n"       \
    "-1 * 3^2 * 5^3 * 7"

TEST_FUNCTION_START(fmpz_factor_fprint, state)
{
    FILE * fs;
    fmpz_t a;
    fmpz_factor_t f;
    int res, tmp;
    char str[sizeof(RESULT_STRING)];

    fs = fopen(TMP_FILENAME, "w+");
    if (fs == NULL)
        flint_throw(FLINT_TEST_FAIL, "Could not open temporary file \"" TMP_FILENAME "\"\n");

    fmpz_init(a);
    fmpz_factor_init(f);
    res = 0;

    /* Print 0 */
    fmpz_zero(a);
    fmpz_factor(f, a);
    tmp = fmpz_factor_fprint(fs, f);
    if (tmp != STRLEN("0"))
        goto wrong_return_value;
    res += tmp;
    res += (fputc('\n', fs) != EOF);

    /* Print 1 */
    fmpz_one(a);
    fmpz_factor(f, a);
    tmp = fmpz_factor_fprint(fs, f);
    if (tmp != STRLEN("1"))
        goto wrong_return_value;
    res += tmp;
    res += (fputc('\n', fs) != EOF);

    /* Print -1 */
    fmpz_set_si(a, -1);
    fmpz_factor(f, a);
    tmp = fmpz_factor_fprint(fs, f);
    if (tmp != STRLEN("-1"))
        goto wrong_return_value;
    res += tmp;
    res += (fputc('\n', fs) != EOF);

    /* Print 3 * 5 * 7 */
    fmpz_set_si(a, 3 * 5 * 7);
    fmpz_factor(f, a);
    tmp = fmpz_factor_fprint(fs, f);
    if (tmp != STRLEN("3 * 5 * 7"))
        goto wrong_return_value;
    res += tmp;
    res += (fputc('\n', fs) != EOF);

    /* Print -1 * 3^2 * 5^3 * 7 */
    fmpz_set_si(a, -1 * 9 * 125 * 7);
    fmpz_factor(f, a);
    tmp = fmpz_factor_fprint(fs, f);
    if (tmp != STRLEN("-1 * 3^2 * 5^3 * 7"))
    {
wrong_return_value:
        fseek(fs, -tmp, SEEK_CUR);
        fread(str, sizeof(char), tmp, fs);
        str[tmp] = '\0';
        flint_throw(FLINT_TEST_FAIL,
                "Wrong return value for a = %{fmpz}\n"
                "Got string: %s\n",
                a, str);
    }
    res += tmp;
    res += (fputc('\0', fs) != EOF);

    fseek(fs, 0, SEEK_SET);
    if (fread(str, sizeof(char), res, fs) != res)
        flint_throw(FLINT_TEST_FAIL, "Could not read %d bytes from filestream.\n", res);

    if (strcmp(str, RESULT_STRING))
        flint_throw(FLINT_TEST_FAIL,
                "Result differed from expected output.\n"
                "Got:\n"
                "%s\n\n"
                "Expected:\n"
                "%s\n",
                str, RESULT_STRING);

    fmpz_clear(a);
    fmpz_factor_clear(f);

    if (fclose(fs))
        flint_throw(FLINT_TEST_FAIL, "Could not close temporary file \"" TMP_FILENAME "\"\n");

    if (remove(TMP_FILENAME))
        flint_throw(FLINT_TEST_FAIL, "Could not remove temporary file \"" TMP_FILENAME "\"\n");

    TEST_FUNCTION_END(state);
}

#undef TMP_FILENAME
#undef STRLEN
#undef RESULT_STRING
