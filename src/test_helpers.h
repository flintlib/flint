/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2023, 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include <time.h> /* clock */
#include <string.h> /* strlen, strcmp, strncmp */
#include <stdlib.h> /* strtol */
#include <limits.h> /* LONG_MIN */
#include "templates.h"
#include "flint.h"

#if FLINT_WANT_PRETTY_TESTS
# define _RED "\033[0;31m"
# define _RED_B "\033[1;31m"
# define _GREEN "\033[0;32m"
# define _GREEN_B "\033[1;32m"
# define _YELLOW "\033[0;33m"
# define _YELLOW_B "\033[1;33m"
# define _RESET "\033[0m"
#else
# define _RED
# define _RED_B
# define _GREEN
# define _GREEN_B
# define _YELLOW
# define _YELLOW_B
# define _RESET
#endif

/* easy way to test a condition in test code */
#define FLINT_TEST(e)                                   \
    do {                                                \
        if (!(e))                                       \
        {                                               \
            flint_throw(FLINT_ERROR, "test %s FAILED\n%s:%d", #e, __FILE__, __LINE__); \
        }                                               \
    } while (0)

/* test function macro *******************************************************/

typedef struct
{
    int (* test_function)(void);
    char * name;
}
test_struct;

#define TEST_FUNCTION(label) { CAT(test, label), #label }

#define TEST_FUNCTION_START(label, state)               \
int CAT(test, label)(void)                              \
{                                                       \
    double FLINT_SET_BUT_UNUSED(_start_time_), FLINT_SET_BUT_UNUSED(_end_time_); \
    char _test_io_string_[128] = #label "                                                        "; \
    const int _label_len_ = sizeof(#label) - 1;         \
    puts(#label "...");                                 \
    FLINT_TEST_INIT(state);                             \
    _start_time_ = clock();

#define TEST_GR_FUNCTION_START(label, state, count_success, count_domain, count_unable) \
int CAT(test, label)(void)                              \
{                                                       \
    slong count_success = 0, count_unable = 0, count_domain = 0; \
    double FLINT_SET_BUT_UNUSED(_start_time_), FLINT_SET_BUT_UNUSED(_end_time_); \
    char _test_io_string_[128] = #label "                                                        "; \
    const int _label_len_ = sizeof(#label) - 1;         \
    puts(#label "...");                                 \
    FLINT_TEST_INIT(state);                             \
    _start_time_ = clock();

#define TEST_TEMPLATE_FUNCTION_START(T, label, state)   \
int TEMPLATE3(test, T, label)(void)                     \
{                                                       \
    double FLINT_SET_BUT_UNUSED(_start_time_), FLINT_SET_BUT_UNUSED(_end_time_); \
    char _test_io_string_[128] = TEMPLATE_STR(T) "_" #label "                                                      "; \
    const int _label_len_ = sizeof(TEMPLATE_STR(T) "_" #label) - 1; \
    puts(TEMPLATE_STR(T) "_" #label "...");             \
    FLINT_TEST_INIT(state);                             \
    _start_time_ = clock();

#define TEST_TEMPLATE2_FUNCTION_START(T, label1, label2, state)\
int TEMPLATE5(test, T, label1, T, label2)(void)         \
{                                                       \
    double FLINT_SET_BUT_UNUSED(_start_time_), FLINT_SET_BUT_UNUSED(_end_time_); \
    char _test_io_string_[128] = TEMPLATE_STR(T) "_" #label1 "_" TEMPLATE_STR(T) "_" #label2 "                                                  "; \
    const int _label_len_ = sizeof(TEMPLATE_STR(T) "_" #label1 "_" TEMPLATE_STR(T) "_" #label2) - 1; \
    puts(TEMPLATE_STR(T) "_" #label1 "_" TEMPLATE_STR(T) "_" #label2 "..."); \
    FLINT_TEST_INIT(state);                             \
    _start_time_ = clock();

#define TEST_FUNCTION_END(state)                        \
    _end_time_ = clock();                               \
    FLINT_TEST_CLEAR(state);                            \
    if (_label_len_ < 48)                               \
        printf("%.48s%6.2f   (" _GREEN_B "PASS" _RESET ")\n", \
                _test_io_string_,                       \
                (_end_time_ - _start_time_) / CLOCKS_PER_SEC); \
    else                                                \
        printf("%.*s\n%54.2f   (" _GREEN_B "PASS" _RESET ")\n", \
                _label_len_, _test_io_string_,          \
                (_end_time_ - _start_time_) / CLOCKS_PER_SEC); \
    return 0;                                           \
}

#define TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable) \
    _end_time_ = clock();                               \
    FLINT_TEST_CLEAR(state);                            \
    printf("%.*s\n  "                                   \
            "%5" _WORD_FMT "d success, "                \
            "%5" _WORD_FMT "d domain, "                 \
            "%5" _WORD_FMT "d unable"                   \
            "%11.2f   "                                 \
            "(" _GREEN_B "PASS" _RESET ")\n",           \
            _label_len_, _test_io_string_,              \
            count_success, count_domain, count_unable,  \
            (_end_time_ - _start_time_) / CLOCKS_PER_SEC); \
    return 0;                                           \
}

#define TEST_FUNCTION_END_SKIPPED(state)                \
    FLINT_TEST_CLEAR(state);                            \
    if (_label_len_ < 54)                               \
        printf("%.*s(" _YELLOW_B "SKIPPED" _RESET ")\n", \
                54, _test_io_string_);                  \
    else                                                \
        printf("%.*s\n%*s(" _YELLOW_B "SKIPPED" _RESET ")\n", \
                _label_len_, _test_io_string_,          \
                54, "");                                \
    return 0;                                           \
}

#define TEST_FUNCTION_FAIL(...)                         \
do                                                      \
{                                                       \
    _end_time_ = clock();                               \
    if (_label_len_ < 48)                               \
        printf("%.48s%6.2f   (" _RED_B "FAIL" _RESET ")\n", \
                _test_io_string_,                       \
                (_end_time_ - _start_time_) / CLOCKS_PER_SEC); \
    else                                                \
        printf("%.*s\n%54.2f   (" _RED_B "FAIL" _RESET ")\n", \
                _label_len_, _test_io_string_,          \
                (_end_time_ - _start_time_) / CLOCKS_PER_SEC); \
    flint_printf(__VA_ARGS__);                          \
    return 1;                                           \
} while (0)

#define TEST_MAIN(tests)                                                    \
int                                                                         \
main(int argc, char ** argv)                                                \
{                                                                           \
    char numthreads_str[] = "--numthreads=";                                \
    char thread_str[] = "--thread=";                                        \
    size_t numthreads_str_len = sizeof(numthreads_str) / sizeof(char) - 1;  \
    size_t thread_str_len = sizeof(thread_str) / sizeof(char) - 1;          \
    long int numthreads, thread;                                            \
    int ix, ix_max, numtests;                                               \
                                                                            \
    numtests = sizeof(tests) / sizeof(test_struct);                         \
                                                                            \
    /* If no arguments was put in, simply run all tests. */                 \
    if (argc < 2)                                                           \
    {                                                                       \
        for (ix = 0; ix < numtests; ix++)                                   \
            if ((tests)[ix].test_function())                                \
                flint_abort();                                              \
                                                                            \
        return 0;                                                           \
    }                                                                       \
                                                                            \
    /* If no options not specified, then run the set of functions specified */\
    if (argv[1][0] != '-' && argv[1][1] != '-')                             \
    {                                                                       \
        int jx;                                                             \
                                                                            \
        for (jx = 1; jx < argc; jx++)                                       \
        {                                                                   \
            for (ix = 0; ix < numtests; ix++)                               \
            {                                                               \
                /* If argument equals to test name, run it */               \
                if (strcmp(argv[jx], (tests)[ix].name) == 0)                \
                {                                                           \
                    if ((tests)[ix].test_function())                        \
                        flint_abort();                                      \
                    break;                                                  \
                }                                                           \
            }                                                               \
                                                                            \
            if (ix == numtests)                                             \
            {                                                               \
                flint_throw(FLINT_ERROR, "Could not find test function for %s\n", argv[jx]); \
            }                                                               \
        }                                                                   \
                                                                            \
        return 0;                                                           \
    }                                                                       \
                                                                            \
    /* If the length of first argument is less or equal to the the length of\
     * "--numthreads=" or that the first argument does not start with       \
     * "--numthreads=", then exit. */                                       \
    if (strlen(argv[1]) <= numthreads_str_len                               \
        || strncmp(argv[1], numthreads_str, numthreads_str_len) != 0)       \
    {                                                                       \
        printf("Invalid first option, must be --numthreads=NUM.\n");        \
        return 1;                                                           \
    }                                                                       \
                                                                            \
    /* Parse the argument. An invalid argument is either returned as LONG_MIN,\
     * LONG_MAX or 0. Moreover, we require the argument to be strictly positive.\
     * If the argument is invalid or non-positive, exit. */                 \
    numthreads = strtol(argv[1] + numthreads_str_len, NULL, 10);            \
    if (numthreads == LONG_MAX || numthreads <= 0)                          \
    {                                                                       \
        printf("Invalid parameter for option --numthreads=NUM.\n");         \
        return 1;                                                           \
    }                                                                       \
                                                                            \
    /* Same goes here as for the first argument */                          \
    if (strlen(argv[2]) <= thread_str_len                                   \
        || strncmp(argv[2], thread_str, thread_str_len) != 0)               \
    {                                                                       \
        printf("Invalid second option, must be --thread=NUM.\n");           \
        return 1;                                                           \
    }                                                                       \
                                                                            \
    /* Same goes here as for the first argument, with the addition that the \
     * parameter for this option must be less or equal to the previous argument. */\
    thread = strtol(argv[2] + thread_str_len, NULL, 10);                    \
    if (thread == LONG_MAX || thread <= 0 || thread > numthreads)           \
    {                                                                       \
        printf("Invalid parameter for option --thread=NUM.\n");             \
        return 1;                                                           \
    }                                                                       \
                                                                            \
    /* We split the tests into `numthreads_par' partitions, and run the     \
     * `thread_par'-th partition. */                                        \
    ix = (thread - 1) * ((numtests + numthreads - 1) / numthreads);         \
    ix_max = FLINT_MIN(thread * ((numtests + numthreads - 1) / numthreads), numtests);\
    for (; ix < ix_max; ix++)                                               \
        if ((tests)[ix].test_function())                                    \
            flint_abort();                                                  \
                                                                            \
    return 0;                                                               \
}

#endif
