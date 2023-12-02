/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include "templates.h"

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
    FLINT_TEST_INIT(state);                             \
    printf(#label "....");                              \
    fflush(stdout);                                     \

#define TEST_GR_FUNCTION_START(label, state, count_success, count_domain, count_unable)\
int CAT(test, label)(void)                              \
{                                                       \
    slong count_success = 0, count_unable = 0, count_domain = 0;\
    FLINT_TEST_INIT(state);                             \
    printf(#label "....");                              \
    fflush(stdout);                                     \

#define TEST_TEMPLATE_FUNCTION_START(T, label, state)   \
int TEMPLATE3(test, T, label)(void)                     \
{                                                       \
    FLINT_TEST_INIT(state);                             \
    printf(TEMPLATE_STR(T) "_" #label "....");          \
    fflush(stdout);                                     \

#define TEST_TEMPLATE2_FUNCTION_START(T, label1, label2, state)\
int TEMPLATE5(test, T, label1, T, label2)(void)         \
{                                                       \
    FLINT_TEST_INIT(state);                             \
    printf(TEMPLATE_STR(T) "_" #label1 "_"              \
           TEMPLATE_STR(T) "_" #label2 "....");         \
    fflush(stdout);                                     \

#define TEST_FUNCTION_END(state)                        \
    FLINT_TEST_CLEANUP(state);                          \
    printf("PASS\n");                                   \
    return 0;                                           \
}

#define TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable)\
    FLINT_TEST_CLEANUP(state);                          \
    printf(" [" WORD_FMT "d success, "                  \
                WORD_FMT "d domain, "                   \
                WORD_FMT "d unable] PASS\n",            \
                count_success, count_domain, count_unable);\
    return 0;                                           \
}

#define TEST_FUNCTION_END_SKIPPED(state)                \
    FLINT_TEST_CLEANUP(state);                          \
    printf("SKIPPED\n");                                \
    return 0;                                           \
}

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
