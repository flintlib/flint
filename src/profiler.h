/*
    Copyright (C) 2007 William Hart
    Copyright (C) 2007 David Harvey
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_PROFILER_H
#define FLINT_PROFILER_H

#if defined(_MSC_VER)
# include <intrin.h>
# include "gettimeofday.h"
# pragma intrinsic(__rdtsc)
#else
# include <sys/time.h>
#endif
#include <time.h>
#include "flint.h"

#if defined(_MSC_VER) || defined(__x86_64__) || defined(__aarch64__)
# define FLINT_HAVE_get_cycle_counter  1
#endif
#if (defined(__unix__) && !defined(__CYGWIN__)) || defined(__APPLE__)
# define FLINT_HAVE_getrusage   1
#endif

#define FLINT_NUM_CLOCKS 20
#if !defined(FLINT_CLOCKSPEED)
# define FLINT_CLOCKSPEED 3100000000.0
#endif
#define FLINT_CLOCK_SCALE_FACTOR (1000000.0 / FLINT_CLOCKSPEED)

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__WIN32) && !defined(__CYGWIN__) && !defined(_MSC_VER)
int gettimeofday(struct timeval * p, void * tz);
#endif

/* memory usage **************************************************************/

#if FLINT_HAVE_getrusage
#if FLINT_HAVE_FILE
void fprint_memory_usage(FILE *);
#endif
void print_memory_usage(void);
#endif

/* timeit ********************************************************************/

typedef struct { slong cpu, wall; } timeit_t[1];

static inline void timeit_start(timeit_t t)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    t->wall = - tv.tv_sec * 1000 - tv.tv_usec / 1000;
    t->cpu = - clock() * 1000 / CLOCKS_PER_SEC;
}

static inline slong timeit_query_wall(timeit_t t)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    return t->wall + tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

static inline void timeit_stop(timeit_t t)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    t->wall += tv.tv_sec * 1000 + tv.tv_usec / 1000;
    t->cpu += clock() * 1000 / CLOCKS_PER_SEC;
}

static inline void timeit_start_us(timeit_t t)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    t->wall = - tv.tv_sec * 1000000 - tv.tv_usec;
}

static inline void timeit_stop_us(timeit_t t)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    t->wall += tv.tv_sec * 1000000 + tv.tv_usec;
}

static inline void timeit_print(timeit_t timer, ulong reps)
{
    flint_printf("cpu/wall(s): %g %g\n", timer->cpu * 0.001 / reps,
                                   timer->wall * 0.001 / reps);
}

#define TIMEIT_REPEAT(__timer, __reps) \
    do \
    { \
        slong __timeit_k; \
        __reps = 1; \
        while (1) \
        { \
            timeit_start(__timer); \
            for (__timeit_k = 0; __timeit_k < __reps; __timeit_k++) \
            {

#define TIMEIT_END_REPEAT(__timer, __reps) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= 100) \
                break; \
            __reps *= 10; \
        } \
    } while (0)

#define TIMEIT_START \
    do { \
        timeit_t __timer; slong __reps; \
        TIMEIT_REPEAT(__timer, __reps)

#define TIMEIT_STOP \
        TIMEIT_END_REPEAT(__timer, __reps); \
        timeit_print(__timer, __reps); \
    } while (0)

#define TIMEIT_STOP_VALUES(tcpu, twall) \
        TIMEIT_END_REPEAT(__timer, __reps); \
        (tcpu) = __timer->cpu*0.001 / __reps; \
        (twall) = __timer->wall*0.001 / __reps; \
    } while (0)

#define TIMEIT_ONCE_START \
    do \
    { \
      timeit_t __timer; \
      timeit_start(__timer); \
      do {

#define TIMEIT_ONCE_STOP \
      } while (0); \
      timeit_stop(__timer); \
      timeit_print(__timer, 1); \
    } while (0)

/* raw clock cycles **********************************************************/

#if FLINT_HAVE_get_cycle_counter
static inline double get_cycle_counter(void)
{
#if defined(_MSC_VER)
    return __rdtsc();
#elif defined(__x86_64__)
    unsigned int hi, lo;
    __asm__ volatile ("rdtsc; movl %%edx,%0; movl %%eax,%1"
            : "=r" (hi), "=r" (lo) : : "%edx", "%eax");
    return hi * (double) (WORD(1) << 32) + lo;
#elif defined(__aarch64__)
    ulong val;
    __asm__ volatile ("mrs %0, cntvct_el0" : "=r" (val));
    return val;
#endif
}

FLINT_DLL extern double clock_last[FLINT_NUM_CLOCKS];
FLINT_DLL extern double clock_accum[FLINT_NUM_CLOCKS];

static inline void init_clock(int n)
{
    clock_accum[n] = 0.0;
}

static inline void init_all_clocks(void)
{
    for (int i = 0; i < FLINT_NUM_CLOCKS; i++)
        init_clock(i);
}

static inline double get_clock(int n)
{
    return clock_accum[n] * FLINT_CLOCK_SCALE_FACTOR;
}

static inline void start_clock(int n)
{
    clock_last[n] = get_cycle_counter();
}

static inline void stop_clock(int n)
{
    double now = get_cycle_counter();
    clock_accum[n] += (now - clock_last[n]);
}

#define DURATION_TARGET 10000.0
#define DURATION_THRESHOLD 5000.0

static inline void prof_start(void) { start_clock(0); }
static inline void prof_stop(void) { stop_clock(0); }

typedef void (*profile_target_t)(void* arg, ulong count);

void prof_repeat(double* min, double* max, profile_target_t target, void* arg);

#endif /* FLINT_HAVE_get_cycle_counter */

#ifdef __cplusplus
}
#endif

#endif /* FLINT_PROFILER_H */
