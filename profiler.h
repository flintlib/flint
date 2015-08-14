/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright 2007 William Hart and David Harvey

******************************************************************************/
#ifndef FLINT_PROFILER_H
#define FLINT_PROFILER_H

#include "flint.h"

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <time.h>
#if defined( _MSC_VER )
#include <intrin.h>
#include "gettimeofday.h"
#pragma intrinsic( __rdtsc )
#else
#include <sys/time.h>
#endif
#if defined (__WIN32) && !defined(__CYGWIN__)
#ifdef __cplusplus
FLINT_DLL void  GetSystemTimeAsFileTime(FILETIME*);

static __inline__ int gettimeofday(struct timeval * p, void * tz)
{
   union {
      slong slong ns100; 
      FILETIME ft;
   } now;

    GetSystemTimeAsFileTime(&(now.ft));
    p->tv_usec=(slong)((now.ns100 / WORD(10)L) % WORD(1000000)L );
    p->tv_sec= (slong)((now.ns100-(WORD(116444736000000000)L))/WORD(10000000)L);
	
    return 0;
}
#else
FLINT_DLL int gettimeofday(struct timeval * p, void * tz);
#endif
#elif !defined(_MSC_VER)
#include <sys/resource.h>
#endif
#undef ulong
#define ulong mp_limb_t

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct  
{
    ulong size;
    ulong peak;
    ulong hwm;
    ulong rss;
} meminfo_t[1];

FLINT_DLL void get_memory_usage(meminfo_t meminfo);

typedef struct
{
    slong cpu;
    slong wall;
} timeit_t[1];

static __inline__
void timeit_start(timeit_t t)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    t->wall = - tv.tv_sec * 1000 - tv.tv_usec / 1000;
    t->cpu = - clock() * 1000 / CLOCKS_PER_SEC;
}

static __inline__
void timeit_stop(timeit_t t)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    t->wall += tv.tv_sec * 1000 + tv.tv_usec / 1000;
    t->cpu += clock() * 1000 / CLOCKS_PER_SEC;
}

/******************************************************************************

    Timer based on the cycle counter
   
******************************************************************************/

#define FLINT_NUM_CLOCKS 20

#define FLINT_CLOCKSPEED 3100000000.0

extern double clock_last[FLINT_NUM_CLOCKS];
extern double clock_accum[FLINT_NUM_CLOCKS];

static __inline__ 
double get_cycle_counter()
{
#if defined( _MSC_VER )
    return (double)__rdtsc();
#else
    unsigned int hi;
   unsigned int lo;

   __asm("rdtsc; movl %%edx,%0; movl %%eax,%1" 
       : "=r" (hi), "=r" (lo)
       : 
       : "%edx", "%eax");

   return (double) hi * (1 << 30) * 4 + lo;
#endif
}

#define FLINT_CLOCK_SCALE_FACTOR (1000000.0 / FLINT_CLOCKSPEED)

static __inline__ 
void init_clock(int n)
{
   clock_accum[n] = 0.0;
}

static __inline__ 
void init_all_clocks()
{
   int i;
   for (i = 0; i < FLINT_NUM_CLOCKS; i++)
      clock_accum[i] = 0.0;
}

static __inline__
double get_clock(int n)
{
   return clock_accum[n] * FLINT_CLOCK_SCALE_FACTOR;
}

static __inline__ 
void start_clock(int n)
{
   clock_last[n] = get_cycle_counter();
}

static __inline__ 
void stop_clock(int n)
{
   double now = get_cycle_counter();
   clock_accum[n] += (now - clock_last[n]);
}

/******************************************************************************

    Framework for repeatedly sampling a single target

******************************************************************************/

static __inline__ 
void prof_start()
{
   start_clock(0);
}

static __inline__ 
void prof_stop()
{
   stop_clock(0);
}

typedef void (*profile_target_t)(void* arg, ulong count);

FLINT_DLL void prof_repeat(double* min, double* max, profile_target_t target, void* arg);

#define DURATION_THRESHOLD 5000.0

#define DURATION_TARGET 10000.0

/******************************************************************************

    Simple timing macros

******************************************************************************/

#define TIMEIT_PRINT(__timer, __reps) \
    flint_printf("cpu/wall(s): %g %g\n", \
        __timer->cpu*0.001/__reps, __timer->wall*0.001 / __reps);

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
    } while (0);

#define TIMEIT_START \
    do { \
        timeit_t __timer; slong __reps; \
        TIMEIT_REPEAT(__timer, __reps)

#define TIMEIT_STOP \
        TIMEIT_END_REPEAT(__timer, __reps) \
        TIMEIT_PRINT(__timer, __reps) \
    } while (0);

#define TIMEIT_ONCE_START \
    do \
    { \
      timeit_t __timer; \
      timeit_start(__timer); \
      do { \

#define TIMEIT_ONCE_STOP \
      } while (0); \
      timeit_stop(__timer); \
      TIMEIT_PRINT(__timer, 1) \
    } while (0); \

#define SHOW_MEMORY_USAGE \
    do { \
        meminfo_t meminfo; \
        get_memory_usage(meminfo); \
        flint_printf("virt/peak/res/peak(MB): %.2f %.2f %.2f %.2f\n", \
            meminfo->size / 1024.0, meminfo->peak / 1024.0, \
            meminfo->rss / 1024.0, meminfo->hwm / 1024.0); \
    } while (0);

#ifdef __cplusplus
}
#endif

#endif

