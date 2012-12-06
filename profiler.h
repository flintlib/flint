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

#undef ulong /* interferes with system includes, redefined by flint.h below */
#include <time.h>
#include <sys/time.h>
#if defined (__WIN32) && !defined(__CYGWIN__)
#ifdef __cplusplus
void  GetSystemTimeAsFileTime(FILETIME*);

static __inline__ int gettimeofday(struct timeval * p, void * tz)
{
   union {
      long long ns100; 
      FILETIME ft;
   } now;

    GetSystemTimeAsFileTime(&(now.ft));
    p->tv_usec=(long)((now.ns100 / 10LL) % 1000000LL );
    p->tv_sec= (long)((now.ns100-(116444736000000000LL))/10000000LL);
	
    return 0;
}
#else
int gettimeofday(struct timeval * p, void * tz);
#endif
#else
#include <sys/resource.h>
#endif
#define ulong unsigned long

#ifndef FLINT_PROFILER_H
#define FLINT_PROFILER_H

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    long cpu;
    long wall;
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

#define FLINT_CLOCKSPEED 2400000000.0

extern double clock_last[FLINT_NUM_CLOCKS];
extern double clock_accum[FLINT_NUM_CLOCKS];

static __inline__ 
double get_cycle_counter()
{
   unsigned int hi;
   unsigned int lo;

   __asm("rdtsc; movl %%edx,%0; movl %%eax,%1" 
       : "=r" (hi), "=r" (lo)
       : 
       : "%edx", "%eax");

   return (double) hi * (1 << 30) * 4 + lo;
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

typedef void (*profile_target_t)(void* arg, unsigned long count);

void prof_repeat(double* min, double* max, profile_target_t target, void* arg);

#define DURATION_THRESHOLD 5000.0

#define DURATION_TARGET 10000.0

#ifdef __cplusplus
}
#endif

#endif

