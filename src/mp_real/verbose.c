/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "flint.h"
#if FLINT_USES_PTHREAD
#include <pthread.h>
#endif
#include "profiler.h"
#include "thread_support.h"
#include "mp_real.h"
#include "impl.h"

/* Progress reports and timings on stderr for huge computations.

   The level is set by mp_real_set_verbose, or by the environment
   variable FLINT_MP_REAL_VERBOSE (read once, at the first message
   check) when nothing has been set:

       0  silent (the default)
       1  progress: one line every MP_REAL_VERBOSE_INTERVAL seconds at
          most from the long loops (the merges of the binary
          splittings, the slices of the bit-burst evaluations, the
          Newton steps), plus the start and end of the phases that
          take longer than the interval
       2  everything: every phase of every evaluation above a few
          thousand limbs, with its size, thread budget and time

   Each line carries the time since the first message of the process
   and the id of the calling thread (0 for the thread that first
   reported, which is normally the main thread).  The messages are
   serialized by a mutex, so lines from different threads never
   interleave. */

static int mp_real_verbose_level = -1;
static double verbose_epoch = -1.0;
static double verbose_last_progress = -1e300;
#if FLINT_USES_PTHREAD
static pthread_mutex_t verbose_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_t verbose_threads[256];
static int verbose_num_threads = 0;
#endif

void
mp_real_set_verbose(int level)
{
    mp_real_verbose_level = level;
}

int
mp_real_get_verbose(void)
{
    if (mp_real_verbose_level < 0)
    {
        const char * s = getenv("FLINT_MP_REAL_VERBOSE");
        mp_real_verbose_level = (s == NULL) ? 0 : atoi(s);
        if (mp_real_verbose_level < 0)
            mp_real_verbose_level = 0;
    }
    return mp_real_verbose_level;
}

double
_mp_real_verbose_time(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double) t.tv_sec + 1e-6 * (double) t.tv_usec;
}

/* the small integer id of the calling thread, in order of first
   appearance */
static int
_thread_id(void)
{
#if FLINT_USES_PTHREAD
    pthread_t me = pthread_self();
    int i;
    for (i = 0; i < verbose_num_threads; i++)
        if (pthread_equal(verbose_threads[i], me))
            return i;
    if (verbose_num_threads < 256)
        verbose_threads[verbose_num_threads++] = me;
    return verbose_num_threads - 1;
#else
    return 0;
#endif
}

static void
_vlog(int progress, const char * fmt, va_list ap)
{
    double now = _mp_real_verbose_time();

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&verbose_mutex);
#endif
    if (verbose_epoch < 0.0)
        verbose_epoch = now;
    if (!progress || now - verbose_last_progress >= MP_REAL_VERBOSE_INTERVAL)
    {
        if (progress)
            verbose_last_progress = now;
        fprintf(stderr, "[mp_real %9.3f s, thread %d, %d avail] ",
            now - verbose_epoch, _thread_id(), flint_get_num_threads());
        flint_vfprintf(stderr, fmt, ap);
        fputc('\n', stderr);
        fflush(stderr);
    }
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&verbose_mutex);
#endif
}

void
_mp_real_log(const char * fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    _vlog(0, fmt, ap);
    va_end(ap);
}

void
_mp_real_progress(const char * fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    _vlog(1, fmt, ap);
    va_end(ap);
}

/* the resident memory and its high-water mark in MB (Linux; -1
   elsewhere), for the phase reports */
static void
_mem_mb(double * cur, double * peak)
{
    *cur = *peak = -1.0;
#ifdef __linux__
    {
        FILE * f = fopen("/proc/self/status", "r");
        char line[256];
        long kb;
        if (f == NULL)
            return;
        while (fgets(line, sizeof(line), f) != NULL)
        {
            if (strncmp(line, "VmRSS:", 6) == 0 && sscanf(line + 6, "%ld", &kb) == 1)
                *cur = kb / 1024.0;
            else if (strncmp(line, "VmHWM:", 6) == 0 && sscanf(line + 6, "%ld", &kb) == 1)
                *peak = kb / 1024.0;
        }
        fclose(f);
    }
#endif
}

/* a phase: start at level 2, end at level 2 or when it took at least
   the interval at level 1 (with the process's resident memory and
   its peak so far, where available) */
void
_mp_real_phase_start(mp_real_phase_t * ph, const char * name, slong n)
{
    ph->t0 = _mp_real_verbose_time();
    ph->name = name;
    ph->n = n;
    if (mp_real_get_verbose() >= 2)
    {
        double cur, peak;
        _mem_mb(&cur, &peak);
        if (cur >= 0.0)
            _mp_real_log("%s: start, %wd limbs, RSS %.0f MB (peak %.0f MB)",
                name, n, cur, peak);
        else
            _mp_real_log("%s: start, %wd limbs", name, n);
    }
}

void
_mp_real_phase_end(mp_real_phase_t * ph)
{
    double dt = _mp_real_verbose_time() - ph->t0;
    if (mp_real_get_verbose() >= 2
        || (mp_real_get_verbose() >= 1 && dt >= MP_REAL_VERBOSE_INTERVAL))
    {
        double cur, peak;
        _mem_mb(&cur, &peak);
        if (cur >= 0.0)
            _mp_real_log("%s: done, %wd limbs, %.3f s, RSS %.0f MB (peak %.0f MB)",
                ph->name, ph->n, dt, cur, peak);
        else
            _mp_real_log("%s: done, %wd limbs, %.3f s", ph->name, ph->n, dt);
    }
}
