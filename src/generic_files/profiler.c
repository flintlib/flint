/*
    Copyright (C) 2007 William Hart
    Copyright (C) 2007 David Harvey
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#if defined(__FreeBSD__) || defined(__DragonFly__) || defined(__NetBSD__) || defined(__OpenBSD__)
# include <sys/types.h>
# include <sys/resource.h>
# include <sys/sysctl.h>
# include <unistd.h>
#elif defined(__APPLE__)
# include <mach/mach.h>
# include <sys/resource.h>
#endif
#if defined(__FreeBSD__) || defined(__DragonFly__)
# include <sys/user.h>
#elif defined(__NetBSD__) || defined(__OpenBSD__)
# include <sys/param.h>
#endif
#include <float.h>
#include <math.h>
#include <stdio.h>
#if defined(__linux__)
# include <stdlib.h>
# include <string.h>
#endif
#include "profiler.h"
#include "longlong.h"

#if FLINT_HAVE_get_cycle_counter
double clock_last[FLINT_NUM_CLOCKS];
double clock_accum[FLINT_NUM_CLOCKS];

void
prof_repeat(double * min, double * max, profile_target_t target, void * arg)
{
    /* Number of timings that were at least DURATION_THRESHOLD microseconds */
    ulong good_count = 0;
    double max_time = DBL_MIN, min_time = DBL_MAX;

    /* First try one loop */
    ulong num_trials = 4;
    double last_time;
    init_clock(0);
    target(arg, num_trials);
    last_time = get_clock(0);

    /* Loop until we have enough good times */
    while (1)
    {
        double per_trial = last_time / num_trials;

        /* If the last recorded time was long enough, record it */
        if (last_time > DURATION_THRESHOLD)
        {
            if (good_count)
            {
                if (per_trial > max_time) max_time = per_trial;
                if (per_trial < min_time) min_time = per_trial;
            }
            else
                max_time = min_time = per_trial;

            if (++good_count == 5) break;
        }

        /* Adjust num_trials so that the elapsed time gravitates towards
           DURATION_TARGET; num_trials can be changed by a factor of
           at most 25%, and must be at least 1 */
        double adjust_ratio;
        if (last_time < 0.0001) last_time = 0.0001;
        adjust_ratio = DURATION_TARGET / last_time;
        if (adjust_ratio > 1.25) adjust_ratio = 1.25;
        if (adjust_ratio < 0.75) adjust_ratio = 0.75;
        num_trials = ceil(adjust_ratio * num_trials);
        /* Just to be safe */
        if (num_trials == 0) num_trials = 1;

        /* Run another trial */
        init_clock(0);
        target(arg, num_trials);
        last_time = get_clock(0);
    }

    /* Store results */
    if (min) *min = min_time;
    if (max) *max = max_time;
}
#endif /* FLINT_HAVE_get_cycle_counter */

#if FLINT_HAVE_getrusage
# define kB (1ULL << 10)
# define MB (1ULL << 20)
# define GB (1ULL << 30)
# define TB (1ULL << 40)
# define PB (1ULL << 50)
static inline void sprint_size(char * str, ulong bytes) {
    double num;
    char fmt[] = "%4.0f";
    if (bytes < kB)
    {
        sprintf(str, "  %3" _WORD_FMT "u", bytes);
        str[5] = ' ', str[6] = 'B';
        return;
    }
    else if (bytes < MB) { num = bytes / (double) kB; str[5] = 'k'; }
    else if (bytes < GB) { num = bytes / (double) MB; str[5] = 'M'; }
    else if (bytes < TB) { num = bytes / (double) GB; str[5] = 'G'; }
    else if (bytes < PB) { num = bytes / (double) TB; str[5] = 'T'; }
    else                 { num = bytes / (double) PB; str[5] = 'P'; }
    fmt[3] = '0' + (num < 10.0) + (num < 100.0);
    sprintf(str, fmt, num);
    str[4] = ' ', str[6] = 'B';
}
#undef kB
#undef MB
#undef GB
#undef TB
#undef PB

void fprint_memory_usage(FILE * fs)
{
#if defined(__linux__)
    FILE * file = fopen("/proc/self/status", "r");
    ulong virt = 0, virtpeak = 0, rss = 0, rsspeak = 0;
    char line[]
        = "virt/peak/rss/peak: 1234567 / 1234567 / 1234567 / 1234567\n";
    char tmp[128];

    while (fgets(tmp, 128, file) != NULL)
    {
        if (strncmp(tmp, "VmSize:", 7) == 0)
            virt = strtoull(tmp + 8, NULL, 10);
        else if (strncmp(tmp, "VmPeak:", 7) == 0)
            virtpeak = strtoull(tmp + 8, NULL, 10);
        else if (strncmp(tmp, "VmRSS:", 6) == 0)
            rss = strtoull(tmp + 7, NULL, 10);
        else if (strncmp(tmp, "VmHWM:", 6) == 0)
            rsspeak = strtoull(tmp + 7, NULL, 10);
    }
    fclose(file);

    sprint_size(line + 20, 1024 * virt);
    sprint_size(line + 30, 1024 * virtpeak);
    sprint_size(line + 40, 1024 * rss);
    sprint_size(line + 50, 1024 * rsspeak);
    fputs(line, fs);
#else
        
    char line[] = "virt/rss/peak: 1234567 / 1234567 / 1234567\n";
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    ulong virt = 0, rss = 0, rsspeak = usage.ru_maxrss;

# if defined(__APPLE__)
    mach_task_basic_info_data_t info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
            (task_info_t) &info, &count);
    virt = info.virtual_size;
    rss = info.resident_size;
# else /* BSD */
    /* ChatJippity generated */
    struct kinfo_proc kp;
    size_t len = sizeof(struct kinfo_proc);
    int mib[4] = {CTL_KERN, KERN_PROC, KERN_PROC_PID, getpid()};
    rsspeak *= 1024; /* given in kilobytes */
    if (sysctl(mib, 4, &kp, &len, NULL, 0) == -1)
        flint_throw(FLINT_ERROR, "report please\n");
#  if defined(__FreeBSD__) || defined(__DragonFly__)
    virt = kp.ki_size;
    rss = kp.ki_rssize * getpagesize();
#  elif defined(__NetBSD__)
    virt = kp.p_vm_msize * getpagesize();
    rss = kp.p_vm_rssize * getpagesize();
#  elif defined(__OpenBSD__)
    virt = kp.p_vmspace.vm_map.size;
    rss = kp.p_vmspace.vm_rssize * getpagesize();
#  endif
# endif /* non-Linux OS */

    sprint_size(line + 15, virt);
    sprint_size(line + 25, rss);
    sprint_size(line + 35, rsspeak);
    fputs(line, fs);
#endif /* Linux or non-Linux  */
}

void print_memory_usage(void) { fprint_memory_usage(stdout); }
#endif /* FLINT_HAVE_getrusage */
