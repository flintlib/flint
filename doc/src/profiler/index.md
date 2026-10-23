---
title: profiler.h
header: profiler.h
description: Timers for quick benchmarks, cycle counters and memory usage reports.
category: Foundations
order: 2
---

This module provides simple tools for timing code and reporting memory
usage. It is used by FLINT's own profiling programs (the `profile`
directories, built with `make profile`) and example programs. It is not
a substitute for a real profiler such as `perf`.

There are three groups of functions. The `timeit` timers and the
`TIMEIT_*` macros measure processor time and wall time in milliseconds
and are available everywhere. The cycle-counter clocks and `prof_repeat`
read the processor's time-stamp counter and are only available on
x86-64 and ARM64 (and with MSVC); the header defines
`FLINT_HAVE_get_cycle_counter` when they are. The memory usage reports
need `getrusage` or `/proc` and are available on Linux, the BSDs and
macOS (`FLINT_HAVE_getrusage`).

## Timers {#timeit}

::: type
timeit_t

A timer, an array of one structure with fields `cpu` and `wall`
(`slong`). After `timeit_stop`, they hold the elapsed processor time and
wall time in milliseconds (in microseconds for `timeit_stop_us`).
:::

::: function
void timeit_start(timeit_t t)
void timeit_stop(timeit_t t)

Start and stop the timer. After `timeit_stop`, `t->cpu` is the processor
time of the process (from `clock`) and `t->wall` the wall time (from
`gettimeofday`) elapsed since `timeit_start`, both in milliseconds.
Processor time adds up over all threads of the process, so for a parallel
computation it can exceed the wall time.
:::

::: function
slong timeit_query_wall(timeit_t t)

Returns the wall time in milliseconds elapsed since `timeit_start(t)`,
without stopping the timer.
:::

::: function
void timeit_start_us(timeit_t t)
void timeit_stop_us(timeit_t t)

Start and stop a wall-time measurement in microseconds, stored in
`t->wall`. The field `t->cpu` is not set.
:::

::: function
void timeit_print(timeit_t t, ulong reps)

Prints `cpu/wall(s): ` followed by the processor time and wall time of
the timer, in seconds, divided by `reps`.

::: details
For example:

```
timeit_t t0;
timeit_start(t0);
/* the code to time */
timeit_stop(t0);
timeit_print(t0, 1);     /* one repetition */
```
:::
:::

::: macro
TIMEIT_REPEAT(timer, reps)
TIMEIT_END_REPEAT(timer, reps)

Run the code between the two macros repeatedly, starting with one
repetition and multiplying the number of repetitions by 10 until the
processor time of a round is at least 100 milliseconds. On exit, the
`slong` variable `reps` holds the number of repetitions of the last round
and the `timeit_t` variable `timer` its time, which the caller divides by
`reps`.
:::

::: macro
TIMEIT_START
TIMEIT_STOP
TIMEIT_STOP_VALUES(tcpu, twall)

Run the code between `TIMEIT_START` and `TIMEIT_STOP` repeatedly, as
`TIMEIT_REPEAT`, and print the processor and wall time of one repetition
with `timeit_print`. `TIMEIT_STOP_VALUES` instead stores the two times,
in seconds, in the `double` variables `tcpu` and `twall`.

::: details
The code is run at least once and typically for 0.1 to 1 second, so it
must be safe to execute repeatedly. For example:

```
TIMEIT_START
fmpz_mul(c, a, b);
TIMEIT_STOP
```
:::
:::

::: macro
TIMEIT_ONCE_START
TIMEIT_ONCE_STOP

Run the code between the two macros exactly once and print the elapsed
processor and wall time with `timeit_print`. The result is imprecise if
the code runs for less than a few milliseconds.
:::

## Cycle counters {#clocks}

These clocks count ticks of the time-stamp counter (`rdtsc` on x86-64,
the virtual counter `cntvct_el0` on ARM64) and convert ticks to
microseconds by dividing by `FLINT_CLOCKSPEED`, which defaults to
$3.1 \cdot 10^9$ ticks per second. The conversion is only right if
`FLINT_CLOCKSPEED` is defined to the actual tick rate before including
`profiler.h`. On current x86-64 processors the counter runs at a fixed
rate close to the nominal frequency, independently of the current clock
speed; on ARM64 it usually runs at a much lower fixed frequency (for
example 24 MHz on Apple processors), so the default gives wrong absolute
times there.

::: function
double get_cycle_counter(void)

Returns the current value of the time-stamp counter.
:::

::: macro
FLINT_NUM_CLOCKS
FLINT_CLOCKSPEED
FLINT_CLOCK_SCALE_FACTOR

`FLINT_NUM_CLOCKS` (20) is the number of clocks. `FLINT_CLOCKSPEED` is
the assumed number of counter ticks per second, and
`FLINT_CLOCK_SCALE_FACTOR` is $10^6$ divided by it, the number of
microseconds per tick.
:::

::: function
void init_clock(int n)
void init_all_clocks(void)
void start_clock(int n)
void stop_clock(int n)
double get_clock(int n)

Accumulating clocks numbered $0 \le n <$ `FLINT_NUM_CLOCKS`.
`init_clock(n)` resets clock $n$ to zero and `init_all_clocks` resets all
of them. Each pair `start_clock(n)`, `stop_clock(n)` adds the ticks
elapsed between the two calls to clock $n$, and `get_clock(n)` returns
the accumulated time in microseconds.

::: details
```
init_all_clocks();
start_clock(1);
/* code to time */
stop_clock(1);
flint_printf("%.3f us\n", get_clock(1));
```
:::
:::

::: variable
double clock_last[FLINT_NUM_CLOCKS]
double clock_accum[FLINT_NUM_CLOCKS]

The state of the clocks: the counter value at the last `start_clock` and
the accumulated ticks.
:::

::: function
void prof_start(void)
void prof_stop(void)

`start_clock(0)` and `stop_clock(0)`, for use in the target functions of
`prof_repeat`.
:::

::: type
profile_target_t

A pointer to a function `void target(void * arg, ulong count)` that runs
the code to be measured `count` times, enclosing each run (or the whole
loop) in `prof_start` and `prof_stop`.
:::

::: function
void prof_repeat(double * min, double * max, profile_target_t target, void * arg)

Measures the time of one run of `target`: calls `target(arg, count)`
repeatedly, adjusting `count` so that each call takes about
`DURATION_TARGET` microseconds of clock 0, until five calls have taken
more than `DURATION_THRESHOLD` microseconds. Sets `*min` and `*max` to the
smallest and largest time per run, in microseconds, of these five calls.
Either pointer may be `NULL`.

::: details
The number of runs per call starts at 4 and changes by at most 25% from
one call to the next. A typical target function is:

```
typedef struct { ulong a, b; } myfunc_t;

void sample_myfunc(void * arg, ulong count)
{
    myfunc_t * params = (myfunc_t *) arg;
    ulong i;

    prof_start();
    for (i = 0; i < count; i++)
        myfunc(params->a, params->b);
    prof_stop();
}
```

which is then measured with
`prof_repeat(&min, &max, sample_myfunc, &params)`. As the times are
derived from clock 0, the remarks above about `FLINT_CLOCKSPEED` apply.
:::
:::

::: macro
DURATION_TARGET
DURATION_THRESHOLD

The target duration of one timing round of `prof_repeat` and the minimum
duration for a round to count, in microseconds (10000 and 5000).
:::

## Memory usage {#memory}

::: function
void fprint_memory_usage(FILE * fs)
void print_memory_usage(void)

Print the memory usage of the current process to `fs` or to standard
output, in one line. On Linux, the line is
`virt/peak/rss/peak: ` followed by the current and peak virtual memory
size and the current and peak resident set size (from
`/proc/self/status`); on macOS and the BSDs it is `virt/rss/peak: `
followed by the virtual size, the resident set size and the peak resident
set size. Sizes are printed with units, as in `24.8 MB`.
:::

## Internal functions {#internal}

::: function
int gettimeofday(struct timeval * p, void * tz)

Replacement for the POSIX function, declared on Windows with compilers
other than MSVC, where the C library does not provide it.
:::
