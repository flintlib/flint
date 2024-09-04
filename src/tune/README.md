# Tuning-suite

Currently working-in-progress, but feedback is much appreciated.

## Usage

Run `make tune` followed by `./build/tuneup`. This pushes all the optimized
parameters into `stdout`. These can then be used to optimize the parameter file
`flint-mparam.h` for your system.

### Set CPU frequency

If tuner is using clock ticks (currently only for x86-64), you can specify your
clock frequency by pushing `export FLINT_CPU_FREQUENCY=3.1e9` to tell the tuner
that your CPU frequency is 3.1 GHz.

### Options

Currently, no command-line options are allowed apart from `-h` and `--help` to
display the usual help message.

However, it would be optimal to be able to specify:

- Functions intended to benchmark (currently does all available)
- Minimum number of runs
- Warmup runs (currently, 10 is the default)
- Minimum amount of time to run each function (?)
- Precision required to terminate successfully (currently 1.25 %)
- Percentage of runs required to be within said precision (currently 13.5 %)

## Issues

Please open up any issues at <https://github.com/flintlib/flint/issues>.

## Requirements

- FLINT was built with Autotools.
- Either that
  * `clock_gettime` is available on the system, or that
  * compiler is GCC compatible and architecture is x86.

## How it works

The program works in the following order:

1. Parses options
2. Sets default values
3. For each function (that is, each variant of each function) tested:
   a. Run a couple of warmups, which are trashed
   b. Run hotlaps, of which the time is saved
   c. Check if there is a smallest time elapsed $t$ for running a function of
      which at least $k$ of the runs have a time in the interval
      $[t, (1 + p) t]$, where $p$ is the precision and $k / n$ is the percentage
      of runs required to be within said precision, where $n$ is the total
      number of runs. If no such $t$ was found, abort.
4. With all $t$ obtained from each families of functions, determine cutoff
   points, methods used, etc.
5. Print the associated `#define` into `stdout`.
