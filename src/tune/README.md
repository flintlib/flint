# Tuning-suite

## Requirements

- Either requires that
  * `clock_gettime` is available on the system, or that
  * compiler is GCC compatible and architecture is x86.

## Usage

### Set CPU frequency

If tuner is using clock ticks (currently only for x86-64), you can specify your
clock frequency by pushing `export FLINT_CPU_FREQUENCY=3.1e9` to tell the tuner
that your CPU frequency is 3.1 GHz.
