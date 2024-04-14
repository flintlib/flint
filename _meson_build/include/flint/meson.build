flint_h = configure_file(
  input: '../../src/flint.h.in',
  output: 'flint.h',
  configuration: {
    'FLINT_MAJOR': FLINT_MAJOR,
    'FLINT_MINOR': FLINT_MINOR,
    'FLINT_PATCH': FLINT_PATCH,
    'FLINT_VERSION_FULL': FLINT_VERSION_FULL,
  }
)

flint_config_h = configure_file(
  output: 'flint-config.h',
  configuration: cfg_data,
)

config_h = configure_file(
  output: 'config.h',
  configuration: cfg_data_internal,
)

if gmp_long_long_limb
  gmpcompat_h_in = 'gmpcompat-longlong.h.in'
else
  gmpcompat_h_in = 'gmpcompat.h.in'
endif

gmpcompat_h = configure_file(
  input: '../../src' / gmpcompat_h_in,
  output: 'gmpcompat.h',
  copy: true,
)

flint_mparam_h_in = '../../src/mpn_extras' / flint_mparam / 'flint-mparam.h'

flint_mparam_h = configure_file(
  input: flint_mparam_h_in,
  output: 'flint-mparam.h',
  copy: true,
)

# Copy all the other headers to the build directory so you can do
#   gcc -Ibuild-dir/include ...
foreach h: headers_all
  fs.copyfile(h)
endforeach

headers_all += [
  flint_h,
  flint_config_h,
  config_h,
  gmpcompat_h,
  flint_mparam_h,
]
