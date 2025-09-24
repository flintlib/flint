header_install_kwargs = {
  'install_dir': get_option('includedir') / 'flint',
  'install': true,
}

flint_h = configure_file(
  input: '../../src/flint.h.in',
  output: 'flint.h',
  configuration: {
    'FLINT_MAJOR': FLINT_MAJOR,
    'FLINT_MINOR': FLINT_MINOR,
    'FLINT_PATCH': FLINT_PATCH,
    'FLINT_VERSION_FULL': FLINT_VERSION_FULL,
  },
  kwargs: header_install_kwargs,
)

flint_config_h = configure_file(
  output: 'flint-config.h',
  configuration: cfg_data,
  kwargs: header_install_kwargs,
)

config_h = configure_file(
  output: 'config.h',
  configuration: cfg_data_internal,
  kwargs: header_install_kwargs,
)

if gmp_long_long_limb
  gmpcompat_h_in = 'gmpcompat-longlong.h.in'
else
  gmpcompat_h_in = 'gmpcompat.h.in'
endif

gmpcompat_h = fs.copyfile(
  '../../src' / gmpcompat_h_in,
  'gmpcompat.h',
  kwargs: header_install_kwargs,
)

flint_mparam_h_in = '../../src/mpn_extras' / flint_mparam / 'flint-mparam.h'

flint_mparam_h = fs.copyfile(
  flint_mparam_h_in,
  'flint-mparam.h',
  kwargs: header_install_kwargs,
)

# Copy all the other headers to the build directory so you can do
#   gcc -Ibuild-dir/include ...
foreach h: headers_all
  fs.copyfile(h, kwargs: header_install_kwargs)
endforeach

headers_all += [
  flint_h,
  flint_config_h,
  config_h,
  gmpcompat_h,
  flint_mparam_h,
]
