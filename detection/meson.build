
fft_small_opt = get_option('fft_small')

# check for headers used by fft_small
has_avx = cc.has_header('immintrin.h')
has_neon = cc.has_header('arm_neon.h')

fft_small_supported = has_avx or has_neon

# Error if fft_small enabled but not possible, otherwise disable if auto
fft_small_opt = fft_small_opt.require(fft_small_supported,
    error_message: 'host CPU does not support AVX/NEON. Set fft_small to disabled or auto')

# If fft_small was auto we enable it if possible
fft_small_opt = fft_small_opt.enable_auto_if(fft_small_supported)
have_fft_small = fft_small_opt.enabled()

message('FFT_SMALL: ', have_fft_small ? 'enabled' : 'disabled')

if have_fft_small
    cfg_data.set('FLINT_HAVE_FFT_SMALL', 1,
        description: 'Define to use the fft_small module')
endif

if have_fft_small
  c_files_all += [
    'src/fft_small/default_ctx.c',
    'src/fft_small/fmpz_poly_mul.c',
    'src/fft_small/mpn_helpers.c',
    'src/fft_small/mpn_mul.c',
    'src/fft_small/mulmod_statisfies_bounds.c',
    'src/fft_small/nmod_poly_mul.c',
    'src/fft_small/sd_fft.c',
    'src/fft_small/sd_fft_ctx.c',
    'src/fft_small/sd_ifft.c',
  ]
  # Maybe fft_small.h should only be added conditionally as well...
  # headers_all += ['src/fft_small.h']
endif

assembly_opt = get_option('assembly')

# Might be better to check for known good operating systems rather than
# exclude known bad ones.
assembly_opt = assembly_opt.require(host_machine.system() != 'windows',
    error_message: 'No assembly support for Windows')

# Checking for assembly support needs to be more complicated than this. ARM
# needs to be v8 and not sure exactly what to check for x86_64. For now
# the checks below would attempt to build assembly on some systems that do not
# support it. The feature can be disabled with -Dassembly=disabled to build on
# those systems.
if host_machine.cpu_family() == 'aarch64'
  assembly_cpu_supported = true
  assembly_dir = 'arm64'
  # FLINT_HAVE_ASSEMBLY_armv8
elif host_machine.cpu_family() == 'x86_64'
  assembly_cpu_supported = true
  assembly_dir = 'x86_64/broadwell'
  # FLINT_HAVE_ASSEMBLY_x86_64_adx
else
  assembly_cpu_supported = false
endif

assembly_opt = assembly_opt.require(assembly_cpu_supported,
    error_message: 'No assembly support for host CPU')

assembly_opt = assembly_opt.enable_auto_if(true)
have_assembly = assembly_opt.enabled()

message('ASSEMBLY: ', have_assembly ? 'enabled' : 'disabled')

if have_assembly

  asm_deps = [
    '../config.m4',
    '../src/mpn_extras/asm-defs.m4',
  ]

  if assembly_dir == 'arm64'
    asm_files = [
      'src/mpn_extras/arm64/mul_2.asm',
      'src/mpn_extras/arm64/mulhigh_hard.asm',
      'src/mpn_extras/arm64/mul_hard.asm',
      'src/mpn_extras/arm64/sqrhigh_hard.asm',
      'src/mpn_extras/arm64/mulhigh_basecase.asm',
      'src/mpn_extras/arm64/sqr_hard.asm',
    ]
    asm_deps += [
      '../src/mpn_extras/arm64/darwin.m4',
      '../src/mpn_extras/arm64/arm64-defs.m4',
    ]
  elif assembly_dir == 'x86_64/broadwell'
    asm_files = [
      'src/mpn_extras/x86_64/broadwell/mullow_hard.asm',
      'src/mpn_extras/x86_64/broadwell/2add_n_inplace.asm',
      'src/mpn_extras/x86_64/broadwell/mulhigh_hard.asm',
      'src/mpn_extras/x86_64/broadwell/sqrhigh_basecase_odd.asm',
      'src/mpn_extras/x86_64/broadwell/mul_hard.asm',
      'src/mpn_extras/x86_64/broadwell/mulhigh_normalised_hard.asm',
      'src/mpn_extras/x86_64/broadwell/sqrhigh_hard.asm',
      'src/mpn_extras/x86_64/broadwell/mullow_basecase.asm',
      'src/mpn_extras/x86_64/broadwell/mulhigh_basecase.asm',
      'src/mpn_extras/x86_64/broadwell/sqrhigh_basecase_even.asm',
      'src/mpn_extras/x86_64/broadwell/sqr_hard.asm',
    ]
    asm_deps += [
      '../src/mpn_extras/x86_64/broadwell/x86_64-defs.m4',
    ]
    if host_machine.system() == 'darwin'
      asm_deps += [
        '../src/mpn_extras/x86_64/broadwell/darwin.m4',
      ]
    endif
  else
    error('Unknown assembly_dir: ', assembly_dir)
  endif

  m4_prog = find_program('m4', native: true)

  s_files = []
  foreach asm_file: asm_files
    s_filename = fs.stem(asm_file) + '.s'
    s_file = custom_target(s_filename,
      input: ['..' / asm_file] + asm_deps,
      output: s_filename,
      command: [m4_prog, '-I..', '@INPUT0@'],
      capture: true,
    )
    s_files += [s_file]
  endforeach

endif
