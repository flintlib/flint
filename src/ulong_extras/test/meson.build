test_exe = executable('main',
  'main.c',
  dependencies: flint_deps,
  link_with: libflint,
  include_directories: src_inc,
  install: false,
)
test('ulong_extras', test_exe)
