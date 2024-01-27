dnl  m4 macros for gmp assembly code, shared by all CPUs.

dnl  Copyright 1999-2006, 2011 Free Software Foundation, Inc.

dnl  This file is part of the GNU MP Library.
dnl
dnl  The GNU MP Library is free software; you can redistribute it and/or modify
dnl  it under the terms of either:
dnl
dnl    * the GNU Lesser General Public License as published by the Free
dnl      Software Foundation; either version 3 of the License, or (at your
dnl      option) any later version.
dnl
dnl  or
dnl
dnl    * the GNU General Public License as published by the Free Software
dnl      Foundation; either version 2 of the License, or (at your option) any
dnl      later version.
dnl
dnl  or both in parallel, as here.
dnl
dnl  The GNU MP Library is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
dnl  for more details.
dnl
dnl  You should have received copies of the GNU General Public License and the
dnl  GNU Lesser General Public License along with the GNU MP Library.  If not,
dnl  see https://www.gnu.org/licenses/.

dnl  m4 macros for amd64 assembler.

dnl  Copyright 1999-2005, 2008, 2009, 2011-2013, 2017 Free Software Foundation,
dnl  Inc.

dnl  This file is part of the GNU MP Library.
dnl
dnl  The GNU MP Library is free software; you can redistribute it and/or modify
dnl  it under the terms of either:
dnl
dnl    * the GNU Lesser General Public License as published by the Free
dnl      Software Foundation; either version 3 of the License, or (at your
dnl      option) any later version.
dnl
dnl  or
dnl
dnl    * the GNU General Public License as published by the Free Software
dnl      Foundation; either version 2 of the License, or (at your option) any
dnl      later version.
dnl
dnl  or both in parallel, as here.
dnl
dnl  The GNU MP Library is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
dnl  for more details.
dnl
dnl  You should have received copies of the GNU General Public License and the
dnl  GNU Lesser General Public License along with the GNU MP Library.  If not,
dnl  see https://www.gnu.org/licenses/.

define(`R32',
	`ifelse($1,`%rax',`%eax',
		$1,`%rbx',`%ebx',
		$1,`%rcx',`%ecx',
		$1,`%rdx',`%edx',
		$1,`%rsi',`%esi',
		$1,`%rdi',`%edi',
		$1,`%rbp',`%ebp',
		$1,`%r8',`%r8d',
		$1,`%r9',`%r9d',
		$1,`%r10',`%r10d',
		$1,`%r11',`%r11d',
		$1,`%r12',`%r12d',
		$1,`%r13',`%r13d',
		$1,`%r14',`%r14d',
		$1,`%r15',`%r15d')')

define(`R8',
	`ifelse($1,`%rax',`%al',
		$1,`%rbx',`%bl',
		$1,`%rcx',`%cl',
		$1,`%rdx',`%dl',
		$1,`%rsi',`%sil',
		$1,`%rdi',`%dil',
		$1,`%rbp',`%bpl',
		$1,`%r8',`%r8b',
		$1,`%r9',`%r9b',
		$1,`%r10',`%r10b',
		$1,`%r11',`%r11b',
		$1,`%r12',`%r12b',
		$1,`%r13',`%r13b',
		$1,`%r14',`%r14b',
		$1,`%r15',`%r15b')')

define(`L', `.L$1')

define(`ALIGN',
	`ifelse($1,`8',`.p2align	3, 0x90',
		$1,`16',`.p2align	4, 0x90',
		$1,`32',`.p2align	5, 0x90')')
