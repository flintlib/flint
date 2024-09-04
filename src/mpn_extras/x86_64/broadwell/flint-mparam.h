/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Tuning for generic CPU. Is probably far from optimal. */

#ifndef FLINT_MPARAM_H
#define FLINT_MPARAM_H

#define FLINT_FFT_SMALL_MUL_THRESHOLD           1000
#define FLINT_FFT_SMALL_SQR_THRESHOLD           1400

#define FLINT_FFT_MUL_THRESHOLD                32000
#define FLINT_FFT_SQR_THRESHOLD                32000

#define FFT_TAB \
   { {4, 4}, {4, 3}, {3, 2}, {2, 1}, {2, 1} }

#define MULMOD_TAB \
   { 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1 }

#define FFT_N_NUM                                 19
#define FFT_MULMOD_2EXPP1_CUTOFF                 128

/* FIXME: This tuning is for x86_64_adx with fft_small */
/* NOTE: We assume that the same cutoff is optimal for both mulhigh and mullow */
#define FLINT_MPN_MULHIGH_MULDERS_CUTOFF 50
#define FLINT_MPN_MULHIGH_MUL_CUTOFF 2000
#define FLINT_MPN_MULHIGH_K_TAB_SIZE 2048

#define FLINT_MPN_SQRHIGH_MULDERS_CUTOFF 90
#define FLINT_MPN_SQRHIGH_SQR_CUTOFF 2000
#define FLINT_MPN_SQRHIGH_K_TAB_SIZE 2048

#define FLINT_MPN_MULHIGH_K_TAB \
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 14, 14, 16, 15, 15, 18, 18, \
  18, 19, 20, 18, 22, 22, 20, 20, 26, 22, 22, 22, 24, 24, 24, 26, 25, 26, 30, 30, 28, 30, 31, 32, 32, 30, 36, 36, 36, 36, \
  38, 39, 39, 38, 39, 40, 40, 40, 44, 40, 44, 44, 40, 44, 44, 48, 44, 48, 44, 48, 48, 52, 52, 52, 44, 52, 52, 52, 52, 56, \
  60, 60, 52, 60, 60, 52, 52, 60, 64, 72, 56, 60, 72, 60, 60, 60, 76, 64, 60, 60, 72, 60, 72, 80, 72, 72, 80, 72, 68, 76, \
  88, 76, 68, 76, 72, 72, 80, 88, 72, 72, 88, 72, 80, 76, 76, 80, 80, 88, 80, 88, 84, 88, 80, 96, 80, 80, 88, 80, 88, 88, \
  80, 88, 96, 96, 88, 96, 92, 96, 96, 92, 100, 88, 96, 104, 88, 108, 96, 104, 104, 104, 112, 112, 108, 104, 104, 112, 112, 120, 104, 112, \
  120, 112, 112, 120, 124, 124, 116, 124, 108, 120, 124, 116, 120, 120, 116, 120, 124, 120, 120, 140, 120, 120, 120, 120, 144, 120, 132, 144, 136, 140, \
  144, 144, 144, 144, 144, 144, 144, 144, 140, 156, 140, 140, 144, 144, 144, 160, 144, 144, 156, 156, 144, 160, 160, 160, 160, 152, 160, 156, 156, 156, \
  160, 160, 144, 160, 164, 156, 156, 156, 172, 156, 156, 160, 176, 160, 160, 164, 176, 156, 160, 160, 156, 156, 160, 160, 156, 160, 172, 160, 188, 172, \
  172, 172, 160, 172, 176, 160, 160, 176, 180, 176, 164, 188, 192, 176, 172, 188, 188, 188, 172, 188, 192, 188, 180, 192, 192, 188, 188, 192, 188, 188, \
  188, 188, 192, 160, 156, 204, 160, 164, 164, 164, 164, 176, 180, 168, 172, 184, 188, 200, 216, 188, 164, 188, 220, 188, 208, 176, 180, 188, 172, 188, \
  184, 188, 204, 208, 220, 196, 220, 196, 208, 212, 188, 220, 176, 176, 184, 192, 208, 184, 188, 196, 204, 244, 208, 212, 212, 228, 256, 188, 204, 196, \
  188, 192, 192, 192, 212, 188, 292, 212, 220, 236, 228, 248, 260, 224, 264, 196, 200, 196, 212, 208, 204, 216, 208, 228, 216, 220, 252, 220, 268, 264, \
  284, 268, 300, 220, 208, 212, 220, 236, 244, 224, 252, 252, 260, 264, 256, 256, 292, 272, 288, 292, 328, 224, 256, 236, 252, 268, 256, 252, 260, 272, \
  284, 296, 300, 280, 300, 284, 252, 236, 328, 324, 264, 264, 256, 264, 280, 268, 284, 284, 292, 304, 260, 304, 264, 256, 328, 328, 260, 276, 328, 284, \
  276, 296, 300, 320, 320, 304, 328, 304, 272, 268, 280, 268, 288, 292, 288, 284, 316, 288, 328, 328, 300, 328, 328, 280, 264, 328, 300, 328, 316, 324, \
  300, 324, 300, 324, 316, 316, 328, 348, 276, 376, 288, 296, 296, 304, 320, 316, 328, 328, 324, 328, 340, 384, 348, 376, 300, 396, 304, 300, 304, 324, \
  300, 324, 328, 328, 328, 440, 448, 384, 376, 456, 464, 384, 376, 472, 480, 376, 352, 328, 376, 352, 376, 392, 392, 384, 456, 456, 480, 448, 456, 456, \
  472, 472, 472, 352, 464, 472, 472, 472, 480, 440, 480, 480, 480, 480, 456, 472, 472, 464, 464, 464, 456, 472, 480, 472, 480, 480, 480, 480, 448, 456, \
  480, 448, 456, 464, 456, 464, 456, 480, 472, 464, 464, 472, 472, 472, 480, 472, 480, 480, 472, 480, 480, 480, 480, 464, 464, 464, 456, 472, 464, 480, \
  472, 472, 480, 472, 480, 480, 464, 464, 472, 464, 472, 472, 480, 464, 480, 472, 480, 480, 576, 576, 560, 480, 472, 480, 568, 480, 480, 464, 480, 472, \
  480, 576, 480, 552, 560, 560, 560, 560, 568, 560, 560, 576, 576, 560, 568, 472, 480, 480, 544, 568, 552, 544, 560, 544, 560, 568, 552, 576, 568, 560, \
  576, 576, 568, 576, 560, 576, 568, 536, 576, 568, 560, 544, 560, 552, 560, 568, 560, 576, 568, 560, 560, 560, 568, 576, 568, 576, 576, 576, 576, 544, \
  576, 576, 568, 576, 560, 576, 576, 576, 544, 552, 568, 576, 552, 560, 576, 560, 568, 560, 576, 560, 544, 576, 576, 576, 576, 568, 576, 568, 560, 576, \
  552, 552, 576, 560, 568, 568, 568, 576, 576, 576, 560, 552, 576, 560, 568, 560, 576, 560, 568, 560, 568, 568, 568, 576, 552, 576, 560, 576, 576, 560, \
  568, 576, 568, 576, 576, 576, 576, 560, 576, 568, 568, 568, 560, 560, 576, 576, 568, 568, 576, 560, 576, 576, 568, 576, 560, 576, 576, 568, 576, 568, \
  576, 568, 576, 576, 568, 576, 576, 576, 576, 568, 576, 576, 568, 568, 576, 576, 784, 576, 776, 576, 568, 576, 576, 576, 576, 576, 576, 576, 776, 776, \
  776, 776, 776, 776, 776, 784, 776, 776, 784, 776, 776, 776, 800, 776, 776, 776, 776, 776, 776, 800, 776, 808, 792, 800, 776, 792, 776, 776, 776, 776, \
  792, 776, 776, 784, 792, 784, 800, 776, 784, 808, 784, 776, 776, 776, 808, 784, 792, 776, 792, 832, 800, 800, 816, 792, 816, 816, 856, 808, 848, 824, \
  870, 832, 792, 776, 784, 784, 784, 784, 800, 792, 800, 792, 784, 800, 800, 800, 816, 824, 824, 824, 832, 816, 816, 832, 824, 824, 848, 832, 856, 856, \
  840, 872, 864, 872, 872, 880, 880, 880, 872, 888, 880, 880, 872, 880, 880, 880, 840, 872, 872, 848, 880, 848, 856, 840, 848, 840, 880, 872, 856, 872, \
  856, 888, 880, 872, 888, 880, 872, 920, 888, 872, 880, 872, 888, 888, 888, 880, 880, 928, 880, 928, 928, 928, 920, 920, 904, 912, 880, 904, 928, 872, \
  872, 880, 888, 880, 896, 880, 872, 896, 888, 896, 896, 928, 904, 896, 896, 912, 904, 904, 920, 880, 912, 920, 928, 928, 880, 920, 920, 880, 888, 904, \
  896, 904, 928, 896, 912, 896, 912, 920, 912, 912, 928, 928, 928, 920, 928, 928, 928, 928, 928, 928, 912, 904, 912, 896, 904, 904, 920, 920, 920, 928, \
  928, 920, 928, 928, 928, 912, 928, 912, 928, 928, 928, 912, 912, 912, 928, 928, 928, 896, 928, 928, 912, 928, 928, 928, 912, 928, 912, 928, 928, 912, \
  928, 912, 928, 928, 928, 928, 928, 912, 928, 928, 928, 928, 912, 912, 928, 912, 928, 1024, 928, 928, 928, 928, 928, 928, 1056, 912, 928, 928, 1024, 1024, \
  928, 928, 1024, 928, 928, 928, 928, 928, 928, 1040, 1040, 928, 1056, 1024, 1072, 1024, 1040, 1040, 1040, 1024, 1088, 1056, 1056, 1088, 1040, 1056, 1072, 1072, 1056, 1056, \
  1024, 1088, 1040, 1024, 1040, 1040, 1024, 1056, 1056, 1056, 1040, 1072, 1056, 1040, 1056, 1056, 1056, 1056, 1056, 1056, 1056, 1120, 1056, 1088, 1056, 1120, 1088, 1072, 1104, 1104, \
  1104, 1120, 1088, 1088, 1072, 1088, 1120, 1104, 1088, 1104, 1088, 1072, 1104, 1088, 1120, 1088, 1072, 1072, 1072, 1088, 1088, 1072, 1072, 1088, 1104, 1152, 1104, 1104, 1088, 1104, \
  1136, 1088, 1104, 1152, 1152, 1152, 1136, 1120, 1136, 1152, 1120, 1152, 1088, 1120, 1104, 1120, 1136, 1104, 1136, 1088, 1136, 1104, 1088, 1104, 1120, 1104, 1104, 1120, 1136, 1136, \
  1120, 1136, 1136, 1136, 1120, 1152, 1152, 1136, 1152, 1136, 1152, 1152, 1104, 1152, 1120, 1136, 1120, 1120, 1152, 1120, 1136, 1152, 1136, 1152, 1120, 1152, 1136, 1136, 1136, 1136, \
  1136, 1152, 1152, 1152, 1152, 1120, 1120, 1152, 1136, 1136, 1136, 1152, 1152, 1120, 1152, 1152, 1152, 1152, 1152, 1104, 1152, 1152, 1120, 1136, 1152, 1120, 1152, 1136, 1152, 1152, \
  1152, 1152, 1152, 1152, 1136, 1152, 1136, 1152, 1152, 1136, 1120, 1136, 1152, 1152, 1152, 1136, 1152, 1152, 1136, 1136, 1152, 1136, 1152, 1136, 1136, 1152, 1152, 1152, 1152, 1152, \
  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1152, 1152, 1152, 1136, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, \
  1136, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1152, 1152, \
  1152, 1152, 1152, 1152, 1152, 1136, 1152, 1136, 1152, 1152, 1152, 1136, 1152, 1152, 1136, 1152, 1152, 1136, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, \
  1136, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, \
  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1152, \
  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, \
  1152, 1344, 1152, 1344, 1344, 1344, 1152, 1152, 1152, 1152, 1344, 1328, 1328, 1328, 1152, 1344, 1152, 1344, 1152, 1344, 1152, 1152, 1328, 1152, 1328, 1344, 1328, 1344, 1328, 1344, \
  1344, 1312, 1328, 1328, 1328, 1344, 1344, 1344, 1328, 1344, 1328, 1344, 1344, 1344, 1344, 1344, 1344, 1328, 1344, 1344, 1344, 1344, 1344, 1344, 1344, 1344, 1344, 1344, 1328, 1344, \
  1328, 1344, 1344, 1344, 1344, 1328, 1344, 1537, 1538, 1539, 1540, 1541, 1542, 1543, 1544, 1545, 1546, 1547, 1548, 1549, 1550, 1551, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, \
  1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1573, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1583, 1568, 1568, 1552, 1568, 1552, 1568, \
  1568, 1568, 1584, 1568, 1584, 1568, 1568, 1568, 1552, 1552, 1552, 1568, 1552, 1584, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1584, 1552, 1584, 1552, 1617, 1568, 1584, \
  1552, 1552, 1584, 1584, 1552, 1625, 1626, 1627, 1628, 1629, 1630, 1631, 1632, 1633, 1634, 1632, 1636, 1632, 1638, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, \
  1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1648, 1669, 1648, 1671, 1632, 1664, 1648, 1632, 1648, 1648, 1648, 1632, \
  1632, 1632, 1632, 1632, 1648, 1664, 1632, 1632, 1632, 1648, 1664, 1632, 1632, 1664, 1632, 1632, 1632, 1632, 1632, 1632, 1680, 1664, 1648, 1664, 1648, 1648, 1648, 1648, 1648, 1648, \
  1680, 1664, 1680, 1696, 1680, 1680, 1664, 1680, 1680, 1648, 1632, 1680, 1696, 1632, 1648, 1648, 1632, 1680, 1680, 1664, 1664, 1664, 1648, 1680, 1664, 1680, 1664, 1680, 1664, 1664, \
  1680, 1696, 1664, 1696, 1712, 1712, 1696, 1680, 1712, 1696, 1728, 1712, 1696, 1728, 1728, 1712, 1728, 1648, 1680, 1696, 1712, 1696, 1712, 1696, 1696, 1680, 1696, 1696, 1696, 1712, \
  1696, 1696, 1696, 1696, 1712, 1728, 1696, 1728, 1696, 1696, 1712, 1728, 1712, 1728, 1712, 1680, 1696, 1728, 1712, 1696, 1696, 1696, 1712, 1712, 1728, 1696, 1728, 1712, 1712, 1728, \
  1696, 1696, 1696, 1712, 1696, 1728, 1712, 1712, 1712, 1728, 1696, 1712, 1728, 1728, 1696, 1728, 1728, 1728, 1728, 1728, 1680, 1712, 1728, 1696, 1728, 1728, 1728, 1728, 1696, 1728, \
  1712, 1728, 1728, 1728, 1728, 1728, 1728, 1728, 1728, 1712, 1712, 1728, 1728, 1728, 1728, 1696, 1728, 1712, 1712, 1712, 1712, 1728, 1712, 1712, 1712, 1728, 1712, 1728, 1728, 1728, \
  1728, 1712, 1728, 1728, 1712, 1728, 1728, 1712, 1712, 1728, 1712, 1712, 1728, 1728, 1712, 1728, 1712, 1728, 1712, 1712, 1728, 1728, 1728, 1712, 1728, 1728, 1728, 1728, 1728, 1712, \
  1728, 1712, 1712, 1712, 1728, 1728, 1728, 1728, 1728, 1728, 1728, 1728, 1728, 1728, 1712, 1712, 1728, 1728, 1728, 1712, 1728, 1728, 1728, 1728, 1728, 1728, 1728, 1728, 1728, 1728, \
  1728, 1728, 1712, 1728, 1728, 1824, 1728, 1728, 1728, 1728, 1824, 1728, 1840, 1728, 1728, 1728, 1728, 1856, 1728, 1840, 1856, 1856, 1840, 1728, 1728, 1856, 1856, 1856, 1856, 1728, \
  1728, 1856, 1856, 1728, 1856, 1728, 1856, 1728, 1840, 1856, 1856, 1840, 1856, 1856, 1856, 1840, 1856, 1856, 1856, 1856, 1856, 1856, 1856, 1840, 1904, 1856, 1856, 1840, 1840, 1856, \
  1856, 1840, 1856, 1840, 1856, 1856, 1856, 1856, 1952, 1856, 1856, 1856, 1856, 1952, 1904, 1904, 1856, 1856, 1856, 1920, 1952, 2001, 1952, 1984, 1952, 1936, 1952, 1904, 1968, 1920, \
  1984, 1920, 1968, 1920, 1936, 1856, 2000, 1920, 1936, 1952, 2000, 1968, 1984, 1968, 1984, 2000, 1952, 2000, 2016, 1984, 2000, 2016, 1984, 1664, 2016, 1984, 2016, 1968, 2016, 2016, \
  1744, 2016, 2016, 1968, 2000, 1728, 1712, 1696

#define FLINT_MPN_SQRHIGH_K_TAB \
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 28, 29, 0, 0, 31, 0, 31, \
  0, 32, 0, 34, 0, 36, 36, 40, 0, 40, 40, 40, 0, 44, 0, 44, 44, 48, 0, 52, 48, 52, 44, 44, 48, 48, 48, 48, 48, 48, \
  48, 48, 48, 48, 52, 52, 52, 52, 52, 56, 56, 52, 56, 56, 56, 56, 56, 56, 56, 56, 60, 60, 68, 64, 60, 60, 60, 60, 64, 64, \
  64, 68, 68, 68, 64, 64, 68, 68, 72, 72, 68, 68, 68, 72, 76, 80, 72, 72, 72, 88, 76, 76, 80, 80, 76, 80, 76, 80, 84, 80, \
  88, 88, 80, 84, 88, 80, 80, 84, 92, 88, 92, 88, 88, 88, 88, 96, 88, 108, 100, 92, 88, 88, 104, 100, 100, 92, 104, 108, 100, 92, \
  104, 100, 104, 96, 108, 104, 96, 96, 104, 100, 100, 104, 112, 116, 108, 104, 104, 116, 108, 104, 104, 120, 116, 104, 108, 132, 116, 108, 120, 108, \
  108, 108, 132, 108, 120, 112, 112, 116, 132, 128, 116, 124, 128, 116, 132, 120, 132, 120, 124, 120, 132, 120, 124, 128, 120, 128, 128, 132, 144, 124, \
  128, 124, 140, 128, 128, 124, 136, 132, 128, 128, 140, 144, 128, 128, 140, 136, 132, 144, 148, 152, 144, 132, 160, 156, 140, 144, 156, 144, 140, 144, \
  140, 156, 156, 156, 140, 144, 168, 156, 156, 164, 168, 156, 156, 160, 144, 144, 180, 156, 152, 168, 156, 160, 156, 148, 180, 168, 180, 156, 164, 156, \
  156, 172, 156, 156, 156, 180, 180, 172, 180, 168, 164, 172, 164, 172, 176, 176, 168, 176, 172, 176, 180, 168, 176, 180, 180, 180, 192, 184, 180, 176, \
  204, 176, 188, 180, 188, 180, 204, 180, 192, 204, 180, 192, 180, 204, 228, 192, 188, 192, 204, 180, 192, 216, 200, 216, 228, 216, 204, 216, 188, 216, \
  216, 204, 216, 192, 204, 212, 228, 204, 228, 216, 204, 216, 192, 204, 204, 212, 204, 216, 204, 216, 228, 228, 224, 212, 204, 216, 216, 224, 228, 212, \
  216, 228, 212, 204, 204, 216, 216, 216, 228, 216, 220, 224, 228, 220, 228, 228, 236, 224, 260, 224, 228, 228, 260, 228, 248, 252, 264, 248, 216, 212, \
  212, 212, 220, 216, 216, 220, 216, 224, 216, 228, 232, 224, 220, 224, 240, 244, 236, 244, 232, 256, 288, 240, 236, 288, 260, 260, 264, 232, 228, 228, \
  228, 236, 236, 228, 236, 236, 236, 248, 256, 260, 232, 236, 264, 256, 260, 252, 284, 252, 264, 276, 244, 296, 244, 240, 244, 248, 240, 252, 256, 260, \
  260, 252, 252, 272, 252, 272, 260, 296, 268, 260, 256, 276, 272, 264, 312, 284, 256, 252, 252, 252, 264, 252, 256, 260, 264, 276, 292, 268, 264, 276, \
  264, 260, 264, 304, 288, 296, 296, 296, 288, 296, 272, 280, 264, 264, 272, 264, 288, 288, 280, 288, 272, 280, 296, 296, 280, 288, 296, 320, 344, 320, \
  344, 272, 344, 344, 304, 288, 280, 280, 312, 280, 280, 304, 304, 312, 304, 296, 288, 328, 320, 352, 320, 320, 328, 360, 344, 344, 360, 288, 288, 304, \
  288, 296, 320, 320, 312, 312, 304, 328, 336, 312, 312, 360, 336, 344, 344, 336, 360, 360, 296, 296, 304, 360, 328, 328, 312, 320, 320, 328, 312, 344, \
  344, 328, 344, 344, 368, 360, 352, 360, 392, 392, 368, 320, 312, 392, 312, 328, 344, 336, 344, 328, 360, 352, 352, 360, 360, 360, 368, 408, 360, 376, \
  392, 392, 376, 336, 344, 352, 360, 352, 344, 344, 384, 344, 360, 376, 392, 368, 360, 408, 408, 448, 432, 384, 392, 336, 360, 344, 360, 360, 368, 376, \
  376, 360, 368, 408, 368, 376, 376, 376, 432, 376, 384, 464, 432, 344, 432, 376, 344, 344, 384, 344, 384, 376, 400, 432, 456, 432, 456, 392, 432, 392, \
  448, 456, 360, 376, 456, 408, 384, 368, 376, 432, 376, 472, 464, 504, 448, 360, 408, 456, 376, 408, 424, 424, 448, 440, 456, 392, 408, 384, 408, 392, \
  384, 416, 424, 432, 400, 472, 480, 408, 432, 432, 464, 456, 504, 464, 456, 472, 496, 416, 392, 424, 504, 400, 440, 432, 472, 448, 456, 456, 432, 448, \
  456, 504, 504, 512, 496, 504, 424, 416, 408, 432, 544, 432, 440, 456, 448, 448, 464, 472, 480, 552, 552, 544, 552, 560, 552, 544, 544, 544, 552, 552, \
  552, 552, 552, 552, 552, 560, 552, 552, 560, 560, 560, 560, 568, 576, 568, 592, 576, 560, 560, 560, 560, 568, 568, 576, 568, 568, 576, 568, 568, 576, \
  568, 568, 600, 592, 576, 608, 576, 576, 576, 576, 584, 584, 576, 592, 592, 600, 584, 592, 584, 600, 600, 600, 584, 624, 584, 592, 600, 592, 616, 592, \
  592, 592, 592, 592, 592, 592, 600, 608, 600, 608, 600, 600, 600, 624, 616, 608, 616, 600, 632, 608, 608, 608, 616, 608, 608, 608, 608, 616, 608, 664, \
  664, 632, 656, 616, 664, 624, 616, 632, 632, 616, 624, 624, 624, 632, 624, 624, 624, 624, 648, 632, 624, 624, 872, 872, 656, 664, 696, 872, 872, 872, \
  632, 872, 872, 896, 872, 896, 872, 896, 880, 872, 896, 872, 888, 872, 872, 872, 872, 872, 888, 872, 880, 824, 872, 856, 880, 888, 800, 848, 800, 880, \
  848, 800, 808, 872, 872, 864, 824, 840, 872, 872, 872, 872, 872, 872, 872, 872, 872, 872, 880, 928, 872, 920, 872, 920, 888, 880, 872, 880, 872, 872, \
  872, 888, 880, 888, 888, 872, 888, 872, 880, 896, 920, 920, 880, 928, 904, 872, 888, 904, 872, 872, 872, 872, 880, 888, 872, 872, 872, 880, 872, 904, \
  920, 888, 888, 872, 880, 880, 888, 896, 896, 880, 872, 904, 880, 912, 896, 872, 904, 904, 912, 880, 880, 920, 912, 896, 928, 928, 872, 872, 872, 872, \
  880, 880, 888, 888, 888, 912, 880, 880, 880, 912, 896, 896, 928, 896, 912, 912, 896, 928, 928, 928, 928, 912, 880, 928, 880, 880, 896, 880, 896, 880, \
  896, 880, 880, 912, 912, 880, 896, 880, 896, 912, 880, 880, 928, 912, 912, 928, 928, 896, 928, 912, 928, 928, 928, 928, 928, 928, 880, 880, 880, 880, \
  880, 912, 896, 912, 912, 896, 928, 896, 928, 896, 928, 928, 928, 928, 912, 880, 928, 912, 896, 928, 896, 880, 912, 896, 912, 896, 896, 896, 928, 896, \
  912, 928, 912, 928, 912, 928, 912, 880, 928, 928, 880, 928, 928, 880, 896, 928, 880, 896, 912, 928, 928, 896, 912, 928, 912, 928, 1024, 928, 928, 1040, \
  912, 928, 1024, 1024, 928, 928, 928, 928, 912, 1040, 1072, 912, 1072, 1072, 1056, 1040, 1072, 1088, 1104, 1040, 1040, 928, 1024, 928, 928, 1056, 928, 928, 1056, 1040, \
  928, 1056, 1024, 1056, 1056, 1056, 1104, 1024, 1056, 1072, 1056, 1056, 1088, 1072, 1120, 1104, 1072, 1072, 1072, 1104, 1120, 1024, 1120, 1120, 1056, 1024, 1040, 1024, 1040, 1040, \
  1040, 1024, 1040, 1056, 1040, 1088, 1040, 1056, 1120, 1040, 1040, 1072, 1072, 1088, 1104, 1088, 1104, 1088, 1088, 1056, 1120, 1104, 1152, 1104, 1088, 1120, 1104, 1120, 1024, 1120, \
  1152, 1136, 1152, 1056, 1056, 1072, 1056, 1152, 1088, 1056, 1104, 1072, 1104, 1088, 1120, 1072, 1104, 1088, 1088, 1120, 1120, 1120, 1152, 1136, 1120, 1120, 1152, 1120, 1136, 1152, \
  1152, 1072, 1136, 1136, 1088, 1152, 1136, 1136, 1104, 1056, 1072, 1088, 1104, 1104, 1120, 1120, 1104, 1104, 1104, 1136, 1152, 1136, 1152, 1136, 1136, 1152, 1152, 1152, 1136, 1072, \
  1152, 1136, 1152, 1152, 1104, 1152, 1104, 1152, 1152, 1088, 1120, 1152, 1136, 1136, 1120, 1136, 1120, 1152, 1136, 1136, 1120, 1136, 1152, 1136, 1152, 1136, 1136, 1120, 1152, 1152, \
  1104, 1152, 1152, 1152, 1104, 1120, 1104, 1120, 1120, 1120, 1120, 1136, 1136, 1136, 1136, 1152, 1120, 1120, 1152, 1152, 1152, 1136, 1152, 1136, 1136, 1152, 1120, 1152, 1152, 1104, \
  1120, 1152, 1136, 1136, 1136, 1136, 1136, 1136, 1120, 1152, 1152, 1152, 1136, 1104, 1152, 1136, 1136, 1152, 1136, 1120, 1136, 1152, 1152, 1136, 1136, 1120, 1136, 1136, 1136, 1120, \
  1152, 1152, 1152, 1136, 1136, 1136, 1152, 1104, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1136, 1136, 1136, 1136, 1152, 1136, 1136, 1152, 1136, 1152, 1152, 1136, 1120, 1152, \
  1136, 1136, 1152, 1136, 1152, 1152, 1152, 1152, 1136, 1136, 1136, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1152, 1136, 1136, 1136, 1136, 1152, 1152, 1152, 1152, 1152, \
  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1152, 1136, 1152, 1136, 1136, 1152, 1120, 1152, 1152, 1152, 1152, 1136, 1136, 1152, 1136, 1152, \
  1120, 1152, 1136, 1136, 1152, 1152, 1152, 1152, 1152, 1152, 1136, 1152, 1136, 1152, 1152, 1136, 1152, 1152, 1152, 1136, 1136, 1152, 1136, 1152, 1152, 1152, 1152, 1152, 1152, 1152, \
  1136, 1152, 1152, 1152, 1152, 1152, 1136, 1136, 1152, 1136, 1152, 1120, 1152, 1152, 1152, 1136, 1152, 1136, 1152, 1152, 1136, 1152, 1136, 1136, 1136, 1136, 1152, 1152, 1152, 1136, \
  1152, 1152, 1152, 1104, 1120, 1152, 1136, 1152, 1152, 1136, 1152, 1136, 1152, 1152, 1152, 1136, 1136, 1152, 1152, 1120, 1152, 1136, 1152, 1152, 1120, 1152, 1120, 1136, 1152, 1152, \
  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1120, 1152, 1152, 1136, 1152, 1136, 1552, 1552, 1136, 1136, 1152, 1552, 1136, 1136, 1584, 1552, 1552, 1552, \
  1152, 1552, 1552, 1152, 1584, 1552, 1152, 1152, 1552, 1584, 1152, 1152, 1552, 1584, 1568, 1552, 1552, 1552, 1152, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1552, \
  1552, 1552, 1552, 1552, 1552, 1625, 1552, 1552, 1552, 1552, 1552, 1552, 1552, 1568, 1632, 1552, 1632, 1568, 1584, 1568, 1568, 1584, 1600, 1552, 1552, 1632, 1600, 1632, 1632, 1632, \
  1632, 1568, 1584, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1552, 1662, 1632, 1632, 1632, 1648, 1632, 1632, 1632, 1568, 1600, 1632, 1552, 1616, 1632, 1632, 1632, 1584, 1632, \
  1632, 1552, 1632, 1632, 1632, 1568, 1600, 1632, 1648, 1616, 1632, 1648, 1632, 1632, 1632, 1632, 1648, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1680, 1664, 1632, 1664, 1632, 1632, \
  1632, 1632, 1632, 1632, 1632, 1648, 1648, 1648, 1616, 1632, 1680, 1632, 1680, 1648, 1632, 1664, 1632, 1648, 1632, 1632, 1632, 1632, 1632, 1664, 1632, 1632, 1632, 1632, 1632, 1632, \
  1648, 1648, 1632, 1632, 1632, 1632, 1632, 1632, 1632, 1664, 1632, 1648, 1664, 1648, 1648, 1632, 1648, 1632, 1680, 1664, 1648, 1680, 1664, 1696, 1632, 1664, 1680, 1632, 1648, 1696, \
  1680, 1632, 1632, 1648, 1648, 1648, 1632, 1664, 1632, 1648, 1648, 1632, 1648, 1648, 1632, 1664, 1664, 1680, 1680, 1632, 1664, 1664, 1648, 1648, 1712, 1728, 1664, 1712, 1696, 1712, \
  1680, 1632, 1728, 1648, 1728, 1728, 1712, 1648, 1632, 1632, 1664, 1648, 1664, 1664, 1680, 1712, 1664, 1632, 1696, 1664, 1712, 1696, 1680, 1680, 1696, 1680, 1680, 1712, 1648, 1712, \
  1632, 1728, 1696, 1648, 1712, 1632, 1712, 1696, 1648, 1680, 1648, 1664, 1648, 1696, 1712, 1648, 1648, 1680, 1664, 1696, 1728, 1696, 1712, 1728, 1696, 1728, 1712, 1664, 1680, 1728, \
  1712, 1712, 1712, 1728, 1728, 1712, 1728, 1696, 1728, 1712, 1696, 1728, 1728, 1680, 1728, 1712, 1664, 1680, 1680, 1696, 1712, 1696, 1696, 1696, 1728, 1696, 1728, 1728, 1696, 1712, \
  1712, 1664, 1712, 1680, 1664, 1728, 1728, 1664, 1696, 1696, 1680, 1712, 1680, 1712, 1696, 1728, 1696, 1696, 1728, 1696, 1712, 1696, 1712, 1712, 1712, 1728, 1712, 1696, 1728, 1680, \
  1696, 1712, 1712, 1728, 1712, 1728, 1696, 1728, 1776, 1728, 1696, 1776, 1728, 1728, 1712, 1712, 1824, 1856, 1728, 1824, 1728, 1728, 1712, 1728, 1728, 1728, 1856, 1696, 1728, 1840, \
  1712, 1712, 1824, 1856, 1792, 1712, 1840, 1728, 1808, 1728, 1824, 1824, 1840, 1824, 1824, 1856, 1856, 1824, 1856, 1840, 1696, 1856, 1840, 1840, 1856, 1856, 1824, 1712, 1792, 1856, \
  1824, 1728, 1808, 1792, 1856, 1728, 1792, 1840, 1808, 1808, 1712, 1808, 1840, 1808, 1824, 1824, 1824, 1840, 1824, 1840, 1856, 1824, 1728, 1856, 1856, 1824, 1856, 1856, 1792, 1792, \
  1856, 1824, 1856, 1824, 1824, 1856, 1808, 1824, 1856, 1856, 1840, 1840, 1840, 1840, 1856, 1840, 1840, 1856, 1824, 1840, 1808, 1824, 1840, 1856, 1856, 1824, 1856, 1840, 1840, 1840, \
  1824, 1824, 1840, 1840, 1840, 1856, 1856, 1856

#endif
