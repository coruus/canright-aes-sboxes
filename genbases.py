"""Script to generate bases for Canright S-box programs."""
from __future__ import division, print_function

one = 0x01
d, d2, d4, d8                 = 0xff, 0x13, 0x1e, 0x4f
d128, d64, d32, d16           = 0x4e, 0x1f, 0x12, 0xfe
L, L2, L4, L8                 = 0xa2, 0xf2, 0x42, 0xaf
L16, L32, L64, L128           = 0xa3, 0xf3, 0x43, 0xae
alpha, alpha2, alpha4, alpha8 = 0xe1, 0x5c, 0xe0, 0x5d
omega, omega2                 = 0xbd, 0xbc

GF2_8 = [\
    [d16, d],
    [d32, d2],
    [d64, d4],
    [d128, d8],
    [L16, L],
    [L32, L2],
    [L64, L4],
    [L128, L8],
    [d, one],
    [d16, one],
    [d2, one],
    [d32, one],
    [d4, one],
    [d64, one],
    [d8, one],
    [d128, one],
    [L, one],
    [L16, one],
    [L2, one],
    [L32, one],
    [L4, one],
    [L64, one],
    [L8, one],
    [L128, one]]
GF2_4 = [\
  [alpha4, alpha],
  [alpha8, alpha2],
  [alpha, one],
  [alpha4, one],
  [alpha2, one],
  [alpha8, one]]
GF2_2 = [\
  [omega2, omega],
  [omega, one],
  [omega2, one]]
