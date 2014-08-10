/* sbox.c
 *
 * by: David Canright
 *
 * illustrates compact implementation of AES S-box via subfield operations
 *   case # 4 : [d^16, d], [alpha^8, alpha^2], [Omega^2, Omega]
 *   nu = beta^8 = N^2*alpha^2, N = w^2
 */

#include <stdio.h>
#include <sys/types.h>

/* to convert between polynomial (A^7...1) basis A & normal basis X */
/* or to basis S which incorporates bit matrix of Sbox */
static int A2X[8] = {0x98, 0xF3, 0xF2, 0x48, 0x09, 0x81, 0xA9, 0xFF},
           X2A[8] = {0x64, 0x78, 0x6E, 0x8C, 0x68, 0x29, 0xDE, 0x60},
           X2S[8] = {0x58, 0x2D, 0x9E, 0x0B, 0xDC, 0x04, 0x03, 0x24},
           S2X[8] = {0x8C, 0x79, 0x05, 0xEB, 0x12, 0x04, 0x51, 0x53};

/* multiply in GF(2^2), using normal basis (Omega^2,Omega) */
int G4_mul(int x, int y) {
  int a, b, c, d, e, p, q;

  a = (x & 0x2) >> 1;
  b = (x & 0x1);
  c = (y & 0x2) >> 1;
  d = (y & 0x1);
  e = (a ^ b) & (c ^ d);
  p = (a & c) ^ e;
  q = (b & d) ^ e;
  return ((p << 1) | q);
}

/* scale by N = Omega^2 in GF(2^2), using normal basis (Omega^2,Omega) */
int G4_scl_N(int x) {
  int a, b, p, q;

  a = (x & 0x2) >> 1;
  b = (x & 0x1);
  p = b;
  q = a ^ b;
  return ((p << 1) | q);
}

/* scale by N^2 = Omega in GF(2^2), using normal basis (Omega^2,Omega) */
int G4_scl_N2(int x) {
  int a, b, p, q;

  a = (x & 0x2) >> 1;
  b = (x & 0x1);
  p = a ^ b;
  q = a;
  return ((p << 1) | q);
}

/* square in GF(2^2), using normal basis (Omega^2,Omega) */
/* NOTE: inverse is identical */
int G4_sq(int x) {
  int a, b;

  a = (x & 0x2) >> 1;
  b = (x & 0x1);
  return ((b << 1) | a);
}

/* multiply in GF(2^4), using normal basis (alpha^8,alpha^2) */
int G16_mul(int x, int y) {
  int a, b, c, d, e, p, q;

  a = (x & 0xC) >> 2;
  b = (x & 0x3);
  c = (y & 0xC) >> 2;
  d = (y & 0x3);
  e = G4_mul(a ^ b, c ^ d);
  e = G4_scl_N(e);
  p = G4_mul(a, c) ^ e;
  q = G4_mul(b, d) ^ e;
  return ((p << 2) | q);
}

/* square & scale by nu in GF(2^4)/GF(2^2), normal basis (alpha^8,alpha^2) */
/*   nu = beta^8 = N^2*alpha^2, N = w^2 */
int G16_sq_scl(int x) {
  int a, b, p, q;

  a = (x & 0xC) >> 2;
  b = (x & 0x3);
  p = G4_sq(a ^ b);
  q = G4_scl_N2(G4_sq(b));
  return ((p << 2) | q);
}

/* inverse in GF(2^4), using normal basis (alpha^8,alpha^2) */
int G16_inv(int x) {
  int a, b, c, d, e, p, q;

  a = (x & 0xC) >> 2;
  b = (x & 0x3);
  c = G4_scl_N(G4_sq(a ^ b));
  d = G4_mul(a, b);
  e = G4_sq(c ^ d);  // really inverse, but same as square
  p = G4_mul(e, b);
  q = G4_mul(e, a);
  return ((p << 2) | q);
}

/* inverse in GF(2^8), using normal basis (d^16,d) */
int G256_inv(int x) {
  int a, b, c, d, e, p, q;

  a = (x & 0xF0) >> 4;
  b = (x & 0x0F);
  c = G16_sq_scl(a ^ b);
  d = G16_mul(a, b);
  e = G16_inv(c ^ d);
  p = G16_mul(e, b);
  q = G16_mul(e, a);
  return ((p << 4) | q);
}

/* convert to new basis in GF(2^8) */
/* i.e., bit matrix multiply */
int G256_newbasis(int x, int b[]) {
  int i, y = 0;

  for (i = 7; i >= 0; i--) {
    if (x & 1)
      y ^= b[i];
    x >>= 1;
  }
  return (y);
}

/* find Sbox of n in GF(2^8) mod POLY */
int Sbox(int n) {
  int t;

  t = G256_newbasis(n, A2X);
  t = G256_inv(t);
  t = G256_newbasis(t, X2S);
  return (t ^ 0x63);
}

/* find inverse Sbox of n in GF(2^8) mod POLY */
int iSbox(int n) {
  int t;

  t = G256_newbasis(n ^ 0x63, S2X);
  t = G256_inv(t);
  t = G256_newbasis(t, X2A);
  return (t);
}

/* compute tables of Sbox & its inverse; print 'em out */
int main() {
  int Sbox_tbl[256], iSbox_tbl[256], i, j;

  for (i = 0; i < 256; i++) {
    Sbox_tbl[i] = Sbox(i);
    iSbox_tbl[i] = iSbox(i);
  }
  printf("char S[256] = {\n");
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      printf("%3d, ", Sbox_tbl[i * 16 + j]);
    }
    printf("\n");
  }
  printf("};\n\n");
  printf("char Si[256] = {\n");
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      printf("%3d, ", iSbox_tbl[i * 16 + j]);
    }
    printf("\n");
  }
  printf("};\n\n");
  return (0);
}
