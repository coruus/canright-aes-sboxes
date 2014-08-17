/* bestboth.c
 *
 * by: David Canright
 *
 * for each input basis, and each of 4 transformation matrices,
 * takes bit matrix and finds equivalent with minimum # of gates
 *  combining both input matrices, and both output matrices
 * NOTE: matrix input order is: [A2X, X2A, X2S, S2X]
 *
 * input should have lines of the form:
 hexstring num
 * where hexstring contains all 4 matrices, num is an ID#, e.g.:
98F3F2480981A9FF64786E8C6829DE60582D9E0BDC0403248C7905EB12045153  4
 * for which the output should be:
basis #  4:
   A2X: 98F3F2480981A9FF   S2X: 8C7905EB12045153
 ncols =  8, gates = 42
  A2Xb: 0000000000012804100810224008808001
  S2Xb: 0028006200000100008800000102044010
 [0,2],  [0,3],  [1,7],  [2,10],  [3,11],  [4,7],  [5,8],  [6,10],  [4,15],
 ncols = 17, gates = 20
   X2S: 582D9E0BDC040324   X2A: 64786E8C6829DE60
 ncols =  8, gates = 38
  X2Sb: 000000000000000040082480180002040100
  X2Ab: 041000800021D00000000000000204080860
 [0,4],  [1,3],  [1,7],  [2,4],  [2,8],  [2,6],  [3,13],  [5,11],  [6,9],  [10,12],
 ncols = 18, gates = 18
***bestgates   4 =    38   =   20 +   18
 * which, for each matrix pair, shows the original versions (8 columns),
 *   the optimized versions, and a list of index pairs for precomputed XORs,
 *   which correspond to new columns. Also shown: # XOR gates required.
 * Note: a "quick" test case is:
F1261450CA86D330C502A8BF412B3590352582D03974323C65C4836C69953380    0
 *
 * uses pruning algorithm to eliminate redundant cases; minimal memory copying
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 8

/* gatematrix is a structure with an array of 16-bit columns,
 list of indices (used in pairs), number of columns, and number of gates*/
typedef struct gatematrix {
  uint32_t mat[128];
  char ind[256];
  int n;
  int g;
} GateMat;

static uint32_t share[65536];
static GateMat test;

/* blockPrint prints columns and index pairs for matrix pair */
void blockPrint(GateMat* p, const char* tag1, const char* tag2) {
  int i;

  printf("%6s: ", tag1);
  for (i = 0; i < p->n; i++) {
    printf("%02X", (p->mat[i]) & 0XFF);
  }
  if ((p->n) > N) {
    printf("\n");
  }
  printf("%6s: ", tag2);
  for (i = 0; i < p->n; i++) {
    printf("%02X", ((p->mat[i]) & 0XFF00) >> 8);
  }
  if ((p->n) > N)
    printf("\n");
  for (i = 0; i < (p->n) - N; i++) {
    printf(" [%1d,%1d], ", p->ind[2 * i], p->ind[2 * i + 1]);
  }
  printf("\n ncols = %2d, gates = %2d\n", p->n, p->g);

} /* end blockPrint */

/* copyMat copies from one to another*/
void copyMat(GateMat* p, GateMat* q) {
  int i, n;

  n = q->n = p->n;
  q->g = p->g;
  memcpy(q->mat, p->mat, n * sizeof(uint32_t));
  memcpy(q->ind, p->ind, (n - N) * 2);
} /* end copyMat */

/*
* bestgates is recursive:
*   takes current matrix, tries all possibilities of adding a gate
*   returns best # of gates
* p points to test matrix on input, and used to store output.
* tree search is pruned if this set of columns previously tried
*/
void bestgates() {
  char indb[256];
  int gb, nb, ci, cj;
  int i, j, n, c, g, io, jo;
  int nm, np, n2, n2p, t;

  gb = 1024; /* best # gates, start high */
  n = test.n;
  g = test.g;
  nm = n - 1;
  np = n + 1;
  n2 = 2 * (n - N);
  n2p = n2 + 1;
  if (n == N) {
    io = jo = 0; /* if orig matrix, no "old" index pair */
  } else {
    io = (test.ind[n2 - 2]);
    jo = (test.ind[n2 - 1]);
  }
  for (i = 0; i < nm; i++) /* for each pair of columns */
    for (j = i + 1; j < n; j++) {
      c = (test.mat[i]) & (test.mat[j]);
      t = share[c];
      if (t) { /* if can share a gate ????????????? */
        if (i < io && j != io && j != jo && j < nm) {
          // if prior, indep. pair, then been there, done that; skip to next j
          continue;
        }
        test.n = np;
        test.g = g - t;
        ci = test.mat[i]; /* save current columns */
        cj = test.mat[j];
        test.mat[i] ^= c; /* update to new columns */
        test.mat[j] ^= c;
        test.mat[n] = c;
        test.ind[n2] = i;
        test.ind[n2p] = j;
        bestgates();      /* recurse with new matrix */
        test.mat[i] = ci; /* restore current columns */
        test.mat[j] = cj;
        if (test.g < gb) { /* if best yet, save data */
          memcpy(indb, test.ind + n2, (test.n - n) * 2);
          nb = test.n;
          gb = test.g;
        }
      }
    }              /* end columns loop */
  if (gb < 1024) { /* if improved, return best data */
    memcpy(test.ind + n2, indb, (nb - n) * 2);
    test.n = nb;
    test.g = gb;
  }
  /*      else {printf("%3d [%2d]",n,g); fflush(stdout);} */
}

/* bestmat reconstructs best matrix */
void bestmat(GateMat* p) {
  int i, j, n, c;
  int nm, np, n2, n2p, t;
  GateMat best;

  p->g = test.g;
  for (i = 0; i < N; i++) test.mat[i] = p->mat[i];
  for (n = 0; n < (test.n - N); n++) {
    i = test.ind[n * 2];
    j = test.ind[n * 2 + 1];
    c = (test.mat[i]) & (test.mat[j]);
    test.mat[i] ^= c;
    test.mat[j] ^= c;
    test.mat[n + N] = c;
  }
}

int main(void) {
  char line[256];
  char name[4][4] = {
      "A2X", "X2S", "S2X", "X2A",
  };
  char bname[4][5] = {
      "A2Xb", "X2Sb", "S2Xb", "X2Ab",
  };
  long int i, j, k, n, nid, gt;
  unsigned u;
  int InitMat[32];
  GateMat orig[2];

  /* share[i] is initialized to 0 if # bits < 2 */
  share[0] = 0;
  for (i = 1; i < 65536; i++) {
    k = 0;
    for (j = i & 0xFFFF; j; j >>= 1) k += j & 1;
    share[i] = k - 1;
  }

  while (fgets(line, 256, stdin) == line) {
    for (i = 0; i < 32; i++) { /* read matrices, ID number */
      sscanf(line + 2 * i, "%02X", &u);
      InitMat[i] = u;
    }
    sscanf(line + 65, "%ld", &nid);
    printf("\nbasis #%3ld:\n", nid);

    /* NOTE: matrix input order is: [A2X, X2A, X2S, S2X] */
    for (i = 0; i < 8; i++) { /* combine input pair; combine output pair */
      (orig[0]).mat[i] = InitMat[8 * 0 + i] | (InitMat[8 * 3 + i] << 8);
      (orig[1]).mat[i] = InitMat[8 * 2 + i] | (InitMat[8 * 1 + i] << 8);
    }

    gt = 0;
    for (k = 0; k < 2; k++) { /* for each matrix pair */
      (orig[k]).n = 8;        /* initialize # columns, # gates */
      for (i = j = 0; i < 8; i++) {
        j += share[(orig[k]).mat[i]];
      }
      (orig[k]).g = j - 8;
      blockPrint(&(orig[k]), name[k], name[k + 2]);
      fflush(stdout);

      copyMat(&(orig[k]), &test);
      bestgates(); /* optimize */
      bestmat(&(orig[k]));

      blockPrint(&test, bname[k], bname[k + 2]);
      fflush(stdout);
      gt += test.g; /* total # gates */
    }
    printf("***bestgates %3ld = %5ld   =%5d +%5d\n", nid, gt, (orig[0]).g, (orig[1]).g);
    fflush(stdout);
  }

  return 0;
}
