// Test Gosper's example: find coth 1/2 via Newton's method.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"
#include "test.h"

static void *cot1_fn(cf_t cf) {
  int n = 1;
  while(cf_wait(cf)) {
    cf_put_int(cf, n);
    n += 2;
  }
  return NULL;
}

int main() {
  mpz_t a[6];
  int i;
  for (i = 0; i < 6; i++) mpz_init(a[i]);
  cf_t x, y;

  // Gosper's example.
  mpz_set_si(a[0], 1);
  mpz_set_si(a[3], -1);
  mpz_set_si(a[5], 1);
  x = cf_new_const(cot1_fn);
  y = cf_new_newton(x, a, a[0]);
  CF_EXPECT_DEC(y, "2.16395341373865284877");
  cf_free(x);
  cf_free(y);
  
  // Confirm sqrt(1-(sin 1)^2) = cos 1
  x = cf_new_sin1();
  cf_t s1[2];
  cf_tee(s1, x);
  mpz_t b[8];
  mpz8_init(b);
  mpz8_set_int(b,
      -1, 0, 0, 1,
       0, 0, 0, 1);

  cf_t bi = cf_new_bihom(s1[0], s1[1], b);
  cf_t n = cf_new_sqrt(bi);

  CF_EXPECT_DEC(n, "0.54030230586813971740");
  cf_free(x);
  cf_free(s1[0]);
  cf_free(s1[1]);
  cf_free(bi);
  cf_free(n);

  x = cf_new_sqrt_int(355, 113);
  CF_EXPECT_DEC(x, "1.77245392615830279609");
  cf_free(x);

  mpz8_clear(b);
  for (i = 0; i < 6; i++) mpz_clear(a[i]);
  return 0;
}
