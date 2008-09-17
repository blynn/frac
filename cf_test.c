// Test demand channel infrastructure.

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"
#include "test.h"

static void *count_fn(cf_t cf) {
  mpz_t z;
  mpz_init(z);
  mpz_set_ui(z, 0);
  while(cf_wait(cf)) {
    cf_put(cf, z);
    mpz_add_ui(z, z, 1);
  }
  mpz_clear(z);
  return NULL;
}

static void *count_int_fn(cf_t cf) {
  int n = 0;
  while(cf_wait(cf)) {
    cf_put_int(cf, n);
    n++;
  }
  return NULL;
}

int main() {
  mpz_t z, z1;
  mpz_init(z);
  mpz_init(z1);
  cf_t a, b;
  a = cf_new_const(count_fn);
  b = cf_new_const(count_int_fn);
  for (int i = 0; i < 100; i++) {
    cf_get(z, a);
    EXPECT(!mpz_cmp_ui(z, i));
    cf_get(z, b);
    EXPECT(!mpz_cmp_ui(z, i));
  }
  cf_free(a);
  cf_free(b);
  mpz_clear(z);
  mpz_clear(z1);
  return 0;
}
