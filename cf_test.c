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
    mpz_add_ui(z, z, 1);
    cf_put(cf, z);
  }
  mpz_clear(z);
  return NULL;
}

int main() {
  mpz_t z;
  mpz_init(z);
  cf_t a;
  a = cf_new(count_fn, NULL);
  for (int i = 1; i <= 10; i++) {
    cf_get(z, a);
    EXPECT(!mpz_cmp_ui(z, i));
  }
  cf_free(a);
  mpz_clear(z);
  return 0;
}
