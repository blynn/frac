#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"
#include "test.h"

static void *count_int_fn(cf_t cf) {
  int n = 1;
  while(cf_wait(cf)) {
    cf_put_int(cf, n);
    n++;
  }
  return NULL;
}

int main() {
  mpz_t z;
  mpz_init(z);
  cf_t x = cf_new_const(count_int_fn);
  cf_t out[2];
  int i[2];
  cf_tee(out, x);

  void get(int k, int n) {
    while(n) {
      cf_get(z, out[k]);
      EXPECT(!mpz_cmp_ui(z, i[k]++));
      n--;
    }
  }

  i[0] = 1;
  i[1] = 1;
  get(0, 5);
  get(1, 10);
  get(0, 20);

  cf_free(out[0]);
  get(1, 100);
  cf_free(out[1]);

  mpz_clear(z);
  return 0;
}
