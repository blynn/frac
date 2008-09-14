#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"

int main() {
  mpz_t z;
  mpz_init(z);
  cf_t x, conv;

  mpz_set_ui(z, 3);
  x = cf_new_tanh(z);

  printf("tanh 3 = ");
  conv = cf_new_cf_to_decimal(x);
  for (int i = 0; i <= 5000; i++) {
    cf_get(z, conv);
    gmp_printf("%Zd", z);
    if (!(i % 5)) putchar(' ');
    if (!(i % 50)) putchar('\n');
  }
  cf_free(conv);
  cf_free(x);

  mpz_set_ui(z, 2);
  x = cf_new_epow(z);

  printf("e^2 = ");
  conv = cf_new_cf_to_decimal(x);
  for (int i = 0; i <= 5000; i++) {
    cf_get(z, conv);
    gmp_printf("%Zd", z);
    if (!(i % 5)) putchar(' ');
    if (!(i % 50)) putchar('\n');
  }
  cf_free(conv);
  cf_free(x);

  mpz_clear(z);
  return 0;
}
