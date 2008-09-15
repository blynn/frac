#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"

void cf_expect_dec(cf_t x, char *result, char *filename, int line) {
  mpz_t z;
  mpz_init(z);
  cf_t conv;
  conv = cf_new_cf_to_decimal(x);
  int len = strlen(result);
  char s[len + 1];
  for (int i = 0; i <= len; i++) {
    cf_get(z, conv);
    s[i] = mpz_get_ui(z) + '0';
  }
  s[len] = 0;
  if (strcmp(s, result)) {
    fprintf(stderr, "\n%s:%d: bad continued fraction decimal expansion\n",
        filename, line);
    fprintf(stderr, "  expected: %s\n", result);
    fprintf(stderr, "    actual: %s\n\n", s);
  }

  cf_free(conv);
  mpz_clear(z);
}

void cf_new_expect_dec(cf_t (*cf_new_fn)(), char *result,
    char *filename, int line) {
  cf_t x = cf_new_fn();
  cf_expect_dec(x, result, filename, line);
  cf_free(x);
}
