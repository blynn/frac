// Solve quadratic equations involving continued fractions.
// A slight misnomer: Gosper describes how to use Newton's method, but
// I use binary search instead, and assume unique roots in the right ranges.

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

// Represents:
//
// a0 xy + a1 x + a2 y + a3
// ------------------------
// a4 xy - a0 x + a5 y - a2

struct newton_data_s {
  cf_t x;
  mpz_t a[6];
  mpz_t lower;
};
typedef struct newton_data_s newton_data_t[1];
typedef struct newton_data_s *newton_data_ptr;

// In the 3D table, we have two sets of equations. Left is a0, b0, c0,
// Right is a1, b1, c1, and these are the terms involving x:
//
//   a0 y + b0     a1 y + b1
//   ---------     ---------
//   c0 y - a0     c1 y - a1
//
// I work left to right, and top to bottom, unlike Gosper.
//
// This notation is different to that in cf_bihom.c, and similar to cf_mobius.

struct abc_s {
  mpz_t a0, a1;
  mpz_t b0, b1;
  mpz_t c0, c1;
};
typedef struct abc_s abc_t[1];

void abc_init(abc_t p) {
  mpz_init(p->a0); mpz_init(p->a1);
  mpz_init(p->b0); mpz_init(p->b1);
  mpz_init(p->c0); mpz_init(p->c1);
}

void abc_clear(abc_t p) {
  mpz_clear(p->a0); mpz_clear(p->a1);
  mpz_clear(p->b0); mpz_clear(p->b1);
  mpz_clear(p->c0); mpz_clear(p->c1);
}

void abc_set_coeff(abc_t p, mpz_t a[6]) {
  mpz_set(p->a0, a[2]); mpz_set(p->b0, a[3]); mpz_set(p->c0, a[5]);
  mpz_set(p->a1, a[0]); mpz_set(p->b1, a[1]); mpz_set(p->c1, a[4]);
}

void abc_print(abc_t p) {
  gmp_printf("%Zd y + %Zd    %Zd y + %Zd\n", p->a0, p->b0, p->a1, p->b1);
  gmp_printf("%Zd y - %Zd    %Zd y - %Zd\n", p->c0, p->a0, p->c1, p->a1);
  gmp_printf("(%Zd+sqrt(%Zd^2+%Zd*%Zd))/%Zd\n", p->a0, p->a0, p->b0, p->c0, p->c0);
  gmp_printf("%Zd*z^2-2*%Zd*z-%Zd\n", p->c0, p->a0, p->b0);
}

// Finds smallest root greater than given lower bound.
// Assumes there exists exactly one root with this property.
// (Impossible to solve equations if roots have same integer part at
// the moment. Don't do that.)
// TODO: Rewrite so you choose if you want smaller or greater root.
// or so it returns both solutions in an array.
static void *newton(cf_t cf) {
  newton_data_ptr nd = cf_data(cf);
  abc_t p;
  abc_init(p);
  abc_set_coeff(p, nd->a);
  cf_t x = nd->x;
  mpz_t z, z0, z1, one, pow2, t0, t1, t2;
  mpz_init(z); mpz_init(z0); mpz_init(z1); mpz_init(pow2);
  mpz_init(t0); mpz_init(t1); mpz_init(t2);
  mpz_init(one);
  mpz_set_ui(one, 1);

  void move_right() {
    cf_get(z, x);
    mpz_mul(t0, z, p->b1);
    mpz_add(t0, t0, p->b0);
    mpz_set(p->b0, p->b1);
    mpz_set(p->b1, t0);

    mpz_mul(t0, z, p->a1);  mpz_mul(t1, z, p->c1);
    mpz_add(t0, t0, p->a0); mpz_add(t1, t1, p->c0);
    mpz_set(p->a0, p->a1);  mpz_set(p->c0, p->c1);
    mpz_set(p->a1, t0);     mpz_set(p->c1, t1);
  }

  void move_down() {
    // Recurrence relation with z, then
    // Subtract and invert z.
    mpz_mul(t0, z, p->a0);
    mpz_add(t0, t0, p->b0);
    mpz_mul(t1, z, p->c0);
    mpz_sub(t1, t1, p->a0);
    mpz_mul(p->a0, z, t1);
    mpz_set(p->b0, p->c0);
    mpz_sub(p->c0, t0, p->a0);
    mpz_set(p->a0, t1);

    mpz_mul(t0, z, p->a1);
    mpz_add(t0, t0, p->b1);
    mpz_mul(t1, z, p->c1);
    mpz_sub(t1, t1, p->a1);
    mpz_mul(p->a1, z, t1);
    mpz_set(p->b1, p->c1);
    mpz_sub(p->c1, t0, p->a1);
    mpz_set(p->a1, t1);
  }

  int sign_quad() {
    // Returns sign of c0 z^2 - 2 a0 z - b0
    mpz_mul_si(t0, p->a0, -2);
    mpz_addmul(t0, p->c0, z);
    mpz_mul(t0, t0, z);
    mpz_sub(t0, t0, p->b0);
    return mpz_sgn(t0);
  }

  int sign_quad1() {
    // Returns sign of c1 z^2 - 2 a1 z - b1
    mpz_mul_si(t0, p->a1, -2);
    mpz_addmul(t0, p->c1, z);
    mpz_mul(t0, t0, z);
    mpz_sub(t0, t0, p->b1);
    return mpz_sgn(t0);
  }

  move_right();  // Get rid of pathological cases.
  move_right();
  move_right();

  // Get integer part, starting search from given lower bound.
  void binary_search(mpz_ptr lower) {
    while (!mpz_sgn(p->c0)) move_right();
    for (;;) {
      mpz_set(z0, lower);
      mpz_set(z, lower);
      int sign = sign_quad();
      mpz_set_ui(pow2, 1);
      for (;;) {
	mpz_add(z, z0, pow2);
	if (sign_quad() != sign) break;
	mpz_mul_2exp(pow2, pow2, 1);
      }
      mpz_set(z1, z);

      for (;;) {
	mpz_add(z, z0, z1);
	mpz_div_2exp(z, z, 1);
	if (!mpz_cmp(z, z0)) break;
	if (sign_quad() == sign) {
	  mpz_set(z0, z);
	} else {
	  mpz_set(z1, z);
	}
      }
      sign = sign_quad1();
      mpz_set(z, z1);
      if (sign_quad1() != sign) break;
      move_right();
    }
  }

  binary_search(nd->lower);
  cf_put(cf, z0);
  mpz_set(z, z0);
  move_down();

  while(cf_wait(cf)) {
    binary_search(one);
    cf_put(cf, z0);
    mpz_set(z, z0);
    move_down();
  }
  abc_clear(p);
  mpz_clear(z); mpz_clear(z0); mpz_clear(z1); mpz_clear(pow2);
  mpz_clear(t0); mpz_clear(t1); mpz_clear(t2);
  for (int i = 0; i < 6; i++) mpz_clear(nd->a[i]);
  free(nd);
  return NULL;
}

cf_t cf_new_newton(cf_t x, mpz_t a[6], mpz_t lower) {
  newton_data_ptr p = malloc(sizeof(*p));
  p->x = x;
  for (int i = 0; i < 6; i++) {
    mpz_init(p->a[i]);
    mpz_set(p->a[i], a[i]);
  }
  mpz_init(p->lower);
  mpz_set(p->lower, lower);
  return cf_new(newton, p);
}

cf_t cf_new_sqrt(cf_t x) {
  newton_data_ptr p = malloc(sizeof(*p));
  p->x = x;
  mpz_init(p->lower);
  for (int i = 0; i < 6; i++) mpz_init(p->a[i]);
  // Solve y = x/y.
  mpz_set_si(p->a[1], 1);
  mpz_set_si(p->a[5], 1);
  return cf_new(newton, p);
}
