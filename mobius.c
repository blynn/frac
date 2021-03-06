#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

// Minds our p's and q's. The two last computed convergents.
struct pqset_s {
  mpz_t pold, p;
  mpz_t qold, q;
};
typedef struct pqset_s *pqset_ptr;
typedef struct pqset_s pqset_t[1];

void pqset_init(pqset_t pq) {
  mpz_init(pq->pold);
  mpz_init(pq->qold);
  mpz_init(pq->p);
  mpz_init(pq->q);
}

void pqset_clear(pqset_t pq) {
  mpz_clear(pq->pold);
  mpz_clear(pq->qold);
  mpz_clear(pq->p);
  mpz_clear(pq->q);
}

void pqset_neg(pqset_t pq) {
  mpz_neg(pq->p, pq->p);
  mpz_neg(pq->q, pq->q);
  mpz_neg(pq->pold, pq->pold);
  mpz_neg(pq->qold, pq->qold);
}

void pqset_print(pqset_t pq) {
  gmp_printf("p's: %Zd %Zd\n", pq->pold, pq->p);
  gmp_printf("q's: %Zd %Zd\n", pq->qold, pq->q);
}

// Compute the next convergent for regular continued fractions.
void pqset_regular_recur(pqset_t pq, mpz_t denom) {
  mpz_addmul(pq->pold, denom, pq->p);
  mpz_swap(pq->pold, pq->p);
  mpz_addmul(pq->qold, denom, pq->q);
  mpz_swap(pq->qold, pq->q);
}

// Compute the next convergent for nonregular continued fractions.
void pqset_nonregular_recur(pqset_t pq, mpz_t num, mpz_t denom) {
  mpz_mul(pq->pold, pq->pold, num);
  mpz_addmul(pq->pold, pq->p, denom);
  mpz_swap(pq->pold, pq->p);
  mpz_mul(pq->qold, pq->qold, num);
  mpz_addmul(pq->qold, pq->q, denom);
  mpz_swap(pq->qold, pq->q);
}

// Get rid of nontrivial GCD for {p, q, pold, qold}.
// t0 and t1 are temporary variables.
void pqset_remove_gcd(pqset_ptr pq, mpz_t t0, mpz_t t1) {
  mpz_gcd(t0, pq->p, pq->q);
  mpz_gcd(t1, pq->pold, pq->qold);
  mpz_gcd(t0, t0, t1);
  if (mpz_cmp_ui(t0, 1)) {
    mpz_divexact(pq->pold, pq->pold, t0);
    mpz_divexact(pq->qold, pq->qold, t0);
    mpz_divexact(pq->p, pq->p, t0);
    mpz_divexact(pq->q, pq->q, t0);
  }
}

// A Mobius transformation: four coefficients and the input.
// TODO: Use an array of size 4.
struct mobius_data_s {
  cf_t input;
  mpz_t a, b, c, d;
};
typedef struct mobius_data_s *mobius_data_ptr;

void pqset_set_mobius(pqset_t pq, mobius_data_ptr md) {
  mpz_set(pq->pold, md->b); mpz_set(pq->p, md->a);
  mpz_set(pq->qold, md->d); mpz_set(pq->q, md->c);
}

// Compute convergents of Mobius function applied to a regular
// continued fraction.
static void *mobius_convergent(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;
  pqset_t pq;
  pqset_init(pq);
  pqset_set_mobius(pq, md);

  mpz_t denom;
  mpz_init(denom);
  while(cf_wait(cf)) {
    cf_get(denom, input);
    pqset_regular_recur(pq, denom);

    cf_put(cf, pq->p);
    cf_put(cf, pq->q);
  }
  mpz_clear(denom);
  pqset_clear(pq);

  mpz_clear(md->a);
  mpz_clear(md->b);
  mpz_clear(md->c);
  mpz_clear(md->d);
  free(md);
  return NULL;
}
// Start a thread that, when signalled, computes the convergents of a Mobius
// transformation of a continued fraction.
cf_t cf_new_mobius_convergent(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a); mpz_init(md->b); mpz_init(md->c); mpz_init(md->d);
  mpz_set(md->a, a); mpz_set(md->b, b); mpz_set(md->c, c); mpz_set(md->d, d);
  md->input = x;
  return cf_new(mobius_convergent, md);
}

// Start a thread that, when signalled, computes the convergents of a continued
// fraction.
cf_t cf_new_cf_convergent(cf_t x) {
  mpz_t one, zero;
  mpz_init(one); mpz_init(zero);
  mpz_set_ui(one, 1); mpz_set_ui(zero, 0);
  cf_t res = cf_new_mobius_convergent(x, one, zero, zero, one);
  mpz_clear(one); mpz_clear(zero);
  return res;
}

// Compute nonregular convergents of a Mobius function applied
// to a nonregular continued fraction.
static void *nonregular_mobius_convergent(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;
  pqset_t pq; pqset_init(pq); pqset_set_mobius(pq, md);
  mpz_t num; mpz_init(num);
  mpz_t denom; mpz_init(denom);
  mpz_t t0, t1; mpz_init(t0); mpz_init(t1);
  void recur() {
    pqset_nonregular_recur(pq, num, denom);
    pqset_remove_gcd(pq, t0, t1);

    cf_put(cf, pq->p);
    cf_put(cf, pq->q);
  }
  mpz_set_ui(num, 1);
  cf_get(denom, input);
  recur();
  while(cf_wait(cf)) {
    cf_get(num, input);
    cf_get(denom, input);
    recur();
  }
  mpz_clear(num);
  mpz_clear(denom);
  pqset_clear(pq);

  mpz_clear(md->a); mpz_clear(md->b); mpz_clear(md->c); mpz_clear(md->d);
  mpz_clear(t0); mpz_clear(t1);
  free(md);
  return NULL;
}

cf_t cf_new_nonregular_mobius_convergent(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a); mpz_init(md->b); mpz_init(md->c); mpz_init(md->d);
  mpz_set(md->a, a); mpz_set(md->b, b); mpz_set(md->c, c); mpz_set(md->d, d);
  md->input = x;
  return cf_new(nonregular_mobius_convergent, md);
}

// Input: Mobius transformation and nonregular continued fraction.
// Output: Regular continued fraction. Assumes input fraction is well-behaved.
static void *mobius_nonregular_throughput(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;
  pqset_t pq; pqset_init(pq); pqset_set_mobius(pq, md);
  mpz_t num; mpz_init(num);
  mpz_t denom; mpz_init(denom);
  mpz_t t0, t1, t2; mpz_init(t2); mpz_init(t1); mpz_init(t0);
  int recur() {
    pqset_nonregular_recur(pq, num, denom);
    pqset_remove_gcd(pq, t0, t1);

    if (mpz_sgn(pq->qold)) {
      mpz_fdiv_qr(t1, t0, pq->pold, pq->qold);
      mpz_mul(t2, t1, pq->q);

      if (mpz_cmp(t2, pq->p) <= 0) {
	mpz_add(t2, t2, pq->q);
	if (mpz_cmp(t2, pq->p) > 0) {
	  // Output continued fraction term.
	  cf_put(cf, t1);
	  // Subtract: remainder of p/q.
	  mpz_sub(t2, t2, pq->p);
	  mpz_sub(t2, pq->q, t2);
	  // Invert
	  mpz_set(pq->pold, pq->qold);
	  mpz_set(pq->qold, t0);
	  mpz_set(pq->p, pq->q);
	  mpz_set(pq->q, t2);
	  return 1;
	}
      }
    }
    return 0;
  }
  mpz_set_ui(num, 1);
  cf_get(denom, input);
  recur();
  while(cf_wait(cf)) {
    do {
      cf_get(num, input);
      cf_get(denom, input);
    } while(!recur());
  }
  mpz_clear(num);
  mpz_clear(denom);
  pqset_clear(pq);

  mpz_clear(md->a); mpz_clear(md->b); mpz_clear(md->c); mpz_clear(md->d);
  free(md);
  mpz_clear(t2); mpz_clear(t1); mpz_clear(t0);
  return NULL;
}

cf_t cf_new_nonregular_to_cf(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a); mpz_init(md->b); mpz_init(md->c); mpz_init(md->d);
  mpz_set(md->a, a); mpz_set(md->b, b); mpz_set(md->c, c); mpz_set(md->d, d);
  md->input = x;
  return cf_new(mobius_nonregular_throughput, md);
}

// This seems to be slower than regularizing the continued fraction
// and then converting to decimal.
static void *nonregular_mobius_decimal(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;
  pqset_t pq; pqset_init(pq); pqset_set_mobius(pq, md);
  mpz_t num; mpz_init(num);
  mpz_t denom; mpz_init(denom);
  mpz_t t0, t1, t2; mpz_init(t2); mpz_init(t1); mpz_init(t0);
  int recur() {
    pqset_nonregular_recur(pq, num, denom);
    pqset_remove_gcd(pq, t0, t1);

    // If the denominator is zero, we can't do anything yet.
    if (mpz_sgn(pq->qold)) {
      mpz_fdiv_qr(t1, t0, pq->pold, pq->qold);
      mpz_mul(t2, t1, pq->q);
      if (mpz_cmp(t2, pq->p) <= 0) {
	mpz_add(t2, t2, pq->q);
	if (mpz_cmp(t2, pq->p) > 0) {
	  // Output a decimal digit.
	  cf_put(cf, t1);
	  // Subtract: remainder of p/q.
	  mpz_sub(t2, t2, pq->p);
	  mpz_sub(t2, pq->q, t2);
	  // Multiply numerator by 10.
	  mpz_mul_ui(pq->pold, t0, 10);
	  mpz_mul_ui(pq->p, t2, 10);
	  return 1;
	}
      }
    }
    return 0;
  }
  mpz_set_ui(num, 1);
  cf_get(denom, input);
  recur();
  while(cf_wait(cf)) {
    do {
      cf_get(num, input);
      cf_get(denom, input);
    } while(!recur());
  }
  mpz_clear(num);
  mpz_clear(denom);
  pqset_clear(pq);
  mpz_clear(t0); mpz_clear(t1); mpz_clear(t2);
  mpz_clear(md->a); mpz_clear(md->b); mpz_clear(md->c); mpz_clear(md->d);
  free(md);
  return NULL;
}

cf_t cf_new_nonregular_mobius_to_decimal(cf_t x, mpz_t a[4]) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a); mpz_init(md->b); mpz_init(md->c); mpz_init(md->d);
  mpz_set(md->a, a[0]); mpz_set(md->b, a[1]);
  mpz_set(md->c, a[2]); mpz_set(md->d, a[3]);
  md->input = x;
  return cf_new(nonregular_mobius_decimal, md);
}

static void determine_sign(cf_t cf, pqset_t pq, mpz_t denom, cf_t input) {
  cf_set_sign(cf, cf_sign(input));
  while (mpz_sgn(pq->pold) != mpz_sgn(pq->p)
      || mpz_sgn(pq->qold) != mpz_sgn(pq->q)) {
    cf_get(denom, input);
    pqset_regular_recur(pq, denom);
  }
  if (mpz_sgn(pq->qold) < 0) {
    mpz_neg(pq->qold, pq->qold);
    mpz_neg(pq->q, pq->q);
    cf_flip_sign(cf);
  }
  if (mpz_sgn(pq->pold) < 0) {
    mpz_neg(pq->pold, pq->pold);
    mpz_neg(pq->p, pq->p);
    cf_flip_sign(cf);
  }
}

// Input: Mobius transformation and regular continued fraction.
// Output: Regular continued fraction.
static void *mobius_throughput(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;
  pqset_t pq; pqset_init(pq); pqset_set_mobius(pq, md);
  mpz_t denom; mpz_init(denom);
  mpz_t t0, t1, t2; mpz_init(t2); mpz_init(t1); mpz_init(t0);

  determine_sign(cf, pq, denom, input);

  int recur() {
    pqset_regular_recur(pq, denom);

    if (mpz_sgn(pq->qold)) {
      mpz_fdiv_qr(t1, t0, pq->pold, pq->qold);
      mpz_mul(t2, t1, pq->q);

      if (mpz_cmp(t2, pq->p) <= 0) {
	mpz_add(t2, t2, pq->q);
	if (mpz_cmp(t2, pq->p) > 0) {
	  // Output continued fraction term.
	  cf_put(cf, t1);
	  // Subtract: remainder of p/q.
	  mpz_sub(t2, t2, pq->p);
	  mpz_sub(t2, pq->q, t2);
	  // Invert
	  mpz_set(pq->pold, pq->qold);
	  mpz_set(pq->qold, t0);
	  mpz_set(pq->p, pq->q);
	  mpz_set(pq->q, t2);
	  return 1;
	}
      }
    }
    return 0;
  }
  while(cf_wait(cf)) {
    do {
      cf_get(denom, input);
    } while(!recur());
  }
  mpz_clear(denom);
  pqset_clear(pq);
  mpz_clear(t0); mpz_clear(t1); mpz_clear(t2);
  mpz_clear(md->a); mpz_clear(md->b); mpz_clear(md->c); mpz_clear(md->d);
  free(md);
  return NULL;
}

cf_t cf_new_mobius_to_cf(cf_t x, mpz_t z[4]) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a); mpz_init(md->b); mpz_init(md->c); mpz_init(md->d);
  mpz_set(md->a, z[0]); mpz_set(md->b, z[1]);
  mpz_set(md->c, z[2]); mpz_set(md->d, z[3]);
  md->input = x;
  return cf_new(mobius_throughput, md);
}

// Input: Mobius transformation and regular continued fraction.
// Output: Decimal representation. The integer part is given first,
// followed by one digit at a time.
static void *mobius_decimal(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;
  pqset_t pq; pqset_init(pq); pqset_set_mobius(pq, md);
  mpz_t denom; mpz_init(denom);
  mpz_t t0, t1, t2; mpz_init(t2); mpz_init(t1); mpz_init(t0);

  determine_sign(cf, pq, denom, input);
  int recur() {
    pqset_regular_recur(pq, denom);

    // If the denominator is zero, we can't do anything yet.
    if (mpz_sgn(pq->qold)) {
      // Each term except possibly the first is one of {0, ..., 9}.
      /* Naive attempt to expoit this didn't seem faster:
       * (and I'd have to handle the first term properly)
      int i;
      mpz_set(t0, pq->q);
      for (i = 0; i <= 9; i++) {
	if (mpz_cmp(t0, pq->p) > 0) break;
	mpz_add(t0, t0, pq->q);
      }
      mpz_set_ui(pq->pold, i);
      mpz_sub(t0, t0, pq->p);
      mpz_sub(t0, pq->q, t0);
      mpz_mul(pq->qold, pq->pold, pq->qnew);
      */

      mpz_fdiv_qr(t1, t0, pq->pold, pq->qold);
      mpz_mul(t2, t1, pq->q);
      if (mpz_cmp(t2, pq->p) <= 0) {
	mpz_add(t2, t2, pq->q);
	if (mpz_cmp(t2, pq->p) > 0) {
	  // Output a decimal digit.
	  cf_put(cf, t1);
	  // Compute t2 = remainder of p/q.
	  mpz_sub(t2, t2, pq->p);
	  mpz_sub(t2, pq->q, t2);
	  // Multiply numerator by 10.
	  mpz_mul_ui(pq->pold, t0, 10);
	  mpz_mul_ui(pq->p, t2, 10);
	  return 1;
	}
      }
    }
    return 0;
  }
  while(cf_wait(cf)) {
    do {
      cf_get(denom, input);
    } while(!recur());
  }
  mpz_clear(denom);
  pqset_clear(pq);
  mpz_clear(t0); mpz_clear(t1); mpz_clear(t2);
  mpz_clear(md->a); mpz_clear(md->b); mpz_clear(md->c); mpz_clear(md->d);
  free(md);
  return NULL;
}

cf_t cf_new_mobius_to_decimal(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a); mpz_init(md->b); mpz_init(md->c); mpz_init(md->d);
  mpz_set(md->a, a); mpz_set(md->b, b); mpz_set(md->c, c); mpz_set(md->d, d);
  md->input = x;
  return cf_new(mobius_decimal, md);
}

cf_t cf_new_cf_to_decimal(cf_t x) {
  mpz_t one, zero;
  mpz_init(one); mpz_init(zero);
  mpz_set_ui(one, 1); mpz_set_ui(zero, 0);
  cf_t res = cf_new_mobius_to_decimal(x, one, zero, zero, one);
  mpz_clear(one); mpz_clear(zero);
  return res;
}

struct funarg_s {
  void *(*fun)(cf_t);
  mpz_t arg;
};
typedef struct funarg_s *funarg_ptr;

static void *one_arg_nonregular(cf_t cf) {
  funarg_ptr p = cf_data(cf);
  mpz_t a, b, c, d;
  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);
  mpz_set_ui(a, 1); mpz_set_ui(b, 0);
  mpz_set_ui(c, 0); mpz_set_ui(d, 1);
  mpz_ptr copy = malloc(sizeof(*copy));
  mpz_init(copy);
  mpz_set(copy, p->arg);
  cf_t nonregular = cf_new(p->fun, copy);
  cf_t conv = cf_new_nonregular_to_cf(nonregular, a, b, c, d);
  mpz_t z;
  mpz_init(z);
  while(cf_wait(cf)) {
    cf_get(z, conv);
    cf_put(cf, z);
  }
  mpz_clear(z);
  cf_free(conv);
  cf_free(nonregular);
  mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);

  mpz_clear(p->arg);
  free(p);
  return NULL;
}

cf_t cf_new_one_arg(void *(*fun)(cf_t), mpz_t z) {
  mpz_ptr p = malloc(sizeof(*p));
  mpz_init(p);
  mpz_set(p, z);
  return cf_new(fun, p);
}

cf_t cf_new_one_arg_nonregular(void *(*fun)(cf_t), mpz_t z) {
  funarg_ptr p = malloc(sizeof(*p));
  p->fun = fun;
  mpz_init(p->arg);
  mpz_set(p->arg, z);
  return cf_new(one_arg_nonregular, p);
}

static void *nonreg_const(cf_t cf) {
  funarg_ptr p = cf_data(cf);
  mpz_t a, b, c, d;
  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);
  mpz_set_ui(a, 1); mpz_set_ui(b, 0);
  mpz_set_ui(c, 0); mpz_set_ui(d, 1);
  cf_t nonregular = cf_new(p->fun, NULL);
  cf_t conv = cf_new_nonregular_to_cf(nonregular, a, b, c, d);
  mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);
  mpz_t z;
  mpz_init(z);
  while(cf_wait(cf)) {
    cf_get(z, conv);
    cf_put(cf, z);
  }
  mpz_clear(z);
  cf_free(conv);
  cf_free(nonregular);
  free(p);
  return NULL;
}

cf_t cf_new_const_nonregular(void *(*fun)(cf_t)) {
  funarg_ptr p = malloc(sizeof(*p));
  p->fun = fun;
  return cf_new(nonreg_const, p);
}
