/*
 * Bignum.h
 *
 *  Created on: 13 sept. 2022
 *  Abdoulkader Moussa
 *  Maria Ukanwa
 *
 */

#ifndef BIGNUM_HPP_
#define BIGNUM_HPP_
#define MAX_DIV_ERR 100
#define CHAR_BITS 8

#include <iostream>
using namespace std;

class Bignum {
private:
  unsigned *tab; // array for the coefficient of our Bignum
  unsigned top;  // the last element used in the table, idx
  unsigned size; // size of the table (numb elems)
  unsigned sign; // sign pos=0; neg=1, default to pos

  /* convenient way to construct a new Bignum  */
  Bignum(unsigned *tab, unsigned top, unsigned size, unsigned sign);

public:
  /* Constructor  */
  Bignum(int = 0); //"int" equiv "signed int"
  Bignum(unsigned = 0u);
  Bignum(Bignum const &); // copy one Bignum to another
  Bignum(Bignum &&);      // copy and steal ptr to tab

  /* Destructor   */
  ~Bignum();

  /* Set functions to change private members  */
  void setTab(unsigned *);
  void setSignPos();
  void setSignNeg();

  /* Get functions to probe private members   */
  unsigned *getTab();
  unsigned getSign();
  bool isNeg();
  bool isPos();

  /* operator overload functions   */
  Bignum &operator=(Bignum const &); // deep copy
  Bignum &operator=(Bignum &&);      // move
  Bignum &operator+=(Bignum const &);
  Bignum &operator-=(Bignum const &);
  Bignum &operator<<=(uint32_t);
  Bignum &operator>>=(uint32_t);
  friend bool operator==(Bignum const &, Bignum const &);
  friend bool operator!=(Bignum const &, Bignum const &);
  friend bool operator>=(Bignum const &, Bignum const &);
  friend bool operator<=(Bignum const &, Bignum const &);
  friend bool operator>(Bignum const &, Bignum const &);
  friend bool operator<(Bignum const &, Bignum const &);
  friend Bignum operator+(Bignum const &, Bignum const &);
  friend Bignum operator-(Bignum const &, Bignum const &);
  friend Bignum operator*(Bignum const &, Bignum const &);
  friend Bignum operator/(Bignum const &, Bignum const &);
  friend Bignum operator%(Bignum const &, Bignum const &);
  friend Bignum operator&(Bignum const &, Bignum const &);
  friend Bignum operator<<(Bignum const &, uint32_t);
  friend Bignum operator>>(Bignum const &, uint32_t);

  /* Addition of unsigned values   */
  friend Bignum addition(Bignum const &, Bignum const &);

  /* Subtraction of unsigned values*/
  friend Bignum subtraction(Bignum const &, Bignum const &);
  friend Bignum signed_subtr(Bignum const &, Bignum const &);

  /*
   * Supporting functions for += and -=
   * similar to addition/subraction functions, except here
   * there is resource reuse (if possible) of the first Bignum
   */
  friend void add_optimized(Bignum &, Bignum const &);
  friend void subtr_optimized(Bignum &, Bignum const &);

  /*
   * multiplication carried out in tradition manner
   * (like hand & paper). Helper function "mult" used to
   * carry out the intermediate multiplies
   */
  friend Bignum multiply(Bignum const &, Bignum const &);
  friend Bignum mult(Bignum const &, uint32_t, unsigned);

  /*
   * Division carried out in tradional mamer (like hand
   * and paper). Algo is called "muliplty precision division".
   * This algo produces as a result both the quotient and
   * the remainder. Division takes the quotient result, while
   * modulus takes the remainder result.
   * the _10 variant is a fast implementation for divisor=10
   */
  friend pair<Bignum, Bignum> division_modulo(Bignum const &, Bignum const &);
  friend pair<Bignum, Bignum> division_modulo_10(Bignum const &);

  /* Display the Bignum as one single string in hexa  */
  friend ostream &operator<<(ostream &os, const Bignum &b);

  /* Debug function that displays the Bignum coefficients. */
  void showInfo();

  /* Return a string in base10 conversion of a Bignum */
  friend string convertInBase10(Bignum x);
  friend string displayBase10(Bignum x);

  /*
   * unsigned comparison function (used in the overloaded funcs <, > , ==)
   * (op1 == op2) return 0;
   * (op1  > op2) return 1;
   * (op1  < op2) return 2
   */
  friend int op1_is_bigger(Bignum const &b1, Bignum const &b2);

  /* These functions compute x^y(mod n) using the algorithm square and multiply  */
  friend Bignum pow(Bignum x, Bignum y, Bignum n);
  friend Bignum square_and_multiply_algorithm(Bignum x, Bignum y, Bignum n);

  /*
   * Generates a random Bignum with specified number of bits.
   * The number of bits must a multiple of 32.
   */
  friend Bignum generate_random_bignum(unsigned nb_of_bits);

  /*
   * Generates a random prime Bignum with specified number of bits.
   * The number of bits must a multiple of 32.
   */
  friend Bignum generate_prime_bignum(unsigned nb_of_bits);
};

/*
 * Generates a random Bignum with specified number of bits.
 * The number of bits must a multiple of 32.
 */
Bignum generate_random_bignum(unsigned nb_of_bits);

/* Generates a random prime Bignum with specified number of bits. */
Bignum generate_prime_bignum(unsigned nb_of_bits);
Bignum generate_prime_bignum2(unsigned nb_of_bits);

/* Convert a bignum to binary */
string convertToBinary(Bignum x);

/*
 * This function implements the extended Euclidean algorithm
 * given two inputs A, M, it calculates: A*x + M*y = 1
 * Returns the tuple = {gcd(A,M), x, y}, where x = modular inverse
 */
tuple<Bignum, Bignum, Bignum> gcdExtended(Bignum const &, Bignum const &);

/*
 * This function checks if the modular inverse exists based on the gcd of
 * the two inputs. If gcd is 1, returns the modular inverse, else returns 0
 */
Bignum modInverse(Bignum const &, Bignum const &);

/*
 * Implementation of the RSA algorithm
 * Encrypts and decrypts the input message
 */
void RSA(string msg, unsigned nb_of_bits_prime_numbers);

/* given a Bignum (representing a plaintext), it is converted to ascii */
string convertPTtoAscii(Bignum const &);

// temp funcs - int versions- use as reference for numbers lt 64bits
int gcdExtended_int(int, int, int *, int *);
void modInverse_int(int, int);

#endif /* BIGNUM_HPP_ */
