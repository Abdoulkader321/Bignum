/*
 * Bignum.cpp
 *
 *  Created on: 13 sept. 2022
 *  Abdoulkader Moussa
 *  Maria Ukanwa
 *
 * Compile:
 *  + make
 *  + ./Main
 *
 */

#include "Bignum.hpp"
#include <algorithm>
#include <climits>
#include <iomanip>
#include <sstream>
#include <tuple>
#include <vector>

using namespace std;
static bool seed_intialized = false;

Bignum::Bignum(int n) : top(1), size(1) {
  tab = new unsigned[1];
  if (n >= 0) {
    sign = 0;
    tab[0] = n;
  } else { // negative, n * -1;
    sign = 1;
    tab[0] = -n;
  }
}

Bignum::Bignum(unsigned n) : top(1), size(1), sign(0) {
  tab = new unsigned[1];
  tab[0] = n;
}

Bignum::Bignum(Bignum const &x) : top(x.top), size(x.size), sign(x.sign) {
  tab = new unsigned[x.size];
  for (unsigned i = 0; i < top; i++)
    tab[i] = x.tab[i];
}

Bignum::Bignum(Bignum &&x)
    : tab(x.tab), top(x.top), size(x.size), sign(x.sign) {
  x.tab = nullptr; // 'this' is new owner of the tab table
}

Bignum::Bignum(unsigned *xtab, unsigned xtop, unsigned xsize, unsigned xsign)
    : tab(xtab), top(xtop), size(xsize), sign(xsign) {}

Bignum::~Bignum() { delete[] tab; }

/*
 * copy overload is a deep copy, the size of the current Bignum
 * needs to be big enough or increased if necessary
 */
Bignum &Bignum::operator=(Bignum const &x) {
  if (x.size > size) {
    delete[] tab;
    tab = new unsigned[x.size];
    size = x.size;
  }

  top = x.top;
  sign = x.sign;
  for (unsigned i = 0; i < x.top; i++) {
    tab[i] = x.tab[i];
  }
  return *this;
}

/*
 * move overload is a shallow copy, ownership of the
 * the tab is transfered from the old Bignum to the new Bignum
 */
Bignum &Bignum::operator=(Bignum &&x) {
  delete[] tab;
  tab = x.tab;
  top = x.top;
  size = x.size;
  sign = x.sign;
  x.tab = nullptr; // insure tab not leaked when x deleted
  return *this;
}

unsigned *Bignum::getTab() { return this->tab; }
unsigned Bignum::getSign() { return this->sign; }
bool Bignum::isNeg() { return (this->sign == 1); }
bool Bignum::isPos() { return (this->sign == 0); }

void Bignum::setTab(unsigned *new_tab) { this->tab = new_tab; }
void Bignum::setSignPos() { this->sign = 0; }
void Bignum::setSignNeg() { this->sign = 1; }

/*
 * Left Shift operator
 * (1) Based on the shift amount, calculate the size of the new tab
 * (2) Data is shifted (moved) from the msb blk of the Bignum input into the
 *     msb blk of the result. This shift is done by simulating the shift
 *     of each blk(32b) into a 64b block and then moving the upper and
 *     lower parts accordingly into the result
 */
Bignum operator<<(Bignum const &x, uint32_t shft) {
  unsigned top;
  unsigned size;
  unsigned sign = x.sign;
  unsigned *tab;

  uint64_t shft_result;
  uint64_t extended_x;
  uint32_t lower_32;
  uint32_t upper_32;
  unsigned last_blk;
  unsigned old_top = x.top;

  /* calculate tab needed to accomodate the shift */
  uint32_t max_shft, shft_blks, shft_bits, shft_bits_fixup, tot_blks,
      shft_pre_fixup;
  max_shft = sizeof(uint32_t) * 8;
  shft_blks = shft / max_shft;  // assume floor fnc
  shft_bits = shft % max_shft;  // shift amt unless shift is multiple max_shft
  tot_blks = shft_blks + x.top; // original tab size + shift blks
  shft_bits = shft % max_shft;

  /* fixup if shift is multiple of max_shift, if shft==0, no fixup */
  shft_pre_fixup = shft_blks ? max_shft : 0;
  shft_bits_fixup = shft_bits ? shft_bits : shft_pre_fixup; // actual shift

  size = shft % max_shft ? tot_blks + 1 : tot_blks;

  top = size;
  tab = new unsigned[size];
  last_blk = shft_bits ? top - (old_top + 1) : top - old_top;

  /* start shifting into the msb blk of the result  */
  for (unsigned s = top; s >= 1; s--) {

    if (old_top > 0) {
      extended_x = 0; // upper 32b must be cleared
      extended_x = x.tab[old_top - 1];
      shft_result = extended_x << shft_bits_fixup;
      lower_32 = (uint32_t)shft_result;
      upper_32 = shft_result >> 32;

      if (s == top) {
        tab[s - 1] = upper_32;
      } else
        tab[s - 1] = tab[s - 1] | upper_32;

      if (s >= 2)
        tab[s - 2] = lower_32;

      old_top--;

    } else {
      if (s <= last_blk) // last blk with valid data
        tab[s - 1] = 0;
    }
  } // end for-loop

  while ((top > 1) && (tab[top - 1] == 0)) // adjust, find the real top
    top--;

  return Bignum(tab, top, size, sign);
}

/*
 * Right Shift operator
 * (1) Based on the shift amount, calculate the size of the new tab
 * (2) Data is shifted (moved) from the msb blk of the Bignum input into what
 *     would be the msb blk of the result. This shift is done
 *     by simulating the shift of each blk(32b) into a 64b block and
 *     then moving the upper and lower parts accordingly into the result.
 *     Bits shifted to the right of the one's place value are dropped
 * (3) Future optimiztion, reuse the input BigNum resource, do not overwrite
 *     the original tab before that data is shifted
 */
Bignum operator>>(Bignum const &x, uint32_t shft) {
  unsigned top;
  unsigned old_top = x.top;
  unsigned size;
  unsigned sign = x.sign;
  unsigned *tab;
  uint64_t shft_result;
  uint64_t extended_x;
  uint32_t lower_32;
  uint32_t upper_32;

  /* calculate reduced new tab   */
  uint32_t blks_lost, s_idx;
  uint32_t max_shft, shft_blks, shft_bits, shft_pre_fixup, shft_bits_fixup,
      shft_blks_fixup;
  max_shft = sizeof(uint32_t) * 8;
  shft_blks = shft / max_shft; // assume floor function
  shft_bits = shft % max_shft; // shift amt unless shift is multiple max_shft
  blks_lost = shft_blks;

  shft_pre_fixup = shft_blks ? max_shft : 0;
  shft_bits_fixup = shft_bits ? shft_bits : shft_pre_fixup; // actual shift
  shft_blks_fixup = (shft_bits == 0) ? shft_blks - 1 : shft_blks;

  size = (old_top <= blks_lost) ? 1 : old_top - blks_lost;
  top = size;
  tab = new unsigned[size];

  /* start shifting at the the lsb blk of the result */
  for (unsigned s = 0; s < old_top; s++) {
    if (((s + shft_blks_fixup) > old_top) || (s > top)) {
      break;
    } else {
      s_idx = s + shft_blks_fixup;
      extended_x = 0; // upper 32b must be cleared
      extended_x = x.tab[s_idx];
      extended_x <<= 32;
      shft_result = extended_x >> shft_bits_fixup;
      lower_32 = (uint32_t)shft_result;
      upper_32 = shft_result >> 32;

      if (s != 0) {
        tab[s - 1] = tab[s - 1] | lower_32;
      }

      tab[s] = upper_32;
    }
  }

  while ((top > 1) && (tab[top - 1] == 0)) // adjust, find the real top
    top--;

  return Bignum(tab, top, size, sign);
}

Bignum &Bignum::operator<<=(uint32_t shft) {
  *this = *this << shft;
  return *this;
}

Bignum &Bignum::operator>>=(uint32_t shft) {
  *this = *this >> shft;
  return *this;
}

bool operator>(Bignum const &x, Bignum const &y) {
  return (op1_is_bigger(x, y) == 1);
}

bool operator<(Bignum const &x, Bignum const &y) {
  return (op1_is_bigger(x, y) == 2);
}

bool operator==(Bignum const &x, Bignum const &y) {
  return (op1_is_bigger(x, y) == 0);
}

bool operator!=(Bignum const &x, Bignum const &y) {
  return (!(op1_is_bigger(x, y) == 0));
}

bool operator>=(Bignum const &x, Bignum const &y) {
  return (op1_is_bigger(x, y) != 2);
}

bool operator<=(Bignum const &x, Bignum const &y) {
  return (op1_is_bigger(x, y) != 1);
}

Bignum operator&(Bignum const &x, Bignum const &y) {

  /* recursive call to insure x > y */
  if (x < y)
    return (y & x);

  unsigned top = x.top;
  unsigned size = x.size;
  unsigned sign = x.sign;
  unsigned *tab;

  tab = new unsigned[size]; // use size of bigger operand

  for (unsigned idx = 0; idx < x.top; idx++)
    tab[idx] = x.tab[idx] & y.tab[idx];

  while ((top > 1) && (tab[top - 1] == 0)) // adjust, find the real top
    top--;

  return Bignum(tab, top, size, sign);
}

Bignum operator+(Bignum const &x, Bignum const &y) {
  Bignum bn(1);

  if (x.sign == y.sign) {
    bn = addition(x, y); // Same sign : addition
  } else {
    bn = subtraction(x, y); // Diff sign : subtraction (bigger-smaller)
  }
  return bn;
}

Bignum operator-(Bignum const &x, Bignum const &y) {
  return (signed_subtr(x, y));
}

/*
 * Addition on unsigned values (coefficients) in the tab.
 * Sign is handled seperately. Addition of coefficients at each
 * place value is done in isolation taking into account the carry
 * from the previous blk, and in turn propagating the carry generated
 * at this current blk to the next place value.
 */
Bignum addition(Bignum const &x, Bignum const &y) {

  /* recursive call to insure x > y */
  if (x.top < y.top)
    return addition(y, x);

  /* mem allocation expensive, insure allocated tab is gt sufficient */
  unsigned size = (x.top == x.size) ? 2 * x.size : x.size;
  unsigned *tab = new unsigned[size];

  uint64_t carry = 0;
  uint64_t sum;

  /* addition of cofficients in the tab for one place value at a time */
  for (unsigned idx = 0; idx < y.top; idx++) {
    // use casting, so that precision is not lost, want uint64_t precision
    sum = (uint64_t)x.tab[idx] + (uint64_t)y.tab[idx] + carry;
    carry = sum >> 32; // check for overflow
    tab[idx] = sum;
  }

  /* finish up the addition of x and carry */
  for (unsigned idx = y.top; idx < x.top; idx++) {
    sum = carry + (uint64_t)x.tab[idx];
    carry = sum >> 32;
    tab[idx] = sum;
  }

  tab[x.top] = carry; // if carry=1, top will advance (below)
  return Bignum(tab, x.top + carry, size, x.sign);
}

Bignum signed_subtr(Bignum const &x, Bignum const &y) {
  Bignum bn(1);

  /* take into account 4 different sign combinations */
  if ((x.sign == 0) && (y.sign == 0)) {
    bn = subtraction(x, y); // x bigger -> pos, b2 bigger -> neg
  }
  if ((x.sign == 1) && (y.sign == 1)) {
    bn = subtraction(y, x); // x bigger -> neg, b2 bigger -> pos
  }

  /* sign adjustment for 2 cases above  */
  bn.sign = (x == y) ? (unsigned)0 : ((x > y) ? x.sign : !y.sign);

  if ((x.sign == 0) && (y.sign == 1)) {
    bn = addition(x, y);
    bn.sign = 0; // fixup for pos result
  }
  if ((x.sign == 1) && (y.sign == 0)) {
    bn = addition(x, y);
    bn.sign = 1; // fixup for neg result
  }

  return bn;
}

/*
 * Signs of input operand determine if the operation is
 * subraction or addition. Seperating the sign allows subtraction on
 * unsigned values in the tab, always insure bigger BN - smaller BN
 * however for some coefficients will require regrouping (small - big):
 * This is handled via the carry propagated one place value to the next
 */
Bignum subtraction(Bignum const &x, Bignum const &y) {

  /* recursive call to insure x > y */
  if (x < y)
    return subtraction(y, x);

  unsigned top = x.top;
  unsigned size = x.size;
  unsigned sign = x.sign;
  unsigned *tab;
  unsigned carry = 0;

  tab = new unsigned[size]; // use size of bigger operand

  for (unsigned idx = 0; idx < y.top; idx++) {    // x > y
    tab[idx] = x.tab[idx] - (y.tab[idx] + carry); // carry=regrouping term
    carry = (x.tab[idx] >= (y.tab[idx] + carry)) &&
                    ((y.tab[idx] != UINT_MAX) || (carry != 1))
                ? (unsigned)0
                : (unsigned)1;
  }

  /*Finish up the subtraction of the bigger number BN and the final carry*/
  for (unsigned i = y.top; i < x.top; i++) {
    tab[i] = x.tab[i] - carry;
    carry = ((x.tab[i] >= carry) && (carry != 1)) ? (unsigned)0 : (unsigned)1;
  }

  while ((top > 1) && (tab[top - 1] == 0)) // adjust, find the real top
    top--;

  return Bignum(tab, top, size, sign);
}

/*
 * add_optimized and sub_optimized are supporting funcs for += and -=
 * similar to addition/subraction functions, except here there is a
 * possible resource optimization: if Bignum x has sufficient space,
 * its tab resource is reused.
 * assumption: |op1| > |op2|, must verify before this funct is called
 */
void add_optimized(Bignum &x, Bignum const &y) {
  uint64_t carry = 0;
  uint64_t sum;

  /* addition of cofficients in the tab for one place value at a time */
  for (unsigned idx = 0; idx < y.top; idx++) {
    sum = x.tab[idx] + y.tab[idx] + carry;
    carry = sum >> 32; // check for overflow
    x.tab[idx] = sum;
  }

  /* finish up the addtion of x and carry */
  for (unsigned idx = y.top; idx < x.top; idx++) {
    sum = carry + x.tab[idx];
    carry = sum >> 32;
    x.tab[idx] = sum;
  }

  x.tab[x.top] = carry; // if carry=1, top needs to advance;
  x.top = x.top + carry;
}

void subtr_optimized(Bignum &x, Bignum const &y) {
  unsigned carry = 0;
  unsigned tmp_sub;

  for (unsigned idx = 0; idx < y.top; idx++) {   // x > y
    tmp_sub = x.tab[idx] - (y.tab[idx] + carry); // carry=regrouping term
    carry = (x.tab[idx] >= (y.tab[idx] + carry)) &&
                    ((y.tab[idx] != UINT_MAX) || (carry != 1))
                ? (unsigned)0
                : (unsigned)1;
    x.tab[idx] = tmp_sub;
  }

  for (unsigned i = y.top; i < x.top; i++) {
    tmp_sub = x.tab[i] - carry;
    carry = ((x.tab[i] >= carry) && (carry != 1)) ? (unsigned)0 : (unsigned)1;
    x.tab[i] = tmp_sub;
  }

  while ((x.top > 1) && (x.tab[x.top - 1] == 0)) // find real top
    x.top--;
}

/*
  unsigned comparator
  if op1 == op2 return 0
  if op1 >  op2 return 1
  if op1 <  op2 return 2
*/
int op1_is_bigger(Bignum const &x, Bignum const &y) {
  unsigned xtop_adj = x.top;
  unsigned ytop_adj = y.top;

  /* find position of msb non-zero coeff*/
  while (xtop_adj > 1 && (x.tab[xtop_adj - 1] == 0))
    xtop_adj--;

  while (ytop_adj > 1 && (y.tab[ytop_adj - 1] == 0))
    ytop_adj--;

  if (xtop_adj == ytop_adj) {
    while (xtop_adj > 0 && (x.tab[xtop_adj - 1] == y.tab[xtop_adj - 1])) {
      xtop_adj--;
    }
    if (xtop_adj == (uint32_t)0)
      return 0;
    else if (x.tab[xtop_adj - 1] > y.tab[xtop_adj - 1])
      return 1;
    else
      return 2;

  } else if (xtop_adj > ytop_adj) {
    return 1;
  } else
    return 2;
}

/* Same as +,- except add_optimization and sub_optimization are called */
Bignum &Bignum::operator+=(Bignum const &y) {

  unsigned min_size = (this->top < y.top) ? y.top + 1 : this->top + 1;

  if (this->size >= min_size) { // reuse 'this' mem resource
    if (this->sign == y.sign) {
      if (*this > y) {
        add_optimized(*this, y); // Same sign do addition, sign comes from *this
      } else {
        *this = *this + y; // normal add
      }
    } else {
      if (*this > y) {             //*this bigger, sign already correct
        subtr_optimized(*this, y); // Diff sign do subtraction (bigger-smaller)
      } else {
        *this = *this - y; // normal subtraction
      }
    }

  } else {
    *this = *this + y; // 'this' not have sufficient resource, normal add
  }

  return *this;
}

Bignum &Bignum::operator-=(Bignum const &y) {

  unsigned min_size = (this->top < y.top) ? y.top : this->top; // is +1 needed?

  if (min_size <= this->size) { // reuse 'this' mem resource
    if (this->sign != y.sign) {
      if (*this > y) {
        add_optimized(*this, y); // Same sign do addition, sign comes from *this
      } else {
        *this = *this + y; // normal add
      }
    } else {
      if (*this > y) {             //*this bigger, sign already correct
        subtr_optimized(*this, y); // Diff sign do subtraction (bigger-smaller)
      } else {
        *this = *this - y; // normal subtraction
      }
    }

  } else {
    *this = *this - y; // 'this' not have sufficient resource, normal add
  }

  return *this;
}

Bignum operator*(Bignum const &x, Bignum const &y) {
  Bignum prod = multiply(x, y);
  return prod;
}

/*
 * multiplication of the multiplicand bignum with one coefficient
 * of the multiplier. The shift adjusts the mult result, so that it
 * is correctly added into the final multiply result.
 */
Bignum mult(Bignum const &x, uint32_t num, unsigned shft) {
  unsigned top = 0;
  unsigned size = x.size + shft + 1; // +1 carry
  unsigned sign = 0;
  unsigned *tab;
  uint64_t coeff = 0;
  uint64_t carry = 0;

  tab = new unsigned[size];

  for (unsigned i = 0; i < shft; i++) {
    tab[top++] = (unsigned)0;
  }

  for (unsigned j = 0; j < x.top; j++) {
    coeff = ((uint64_t)x.tab[j] * (uint64_t)num) + carry;
    tab[top++] = (uint32_t)coeff; // lower 32 bits
    carry = coeff >> 32;          // upper 32 bits
  }

  if (carry)
    tab[top++] = carry;
  return (Bignum(tab, top, size, sign));
}

/*
 * Implements classic multiply
 * (1) multiplicand is multiplied with each coefficient in the multiplier
 * (2) these results are shifted based on the coefficient place value
 * (3) these small multiplications are added up to give the result
 */
Bignum multiply(Bignum const &x, Bignum const &y) {

  /* recursive call to insure x > y */
  if (x < y)
    return multiply(y, x);

  Bignum sum_of_prod(0);
  unsigned rows = y.top;

  /* start with msb of multiplier and work backwards */
  for (unsigned r = rows; r > 0; r--) {
    sum_of_prod = sum_of_prod + mult(x, y.tab[r - 1], r - 1);
  }

  sum_of_prod.sign = (x.sign == y.sign) ? 0 : 1;
  return (sum_of_prod);
}

/*
 The divide Algo is called "multiple-precision division", it gives
 both the quotient and the remainder
 Here is the algorithm as outlined in : Handbook of Applied Cryptography
 1. For j from 0 to (n-t) inclusive, do: qj <- 0

 2. While [x >= y*b^(n-t)] do
      q@(n-t) <- q@(n-t) + 1
            x <- x - y*b^(n-t)

 3. for (i = n; i >= (t+1); i--) // inclusive
    3.1 q@(i-t-1) = (If xi == yt) ? q@(i-t-1) ? b-1 : [xi*b + x@(i-1)]//yt
        whereby b-1 = 2^32 - 1 = max value represented in 32bit space
        whereby [xi*b + x@(i-1)] is formed by glueing the two 32bit values
        into one uint64_t
        whereby the floor function is handled by the integer division
        i.e. integer division drops the decimal portion
    3.2 While ( q@(i-t-1)[yt * b + y@(t-1)] >
                  xi * b^2 + x@(i-1) * b + x@(i-2) )
          q@(i-t-1) <-  q@(i-t-1) - 1
    3.3 x <- x - q@a(i-t-1) * y * b@(i-t-1)
    3.4 If x < 0
        x <- x + y * b@(i-t-1)
        q@(i-t-1) <- q@(i-t-1) - 1

 Implement an alternative in steps 3.3/3.4
    in unsigned math can detect if
    x<0 by looking at msb of the x

    change the code to implement the following
    if x < q@a(i-t-1) * y * b@(i-t-1)
      execute 3.4
      x <- x - q@a(i-t-1) * y * b@(i-t-1) + y * b@(i-t-1)
    else
      execute 3.3

 Note: n = number of coefficients in x, first coeff is n0
       t = number of coefficients in y, first coeff is t0
     n-t = number of coefficients in

 Special case: if t=0, this algo does not work and a fixup is required
 (1) change to t=1 and left shift dividend, divisor, and remainder by 32 (1 blk)
     then at the end do a fixup and divide the q and r by 32
 or
 (2) access [t-1] elem only if t!=0

 Optimization: in step (1) instead of blindly loading q with zeros
 we load it with the init value x//y, this speeds up the division since
 the quotient starts closer to the result rather than zeros.
 init q with max[(x//y)-1, 0]
*/
pair<Bignum, Bignum> division_modulo(Bignum const &x, Bignum const &y) {
  unsigned q_top;
  unsigned q_size = x.top - y.top + 1;
  unsigned q_sign = 0;
  unsigned *q_tab;
  bool r_fixup = false; // if t need fixup, so does r

  unsigned n = x.top - 1;
  unsigned t = y.top - 1;

  Bignum x_cpy = x, y_cpy = y;

  /* For n<t, this case not handled by the algo, result: q=0, r=x */
  if (n < t) {
    Bignum q(0);
    pair<Bignum, Bignum> result(q, x_cpy);
    return result;
  }

  x_cpy.setSignPos(); // force unsigned div
  y_cpy.setSignPos();

  /*algo not work for t=0, change it to t=1, shift operand & divide by 32*/
  /*
  if (t == 0) {
    x_cpy = x << 32;
    y_cpy = y << 32;
    n++;
    t++;
    r_fixup = true;
  }
  */

  q_tab = new unsigned[q_size];

  /*********** formulate y*b^(n-t)  ************/
  Bignum yb(y_cpy);
  unsigned q_msb = 0;

  for (unsigned k = 0; k < (n - t); k++)
    yb <<= 32;

  /* optimization [(xn/yt) - 1] hint to the dividend */
  if (x_cpy.tab[n] > y_cpy.tab[t])
    q_msb = x_cpy.tab[n] / y_cpy.tab[t]; // need to calculate (floor - 1)

  /********************************  Step 1.   ********************/
  for (unsigned j = 0; j <= (n - t); j++) {
    q_tab[j] = 0;
    if (j == (n - t)) {
      q_tab[j] = q_msb;
      x_cpy = x_cpy - (yb * q_msb);
    }
  }

  /*********************************  Step 2.    ******************/
  while (x_cpy >= yb) {
    q_tab[n - t] = q_tab[n - t] + 1;
    x_cpy = x_cpy - yb;
  }

  /*********************************  Step 3.    ******************/
  uint64_t x_long, y_long, floor_x_y;
  uint32_t x_upper, x_lower, floor_lower;
  Bignum x_x_x(1);
  Bignum q_y_y(1);

  for (unsigned i = n; i >= (t + 1); i--) {

    /****************** formulate [xi*b + x@(i-1)]    *************/
    x_upper = x_cpy.tab[i];
    x_lower = x_cpy.tab[i - 1];

    x_long = 0; // insure upper 32b are cleared
    x_long = x_upper;
    x_long <<= 32;     // wr to the upper register
    x_long += x_lower; // wr to the lower register

    /********************* formulate floor term *******************/
    y_long = 0;
    y_long = y_cpy.tab[t];
    floor_x_y = x_long / y_long;
    floor_lower = (uint32_t)floor_x_y; // get lower 32 bits

    /************************* execute 3.1    *********************/
    q_tab[i - t - 1] = (x_cpy.tab[i] == y_cpy.tab[t]) ? UINT_MAX : floor_lower;

    /*****************  q@(i-t-1)[yt * b + y@(t-1)] ***************/
    Bignum q_i_t_minus_1(q_tab[i - t - 1]);
    Bignum y_t(y_cpy.tab[t]);

    /* FIX for t=0, do this operation only if t-1 place value exist*/
    if (t != 0) {
      y_t <<= 32; // move to next place value, ie *b
      y_t = y_t + y_cpy.tab[t - 1];
    }

    q_y_y = q_i_t_minus_1 * y_t;

    /***************** xi * b^2 + x@(i-1) * b + x@(i-2) ***********/
    Bignum x_x_x(x_cpy.tab[i]);
    x_x_x <<= 32;
    x_x_x = x_x_x + x_cpy.tab[i - 1];

    /*FIX FOR ALGO: t=0, y op not have t-1 place value, fix x as well */
    if (t != 0) {
      x_x_x <<= 32;
      x_x_x = x_x_x + x_cpy.tab[i - 2];
    }

    /************************** execute 3.2  **********************/
    while (q_y_y > x_x_x) {
      q_tab[i - t - 1] = q_tab[i - t - 1] - 1;

      q_i_t_minus_1 = q_tab[i - t - 1]; // update q_y_y for the while cond
      q_y_y = q_i_t_minus_1 * y_t;
    }

    /* Alternative 3.3/3.4 */
    /******************* formulate y*b^(i-t-1)  *******************/
    yb = y_cpy; // reuse this var
    for (unsigned k = 0; k < (i - t - 1); k++) {
      yb = yb << 32;
    }

    /************************* execute 3.3  ***********************/
    Bignum q_yb(1);
    q_yb = q_tab[i - t - 1] * yb;
    if (x_cpy < q_yb) {
      x_cpy = x_cpy - q_yb + yb; // x = x - y(b@i-t-1)[(q@i-t-1) - 1]
      q_tab[i - t - 1] = q_tab[i - t - 1] - 1;
    } else {
      x_cpy = x_cpy - q_yb;
    }
  }

  q_top = q_size;
  while ((x_cpy.top > 1) && (x_cpy.tab[x_cpy.top - 1] == 0)) // find real top
    x_cpy.top--;
  while ((q_top > 1) && (q_tab[q_top - 1] == 0)) // find real top
    q_top--;

  if (r_fixup)
    x_cpy = x_cpy >> 32;

  Bignum quotient = Bignum(q_tab, q_top, q_size, q_sign);
  Bignum r = x_cpy;

  pair<Bignum, Bignum> result(quotient, r);

  return result;
}

/* take the quotient from the division_module fnc  */
Bignum operator/(Bignum const &x, Bignum const &y) {

  pair<Bignum, Bignum> result = division_modulo(x, y);

  return result.first;
}

/* take the remainder from the division_module fnc */
Bignum operator%(Bignum const &x, Bignum const &y) {

  pair<Bignum, Bignum> result = division_modulo(x, y);

  return result.second;
}

/* For Debug: Displays the coefficients of the BN seperately */
void Bignum::showInfo() {
  char mysign = '+';
  if (this->sign == 1)
    mysign = '-';

  cout << "top:" << this->top << " size:" << this->size;
  cout << " val: (" << mysign << ") ";

  for (unsigned i = 0; i < top; i++) { // Little Endian
    if ((i + 1) == top)
      cout << tab[i] << "*(2^" << 32 * (i) << ") ";
    else
      cout << tab[i] << "*(2^" << 32 * (i) << ") + ";
  }

  cout << "\n";
}

ostream &operator<<(ostream &os, const Bignum &x) {

  ostringstream sstream;
  string result;
  // pad_len = 0;

  (x.sign == 0) ? os << "+" : os << "- ";

  for (unsigned i = x.top; i > 0; i--) { // big endian
    if (i == x.top)
      os << "0x";

    sstream << setfill('0') << setw(8) << hex << x.tab[i - 1];
    result = sstream.str();
    sstream.str(string()); // flush the sstreamn for next iter

    if (i - 1 == 0)
      os << result;
    else
      os << result << "_";
  }
  return os;
}

/*
 * Divide by 10 Algo:
 * q = q * .8/8 = q * .1 = q/10
 */
pair<Bignum, Bignum> division_modulo_10(Bignum const &x) {

  Bignum q(1), r(1);
  q = (x >> 1) + (x >> 2);       // n * .75                  = n*.75
  q = q + (q >> 4);              // n * 51/64                = n*.797
  q = q + (q >> 8);              // n * 13107/16384          = n*.799
  q = q + (q >> 16);             // n * 858993458/1073741824 = n*.8
  q = q >> 3;                    //(n *.8)/ 8                = n/10
  r = x - (((q << 2) + q) << 1); // n- 2*[4*n/10 + n/10]  = n - 10n/10

  r = x - q * 10; // calculate r accurately

  /* the error term is too big, make it smaller */
  if (r > MAX_DIV_ERR) {
    pair<Bignum, Bignum> result = division_modulo_10(r);
    q = q + result.first;
    r = result.second;
  } else {
    Bignum fixup = r / 10; // floor
    if (r.isNeg()) {
      fixup.setSignNeg();
    }

    q = q + fixup;
    r = r - (fixup * 10);
  }

  if (r >= 10)
    cout << "Bad division, remainder ge 10; r = " << r << "\n";

  pair<Bignum, Bignum> result(q, r);
  return result;
}

/*
 * Convert the input to diplay in base10, by using the fast
 * divide-by-10 algo
 */
string convertInBase10(Bignum x) {
  Bignum r(1);

  if (x >= 10) {

    pair<Bignum, Bignum> result = division_modulo_10(x);
    x = result.first;
    r = result.second;
    return convertInBase10(x) + std::to_string(r.tab[0]);
  }

  char sign = x.isNeg() ? '-' : '+';
  pair<Bignum, Bignum> result = division_modulo_10(x);
  return sign + std::to_string(result.second.tab[0]); // called once
}

/* same as above, except use the generic divide function */
string displayBase10(Bignum x) {

  Bignum r(1);

  if (x >= 10) {
    pair<Bignum, Bignum> result = division_modulo(x, 10);
    x = result.first;
    r = result.second;
    return displayBase10(x) + std::to_string(r.tab[0]);
  }

  char sign = x.isNeg() ? '-' : '+';
  pair<Bignum, Bignum> result = division_modulo(x, 10);
  return sign + std::to_string(result.second.tab[0]); // called once
}

Bignum generate_random_bignum(unsigned nb_of_bits) {

  Bignum random_bignum(0);

  if (!seed_intialized) {
    srand(time(NULL));
    seed_intialized = true;
  }
  unsigned numberOfBlocks = nb_of_bits / 32;
  unsigned *tab = new unsigned[numberOfBlocks];

  for (unsigned i = 0; i < numberOfBlocks; i++) {

    unsigned random_number = rand();
    tab[i] = random_number;
  }

  return Bignum(tab, numberOfBlocks, numberOfBlocks, 0);
}

/**
 * Return True if the given bignum is a prime, False otherwise.
 * Checks if the bignum is a prime or not using Fermat's little theorem.
 */
static bool is_a_prime_bignum(Bignum x) {

  unsigned tab[4] = {2, 3, 5, 7};

  for (unsigned i = 0; i < 4; i++) {

    if (pow(Bignum(tab[i]), x, x) == Bignum(tab[i])) {

    } else {
      return false;
    }
  }
  return true;
}

Bignum generate_prime_bignum(unsigned nb_of_bits) {
  int iterations = 0;
  bool is_finished = false;

  Bignum x(0);
  while (!is_finished) {
    x = generate_random_bignum(nb_of_bits);

    if (x % 2 == Bignum(0)) {
      x = x + 1;
    }
    is_finished = is_a_prime_bignum(x);

    if (!is_finished)
      iterations++;
  }

  return x;
}

string convertToBinary(Bignum x) {

  char binary = convertInBase10(x % 2)[1];
  string binary_string;
  binary_string += binary;

  x = x >> 1;

  return x == Bignum(0) ? binary_string : convertToBinary(x) + binary_string;
}

/**
 * Compute x^y(mod n) using the square and multiply algorithm.
 */
Bignum square_and_multiply_algorithm(Bignum x, Bignum y, Bignum n) {

  Bignum b(1);
  string y_binary = convertToBinary(y);

  for (int i = 0; i < y_binary.length(); i++) {
    b = (b * b) % n;
    if (y_binary[i] == '1') {
      b = (b * x) % n;
    }
  }

  return b;
}

/**
 * Compute x^y(mod n) using the square and multiply algorithm.
 */
Bignum pow(Bignum x, Bignum y, Bignum n) {
  return square_and_multiply_algorithm(x, y, n);
}

/* A*x + M*y = 1 calculate x and y terms, x is the modular inverse */
tuple<Bignum, Bignum, Bignum> gcdExtended(Bignum const &A, Bignum const &M) {

  if (M == 0) {
    tuple<Bignum, Bignum, Bignum> gcd_result(A, 1, 0); // return {a, 1, 0}
    return gcd_result;

  } else {
    Bignum gcd(1), p0(1), p1(1);
    tie(gcd, p1, p0) = gcdExtended(M, Bignum(A % M));

    p1 = p1 - (A / M) * p0;

    tuple<Bignum, Bignum, Bignum> gcd_result(gcd, p0, p1);

    return gcd_result;
  }
}

/*
 *  if gcd(e,n)  = 1 return the modular inverse
 *  if gcd(e,n) != 1 return 0, modular inverse not exist
 *  given e and n, modular inverse solves for d:
 *  e*d = 1 mod n
 *
 */
Bignum modInverse(Bignum const &e, Bignum const &n) {
  Bignum gcd(1), x(1), y(1);

  tie(gcd, x, y) = gcdExtended(e, n);

  if (gcd == 1) {
    if (x.isNeg())
      x = n + x;
    return x;

  } else {
    return 0;
  }
}

/*
 * Low-level Primality Test : check if the candidate is
 * divisible by a small set of prime numbers
 */
bool is_a_prime_pass1(Bignum x) {

  unsigned tab[5] = {2, 3, 5, 7, 11};

  for (unsigned i = 0; i < 5; i++) {
    if (x % tab[i] == 0)
      return false;
  }

  return true;
}

/* Prevent overflow of check_composite function */
static Bignum mul(Bignum a, Bignum b, Bignum mod) {
  Bignum res(0);

  while (b != 0) {
    if ((b % 2) == 1)
      res = (res + a) % mod;
    a = (a << 1) % mod;
    b >>= 1;
  }

  return res;
}

/*
 * Verify that the prime candidate is not a composite, ie it only has
 * factors 1 and itself.
 */
static bool check_composite(Bignum p_cand, Bignum p_base, Bignum p_pow,
                            Bignum two) {

  Bignum x = pow(p_base, p_pow, p_cand); //(p_base^p_pow) % p_cand

  if ((x == 1) || (x == (p_cand - 1)))
    return false;

  for (int k = 1; k < two; k++) {
    x = mul(x, x, p_cand);

    if (x == (p_cand - 1))
      return false;
  }

  return true;
}

/*
 * High-level Primality Test: Rabin Miller Primality Test
 */
bool is_a_prime_pass2(Bignum x) {

  Bignum tab[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

  Bignum s(0);
  Bignum d = x - 1;

  while ((d % 2) == 0) { // check lsb until 1 is found
    d >>= 1;
    s = s + 1;
  }

  for (int i = 0; i < 12; i++) {
    if (check_composite(x, tab[i], d, s))
      return false;
  }

  return true;
}

/*
 * Takes a "prime candidate" and runs a low-level primality test and
 * then a high-level primality test.
 */
Bignum generate_prime_bignum2(unsigned nb_of_bits) {
  int iterations = 0;
  bool is_finished = false;
  uint32_t nb_minus1 = nb_of_bits - 1;
  uint32_t nb_minus2 = nb_of_bits - 2;

  Bignum x(0);
  Bignum value_one(1);

  while (!is_finished) {
    x = generate_random_bignum(nb_of_bits);

    if ((x % 2) == 0)
      x = x + 1;

    if ((x >> nb_minus2) == 0)
      x = x + (value_one << nb_minus2);

    value_one = 1;
    if ((x >> nb_minus1) == 0)
      x = x + (value_one << nb_minus1);

    is_finished = is_a_prime_pass1(x) && is_a_prime_pass2(x);
    if (!is_finished) {
      iterations++;
    }
  }
  cout << "iterations to get prime number: " << iterations << "\n";
  if (is_a_prime_bignum(x))
    cout << "this prime passes Abdul's test!!\n";
  else
    cout << "this prime DOES NOT PASS Abdul's test!!\n";

  return x;
}

/*
 * Implement the RSA algo as follows
 * (1) request user for a 16 byte secret message (16 chars)
 * (2) Generate 2 prime numbers, p and q, where the product p*q is
 *     approximately the size of the required bit length (1024 bit=32 blks)
 *     i.e. generate primes with p_len = |k/2| = 512
 * (3) choose e = 65537 as it has two bits of 1 in the binary
 *     representation making for a faster modular exponentiation operation
 *     e = 1_0000_0000_0000_0001 = 0x10001
 * (4) Euler's totient is calculated as : phi = (p-1) * (q-1)
 * (5) modular n value is calculated as : n = pq
 * (6) inverse modular to find the private key d:
 *     e*d =1 mod phi or [d^-1]mod phi = 1
 * (7) public key (e,n) ; private key (d,n) used to encrypt/decrypt where
 *     pt = plaintext and ct = ciphertext
 *     Encryption(pt,e,n): ct = (pt^e) mod n
 *     Decryption(ct,d,n): pt = (ct^d) mod n
 */
void RSA(string msg, unsigned nb_of_bits_prime_numbers) {

  /* msg goes into a bignum, take each char and store as an int */
  Bignum pt(0);
  vector<char> msg_bytes_v(msg.begin(), msg.end());
  int last_elem = msg_bytes_v.size() - 1;

  for (int i = 0; i < msg_bytes_v.size(); i++) {
    pt += msg_bytes_v[i];
    if (i != last_elem)
      pt <<= CHAR_BITS;
  }

  cout << "-------------- Welcome! Text to encrypt and decrypt: ";
  cout << "\n-------------- \"" << msg << "\"\n";

  cout << "************ Let's generate RSA parameters ************\n";
  cout << ">> Generating prime numbers may take a few seconds!<<\n";

  cout << "+ Generating p ..." << flush;
  Bignum p = generate_prime_bignum(nb_of_bits_prime_numbers);
  cout << "OK!\n";
  cout << "+ Generating q ..." << flush;
  Bignum q = generate_prime_bignum(nb_of_bits_prime_numbers);
  cout << "OK!\n";
  Bignum e(65537);
  cout << "+ Selecting e ...OK!\n";

  cout << "+ Computing N = p * q ..." << flush;
  Bignum n = p * q;
  cout << "OK!\n";
  cout << "+ Computing phi = (p-1) * (q-1) ..." << flush;
  Bignum phi = (p - 1) * (q - 1); // Euler's totient
  cout << "OK!\n";

  cout << "+ Computing d such that e * d ≡ 1 mod [phi]  ..." << flush;
  Bignum d = modInverse(e, phi);
  cout << "OK!\n";

  cout << "************ RSA parameters ***********\n";

  cout << "p=  " << convertInBase10(p) << "\n\n";
  cout << "q=  " << convertInBase10(q) << "\n\n";
  cout << "e=  " << convertInBase10(e) << "\n\n";
  cout << "N=  " << convertInBase10(n) << "\n\n";
  cout << "phi= " << convertInBase10(phi) << "\n\n";
  cout << "d=  " << convertInBase10(d) << "\n\n";

  cout << "************ Last step ***********\n";

  cout << "+ Converting "
       << "\"" << msg << "\""
       << " to bignum ...OK!:\n"
       << "pt= " << convertInBase10(pt) << "\n\n";

  /*
   * ct = pt^e mod n
   * pt = ct^d mod n
   */
  Bignum ct(1), pt_(1);
  ct = Bignum(pow(pt, e, n));
  pt_ = Bignum(pow(ct, d, n));

  cout << "+ Encrypting the msg...";

  cout << "OK!:\nEncrypted msg: ct= " << convertInBase10(ct) << "\n\n";

  cout << "+ Decrypting the msg...OK!:\n";

  cout << "Decrypted msg: pt_= " << convertInBase10(pt_) << "\n\n";

  cout << "+ Converting pt_ to ascii...OK!:\n";

  cout << "Decrypted msg: \"" << convertPTtoAscii(pt_) << "\"\n";
}

string convertPTtoAscii(Bignum const &pt) {
  uint64_t binary_length = 0, num_chars_pt = 0, msb_char_bits = 0;
  Bignum pt_cpy = pt;

  while (pt_cpy != 0) {
    pt_cpy >>= 1;
    binary_length++;
  }

  num_chars_pt = binary_length / CHAR_BITS;
  msb_char_bits = binary_length % CHAR_BITS;
  if (msb_char_bits)
    num_chars_pt++;

  Bignum char_msk(255);
  unsigned lsb_char;
  unsigned char_count = 0;

  pt_cpy = pt;
  vector<char> pt_v(num_chars_pt, '-');

  while (num_chars_pt != 0) {

    lsb_char = (pt_cpy & char_msk).getTab()[0];
    pt_v[char_count] = lsb_char;
    if (num_chars_pt != 1)
      pt_cpy >>= 8;
    num_chars_pt--;
    char_count++;
  }

  reverse(pt_v.begin(), pt_v.end());

  string ascii_msg = "";

  for (int x : pt_v)
    ascii_msg += (char)x;

  return ascii_msg;
}

// Function for extended Euclidean Algorithm
int gcdExtended_int(int e, int n, int *x, int *y) {
  // Base Case
  if (e == 0) {
    *x = 0;
    *y = 1;
    return n;
  }

  // store results of recursive call
  int x1, y1;
  int gcd = gcdExtended_int(n % e, e, &x1, &y1);

  // Update x and y using results of recursive call
  *x = y1 - (n / e) * x1;
  *y = x1;

  return gcd;
}

void modInverse_int(int e, int n) {
  int x, y, res = 0;

  int g = gcdExtended_int(e, n, &x, &y);
  cout << "g= " << g << " x= " << x << " e= " << e << " n = " << n << "\n";
  if (g != 1)
    cout << "Inverse does not exist for (e^-1)mod n\n";
  else {
    // res = (x % n + n) % n; //+n for the case of neg x
    if (x < 0) {
      res = (x + n);
    } else {
      res = x;
    }
    cout << "Modular multiplicative inverse is " << res << "\n";
  }
}
