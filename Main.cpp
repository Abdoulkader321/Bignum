/*
 * Main.cpp
 *
 *  Created on: 7 nov. 2022
 *  Maria Ukanwa
 *  Abdoulkader Moussa
 *
 * Compile:
 *  + make
 *  + ./Main
 *
 */

#include "Bignum.hpp"

//#include "math.h"
//#include <limits.h>
//#include <cstdint>
//#include <iomanip>
//#include <iostream>
//#include <sstream>
#include <climits>
#include <vector>


void TEST(string expected, string got) {
  cout << "EXPECTED: " << expected;
  cout << " GOT: " << got;
  bool is_equal = expected == got;
  cout << (is_equal ? " (passed)" : " (failed!)") << "\n";
}

void TEST_MOVE() {
  Bignum b1(UINT_MAX - 1);
  cout << "b1= " << b1 << '\n';
  Bignum b2(b1);        // testing deep copy
  
  if ((b2 == b1) && (b2.getSign() == b1.getSign()))
    cout << "deep copy passed\n";
  else
    cout << "deep copy failed\n";
  
  Bignum b3 = move(b1); // testing the move, b11 no longer accessible
  cout << "b2= " << b2 << "\nb3= " << b3 << '\n';

}

int main() {

  TEST_MOVE();
  cout << "\n\n***********TESTING ADDITION*************\n";
  
  Bignum b1(3);
  Bignum b2(-5);
  Bignum b4(UINT_MAX); // 4 294 967 295
  Bignum b5(UINT_MAX);
  Bignum b6(5);

  b1-=b4;
  TEST("-4294967292", convertInBase10(b1));
  
  Bignum b44(UINT_MAX);
  b44 += b44;
  TEST("+8589934590", convertInBase10(b44));
  
  b44 -= b44;
  TEST("+0", convertInBase10(b44));
  TEST("+4294967298", convertInBase10(Bignum(3)+Bignum(UINT_MAX)));

  Bignum b7(0xFFFFFFFC); b7.setSignNeg();
  TEST("-4294967297", convertInBase10(b7+Bignum(-5)));
  
  TEST("+8589934590", convertInBase10(b4+b4));
  Bignum b3 = b4 + b4;

  b3 += b3;
  TEST("+17179869180", convertInBase10(b3));
  
  // test recursive call, switch ops
  // 4294967295 + 17179869180
  TEST("+21474836475", convertInBase10(b4 + b3));
  
  TEST("-2", convertInBase10(Bignum(3)+Bignum(-5)));
  
  cout << "\n\n***********TESTING SUBTRACTION*************\n";
  Bignum s1((uint32_t)1 << 31);
  Bignum s2((uint32_t)1 << 31);
  Bignum s3(1);

  s1 += s1;
  for (size_t i = 0; i < 5; i++) {
    s1 = s1 * s1;
  }
  
  for (size_t i = 0; i < 1; i++) {
    s2 = s2 * s2;
  }

TEST("-539307940458694772318791557236707420085393093682691971820290243473198027416502889398125431967222608063360341639614180072976369306443249867478542291918422373133303680274596455828906658803738282358359248856255017306514452047027388644421739331622481711490051532053758894719841737815439148914506068988872672411648", convertInBase10(s1*Bignum(-3)));
  
  s2 -= s1;
  TEST("-179769313486231590772930519078902473361797697894230657273430081157732675805500963132708477322407536021120113879871393357658789768814416622492847430639474124377767893424865485276302219601246094119453082952085005768838150682342462881473913110540827237163350510684586298239947245938479716304830744643605796749312", convertInBase10(s2));

  cout << "\n\n***********TESTING ADDITION*************\n";
  Bignum m1((uint32_t)1 << 31);
  cout << m1 << "\n";
  for (size_t i = 0; i < 10; i++) {
    m1 += m1;
  }

  TEST("+2199023255552", convertInBase10(m1));
  
  cout << "\n\n***********TESTING MULTIPLY*************\n";
  TEST("+21474836475", convertInBase10(b4*b6));
  Bignum m2((uint32_t)1 << 31);
  m2 += m2;
  for (size_t i = 0; i < 4; i++) {
    m2 = m2 * m2;
  }

  TEST("+13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084096", convertInBase10(m2));
  
  cout << "\n\n********TESTING LEFT/RIGHT SHFT *********\n";
  Bignum e11(0xfdbc593b);
  Bignum e22(0xffff000f);
  Bignum e33(0x9ace369b);
  uint32_t shfte = 32;

  e11 = e11 << shfte;
  e11 = e11 + e22;
  e11 = e11 << shfte;
  e11 = e11 + e33;
  e11 = e11 << 65; // this is the real test

  TEST("+2897149957052279815278654728864757670563488464896", convertInBase10(e11));
  
  cout << "\n\n********MORE LEFT/RIGHT SHFT *********\n";
  Bignum d1((uint32_t)1 << 31);
  Bignum d22(0xffff000f);
  Bignum d33(0x12345678);
  Bignum d55(0xababcdcd);

  uint32_t shft = 32;

  d1 = d1 << shft;
  d1 = d1 + d22;
  d1 = d1 << shft;
  d1 = d1 + d33;
  d1 = d1 >> shft;
  d1 = d1 << 50;
  d1 = d1 + d33;
  d1 = d1 << 35;
  d1 = d1 + d55;
  d1 = d1 << 31;
  d1 = d1 + d33;

  TEST("+766247770789750909083944017408679933334397452203349624", convertInBase10(d1));
  
  d1 = d1 >> 64;

  TEST("+41538374887621139061726853708470752", convertInBase10(d1));
  
  d1 = d1 >> 28;

  TEST("+154742504982729029140348932", convertInBase10(d1));
  
  d1 = d22;
  d1 = d1 >> 31;
  TEST("+1", convertInBase10(d1));
  
  cout << "\n\n************ TESTING DIVISON **************\n";
  Bignum d2((uint32_t)1 << 31);
  Bignum d3((uint32_t)1 << 31);
  Bignum d4((uint32_t)1 << 31);
  Bignum d5((uint32_t)1 << 31);
  Bignum d6((uint32_t)1 << 31);
  d2 = d2 << shft;
  d2 = d2 << shft;
  d2 = d2 << shft;
  d3 = d3 << shft;
  d3 = d3 << 28;

  cout << "d2 " << d2 << "\n";
  cout << "d3 " << d3 << "\n";
  cout << "d2/d3\n";
  
  d2 = d2 / d3;
  TEST("+68719476736", convertInBase10(d2));

  d4 = d4 << shft;
  d4 = d4 << shft;
  d4 = d4 << shft;
  d5 = d5 << shft;
  d5 = d5 + d6;

  cout << "d4 " << d4 << "\n";
  cout << "d5 " << d5 << "\n";
  cout << "d4/d5\n";
  
  d4 = d4 / d5;
  
  TEST("+18446744069414584320", convertInBase10(d4));
  
  Bignum a40(0x00000028);  // int = 40
  a40 = a40 << shft;
  Bignum a23(0x00000017);  // int = 23
  a23 = a23 << shft;
  a23 = a23 + Bignum(0x30000000);
  cout << "a40: " << a40 << "\n";
  cout << "a23: " << a23 << "\n";  //DEBUG, CORRECT ANS, WHY DIV BUG NOT KICK IN
  a40 = a40 / a23;   // 0x28_00000000 div 0x17_30000000
  cout << "a40/a23\n";

  TEST("+1", convertInBase10(a40));
  
  cout << "\n********TESTING Add op *********\n";
  // trace e1 - e2
  Bignum e1(0xfabc1234);
  Bignum e2(0xff32dd32);
  Bignum e3(0x12345678);
  Bignum a2(0x0476cafd); // e1 - e2:   0x0476cafd_fb893502_00000000

  e1 = e1 << shft;
  e1 = e1 + e2;
  e1 = e1 << shft;
  e1 = e1 + e3; // e1 = 0xfabc1234_ff32dd32_12345678

  e2 = e2 << shft;
  e2 = e2 + 0xfabc1234;
  e2 = e2 << shft;
  e2 = e2 + e3; // e2 = 0xff32dd32_fabc1234_12345678

  //answer to e1 + e2, check this answer
  a2 = a2 << shft;
  a2 = a2 + 0xfb893502;
  a2 = a2 << shft;

  e2.setSignNeg();
  cout << "e1:" << e1 << "\n";
  cout << "e2:" << e2 << "\n";
  e1 = e1 + e2;

  cout << "e1 + e2 :" << e1 << "\n";
  TEST("-1381551889180773292399656960", convertInBase10(e1));

  cout << "\n********TESTING Shift op *********\n";
  uint32_t my_shft = 65;
  Bignum f1(UINT_MAX);
  f1 = f1 << shft;
  f1 = f1 + 0xfabc1234;
  f1 = f1 << 1; // use 3 boxes

  cout << f1 << "\n";

  f1 = f1 << my_shft;
  cout << f1 << "\n";
  TEST("+1361129467677235669681445883163099267072",convertInBase10(f1));
  
  cout << "\n********TESTING Modulo *********\n";
  cout << "(UINT32_MAX+1) % UINT32_MAX: ";
  Bignum r1(UINT32_MAX);
  Bignum r2(UINT32_MAX);
  r1 += 2;

  Bignum r3(10);
  Bignum r4(9);
  cout << "\n" << r3 << "/" << r4;
  r3 = r3 / r4;
  cout << "= " << r3 << "\n";

  cout << "\n" << r1 << "/" << r2;
  // r1 = r1/r2;  //causes seg fault, prob with div algo
  cout << "= " << r1 << "\n";

  cout << "\n" << r1 << "%" << r2;
  // r1 = r1 % r2; //causes seg fault, prob with div algo
  cout << "= " << r1 << "\n";
  r1.showInfo();

  cout << "\n********TESTING FIX SUBTRACTION BUG *********\n";
  Bignum t1(0x00000003);
  Bignum t2(0x00000003);
  t1 = t1 << 32;
  cout << "t1: " << t1 << "\n";
  t1 = t1 + 0xffffffff;
  t1 = t1 << 32;
  cout << "t1: " << t1 << "\n";
  t1 = t1 + 0xeaf048d0;
  t1 = t1 << 64;
  cout << "t1: " << t1 << "\n";

  t2 = t2 << 32;
  t2 = t2 + 0xfffffffb;
  t2 = t2 << 32;
  t2 = t2 + 0xeaf048d0;
  t2 = t2 << 32;
  t2 = t2 + 0x150fb730;
  t2 = t2 << 32;

  Bignum t3 = t1 - t2;
  cout << "t1: " << t1 << "\n";
  cout << "t2: " << t2 << "\n";
  cout << "t3: " << t3 << "\n";

  cout << "\n********TESTING Display Base10 *********\n";
  cout << f1 << "\n";
  cout << convertInBase10(f1);
  cout << "\n\n" << Bignum(UINT32_MAX) << "\n";

  TEST("+302", displayBase10(Bignum(302)));
  TEST("+4294967295", displayBase10(Bignum(UINT32_MAX)));

  cout << "\n*********** TESTING MULTIPLICATION 2*********\n";

  TEST("+0",
       convertInBase10(Bignum(UINT32_MAX) * Bignum(0))); // Multiply By 0
  
  Bignum A(UINT32_MAX);
  A += 1;
  A = A * Bignum(-1);
  TEST("-4294967296", convertInBase10(A)); // Multiply By -1;

  TEST("+18446744073709551616",
       convertInBase10(A * A)); // Multiply 2 negative numbers

  Bignum res(UINT32_MAX);
  res += 1;

  res = res * Bignum(2);
  TEST("+8589934592", convertInBase10(res));

  res = res * Bignum(1000);
  TEST("+8589934592000", convertInBase10(res));

  cout << "\n*********** TESTING SUBTRACTION 2*********\n";

  A = Bignum(UINT32_MAX);
  A += 1;
  A = A * 3;

  TEST("+0", convertInBase10(A - A));
  TEST("+12872556210", convertInBase10(A - Bignum(12345678)));
  TEST("+11650333997", convertInBase10(A - Bignum(1234567891)));

  // Sign test
  Bignum B(UINT32_MAX);
  B += 1;
  B = B * 5;

  TEST("-8589934592", convertInBase10(A - B));
  TEST("+8589934592", convertInBase10(B - A));

  A = A - B; // Neg
  TEST("+0", convertInBase10(A - A));

  cout << "\n*********** Generate random bignum *********\n";
  unsigned nb_blocks = 32;

  Bignum random1 = generate_random_bignum(nb_blocks);
  cout << "Random1:\n" << convertInBase10(random1) << "\n";

  Bignum random2 = generate_random_bignum(nb_blocks);
  cout << "Random2:\n" << convertInBase10(random2) << "\n";

  cout << "Random1 * Random2:\n" << convertInBase10(random1 * random2) << "\n";

  cout << "\n************* (A^B) % N **************\n";

  A = Bignum(987654321);
  B = Bignum(123456789);
  Bignum N = Bignum(70);

  TEST("+41", convertInBase10(pow(A, B, N)));

  cout << "\n*************** MODULO INVERSE ************\n";
  modInverse_int(15, 26);
  modInverse_int(3, 13);
  modInverse_int(0x0fff6753, 0xffff8971);
  
  TEST("+7", convertInBase10(modInverse(Bignum(15), Bignum(26))));
  TEST("+9", convertInBase10(modInverse(Bignum(3), Bignum(13))));
  TEST("+636088586", convertInBase10(modInverse(Bignum(0x0fff6753),
                                                Bignum(0xffff8971))));
  
  cout << "\n*************** Shift_function ************\n";
  Bignum t11(0x00000003);
  Bignum t22(0x00000003);
  Bignum t_shft(32);
  t11 <<= 32;
  t11 = t11 + 0xffffffff;
  t11 <<= 32;
  t11 = t11 + 0xeaf048d0;
  t11 <<= 64;

  t22 <<= 32;
  t22 = t22 + 0xfffffffb;
  t22 <<= 32;
  t22 = t22 + 0xeaf048d0;
  t22 <<= 32;
  t22 = t22 + 0x150fb730;
  t22 <<= 32;

  Bignum t33 = t11 - t22;
  cout << "t33: " << t33 << "\n";
  cout << "t3:  " << t3 << "\n";
  
  cout << "t11: " << t11 << "\n";
  cout << "t22: " << t22 << "\n";
  if (t22 < t11) {
    cout<< "t22 is lt t11" << "\n";
  }
  
  cout << "\n*************** Compare *******************\n";
  if (t11 > t22) {
    cout<< "t11 is gt t22" << "\n";
  }
  
  if (t11 == t11) {
    cout<< "t11 is eq t11" << "\n";
  }
  
  cout << "\n*************** Modulo ********************\n";
  Bignum T = Bignum(1) * Bignum(1);
  T = T % Bignum(5);
  cout << "T= " << T % Bignum(5) << "\n";
  
  cout << "\n************* RSA example ******************\n";
  
  string msg;
  cout << "Please enter a secret message:\n";
  getline(cin, msg);
  
  RSA(msg, 512);
  
  return 0;
}
