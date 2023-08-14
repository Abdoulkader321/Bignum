/*
 * Main.cpp
 *
 *  Created on: 25 dec. 2022
 *  Maria Ukanwa
 *  Abdoulkader Moussa
 *
 * Compile:
 *  + make Rsa
 *  + ./Rsa 
 *
 */

#include "Bignum.hpp"

int main(int argc, char *argv[]) {

  if (argc != 3) {
    const string usage_msg = "Usage:    ./Rsa <msg> <key_nb_of_bits>\n"
                             "Example:  ./Rsa SecretEncryptRSA 512\n\n"
                             "<msg>: the message to encrypt and decrypt\n"
                             "<key_nb_of_bits>: the number of bits of prime numbers\n";
    cout << usage_msg;
    exit(EXIT_FAILURE);
  }

  unsigned nb_of_bits = atoi(argv[2]);
  string msg = argv[1];

  RSA(msg, nb_of_bits);
}