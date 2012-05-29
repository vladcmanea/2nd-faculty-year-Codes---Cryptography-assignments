#ifndef _DES
#define _DES

#include <ctime>
#include <NTL/ZZ.h>
#include <fstream>
#include <iostream>
using namespace std;

NTL_CLIENT

// RSA class
class RSA
{
	// computes x ^ p mod n
	ZZ RSA::ModPower(const ZZ &x, const ZZ &p, const ZZ &n);

	// finds if l, d can be the l and d in e * d - 1 = l * phi(n)
	bool criteria(ZZ l, ZZ d, ZZ n, ZZ e);

public:

	// RSA constructors
	RSA::RSA(istream &fp, istream &fs, bool copeWiener = 1);
	RSA::RSA(ostream &fp, ostream &fs, bool copeWiener = 1);

	// encrypts to fa
	void applyEncrypt(istream &fc, istream &fp, ostream &fa);

	// decrypts to fa
	void applyDecrypt(istream &fc, istream &fp, istream &fs, ostream &fa);

	// applies Wiener decryption
	void applyWienerDecrypt(istream &fc, istream &fp, ostream &fa);

	// divides into factors
	void applyFactors(istream &fp, istream &fs, ostream &fa);

	// RSA destructor
	~RSA();
};

#endif
