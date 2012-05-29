#include "DES.h"
#include <algorithm>

using namespace std;

// DES constructor
DES::DES()
{
	// calculate pow normal
	pow64[0] = 1;
	for (int i = 1; i < 64; ++i)
		pow64[i] = pow64[i-1] << 1;

	// calculez pow64[64]-1
	for (int i = 0; i < 64; ++i)
		pow64[64] |= pow64[i];

	// open P0 file
	ifstream fp0("P0.txt");

	// read P0 from P0 file
	for (int i = 1; i <= 64; ++i)
		if (!(fp0 >> P0[i]))
			throw("[DES::DES] Could not read P0");

	// close P0 file
	fp0.close();

	// calculate P0i (P0's inverse)
	for (int i = 1; i <= 64; ++i)
		P0i[P0[i]] = i;

	// open E file
	ifstream fe("E.txt");

	// read E from E file
	for (int i = 1; i <= 48; ++i)
		if (!(fe >> E[i]))
			throw("[DES::DES] Could not read E");

	// close E file
	fe.close();

	// open P file
	ifstream fp("P.txt");

	// read P from P file
	for (int i = 1; i <= 32; ++i)
		if (!(fp >> P[i]))
			throw("[DES::DES] Could not read P");

	// calculate Pi (P's inverse)
	for (int i = 1; i <= 32; ++i)
		Pi[P[i]] = i;

	// close P file
	fp.close();

	// open M file
	ifstream fm("M.txt");

	// read M from M file
	for (int k = 1; k <= 8; ++k)
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 16; ++j)
				if (!(fm >> M[k][i][j]))
					throw("[DES::DES] Could not read M");

	// close M file
	fm.close();

	// open PC1 file
	ifstream fpc1("PC1.txt");

	// read PC1 from PC1 file
	for (int i = 1; i <= 56; ++i)
		if (!(fpc1 >> PC1[i]))
			throw("[DES::DES] Could not read PC1");

	// close PC1 file
	fpc1.close();

	// open PC2 file
	ifstream fpc2("PC2.txt");

	// read PC2 from PC2 file
	for (int i = 1; i <= 48; ++i)
		if (!(fpc2 >> PC2[i]))
			throw("[DES::DES] Could not read PC2");

	// close PC2 file
	fpc2.close();
}

// prepares a brute text
void DES::prepare(ifstream &fb, ofstream &fp)
{
	char c;
	int64 n;

	// read from brute file
	while (fb.get(c))
	{
		// n is a uint64 number
		n = c;

		// write number n
		if (!(fp << n))
			throw("[DES::prepare] Could not write prepared");

		if (!(fp << ' '))
			throw("[DES::prepare] Could not write prepared");
	}
}

// unprepares to a brute text
void DES::unprepare(ifstream &fd, ofstream &fu)
{
	char c;
	int64 n;

	// read from brute file
	while (fd >> n)
	{
		// decode character c
		c = char(n);

		// write character c
		if (!(fu << c))
				throw("[DES::unprepare] Could not write unprepared");
	}
}


// applies P0 to an unsigned long argument
int64 DES::applyP0(int power, int64 x)
{
	// ans will be returned
	int64 ans = 0;

	// create bits in ans
	for (int i = 1, j = 64; i <= 64 && power == 1; ++i, --j)
	{
		// ans is added with the bit at the position
		ans |= ((x & int64(pow64[64 - P0[i]])) >> (64 - P0[i])) << (j - 1);
	}

	// create bits in ans
	for (int i = 1, j = 64; i <= 64 && power == -1; ++i, --j)
	{
		// ans is added with the bit at the position
		ans |= ((x & int64(pow64[64 - P0i[i]])) >> (64 - P0i[i])) << (j - 1);
	}

	// return result in ans
	return ans;
}

// applies REV to x
int64 DES::applyREV(int64 x)
{
	// constants
	int64 d32 = int64(4294967295), d64 = int64(18446744073709551615);

	// return result in x
	x = ((x & d32) << 32) | ((x & d64) >> 32);
	return x;
}

// applies T[k[i]] to x
int64 DES::applyT(int64 f, int64 x)
{
	// constants
	int64 d32 = int64(4294967295), d64 = int64(18446744073709551615), y;

	// compute in y the function f
	y = applyF(f, (x & d32));

	// compute in y the xor y
	y ^= (x & d64) >> 32;

	// return result in x
	x = ((x & d32) << 32) | y;
	return x;
}

// calculates function F(v, K[i]) =	P(S(E[v] xor K[i]))
int64 DES::applyF(int64 f, int64 x)
{
	// apply E
	x = applyE(x);

	// apply xor Ki
	x ^= f;

	// apply S
	x = applyS(x);

	// apply P
	x = applyP(1, x);

	// return answer in x
	return x;
}

// calculate E using E expansion permutation
int64 DES::applyE(int64 x)
{
	// ans will be returned
	int64 ans = 0;

	// create 48 bits in ans
	for (int i = 1, j = 48; i <= 48; ++i, --j)
	{
		// ans is added with the bit at the position
		ans |= ((x & int64(pow64[32 - E[i]])) >> (32 - E[i])) << (j - 1);
	}

	// return result in ans
	return ans;
}

// applies P to x
int64 DES::applyP(int f, int64 x)
{
	// ans will be returned
	int64 ans = 0;

	// create 32 bits in ans
	for (int i = 1, j = 32; i <= 32 && f == 1; ++i, --j)
	{
		// ans is added with the bit at the position
		ans |= ((x & int64(pow64[32 - P[i]])) >> (32 - P[i])) << (j - 1);
	}

	// create 32 bits in ans
	for (int i = 1, j = 32; i <= 32 && f == -1; ++i, --j)
	{
		// ans is added with the bit at the position
		ans |= ((x & int64(pow64[32 - Pi[i]])) >> (32 - Pi[i])) << (j - 1);
	}

	// return result in ans
	return ans;
}

// calculates the substitution S
int64 DES::applyS(int64 x)
{
	// ans will be returned
	int64 ans = 0, wi = 0, d6 = 63, d5 = 32, d4 = 16, d3 = 8, d2 = 4, d1 = 2, byte1 = 0, byte2 = 0;

	// calculate Si
	for (int i = 1; i <= 8; ++i)
	{
		// find the leftmost 6 bits in i
		wi = (x & (d6 << ((8 - i) * 6))) >> ((8 - i) * 6);

		// calculate byte1
		byte1 = (((wi & d5) >> 5) << 1) | (wi & 1);

		// calculate byte2
		byte2 = (((wi & d4) >> 4) << 3) | (((wi & d3) >> 3) << 2) | (((wi & d2) >> 2) << 1) | ((wi & d1) >> 1);

		// add to ans
		ans |= M[i][byte1][byte2] << ((8 - i) * 4);
	}

	// return result in ans
	return ans;
}

// applies PC1 to x
int64 DES::applyPC1(int64 x)
{
	// ans will be returned
	int64 ans = 0;

	// create 56 bits in ans
	for (int i = 1, j = 56; i <= 56; ++i, --j)
	{
		// ans is added with the bit at the position
		ans |= ((x & int64(pow64[64 - PC1[i]])) >> (64 - PC1[i])) << (j - 1);
	}

	// return result in ans
	return ans;
}

// applies PC2 to x
int64 DES::applyPC2(int64 x)
{
	// ans will be returned
	int64 ans = 0;

	// create 56 bits in ans
	for (int i = 1, j = 48; i <= 48; ++i, --j)
	{
		// ans is added with the bit at the position
		ans |= ((x & int64(pow64[56 - PC2[i]])) >> (56 - PC2[i])) << (j - 1);
	}

	// return result in ans
	return ans;
}

// applies LSi to key x
int64 DES::applyLS(int i, int64 x)
{
	// necessary variables
	int64 left = 0, right = 0, rmask = int64(268435455), lmask = int64(72057594037927935), bit, d27 = int64(134217728);

	// calculate left part
	left = (x & lmask) >> 28;

	// calculate right part
	right = x & rmask;

	// rotate one bit only in left part
	bit = (left & d27) >> 27;
	left &= (d27 - 1);
	left <<= 1;
	left |= bit;

	// rotate one bit only in right part
	bit = (right & d27) >> 27;
	right &= (d27 - 1);
	right <<= 1;
	right |= bit;

	if (!(i == 1 || i == 2 || i == 9 || i == 16))
	{
		// rotate one bit only in left part
		bit = (left & d27) >> 27;
		left &= (d27 - 1);
		left <<= 1;
		left |= bit;

		// rotate one bit only in right part
		bit = (right & d27) >> 27;
		right &= (d27 - 1);
		right <<= 1;
		right |= bit;
	}

	// answer
	x = (left << 28) | right;
	return x;
}

// crypts x
int64 DES::applyCrypt(int mode, int rds, int64 key, int64 x)
{
	// apply K0
	K[0] = applyPC1(key);

	// generate K1, ..., K16 keys where Ki is 48-bit
	for (int i = 1; i <= rds; ++i)
		K[i] = applyLS(i, K[i-1]);

	// apply PC2 to keys
	for (int i = 1; i <= rds; ++i)
		K[i] = applyPC2(K[i]);
	
	// apply P0 to x
	x = applyP0(1, x);

	// apply TKs to x in encrypt mode
	for (int i = 1; i <= rds && mode == 1; ++i)
	x = applyT(K[i], x);

	// apply TKs to x in decrypt mode
	for (int i = rds; i >= 1 && mode == -1; --i)
	x = applyT(K[i], x);

	// apply REV
	x = applyREV(x);

	// apply P0-1 to x
	x = applyP0(-1, x);

	// return answer
	return x;
}

// crypts crypted text
void DES::crypt(int mode, int rds, ifstream &fe, ifstream &fk, ofstream &fd)
{
	int64 key, x;

	// read key from file
	if (!(fk >> key))
		throw("[DES::crypt] Could not read key");

	// read from prepared file to decrypt each uint64 number
	while (fe >> x)
	{
		x = applyCrypt(mode, rds, key, x);

		// write to encrypted file
		if (!(fd << x))
			throw("[DES::crypt] Could not write decrypted");

		if (!(fd << ' '))
			throw("[DES::crypt] Could not write prepared");
	}
}

// generates all possibilities
void DES::genAll(int f, int w, int64 en[], int64 es[], int64 s[])
{
	// initialize exists
	for (int k = 0; k <= 8; ++k)
		for (int i = 0; i < 64; ++i)
			possible[k][i] = poss[k][i] = 0;

	// for each witness k
	for (int k = 0; k < w; ++k)
	{
		// for each s1, s2, ..., s7, s8 
		for (int i = 1; i <= 8; ++i)
		{
			// si is the string on 4 bytes
			int64 si = (s[k] & (pow64[(9 - i) * 4] - 1)) >> ((8 - i) * 4), res;

			// ei is e' = e xor e* sliced to i
			int64 ei = ((en[k] ^ es[k]) & (pow64[(9 - i) * 6] - 1)) >> ((8 - i) * 6);
			
			// test array
			int test[64];

			// initialize test
			for (int j = 0; j < 64; test[j++] = 0);

			// for any possible B
			for (int64 j = 0; j < pow64[6]; ++j)
			{
				// calculate xor value of S(w) (+) S(w (+) e')
				res = applyS(j << ((8 - i) * 6)) ^ applyS((j ^ ei) << ((8 - i) * 6));

				// take just the necessary bits into account
				res = (res >> ((8 - i) * 4)) & (pow64[4] - 1);

				// INj(Bi, Ci) = {Bj E (2~)~ : Sj(Bj) @ Sj(Bj $ Bi) = Cj)}

				// if this j is correct, put it in test
				if (res == si)
					test[j ^ ((en[k] >> (8 - i) * 6) & (pow64[6] - 1))] = 1;
			}

			// compute possibilities
			for (int j = 0; j < pow64[6]; ++j)
			{
				// if there is in the test
				if (test[j])
					poss[i][j]++;

				// add to correct possible if necessary
				if (poss[i][j] == w)
					possible[i][++possible[i][0]] = j;
			}
		}
	}

	// compute all possibilities
	possib[f][0] = 0;
	for (int i1 = 1; i1 <= possible[1][0]; ++i1)
		for (int i2 = 1; i2 <= possible[2][0]; ++i2)
			for (int i3 = 1; i3 <= possible[3][0]; ++i3)
				for (int i4 = 1; i4 <= possible[4][0]; ++i4)
					for (int i5 = 1; i5 <= possible[5][0]; ++i5)
						for (int i6 = 1; i6 <= possible[6][0]; ++i6)
							for (int i7 = 1; i7 <= possible[7][0]; ++i7)
								for (int i8 = 1; i8 <= possible[8][0]; ++i8)
								{
									// found another possibility
									if (possib[f][0] < 100)
										possib[f][++possib[f][0]] = 0;

									// add a possibility
									possib[f][possib[f][0]] |= (possible[1][i1] << 42);
									possib[f][possib[f][0]] |= (possible[2][i2] << 36);
									possib[f][possib[f][0]] |= (possible[3][i3] << 30);
									possib[f][possib[f][0]] |= (possible[4][i4] << 24);
									possib[f][possib[f][0]] |= (possible[5][i5] << 18);
									possib[f][possib[f][0]] |= (possible[6][i6] << 12);
									possib[f][possib[f][0]] |= (possible[7][i7] << 6);
									possib[f][possib[f][0]] |= (possible[8][i8] << 0);
								}
}

// tests a key solution against plaintext and cryptotext
int DES::testSol(int64 K)
{
	int64 x, y;

	// create a new 3-round DES cryptosystem
	DES *cryptoSystem = new DES();

	ifstream fpl("des3plain.txt");
	ifstream fen("des3encrypted.txt");

	// read all numbers from fpl
	while (fpl >> x)
	{
		// read corresponding number from fen
		if (!(fen >> y))
			throw("[DES::testSol] Could not read number");

		// test if they do not match by using key K
		if ((x = applyCrypt(1, 3, K, x)) != y)
		{
			fpl.close();
			fen.close();

			return 0;
		}
	}

	fpl.close();
	fen.close();

	// they match, hurray :)
	return 1;
}

// decrypts DES in 3 rounds and compares with the plaintext/cryptotext
void DES::dec3(ifstream &fp, ifstream &f3, ofstream &fe)
{
	// u0, v0, u0*, v0* values
	int64 u0n[128], u0s[128], v0n[128], v0s[128], val0n[128], val0s[128];
	int64 u3n[128], u3s[128], v3n[128], v3s[128], val3n[128], val3s[128], e3n[128], e3s[128], s3p[128], K3;
	int64 pk, pf;
	int k3, wits;

	// read the third round values
	for (int i = 0; i < 8; ++i)
	{
		// read the first plaintext
		if (!(f3 >> des3i[i]))
			throw("[DES::dec3] Could not read number");
	}

	for (int i = 0; i < 48; ++i)
	{
		// read the first plaintext
		if (!(f3 >> des3v[i]))
			throw("[DES::dec3] Could not read number");
	}

	// read the number of witnesses
	if (!(fp >> wits))
		throw("[DES::dec3] Could not read number");

	// for each witness
	for (int i = 0; i < wits; ++i)
	{
		// read the first plaintext
		if (!(fp >> val0n[i]))
			throw("[DES::dec3] Could not read number");

		// read the first cryptotext
		if (!(fp >> val3n[i]))
			throw("[DES::dec3] Could not read number");

		// read the second plaintext
		if (!(fp >> val0s[i]))
			throw("[DES::dec3] Could not read number");

		// read the second cryptotext
		if (!(fp >> val3s[i]))
			throw("[DES::dec3] Could not read number");

		// compute u0
		u0n[i] = (val0n[i] & pow64[64]) >> 32;

		// compute v0
		v0n[i] = (val0n[i] & (pow64[32] - 1));

		// compute u0
		u0s[i] = (val0s[i] & pow64[64]) >> 32;

		// compute v0
		v0s[i] = (val0s[i] & (pow64[32] - 1));

		// compute u3
		u3n[i] = (val3n[i] & pow64[64]) >> 32;

		// compute v3
		v3n[i] = (val3n[i] & (pow64[32] - 1));

		// compute u3
		u3s[i] = (val3s[i] & pow64[64]) >> 32;

		// compute v3
		v3s[i] = (val3s[i] & (pow64[32] - 1));

		// calculate e for K3
		e3n[i] = applyE(u3n[i]);
		e3s[i] = applyE(u3s[i]);
		s3p[i] = applyP(-1, (u0n[i] ^ u0s[i]) ^ (v3n[i] ^ v3s[i]));
	}

	// generate all possibilities for k3
	genAll(3, wits, e3n, e3s, s3p);

	// for each possible key k3
	for (k3 = 1; k3 <= possib[3][0]; ++k3)
	{
		// compute K3
		K3 = possib[3][k3];

		// initialize the sure key
		pk = 0;

		// construct the sure key K
		for (int i = 0, j = 48; i < 48; ++i, j--)
		{
			pk |= ((K3 >> (j - 1)) & 1) << (64 - des3v[i]);
		}

		int d = 0;

		// generate all 2^8 = 256 possibilities for the bytes that don't match
		for (int i1 = 0; i1 <= 1; ++i1)
			for (int i2 = 0; i2 <= 1; ++i2)
				for (int i3 = 0; i3 <= 1; ++i3)
					for (int i4 = 0; i4 <= 1; ++i4)
						for (int i5 = 0; i5 <= 1; ++i5)
							for (int i6 = 0; i6 <= 1; ++i6)
								for (int i7 = 0; i7 <= 1; ++i7)
									for (int i8 = 0; i8 <= 1; ++i8)
									{
										pf = pk;

										// add values
										if (i1)
											pf |= pow64[64 - des3i[0]];
										if (i2)
											pf |= pow64[64 - des3i[1]];
										if (i3)
											pf |= pow64[64 - des3i[2]];
										if (i4)
											pf |= pow64[64 - des3i[3]];
										if (i5)
											pf |= pow64[64 - des3i[4]];
										if (i6)
											pf |= pow64[64 - des3i[5]];
										if (i7)
											pf |= pow64[64 - des3i[6]];
										if (i8)
											pf |= pow64[64 - des3i[7]];

										// add parity bits
										for (int v = 64; v; v-= 8)
										{
											// calculate sum of bits
											int sum = 0;
											for (int j = 1; j <= 7; ++j)
												sum += (pf >> (v - j)) & 1;

											// add parity bit so that there is an odd number of bits
											if (!(sum & 1))
												pf |= pow64[v - 8];
										}

										// tests this key against some enc-dec files
										if (testSol(pf))
										{
											// print solution in fe
											if (!(fe << pf))
												throw("[DES::dec3] Could not write K1");

											// print solution in fe
											if (!(fe << '\n'))
												throw("[DES::dec3] Could not write space");
										}
									}
	}
}

// DES destructor
DES::~DES()
{

}
