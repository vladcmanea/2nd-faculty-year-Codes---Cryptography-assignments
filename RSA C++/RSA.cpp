#include "RSA.h"

// computes x ^ p mod n
ZZ RSA::ModPower(const ZZ &x, const ZZ &p, const ZZ &n)
{
	// answer
	ZZ a = to_ZZ(1);

	// count power logarithm
	long k = NumBits(p);

	// logarithmic exponentiation
	for (long i = k - 1; i >= 0; i--) 
	{
		// compute squared for ans
		a = SqrMod(a, n);

		// if power has ith bit
		if (bit(p, i) == 1) 
		{
			// multiply a
			a = MulMod(a, x, n);
		}
	}
	
	// return answer
	return a;
}

// RSA constructor testing input files correctness
RSA::RSA(istream &fp, istream &fs, bool copeWiener)
{
	ZZ p, q, n, e, d, f, w;	
	
	// read in public n
	if (!(fp >> n))
		throw("[RSA::RSA] Could not read public n");

	// read in public e
	if (!(fp >> e))
		throw("[RSA::RSA] Could not read public e");

	// read in secret p
	if (!(fs >> p))
		throw("[RSA::RSA] Could not read secret p");

	// read in secret q
	if (!(fs >> q))
		throw("[RSA::RSA] Could not read secret q");

	// read in secret d
	if (!(fs >> d))
		throw("[RSA::RSA] Could not read secret d");

	// test condition for p
	if (NumBits(p) < 512)
		throw("[RSA::RSA] Secret p is too small");

	// test condition for p
	if (ProbPrime(p, 30) == 0)
		throw("[RSA::RSA] Secret p is not prime");

	// test condition for q
	if (NumBits(q) < 512)
		throw("[RSA::RSA] Secret q is too small");

	// test condition for q
	if (ProbPrime(q, 30) == 0)
		throw("[RSA::RSA] Secret q is not prime");

	// test condition for n
	if (n != p * q)
		throw("[RSA::RSA] Public n is not p x q");

	// compute phi(n)
	f = (p - 1) * (q - 1);

	// test condition for e
	if (e <= 1)
		throw("[RSA::RSA] Public e is too small");

	// test condition for e
	if (e >= f)
		throw("[RSA::RSA] Public e is too large");

	// test condition for e
	if (GCD(e, f) != 1)
		throw("[RSA::RSA] Public e is not in Zf*");

	// test condition for e
	if (e >= f)
		throw("[RSA::RSA] Public e is too large");

	// test condition for e
	if (d <= 0)
		throw("[RSA::RSA] Public d is too small");

	// test condition for f
	if (MulMod(e, d, f) != 1)
		throw("[RSA::RSA] Private d is not the inverse of e in Zf");

	// test cope Wiener
	div(w, SqrRoot(SqrRoot(n)), to_ZZ(3));
	if (copeWiener == true && d < w)
		throw("[RSA::RSA] Wiener attack can be performed although copeWiener set to 1");
	else if (copeWiener == false && d >= w)
		throw("[RSA::RSA] Wiener attack cannot be performed although copeWiener set to 0");
}

// RSA constructor generating output files
RSA::RSA(ostream &fp, ostream &fs, bool copeWiener)
{
	// variables
	ZZ p, q, n, e, d, f, g, h, w;

	// initialize random number generator
	SetSeed(to_ZZ((unsigned long)time(0)));

	do 
	{
		// generate random primes
		p = RandomPrime_ZZ(512, 30);
		q = RandomPrime_ZZ(512, 30);
	} 
	while (NumBits(abs(p - q)) < 128 || NumBits(GCD(p - 1, q - 1)) > 128);

	// compute value of n
	n = p * q;

	// compute value of phi(n)
	f = (p - 1) * (q - 1);

	// generate for cope wiener cases
	if (copeWiener == true)
	{
		do 
		{
			do 
			{
				// generate e
				e = RandomBnd(f);
			} 
			// while e not in Zn*
			while (e < 1 || GCD(e, f) != 1);

			// compute d = e ^ (-1) mod f
			d = InvMod(e, f);

			// compute 1 / 3 root4(n)
			div(w, SqrRoot(SqrRoot(n)), to_ZZ(3));
		}
		// while it does not copeWiener
		while (d < w);
	}
	else
	{
		// compute 1 / 3 root4(n)
		div(w, SqrRoot(SqrRoot(n)), to_ZZ(3));

		do 
		{
			// generate d between 1 and w = 1 / 3 root4(n)
			d = RandomBnd(w);
		} 
		// while d not in Zn*
		while (d < 1 || GCD(d, f) != 1);

		// compute e = d ^ (-1) mod f
		e = InvMod(d, f);
	}

	// write in public n
	if (!(fp << n))
		throw("[RSA::RSA] Could not write public n");

	// write in public enter
	if (!(fp << '\n'))
		throw("[RSA::RSA] Could not write public enter");

	// write in public e
	if (!(fp << e))
		throw("[RSA::RSA] Could not write public e");

	// write in secret p
	if (!(fs << p))
		throw("[RSA::RSA] Could not write secret p");

	// write in secret enter
	if (!(fs << '\n'))
		throw("[RSA::RSA] Could not write secret enter");

	// write in secret q
	if (!(fs << q))
		throw("[RSA::RSA] Could not write secret q");

	// write in secret enter
	if (!(fs << '\n'))
		throw("[RSA::RSA] Could not write secret enter");

	// write in secret d
	if (!(fs << d))
		throw("[RSA::RSA] Could not write secret d");
}

// returns x ^ e mod n
void RSA::applyEncrypt(istream &fc, istream &fp, ostream &fa)
{
	// variables
	ZZ n, e, a, x;	

	// read in public n
	if (!(fp >> n))
		throw("[RSA::applyEncrypt] Could not read public n");

	// read in public e
	if (!(fp >> e))
		throw("[RSA::applyEncrypt] Could not read public e");

	// read in plain x
	if (!(fc >> x))
		throw("[RSA::applyEncrypt] Could not read plain x");

	// test x
	if (x <= 0)
		throw("[RSA::applyEncrypt] Plain x is too small");

	// compute x ^ e mod n in answer
	// PowerMod exists and can be used like this too
	a = ModPower(x, e, n);

	// write answer
	if (!(fa << a))
		throw("[RSA::applyEncrypt] Could not write encrypted a");
}

// returns y ^ d mod n
void RSA::applyDecrypt(istream &fc, istream &fp, istream &fs, ostream &fa)
{
	ZZ p, q, n, e, d, f, a, y, d1, d2, a1, a2, x1, x2, v, w, x;	

	// read in public n
	if (!(fp >> n))
		throw("[RSA::applyDecrypt] Could not read public n");

	// read in public e
	if (!(fp >> e))
		throw("[RSA::applyDecrypt] Could not read public e");

	// read in secret p
	if (!(fs >> p))
		throw("[RSA::applyDecrypt] Could not read secret p");

	// read in secret q
	if (!(fs >> q))
		throw("[RSA::applyDecrypt] Could not read secret q");

	// read in secret d
	if (!(fs >> d))
		throw("[RSA::applyDecrypt] Could not read secret d");

	// read in encrypted y
	if (!(fc >> y))
		throw("[RSA::applyDecrypt] Could not read encrypted y");

	// test y
	if (y <= 0)
		throw("[RSA::applyDecrypt] Encrypted y is too small");

	// standard decryption would have been like this:
	// a = ModPower(y, d, n);
	// but now compute y ^ d mod n by using the Chinese Remainder Theorem

	// compute n1 = d mod f(p) = d mod (p - 1)
	rem(d1, d, p - 1);

	// compute n2 = d mod f(p) = d mod (q - 1)
	rem(d2, d, q - 1);

	// compute a1 = y mod p
	rem(a1, y, p);

	// compute a2 = y mod q
	rem(a2, y, q);

	// compute x1 = a1 ^ n1 mod p
	x1 = ModPower(a1, d1, p);

	// compute x2 = a2 ^ n2 mod q
	x2 = ModPower(a2, d2, q);

	// for not to have errors (?!)
	if (p < q)
	{
		// v = p ^ -1 mod q
		v = InvMod(p, q);

		// w = v mod q
		rem(w, v, q);

		// x = diff * w mod q
		rem(x, (x2 - x1) * w, q);

		// compute a
		a = x1 + p * x;
	}
	else
	{
		// v = q ^ -1 mod p
		v = InvMod(q, p);

		// w = v mod p
		rem(w, v, p);

		// x = diff * w mod p
		rem(x, (x1 - x2) * w, p);

		// compute a
		a = x2 + p * x;
	}

	// write answer
	if (!(fa << a))
		throw("[RSA::applyDecrypt] Could not write plain a");
}

// finds if l, d can be the l and d in e * d - 1 = l * phi(n)
bool RSA::criteria(ZZ l, ZZ d, ZZ n, ZZ e)
{
	ZZ a, b, r, t, r1, r2, q1, q2;

	// a is n
	a = n;

	// compute remainder of (e * d - 1) / l
	if (l != 0)
		rem(r, e * d - 1, l);
	else
		return false;

	// test if r is zero
	if (r != 0)
		return false;

	//cout << "rest = " << r << '\n';

	// compute quotient of (e * d - 1) / l
	if (l != 0)
		div(b, e * d - 1, l);
	else
		return false;

	//cout << "ed-1/l = " << b << '\n';


	// compute delta = b ^ 2
	if (b >= a + 1)
		t = (b - a - 1) * (b - a - 1);
	else
		t = (a + 1 - b) * (a + 1 - b);

	//cout << "b2 = " << t << '\n';

	// compute delta = delta - 4 * a * c
	if (t >= 4 * a)
		t = t - 4 * a;
	else
		return false;

	//cout << "b2-4ac = " << t << '\n';

	// test if delta is something squared
	if (SqrRoot(t) * SqrRoot(t) != t)
		return false;

	//cout << "rad(delta) = " << SqrRoot(t) << '\n';

	// compute quotient
	rem(r1, a + 1 - b + SqrRoot(t), to_ZZ(2));
	rem(r2, a + 1 - b - SqrRoot(t), to_ZZ(2));

	// compute remainder
	div(q1, a + 1 - b + SqrRoot(t), to_ZZ(2));
	div(q2, a + 1 - b - SqrRoot(t), to_ZZ(2));

	//cout << "rest(-b+rad(delta)/2a) = " << r1 << "; cat(-b+rad(delta)/2a) = " << q1 << '\n';
	//cout << "rest(-b-rad(delta)/2a) = " << r2 << "; cat(-b-rad(delta)/2a) = " << q2 << '\n';

	// test if at least one solution exists
	if (r1 == 0 && q1 >= 0 || r2 == 0 && q2 >= 0)
		return true;

	//cout << "rest(-b+rad(delta)/2a) = " << r << '\n';

	return true;
}

// applies Wiener decryption
void RSA::applyWienerDecrypt(istream &fc, istream &fp, ostream &fa)
{
	ZZ n, e, d, f, y, i, l, a[3], b[3], q[3], r[3];

	// read in public n
	if (!(fp >> n))
		throw("[RSA::applyWienerDecrypt] Could not read public n");

	// read in public e
	if (!(fp >> e))
		throw("[RSA::applyWienerDecrypt] Could not read public e");

	// read in encrypted y
	if (!(fc >> y))
		throw("[RSA::applyWienerDecrypt] Could not read encrypted y");

	// test y
	if (y <= 0)
		throw("[RSA::applyWienerDecrypt] Encrypted y is too small");

	//n = 8999999;
	//e = 3174353;

	//cout << "e = " << e << '\n';

	// set i to 0
	i = 0;

	do
	{
		// increment i
		i = i + 1;

		// update all
		a[0] = a[1], a[1] = a[2], b[0] = b[1], b[1] = b[2], q[0] = q[1], q[1] = q[2], r[0] = r[1], r[1] = r[2];

		// compute q i, r i, alpha i, beta i for different cases
		if (i == 1)
		{
			div(q[2], e, n);
			rem(r[2], e, n);
			a[2] = q[2];
			b[2] = to_ZZ(1);
		}
		else if (i == 2)
		{
			if (r[1] > 0)
			{
				div(q[2], n, r[1]);
				rem(r[2], n, r[1]);
			}
			else
				break;

			a[2] = q[0] * q[1] + 1;
			b[2] = q[2];
		}
		else
		{
			if (r[1] > 0)
			{
				div(q[2], r[0], r[1]);
				rem(r[2], r[0], r[1]);
			}
			else
				break;

			a[2] = q[2] * a[1] + a[0];
			b[2] = q[2] * b[1] + b[0];
		}

		//cout << "i = " << i << '\n'; 
	}
	while(criteria(a[2], b[2], n, e) == false);

	//cout << "d = " << b[2] << '\n';

	// compute answer
	f = ModPower(y, b[2], n);

	// write answer
	if (!(fa << f))
		throw("[RSA::applyWienerDecrypt] Could not write plain a");

	// write answer
	if (!(fa << '\n'))
		throw("[RSA::applyWienerDecrypt] Could not write enter");

	// write answer
	if (!(fa << b[2]))
		throw("[RSA::applyWienerDecrypt] Could not write secret d");
}

// divides into factors
void RSA::applyFactors(istream &fp, istream &fs, ostream &fa)
{
	ZZ p, q, n, e, d, f, r, v, x, u, w;
	long m, s;

	// read in public n
	if (!(fp >> n))
		throw("[RSA::applyFactors] Could not read public n");

	// read in public e
	if (!(fp >> e))
		throw("[RSA::applyFactors] Could not read public e");

	// read in secret p but skip it
	if (!(fs >> f))
		throw("[RSA::applyFactors] Could not read secret p");

	// read in secret q but skip it
	if (!(fs >> f))
		throw("[RSA::applyFactors] Could not read secret q");

	// read in secret d
	if (!(fs >> d))
		throw("[RSA::applyFactors] Could not read secret d");

	// compute e * d - 1 and redundance
	r = e * d - 1;

	// compute mod 2 remainder
	rem(v, r, to_ZZ(2));

	// while v is even
	for (s = 0; v == 0; ++s)
	{
		// divide v
		div(r, r, to_ZZ(2));

		// compute mod 2 remainder
		rem(v, r, to_ZZ(2));
	}

	// redundancy
	m = s;

	// generate random values for x
	while (1)
	{
		// reload s
		s = m;

		// generate random x from Zn*
		x = RandomBnd(n);
		while (GCD(x, n) != 1)
			x = RandomBnd(n);

		// compute 2 ^ s in w
		power(w, 2, s);

		// iterate through s
		for (--s; s >= 0; --s)
		{
			// divide w
			div(w, w, to_ZZ(2));
			
			// compute x ^ (2 ^ s * r) mod n in u
			u = ModPower(x, w * r, n);

			//non-trivial
			if (u != 1 && u != n - 1)
			{
				// compute p
				p = GCD(u + 1, n);

				// compute q
				div(q, n, GCD(u + 1, n));

				// write answer
				if (!(fa << p))
					throw("[RSA::applyFactors] Could not write factor p");

				// write answer
				if (!(fa << '\n'))
					throw("[RSA::applyFactors] Could not write enter");

				// write answer
				if (!(fa << q))
					throw("[RSA::applyFactors] Could not write factor q");

				// finished
				return;
			}
			
			// otherwise
			if (u == n - 1)
				break;
		}
	}
}

// RSA destructor
RSA::~RSA()
{

}
