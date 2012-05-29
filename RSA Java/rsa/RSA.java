/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package rsa;

import java.math.BigInteger;
import java.util.Random;

/**
 *
 * @author vladm
 */
public class RSA {

    protected BigInteger p, q, n, f, e, d;

    public RSA() {
        Random generator = new Random();
        p = BigInteger.probablePrime(1024, generator);
        q = BigInteger.probablePrime(1024, generator);
        while (q.compareTo(p) == 0) {
            q = BigInteger.probablePrime(1024, generator);
        }
        n = p.multiply(q);
        f = p.subtract(BigInteger.ONE).multiply(q.subtract(BigInteger.ONE));
        e = BigInteger.valueOf(3 + generator.nextInt(Integer.MAX_VALUE));
        while (e.compareTo(f) >= 0 || !e.gcd(f).equals(BigInteger.ONE)) {
            e = BigInteger.valueOf(3 + generator.nextInt(Integer.MAX_VALUE));
        }
        d = e.modInverse(f);
        while (d.compareTo(BigInteger.ZERO) < 0) {
            d = d.add(f);
        }
    }

    public BigInteger encrypt(BigInteger m) {
        if (m.compareTo(BigInteger.ZERO) > 0 && m.compareTo(n) < 0) {
            return m.modPow(e, n);
        } else {
            return BigInteger.ZERO;
        }
    }

    public BigInteger decrypt(BigInteger m) {
        if (m.compareTo(BigInteger.ZERO) > 0 && m.compareTo(n) < 0) {
            return m.modPow(d, n);
        } else {
            return BigInteger.ZERO;
        }
    }

    public BigInteger sign(BigInteger m) {
        return decrypt(m);
    }

    public BigInteger verify(BigInteger m) {
        return encrypt(m);
    }
}
