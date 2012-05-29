package elgamal;

import java.math.BigInteger;
import java.util.Random;

/**
 * Class ElGamal
 * @author Vlad Manea
 */
public class ElGamal {

    // ElGamal parameters
    protected BigInteger p, alpha, alphaA;
    private BigInteger a;
    private int bits;
    protected Random random = new Random();

    /**
     * generateQ generates prime q
     * @return q
     */
    private BigInteger generateQ() {
        // generate q as a prime
        BigInteger q = BigInteger.probablePrime(bits - 1, random);
        // while 2*q+1 is not a prime in itself
        while (!q.multiply(BigInteger.valueOf(2)).
                add(BigInteger.valueOf(1)).isProbablePrime(10)) {
            // generate q again
            q = BigInteger.probablePrime(bits - 1, random);
        }
        // return a prime q
        return q;
    }

    private void generateAlpha() {
        // generate new random element
        alpha = new BigInteger(p.bitLength(), random);
        // while alpha is not in Zp*
        while (p.compareTo(alpha) >= 0
                || p.gcd(alpha).compareTo(BigInteger.valueOf(1)) > 0
                || alpha.modPow(BigInteger.valueOf(2), p).
                compareTo(BigInteger.valueOf(1)) == 0
                || alpha.modPow(p.subtract(BigInteger.valueOf(1).
                divide(BigInteger.valueOf(2))).
                mod(p.subtract(BigInteger.valueOf(1))), p).
                compareTo(BigInteger.valueOf(1)) == 0) {
            // generate new random element
            alpha = new BigInteger(p.bitLength(), random);
        }
    }

    /**
     * generatePAlpha generates parameters p and alpha
     */
    private void generatePAlpha() {
        // generate a random prime q = ans[1] of 1023 bits
        BigInteger q = generateQ();
        // now you know p = 2*q+1 = 2*ans[1]+1 where p-1 = 2^1*q^1
        p = q.multiply(BigInteger.valueOf(2)).add(BigInteger.valueOf(1));
        // choose a random element in Zp*
        generateAlpha();
    }

    /**
     * generateA generates parameter a
     */
    private void generateA() {
        // generate a random element
        a = new BigInteger(p.bitLength(), random);
        // while it is not in [1, p-2]
        while (BigInteger.ONE.compareTo(a) > 0
                || p.subtract(BigInteger.valueOf(2)).compareTo(a) < 0) {
            // generate new random element
            a = new BigInteger(p.bitLength(), random);
        }
    }

    /**
     * generateAlphaA generates alpha^a mod p
     */
    private void generateAlphaA() {
        alphaA = alpha.modPow(a.mod(p.subtract(BigInteger.valueOf(1))), p);
    }

    /**
     * decrypt decrypts encrypted text
     * @param gamma part of encrypted text
     * @param delta part of encrypted text
     * @return decrypted text
     */
    protected BigInteger decrypt(BigInteger[] encrypted) throws Exception {
        // test encrypted text validity
        if (encrypted.length != 2) {
            // throw invalid exception
            throw new Exception("encrypted text is invalid");
        }
        // return decrypted text
        return encrypted[1].multiply(encrypted[0].modPow(
                p.subtract(BigInteger.ONE).subtract(a), p)).mod(p);
    }

    /**
     * ElGamal constructor
     */
    protected ElGamal(int bits) throws Exception {
        // are bits enough?
        if (bits < 128) {
            // throw invalid exception
            throw new Exception("Number of bits is too small, try "
                    + "[128, infinity]");
        }
        // set number of bits
        this.bits = bits;
        // generate p and alpha
        generatePAlpha();
        // generate a
        generateA();
        // compute alpha^a
        generateAlphaA();
    }
}
