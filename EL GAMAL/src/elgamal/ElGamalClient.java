/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package elgamal;

import java.math.BigInteger;

/**
 * ElGamalUsage class
 * @author Vlad Manea
 */
public class ElGamalClient extends ElGamal {
    /**
     * ElGamal constructor
     */
    public ElGamalClient(int bits) throws Exception {
        super(bits);
    }

    /**
     * decrypt decrypts encrypted text
     * @param gamma part of encrypted text
     * @param delta part of encrypted text
     * @return decrypted text
     */
    @Override
    public BigInteger decrypt(BigInteger[] encrypted) throws Exception {
        // use decrypt in parent
        return super.decrypt(encrypted);
    }

    /**
     * generateK generates parameter k
     * @return k
     */
    private BigInteger generateK() {
        // generate a random element
        BigInteger k = new BigInteger(p.bitLength(), random);
        // while it is not in [1, p-2]
        while (BigInteger.ONE.compareTo(k) > 0
                || p.subtract(BigInteger.valueOf(2)).compareTo(k) < 0) {
            // generate new random element
            k = new BigInteger(p.bitLength(), random);
        }
        // return it
        return k;
    }

    /**
     * encrypt encrypts plain text
     * @param message message to be encrypted
     * @return encrypted text
     * @throws Exception
     */
    public BigInteger[] encrypt(BigInteger m) throws Exception {
        // test plain text validity
        if (m.compareTo(BigInteger.ZERO) < 0 || m.compareTo(p) >= 0) {
            // throw invalid exception
            throw new Exception("message is not in the interval [0, p-1]");
        }
        // generate k as a random in [1, p-2]
        BigInteger k = generateK();
        // set gamma
        BigInteger gamma = alpha.modPow(k, p);
        // set delta
        BigInteger delta = m.multiply(alphaA.modPow(k, p)).mod(p);
        // return encrypted text
        return new BigInteger[]{gamma, delta};
    }

    /**
     * change changes message to be ok with this ElGamal's p
     * @param m message to be tested
     * @return new message
     */
    public BigInteger change(BigInteger m) {
        // return value
        return m.mod(p);
    }
}
