/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package elgamal;

import java.math.BigInteger;
import java.util.Random;

/**
 *
 * @author vladm
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
        ElGamalClient g = new ElGamalClient(128);
        BigInteger x = g.change(new BigInteger(128, new Random()));
        BigInteger[] y = g.encrypt(x);
        System.out.println(x.toString());
        System.out.println(g.decrypt(y).toString());
    }
}
