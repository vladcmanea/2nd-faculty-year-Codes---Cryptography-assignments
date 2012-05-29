/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package rsa;

import java.math.BigInteger;

/**
 *
 * @author vladm
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        RSA rsa = new RSA();
        BigInteger m = BigInteger.valueOf(1234567890);
        BigInteger em = rsa.encrypt(m);
        BigInteger dem = rsa.decrypt(em);
        System.out.println(m.equals(dem));        
    }
}
