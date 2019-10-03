package edu.scripps.yates;

import java.io.BufferedReader;
import java.io.FileReader;

public class PwMsScanFixer {

    public static void main(String[] args) throws Exception {
        // write your code here

        BufferedReader br = new BufferedReader(new FileReader(args[0]));
        String eachline = null;

        int scanCount=1;
        while ((eachline = br.readLine()) != null) {

            if (eachline.startsWith("S\t")) {
               // System.out.println(eachline);
                String[] arr = eachline.split("\t");
                System.out.print("S\t" + scanCount + "\t" + scanCount + "\t");

  /*
                for(int i=1;i<arr.length;i++) {
                    System.out.print("\t" + arr[i]);
                }
*/
                System.out.println(arr[arr.length-1]);

                scanCount++;
            } else {
                System.out.println(eachline);
            }



        }
    }
}
