# JavaPrimer3

A Java port of Primer3's "thal" function.  Maybe later I will add the "oligotm" function.

To build the JavaPrimer3 jar file run:

    ./gradlew assemble

To run a command-line test do this:

    java -jar ./build/libs/JavaPrimer3.jar oligo1 [oligo2]

The original C version of Primer3 is [here](https://github.com/primer3-org/primer3).

## TODO

Add oligotm() function.

Put thermodynamic parameters in Java property files instead of ThalParameters.java.
