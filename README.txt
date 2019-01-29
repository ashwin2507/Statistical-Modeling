Statistical-Modeling

Name: Ashwin Babu
ID: 1001392860
Class: CSE 5301

Programming Language Used: Java 8, Inc.
--------------------------

----------------
CODE STRUCTURE
-----------------
1. All the distributions are implemented as different functions.
2. Make sure to name the correct distributions before running.
3. Bernoulli and Arbitrary discrete distributions only prints the final output of selecting the elements with the assigned probability. 


--------------------
HOW TO RUN THE CODE
--------------------

Command
-------
javac simulateDist.java to execute java file
java simulateDist Dist <number-of-samples> <distribution> <parameters>

For example:
java simulateDist Dist 10 binomial 10 0.6

NOTE
----
The program generates simulated samples from a given distribution by using samples generated from a Standard Uniform Distribution (Random Number Generator). 
The type of distribution and the parameters are given as command line arguments. Command line format is simulateDist <number-of-samples> <distribution> <parameters>
The simulateDist should be in the same directory as at the terminal.

Distribution Implemented are as follows:
Binomial(n, p);
bernoulli_distribution(p, n);
poisson(n, lambda);
Uniform(a, b, n);
Normal(mu, sigma, n);
exponential(lambda, n);
geometric(p, n);
negativebinomial(p, k, n);
Gamma(alpha, lambda, n);
Arb_disc(p, n);

