
import java.util.*;
import java.util.Random;
public class simulateDist {

    public static void main(String[] args) {
        String x = args[2];


        if(x.compareTo("binomial") == 0)
        {
            int n = Integer.parseInt(args[1]);
            //System.out.println(n);
            double p = Double.parseDouble(args[3]);
            //System.out.println(p);
            Binomial(n, p);
        }
        else if(x.compareTo("bernoulli") == 0)
        {
            int n = Integer.parseInt(args[1]);
            double p = Double.parseDouble(args[3]);
            System.out.println(" The randomnly generated sample are: ");
            bernoulli_distribution(p, n);

        }
        else if(x.compareTo("poisson") == 0)
        {
            int n = Integer.parseInt(args[1]);
            double lambda = Double.parseDouble(args[3]);
            poisson(n, lambda);

        }

        else if(x.compareTo("uniform") == 0)
        {
            int n = Integer.parseInt(args[1]);
            float a = Float.parseFloat(args[3]);
            float b = Float.parseFloat(args[4]);
            Uniform(a, b, n);

        }
        else if(x.compareTo("normal") == 0)
        {
            int n = Integer.parseInt(args[1]);
            double mu = Double.parseDouble(args[3]);
            double sigma = Double.parseDouble(args[4]);
            Normal(mu, sigma, n);

        }
        else if(x.compareTo("exponential") == 0)
        {
            int n = Integer.parseInt(args[1]);
            double lambda = Double.parseDouble(args[3]);
            exponential(lambda, n);

        }

        else if(x.compareTo("geometric") == 0)
        {
            int n = Integer.parseInt(args[1]);
            double p = Double.parseDouble(args[3]);
            geometric(p, n);

        }

        else if(x.compareTo("neg-binomial") == 0)
        {
            int n = Integer.parseInt(args[1]);
            int k = Integer.parseInt(args[3]);
            double p = Double.parseDouble(args[4]);
            negativebinomial(p, k, n);

        }
        else if(x.compareTo("gamma") == 0)
        {
            int n = Integer.parseInt(args[1]);
            double alpha = Double.parseDouble(args[3]);
            double lambda = Double.parseDouble(args[4]);
            Gamma(alpha, lambda, n);

        }

        else if(x.compareTo("arb-discrete") == 0)
        {
            int n = Integer.parseInt(args[1]);
            String temp = " ";
            for(int i=3; i<args.length; i++)
            {
                temp = temp +" "+ args[i];
            }

            String[] temp1 = temp.split(" ");
            double[] pr = new double[temp1.length];

            int count =0;
            int iterator =0;


            for(int i=0; i<temp1.length;i++)
            {
                if (temp1[i].isEmpty())
                {
                    pr[i] = -1.0;
                    count++;
                }
                else {
                    pr[i] = Double.parseDouble(temp1[i]);
                }
            }

            double[] p = new double[temp1.length - count];

            for(int j=0; j<pr.length; j++)
            {
                if(pr[j] == -1)
                {

                }
                else
                {
                    p[iterator] = pr[j];
                    iterator++;
                }
            }
            for(int z=0; z<p.length; z++)
            {
                System.out.println(p[z]);
            }

            Arb_disc(p, n);


        }


    }



    // binomial starts here

    public static void Binomial (int trials, double p)
    {

        //Binomial distribution is the probability of exactly x success in n trials
    // The range of each sample element is taken as n
    // elements in the random sample range from 0 to n

    // Function to calculate nCr
        Random rn = new Random();
        int[] succs = new int[trials];
        int min = 0;
        int max = trials;
        System.out.println(" The random generated numbers are: ");
        for (int i = 0; i < (succs.length); i++) {
            succs[i] = rn.nextInt(max - min + 1) + min;
            System.out.print(succs[i]);
            System.out.print(" ");
        }

        System.out.println(" ");
        generateCombo(trials, succs, p);
    }




    /* Binomial */
    public static void generateCombo(int trials, int[] succ, double probsucc){
        int printtrials = trials;//vi p
        int[] succfacts = new int[succ.length];
        int[] diffFacts = new int[succ.length];
        int [] diff = new int[succ.length];
        double[] combos = new double[succ.length];
        int trialfact = 1;
        int succfact = 1;
        int diffFact =1;

        // generates factorial for total number of trials
        for (int i=trials; i>0; i--){
            trialfact = trialfact * i;
        }


        // generate factorial for each possible success
        for (int i=0; i<succ.length; i++){
            if (succ[i] == 0){
                succfacts[i] = 1;
            }
            else{
                succfact =1;
                for(int j=1; j<=succ[i]; j++)
                {
                    succfact = succfact * j;
                    succfacts[i] = succfact;
                }
            }
        }

        for(int i=0; i<diff.length; i++)
        {
            diff[i] = trials - succ[i];
        }
        for (int i=0; i<diff.length; i++){
            if (diff[i] == 0){
                diffFacts[i] = 1;
            }
            else{
                diffFact =1;
                for(int j=1; j<=diff[i]; j++)
                {
                    diffFact = diffFact * j;
                    diffFacts[i] = diffFact;
                }
            }
        }




        // generate combinations for each success and trial number
        for (int i=0; i<succ.length; i++){
            double combo = (trialfact / (succfacts[i] * diffFacts[i])); // added -1
            combos[i] = combo;
        }
        prob(combos, probsucc, trials, succ);

    }

    /* Binomial */
    public static void prob(double[] combos, double probsucc, int trials, int[] succ){
        int printtrials = trials;
        double prob;
        double[] probs = new double[combos.length];
        System.out.println(" ");
        //generates probabilities for each success
        System.out.println(" probability of each are:");
        for (int i=0; i<combos.length; i++){
            prob = (combos[i] * (Math.pow(probsucc, succ[i])) * (Math.pow((1 - probsucc), ( trials - succ[i]))));
            probs[i] = prob;
            System.out.print(probs[i]);
            System.out.print(" ");

        }



    }
// Binomial ends here


    // Bernoulli starts here
    public static void bernoulli_distribution(double p, int n)
    {
        long [] result = new long[n];
        double count =0;
        // Generating random samples of successes and faliures
        for(int i=0; i<n; i++)
        {
            result[i] = Math.round(Math.random());
            if (result[i] == 1)
            {
                count ++;
            }
        }

        

        // distribution list maintains p% 1s and (1-p)% 0s

        if(p>0.1 || p<0.9) {

            if ((count / n) >= (p - 0.1) && (count / n) <= (p + 0.1)) {
                for (int i = 0; i < n; i++) {
                    System.out.print(result[i]);
                    System.out.print(" ");
                }
            } else {
                bernoulli_distribution(p, n);
            }
        }

        else
        {
            if ((count / n) >= (p - 0.01) && (count / n) <= (p + 0.01)) {
                for (int i = 0; i < n; i++) {
                    System.out.print(result[i]);
                    System.out.print(" ");
                }
            } else {
                bernoulli_distribution(p, n);
            }
        }

    }
    // Bernoulli ends here


    /* Poisson */
    public static void poisson(int n, double lambda) {
//        Number of x rare events happen in lambda_ time
//     The range of each sample element is taken as 10
//     elements in the random sample range from 1 to 10, since the range is not specified for each sample
        double[] sample = new double[n];
        Random rn = new Random();
        int max=n;
        int min =0;
        System.out.println("The randomly generated numbers are");
        for(int i=0; i<n; i++)
        {
            sample[i] = rn.nextInt(max - min + 1) + min;
            System.out.print(sample[i]);
            System.out.print(" ");
        }
        System.out.println(" ");
        System.out.println("The probabilities of each are:");

        double[] result = new double[n];
        double exp =0;
        double lamba_power =0;
        double fact=0;
        for(int j=0; j<n; j++)
        {
            exp = Math.exp(-(lambda));
            lamba_power = Math.pow(lambda, sample[j]);
            fact = fact_poi(sample[j]);
            result[j] = ((exp * lamba_power)/ fact);
            System.out.print(result[j]);
            System.out.print(" ");
        }

    }
    public static double fact_poi(double n)
    {
        double res = 1;
        for (int i = 2; i <= n; i++)
            res = res * i;
        return res;
    }

    // Poisson ends here

    // Uniform starts here
    public static void Uniform(float a, float b, int n)
    {
        //Number of x rare events happen in lambda_ time
    // The range of each sample element is taken as 10
    // elements in the random sample range from 1 to 10, since the range is not specified for each sample
        Random rn = new Random();
        double[] sample = new double[n];
        int max = Math.round(b-1);
        int min = Math.round(a+1);
        System.out.println("The randomly generated numbers are:");
        for(int i=0; i<n; i++)
        {
            sample[i] = rn.nextInt(max - min + 1) + min;
            System.out.print(sample[i]);
            System.out.print(" ");
        }

        double[] result = new double[n];
        double f_x =0;

        System.out.println(" ");
        System.out.println(" The probability of each are:");
        for(int j=0; j<n; j++)
        {
            f_x=0;
            f_x = sample[j]/(b-a);
            result[j] = f_x;
            System.out.print(result[j]);
            System.out.print(" ");

        }


    }
    // Uniform ends here

    // exponential starts here
    public static void exponential(double lambda, int n)
    {
        //Exponential dsitributon is F(x) = (1-e^(-lambda x) for x > 0
    // The range of each sample element is taken as n
        // elements in the random sample range from 1 to n, since the range is not specified for each sample
        //Calculates P(X<x)
        Random rn = new Random();
        double[] sample = new double[n];
        int max = n;
        int min = 1;
        System.out.println(" The randomnly generated samples are: ");
        for(int i=0; i<n; i++)
        {
            sample[i] = rn.nextInt(max - min + 1) + min;
            System.out.print(sample[i]);
            System.out.print(" ");
        }
        System.out.println(" ");
        System.out.println(" Their Probabilities are: ");
        double[] result = new double[n];
        double exp =0;
        for(int j=0; j<n; j++)
        {
            exp =0;
            exp = 1- Math.exp((-lambda) * sample[j]);
            result[j] = exp;
            System.out.print(result[j]);
            System.out.print(" ");

        }

    }

    // exponential ends here

    // Geometric starts here
    public static void geometric(double p, int n)
    {

        //Geometric distribution calculates the probability of first success on the xth trial
    // The range of each sample element is taken as 10
    // elements in the random sample range from 1 to 10, since the range is not specified for each sample
        Random rn = new Random();
        int max = n;
        int min =0;
        double[] samples = new double[n];
        System.out.println(" The randomnly generated samples are: ");
        for(int i=0; i<n; i++)
        {
            samples[i] = rn.nextInt(max - min + 1) + min;
            System.out.print(samples[i]);
            System.out.print(" ");
        }

        double[] result = new double[n];
        System.out.println(" ");
        System.out.println(" Their Probabilities");
        for(int i=0; i<n; i++)
        {
            result[i] = Math.pow((1-p), (samples[i] -1)) * p;
        }

        for(int z=0; z<n; z++)
        {
            System.out.print(result[z]);
            System.out.print(" ");
        }

    }

    // Geometric ends here

    // Negativebinomial starts here

    public static void negativebinomial(double p, int k, int n)
    {

        //The xth trial result in the kth success
        //x can vary from k to some number, x>k
    // The range of each sample element is taken as 10
    // elements in the random sample range from k to 10, since the range is not specified for each sample
        int[] sample = new int[n];
        Random rn = new Random();
        int max =n;
        int min=k;
        System.out.println(" The randomnly generated samples are: ");
        for(int i=0; i<n; i++)
        {
            sample[i] = rn.nextInt(max - min + 1) + min;
            System.out.print(sample[i]);
            System.out.print(" ");
        }
        System.out.println(" ");
        double[] result = new double[n];
        int combination =0;
        double failure =0;
        double success =0;
        double z = Double.valueOf(k);
        System.out.println(" Their Probability: ");
        for(int j=0; j<n; j++)
        {
            combination = nCr(sample[j]-1, k-1);
            failure = Math.pow((1-p), (sample[j] -z));
            success = Math.pow(p,z);
            result[j] = (combination * failure * success);
            System.out.print(result[j]);
            System.out.print(" ");
        }

    }

    static int nCr(int n, int r)
    {
        return fact_CR(n) / (fact_CR(r) *
                fact_CR(n - r));
    }

    // Returns factorial of n
    static int fact_CR(int n)
    {
        int res = 1;
        for (int i = 2; i <= n; i++)
            res = res * i;
        return res;
    }

    // Negative Binomial ends here

    // Normal starts here

    public static void Normal(double mu, double sigma, int n)
    {
        //Sample size is given as sample
    // elements in the random sample range from -2 to 2, since the range is not specified for each sample
        double[] samples = new double[n];

        double max =2.0;
        double min = -2.0;
        System.out.println(" The randomnly generated samples are:");
        for(int i=0; i<n; i++)
        {
            samples[i] = (Math.random()*((max-min)+1))+min;
            System.out.print(samples[i]);
            System.out.print(" ");
        }
        System.out.println(" ");
        double[] result = new double[n];
        double z=0;
        System.out.println(" Their Probabilities: ");
        for(int j=0; j<n; j++)
        {
            z = (samples[j] - mu) / sigma;
            result[j] = ((1.0 + erf((z/Math.sqrt(2.0)))) / 2.0);
            System.out.print(result[j]);
            System.out.print(" ");

        }
    }

    public static double erf(double z) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

        // use Horner's method
        double ans = 1 - t * Math.exp( -z*z   -   1.26551223 + t * ( 1.00002368 + t * ( 0.37409196 + t * ( 0.09678418 + t * (-0.18628806 + t * ( 0.27886807 + t * (-1.13520398 + t * ( 1.48851587 + t * (-0.82215223 + t * ( 0.17087277))))))))));
        if (z >= 0) return  ans;
        else        return -ans;
    }

    // Normal ends here

    // Gamma starts here
    public static void Gamma(double alpha, double lambda, int n)
    {
        double[] samples = new double[n];
        Random rn = new Random();
        int max = 10;
        int min = 0;
        System.out.println(" The randomnly generated samples are: ");
        for(int i=0; i<n; i++)
        {
            samples[i] = rn.nextInt(max - min + 1) + min;
            System.out.print(samples[i]);
            System.out.print(" ");
        }
        System.out.println(" ");
        System.out.println(" Their Probabilities: ");
        double[] result = new double[n];
        double temp =0;
        double temp1 =0;
        for(int j=0; j<n; j++)
        {
            temp = ((Math.pow(lambda, alpha)) / gamma(alpha));
            temp1 = Math.pow(samples[j], (alpha-1)) * Math.exp((-lambda) * samples[j]);
            result[j] = temp * temp1;
            System.out.print(result[j]);
            System.out.print(" ");

        }
    }
    public static double logGamma(double x) {
        double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
        double ser = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
                + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
                +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
        return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
    }
    public static double gamma(double x) { return Math.exp(logGamma(x)); }

    // Gamma ends here

    // Arb_disc starts Here
    public static void Arb_disc(double[] p, int n)
    {
        Random rn = new Random();
        double[] samples = new double[n];
        double[] count = new double[n];
        //int max = p.length-1;
        int min =0;
        double elements =0;
        for(int i=0; i<p.length; i++)
        {
            elements = (p[i] * n);
            if(elements<1)
            {
                elements = Math.ceil(elements);
            }
            else {
                elements = Math.floor(elements);
            }

            for(int j=0; j< elements; j++)
            {
                count[j] = min;
            }
            min++;
        }
        System.out.println(" ");

        RandomizeArray(count);

        double[] result = new double[n];
        System.out.println(" The randomnly generated sample is: ");
        for(int z=0; z<result.length; z++)
        {
            getRandom(count);
        }
    }
    public static double[] RandomizeArray(double[] array){
        Random rgen = new Random();  // Random number generator

        for (int i=0; i<array.length; i++) {
            int randomPosition = rgen.nextInt(array.length);
            double temp = array[i];
            array[i] = array[randomPosition];
            array[randomPosition] = temp;
        }

        return array;
    }
    public static void getRandom(double[] array) {
        int rnd = new Random().nextInt(array.length);
        System.out.print(array[rnd]);
        System.out.print(" ");
    }

    // Arb_disc ends here



}
