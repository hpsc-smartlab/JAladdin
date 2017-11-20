package it.uniparthenope;

import java.util.ArrayList;

public class Utility {

    public static ArrayList<Double> linspace(double min, double max, long n){
        ArrayList<Double> v = new ArrayList<Double>();
        double delta = (max - min) / 2;
        double accDelta = 0.0;
        for(long i = 0; i<n; i++){
            v.add(min+accDelta);
            accDelta +=delta;
        }
        return v;
    }

    public static ArrayList<Double> logspace(double min, double max, long n){
        //generate an arraylist of n logaritmically spaced elements between min and max (base 10)
        double logBase = 10;
        double logMin = Math.log10(min);
        double logMax = Math.log10(max);
        double delta = (logMax - logMin) / n;
        double accDelta = 0;
        ArrayList<Double> v = new ArrayList<Double>();
        for(long i = 0; i <= n; ++i){
            v.add(Math.pow(logBase, logMin+accDelta));
            accDelta += delta;
        }
        return v;
    }

    public static Double[][] zeros(int rows, int cols){
        Double[][] matrix = new Double[rows][cols];
        for(int i = 0; i< rows; i++){
            for(int j = 0; j< cols; j++){
                matrix[i][j] = 0.0;
            }
        }
        return matrix;
    }

    public static Double linfit_origin(ArrayList<Double> x, ArrayList<Double> y){//From lingit_origin.m
        // Least square fit of linear funztion through origin.
        // cp. Numerical Recipes,
        // http://hebb.mit.edu/courses/9.29/2002/readings/c15-2.pdf
        // set a=0 in Eq.(15.2.6), with sigma=1 in definitions Eq.(15.2.4)
        // and post-processing after setting a=0
        int S = x.size();
        Double Sx = 0.0;
        ArrayList<Double> Sx2 = new ArrayList<Double>();
        for(Double element : x){
            Sx += element;
            Sx2.add(Math.pow(element,2));
        }
        Double Sy = 0.0;
        for(Double element : y){
            Sy += element;
        }
        Double Sxx = 0.0;
        for(Double element : Sx2){
            Sxx += element;
        }
        Double Delta = S*Sxx - Math.pow(Sx,2);
        return (S*Sxx*Sy/Sx - Sx*Sy)/Delta;
    }

    public static Double fzero(Double k3, Double k2, Double k0, long n_exp, Double v_max_ms, Double x0){
        //This function return zeros of the function defined in ship_resitance.m file
        //k3*x^(3+n_exp)/v_max_ms^n_exp + k2 * x^2 +k0
        //The function domain is defined between -2*x0 and +2*x0 where x0 is the starting point (v_search)

        Double deltaX = 0.001;
        Double minD = -2.0*x0;
        Double maxD = 2.0*x0;
        //Domain definition
        ArrayList<Double> domain = new ArrayList<Double>();
        for(Double xval = minD; xval <= maxD; xval+=deltaX){
            domain.add(xval);
        }
        //Function evaluation
        ArrayList<Double> values = new ArrayList<Double>();
        for(Double element : domain){
            Double value = k3 * Math.pow(element, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(element, 2) + k0;
            values.add(value);
        }
        //function zero
        Double zero = Double.NaN;
        boolean found = false;
        int i = 0;
        while((!found) && (i<values.size()-1)){
            if(values.get(i)==0){
                found = true;
                zero = domain.get(i);
            }
            i++;
        }
        System.out.println("zero="+zero);
        return zero;
    }

}
