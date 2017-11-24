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
        Double step = (Math.abs(max)-Math.abs(min))/(n-1);
        Double base = 10.0;
        Double accDelta = step;
        ArrayList<Double> log = new ArrayList<>();
        for(int i =0; i<n; i++){
            Double val = Math.pow(base,accDelta);
            log.add(val);
            accDelta+=step;
        }
        return log;
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
        //The function domain is defined between -10*x0 and +10*x0 where x0 is the starting point (v_search)
        //This function uses bisection method to find function root. It was used because if initial
        //condition are satisfied, always converge (global method)

        //Initial conditions: f(x) is a continous function in [a,b]
        //sign(f(a)) != sign(f(b))
        Double a = -10*x0;
        Double b= 10*x0;
        //Checking initial conditions:
        //Function evaluations in a and b
        Double fa = k3 * Math.pow(a, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(a, 2) + k0;//Function evaluation in 'a'
        Double fb = k3 * Math.pow(b, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(b, 2) + k0;//Function evaluation in 'b'
        Double zero = Double.NaN;
        if((fa*fb)<0){//if sign(fa) != sign(fb)
            System.out.println("OK! Initial conditions satisfied! fa="+fa+", fb= "+fb);
            zero = bisection(a,b,k3,k2,k0,n_exp,v_max_ms);
        } else{
            System.out.println("Initial conditions not satisfied! fa="+fa+", fb= "+fb);
            System.exit(-1);
        }

        return zero;
    }

    public static Double fzero_secant(Double k3, Double k2, Double k0, long n_exp, Double v_max_ms, Double x0){
        //This function return zeros of the function defined in ship_resitance.m file
        //k3*x^(3+n_exp)/v_max_ms^n_exp + k2 * x^2 +k0
        //This function uses secant method to find function root.
        return secant(x0, (x0+ 0.5), k3, k2, k0, n_exp, v_max_ms);
    }

    private static Double bisection(Double a, Double b, Double k3, Double k2, Double k0, long n_exp, Double v_max_ms){
        //Solve f(x) = 0 with recoursive bisection method.
        Double delta_ass = 0.001;//Max absolute error threshold
        if(Math.abs((b-a)) <= (delta_ass + Math.max(Math.abs(a),Math.abs(b)))){//base case
            return (a+b)/2;//Function root is approssimated with the middle point of the range
        }
        else{//recoursive function calls
            Double fa = k3 * Math.pow(a, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(a, 2) + k0;//Function evaluation in 'a'
            Double meanPoint = (a+b)/2;//Mean point of the range [a, b]
            Double fmean = k3 * Math.pow(meanPoint, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(meanPoint, 2) + k0; // Function evaluation in 'meanPoint'
            if((fa*fmean)<0){//if sign(fa) != sign(fmean)
                return bisection(a, meanPoint, k3, k2, k0, n_exp, v_max_ms);//Root is in the range [a, meanPoint]
            }else {
                return bisection(meanPoint, b, k3, k2, k0, n_exp, v_max_ms);//Root is in the range [meanPoint, b]
            }
        }
    }

    private static Double secant(Double x0, Double x1, Double k3, Double k2, Double k0, long n_exp, Double v_max_ms){
        //Solve f(x) = 0 with secant method.
        Double delta_ass = 0.0000001;//Max absolute error threshold
        int kmax = 10000000; //Max iterations number
        int k = 1;
        //Initial xk and xk+1
        Double xk = x0;
        Double xk1 = x1;
        Double root = xk;
        //Function evaluation in xk and xk+1
        Double fxk = k3 * Math.pow(xk, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk, 2) + k0;
        Double fxk1 = k3 * Math.pow(xk1, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk1, 2) + k0;
        //line gradient approssimation
        Double pxk = (fxk-fxk1)/(xk-xk1);
        //k-step correction
        Double ck = -(fxk/pxk);
        while( (ck > delta_ass*Math.ulp(1.0)) && (k<=kmax) ){//stop criteria
            xk1 = xk;
            fxk1 = fxk;
            xk = xk + ck;
            fxk = k3 * Math.pow(xk, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk, 2) + k0;
            pxk = (fxk-fxk1)/(xk-xk1);
            ck = -(fxk/pxk);
            k++;
        }
        return root;
    }

    public static ArrayList<Object> meshgrid(ArrayList<Double> x){
        return meshgrid(x,x);
    }

    public static ArrayList<Object> meshgrid(ArrayList<Double> x, ArrayList<Double> y){
        //[X, Y] =meshgrid(x,y) returns 2-D grid coordinates based on the coordinates contained
        // in vectors x and y. X is a matrix where each row is a copy of x,
        // and Y is a matrix where each column is a copy of y.
        // The grid represented by the coordinates X and Y has length(y) rows and length(x) columns.
        int nRows = y.size();
        int nCols = x.size();
        Double[][] X = new Double[nRows][nCols];
        Double[][] Y = new Double[nRows][nCols];
        for(int i =0;i<nRows;i++){
            for(int j=0;j<nCols; j++){
                X[i][j] = x.get(j);
                Y[i][j] = y.get(j);
            }
        }

        ArrayList<Object> output = new ArrayList<>();
        output.add((Object) X);
        output.add((Object) Y);
        return output;
    }
}
