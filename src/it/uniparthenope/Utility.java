package it.uniparthenope;


import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import it.uniparthenope.Boxing.meshgridResults;

public class Utility {

    public static ArrayList<Double> linspace(double min, double max, long n){
        ArrayList<Double> v = new ArrayList<Double>();
        double delta = (max - min) / (n-1);
        double accDelta = 0.0;
        for(long i = 0; i<n; i++){
            v.add(min+accDelta);
            accDelta +=delta;
        }
        return v;
    }

    public static ArrayList<Double> logspace(double min, double max, long n){
        //generate an arraylist of n logaritmically spaced elements between min and max (base 10)
        ArrayList<Double> log = new ArrayList<>();
        ArrayList<Double> exps = linspace(min, max, n);
        Double base = 10.0;
        for(Double exp: exps){
            log.add(Math.pow(base,exp));
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

    public static Double[][] ones(int size){
        return ones(size, size);
    }

    public static Double[][] ones(int rows, int cols){
        Double[][] matrix = new Double[rows][cols];
        for(int i = 0; i< rows; i++){
            for(int j = 0; j< cols; j++){
                matrix[i][j] = 1.0;
            }
        }
        return matrix;
    }

    public static Double[][] NaNmatrix(int rows, int cols){
        Double[][] matrix = new Double[rows][cols];
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                matrix[i][j]=Double.NaN;
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


    //This function return zero of the function defined in ship_resitance.m file
    //k3*x^(3+n_exp)/v_max_ms^n_exp + k2 * x^2 +k0
    //it's an approssimation of fzero matlab function.
    public static Double Newton(Double k3, Double k2, Double k0, long n_exp, Double v_max_ms, Double x0){
        Double xk = x0;
        int k=0;
        Double delta_ass = 0.0000001;//Max absolute error threshold
        int kmax = 1000000; //Max iterations number
        Double fxk = k3 * Math.pow(xk, 3+n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk, 2) + k0;
        Double fprimoxk = k3*(3+n_exp)*Math.pow(xk,(3+n_exp)-1)/Math.pow(v_max_ms, n_exp) + k2*2*xk;
        Double ck= -fxk/fprimoxk;
        while( (Math.abs(ck) > delta_ass*Math.ulp(1.0)) && (k<kmax) ){
            xk+=ck;
            fxk = k3 * Math.pow(xk, 3+n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk, 2) + k0;
            fprimoxk = k3*(3+n_exp)*Math.pow(xk,(3+n_exp)-1)/Math.pow(v_max_ms, n_exp) + k2*2*xk;
            ck=-fxk/fprimoxk;
            k++;
        }
        return xk;
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

    public static Double fzero_secant(Double k3, Double k2, Double k0, long n_exp, Double v_max_ms, Double x0){
        //This function return zeros of the function defined in ship_resitance.m file
        //k3*x^(3+n_exp)/v_max_ms^n_exp + k2 * x^2 +k0
        //This function uses secant method to find function root.
        return secant(x0, (x0+ 0.5), k3, k2, k0, n_exp, v_max_ms);
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
        Double fxk = k3 * Math.pow(xk, 3+n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk, 2) + k0;
        Double fxk1 = k3 * Math.pow(xk1, 3+n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk1, 2) + k0;
        //line gradient approssimation
        Double pxk = (fxk-fxk1)/(xk-xk1);
        //k-step correction
        Double ck = -(fxk/pxk);
        while( (Math.abs(ck) > delta_ass*Math.ulp(1.0)) && (k<=kmax) ){//stop criteria
            xk1 = xk;
            fxk1 = fxk;
            xk = xk + ck;
            fxk = k3 * Math.pow(xk, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk, 2) + k0;
            pxk = (fxk-fxk1)/(xk-xk1);
            ck = -(fxk/pxk);
            k++;
        }
        root=xk;
        return root;
    }

    public static long getUnsignedInt(int x) {
        return x & 0x00000000ffffffffL;
    }

    public static int BigToLittleEndian(int x){
       return ByteBuffer.allocate(4).order(ByteOrder.BIG_ENDIAN).putInt(x).order(ByteOrder.LITTLE_ENDIAN).getInt(0);
    }


    public static meshgridResults meshgrid(ArrayList<Double> x){
        return meshgrid(x,x);
    }

    public static meshgridResults meshgrid(ArrayList<Double> x, ArrayList<Double> y){
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
                Y[i][j] = y.get(i);
            }
        }
        return new meshgridResults(X,Y);
    }


    public static Double[] reshape(Double[][] A, int dim){// Reshape 1
        if(A.length*A[0].length != dim){
            System.out.println("size(A) must be = dim!");
            System.exit(0);
        }
        Double[] output = new Double[dim];
        int i=0;
        int j=0;
        for(int k=0;k<dim;k++){
            output[k]=A[i][j];
            i++;
            if(i==A.length){
                i=0;
                j++;
            }
        }

        return output;
    }

    public static Double[][] reshape(ArrayList<Double> A, int[] dim){//reshape 2
        int nRows = dim[0];
        int nCols = dim[1];
        if(A.size() != nRows*nCols){
            System.out.println("size(A) must be = to nRows*nCols!");
            System.exit(0);
        }
        Double[][] out = new Double[nRows][nCols];
        int i=0;
        int j=0;
        for(int k=0;k<A.size();k++){
            out[i][j]=A.get(k);
            i++;
            if(i==nRows){
                i=0;
                j++;
            }

        }
        return out;
    }

    public static Double[][] reshape(Double[] A, int[] dim){//reshape 3
        int nRows = dim[0];
        int nCols = dim[1];
        if(A.length != nRows*nCols){
            System.out.println("A.length must be = to nRows*nCols!");
            System.exit(0);
        }
        Double[][] out = new Double[nRows][nCols];
        int i=0;
        int j=0;
        int n_element = 0;
        for(int k=0;k<A.length;k++){
            out[i][j]=A[n_element];
            n_element++;
            i++;
            if(i==nRows){
                i=0;
                j++;
            }

        }
        return out;
    }

    public static Double[][] reshape(Double[][] A, int[] dim){//Reshape 4
        int nRows = dim[0];
        int nCols = dim[1];
        if(A.length*A[0].length != nRows*nCols){
            System.out.println("size(A) must be = to nRows*nCols!");
            System.exit(0);
        }
        Double[][] out = new Double[nRows][nCols];
        int i=0;
        int j=0;
        int k=0;
        int l=0;
        for(int u=0;u<nRows*nCols;u++){
            out[k][l]=A[i][j];
            i++;
            k++;
            if(i==A.length){
                i=0;
                j++;
            }
            if(k==nRows){
                k=0;
                l++;
            }
        }
        return out;
    }

    public static boolean any(long[] A, String condition, long threshold){
        boolean output = false;
        switch(condition){
            case(">"):
                for(int i=0;i<A.length;i++){
                    if(A[i]>threshold){
                        output = true;
                    }
                }
                break;
            case("<"):
                for(int i=0;i<A.length;i++){
                    if(A[i]<threshold){
                        output = true;
                    }
                }
                break;
            case("=="):
                for(int i=0;i<A.length;i++){
                    if(A[i]==threshold){
                        output = true;
                    }
                }
                break;
            case(">="):
                for(int i=0;i<A.length;i++){
                    if(A[i]>=threshold){
                        output = true;
                    }
                }
                break;
            case("<="):
                for(int i=0;i<A.length;i++){
                    if(A[i]<=threshold){
                        output = true;
                    }
                }
                break;
        }
        return output;
    }


    public static Double[][] transposeMatrix(Double[][] matrix){
        Double[][] output = new Double[matrix[0].length][matrix.length];
        for(int i=0;i< matrix.length;i++){
            for(int j=0;j<matrix[0].length;j++){
                output[j][i]=matrix[i][j];
            }
        }
        return output;
    }

    public static Double[][] MatrixComponentXcomponent(Double[][] a, Double[][] b){
        int nRows = a.length;
        int nCols = a[0].length;
        if(a.length!=b.length){
            System.out.println("a.length != b.length");
            System.exit(0);
        }
        if(a[0].length!=b[0].length){
            System.out.println("a[0].length != b[0].length");
            System.exit(0);
        }
        Double[][] output = new Double[nRows][nCols];
        for(int i=0;i<nRows;i++){
            for(int j=0;j<nCols;j++){
                output[i][j]=a[i][j]*b[i][j];
            }
        }
        return  output;
    }

    public static Double min(Double[] array){
        Double min = array[0];
        for(int i=1;i<array.length;i++){
            if(array[i]<min){
                min = array[i];
            }
        }
        return min;
    }

    public static Double[] min(Double[][] matrix, int flag){
        //MATLAB min(x,[],y) implementation
        if(flag==1){
            //min(matrix,[],1) (min cols)
            Double[] out = new Double[matrix[0].length];
            for(int i=0;i<matrix[0].length;i++){
                Double minTmp = matrix[0][i];
                for(int j=1;j<matrix.length;j++){
                    if(matrix[j][i]<minTmp){
                        minTmp = matrix[j][i];
                    }
                }
                out[i]=minTmp;
            }
            return  out;
        } else {
            //min(matrix,[],2) (min rows)
            Double[] out = new Double[matrix.length];
            for(int i=0;i<matrix.length;i++){
                Double minTmp = matrix[i][0];
                for(int j=1;j<matrix[0].length;j++){
                    if(matrix[i][j]<minTmp){
                        minTmp=matrix[i][j];
                    }
                }
                out[i]=minTmp;
            }
            return out;
        }
    }

    public static int sign(long x){
        if(x>0){
            return 1;
        } else if(x==0){
            return 0;
        } else {
            return -1;
        }
    }

    public static int sign(double x){
        if(x>0){
            return 1;
        } else if(x==0){
            return 0;
        } else {
            return -1;
        }
    }

}
