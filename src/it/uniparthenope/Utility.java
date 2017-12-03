package it.uniparthenope;


import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Collections;

import it.uniparthenope.Boxing.meshgridResults;
import it.uniparthenope.Debug.MyFileWriter;
import org.apache.commons.math3.complex.Complex;

public class Utility {

    public static long Tic(){
        return System.nanoTime();
    }

    public static long Toc(long startTime){
        return  System.nanoTime() - startTime;
    }

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


    public static double[][] zeros(int rows, int cols){
        double[][] matrix = new double[rows][cols];
        for(int i = 0; i< rows; i++){
            for(int j = 0; j< cols; j++){
                matrix[i][j] = 0.0;
            }
        }
        return matrix;
    }



    public static double[][] ones(int size){
        return ones(size,size);
    }

    public static double[][] ones(int rows, int cols){
        double[][] matrix = new double[rows][cols];
        for(int i = 0; i< rows; i++){
            for(int j = 0; j< cols; j++){
                matrix[i][j] = 1.0;
            }
        }
        return matrix;
    }


    public static double[][] NaNmatrix(int rows, int cols){
        double[][] matrix = new double[rows][cols];
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                matrix[i][j]=Double.NaN;
            }
        }
        return matrix;
    }

    public static double[] NaNarray(int dim){
        double[] array = new double[dim];
        for(int i=0;i<dim;i++){
            array[i]=Double.NaN;
        }
        return array;
    }

    public static double nr_cubic(double a, double b, double c){
        //analytical search of roots: return max root only
        //http://en.wikipedia.org/wiki/Cubic_function#General_formula_for_roots
        double QQ = (Math.pow(a,2) -3*b)/9;
        double RR = (2*Math.pow(a,3) -9*a*b + 27*c)/54;
        double DD = Math.pow(RR,2) - Math.pow(QQ,3);

        if(DD<0){//3 real solution (large wave height) ;
            double theta = Math.acos( RR/Math.sqrt(Math.pow(QQ,3)) );

            double x1 = -2*Math.sqrt(QQ)* Math.cos(theta/3) - a/3;
            double x2 = -2*Math.sqrt(QQ)* Math.cos((theta+2*Math.PI)/3) - a/3;
            double x3 = -2*Math.sqrt(QQ)* Math.cos((theta-2*Math.PI)/3) - a/3;
            ArrayList<Double> roots = new ArrayList<>();
            roots.add(x1);
            roots.add(x2);
            roots.add(x3);
            return Collections.max(roots);

        } else {//1 real solution (small wave height)
            double AA = -1*sign(RR)*Math.pow((Math.abs(RR) + Math.sqrt(DD)),(1/3.0));
            double BB = 0;
            if(AA != 0){
                BB = QQ/AA;
            }

            double x1 = AA+BB - a/3;
            Complex tmp = new Complex(0,1).multiply(Math.sqrt(3)).divide(2).multiply((AA-BB));
            Complex x2 = new Complex(-1*(AA+BB)/2 -a/3,0).add(tmp);
            Complex x3 = x2.conjugate();
            ArrayList<Double> roots = new ArrayList<>();
            roots.add(x1);
            roots.add(x2.getReal());
            roots.add(x3.getReal());
            return Collections.max(roots);
        }
    }

    //This function return zero of the function defined in ship_resitance.m file
    //k3*x^(3+n_exp)/v_max_ms^n_exp + k2 * x^2 +k0
    //it's an approssimation of fzero matlab function.
    public static double Newton(double k3, double k2, double k0, long n_exp, double v_max_ms, double x0){
        double xk = x0;
        int k=0;
        double delta_ass = 0.0000001;//Max absolute error threshold
        int kmax = 1000000; //Max iterations number
        double fxk = k3 * Math.pow(xk, 3+n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk, 2) + k0;
        //Double fprimoxk = k3*(3+n_exp)*Math.pow(xk,(3+n_exp)-1)/Math.pow(v_max_ms, n_exp) + k2*2*xk;
        double fprimoxk=2*k2*xk+(3+n_exp)*(1/Math.pow(v_max_ms,n_exp))*k3*Math.pow(xk,(2+n_exp));
        double ck= -fxk/fprimoxk;
        while( (Math.abs(ck) > delta_ass*Math.ulp(1.0)) && (k<kmax) ){
            xk+=ck;
            fxk = k3 * Math.pow(xk, 3+n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk, 2) + k0;
            fprimoxk = k3*(3+n_exp)*Math.pow(xk,(3+n_exp)-1)/Math.pow(v_max_ms, n_exp) + k2*2*xk;
            ck=-fxk/fprimoxk;
            k++;
        }
        return xk;
    }


    private static double bisection(double a, double b, double k3, double k2, double k0, long n_exp, double v_max_ms){
        //Solve f(x) = 0 with recoursive bisection method.
        double delta_ass = 0.001;//Max absolute error threshold
        if(Math.abs((b-a)) <= (delta_ass + Math.max(Math.abs(a),Math.abs(b)))){//base case
            return (a+b)/2;//Function root is approssimated with the middle point of the range
        }
        else{//recoursive function calls
            double fa = k3 * Math.pow(a, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(a, 2) + k0;//Function evaluation in 'a'
            double meanPoint = (a+b)/2;//Mean point of the range [a, b]
            double fmean = k3 * Math.pow(meanPoint, n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(meanPoint, 2) + k0; // Function evaluation in 'meanPoint'
            if((fa*fmean)<0){//if sign(fa) != sign(fmean)
                return bisection(a, meanPoint, k3, k2, k0, n_exp, v_max_ms);//Root is in the range [a, meanPoint]
            }else {
                return bisection(meanPoint, b, k3, k2, k0, n_exp, v_max_ms);//Root is in the range [meanPoint, b]
            }
        }
    }

    public static double fzero_secant(double k3, double k2, double k0, long n_exp, double v_max_ms, double x0){
        //This function return zeros of the function defined in ship_resitance.m file
        //k3*x^(3+n_exp)/v_max_ms^n_exp + k2 * x^2 +k0
        //This function uses secant method to find function root.
        return secant(x0, (x0+ 0.5), k3, k2, k0, n_exp, v_max_ms);
    }

    private static double secant(double x0, double x1, double k3, double k2, double k0, long n_exp, double v_max_ms){
        //Solve f(x) = 0 with secant method.
        double delta_ass = 0.0000001;//Max absolute error threshold
        int kmax = 10000000; //Max iterations number
        int k = 1;
        //Initial xk and xk+1
        double xk = x0;
        double xk1 = x1;
        double root = xk;
        //Function evaluation in xk and xk+1
        double fxk = k3 * Math.pow(xk, 3+n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk, 2) + k0;
        double fxk1 = k3 * Math.pow(xk1, 3+n_exp)/Math.pow(v_max_ms, n_exp) + k2 * Math.pow(xk1, 2) + k0;
        //line gradient approssimation
        double pxk = (fxk-fxk1)/(xk-xk1);
        //k-step correction
        double ck = -(fxk/pxk);
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
        double[][] X = new double[nRows][nCols];
        double[][] Y = new double[nRows][nCols];
        for(int i =0;i<nRows;i++){
            for(int j=0;j<nCols; j++){
                X[i][j] = x.get(j);
                Y[i][j] = y.get(i);
            }
        }
        return new meshgridResults(X,Y);
    }


    public static double[] reshape(double[][] A, int dim){// Reshape 1
        if(A.length*A[0].length != dim){
            System.out.println("size(A) must be = dim!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("reshape: size(A) must be = dim!");
            debug.CloseFile();
            System.exit(0);
        }
        double[] output = new double[dim];
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

    public static double[][] reshape(ArrayList<Double> A, int[] dim){//reshape 2
        int nRows = dim[0];
        int nCols = dim[1];
        if(A.size() != nRows*nCols){
            System.out.println("size(A) must be = to nRows*nCols!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("reshape: size(A) must be = to nRows*nCols!");
            debug.CloseFile();
            System.exit(0);
        }
        double[][] out = new double[nRows][nCols];
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

    public static double[][] reshape(double[] A, int[] dim){//reshape 3
        int nRows = dim[0];
        int nCols = dim[1];
        if(A.length != nRows*nCols){
            System.out.println("A.length must be = to nRows*nCols!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("reshape: A.length must be = to nRows*nCols!");
            debug.CloseFile();
            System.exit(0);
        }
        double[][] out = new double[nRows][nCols];
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

    public static double[][] reshape(double[][] A, int[] dim){//Reshape 4
        int nRows = dim[0];
        int nCols = dim[1];
        if(A.length*A[0].length != nRows*nCols){
            System.out.println("size(A) must be = to nRows*nCols!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("reshape: size(A) must be = to nRows*nCols!");
            debug.CloseFile();
            System.exit(0);
        }
        double[][] out = new double[nRows][nCols];
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


    public static double[][] transposeMatrix(double[][] matrix){
        double[][] output = new double[matrix[0].length][matrix.length];
        for(int i=0;i< matrix.length;i++){
            for(int j=0;j<matrix[0].length;j++){
                output[j][i]=matrix[i][j];
            }
        }
        return output;
    }

    public static double[][] MatrixComponentXcomponent(double[][] a, double[][] b){
        int nRows = a.length;
        int nCols = a[0].length;
        if(a.length!=b.length){
            System.out.println("a.length != b.length");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("MatrixComponentXcomponent: a.length != b.length");
            debug.CloseFile();
            System.exit(0);
        }
        if(a[0].length!=b[0].length){
            System.out.println("a[0].length != b[0].length");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("MatrixComponentXcomponent: a[0].length != b[0].length");
            debug.CloseFile();
            System.exit(0);
        }
        double[][] output = new double[nRows][nCols];
        for(int i=0;i<nRows;i++){
            for(int j=0;j<nCols;j++){
                output[i][j]=a[i][j]*b[i][j];
            }
        }
        return  output;
    }

    public static double min(double[] array){
        double[] notNaNs = removeNaNs(array);
        if(notNaNs.length == 0){
            return Double.NaN;
        }
        else{
            double min=notNaNs[0];
            for(int i=1;i<notNaNs.length-1;i++){
                if(notNaNs[i]<min){
                    min = notNaNs[i];
                }
            }
            return min;
        }
    }

    private static double[] removeNaNs(double[] in){
        ArrayList<Double> notNaNs = new ArrayList<>();
        for(int i=0;i<in.length;i++){
            if(!Double.isNaN(in[i])){
                notNaNs.add(in[i]);
            }
        }
        double[] out = new double[notNaNs.size()];
        for(int i=0;i<notNaNs.size();i++){
            out[i] = notNaNs.get(i);
        }
        return out;
    }

    public static double[] min(double[][] matrix, int flag){
        //MATLAB min(x,[],y) implementation
        if(flag==1){
            //min(matrix,[],1) (min cols)
            double[] out = new double[matrix[0].length];
            for(int i=0;i<matrix[0].length;i++){
                double minTmp = matrix[0][i];
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
            double[] out = new double[matrix.length];
            for(int i=0;i<matrix.length;i++){
                double minTmp = matrix[i][0];
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
