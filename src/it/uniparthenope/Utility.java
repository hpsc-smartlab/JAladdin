package it.uniparthenope;


import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;

import flanagan.interpolation.BiCubicSplineFast;
import it.uniparthenope.Boxing.*;
import it.uniparthenope.Debug.MyFileWriter;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.joda.time.DateTime;
import org.joda.time.Days;
import org.joda.time.Period;

public class Utility {

    public static double[] arrayfy(ArrayList<Double> in){
        double[] out = new double[in.size()];
        for(int i=0;i<out.length;i++)
            out[i] = in.get(i);
        return out;
    }

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
        double base = 10.0;
        for(double exp: exps){
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

    public static double[] NaN(int dim){
        double[] out = new double[dim];
        for(int i=0;i<dim;i++)
            out[i]=Double.NaN;
        return out;
    }

    public static double[][] NaN(int rows, int cols){
        double[][] out = new double[rows][cols];
        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                out[i][j]=Double.NaN;
        return out;
    }

    public static double[][][] NaN(int rows, int cols, int zDim){
        double[][][] out = new double[zDim][rows][cols];
        for(int k=0;k<zDim;k++)
            for(int i=0;i<rows;i++)
                for(int j=0;j<cols;j++)
                    out[k][i][j] = Double.NaN;
        return out;
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

    public static double[][][] NaN3Dmatrix(int x, int y, int z){
        double[][][] matrix3d = new double[x][y][z];
        for(int i=0;i<x;i++){
            for(int j=0;j<y;j++){
                for(int k=0; k<z; k++){
                    matrix3d[i][j][k] = Double.NaN;
                }
            }
        }
        return matrix3d;
    }

    public static double[][] NaN2Dmatrix(int x, int y){
        double[][] out = new double[x][y];
        for(int i=0;i<x;i++){
            for(int j=0;j<y;j++){
                out[i][j] = Double.NaN;
            }
        }
        return out;
    }

    public static double[] nthroot(double[] x, int n){
        double[] out = new double[x.length];
        for(int i=0;i<out.length;i++)
            out[i] = nthroot(x[i], n);
        return out;
    }

    public static double nthroot(double x, int n)
    {
        double threshold = 0.0001;
        if(x < 0)
        {
            System.err.println("nthroot: Negative!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("nthroot: Negative!");
            debug.CloseFile();
            System.exit(0);
            return -1;
        }
        if(x == 0)
            return 0;
        double x1 = x;
        double x2 = x / n;
        while (Math.abs(x1 - x2) > threshold)
        {
            x1 = x2;
            x2 = ((n - 1.0) * x2 + x / Math.pow(x2, n - 1.0)) / n;
        }
        return x2;
    }

    public static double[][][] ones3Dmatrix(int x, int y, int z){
        double[][][] matrix3d = new double[x][y][z];
        for(int i=0;i<x;i++){
            for(int j=0;j<y;j++){
                for(int k=0; k<z; k++){
                    matrix3d[i][j][k] = 1.0;
                }
            }
        }
        return matrix3d;
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

    public static double[] abs(double[] in){
        double[] out = new double[in.length];
        for(int i=0;i<out.length; ++i)
            out[i] = Math.abs(in[i]);
        return out;
    }

    public static meshgridResults meshgrid(double[][] x, double[][] y){
        //[X, Y] =meshgrid(x,y) returns 2-D grid coordinates based on the coordinates contained
        // in vectors x and y. X is a matrix where each row is a copy of x,
        // and Y is a matrix where each column is a copy of y.
        // The grid represented by the coordinates X and Y has length(y) rows and length(x) columns.
        int nRows = y.length*y[0].length;
        int nCols = x.length*x[0].length;
        double[][] X = new double[nRows][nCols];
        double[][] Y = new double[nRows][nCols];
        int nElem =0;
        for(int j=0;j<X[0].length;j++){
            for(int i=0;i<X.length;i++){
                X[i][j] = getElementById(x, nElem);
            }
            nElem++;
        }
        nElem=0;
        for(int i=0;i<Y.length;i++){
            for(int j=0;j<Y[0].length;j++){
                Y[i][j] = getElementById(y, nElem);
            }
            nElem++;
        }
        return new meshgridResults(X,Y);
    }

    private static double getElementById(double[][] matrix, int id){
        int nRows=matrix.length;
        int colIndex = ((int) Math.floor(id/nRows));
        int rowIndex = (id%nRows);
        return matrix[rowIndex][colIndex];
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

    public static meshgridResults meshgrid(double[] x, double[] y){
        //[X, Y] =meshgrid(x,y) returns 2-D grid coordinates based on the coordinates contained
        // in vectors x and y. X is a matrix where each row is a copy of x,
        // and Y is a matrix where each column is a copy of y.
        // The grid represented by the coordinates X and Y has length(y) rows and length(x) columns.
        int nRows = y.length;
        int nCols = x.length;
        double[][] X = new double[nRows][nCols];
        double[][] Y = new double[nRows][nCols];
        for(int i =0;i<nRows;i++){
            for(int j=0;j<nCols; j++){
                X[i][j] = x[j];
                Y[i][j] = y[i];
            }
        }
        return new meshgridResults(X,Y);
    }

    public static double[][] reshape(double[][][] A, int[] dim){
        int rows = dim[0];
        int cols = dim[1];
        if(A.length*A[0].length*A[0][0].length != rows*cols){
            System.out.println("size(A) must be = dim!");
            System.exit(0);
        }
        int contaRighe2D = 0;
        int contaCol2D = 0;
        double[][] out = new double[rows][cols];
        for(int k=0;k<A.length;k++){
            int contaRighe3D = 0;
            int contaCol3D = 0;
            for(int u=0;u<(A[0].length*A[0][0].length);u++){
                out[contaRighe2D][contaCol2D] = A[k][contaRighe3D][contaCol3D];
                contaRighe2D++;
                contaRighe3D++;
                if(contaRighe2D==rows){
                    contaRighe2D=0;
                    contaCol2D++;
                }
                if(contaRighe3D==A[0].length){
                    contaRighe3D=0;
                    contaCol3D++;
                }
            }
        }
        return out;
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

    public static int[][] reshape(int[] A, int[] dim){//reshape 3 for Long Int
        int nRows = dim[0];
        int nCols = dim[1];
        if(A.length != nRows*nCols){
            System.out.println("A.length must be = to nRows*nCols!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("reshape: A.length must be = to nRows*nCols!");
            debug.CloseFile();
            System.exit(0);
        }
        int[][] out = new int[nRows][nCols];
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

    public static findResults find(double[] A, String condition1, double b, String condition2, double c){//TODO: implement other cases
        ArrayList<Integer> indexes = new ArrayList<>();
        ArrayList<Double> elements = new ArrayList<>();

        switch(condition1){
            case("<="):
                switch(condition2){
                    case("<="):
                        for(int i=0;i<A.length;i++){
                            if(b <= A[i] && A[i] <= c){
                                elements.add(A[i]);
                                indexes.add(i);
                            }
                        }
                        break;
                }
                break;
            case(">="):
                switch(condition2){
                    case("<="):
                        for(int i=0;i<A.length;i++){
                            if(b >= A[i] && A[i] <= c){
                                elements.add(A[i]);
                                indexes.add(i);
                            }
                        }
                        break;
                }
                break;
        }

        int[] ind = new int[indexes.size()];
        double[] elem = new double[elements.size()];
        for(int i=0;i<indexes.size();i++){
            ind[i] = indexes.get(i);
            elem[i] = elements.get(i);
        }
        return new findResults(indexes.size(), elem, ind);

    }

    public static ArrayList<Integer> find(boolean[] in){
        ArrayList<Integer> _indexes = new ArrayList<>();
        for(int i=0;i<in.length;++i){
            if(in[i])
                _indexes.add(i);
        }
        return _indexes;
    }

    public static ArrayList<Integer> findNaNs(double[][] A){
        ArrayList<Integer> _indexes = new ArrayList<>();//if i%2 =0 row index, else col index
        for(int i=0;i<A.length;i++){
            for(int j=0;j<A[0].length;j++){
                if(Double.isNaN(A[i][j])){
                    _indexes.add(i);
                    _indexes.add(j);
                }
            }
        }
        return _indexes;
    }

    public static ArrayList<Integer> findNaNs(double[] A){
        ArrayList<Integer> _indexes = new ArrayList<>();
        for(int i=0;i<A.length;i++){
            if(Double.isNaN(A[i]))
                _indexes.add(i);
        }
        return _indexes;
    }

    public static findResults find(double[] A, String condition, double b){
        ArrayList<Integer> indexes = new ArrayList<>();
        ArrayList<Double> elements = new ArrayList<>();
        switch (condition){
            case(">"):
                for(int i=0;i<A.length;i++){
                    if(A[i]>b){
                        elements.add(A[i]);
                        indexes.add(i);
                    }
                }
            break;
            case("<"):
                for(int i=0;i<A.length;i++){
                    if(A[i]<b){
                        elements.add(A[i]);
                        indexes.add(i);
                    }
                }
                break;
            case("=="):
                for(int i=0;i<A.length;i++){
                    if(A[i]==b){
                        elements.add(A[i]);
                        indexes.add(i);
                    }
                }
                break;
            case(">="):
                for(int i=0;i<A.length;i++){
                    if(A[i]>=b){
                        elements.add(A[i]);
                        indexes.add(i);
                    }
                }
                break;
            case("<="):
                for(int i=0;i<A.length;i++){
                    if(A[i]<=b){
                        elements.add(A[i]);
                        indexes.add(i);
                    }
                }
                break;
        }

        int[] ind = new int[indexes.size()];
        double[] elem = new double[elements.size()];
        for(int i=0;i<indexes.size();i++){
            ind[i] = indexes.get(i);
            elem[i] = elements.get(i);
        }
        return new findResults(indexes.size(), elem, ind);
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

    public static boolean any(int[][] A, String condition, int threshold){
        boolean output = false;
        switch(condition){
            case(">"):
                for(int i=0;i<A.length;i++){
                    for(int j=0;j<A[0].length;j++){
                        if(A[i][j]>threshold){
                            output = true;
                        }
                    }
                }
                break;
            case("<"):
                for(int i=0;i<A.length;i++){
                    for(int j=0;j<A[0].length;j++){
                        if(A[i][j]<threshold){
                            output = true;
                        }
                    }
                }
                break;
            case("=="):
                for(int i=0;i<A.length;i++){
                    for(int j=0;j<A[0].length;j++){
                        if(A[i][j]==threshold){
                            output = true;
                        }
                    }
                }
                break;
            case(">="):
                for(int i=0;i<A.length;i++){
                    for(int j=0;j<A[0].length;j++){
                        if(A[i][j]>=threshold){
                            output = true;
                        }
                    }
                }
                break;
            case("<="):
                for(int i=0;i<A.length;i++){
                    for(int j=0;j<A[0].length;j++){
                        if(A[i][j]<=threshold){
                            output = true;
                        }
                    }
                }
                break;
        }
        return output;
    }

    public static double[][] squeeze(double[][][] mat3d, int dimToSqueeze, int index){
        int z_dim = mat3d.length;
        int rows = mat3d[0].length;
        int cols = mat3d[0][0].length;
        double[][] out = new double[0][0];
        switch (dimToSqueeze){
            case 0://Z
                //squeezing Zs
                if(index>=z_dim){
                    System.out.println("squeeze: index "+index+" exceed Z dimension ("+z_dim+")");
                    MyFileWriter debug = new MyFileWriter("","debug",false);
                    debug.WriteLog("squeeze: index "+index+" exceed Z dimension ("+z_dim+")");
                    debug.CloseFile();
                    System.exit(0);
                }
                out = new double[rows][cols];
                for(int i=0; i<mat3d[0].length; i++){
                    for(int j=0; j<mat3d[0][0].length; j++){
                        out[i][j] = mat3d[index][i][j];
                    }
                }
                break;
            case 1://rows:
                //squeezing rows:
                if(index>=rows){
                    System.out.println("squeeze: index "+index+" exceed rows dimension ("+rows+")");
                    MyFileWriter debug = new MyFileWriter("","debug",false);
                    debug.WriteLog("squeeze: index "+index+" exceed rows dimension ("+rows+")");
                    debug.CloseFile();
                    System.exit(0);
                }
                out = new double[cols][z_dim];
                for(int i = 0; i<cols; i++){
                    for(int j = 0; j<z_dim; j++){
                        out[i][j] = mat3d[j][index][i];
                    }
                }
                break;
            case 2://cols
                //squeezing cols:
                if(index>=cols){
                    System.out.println("squeeze: index "+index+" exceed cols dimension ("+cols+")");
                    MyFileWriter debug = new MyFileWriter("","debug",false);
                    debug.WriteLog("squeeze: index "+index+" exceed cols dimension ("+cols+")");
                    debug.CloseFile();
                    System.exit(0);
                }
                out = new double[rows][z_dim];
                for(int i=0;i<rows;i++){
                    for(int j=0;j<z_dim;j++){
                        out[i][j] = mat3d[j][i][index];
                    }
                }
                break;
            default:
                System.out.println("squeeze: "+dimToSqueeze+" is not a valid dimension");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("squeeze: "+dimToSqueeze+" is not a valid dimension");
                debug.CloseFile();
                System.exit(0);
                break;
        }
        return out;
    }


    public static double[][] squeeze(double[][][] mat3d){
        int z_dim = mat3d.length;
        int rows = mat3d[0].length;
        int cols = mat3d[0][0].length;
        double[][] out = new double[0][0];
        //check for singleton dimension:
        if(z_dim==1){
            //3rd dimension singleton
            out = new double[rows][cols];
            for(int i=0; i<mat3d[0].length; i++){
                for(int j=0; j<mat3d[0][0].length; j++){
                    out[i][j] = mat3d[0][i][j];
                }
            }
        } else {
            //rows singleton
            if(rows == 1){
                out = new double[cols][z_dim];
                for(int i = 0; i<cols; i++){
                    for(int j = 0; j<z_dim; j++){
                        out[i][j] = mat3d[j][0][i];
                    }
                }
            } else {
                //cols singleton
                if(cols == 1){
                    out = new double[rows][z_dim];
                    for(int i=0;i<rows;i++){
                        for(int j=0;j<z_dim;j++){
                            out[i][j] = mat3d[j][i][0];
                        }
                    }
                } else {
                    System.out.println("squeeze: no singleton dimension found!");
                    MyFileWriter debug = new MyFileWriter("","debug",false);
                    debug.WriteLog("squeeze: no singleton dimension found!");
                    debug.CloseFile();
                    System.exit(0);
                }
            }
        }
        return out;
    }

    public static int rank(double[][] matrix){
        RealMatrix m = MatrixUtils.createRealMatrix(matrix);
        SingularValueDecomposition x = new SingularValueDecomposition(m);
        return x.getRank();
    }


    public static boolean[][] isnan(double[][] in){
        boolean[][] out = new boolean[in.length][in[0].length];
        for(int i=0;i<in.length;i++){
            for(int j=0;j<in[0].length; j++){
                if(Double.isNaN(in[i][j])){
                    out[i][j] = true;
                } else{
                    out[i][j] = false;
                }
            }
        }
        return out;
    }

    public static double[][] convertDouble(boolean[][] in){
        double[][] out = new double[in.length][in[0].length];
        for(int i=0;i<in.length;i++){
            for(int j=0;j<in[0].length; j++){
                if(in[i][j]){
                    out[i][j] = 1.0;
                } else{
                    out[i][j] = 0.0;
                }
            }
        }
        return out;
    }

    public static double[] sum(double[][] in){
        double[] out = new double[in[0].length];
        for(int i=0;i<in.length; i++){
            for(int j=0;j<in[0].length;j++){
                out[j]+=in[i][j];
            }
        }
        return out;
    }


    public static ArrayList<Integer> find(double[][] mat, double element){
        ArrayList<Integer> indexes = new ArrayList<>();
        int count=0;
        boolean checkNaN = Double.isNaN(element);
        for(int i=0;i<mat[0].length;i++){
            for(int j=0;j<mat.length;j++){
                if(checkNaN){
                    if(Double.isNaN(mat[j][i])) {
                        indexes.add(count);
                    }
                } else{
                    if(mat[j][i] == element) {
                        indexes.add(count);
                    }
                }
                count++;
            }
        }
        return  indexes;
    }

    public static int[] removeElements(int[] array, int element){
        ArrayList<Integer> indexes = find(array,element);
        int[] out = new int[indexes.size()];
        for(int i=0;i<indexes.size();i++){
            out[i]=array[indexes.get(i)];
        }
        return out;
    }

    private static ArrayList<Integer> find(int[] array, int NotElement){
        ArrayList<Integer> indexes = new ArrayList<>();
        for(int i=0;i<array.length;i++){
            if(array[i]!=NotElement){
                indexes.add(i);
            }
        }
        return indexes;
    }


    public static double[][] interp1(double[] x, double[][] y, double[] xi){
        //wrapper for the interp1 method
        double[][] yi = new double[y.length][y[0].length];
        for(int row=0;row<y.length; row++){
            //splitting y in nRows arrays
            double[] Y = new double[y[0].length];
            for(int i=0;i<y[0].length;i++){
                Y[i]=y[row][i];
            }
            //calculating y-row and adding it to result
            double[] yirow = interp1(x,Y,xi);
            for(int i=0;i<y[0].length;i++){
                yi[row][i]=yirow[i];
            }
        }
        return yi;
    }

    public static double[][][] interp1(double[] x, double[][][] y, int[] xi){
        //wrapper for the interp1 method
        double[] xiDouble = new double[xi.length];
        for(int i=0;i<xiDouble.length;i++)
            xiDouble[i]=xi[i]+0.0;
        return interp1(x, y, xiDouble);
    }

    public static double interp1(ArrayList<Double> x, double[] y, double xi){
        //wrapper for the interp1 method
        double[] xArray = new double[x.size()];
        for(int i=0;i<x.size();i++){
            xArray[i]=x.get(i);
        }
        double[] xiArray = new double[1];
        xiArray[0]=xi;
        double[] yi = interp1(xArray, y, xiArray);
        return yi[0];
    }

    public static double[][][] interp1(double[] x, double[][][] y, double[] xi){
        //wrapper for the interp1 method
        int Z = y.length;
        int rows = y[0].length;
        int cols = y[0][0].length;
        double[][][] out = new double[Z][xi.length][cols];
        //cols extraction:
        for(int z=0;z<Z;z++){
            for(int j=0;j<cols;j++){
                double[] row = new double[rows];
                for(int i=0;i<rows;i++){
                    row[i]=y[z][i][j];
                }
                double[] rowOut = interp1(x,row,xi);
                for(int i=0;i<rowOut.length;i++){
                    out[z][i][j] = rowOut[i];
                }
            }
        }
        return out;
    }


    public static double[] interp1(double[] x, double[] y, double[] xi, double extrapolation){
        //interp1 MATLAB function porting
        if(x.length != y.length){
            System.out.println("interp1: X and Y must be the same length!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("interp1: X and Y must be the same length!");
            debug.CloseFile();
            System.exit(0);
        }
        if(x.length == 1){
            System.out.println("interp1: X must contain more than one value");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("interp1: X must contain more than one value");
            debug.CloseFile();
            System.exit(0);
        }

        double[] dx = new double[x.length - 1];
        double[] dy = new double[x.length - 1];
        double[] slope = new double[x.length - 1];
        double[] intercept = new double[x.length - 1];

        // Calculate the line equation (i.e. slope and intercept) between each point
        for (int i = 0; i < x.length - 1; i++) {
            dx[i] = x[i + 1] - x[i];
            if (dx[i] == 0) {
                System.out.println("interp1: X must be montotonic. A duplicate \" + \"x-value was found");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("interp1: X must be montotonic. A duplicate \" + \"x-value was found");
                debug.CloseFile();
                System.exit(0);
            }
            if (dx[i] < 0) {
                System.out.println("interp1: X must be sorted");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("interp1: X must be sorted");
                debug.CloseFile();
                System.exit(0);
            }
            dy[i] = y[i + 1] - y[i];
            slope[i] = dy[i] / dx[i];
            intercept[i] = y[i] - x[i] * slope[i];
        }

        // Perform the interpolation here
        double[] yi = new double[xi.length];
        for (int i = 0; i < xi.length; i++) {
            if ((xi[i] > x[x.length - 1]) || (xi[i] < x[0])) {
                yi[i] = extrapolation;
            }
            else {
                if(Double.isNaN(xi[i])){
                    yi[i] = Double.NaN;
                } else {
                    int loc = Arrays.binarySearch(x, xi[i]);
                    if (loc < -1) {
                        loc = -loc - 2;
                        try {
                            yi[i] = slope[loc] * xi[i] + intercept[loc];
                        } catch (Exception e) {
                            System.exit(-1);
                        }
                    } else {
                        yi[i] = y[loc];
                    }
                }
            }
        }
        return yi;
    }

    /**NOTE:
     * Arrays.binarySearch method returns:
     * index -> if the item is present in the collection
     * -(index) -1 -> otherwise*/

    public static double[] interp1(double[] x, double[][] y,int yColIdx, double[][] xi, int xiColIdx){
        //interp1 MATLAB function porting
        if(x.length != y.length){
            System.out.println("interp1: X and Y must be the same length!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("interp1: X and Y must be the same length!");
            debug.CloseFile();
            System.exit(0);
        }
        if(x.length == 1){
            System.out.println("interp1: X must contain more than one value");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("interp1: X must contain more than one value");
            debug.CloseFile();
            System.exit(0);
        }

        double[] dx = new double[x.length - 1];
        double[] dy = new double[x.length - 1];
        double[] slope = new double[x.length - 1];
        double[] intercept = new double[x.length - 1];

        // Calculate the line equation (i.e. slope and intercept) between each point
        for (int i = 0; i < x.length - 1; i++) {
            dx[i] = x[i + 1] - x[i];
            if (dx[i] == 0) {
                System.out.println("interp1: X must be montotonic. A duplicate \" + \"x-value was found");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("interp1: X must be montotonic. A duplicate \" + \"x-value was found");
                debug.CloseFile();
                System.exit(0);
            }
            if (dx[i] < 0) {
                System.out.println("interp1: X must be sorted");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("interp1: X must be sorted");
                debug.CloseFile();
                System.exit(0);
            }
            dy[i] = y[i + 1][yColIdx] - y[i][yColIdx];
            slope[i] = dy[i] / dx[i];
            intercept[i] = y[i][yColIdx] - x[i] * slope[i];
        }

        // Perform the interpolation here
        double[] yi = new double[xi.length];
        for (int i = 0; i < xi.length; i++) {
            if ((xi[i][xiColIdx] > x[x.length - 1]) || (xi[i][xiColIdx] < x[0])) {
                yi[i] = Double.NaN;
            }
            else {
                if(Double.isNaN(xi[i][xiColIdx])){
                    yi[i] = Double.NaN;
                } else {
                    int loc = Arrays.binarySearch(x, xi[i][xiColIdx]);
                    if (loc < -1) {
                        loc = -loc - 2;
                        yi[i] = slope[loc] * xi[i][xiColIdx] + intercept[loc];
                    } else {
                        yi[i] = y[loc][yColIdx];
                    }
                }
            }
        }
        return yi;
    }


    public static double[] interp1(double[] x, double[] y, double[] xi){
        //interp1 MATLAB function porting
        if(x.length != y.length){
            System.out.println("interp1: X and Y must be the same length!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("interp1: X and Y must be the same length!");
            debug.CloseFile();
            System.exit(0);
        }
        if(x.length == 1){
            System.out.println("interp1: X must contain more than one value");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("interp1: X must contain more than one value");
            debug.CloseFile();
            System.exit(0);
        }

        double[] dx = new double[x.length - 1];
        double[] dy = new double[x.length - 1];
        double[] slope = new double[x.length - 1];
        double[] intercept = new double[x.length - 1];

        // Calculate the line equation (i.e. slope and intercept) between each point
        for (int i = 0; i < x.length - 1; i++) {
            dx[i] = x[i + 1] - x[i];
            if (dx[i] == 0) {
                System.out.println("interp1: X must be montotonic. A duplicate \" + \"x-value was found");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("interp1: X must be montotonic. A duplicate \" + \"x-value was found");
                debug.CloseFile();
                System.exit(0);
            }
            if (dx[i] < 0) {
                System.out.println("interp1: X must be sorted");
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("interp1: X must be sorted");
                debug.CloseFile();
                System.exit(0);
            }
            dy[i] = y[i + 1] - y[i];
            slope[i] = dy[i] / dx[i];
            intercept[i] = y[i] - x[i] * slope[i];
        }

        // Perform the interpolation here
        double[] yi = new double[xi.length];
        for (int i = 0; i < xi.length; i++) {
            if ((xi[i] > x[x.length - 1]) || (xi[i] < x[0])) {
                yi[i] = Double.NaN;
            }
            else {
                if(Double.isNaN(xi[i])){
                    yi[i] = Double.NaN;
                } else {
                    int loc = Arrays.binarySearch(x, xi[i]);
                    if (loc < -1) {
                        loc = -loc - 2;
                        yi[i] = slope[loc] * xi[i] + intercept[loc];
                    } else {
                        yi[i] = y[loc];
                    }
                }
            }
        }
        return yi;
    }


    public static double[][] interp2(double[][] X, double[][] Y, double[][] Z, double[] Xq, double[] Yq){//wrapper for interp2 method
        double[] Xarray = new double[X[0].length];
        for(int i=0;i<X[0].length;i++)
            Xarray[i]=X[0][i];
        double[] Yarray = new double[Y.length];
        for(int i=0;i<Y.length;i++)
            Yarray[i]=Y[i][0];
        return interp2(Xarray, Yarray, Z, Xq, Yq);
    }


    public static double[][] interp2(double[] X, double[] Y, double[][] Z, double[] Xq, double[] Yq){//interp2 method implemented with 'linear'

        double[][] out = new double[Yq.length][Xq.length];
        double[][] Zt = Utility.transposeMatrix(Z);
        for(int i=0;i<Yq.length;++i){
            for(int j=0;j<Xq.length;++j){
                int x_idx=0, x1_idx=0, x2_idx=0 , y_idx=0, y1_idx=0, y2_idx=0;
                if( (Xq[j] > X[X.length-1]) || (Xq[j] < X[0]) || (Yq[i] > Y[Y.length-1]) || (Yq[i] < Y[0]) ){ //Xq o Yq > max o < min
                    out[i][j] = Double.NaN;
                } else {
                    if(Double.isNaN(Xq[j]) || Double.isNaN(Yq[i])){
                        out[i][j] = Double.NaN;
                    } else {
                        x_idx = Arrays.binarySearch(X, Xq[j]);
                        if (x_idx < -1) { //element not found
                            x1_idx = -x_idx - 2;
                            x2_idx = -x_idx - 1;
                        } else { //element found
                            if (x_idx == 0) {//1st element of the array
                                x1_idx = x_idx;
                                x2_idx = x_idx + 1;
                            } else {//last element of the array (or between)
                                x1_idx = x_idx - 1;
                                x2_idx = x_idx;
                            }
                        }
                        y_idx = Arrays.binarySearch(Y, Yq[i]);
                        if (y_idx < -1) { //element not found
                            y1_idx = -y_idx - 2;
                            y2_idx = -y_idx - 1;
                        } else { //element found
                            if (y_idx == 0) { //1st element of the array
                                y1_idx = y_idx;
                                y2_idx = y_idx + 1;
                                //last element of the array (or between)
                            } else {
                                y1_idx = y_idx - 1;
                                y2_idx = y_idx;
                            }
                        }

                        double fx_y1 = (((X[x2_idx] - Xq[j]) / (X[x2_idx] - X[x1_idx])) * Zt[x1_idx][y1_idx]) + (((Xq[j] - X[x1_idx]) / (X[x2_idx] - X[x1_idx])) * Zt[x2_idx][y1_idx]);
                        double fx_y2 = (((X[x2_idx] - Xq[j]) / (X[x2_idx] - X[x1_idx])) * Zt[x1_idx][y2_idx]) + (((Xq[j] - X[x1_idx]) / (X[x2_idx] - X[x1_idx])) * Zt[x2_idx][y2_idx]);

                        out[i][j] = ((Y[y2_idx] - Yq[i]) / (Y[y2_idx] - Y[y1_idx]) * fx_y1) + ((Yq[i] - Y[y1_idx]) / (Y[y2_idx] - Y[y1_idx]) * fx_y2);
                    }
                }
            }
        }
        return out;
    }


    //OLD INTERP2, 'cubic' implemented with Flangan's Math library

//    public static double[][] interp2(double[][] X, double[][] Y, double[][] Z, double[] Xq, double[] Yq){
//        double[] Xarray = new double[X[0].length];
//        for(int i=0;i<X[0].length;i++)
//            Xarray[i]=X[0][i];
//        double[] Yarray = new double[Y.length];
//        for(int i=0;i<Y.length;i++)
//            Yarray[i]=Y[i][0];
//
//        double[][] Zt = transposeMatrix(Z);
//        BiCubicSplineFast bcs = new BiCubicSplineFast(Xarray,Yarray, Zt);
//        double[][] out = new double[Yq.length][Xq.length];
//        for(int i=0;i<Yq.length;i++){
//            for(int j=0;j<Xq.length;j++){
//                out[i][j] = bcs.interpolate(Xq[j],Yq[i]);
//            }
//        }
//        return out;
//    }

    public static double[][] transposeMatrix(double[][] matrix){
        double[][] output = new double[matrix[0].length][matrix.length];
        for(int i=0;i< matrix.length;i++){
            for(int j=0;j<matrix[0].length;j++){
                output[j][i]=matrix[i][j];
            }
        }
        return output;
    }

    public static double[][][] transpose3DMatrix(double[][][] in){//inverts rows with cols
        double[][][] out = new double[in.length][in[0][0].length][in[0].length];
        for(int i=0;i<in.length;i++){
            out[i]=transposeMatrix(in[i]);
        }
        return out;
    }

    public static double[][][] permute3DColsWithZdim(double[][][] in){
        double[][][] out = new double[in[0][0].length][in[0].length][in.length];
        for(int i=0;i<out[0].length;++i){
            for(int z=0;z<out.length;++z){
                for(int j=0;j<out[0][0].length;++j){
                    out[z][i][j] = in[j][i][z];
                }
            }
        }
        return out;
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

    public static double min(double a, double b){
        if(Double.isNaN(a) && Double.isNaN(b)){
            return Double.NaN;
        }else{
            if(Double.isNaN(a)){
                return b;
            } else {
                if(Double.isNaN(b)){
                    return a;
                } else{
                    if(a<b){
                        return a;
                    } else {
                        return b;
                    }
                }
            }
        }
    }

    public static int max(int[][] matrix, int colID){
        //returns the max element of the column colID
        int max = matrix[0][colID];
        for(int i=1;i<matrix.length;i++){
            if(matrix[i][colID]>max)
                max = matrix[i][colID];
        }
        return max;
    }

    public static minResults minWithIndex(double[] array){
        //check for NaNs
        for(int i=0;i<array.length; i++){
            if(Double.isNaN(array[i])){
                array[i] = Double.MAX_VALUE;
            }
        }
        //found min and indexMin
        int index = 0;
        double min = array[0];
        for(int i=1;i<array.length-1;i++){
            if(array[i]<array[i+1]){
                index = i;
                min = array[i];
            }
        }
        return new minResults(min, index);
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

    public static double[][] mean3d(double[][][] mat3d){
        int dim0 = mat3d.length;
        int dim1 = mat3d[0].length;
        int dim2 = mat3d[0][0].length;
        double[][] output = new double[dim0][dim1];
        for(int i=0;i<mat3d.length;i++){
            for(int j=0;j<mat3d[0].length; j++){
                double mean = 0.0;
                for(int k=0;k<mat3d[0][0].length; k++){
                    mean += mat3d[i][j][k];
                }
                mean = mean/dim2;
                output[i][j] = mean;
            }
        }
        return output;
    }

    public static double min3d(double[][][] mat3d){
        double min = Double.MAX_VALUE;
        for(int i=0;i<mat3d.length;i++){
            for(int j=0;j<mat3d[0].length;j++){
                for(int k=0; k<mat3d[0][0].length;k++){
                    if(mat3d[i][j][k]<min && !Double.isNaN(mat3d[i][j][k])){
                        min = mat3d[i][j][k];
                    }
                }
            }
        }
        return min;
    }

    public static double min2d(double[][] mat){
        double min = Double.MAX_VALUE;
        for(int i=0;i<mat.length;i++){
            for(int j=0;j<mat[0].length;j++){
                if(mat[i][j]<min && !Double.isNaN(mat[i][j])){
                    min = mat[i][j];
                }
            }
        }
        return min;
    }

    public static double max3d(double[][][] mat3d){
        double max = Double.MIN_VALUE;
        for(int i=0;i<mat3d.length;i++){
            for(int j=0;j<mat3d[0].length;j++){
                for(int k=0; k<mat3d[0][0].length;k++){
                    if(mat3d[i][j][k]>max && !Double.isNaN(mat3d[i][j][k])){
                        max = mat3d[i][j][k];
                    }
                }
            }
        }
        return max;
    }

    public static double max2d(double[][] mat){
        double max = Double.MIN_VALUE;
        for(int i=0;i<mat.length;i++){
            for(int j=0;j<mat[0].length;j++){
                if(mat[i][j]>max && !Double.isNaN(mat[i][j])){
                    max = mat[i][j];
                }
            }
        }
        return max;
    }

    public static double max(double[][] matrix){
        double max = matrix[0][0];
        for(int i=0;i<matrix.length;i++){
            for(int j=0;j<matrix[0].length; j++){
                if(matrix[i][j]>max){
                    max = matrix[i][j];
                }
            }
        }
        return max;
    }

    public static double maxNotNaN(double[][] matrix){
        double max = Double.MIN_VALUE;
        for(int i=0;i<matrix.length;i++){
            for(int j=0;j<matrix[0].length; j++){
                if(matrix[i][j]>max && !Double.isNaN(matrix[i][j])){
                    max = matrix[i][j];
                }
            }
        }
        return max;
    }

    public static int max(int[] array){
        int max = array[0];
        for(int i=1;i<array.length;i++){
            if(array[i]> max){
                max = array[i];
            }
        }
        return max;
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



    public static double min(double[][] matrix){
        double min = matrix[0][0];
        for(int i=0;i<matrix.length;i++){
            for(int j=0;j<matrix[0].length; j++){
                if(matrix[i][j]<min){
                    min = matrix[i][j];
                }
            }
        }
        return min;
    }
    public static double minNotNaN(double[][] matrix){
        double min = Double.MAX_VALUE;
        for(int i=0;i<matrix.length;i++){
            for(int j=0;j<matrix[0].length; j++){
                if(matrix[i][j]<min && !Double.isNaN(matrix[i][j])){
                    min = matrix[i][j];
                }
            }
        }
        return min;
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

    public static boolean[] ismember(int[] elements, int[] array){
        boolean[] out = new boolean[elements.length];
        for(int i=0;i<out.length;i++)
            out[i]=ismebmer(elements[i],array);
        return out;
    }

    public static boolean ismebmer(int element, int[] array){
        boolean found = false;
        int i=0;
        while(!found && i<array.length){
            if(array[i]==element)
                found=true;
            i++;
        }
        return found;
    }

    public static int[] cumsum(boolean[] array){
        int value = 0;
        int[] out = new int[array.length];
        for(int i=0;i<out.length;i++){
            if(array[i])
                value++;
            out[i]=value;
        }
        return out;
    }

    public static long numel(double[][][] mat){
        return mat.length*mat[0].length*mat[0][0].length;
    }

    public static int numel(double[][] mat){ return mat.length* mat[0].length;}

    //MATLAB inpolygon implementation
    // Given three colinear points p, q, r, the function checks if
    // point q lies on line segment 'pr'
    private static boolean onSegment(Point p, Point q, Point r)
    {
        if (q.getX() <= Math.max(p.getX(), r.getX()) && q.getX() >= Math.min(p.getX(), r.getX()) &&
                q.getY() <= Math.max(p.getY(), r.getY()) && q.getY() >= Math.min(p.getY(), r.getY()))
            return true;
        return false;
    }

    // To find orientation of ordered triplet (p, q, r).
    // The function returns following values
    // 0 --> p, q and r are colinear
    // 1 --> Clockwise
    // 2 --> Counterclockwise
    private static int orientation(Point p, Point q, Point r)
    {
        int val = new Double(Math.floor((q.getY() - p.getY()) * (r.getX() - q.getX()) - (q.getX() - p.getX()) * (r.getY() - q.getY()))).intValue();

        if (val == 0) return 0;  // colinear
        return (val > 0)? 1: 2; // clock or counterclock wise
    }

    // The function that returns true if line segment 'p1q1'
    // and 'p2q2' intersect.
    private static boolean doIntersect(Point p1, Point q1, Point p2, Point q2)
    {
        // Find the four orientations needed for general and
        // special cases
        int o1 = orientation(p1, q1, p2);
        int o2 = orientation(p1, q1, q2);
        int o3 = orientation(p2, q2, p1);
        int o4 = orientation(p2, q2, q1);

        // General case
        if (o1 != o2 && o3 != o4)
            return true;

        // Special Cases
        // p1, q1 and p2 are colinear and p2 lies on segment p1q1
        if (o1 == 0 && onSegment(p1, p2, q1)) return true;

        // p1, q1 and p2 are colinear and q2 lies on segment p1q1
        if (o2 == 0 && onSegment(p1, q2, q1)) return true;

        // p2, q2 and p1 are colinear and p1 lies on segment p2q2
        if (o3 == 0 && onSegment(p2, p1, q2)) return true;

        // p2, q2 and q1 are colinear and q1 lies on segment p2q2
        if (o4 == 0 && onSegment(p2, q1, q2)) return true;

        return false; // Doesn't fall in any of the above cases
    }

    public static inpolygonResults inpolygon(double xP, double yP, ArrayList<Double> xPolygon, ArrayList<Double> yPolygon){
        inpolygonResults result = new inpolygonResults();
        if(xPolygon.size() != yPolygon.size()){
            System.out.println("xPolygon.size != yPolygon.size");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("xPolygon.size != yPolygon.size");
            debug.CloseFile();
            System.exit(0);
        }
        if(xPolygon.size() < 3){
            System.out.println("there must be at least 3 vertices in a polygon!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("there must be at least 3 vertices in a polygon!");
            debug.CloseFile();
            System.exit(0);
        }
        Point p = new Point(xP,yP);
        Point extreme = new Point(10000000000.0, p.getY());
        Point[] polygon = new Point[xPolygon.size()];
        for(int i=0;i<polygon.length;i++){
            polygon[i] = new Point(xPolygon.get(i),yPolygon.get(i));
        }
        int count = 0;
        int i = 0;
        do{
            int next = (i+1)%xPolygon.size();
            if(doIntersect(polygon[i], polygon[next], p, extreme)){
                if (orientation(polygon[i], p, polygon[next]) == 0){
                    if(onSegment(polygon[i], p, polygon[next])){
                        result.setOn(true);
                        i=0;
                    }
                }
                count++;
            }
            i = next;
        }while(i !=0 );
        if(!result.getOn()){
            if((count & 1)==1){
                result.setIn(true);
            }
        }
        return  result;
    }

    public static double[][] atan2(double[][] a, double[][] b){
        if((a.length != b.length) || (a[0].length != b[0].length)){
            System.out.println("Utility.atan2: A and B must have same size!");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("Utility.atan2: A and B must have same size!");
            debug.CloseFile();
            System.exit(0);
        }
        double[][] out = new double[a.length][a[0].length];
        for(int i =0 ;i<a.length;i++){
            for(int j=0; j<a[0].length;j++){
                out[i][j] = Math.atan2(a[i][j],b[i][j]);
            }
        }
        return out;
    }

    public static long datenum(String dateString, String format){
        SimpleDateFormat sdf = new SimpleDateFormat(format);
        try {
            Date data = sdf.parse("000001010000");
            DateTime from = new DateTime(data.getTime());
            data = sdf.parse(dateString);
            DateTime to = new DateTime(data.getTime());
            return Days.daysBetween(from,to).getDays();
        } catch (Exception ex){
            ex.printStackTrace();
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("datenum: "+ex.getMessage());
            debug.CloseFile();
            return 0;
        }
    }

    public static long hoursBetween(String format, String fromString, String toString){
        SimpleDateFormat sdf = new SimpleDateFormat(format);
        long hoursBetween = 0;
        try{
            Date data = sdf.parse(fromString);
            DateTime fromTime = new DateTime(data.getTime());
            data = sdf.parse(toString);
            DateTime toTime = new DateTime(data.getTime());
            Period p = new Period(fromTime, toTime);
            long days = p.getDays();
            long hours = p.getHours();
            long minutes = p.getMinutes();
            hoursBetween = (days*24)+hours+minutes;
        } catch (Exception ex) {
            ex.printStackTrace();
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("hoursBetween: "+ex.getMessage());
            debug.CloseFile();
        }
        return  hoursBetween;
    }

    public static double[] min(double[] x, double[] y){//excluding NaNs
        if(x.length!=y.length){
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("Utility.min(): x.length and y.lenght must be equals!");
            debug.CloseFile();
            System.exit(0);
        }
        double[] out = new double[x.length];
        for(int i=0;i<out.length;i++){
            if(Double.isNaN(x[i])){
                if(Double.isNaN(y[i]))
                    out[i]=Double.NaN;
                else
                    out[i]= y[i];
            } else {
                if(Double.isNaN(y[i]))
                    out[i]=x[i];
                else{
                    out[i] = (x[i]<y[i]) ? x[i] : y[i];
                }
            }
        }
        return out;
    }
    public static double[] mean(double[] x, double[] y){
        double[] out = new double[x.length];
        if(x.length!=y.length){
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("Utility.mean(): x.length and y.lenght must be equals!");
            debug.CloseFile();
            System.exit(0);
        }
        for(int i=0;i<x.length;i++)
            out[i]=(x[i]+y[i])/2;
        return out;
    }


    public static long datenum(long year, long month, long day, long hour, long min){
        SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMddHHmm");
        try {
            Date data = sdf.parse("000001000000");
            DateTime from = new DateTime(data.getTime());

            String stringYear = ""+year;
            if((year-100)<0){//check year format (2015 or 15)
                stringYear="20"+year;
            }
            String stringMonth = ""+month;
            if(month<10){
                stringMonth = "0"+month;
            }
            String stringday = ""+day;
            if(day<10){
                stringday = "0"+day;
            }
            String stringhour=""+hour;
            if(hour<10){
                stringhour = "0"+hour;
            }
            String stringMin=""+min;
            if(min<10){
                stringMin="0"+min;
            }

            String toString = ""+year+""+stringMonth+""+stringday+""+stringhour+""+stringMin;
            data = sdf.parse(toString);
            DateTime to = new DateTime(data.getTime());
            return Days.daysBetween(from,to).getDays();
        } catch (Exception ex){
            ex.printStackTrace();
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("datenum: "+ex.getMessage());
            debug.CloseFile();
            return 0;
        }
    }


    public static String datestr(long datenum, long hours, long minutes){
        SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMddHHmm");
        try {
            Date data = sdf.parse("000000000000");
            Calendar c = Calendar.getInstance();
            c.setTime(data);
            c.add(Calendar.DATE, (int) datenum);
            SimpleDateFormat outDate = new SimpleDateFormat("dd-MM-yyyy");
            String out = outDate.format(c.getTime());
            if(hours<10)
                out+=" 0"+hours;
            else
                out+=" "+hours;
            if(minutes<10)
                out+=":0"+minutes;
            else
                out+=":"+minutes;
            return out;
        } catch (Exception ex) {
            ex.printStackTrace();
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("datestr: "+ex.getMessage());
            debug.CloseFile();
            return "";
        }
    }

    public static int fix(double a){
        return (int) a;
    }


    public static int[] deepCopy(int[] in, int until){
        int[] out = new int[in.length+until];
        for(int i=0;i<out.length;i++){
            out[i]=in[i];
        }
        return out;
    }

    public static double[] deepCopy(double[] in){
        double[] out = new double[in.length];
        for(int i=0;i<in.length;i++){
            out[i]=in[i];
        }
        return  out;
    }

    public static double[][] deepCopy(double[][] in){
        double[][] out = new double[in.length][in[0].length];
        for(int i=0;i<in.length;i++){
            for(int j=0;j<in[0].length;j++){
                out[i][j] = in[i][j];
            }
        }
        return out;
    }

    public static double[][][] deepCopy(double[][][] in){
        double[][][] out = new double[in.length][in[0].length][in[0][0].length];
        for(int k=0;k<in.length;k++){
            for(int i=0;i<in[0].length;i++){
                for(int j=0;j<in[0][0].length;j++){
                    out[k][i][j] = in[k][i][j];
                }
            }
        }
        return out;
    }
}

