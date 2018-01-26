package it.uniparthenope;

import java.util.stream.DoubleStream;

public class MemorySaver {

    public static double[][][] getXInset(double[][][] A, double deg2rad){
        double[][][] B = new double[A.length][A[0].length][A[0][0].length];
        for(int k=0; k<B.length; ++k){
            for(int i=0;i<B[0].length;++i){
                for(int j=0;j<B[0][0].length;++j){
                    B[k][i][j] = Math.sin(deg2rad * A[k][i][j]);
                }
            }
        }
        return B;
    }

    public static double[][][] getYInset(double[][][] A, double deg2rad){
        double[][][] B = new double[A.length][A[0].length][A[0][0].length];
        for(int k=0; k<B.length; ++k){
            for(int i=0;i<B[0].length;++i){
                for(int j=0;j<B[0][0].length;++j){
                    B[k][i][j] = Math.cos(deg2rad * A[k][i][j]);
                }
            }
        }
        return B;
    }

    public static int getIdx(int[][] free_edges, int valueA, int valueB){
        boolean found = false;
        boolean overIndex = false;
        int i=0;
        while(!found && !overIndex){
            if( (free_edges[i][0] == valueA )&& (free_edges[i][1] == valueB) ){
                found = true;
            } else{
                ++i;
                if(i==free_edges.length){
                    overIndex = true;
                }
            }

        }
        if(found) return i;
        else return -1;
    }

    public static double[] addNaNAtTheEnd(double[] A){
        double[] B = new double[A.length+1];
        for(int i=0;i<A.length;++i)
            B[i] = A[i];
        B[B.length-1] = Double.NaN;
        return B;
    }
    
}
