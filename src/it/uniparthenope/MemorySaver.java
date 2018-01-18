package it.uniparthenope;

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
    
}
