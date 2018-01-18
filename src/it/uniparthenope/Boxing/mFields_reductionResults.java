package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class mFields_reductionResults {
    private double[] lat_red;
    private double[] lon_red;
    private ArrayList<double[][][]> outs = new ArrayList<>();


    public double[] getLat_red() {
        return lat_red;
    }

    public void setLat_red(double[] lat_red) {
        this.lat_red = lat_red;
    }

    public double[] getLon_red() {
        return lon_red;
    }

    public void setLon_red(double[] lon_red) {
        this.lon_red = lon_red;
    }

//    public double[][][] getOut1() {
//        return out1;
//    }
//
//    public void setOut1(double[][][] out1) {
//        this.out1 = out1;
//    }
//
//    public double[][][] getOut2() {
//        return out2;
//    }
//
//    public void setOut2(double[][][] out2) {
//        this.out2 = out2;
//    }
//
//    public double[][][] getOut3() {
//        return out3;
//    }
//
//    public void setOut3(double[][][] out3) {
//        this.out3 = out3;
//    }

    public void addOut(double[][][] out){ this.outs.add(out); }

    public double[][][] getOut(int index){
        return this.outs.get(index);
    }

    public void setOut(int index, double[][][] mat){
        this.outs.set(index, mat);
    }

    public void multiplyElementFor(int index, double factor){
        for(int k=0;k<this.outs.get(index).length;++k){
            for(int i=0;i<this.outs.get(index)[0].length;++i){
                for(int j=0;j<this.outs.get(index)[0][0].length;++j){
                    this.outs.get(index)[k][i][j] = this.outs.get(index)[k][i][j] * factor;
                }
            }
        }
    }
}
