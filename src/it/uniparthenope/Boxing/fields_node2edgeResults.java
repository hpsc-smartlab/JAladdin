package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class fields_node2edgeResults {
    double[][] waveHeight_edges;
    double[][] wavePeriod_edges;
    double[][] Ywave_edges;
    double[][] Xwave_edges;
    double[][] windMAGN_edges;
    double[][] Ywind_edges;
    double[][] Xwind_edges;
    double[][] waveLenght_edges;
    double[] bathy_edges;
    double[] J_edges;

    public fields_node2edgeResults(){
        waveHeight_edges = null;
        wavePeriod_edges = null;
        Ywave_edges = null;
        Xwave_edges = null;
        windMAGN_edges = null;
        Ywind_edges = null;
        Xwind_edges = null;
        waveLenght_edges = null;
        bathy_edges = null;
    }

    public void setValues(ArrayList<double[][]> mDynamicF_EWeightsResults){
        if(mDynamicF_EWeightsResults.size()==7){//wave and wind fields
            waveHeight_edges = mDynamicF_EWeightsResults.get(0);
            wavePeriod_edges = mDynamicF_EWeightsResults.get(1);
            Ywave_edges = mDynamicF_EWeightsResults.get(2);
            Xwave_edges = mDynamicF_EWeightsResults.get(3);
            windMAGN_edges = mDynamicF_EWeightsResults.get(4);
            Ywind_edges = mDynamicF_EWeightsResults.get(5);
            Xwind_edges = mDynamicF_EWeightsResults.get(6);
        } else { //just wave fields
            waveHeight_edges = mDynamicF_EWeightsResults.get(0);
            Ywave_edges = mDynamicF_EWeightsResults.get(1);
            Xwave_edges = mDynamicF_EWeightsResults.get(2);
        }
    }

    public double[][] getWaveHeight_edges() {
        return waveHeight_edges;
    }

    public void setWaveHeight_edges(double[][] waveHeight_edges) {
        this.waveHeight_edges = waveHeight_edges;
    }

    public double[][] getWavePeriod_edges() {
        return wavePeriod_edges;
    }

    public void setWavePeriod_edges(double[][] wavePeriod_edges) {
        this.wavePeriod_edges = wavePeriod_edges;
    }

    public double[][] getYwave_edges() {
        return Ywave_edges;
    }

    public void setYwave_edges(double[][] ywave_edges) {
        Ywave_edges = ywave_edges;
    }

    public double[][] getXwave_edges() {
        return Xwave_edges;
    }

    public void setXwave_edges(double[][] xwave_edges) {
        Xwave_edges = xwave_edges;
    }

    public double[][] getWindMAGN_edges() {
        return windMAGN_edges;
    }

    public void setWindMAGN_edges(double[][] windMAGN_edges) {
        this.windMAGN_edges = windMAGN_edges;
    }

    public double[][] getYwind_edges() {
        return Ywind_edges;
    }

    public void setYwind_edges(double[][] ywind_edges) {
        Ywind_edges = ywind_edges;
    }

    public double[][] getXwind_edges() {
        return Xwind_edges;
    }

    public void setXwind_edges(double[][] xwind_edges) {
        Xwind_edges = xwind_edges;
    }

    public double[][] getWaveLenght_edges() {
        return waveLenght_edges;
    }

    public void setWaveLenght_edges(double[][] waveLenght_edges) {
        this.waveLenght_edges = waveLenght_edges;
    }

    public double[] getBathy_edges() {
        return bathy_edges;
    }

    public void setBathy_edges(double[] bathy_edges) {
        this.bathy_edges = bathy_edges;
    }

    public double[] getJ_edges() {
        return J_edges;
    }

    public void setJ_edges(double[] j_edges) {
        J_edges = j_edges;
    }
}
