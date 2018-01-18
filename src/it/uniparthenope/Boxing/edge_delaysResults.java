package it.uniparthenope.Boxing;

import it.uniparthenope.DangerIndexes;

import java.util.ArrayList;

public class edge_delaysResults {
    private double[][] sh_delay;
    private int[][][] safe_indexes;
    private ArrayList<DangerIndexes> tdep_danger_idx;
    private ArrayList<int[][]> varargout;

    public edge_delaysResults(){
        this.varargout = new ArrayList<>();
    }

    public edge_delaysResults(double[][] sh_delay, int[][][] safe_indexes, ArrayList<DangerIndexes> tdep_danger_idx) {
        this.sh_delay = sh_delay;
        this.safe_indexes = safe_indexes;
        this.tdep_danger_idx = tdep_danger_idx;
        this.varargout = new ArrayList<>();
    }

    public double[][] getSh_delay() {
        return sh_delay;
    }

    public void setSh_delay(double[][] sh_delay) {
        this.sh_delay = sh_delay;
    }

    public int[][][] getSafe_indexes() {
        return safe_indexes;
    }

    public void setSafe_indexes(int[][][] safe_indexes) {
        this.safe_indexes = safe_indexes;
    }

    public ArrayList<DangerIndexes> getTdep_danger_idx() {
        return tdep_danger_idx;
    }

    public void setTdep_danger_idx(ArrayList<DangerIndexes> tdep_danger_idx) {
        this.tdep_danger_idx = tdep_danger_idx;
    }

    public ArrayList<int[][]> getVarargout() {
        return varargout;
    }

    public void AddVarargout(int[][] arg){
        varargout.add(arg);
    }
}
