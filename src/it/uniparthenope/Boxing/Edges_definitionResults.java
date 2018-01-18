package it.uniparthenope.Boxing;

import it.uniparthenope.DangerIndexes;

import java.util.ArrayList;

public class Edges_definitionResults {
    int[][] free_edges;
    ArrayList<Integer> nogo_edges;
    double[] theta_grid;
    double[] edge_lenght;
    double[][] sh_delay;
    int[][] gear_idx;
    int[][][] safe_indexes;
    ArrayList<DangerIndexes> tdep_danger_idx;
    double[][] waveHeight_edges;
    double[][] wavePeriod_edges;
    double[][] waveLength_edges;
    double[][] waveDir_edges;
    double[][] windMAGN_edges;
    double[][] windDIR_edges;
    double[] bathy_edges;
    boolean[] I_bool;
    int[] I_ord;
    int[] I_point;

    public Edges_definitionResults(int[][] free_edges, ArrayList<Integer> nogo_edges, double[] theta_grid, double[] edge_lenght,
                                   double[][] sh_delay, int[][] gear_idx, int[][][] safe_indexes,
                                   ArrayList<DangerIndexes> tdep_danger_idx, double[][] waveHeight_edges,
                                   double[][] wavePeriod_edges, double[][] waveLength_edges, double[][] waveDir_edges,
                                   double[][] windMAGN_edges, double[][] windDIR_edges, double[] bathy_edges,
                                   boolean[] i_bool, int[] i_ord, int[] i_point) {
        this.free_edges = free_edges;
        this.nogo_edges = nogo_edges;
        this.theta_grid = theta_grid;
        this.edge_lenght = edge_lenght;
        this.sh_delay = sh_delay;
        this.gear_idx = gear_idx;
        this.safe_indexes = safe_indexes;
        this.tdep_danger_idx = tdep_danger_idx;
        this.waveHeight_edges = waveHeight_edges;
        this.wavePeriod_edges = wavePeriod_edges;
        this.waveLength_edges = waveLength_edges;
        this.waveDir_edges = waveDir_edges;
        this.windMAGN_edges = windMAGN_edges;
        this.windDIR_edges = windDIR_edges;
        this.bathy_edges = bathy_edges;
        I_bool = i_bool;
        I_ord = i_ord;
        I_point = i_point;
    }
}
