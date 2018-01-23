package it.uniparthenope.Boxing;

import it.uniparthenope.DangerIndexes;
import it.uniparthenope.Parser.JSONManager;
import org.json.simple.parser.ParseException;

import java.io.IOException;
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

    public Edges_definitionResults(boolean flag) throws IOException, ParseException {
        if(flag == true){
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/Edges_Definition.json");
            free_edges = reader.retrieveInt2D("free_edges");
            nogo_edges = reader.retrieveIntegerArrayList("nogo_edges");
            theta_grid = reader.retrieveDoubleArray("theta_grid");
            edge_lenght = reader.retrieveDoubleArray("edge_lenght");
            sh_delay = reader.retrieveDouble2D("sh_delay");
            gear_idx = reader.retrieveInt2D("gear_idx");
            safe_indexes = reader.retrieveInt3D("safe_indexes");

            tdep_danger_idx = new ArrayList<>();
            int n_danger_idxs = reader.retrieveInteger("DangerIndexesNumber");
            for(int i=0;i<n_danger_idxs;++i)
                tdep_danger_idx.add(new DangerIndexes(i));

            waveHeight_edges = reader.retrieveDouble2D("waveHeight_edges");
            wavePeriod_edges = reader.retrieveDouble2D("wavePeriod_edges");
            waveLength_edges = reader.retrieveDouble2D("waveLength_edges");
            waveDir_edges = reader.retrieveDouble2D("waveDir_edges");
            windMAGN_edges = reader.retrieveDouble2D("windMAGN_edges");
            windDIR_edges = reader.retrieveDouble2D("windDIR_edges");
            bathy_edges = reader.retrieveDoubleArray("bathy_edges");
            I_bool = reader.retrieveBooleanArray("I_bool");
            I_ord = reader.retrieveIntArray("I_ord");
            I_point = reader.retrieveIntArray("I_point");

            reader.dispose();
        }
    }

    public void saveState() throws IOException{
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/Edges_definition.json");
        writer.putInt2D("free_edges", free_edges);
        writer.putIntegerArrayList("nogo_edges", nogo_edges);
        writer.putDoubleArray("theta_grid", theta_grid);
        writer.putDoubleArray("edge_lenght", edge_lenght);
        writer.putDouble2D("sh_delay", sh_delay);
        writer.putInt2D("gear_idx", gear_idx);
        writer.putInt3D("safe_indexes", safe_indexes);
        writer.putInteger("DangerIndexesNumber",tdep_danger_idx.size());
        for(int i=0;i<tdep_danger_idx.size();++i)
            tdep_danger_idx.get(i).saveState(i);
        writer.putDouble2D("waveHeight_edges", waveHeight_edges);
        writer.putDouble2D("wavePeriod_edges", wavePeriod_edges);
        writer.putDouble2D("waveLength_edges", waveLength_edges);
        writer.putDouble2D("waveDir_edges", waveDir_edges);
        writer.putDouble2D("windMAGN_edges", windMAGN_edges);
        writer.putDouble2D("windDIR_edges", windDIR_edges);
        writer.putDoubleArray("bathy_edges", bathy_edges);
        writer.putBooleanArray("I_bool", I_bool);
        writer.putIntArray("I_ord", I_ord);
        writer.putIntArray("I_point", I_point);
        writer.dispose();
    }


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



    public int[][] getFree_edges() {
        return free_edges;
    }

    public ArrayList<Integer> getNogo_edges() {
        return nogo_edges;
    }

    public double[] getTheta_grid() {
        return theta_grid;
    }

    public double[] getEdge_lenght() {
        return edge_lenght;
    }

    public double[][] getSh_delay() {
        return sh_delay;
    }

    public int[][] getGear_idx() {
        return gear_idx;
    }

    public int[][][] getSafe_indexes() {
        return safe_indexes;
    }

    public ArrayList<DangerIndexes> getTdep_danger_idx() {
        return tdep_danger_idx;
    }

    public double[][] getWaveHeight_edges() {
        return waveHeight_edges;
    }

    public double[][] getWavePeriod_edges() {
        return wavePeriod_edges;
    }

    public double[][] getWaveLength_edges() {
        return waveLength_edges;
    }

    public double[][] getWaveDir_edges() {
        return waveDir_edges;
    }

    public double[][] getWindMAGN_edges() {
        return windMAGN_edges;
    }

    public double[][] getWindDIR_edges() {
        return windDIR_edges;
    }

    public double[] getBathy_edges() {
        return bathy_edges;
    }

    public boolean[] getI_bool() {
        return I_bool;
    }

    public int[] getI_ord() {
        return I_ord;
    }

    public int[] getI_point() {
        return I_point;
    }
}
