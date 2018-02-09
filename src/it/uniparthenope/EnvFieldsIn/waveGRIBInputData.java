package it.uniparthenope.EnvFieldsIn;

import it.uniparthenope.Debug.MyFileWriter;
import net.sourceforge.jgrib.GribFile;

import java.io.IOException;

public class waveGRIBInputData implements waveInputData{
    String filename;
    GribFile file = null;
    private String outDir;
    double minLat, minLon, deltaLat, deltaLon;
    int nLat, nLon, offset, timestep;
    double[][][] VTDH = null;
    double[][][] VTPK = null;
    double[][][] VDIR = null;

    public waveGRIBInputData(String outDir){
        this.outDir = outDir;
        org.apache.log4j.BasicConfigurator.configure();
    }

    @Override
    public boolean open(String filename) throws IOException {
        this.filename = filename;
        try{
            this.file = new GribFile(filename);
            this.minLon = file.getLightRecords()[0].getGDS().getGridLon1();
            this.minLat = file.getLightRecords()[0].getGDS().getGridLat1();
            this.deltaLon = file.getLightRecords()[0].getGDS().getGridDX();
            this.deltaLat = file.getLightRecords()[0].getGDS().getGridDY();
            this.nLon = file.getLightRecords()[0].getGDS().getGridNX();
            this.nLat = file.getLightRecords()[0].getGDS().getGridNY();
            this.offset = (int) Math.floor(file.getRecordCount()/2);
            this.timestep = (int) Math.floor(offset/(file.getTypeNames().length/2))-1;//time step for each var
            return true;
        } catch (Exception e){
            e.printStackTrace();
            return false;
        }
    }

    @Override
    public int[] getTimestep() {
        return new int[0];
    }

    @Override
    public double[] getLongitude() {
        //getting longitude array
        double[] longitude = new double[nLon];
        double nextVal = minLon;
        for(int i=0; i<nLon; ++i){
            longitude[i] = nextVal;
            nextVal += deltaLon;
        }
        return longitude;
    }

    @Override
    public double[] getLatitude() {
        //getting latitude array
        double[] latitude = new double[nLat];
        double nextVal = minLat;
        for(int i=0; i<nLat; ++i){
            latitude[i] = nextVal;
            nextVal += deltaLat;
        }
        return latitude;
    }

    @Override
    public double[][][] getDirection() {
        if(this.VDIR == null)
            gettingData();
        return this.VDIR;
    }

    @Override
    public double[][][] getHeight() {
        if(this.VTDH == null)
            gettingData();
        return this.VTDH;
    }

    @Override
    public double[][][] getPeriod() {
        if(this.VTPK == null)
            gettingData();
        return this.VTPK;
    }

    @Override
    public void dispatch() {
        //freeing memory
        this.VDIR = null;
        this.VTPK = null;
        this.VTDH = null;
    }

    private void errLog(Exception ex){
        MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
        debug.WriteLog("waveGRIBInputData: "+ex.getMessage());
        debug.CloseFile();
        System.exit(0);
    }

    private void gettingData(){
        try{
            int wDirIDX = 6;
            int wPeriodIDX = 4;
            int wHeigthIDX = 1;
            int myIDX=0;//every 3 hs
            for(int k=0; k<this.timestep;++k){
                int currentwHeigthIDX = ((wHeigthIDX+(k*7))+offset);
                int currentwPeriodIDX = ((wPeriodIDX+(k*7))+offset);
                int currentwDirIDX = ((wDirIDX+(k*7))+offset);
                for(int i=0; i<nLon; ++i){
                    for(int j=0;j<nLat; ++j){
                        this.VTDH[i][myIDX][j] = file.getRecord(currentwHeigthIDX).getBDS().getValue(i + j * nLon);
                        if(this.VTDH[i][myIDX][j] >= 9.9E24)
                            this.VTDH[i][myIDX][j] = Double.NaN;
                        this.VTDH[i][myIDX+1][j] = this.VTDH[i][myIDX][j];
                        this.VTDH[i][myIDX+2][j] = this.VTDH[i][myIDX][j];

                        this.VTPK[i][myIDX][j] = file.getRecord(currentwPeriodIDX).getBDS().getValue(i + j * nLon);
                        if(this.VTPK[i][myIDX][j] >= 9.9E24)
                            this.VTPK[i][myIDX][j] = Double.NaN;
                        this.VTPK[i][myIDX+1][j] = this.VTPK[i][myIDX][j];
                        this.VTPK[i][myIDX+2][j] = this.VTPK[i][myIDX][j];

                        this.VDIR[i][myIDX][j] = file.getRecord(currentwDirIDX).getBDS().getValue(i + j * nLon);
                        if(this.VDIR[i][myIDX][j] >= 9.9E24)
                            this.VDIR[i][myIDX][j] = Double.NaN;
                        this.VDIR[i][myIDX+1][j] = this.VDIR[i][myIDX][j];
                        this.VDIR[i][myIDX+2][j] = this.VDIR[i][myIDX][j];
                    }
                }
                myIDX=myIDX+3;
            }
        } catch (Exception ex){
            errLog(ex);
            System.exit(0);
        }
    }

}
