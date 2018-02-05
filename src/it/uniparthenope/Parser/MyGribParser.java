package it.uniparthenope.Parser;

import it.uniparthenope.Boxing.waveForecastResults;
import net.sourceforge.jgrib.GribFile;
import net.sourceforge.jgrib.NoValidGribException;
import net.sourceforge.jgrib.NotSupportedException;

import java.io.IOException;

public class MyGribParser {
    public static waveForecastResults parseGrib(String filename, String outDir) throws NotSupportedException, NoValidGribException, IOException {
        org.apache.log4j.BasicConfigurator.configure();
        GribFile file = new GribFile(filename);
        double minLon = file.getLightRecords()[0].getGDS().getGridLon1();
        double minLat = file.getLightRecords()[0].getGDS().getGridLat1();
        double deltaLon = file.getLightRecords()[0].getGDS().getGridDX();
        double deltaLat = file.getLightRecords()[0].getGDS().getGridDY();
        int nLon = file.getLightRecords()[0].getGDS().getGridNX();
        int nLat = file.getLightRecords()[0].getGDS().getGridNY();

        //getting longitude array
        double[] longitude = new double[nLon];
        double nextVal = minLon;
        for(int i=0; i<nLon; ++i){
            longitude[i] = nextVal;
            nextVal += deltaLon;
        }
        //getting latitude array
        double[] latitude = new double[nLat];
        nextVal = minLat;
        for(int i=0; i<nLat; ++i){
            latitude[i] = nextVal;
            nextVal += deltaLat;
        }

        //getting timesteps:
        int offset = (int) Math.floor(file.getRecordCount()/2);
        int timeStep = (int) Math.floor(offset/(file.getTypeNames().length/2))-1;//time step for each var
        //Wave Direction, mean period and wave heigth
        double[][][] waveDir = new double[nLon][timeStep*3][nLat];
        int wDirIDX = 6;
        double[][][] wavePeriod = new double[nLon][timeStep*3][nLat];
        int wPeriodIDX = 4;
        double[][][] waveHeigth = new double[nLon][timeStep*3][nLat];
        int wHeigthIDX = 1;
        int myIDX=0;//every 3 hs
        for(int k=0; k<timeStep;++k){
            int currentwHeigthIDX = ((wHeigthIDX+(k*7))+offset);
            int currentwPeriodIDX = ((wPeriodIDX+(k*7))+offset);
            int currentwDirIDX = ((wDirIDX+(k*7))+offset);
            for(int i=0; i<nLon; ++i){
                for(int j=0;j<nLat; ++j){
                    waveHeigth[i][myIDX][j] = file.getRecord(currentwHeigthIDX).getBDS().getValue(i + j * nLon);
                    if(waveHeigth[i][myIDX][j] >= 9.9E24)
                        waveHeigth[i][myIDX][j] = Double.NaN;
                    waveHeigth[i][myIDX+1][j] = waveHeigth[i][myIDX][j];
                    waveHeigth[i][myIDX+2][j] = waveHeigth[i][myIDX][j];

                    wavePeriod[i][myIDX][j] = file.getRecord(currentwPeriodIDX).getBDS().getValue(i + j * nLon);
                    if(wavePeriod[i][myIDX][j] >= 9.9E24)
                        wavePeriod[i][myIDX][j] = Double.NaN;
                    wavePeriod[i][myIDX+1][j] = wavePeriod[i][myIDX][j];
                    wavePeriod[i][myIDX+2][j] = wavePeriod[i][myIDX][j];

                    waveDir[i][myIDX][j] = file.getRecord(currentwDirIDX).getBDS().getValue(i + j * nLon);
                    if(waveDir[i][myIDX][j] >= 9.9E24)
                        waveDir[i][myIDX][j] = Double.NaN;
                    waveDir[i][myIDX+1][j] = waveDir[i][myIDX][j];
                    waveDir[i][myIDX+2][j] = waveDir[i][myIDX][j];
//                    waveHeigth[i][k][j] = file.getRecord(currentwHeigthIDX).getBDS().getValue(i + j * nLon);
//                    wavePeriod[i][k][j] = file.getRecord(currentwPeriodIDX).getBDS().getValue(i + j * nLon);
//                    waveDir[i][k][j] = file.getRecord(currentwDirIDX).getBDS().getValue(i + j * nLon);
//                    if(waveHeigth[i][k][j] >= 9.9E24)
//                        waveHeigth[i][k][j] = Double.NaN;
//                    if(wavePeriod[i][k][j] >= 9.9E24)
//                        waveHeigth[i][k][j] = Double.NaN;
//                    if(waveDir[i][k][j] >= 9.9E24)
//                        waveHeigth[i][k][j] = Double.NaN;
                }
            }
            myIDX=myIDX+3;
        }
       return new waveForecastResults(waveDir, waveHeigth, wavePeriod, null, latitude, longitude);
    }
}
