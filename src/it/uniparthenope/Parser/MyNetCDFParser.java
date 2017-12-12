package it.uniparthenope.Parser;

import it.uniparthenope.Boxing.waveForecastResults;
import it.uniparthenope.Debug.MyFileWriter;
import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayInt;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;
import java.util.ArrayList;

import it.uniparthenope.Boxing.parseMedOneMinResults;

public class MyNetCDFParser {
    private String filePath;
    private NetcdfFile dataFile = null;
    private boolean fileExists;

    public MyNetCDFParser(String filePath){
        this.filePath = filePath;
        try {
            this.dataFile = NetcdfFile.open(filePath, null);
            this.fileExists = true;
        } catch (Exception e){
            e.printStackTrace();
            this.fileExists = false;
        }
    }

    public boolean isFileExists() {
        return fileExists;
    }

    public parseMedOneMinResults parseMedOneMin(){
        //Retrive the variable named "latitude"
        Variable latitude = dataFile.findVariable("latitude");
        if(latitude == null){
            System.out.println("Can't find the variable named latitude");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedOneMinResults: Can't find the variable named latitude");
            debug.CloseFile();
            System.exit(0);
        }
        int[] latitudeShape = latitude.getShape();
        int[] latitudeOrigin = new int[1];
        //Retrive the variable named "longitude"
        Variable longitude = dataFile.findVariable("longitude");
        if(longitude == null){
            System.out.println("Can't find the variable named longitude");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedOneMinResults: Can't find the variable named longitude");
            debug.CloseFile();
            System.exit(0);
        }
        int[] longitudeShape = longitude.getShape();
        int[] longitudeOrigin = new int[1];
        //Retrive the variable named "z" (depth)
        Variable z = dataFile.findVariable("z");
        if(z == null){
            System.out.println("Can't find the variable named z");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedOneMinResults: Can't find the variable named z");
            debug.CloseFile();
            System.exit(0);
        }
        int[] zShape = z.getShape();
        int[] zOrigin = new int[2];

        try{
            //Loading data from NetCDF file
            ArrayDouble.D1 latitudeArray = (ArrayDouble.D1) latitude.read(latitudeOrigin,latitudeShape);
            ArrayDouble.D1 longitudeArray = (ArrayDouble.D1) longitude.read(longitudeOrigin, longitudeShape);
            ArrayDouble.D2 zArray = (ArrayDouble.D2) z.read(zOrigin, zShape);

            //Passing from ArrayDouble to ArrayList<Double> and Double[][]
            ArrayList<Double> lat = new ArrayList<>();
            ArrayList<Double> lon = new ArrayList<>();
            double[][] depth = new double[zShape[0]][zShape[1]];

            for(int i=0; i<latitudeShape[0]; i++){
                lat.add(latitudeArray.get(i));
            }
            for(int i=0; i<longitudeShape[0];i++){
                lon.add(longitudeArray.get(i));
            }
            for(int i=0;i<zShape[0];i++){
                for(int j=0;j<zShape[1];j++){
                    depth[i][j] = zArray.get(i, j);
                }
            }

            return new parseMedOneMinResults(lat, lon, depth);

        } catch(Exception e){
            e.printStackTrace();
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedOneMinResults: "+e.getMessage());
            debug.CloseFile();
            System.exit(0);
        } finally {
            try {
                this.dataFile.close();
            } catch (Exception e){
                e.printStackTrace();
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("parseMedOneMinResults: "+e.getMessage());
                debug.CloseFile();
                System.exit(0);

            }
        }
        return null;
    }

    public waveForecastResults parseWaveForecastData(){
        //Retrive the variable named "time"
        Variable time = dataFile.findVariable("time");
        if(time==null){
            System.out.println("Can't find the variable named time");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedWaveForecastData: Can't find the variable named time");
            debug.CloseFile();
            System.exit(0);
        }
        int[] timeShape = time.getShape();
        int[] timeOrigin = new int[1];
        //Retrive the variable named "latitude"
        Variable latitude = dataFile.findVariable("latitude");
        if(latitude == null){
            System.out.println("Can't find the variable named latitude");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedWaveForecastData: Can't find the variable named latitude");
            debug.CloseFile();
            System.exit(0);
        }
        int[] latitudeShape = latitude.getShape();
        int[] latitudeOrigin = new int[1];
        //Retrive the variable named "longitude"
        Variable longitude = dataFile.findVariable("longitude");
        if(longitude == null){
            System.out.println("Can't find the variable named longitude");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedWaveForecastData: Can't find the variable named longitude");
            debug.CloseFile();
            System.exit(0);
        }
        int[] longitudeShape = longitude.getShape();
        int[] longitudeOrigin = new int[1];
        //Retrive the variable named "VDIR"
        Variable VDIR = dataFile.findVariable("VDIR");
        if(VDIR==null){
            System.out.println("Can't find the variable named VDIR");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedWaveForecastData: Can't find the variable named VDIR");
            debug.CloseFile();
            System.exit(0);
        }
        int[] VDIRShape = VDIR.getShape();
        int[] VDIROrigin = new int[3];
        //Retrive the variable named "VTDH"
        Variable VTDH = dataFile.findVariable("VTDH");
        if(VTDH==null){
            System.out.println("Can't find the variable named VTDH");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedWaveForecastData: Can't find the variable named VTDH");
            debug.CloseFile();
            System.exit(0);
        }
        int[] VTDHShape = VDIR.getShape();
        int[] VTDHOrigin = new int[3];
        //Retrive the variable named "VTPK"
        Variable VTPK = dataFile.findVariable("VTPK");
        if(VTDH==null){
            System.out.println("Can't find the variable named VTPK");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedWaveForecastData: Can't find the variable named VTPK");
            debug.CloseFile();
            System.exit(0);
        }
        int[] VTPKShape = VDIR.getShape();
        int[] VTPKOrigin = new int[3];

        try {
            //Loading data from NetCDF file
            ArrayInt.D1 timeArray = (ArrayInt.D1) time.read(timeOrigin,timeShape);
            ArrayDouble.D1 latitudeArray = (ArrayDouble.D1) latitude.read(latitudeOrigin,latitudeShape);
            ArrayDouble.D1 longitudeArray = (ArrayDouble.D1) longitude.read(longitudeOrigin, longitudeShape);
            ArrayDouble.D3 vdir3DMatrix = (ArrayDouble.D3) VDIR.read(VDIROrigin, VDIRShape);
            ArrayDouble.D3 vtdh3DMatrix = (ArrayDouble.D3) VTDH.read(VTDHOrigin, VTDHShape);
            ArrayDouble.D3 vtpk3DMatrix = (ArrayDouble.D3) VTPK.read(VTPKOrigin, VTPKShape);
            //passing to parseWaveForecastDataResults data types
            int[] parsedTime = new int[timeShape[0]];
            for(int i=0;i<timeShape[0];i++){
                parsedTime[i] = timeArray.get(i);
            }
            double[] parsedLatitude = new double[latitudeShape[0]];
            for(int i=0;i<latitudeShape[0];i++){
                parsedLatitude[i] = latitudeArray.get(i);
            }
            double[] parsedLongitude = new double[longitudeShape[0]];
            for(int i=0;i<longitudeShape[0];i++){
                parsedLongitude[i] = longitudeArray.get(i);
            }
            double[][][] parsedVDIR = new double[VDIRShape[0]][VDIRShape[1]][VDIRShape[2]];
            for(int i=0;i<VDIRShape[0];i++){
                for(int j=0;j<VDIRShape[1];j++){
                    for(int k=0;k<VDIRShape[2];k++){
                        parsedVDIR[i][j][k] = vdir3DMatrix.get(i, j, k);
                    }
                }
            }
            double[][][] parsedVTDH = new double[VTDHShape[0]][VTDHShape[1]][VTDHShape[2]];
            for(int i=0;i<VTDHShape[0];i++){
                for(int j=0;j<VTDHShape[1];j++){
                    for(int k=0;k<VTDHShape[2];k++){
                        parsedVTDH[i][j][k] = vtdh3DMatrix.get(i, j, k);
                    }
                }
            }
            double[][][] parsedVTPK = new double[VTPKShape[0]][VTPKShape[1]][VTPKShape[2]];
            for(int i=0;i<VTPKShape[0];i++){
                for(int j=0;j<VTPKShape[1];j++){
                    for(int k=0;k<VTPKShape[2];k++){
                        parsedVTPK[i][j][k] = vtpk3DMatrix.get(i, j, k);
                    }
                }
            }
            return new waveForecastResults(parsedVDIR, parsedVTDH, parsedVTPK, parsedTime, parsedLatitude, parsedLongitude);

        } catch (Exception e){
            e.printStackTrace();
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("parseMedWaveForecastData: "+e.getMessage());
            debug.CloseFile();
            System.exit(0);
        } finally {
            try {
                this.dataFile.close();
            } catch (Exception e){
                e.printStackTrace();
                MyFileWriter debug = new MyFileWriter("","debug",false);
                debug.WriteLog("parseMedWaveForecastData: "+e.getMessage());
                debug.CloseFile();
                System.exit(0);

            }
        }
        return null;
    }

}
