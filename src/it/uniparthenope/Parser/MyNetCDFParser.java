package it.uniparthenope.Parser;

import it.uniparthenope.Boxing.waveForecastResults;
import it.uniparthenope.Debug.MyFileWriter;
import it.uniparthenope.Utility;
import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayFloat;
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

    public parseMedOneMinResults parseMedOneMin(String outdir){//Bathy
        //Retrive the variable named "latitude"
        Variable latitude = dataFile.findVariable("latitude");
        if(latitude == null){
            System.out.println("Can't find the variable named latitude");
            MyFileWriter debug = new MyFileWriter("","debug",false, outdir);
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
            MyFileWriter debug = new MyFileWriter("","debug",false, outdir);
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
            MyFileWriter debug = new MyFileWriter("","debug",false, outdir);
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
            errLog(e, outdir);
        } finally {
            try {
                this.dataFile.close();
            } catch (Exception e){
                errLog(e, outdir);

            }
        }
        return null;
    }

    public waveForecastResults parseWaveForecastData(String outdir){
        boolean flag = false;
        //Retrive the variable named "time"
        Variable time = dataFile.findVariable("time");
        if(time==null){
            System.out.println("Can't find the variable named time");
            MyFileWriter debug = new MyFileWriter("","debug",false, outdir);
            debug.WriteLog("parseMedWaveForecastData: Can't find the variable named time");
            debug.CloseFile();
            System.exit(0);
        }
        int[] timeShape = time.getShape();
        int[] timeOrigin = new int[1];
        //Retrive the variable named "latitude"
        Variable latitude = dataFile.findVariable("lat");
        if(latitude == null){
            System.out.println("Can't find the variable named lat");
            MyFileWriter debug = new MyFileWriter("","debug",false, outdir);
            debug.WriteLog("parseMedWaveForecastData: Can't find the variable named lat");
            debug.CloseFile();
            System.exit(0);
        }
        int[] latitudeShape = latitude.getShape();
        int[] latitudeOrigin = new int[1];
        //Retrive the variable named "longitude"
        Variable longitude = dataFile.findVariable("lon");
        if(longitude == null){
            System.out.println("Can't find the variable named lon");
            MyFileWriter debug = new MyFileWriter("","debug",false, outdir);
            debug.WriteLog("parseMedWaveForecastData: Can't find the variable named lon");
            debug.CloseFile();
            System.exit(0);
        }
        int[] longitudeShape = longitude.getShape();
        int[] longitudeOrigin = new int[1];
        //Retrive the variable named "VDIR"
        Variable VDIR = dataFile.findVariable("VDIR");
        if(VDIR==null){
            flag=true;
            VDIR = dataFile.findVariable("dirmn");
            if(VDIR==null){
                System.out.println("Can't find the variable named VDIR or dirmn (wave direction)");
                MyFileWriter debug = new MyFileWriter("","debug",false, outdir);
                debug.WriteLog("parseMedWaveForecastData: Can't find the variable named VDIR");
                debug.CloseFile();
                System.exit(0);
            }
        }
        int[] VDIRShape = VDIR.getShape();
        int[] VDIROrigin = new int[3];
        //Retrive the variable named "VTDH"
        Variable VTDH = dataFile.findVariable("VTDH");
        if(VTDH==null){
            VTDH = dataFile.findVariable("hs");
            if(VTDH == null){
                System.out.println("Can't find the variable named VTDH or hs (wave height)");
                MyFileWriter debug = new MyFileWriter("","debug",false, outdir);
                debug.WriteLog("parseMedWaveForecastData: Can't find the variable named VTDH");
                debug.CloseFile();
                System.exit(0);
            }
        }
        int[] VTDHShape = VDIR.getShape();
        int[] VTDHOrigin = new int[3];
        //Retrive the variable named "VTPK"
        Variable VTPK = dataFile.findVariable("VTPK");
        if(VTPK==null){
            VTPK = dataFile.findVariable("peakp");
            if(VTPK==null) {
                System.out.println("Can't find the variable named VTPK or tmn (wave period)");
                MyFileWriter debug = new MyFileWriter("", "debug", false, outdir);
                debug.WriteLog("parseMedWaveForecastData: Can't find the variable named VTPK");
                debug.CloseFile();
                System.exit(0);
            }
        }
        int[] VTPKShape = VDIR.getShape();
        int[] VTPKOrigin = new int[3];
        try {
            //Loading data from NetCDF file and passing them to parseWaveForecastDataResults data types
            int[] parsedTime = new int[timeShape[0]];
            double[] parsedLatitude = new double[latitudeShape[0]];
            double[] parsedLongitude = new double[longitudeShape[0]];
            double[][][] parsedVDIR;
            double[][][] parsedVTDH;
            double[][][] parsedVTPK;

            //date extraction from namefile
            String noExtension = this.filePath.substring(0, this.filePath.length()-3);
            String baseDate = noExtension.substring(noExtension.length()-8);
            int base_year = Integer.parseInt(baseDate.substring(0,4));
            int base_month = Integer.parseInt(baseDate.substring(4,6));
            int base_day = Integer.parseInt(baseDate.substring(6,8));
            int base_hours = 0;

            if(!flag){
                ArrayInt.D1 timeArray = (ArrayInt.D1) time.read(timeOrigin,timeShape);
                ArrayFloat.D1 latitudeArray = (ArrayFloat.D1) latitude.read(latitudeOrigin, latitudeShape);
                ArrayFloat.D1 longitudeArray = (ArrayFloat.D1) longitude.read(longitudeOrigin, longitudeShape);
                ArrayFloat.D3 vdir3DMatrix = (ArrayFloat.D3) VDIR.read(VDIROrigin, VDIRShape);
                ArrayFloat.D3 vtdh3DMatrix = (ArrayFloat.D3) VTDH.read(VTDHOrigin, VTDHShape);
                ArrayFloat.D3 vtpk3DMatrix = (ArrayFloat.D3) VTPK.read(VTPKOrigin, VTPKShape);
                for(int i=0;i<timeShape[0];i++){
                    parsedTime[i] = timeArray.get(i);
                }
                for(int i=0;i<latitudeShape[0];i++){
                    parsedLatitude[i] = latitudeArray.get(i);
                }
                for(int i=0;i<longitudeShape[0];i++){
                    parsedLongitude[i] = longitudeArray.get(i);
                }
                parsedVDIR = ArrayFloatD3ToDouble3D(vdir3DMatrix, VDIRShape);
                parsedVTDH = ArrayFloatD3ToDouble3D(vtdh3DMatrix, VTDHShape);
                parsedVTPK = ArrayFloatD3ToDouble3D(vtpk3DMatrix, VTPKShape);
            } else{
                ArrayDouble.D1 tmpArray = (ArrayDouble.D1) time.read(timeOrigin, timeShape);
                for(int i=0;i<timeShape[0];i++){
                    //ADDED
                    int hoursSinceBaseDate =(int) Math.floor(tmpArray.get(i));
                    int daysToAdd = (int) Math.floor(hoursSinceBaseDate/23);
                    int hoursToAdd = hoursSinceBaseDate%23;
                    int yearsToAdd = 0;

                    long hh = base_hours + hoursToAdd;
                    long dd = base_day + daysToAdd;
                    long MM = base_month;
                    if(switchMonth(dd,base_month)){
                        MM++;
                        dd=1;
                        if(MM>=12){
                            yearsToAdd++;
                            MM=1;
                        }
                    }
                    //parsedTime[i]=(int) Utility.datenum(base_year+yearsToAdd, MM,dd,hh,0,outdir);
                    //getting string
                    String finaldate = ""+(base_year+yearsToAdd);
                    if(MM<10)
                        finaldate+="0"+MM;
                    else
                        finaldate+=MM;
                    if(dd<10)
                        finaldate+="0"+dd;
                    else
                        finaldate+=dd;
                    if(hh<10)
                        finaldate+="0"+hh;
                    else
                        finaldate+=hh;
                    finaldate+="00";

                    //getting time
                    parsedTime[i] = (int) Utility.getTimeStamp(finaldate, outdir);

                }
                tmpArray = (ArrayDouble.D1) latitude.read(latitudeOrigin, latitudeShape);
                for(int i=0;i<latitudeShape[0];i++){
                    parsedLatitude[i] = tmpArray.get(i);
                }
                tmpArray = (ArrayDouble.D1) longitude.read(longitudeOrigin, longitudeShape);
                for(int i=0;i<longitudeShape[0];i++){
                    parsedLongitude[i] = tmpArray.get(i);
                }
                ArrayFloat.D3 vdir3DMatrix = (ArrayFloat.D3) VDIR.read(VDIROrigin, VDIRShape);
                ArrayFloat.D3 vtdh3DMatrix = (ArrayFloat.D3) VTDH.read(VTDHOrigin, VTDHShape);
                ArrayFloat.D3 vtpk3DMatrix = (ArrayFloat.D3) VTPK.read(VTPKOrigin, VTPKShape);
                parsedVDIR = ArrayFloatD3ToDouble3D(vdir3DMatrix, VDIRShape, true);
                parsedVTDH = ArrayFloatD3ToDouble3D(vtdh3DMatrix, VTDHShape);
                parsedVTPK = ArrayFloatD3ToDouble3D(vtpk3DMatrix, VTPKShape);
            }

            return new waveForecastResults(parsedVDIR, parsedVTDH, parsedVTPK, parsedTime, parsedLatitude, parsedLongitude);

        } catch (Exception e){
            errLog(e, outdir);
        } finally {
            try {
                this.dataFile.close();
            } catch (Exception e){
                errLog(e, outdir);
            }
        }
        return null;
    }

    private boolean switchMonth(long day, long month){
        if(month==1 && day>=31)
            return true;
        if(month==2 && day>=28)
            return true;
        if(month==3 && day>=31)
            return true;
        if(month==4 && day>=30)
            return true;
        if(month==5 && day>=31)
            return true;
        if(month==6 && day>=30)
            return true;
        if(month==7 && day>=31)
            return true;
        if(month==8 && day>=31)
            return true;
        if(month==9 && day>=30)
            return true;
        if(month==10 && day>=31)
            return true;
        if(month==11 && day>=30)
            return true;
        if(month==12 && day>=31)
            return true;
        return false;

    }

    private double[][][] ArrayFloatD3ToDouble3D(ArrayFloat.D3 Matrix, int[] shape){
        double[][][] out = new double[shape[2]][shape[0]][shape[1]];
        for(int k=0;k<shape[2];k++){
            for(int i=0;i<shape[0];i++){
                for(int j=0;j<shape[1];j++){
                    if((Matrix.get(i,j,k) >= 1.0E20) || Matrix.get(i,j,k) <= -999.9/*-999.9f*/){
                        out[k][i][j] = Double.NaN;
                    } else {
                        out[k][i][j] = Matrix.get(i,j,k);
                    }
                }
            }
        }
        return out;
    }


    private double[][][] ArrayFloatD3ToDouble3D(ArrayFloat.D3 Matrix, int[] shape, boolean flag){
        //Converting radiants to degree and degree to degree north
        double[][][] out = new double[shape[2]][shape[0]][shape[1]];
        double radToDegreeFactor = 180/Math.PI;
        for(int k=0;k<shape[2];k++){
            for(int i=0;i<shape[0];i++){
                for(int j=0;j<shape[1];j++){
                    if((Matrix.get(i,j,k) >= 1.0E20) || Matrix.get(i,j,k) <= -999.9/*-999.9f*/){
                        out[k][i][j] = Double.NaN;
                    } else {
                        out[k][i][j] = Matrix.get(i,j,k)*radToDegreeFactor;
                        //conversion to degree north
                        if((out[k][i][j]>=0.0) && (out[k][i][j]<=90.0)){
                            //1st quarter
                            out[k][i][j] = 90.0-out[k][i][j];
                        } else{
                            if((out[k][i][j]>90.0) && (out[k][i][j]<180.0)){
                                //2-nd quarter
                                out[k][i][j] = 90.0+out[k][i][j];
                            } else{
                                //3-rd and 4-th quarter
                                out[k][i][j] = 360.0-out[k][i][j]+90.0;
                            }
                        }
                    }
                }
            }
        }
        return out;
    }

    private void errLog(Exception e, String outdir){
        e.printStackTrace();
        MyFileWriter debug = new MyFileWriter("","debug",false, outdir);
        debug.WriteLog("parseMedWaveForecastData: "+e.getMessage());
        debug.CloseFile();
        System.exit(0);
    }
}
