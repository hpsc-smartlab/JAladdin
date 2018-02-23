package it.uniparthenope.EnvFieldsIn;

import it.uniparthenope.Debug.MyFileWriter;
import it.uniparthenope.Utility;
import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayFloat;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;

public class waveCustomNetCDFInputData implements waveInputData {
    private String filename;
    private NetcdfFile dataFile = null;
    private String outDir;

    public waveCustomNetCDFInputData(String outDir){ this.outDir = outDir;}

    public boolean open(String filename) throws IOException {
        this.filename = filename;
        try{
            this.dataFile = NetcdfFile.open(this.filename, null);
            return true;
        } catch (Exception e){
            e.printStackTrace();
            return false;
        }
    }

    @Override
    public double[] getLongitude() {
        //Retrive the variable named "longitude"
        Variable longitude = dataFile.findVariable("lon");
        if(longitude == null){
            System.out.println("Can't find the variable named lon");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveCustomNetCDFInputData: Can't find the variable named lon");
            debug.CloseFile();
            System.exit(0);
        }
        int[] longitudeShape = longitude.getShape();
        int[] longitudeOrigin = new int[1];
        try{
            double[] parsedLongitude = new double[longitudeShape[0]];
            ArrayDouble.D1 longitudeArray = (ArrayDouble.D1) longitude.read(longitudeOrigin, longitudeShape);
            for(int i=0;i<longitudeShape[0];++i)
                parsedLongitude[i] = longitudeArray.get(i);
            return parsedLongitude;
        }catch (Exception ex){
            errLog(ex);
            return new double[0];
        }
    }

    @Override
    public double[] getLatitude() {
        //Retrive the variable named "latitude"
        Variable latitude = dataFile.findVariable("lat");
        if(latitude == null){
            System.out.println("Can't find the variable named lat");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveCustomNetCDFInputData: Can't find the variable named lat");
            debug.CloseFile();
            System.exit(0);
        }
        int[] latitudeShape = latitude.getShape();
        int[] latitudeOrigin = new int[1];
        try {
            double[] parsedLatitude = new double[latitudeShape[0]];
            ArrayDouble.D1 latitudeArray = (ArrayDouble.D1) latitude.read(latitudeOrigin, latitudeShape);
            for(int i=0;i<latitudeShape[0];++i){
                parsedLatitude[i] = latitudeArray.get(i);
            }
            return parsedLatitude;
        } catch (Exception ex){
            errLog(ex);
            return new double[0];
        }
    }

    @Override
    public int[] getTimestep() {
        Variable time = dataFile.findVariable("time");
        if(time==null){
            System.out.println("Can't find the variable named time");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveCustomNetCDFInputData: Can't find the variable named time");
            debug.CloseFile();
            System.exit(0);
        }
        int[] timeShape = time.getShape();
        int[] timeOrigin = new int[1];
        try{
            //Loading data from NetCDF file and passing them to parseWaveForecastDataResults data types
            ArrayDouble.D1 tmpArray = (ArrayDouble.D1) time.read(timeOrigin, timeShape);
            int[] parsedTime = new int[timeShape[0]];
            //date extraction from namefile
            String noExtension = this.filename.substring(0, this.filename.length()-3);
            String baseDate = noExtension.substring(noExtension.length()-8);
            int base_year = Integer.parseInt(baseDate.substring(0,4));
            int base_month = Integer.parseInt(baseDate.substring(4,6));
            int base_day = Integer.parseInt(baseDate.substring(6,8));
            int base_hours = 0;
            for(int i=0;i<timeShape[0];++i){
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
                parsedTime[i] = (int) Utility.getTimeStamp(finaldate, this.outDir);
            }
            return parsedTime;
        } catch (Exception ex){
            errLog(ex);
            return new int[0];
        }
    }


    @Override
    public double[][][] getDirection() {
        Variable VDIR = dataFile.findVariable("dirmn");
        if(VDIR==null){
            System.out.println("Can't find the variable named dirmn (wave direction)");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveCustomNetCDFInputData: Can't find the variable named dirmn");
            debug.CloseFile();
            System.exit(0);
        }
        int[] VDIRShape = VDIR.getShape();
        int[] VDIROrigin = new int[3];
        try {
            ArrayFloat.D3 vdir3DMatrix = (ArrayFloat.D3) VDIR.read(VDIROrigin, VDIRShape);
            return ArrayFloatD3ToDouble3DAndConvert(vdir3DMatrix, VDIRShape);
        } catch (Exception ex){
            errLog(ex);
            return new double[0][][];
        }
    }

    @Override
    public double[][][] getHeight() {
        //Retrive the variable named "VTDH"
        Variable VTDH = dataFile.findVariable("hs");
        if(VTDH == null){
            System.out.println("Can't find the variable named hs (wave height)");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveCustomNetCDFInputData: Can't find the variable named hr");
            debug.CloseFile();
        }
        int[] VTDHShape = VTDH.getShape();
        int[] VTDHOrigin = new int[3];
        try {
            ArrayFloat.D3 vtdh3DMatrix = (ArrayFloat.D3) VTDH.read(VTDHOrigin, VTDHShape);
            return ArrayFloatD3ToDouble3D(vtdh3DMatrix, VTDHShape);
        } catch (Exception ex){
            errLog(ex);
            return new double[0][][];
        }
    }

    @Override
    public double[][][] getPeriod() {
        //Retrive the variable named "VTPK"
        Variable VTPK = dataFile.findVariable("peakp");
        if(VTPK==null){
            System.out.println("Can't find the variable named peakp (wave period)");
            MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
            debug.WriteLog("waveCustomNetCDFInputData: Can't find the variable named peakp");
            debug.CloseFile();
        }
        int[] VTPKShape = VTPK.getShape();
        int[] VTPKOrigin = new int[3];
        try {
            ArrayFloat.D3 vtpk3DMatrix = (ArrayFloat.D3) VTPK.read(VTPKOrigin, VTPKShape);
            return ArrayFloatD3ToDouble3D(vtpk3DMatrix, VTPKShape);
        } catch (Exception ex){
            errLog(ex);
            return new double[0][][];
        }
    }

    @Override
    public void dispatch() {
        if(this.dataFile != null){
            try {
                this.dataFile.close();
            } catch (Exception ex){
                ex.printStackTrace();
            }
        }
    }

    private void errLog(Exception ex){
        MyFileWriter debug = new MyFileWriter("","debug",false, this.outDir);
        debug.WriteLog("waveCustomNetCDFInputData: "+ex.getMessage());
        debug.CloseFile();
        System.exit(0);
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

    private double radToNorthDeg(double rad){
        //conversion factor = 180/Math.PI
        double deg = rad*(180/Math.PI);
        if(deg>90 && deg<360)
            deg=(360-deg+90);
        else{
            if(deg>=0 && deg<=90)
                deg=deg*(-1)+90;
        }
        return deg;
    }

    private double[][][] ArrayFloatD3ToDouble3DAndConvert(ArrayFloat.D3 Matrix, int[] shape){
        //Converting radiants to degree and degree to degree north
        double[][][] out = new double[shape[2]][shape[0]][shape[1]];
        for(int k=0;k<shape[2];k++){
            for(int i=0;i<shape[0];i++){
                for(int j=0;j<shape[1];j++){
                    if(Matrix.get(i,j,k) <= -999.9/*-999.9f*/){
                        out[k][i][j] = Double.NaN;
                    } else {
                        //Converting from radiant to north degree
                        out[k][i][j] = radToNorthDeg(Matrix.get(i, j, k));
                    }
                }
            }
        }
        return out;
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
}
