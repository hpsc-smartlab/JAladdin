package it.uniparthenope.Parser;

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

    public MyNetCDFParser(String filePath){
        this.filePath = filePath;
        try {
            this.dataFile = NetcdfFile.open(filePath, null);
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    public parseMedOneMinResults parseMedOneMin(){
        //Retrive the variable named "latitude"
        Variable latitude = dataFile.findVariable("latitude");
        if(latitude == null){
            System.out.println("Can't find the variable named latitude");
            return null;
        }
        int[] latitudeShape = latitude.getShape();
        int[] latitudeOrigin = new int[1];
        //Retrive the variable named "longitude"
        Variable longitude = dataFile.findVariable("longitude");
        if(longitude == null){
            System.out.println("Can't find the variable named longitude");
            return null;
        }
        int[] longitudeShape = longitude.getShape();
        int[] longitudeOrigin = new int[1];
        //Retrive the variable named "z" (depth)
        Variable z = dataFile.findVariable("z");
        if(z == null){
            System.out.println("Can't find the variable named z");
            return null;
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
            Double[][] depth = new Double[zShape[0]][zShape[1]];

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
            return null;
        } finally {
            try {
                this.dataFile.close();
            } catch (Exception e){
                e.printStackTrace();
            }
        }

    }

}
