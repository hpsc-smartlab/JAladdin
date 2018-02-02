package it.uniparthenope;

import it.uniparthenope.Boxing.Point;

import java.util.ArrayList;

public class Polygon {
    private ArrayList<Point> corners;

    public Polygon(){
        this.corners = new ArrayList<>();
    }

    public void addCorner(double lon, double lat){
        this.corners.add(new Point(lon, lat));
    }

    public int getCornersNumber(){
        return this.corners.size();
    }

    public Point getCorner(int i){
        return this.corners.get(i);
    }
}

