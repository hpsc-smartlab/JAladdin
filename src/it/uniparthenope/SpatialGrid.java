package it.uniparthenope;

public class SpatialGrid {
    private long numberOfNeighbors;
    private long numberOfShortestPaths;
    private double nodesLargeN;
    private long invStepFields;
    private long minNoGridPoints;

    public SpatialGrid(){
        this.numberOfNeighbors = 24;
        this.numberOfShortestPaths = 1;
        this.nodesLargeN = 1.0 * Math.pow(10,5);//Limited by ram availability
        this.invStepFields = 4;//inverse grid size [deg^{-1}] of original wind forecast fields
        this.minNoGridPoints = 4; //minimum number of grid points along either direction (meridional or zonal) - TWA inconsistencies may result from tto small such a value
    }

    public double getNodesLargeN() {
        return nodesLargeN;
    }

    public long getInvStepFields() {
        return invStepFields;
    }

    public long getMinNoGridPoints() {
        return minNoGridPoints;
    }

    public long getNumberOfNeighbors() {
        return numberOfNeighbors;
    }

    public long getNumberOfShortestPaths() {
        return numberOfShortestPaths;
    }

}
