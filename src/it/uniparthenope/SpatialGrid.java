package it.uniparthenope;

public class SpatialGrid {
    private int numberOfNeighbors;
    private int numberOfShortestPaths;
    private double nodesLargeN;
    private int invStepFields;
    private int minNoGridPoints;

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

    public int getInvStepFields() {
        return invStepFields;
    }

    public int getMinNoGridPoints() {
        return minNoGridPoints;
    }

    public int getNumberOfNeighbors() {
        return numberOfNeighbors;
    }

    public int getNumberOfShortestPaths() {
        return numberOfShortestPaths;
    }

    public void setInvStepFields(int invStepFields) {
        this.invStepFields = invStepFields;
    }

    public void setMinNoGridPoints(int minNoGridPoints) {
        this.minNoGridPoints = minNoGridPoints;
    }

    public void setNodesLargeN(double nodesLargeN) {
        this.nodesLargeN = nodesLargeN;
    }

    public void setNumberOfNeighbors(int numberOfNeighbors) {
        this.numberOfNeighbors = numberOfNeighbors;
    }

    public void setNumberOfShortestPaths(int numberOfShortestPaths) {
        this.numberOfShortestPaths = numberOfShortestPaths;
    }
}
