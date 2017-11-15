package it.uniparthenope;

public class VisirModel {
    private Const constants;
    private Optim optim;
    private Ship ship;
    private SpatialGrid sGrid;
    private TemporalGrid tGrid;

    public VisirModel(){//Initialize with standard values
        this.constants = new Const();
        this.optim = new Optim();
        this.ship = new Ship();
        this.sGrid = new SpatialGrid();
        this.tGrid = new TemporalGrid();
    }

    public Const getConstants() {
        return this.constants;
    }

    public Optim getOptim() {
        return this.optim;
    }

    public Ship getShip() {
        return this.ship;
    }

    public SpatialGrid getsGrid() {
        return this.sGrid;
    }

    public TemporalGrid gettGrid() {
        return this.tGrid;
    }
}
