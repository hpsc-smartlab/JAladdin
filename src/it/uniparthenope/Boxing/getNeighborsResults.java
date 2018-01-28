package it.uniparthenope.Boxing;

public class getNeighborsResults {
    private int[] rows;
    private int[] neighbors;

    public getNeighborsResults(int[] rows, int[] neighbors) {
        this.rows = rows;
        this.neighbors = neighbors;
    }

    public int[] getRows() {
        return rows;
    }

    public int[] getNeighbors() {
        return neighbors;
    }
}
