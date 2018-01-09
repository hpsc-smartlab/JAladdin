package it.uniparthenope.Boxing;

public class idx_ref2inset_gridV2Results {
    int[][] idx;
    int newDim;

    public idx_ref2inset_gridV2Results(int[][] idx, int newDim) {
        this.idx = idx;
        this.newDim = newDim;//n rows without -1 elements in the column
    }

    public int[][] getIdx() {
        return idx;
    }

    public int getNewDim() {
        return newDim;
    }
}
