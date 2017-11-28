package it.uniparthenope.Boxing;

public class idx_ref2inset_gridResults {
    long[] row;
    long[] col;
    long[] idx;

    public idx_ref2inset_gridResults(long[] row, long[] col, long[] idx){
        this.row = row;
        this.col = col;
        this.idx = idx;
    }

    public long[] getRow() {
        return row;
    }

    public long[] getCol() {
        return col;
    }

    public long[] getIdx() {
        return idx;
    }
}
