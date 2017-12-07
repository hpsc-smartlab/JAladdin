package it.uniparthenope.Boxing;

public class findResults {
    int size;
    double[] elements;
    int[] indexes;

    public findResults(int size, double[] elements, int[] indexes){
        this.size = size;
        this.elements = elements;
        this.indexes = indexes;
    }

    public int getSize() {
        return size;
    }

    public double[] getElements() {
        return elements;
    }

    public int[] getIndexes() {
        return indexes;
    }
}
