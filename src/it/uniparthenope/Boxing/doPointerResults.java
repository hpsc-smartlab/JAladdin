package it.uniparthenope.Boxing;

public class doPointerResults {
    boolean[] head_bool;
    int[] head_ord;
    int[] pointer;

    public doPointerResults(boolean[] head_bool, int[] head_ord, int[] pointer) {
        this.head_bool = head_bool;
        this.head_ord = head_ord;
        this.pointer = pointer;
    }

    public boolean[] getHead_bool() {
        return head_bool;
    }

    public int[] getHead_ord() {
        return head_ord;
    }

    public int[] getPointer() {
        return pointer;
    }
}
