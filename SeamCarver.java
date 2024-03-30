
/******************************************************************************
 *  Compilation:  javac SeamCarver.java
 *  Execution:    java SeamCarver input.png columnsToRemove rowsToRemove
 *
 *                
 *
 *  Read image from file specified as command line argument. Use SeamCarver
 *  to remove number of rows and columns specified as command line arguments.
 *  Show the images and print time elapsed to screen.
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.Stopwatch;

public class SeamCarver {
    private double[][] energy;

    private int[][] rgb;
    private int height;
    private int width;
    private boolean tpose;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null) {
            throw new IllegalArgumentException("picture is null");
        }
        this.height = picture.height();
        this.width = picture.width();
        this.energy = new double[width][height];
        this.rgb = new int[width][height];
        findEnergy(picture);
    }

    // finding energy of all pixles
    private void findEnergy(Picture picture) {
        for (int col = 0; col < width; col++) {
            for (int row = 0; row < height; row++) {
                if (row == 0 || row == height - 1 || col == 0 || col == width - 1) {
                    energy[col][row] = 1000.00;
                    rgb[col][row] = picture.getRGB(col, row);
                    continue;
                }
                int left = picture.getRGB(col - 1, row);
                int right = picture.getRGB(col + 1, row);
                int up = picture.getRGB(col, row - 1);
                int down = picture.getRGB(col, row + 1);
                rgb[col][row] = picture.getRGB(col, row);
                energy[col][row] = Math.sqrt(energyFunction(left, right) + energyFunction(up, down));
            }
        }

    }

    // finding the energy gradient
    private double energyFunction(int first, int second) {
        int r1 = (first >> 16) & 0xFF, g1 = (first >> 8) & 0xFF, b1 = (first >> 0) & 0xFF;
        int r2 = (second >> 16) & 0xFF, g2 = (second >> 8) & 0xFF, b2 = (second >> 0) & 0xFF;
        return Math.pow((r1 - r2), 2) + Math.pow((g1 - g2), 2) + Math.pow((b1 - b2), 2);
    }

    // current picture
    public Picture picture() {
        if (tpose) {
            transposeEnergy();
            tpose = false;
        }
        Picture picture = new Picture(width, height);
        for (int col = 0; col < width; col++) {
            for (int row = 0; row < height; row++) {
                picture.setRGB(col, row, rgb[col][row]);
            }
        }
        return picture;
    }

    // width of current picture
    public int width() {
        return tpose ? height : width;
    }

    // height of current picture
    public int height() {
        return tpose ? width : height;
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || x >= (!tpose ? width : height) || y < 0 || y >= (!tpose ? height : width)) {
            throw new IllegalArgumentException("arguments are out of range");
        }
        return tpose ? energy[y][x] : energy[x][y];
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        if (!tpose) {
            transposeEnergy();
            tpose = true;
        }
        if (energy[0].length == 1) {
            int[] arr = new int[1];
            arr[0] = energy.length / 2;
            return arr;
        }
        return findVerticalSP();
    }

    // transposing the energy matrix
    private void transposeEnergy() {
        double[][] transpose = new double[height][width];
        int[][] transpose2 = new int[height][width];
        for (int col = 0; col < width; col++) {
            for (int row = 0; row < height; row++) {
                transpose[row][col] = energy[col][row];
                transpose2[row][col] = rgb[col][row];
            }
        }
        energy = new double[height][width];
        rgb = new int[height][width];
        for (int i = 0; i < height; i++) {
            System.arraycopy(transpose[i], 0, energy[i], 0, width);
            System.arraycopy(transpose2[i], 0, rgb[i], 0, width);
        }
        height = energy[0].length;
        width = energy.length;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        if (tpose) {
            transposeEnergy();
            tpose = false;
        }
        if (energy[0].length == 1) {
            int[] arr = new int[1];
            arr[0] = energy.length / 2;
            return arr;
        }
        return findVerticalSP();
    }

    // topological sort helper method
    private int[] findVerticalSP() {
        double distTo[] = new double[width * height];
        int row = 0;
        int edgeTo[] = new int[width * height];
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                distTo[j * width + i] = Double.POSITIVE_INFINITY;
            }
        }
        for (int c = 0; c < width; c++) {
            distTo[width * row + c] = 1000.0;
            edgeTo[width * row + c] = -1;
        }
        while (row < height) {
            int s = 0;
            while (s < width) {
                for (int i = Math.max(0, s - 1); i <= Math.min(width - 1, s + 1); i++)
                    relax(distTo, edgeTo, i, row + 1, s);
                s++;
            }
            row++;
        }
        int c2 = findMinimumCol(distTo);
        int[] sp = new int[height];
        int i = height - 1;
        for (int x = width * (height - 1) + c2; x != -1 && i >= 0; x = edgeTo[x]) {
            sp[i--] = x % width;
        }
        return sp;
    }

    // edge relaxation
    private void relax(double[] distTo, int[] edgeTo, int c, int r, int col) {
        if (c < 0 || c > width - 1 || r > height - 1) {
            return;
        }
        if (distTo[width * r + c] > distTo[width * (r - 1) + col] + energy[c][r]) {
            distTo[width * r + c] = distTo[width * (r - 1) + col] + energy[c][r];
            edgeTo[width * r + c] = width * (r - 1) + col;
        }
    }

    private int findMinimumCol(double[] distTo) {
        int min = 0;
        double mindist = Double.POSITIVE_INFINITY;
        for (int c = 0; c < width; c++) {
            if (distTo[width * (height - 1) + c] < mindist) {
                min = c;
                mindist = distTo[width * (height - 1) + c];
            }
        }
        return min;
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null) {
            throw new IllegalArgumentException("seam is null");
        }
        if (tpose) {
            transposeEnergy();
            tpose = false;
        }
        if (height <= 1) {
            throw new IllegalArgumentException("height is lenght 1");
        }
        if (seam.length != width) {
            throw new IllegalArgumentException("seam is incorrect" + seam.length + width);
        }
        for (int i = 1; i < seam.length; i++) {
            if (seam[i - 1] < 0 || seam[i] < 0 || seam[i - 1] >= height || seam[i] >= height
                    || seam[i - 1] - seam[i] < -1
                    || seam[i - 1] - seam[i] > 1) {
                throw new IllegalArgumentException("seam is incorrect");
            }
        }
        removeSeamHori(seam);
    }

    // helper method to delete the seam
    private void removeSeamHori(int[] seam) {
        int j = 0;
        for (int i = 0; i < seam.length; i++) {
            energy[j][Math.max(seam[i] - 1, 0)] = -1;
            energy[j][Math.min(seam[i] + 1, height - 1)] = -1;
            energy[Math.max(j - 1, 0)][seam[i]] = -1;
            energy[Math.min(j + 1, width - 1)][seam[i]] = -1;
            j++;
        }
        height--;
        double[][] arr = new double[width][height];
        int[][] rgbarr = new int[width][height];
        for (int i = 0; i < width; i++) {
            System.arraycopy(energy[i], 0, arr[i], 0, seam[i]);
            System.arraycopy(energy[i], seam[i] + 1, arr[i], seam[i], height - seam[i]);
            System.arraycopy(rgb[i], 0, rgbarr[i], 0, seam[i]);
            System.arraycopy(rgb[i], seam[i] + 1, rgbarr[i], seam[i], height - seam[i]);
        }
        energy = new double[width][height];
        rgb = new int[width][height];
        for (int i = 0; i < width; i++) {
            System.arraycopy(arr[i], 0, energy[i], 0, height);
            System.arraycopy(rgbarr[i], 0, rgb[i], 0, height);
        }
        recalculateEnergy();
    }

    // recalculate the energy
    private void recalculateEnergy() {
        for (int col = 0; col < width; col++) {
            for (int row = 0; row < height; row++) {
                if (energy[col][row] == -1) {
                    if (row == 0 || row == height - 1 || col == 0 || col == width - 1) {
                        energy[col][row] = 1000.00;
                        continue;
                    }
                    int left = rgb[col - 1][row];
                    int right = rgb[col + 1][row];
                    int up = rgb[col][row - 1];
                    int down = rgb[col][row + 1];
                    energy[col][row] = Math.sqrt(energyFunction(left, right) + energyFunction(up, down));
                }
            }
        }
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (seam == null) {
            throw new IllegalArgumentException("seam is null");
        }
        if (!tpose) {
            transposeEnergy();
            tpose = true;
        }
        if (height <= 1) {
            throw new IllegalArgumentException("height is lenght 1");
        }
        if (seam.length != width) {
            throw new IllegalArgumentException("seam is incorrect");
        }
        for (int i = 1; i < seam.length; i++) {
            if (seam[i - 1] < 0 || seam[i] < 0 || seam[i - 1] >= height || seam[i] >= height
                    || seam[i - 1] - seam[i] < -1
                    || seam[i - 1] - seam[i] > 1) {
                throw new IllegalArgumentException("seam is incorrect");
            }
        }
        removeSeamHori(seam);
    }

    // unit testing (optional)
    public static void main(String[] args) {
        if (args.length != 3) {
            StdOut.println("Usage:\njava SeamCaver [image filename] [num cols to remove] [num rows to remove]");
            return;
        }

        Picture inputImg = new Picture(args[0]);
        int removeColumns = Integer.parseInt(args[1]);
        int removeRows = Integer.parseInt(args[2]);

        StdOut.printf("image is %d columns by %d rows\n", inputImg.width(), inputImg.height());
        SeamCarver sc = new SeamCarver(inputImg);

        Stopwatch sw = new Stopwatch();

        for (int i = 0; i < removeRows; i++) {
            int[] horizontalSeam = sc.findHorizontalSeam();
            sc.removeHorizontalSeam(horizontalSeam);
        }

        for (int i = 0; i < removeColumns; i++) {
            int[] verticalSeam = sc.findVerticalSeam();
            sc.removeVerticalSeam(verticalSeam);
        }
        Picture outputImg = sc.picture();

        StdOut.printf("new image size is %d columns by %d rows\n", sc.width(), sc.height());

        StdOut.println("Resizing time: " + sw.elapsedTime() + " seconds.");
        inputImg.show();
        outputImg.show();
    }

}