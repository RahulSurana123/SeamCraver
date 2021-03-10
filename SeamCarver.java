import java.awt.Color;
import edu.princeton.cs.algs4.Picture;
import jdk.jfr.internal.PrivateAccess;

public class SeamCarver {

	Private Picture pic;
	
	private double[][] energy;
	
	private int[][] parent;
	
	// create a seam carver object based on the given picture
	public SeamCarver(Picture picture) {
		this.pic = picture;
		energy = new double[width()][height()];
		parent = new int[width()][height()];
		for (int y = 0; y < height(); y++) {
            for (int x = 0; x < width(); x++) {
                energy[x][y] = energy(x, y);
            }
        }	
	}

	// current picture
	public Picture picture() {
		return this.pic;
	}

	// width of current picture
	public int width() {
		return pic.width();
	}

	// height of current picture
	public int height() {
		return pic.height();
	}

	// energy of pixel at column x and row y
	public double energy(int x, int y) {
		
		if (x < 0 || x > width() - 1 || y < 0 || y > height() - 1) {
            throw new IndexOutOfBoundsException();
        }
        if (x == 0 || x == width() - 1 || y == 0 || y == height() - 1) {
            return 1000;
        }
		Color cx1 = pic.get(x - 1, y);
		Color cx2 = pic.get(x + 1, y);
		Color cy1 = pic.get(x, y - 1);
		Color cy2 = pic.get(x, y + 1);

		double ans = Math.sqrt(Math.pow(cx1.getRed() - cx2.getRed(), 2)
				+ Math.pow(cx1.getGreen() - cx2.getGreen(), 2)
				+ Math.pow(cx1.getBlue() - cx2.getBlue(), 2)
				+ Math.pow(cy1.getRed() - cy2.getRed(), 2)
				+ Math.pow(cy1.getGreen() - cy2.getGreen(), 2)
				+ Math.pow(cy1.getBlue() - cy2.getBlue(), 2));

		return ans;
	}

	// sequence of indices for horizontal seam
	public int[] findHorizontalSeam() {
		int[] removing_index = new int[height()];
		double[] distold = new double[width()];
		double[] distlatest = new double[width()];
		for (int i = 0; i < pic.width(); i++) {
			for (int j = 0; j < pic.height(); j++) {
					shortestPath(i, j, distlatest, distold);
			}
			System.arraycopy(distlatest, 0, distold, 0, width());
		}
		int best_index = 0;
		double least_energy = distlatest[0];
		for (int i = 0; i < width(); i++) {
			if(least_energy > distlatest[i]) {
				best_index = i;
				least_energy = distlatest[i];
			}
		}
		
		removing_index[height() - 1] = best_index;
        for (int i = height() - 2; i >= 0; i--) {
            removing_index[i] = parent[best_index][i + 1];
            best_index = parent[best_index][i + 1];
        }
		return removing_index;
	}
	
	
	
	private void shortestPath(int x, int y, double[] distTo, double[] oldDistTo) {
		if(x == 0) {
			double a = oldDistTo[x];
			double b = oldDistTo[x + 1];
			double min = Math.min(a, b);
			distTo[x] = min + energy[x][y];
			if (a == min) {
				parent[x][y] = x;
			}
			else {
				parent[x][y] = x + 1;
			}
		}
		if(x == width()-1) {
			double a = oldDistTo[x];
			double b = oldDistTo[x - 1];
			double min = Math.min(a, b);
			distTo[x] = min + energy[x][y];
			if (a == min) {
				parent[x][y] = x;
			}
			else {
				parent[x][y] = x - 1;
			}
		}
		
		double a = oldDistTo[x];
		double b = oldDistTo[x - 1];
		double c = oldDistTo[x + 1];
		double min = Math.min(Math.min(a, b), c);
		
		distTo[x] = min + energy[x][y];
		
		if (a == min) {
			parent[x][y] = x;
		}
		else if(c == min) {
			parent[x][y] = x + 1;
		}
		else {
			parent[x][y] = x - 1;
		}
	}

	
	// sequence of indices for vertical seam
	public int[] findVerticalSeam() {
		rotate();
		int[] removing_index = findHorizontalSeam();
		rotate();
		return removing_index;
	}


	private void rotate() {
		Picture Picture = new Picture(pic.height(), pic.width());
        double[][] newEnergy = new double[pic.height()][pic.width()];
        for (int i = 0; i < pic.width(); i++)
            for (int k = 0; k < pic.height(); k++) {
                Picture.set(k, i, pic.get(i, k));
                newEnergy[k][i] = energy[i][k];
            }
        energy = newEnergy;
        pic = Picture;
        parent = new int[pic.width()][pic.height()];
	}
	
	private void check(int[] seam) {
        if (width() <= 1 || height() <= 1) {
            throw new IllegalArgumentException("The width and height of the picture must be greatern than 1");
        }
        if (seam.length <= 1) {
            throw new IllegalArgumentException("The seam size must be greater than 1.");
        }

        for (int i = 0; i < seam.length - 1; i++) {
            if (Math.abs(seam[i] - seam[i + 1]) > 1) {
                throw new IllegalArgumentException();
            }
        }
    }
	
	// remove horizontal seam from current picture
	public void removeHorizontalSeam(int[] seam) {
		check(seam);
        if (seam.length > height()) {
            throw new IllegalArgumentException("The seam must not be greater than the image height!");
        }

        pic = removeSeam(seam, false);
        energy = new double[width()][height()];

        for (int y = 0; y < height(); y++) {
            for (int x = 0; x < width(); x++) {
                energy[x][y] = energy(x, y);
            }
        }
	}

	// remove vertical seam from current picture
	public void removeVerticalSeam(int[] seam) {
		check(seam);
        if (seam.length > height()) {
            throw new IllegalArgumentException("The seam must not be greater than the image height!");
        }

        pic = removeSeam(seam, true);
        energy = new double[width()][height()];

        for (int y = 0; y < height(); y++) {
            for (int x = 0; x < width(); x++) {
                energy[x][y] = energy(x, y);
            }
        }
	}

	
	private Picture removeSeam(int[] seam, boolean vertical) {
        if (vertical) {
            Picture p = new Picture(width() - 1, height());
            for (int y = 0; y < height(); y++) {
                int k = 0;
                for (int x = 0; x < width(); x++) {
                    if (x != seam[y]) {
                        p.set(k, y, pic.get(x, y));
                        k++;
                    }
                }
            }
            return p;
        }

        Picture p = new Picture(width(), height() - 1);
        for (int y = 0; y < width(); y++) {
            int k = 0;
            for (int x = 0; x < height(); x++) {
                if (x != seam[y]) {
                    p.set(y, k, pic.get(y, x));
                    k++;
                }
            }
        }
        return p;
    }
	
	// unit testing (optional)
	public static void main(String[] args) {
	}

}