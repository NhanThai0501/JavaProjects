/**
 * 
 */
// package cs3310;

/**
 * @author NhanThai
 *
 */

// IMPORTANT NOTES: I created this file in Java Eclipse IDE
import java.util.Random;
import java.util.Scanner;
import java.util.function.BiFunction;

public class Project_One {
	public static void main(String[] args) {
		
		int min = 0;
		int max = 10;
		Scanner scan = new Scanner(System.in);

		System.out.println("Enter the size:");
		int size = scan.nextInt();
		
		matrixCalculation(size, min, max);
		timeCalculation(size);

	} // end main
	
	public static void matrixCalculation(int size, int min, int max) {
		Random rand = new Random();
		Scanner scan = new Scanner(System.in);
		Matrix cl = new Matrix(size);
		Matrix dc = new Matrix(size);
		Matrix st = new Matrix(size);
		
		System.out.println("Enter number of matrix 1:");
		int[][] A = new int[size][size];
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			//	A[i][j] = scan.nextInt();
			//	A[i][j] = rand.nextInt(min, max);
				A[i][j] = rand.nextInt(max);

		cl.printMatrix(A); System.out.println();
		cl.printMatrix_standard(A);
		

		System.out.println("Enter number of matrix 2:");
		int[][] B = new int[size][size];
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			//	B[i][j] = scan.nextInt();
			//	B[i][j] = rand.nextInt(min, max);
			B[i][j] = rand.nextInt(max);
		
		cl.printMatrix(B); System.out.println();
		cl.printMatrix_standard(B);

		int[][] classic = cl.classicMultiplication(A, B);

		System.out.println("\nResult of Classical Matrix Multiplication : ");
		Matrix.printMatrix(classic); System.out.println();
		cl.printMatrix_standard(classic);

		int[][] divi = dc.divideConquerMultiplication(A, B);
		System.out.println("\nResult of Divide and Conquer Matrix Multiplication : ");
		Matrix.printMatrix(divi); System.out.println();
		dc.printMatrix_standard(divi);

		int[][] strn = st.strassenMultiplication(A, B);
		System.out.println("\nResult of Strassen Matrix Multiplication : ");
		Matrix.printMatrix(strn); System.out.println();
		st.printMatrix_standard(strn);
	}


	public static void timeCalculation(int size) {
	      
        int trials = 15;
        int min = 1;
        int max = 9;

        // Use the same random matrices for all tests
        Random rand = new Random();
        int[][] A = generateRandomMatrix(size, min, max, rand);
        int[][] B = generateRandomMatrix(size, min, max, rand);

        System.out.print("Number of trials: " + trials + "\n");

        Matrix cl = new Matrix(size);
        Matrix dc = new Matrix(size);
        Matrix st = new Matrix(size);
		
        // Running time for each multiplication method
        Statistics classicalResult = statisticsMatrixMultiplication(A, B, trials, cl::classicMultiplication);
        Statistics divconResult = statisticsMatrixMultiplication(A, B, trials, dc::divideConquerMultiplication);
        Statistics strassenResult = statisticsMatrixMultiplication(A, B, trials, st::strassenMultiplication);

        // Print results
        printResults("Classical", classicalResult);
        printResults("Divide and Conquer", divconResult);
        printResults("strassenMultiplication", strassenResult);
    }

	
	//This function gathers statistics on matrix multiplication performance for a given multiplication method. 
	private static Statistics statisticsMatrixMultiplication(int[][] A, int[][] B, int trials, BiFunction<int[][], int[][], int[][]> multiplicationMethod) {
	    double bestTime = Double.MAX_VALUE; // Initialize with maximum value to easily find the minimum time.
	    double worstTime = Double.MIN_VALUE; // Initialize with minimum value to easily find the maximum time.
	    double sumTime = 0.0;

	    // Loop through the number of trials to perform multiplication and gather performance data.
	    for (int i = 0; i < trials; i++) {
	        long startTime = System.nanoTime(); // Capture the start time.
	        multiplicationMethod.apply(A, B);   // Execute the multiplication.
	        long endTime = System.nanoTime();   // Capture the end time.

	        double timeTaken = endTime - startTime; // Compute the time taken for this trial.

	        // Update best and worst times.
	        bestTime = Math.min(bestTime, timeTaken);
	        worstTime = Math.max(worstTime, timeTaken);

	        sumTime += timeTaken; // Sum up all the times for later calculation of average.
	    }

	    // Calculate average time taken, excluding the best and worst times.
	    // If there are less than 3 trials, then we simply compute the average of all trials.
	    double averageTime = (sumTime - bestTime - worstTime) / (trials - 2);

	    return new Statistics(bestTime, worstTime, averageTime); // Return the computed statistics.
	}

	
	// Generates a square matrix of a given size with random integer values between specified minimum and maximum.
	private static int[][] generateRandomMatrix(int size, int min, int max, Random rand) {
	    int[][] matrix = new int[size][size]; // Initialize the matrix.

	    // Loop through the matrix's rows.
	    for (int i = 0; i < size; i++) {
	        // Loop through the matrix's columns.
	        for (int j = 0; j < size; j++) {
	            // matrix[i][j] = rand.nextInt(min, max + 1); // Assign a random number. Note: max + 1 since Random.nextInt is exclusive of the upper bound.
				matrix[i][j] = rand.nextInt(max + 1);
	        }
	    }
	    return matrix; // Return the populated matrix.
	}


    private static void printResults(String method, Statistics result) {
        System.out.println("\nResult of " + method + " Matrix Multiplication:");
        System.out.println("Best case: " + result.bestTime + " ns");
        System.out.println("Worst case: " + result.worstTime + " ns");
        System.out.println("Average time: " + result.averageTime + " ns (excluding best and worst case)");
    }

    private static class Statistics {
        double bestTime;
        double worstTime;
        double averageTime;

        Statistics(double bestTime, double worstTime, double averageTime) {
            this.bestTime = bestTime;
            this.worstTime = worstTime;
            this.averageTime = averageTime;
        }
    }
    
    
//begin matrix class 
	public static class Matrix {
// start from 0 to count time 
		int classicCount = 0;
		int diviCount = 0;
		int straCount = 0;
		int size;
//create array for Matrix 
		int m1[][];
		int m2[][];

		public Matrix(int size) {
			this.size = size;
			m1 = new int[size][size];
			m2 = new int[size][size];
		}

		public void initializeMatrix() {
// using random for Matrix number 
			Random r1 = new Random();
			for (int a = 0; a < size; a++) {
				for (int b = 0; b < size; b++)

				{
					m1[a][b] = r1.nextInt(9);
					m2[a][b] = r1.nextInt(9);
				}
			}
		}

		public void initializeMatrix(int[][] copyFrom, int[][] copyTo, int startRow, int startCol) {
			for (int a = 0, a2 = startRow; a < copyTo.length; a++, a2++) { // for row
				for (int b = 0, b2 = startCol; b < copyTo.length; b++, b2++) { // for column
					copyTo[a][b] = copyFrom[a2][b2];
				}
			}
		}


		// Standard output
		public void printMatrix_standard(int matrix[][]) {
			for (int a = 0; a < matrix.length; a++) {
				for (int b = 0; b < matrix.length; b++) {

					System.out.print(matrix[a][b] + " ");
				}
				System.out.println();
			}
		}
		
		// Decorative output
	    public static void printMatrix(int[][] matrix) {
	        int maxLength = findLongestNumberLength(matrix);

	        for (int[] row : matrix) {
	            System.out.print("{");
	            for (int i = 0; i < row.length; i++) {
	                // Use the format string to align the number to the right
	                System.out.printf("%" + maxLength + "d", row[i]);

	                // Append a comma if it's not the last number in the row
	                if (i < row.length - 1) {
	                    System.out.print(", ");
	                }
	            }
	            System.out.println("}");
	        }
	    }

	    // Helper method to find the length of the longest number in the matrix
	    private static int findLongestNumberLength(int[][] matrix) {
	        int maxLength = 0;
	        for (int[] row : matrix) {
	            for (int number : row) {
	                maxLength = Math.max(maxLength, Integer.toString(number).length());
	            }
	        }
	        // Add space for comma and space
	        return maxLength; // one for the comma, one for the space after the comma
	    }
	    
	    // Classical multiplication
		public int[][] classicMultiplication(int[][] m1, int[][] m2) {
			classicCount++;
			int result[][] = new int[m1.length][m1.length];
			for (int i = 0; i < m1.length; i++)
				for (int j = 0; j < m1.length; j++) {
					result[i][j] = 0;
					for (int k = 0; k < m1.length; k++)
						result[i][j] = result[i][j] + m1[i][k] * m2[k][j];
				}
			return result;
		}

	    // divide and conquer multiplication
		public int[][] divideConquerMultiplication(int[][] a, int[][] b) {
			diviCount++;
			int n = a.length;
			int[][] result = new int[n][n];
			if (n == 1) {
				result[0][0] = a[0][0] * b[0][0];
			} else {
				int[][] a11 = new int[n / 2][n / 2];
				int[][] a12 = new int[n / 2][n / 2];
				int[][] a21 = new int[n / 2][n / 2];
				int[][] a22 = new int[n / 2][n / 2];
				int[][] b11 = new int[n / 2][n / 2];
				int[][] b12 = new int[n / 2][n / 2];
				int[][] b21 = new int[n / 2][n / 2];
				int[][] b22 = new int[n / 2][n / 2];

				initializeMatrix(a, a11, 0, 0);
				initializeMatrix(a, a12, 0, n / 2);
				initializeMatrix(a, a21, n / 2, 0);
				initializeMatrix(a, a22, n / 2, n / 2);
				initializeMatrix(b, b11, 0, 0);
				initializeMatrix(b, b12, 0, n / 2);
				initializeMatrix(b, b21, n / 2, 0);
				initializeMatrix(b, b22, n / 2, n / 2);

				int[][] c11 = addMatrix(divideConquerMultiplication(a11, b11), divideConquerMultiplication(a12, b21));
				int[][] c12 = addMatrix(divideConquerMultiplication(a11, b12), divideConquerMultiplication(a12, b22));
				int[][] c21 = addMatrix(divideConquerMultiplication(a21, b11), divideConquerMultiplication(a22, b21));
				int[][] c22 = addMatrix(divideConquerMultiplication(a21, b12), divideConquerMultiplication(a22, b22));
				attachMatrix(c11, result, 0, 0);
				attachMatrix(c12, result, 0, n / 2);
				attachMatrix(c21, result, n / 2, 0);
				attachMatrix(c22, result, n / 2, n / 2);
			}
			return result;

		}

		// strassen multiplication 
		public int[][] strassenMultiplication(int[][] a, int[][] b) {
			straCount++;
			int n = a.length;
			int[][] result = new int[n][n];
			if (n == 1) {
				result[0][0] = a[0][0] * b[0][0];
			} else {
				int[][] a11 = new int[n / 2][n / 2];
				int[][] a12 = new int[n / 2][n / 2];
				int[][] a21 = new int[n / 2][n / 2];
				int[][] a22 = new int[n / 2][n / 2];
				int[][] b11 = new int[n / 2][n / 2];
				int[][] b12 = new int[n / 2][n / 2];
				int[][] b21 = new int[n / 2][n / 2];
				int[][] b22 = new int[n / 2][n / 2];

				initializeMatrix(a, a11, 0, 0);
				initializeMatrix(a, a12, 0, n / 2);
				initializeMatrix(a, a21, n / 2, 0);
				initializeMatrix(a, a22, n / 2, n / 2);
				initializeMatrix(b, b11, 0, 0);
				initializeMatrix(b, b12, 0, n / 2);
				initializeMatrix(b, b21, n / 2, 0);
				initializeMatrix(b, b22, n / 2, n / 2);

				int[][] m1 = strassenMultiplication(addMatrix(a11, a22), addMatrix(b11, b22));
				int[][] m2 = strassenMultiplication(addMatrix(a21, a22), b11);
				int[][] m3 = strassenMultiplication(a11, subtractMatrix(b12, b22));
				int[][] m4 = strassenMultiplication(a22, subtractMatrix(b21, b11));
				int[][] m5 = strassenMultiplication(addMatrix(a11, a12), b22);
				int[][] m6 = strassenMultiplication(subtractMatrix(a21, a11), addMatrix(b11, b12));
				int[][] m7 = strassenMultiplication(subtractMatrix(a12, a22), addMatrix(b21, b22));

				int[][] c11 = addMatrix(subtractMatrix(addMatrix(m1, m4), m5), m7);
				int[][] c12 = addMatrix(m3, m5);
				int[][] c21 = addMatrix(m2, m4);
				int[][] c22 = addMatrix(subtractMatrix(addMatrix(m1, m3), m2), m6);

				attachMatrix(c11, result, 0, 0);
				attachMatrix(c12, result, 0, n / 2);
				attachMatrix(c21, result, n / 2, 0);
				attachMatrix(c22, result, n / 2, n / 2);

			}
			return result;

		}
		
		//This function attaches (or copies) a smaller matrix (child) into a specific location 
		// of a larger matrix (parent) starting from the specified row and column.
		public void attachMatrix(int[][] child, int[][] parent, int startRow, int startCol) {
			for (int i = 0, i2 = startRow; i < child.length; i++, i2++) {
				for (int j = 0, j2 = startCol; j < child.length; j++, j2++) {
					parent[i2][j2] = child[i][j];
				}
			}
		}

		// Add two matrices
		public int[][] addMatrix(int[][] m1, int[][] m2) {
			int result[][] = new int[m1.length][m1.length];
			for (int a = 0; a < m1.length; a++) {
				for (int b = 0; b < m1.length; b++) {
					result[a][b] = m1[a][b] + m2[a][b];
				}
			}
			return result;
		}

		// Subtract two matrices
		public int[][] subtractMatrix(int[][] m1, int[][] m2) {
			int result[][] = new int[m1.length][m1.length];
			for (int a = 0; a < m1.length; a++) {
				for (int b = 0; b < m1.length; b++) {
					result[a][b] = m1[a][b] - m2[a][b];
				}
			}
			return result;
		}
	}// end class Matrix

} // end class Project_One
