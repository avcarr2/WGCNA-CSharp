using System;
using System.Data;
using System.Text;
using System.Diagnostics;
using System.Linq;
using System.Reflection;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;
using System.Collections.Generic;
using System.Collections;
using System.Threading.Tasks;
using System.Threading; 
using Microsoft.VisualBasic;

namespace WGCNA_Basic
{
	public class Program
	{
		public double[,] TestAdjMatrix = { { 0, 1, 1, 1 },
										{ 1, 0, 0, 1 },
										{ 1, 0, 0, 1 },
										{ 1, 1, 1, 0 } };
		static void Main()
		{
			string[] accessions = { "P0001", "P0002", "P0003" }; 
			double[] intensities1 = { 15.0, 14.0, 13.0 };
			double[] intensities2 = { 16.0, 14.0, 12.0 };
			double[] intensities3 = { 17.0, 12, 11 };
			WGCNA sample1 = new WGCNA(accessions, intensities1); sample1.Sample = "Sample1"; sample1.ExperimentalGroup = "Group1"; 
			WGCNA sample2 = new WGCNA(accessions, intensities2); sample2.Sample = "Sample2"; sample2.ExperimentalGroup = "Group1"; 
			WGCNA sample3 = new WGCNA(accessions, intensities3); sample3.Sample = "Sample3"; sample3.ExperimentalGroup = "Group1";
			WGCNA sample4 = new WGCNA(accessions, intensities1); sample3.Sample = "Sample4"; sample3.ExperimentalGroup = "Group2";

			Experiment Experiment1 = new Experiment();
			Experiment1.AddSample(sample1);
			Experiment1.AddSample(sample2);
			Experiment1.AddSample(sample3);
			int beta = 12; 
			Experiment1.PerformWGCNA(beta); 
		}
		
	}
	public class WGCNA 
	{
		public string Sample { get; set; }
		public string ExperimentalGroup { get; set; }

		public WGCNA(string[] Accessions, double[] Intensities)
		{
			ExpData = new Dictionary<string, double>(); 
			for (int i = 0; i < Accessions.Length; i++)
			{
				ExpData.Add(Accessions[i], Intensities[i]);
			}
		}
		public Dictionary<string, double> ExpData { get; set; }
		public static void PrintDictionaryContentsToConsole(WGCNA wgcnaobject)
		{
			foreach (KeyValuePair<string, double> entry in wgcnaobject.ExpData)
				Console.WriteLine("Accession = {0}; Intensity = {1}", entry.Key, entry.Value);
		}
		public static double PearsonCorrelation(double[] dataSet1, double[] dataSet2)
		{
			double meanX = dataSet1.Mean();
			double meanY = dataSet2.Mean();
			double numerator = new double();
			double denominator_xix = new double();
			double denominator_yiy = new double();
			double correlation = new double();
			for (int i = 0; i < dataSet1.Length; i++)
			{
				numerator += (dataSet1[i] - meanX) * (dataSet2[i] - meanY);
				denominator_xix += Math.Pow(dataSet1[i] - meanX, 2);
				denominator_yiy += Math.Pow(dataSet2[i] - meanY, 2);
			}
			return correlation = numerator / (Math.Sqrt(denominator_xix * denominator_yiy));
		}
		public static double PearsonCorrelation(double[][] dataset)
		{
			double[] dataset1 = dataset[0];
			double[] dataset2 = dataset[1];

			double meanX = dataset1.Mean();
			double meanY = dataset2.Mean();
			double numerator = new double();
			double denominator_xix = new double();
			double denominator_yiy = new double();
			double correlation = new double();
			for (int i = 0; i < dataset1.Length; i++)
			{
				numerator += (dataset1[i] - meanX) * (dataset2[i] - meanY);
				denominator_xix += Math.Pow(dataset1[i] - meanX, 2);
				denominator_yiy += Math.Pow(dataset2[i] - meanY, 2);
			}
			return correlation = numerator / (Math.Sqrt(denominator_xix * denominator_yiy));
		}
		public static double[] CorrelationFromTwoDictionaries(Dictionary<string, double[]> dict1, Dictionary<string, double[]> dict2)
		{
			// Takes a dictionary of the format <"accession", "intensities"> between two experimental conditions and returns the pearson correlation
			double[] result = new double[dict1.Count];
			double[][] entry = ConvertDictValuesToJaggedArray(dict1);
			double[][] entry2 = ConvertDictValuesToJaggedArray(dict2);

			for (int i = 0; i < result.Length; i++)
			{
				result[i] = PearsonCorrelation(entry[i], entry2[i]);
			}

			return result;
		}
		public static double[,] GenerateCorrelationMatrix(Dictionary<string, double[]> dict1)
		{
			double[,] CorrelationMatrix = new double[dict1.Count, dict1.Count];
			double[][] jaggedArray = ConvertDictValuesToJaggedArray(dict1);
			for (int i = 0; i < CorrelationMatrix.GetLength(0); i++)
			{
				for (int j = 0; j < CorrelationMatrix.GetLength(1); j++)
				{
					CorrelationMatrix[i, j] = PearsonCorrelation(jaggedArray[i], jaggedArray[j]);
				}
			}
			return CorrelationMatrix;
		}
		public static double[,] CreateAdjacencyMatrix(double[,] correlationMatrix, int beta)
		{
			// Creates a signed adjacency matrix
			double[,] result = new double[correlationMatrix.GetLength(0), correlationMatrix.GetLength(1)];
			for (int i = 0; i < correlationMatrix.GetLength(0); i++)
			{
				for (int j = 0; j < correlationMatrix.GetLength(1); j++)
				{
					result[i, j] = CalculateSignedAdjacency(correlationMatrix[i, j], beta);
				}
			}
			return result;
		}
		public static double CalculateSignedAdjacency(double corr, int beta)
		{
			double adj = Math.Pow((0.5 * (1 + corr)), beta);
			return adj;
		}
		public static double[][] ConvertDictValuesToJaggedArray(Dictionary<string, double[]> dict1)
		{
			double[][] entry = (new List<double[]>(dict1.Values)).ToArray();
			return entry;
		}
		public static double CalculatePearsonCorrelation(List<double> list1, List<double> list2)
		{
			double meanX = list1.Mean();
			double meanY = list2.Mean();
			double numerator = new double();
			double denominator_xix = new double();
			double denominator_yiy = new double(); 
			for(int i = 0; i < list1.Count; i++)
			{
				numerator += (list1.ElementAt(i) - meanX) * (list2.ElementAt(i) - meanY);
				denominator_xix += Math.Pow(list1.ElementAt(i) - meanX, 2);
				denominator_yiy += Math.Pow(list2.ElementAt(i) - meanY, 2); 
			}
			return numerator / (Math.Sqrt(denominator_xix * denominator_yiy)); 
		}



	}
	public class Experiment : IEnumerable<WGCNA> // An experiment is a collection of WGCNA objects
	{
		readonly IList<WGCNA> _wgcnaSamples;
		public double[,] CorrelationMatrix;
		public double[,] AdjacencyMatrix;
		public double[,] TOMMatrix; 
		public Experiment()
		{
			_wgcnaSamples = new List<WGCNA>(); 
		}
		public IEnumerator<WGCNA> GetEnumerator()
		{
			foreach(WGCNA sample in _wgcnaSamples)
			{
				yield return sample; 
			}
		}
		IEnumerator IEnumerable.GetEnumerator()
		{
			return GetEnumerator(); 
		}
		public void AddSample(WGCNA sample)
		{
			_wgcnaSamples.Add(sample); 
		}
	}
	public static class ExperimentHelpers 
	{ 
		public static List<List<double>> CombineSamplesFromExperiment(this Experiment exp)
		{
			List<List<double>> results = new List<List<double>>();
			foreach (WGCNA sample in exp)
			{
				results.Add(sample.ExpData.Values.ToList());
			}
			return results;
		}
		public static void CalculateCorrelationMatrix(this Experiment exp)
		{
			List<List<double>> intensitiesList = exp.CombineSamplesFromExperiment().CreateCorrelationMatrixPrecursor();
			int[] countIntensitiesListElements = new int[intensitiesList.Count]; 
			for(int i = 0; i < countIntensitiesListElements[0]; i++)
			{
				countIntensitiesListElements[i] = intensitiesList[i].Count; 
			}
			double[,] correlationMatrix = new double[countIntensitiesListElements[0], countIntensitiesListElements[0]];
			
			
		}
		public static List<List<double>> CreateCorrelationMatrixPrecursor(this List<List<double>> combinedSamples)
		{
			List<List<double>> result = new List<List<double>>();
			int[] count = new int[combinedSamples.Count]; 
			for(int i = 0; i < combinedSamples.Count; i++)
			{
				count[i] = combinedSamples[i].Count; 
			}
			// Collect the values from each sample in the collection
				for(int k = 0; k < count[0]; k++)
				{
					List<double> intermediateList = new List<double>();
					foreach (List<double> list in combinedSamples)
						{
							intermediateList.Add(list.ElementAt(k)); 
						}
					result.Add(intermediateList);
				}
			return result; 
		}
		public static int[] CountListsReturnInt(this List<List<double>> listsOfIntensities)
		{
			List<int> results = new List<int>(); 
			foreach(List<double> list in listsOfIntensities)
			{
				results.Add(list.Count); 
			}
			return results.ToArray(); 
		}
		public static double[,] CreatePearsonMatrix(this List<List<double>> intensitiesList)
		{
			int[] counts = intensitiesList.CountListsReturnInt();
			double[,] corrMatrix = new double[counts[0], counts[0]]; 
			
			for(int i = 0; i < counts[0]; i++)
			{
				for(int j = 0; j < counts[0]; j++)
				{
					corrMatrix[i,j] = WGCNA.CalculatePearsonCorrelation(intensitiesList.ElementAt(i), intensitiesList.ElementAt(j)); 
				}
			}
			return corrMatrix; 
		}
		public static double[,] CreateAdjacencyMatrix(this double[,] corrMatrix, int beta)
		{
			return WGCNA.CreateAdjacencyMatrix(corrMatrix, beta); 
		}
		public static double[,] CreateTOMMatrix(this double[,] adjMatrix)
		{
			return MatrixMath.CalculateTOM(adjMatrix); 
		}
		public static double[,] CreateDissTOMMatrix(this double[,] TOMMatrix)
		{
			return TOMMatrix.CalculateDissTOM(); 
		}
		public static void PerformWGCNA(this Experiment exp, int beta)
		{
			exp.CorrelationMatrix = exp.CombineSamplesFromExperiment().CreateCorrelationMatrixPrecursor().CreatePearsonMatrix();
			exp.AdjacencyMatrix = exp.CorrelationMatrix.CreateAdjacencyMatrix(beta);
			exp.TOMMatrix = exp.AdjacencyMatrix.CreateTOMMatrix().CreateDissTOMMatrix(); 
		}
	}
	public class WGCNAHelpers 
	{
		public static double[] RandomDoubleArray(int length, double[] range) //range needs to be of length 2
		{
			if (range.Length > 2)
			{
				throw new ArgumentException("range must be of length 2");
			}
			double[] resultsDouble = new double[length];
			for (int i = 0; i < length; i++)
			{
				resultsDouble[i] = RandomDoubleGenerator(range[0], range[1]);
			}
			return resultsDouble;
		}
		private static double RandomDoubleGenerator(double minValue, double maxValue)
		{
			Random random = new Random();
			double n = random.NextDouble() * (maxValue - minValue) + minValue;
			return n;
		}
		public static bool KeysEqual<TKey, TValue>(IDictionary<TKey, TValue> first, IDictionary<TKey, TValue> second)
		{
			if (first.Count != second.Count)
			{
				return false;
			}
			foreach (var key in first.Keys)
			{
				if (!second.ContainsKey(key))
				{
					return false;
				}
			}
			return true;

		}
		public static double[,] MakeRandomBinarySymmetricMatrix(int rows)
		{
			Random rng = new Random();
			double[,] results = new double[rows, rows]; 
			for(int i = 0; i < rows; i++)
			{
				for(int j = 0; j < rows; j++)
				{
					results[i, j] = rng.Next(0, 2);
					results[j, i] = results[i, j]; 
				}
			}
			return results; 
		}
	}
	public static class MatrixMath 
	{
		public static double[,] MatrixMultiply(double[,] mat, double number)
		{
			double[,] result = new double[mat.GetLength(0), mat.GetLength(1)];
			for (int i = 0; i < mat.GetLength(0); i++)
			{
				for (int j = 0; j < mat.GetLength(1); j++)
				{
					result[i, j] = number * mat[i, j];
				}
			}
			return result;
		}
		public static double DotProduct(double[] vec1, double[] vec2)
		{
			double result = new double();
			if (vec1.Length == vec2.Length)
			{
				for (int i = 0; i < vec1.Length; i++)
				{
					result += vec1[i] * vec2[i];
				}
			}
			return result;
		}
		public static double[,] MatrixMultiplication(double[,] mat1, double[,] mat2)
		{
			// For matrix multiplication, the number of columns in first matrix must equal 
			// the number of rows in the second matrix. The resulting matrix
			// has the number of rows of the first matrix, and columns of the second matrix
			int numberRows = mat1.GetLength(0);
			int numberCols = mat2.GetLength(1);
			bool checkEqualityOfRowsAndColumns = numberRows.Equals(numberCols);
			if (checkEqualityOfRowsAndColumns != true)
			{
				throw new ArgumentException();
			}

			double[,] result = new double[numberRows, numberCols];
			for (int i = 0; i < numberCols; i++)
			{
				for (int j = 0; j < numberRows; j++)
				{
					result[i, j] = DotProduct(AccessRow(mat1, i), AccessColumn(mat2, j));
				}
			}
			return result;
		}
		public static double[] RowSum(double[,] matrix)
		{
			double[] result = new double[matrix.GetLength(0)];

			for (int i = 0; i < matrix.GetLength(0); i++)
			{
				result[i] = AccessRow(matrix, i).Sum();
			}
			return result;
		}
		public static double[] ColSum(double[,] matrix)
		{
			double[] result = new double[matrix.GetLength(1)];

			for (int i = 0; i < matrix.GetLength(1); i++)
			{
				result[i] = AccessColumn(matrix, i).Sum();
			}
			return result;
		}
		public static double[,] ReplaceDiagonal(this double[,] matrix, double replacement)
		{
			double[,] result = new double[matrix.GetLength(0), matrix.GetLength(1)];
			for (int i = 0; i < matrix.GetLength(0); i++)
			{
				for (int j = 0; j < matrix.GetLength(1); j++)
				{
					if (i != j)
					{
						result[i, j] = matrix[i, j];

					}
					else result[i, j] = replacement;
				}
			}
			return result;
		}
		public static double[,] AddSingleNumberToMatrix(double[,] matrix, double number)
		{
			double[,] result = new double[matrix.GetLength(0), matrix.GetLength(1)];
			for (int i = 0; i < matrix.GetLength(0); i++)
			{
				for (int j = 0; j < matrix.GetLength(1); j++)
				{
					result[i, j] = matrix[i, j] + number;
				}
			}
			return result;
		}
		public static double[,] AddTwoMatrices(double[,] matrix1, double[,] matrix2)
		{
			bool equalityCheck = CheckDimensionEquality_Addition(matrix1, matrix2);
			if (equalityCheck == false)
			{
				throw new ArgumentException();
			}

			double[,] result = new double[matrix1.GetLength(0), matrix1.GetLength(0)];
			for (int i = 0; i < matrix1.GetLength(0); i++)
			{
				for (int j = 0; j < matrix1.GetLength(1); j++)
				{
					result[i, j] = matrix1[i, j] + matrix2[i, j];
				}
			}
			return result;
		}
		public static double[,] SubtractTwoMatrices(double[,] matrix1, double[,] matrix2)
		{
			bool equalityCheck = CheckDimensionEquality_Multiplication(matrix1, matrix2);
			if (equalityCheck == false)
			{
				throw new ArgumentException();
			}
			double[,] result = new double[matrix1.GetLength(0), matrix1.GetLength(0)];
			for (int i = 0; i < matrix1.GetLength(0); i++)
			{
				for (int j = 0; j < matrix1.GetLength(1); j++)
				{
					result[i, j] = matrix1[i, j] - matrix2[i, j];
				}
			}
			return result;
		}
		public static bool CheckDimensionEquality_Multiplication(double[,] matrix1, double[,] matrix2)
		{
			int dim1_mat1 = matrix1.GetLength(0);
			int dim2_mat2 = matrix2.GetLength(1);
			if (dim1_mat1 == dim2_mat2)
			{
				return true;
			}
			return false;
		}
		public static bool CheckDimensionEquality_Addition(double[,] matrix1, double[,] matrix2)
		{
			bool rowsEquality = matrix1.GetLength(0).Equals(matrix2.GetLength(0));
			bool colsEquality = matrix1.GetLength(1).Equals(matrix2.GetLength(1));
			if (rowsEquality == true && colsEquality == true)
			{
				return true;
			}
			else
			{
				return false;
			}

		}
		public static double[] AccessRow(double[,] matrix, int rowNumber)
		{
			return Enumerable.Range(0, matrix.GetLength(1))
				.Select(x => matrix[rowNumber, x])
				.ToArray();
		}
		public static double[] AccessColumn(double[,] matrix, int colNumber)
		{
			return Enumerable.Range(0, matrix.GetLength(1))
				.Select(x => matrix[x, colNumber])
				.ToArray();
		}
		public static double[,] FindTOMMinimums(double[,] adjacencyMatrix)
		{
			double[] k_j = RowSum(adjacencyMatrix);
			double[] k_i = ColSum(adjacencyMatrix);
			double[,] minimums = new double[k_i.Length, k_j.Length];
			for (int i = 0; i < k_i.Length; i++)
			{
				for (int j = 0; j < k_j.Length; j++)
				{
					double[] int_mins = { k_i[i], k_j[j] };
					minimums[i, j] = int_mins.Minimum();
				}
			}
			return minimums;
		}
		public static double[,] ElementWiseMatrixMultiplication(double[,] matrix1, double[,] matrix2)
		{
			bool dimensionCheck = CheckDimensionEquality_Multiplication(matrix1, matrix2);
			if (dimensionCheck == false)
			{
				throw new ArgumentException();
			}
			double[,] results = new double[matrix1.GetLength(0), matrix2.GetLength(1)];
			int rowsMatrix = matrix1.GetLength(0);
			int colsMatrix = matrix1.GetLength(1);
			for (int i = 0; i < rowsMatrix; i++)
			{
				for (int j = 0; j < colsMatrix; j++)
				{
					results[i, j] = matrix1[i, j] * matrix2[i, j];
				}
			}
			return results;
		}
		public static double[,] MatrixDivision(double[,] matrix1, double[,] matrix2)
		{
			bool dimensionCheck = CheckDimensionEquality_Multiplication(matrix1, matrix2);
			if (dimensionCheck == false)
			{
				throw new ArgumentException();
			}
			int rowsMatrix = matrix1.GetLength(0);
			int colsMatrix = matrix1.GetLength(1);
			double[,] results = new double[rowsMatrix, colsMatrix];
			for (int i = 0; i < rowsMatrix; i++)
			{
				for (int j = 0; j < colsMatrix; j++)
				{
					results[i, j] = matrix1[i, j] / matrix2[i, j];
				}
			}
			return results;
		}
		public static double[,] TransposeMatrix(this double[,] matrix)
		{
			int rows = matrix.GetLength(0);
			int cols = matrix.GetLength(1);

			double[,] result = new double[rows, cols];
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
				{
					result[j, i] = matrix[i, j];
				}
			}
			return result;
		}
		public static double[,] CalculateTOM(double[,] adjacencyMatrix)
		{
			double[,] L = MatrixMultiplication(adjacencyMatrix, adjacencyMatrix);
			double[] K = ColSum(adjacencyMatrix);
			double[,] _TOM = new double[adjacencyMatrix.GetLength(0), adjacencyMatrix.GetLength(1)];
			for (int i = 0; i < adjacencyMatrix.GetLength(0); i++)
			{
				double[] int_min = new double[2];
				for (int j = 0; j < adjacencyMatrix.GetLength(0); j++)
				{
					int_min[0] = K[i];
					int_min[1] = K[j];
					double numerator = L[i, j] + adjacencyMatrix[i, j];
					double denominator = int_min.Minimum() + 1 - adjacencyMatrix[i, j];
					_TOM[i, j] = numerator / denominator;

				}
			}
			// return AddTwoMatrices(_TOM, adjacencyMatrix.TransposeMatrix()).ReplaceDiagonal(1); 
			return _TOM.ReplaceDiagonal(1);
		}
		public static double[,] CalculateDissTOM(this double[,] TOMsimilarity)
		{
			return AddSingleNumberToMatrix(TOMsimilarity, -1);
		}
		public static double[,] MatrixMultiplicationParallel(double[,] matrix1, double[,] matrix2)
		{
			int rowsMatrix1 = matrix1.GetLength(0);
			int colsMatrix1 = matrix1.GetLength(1);
			int colsMatrix2 = matrix2.GetLength(1);
			double[,] result = new double[rowsMatrix1, colsMatrix2];

			if (CheckDimensionEquality_Multiplication(matrix1, matrix2) != true)
				throw new ArgumentException();

			Parallel.For(0, rowsMatrix1, i =>
			{
				for (int j = 0; j < colsMatrix2; j++)
				{
					result[i, j] = DotProduct(AccessRow(matrix1, i), AccessColumn(matrix2, j));
				}
			});
			return result;
		}
		public static double[,] ElementWiseMatrixMultiplicationParallel(double[,] matrix1, double[,] matrix2)
		{
			bool dimensionCheck = CheckDimensionEquality_Multiplication(matrix1, matrix2);
			if (dimensionCheck == false)
			{
				throw new ArgumentException();
			}
			double[,] results = new double[matrix1.GetLength(0), matrix2.GetLength(1)];
			int rowsMatrix = matrix1.GetLength(0);
			int colsMatrix = matrix1.GetLength(1);

			Parallel.For(0, rowsMatrix, i =>
			{
				for (int j = 0; j < colsMatrix; j++)
				{
					results[i, j] = matrix1[i, j] + matrix2[i, j];
				}
			});
			return results;
		}
		public static double[,] ReplaceDiagonalParallel(this double[,] matrix1, double replacement)
		{
			int rowsMatrix = matrix1.GetLength(0);
			int colsMatrix = matrix1.GetLength(1);
			double[,] result = new double[rowsMatrix, colsMatrix];
			Parallel.For(0, rowsMatrix, i =>
			{
				for (int j = 0; j < colsMatrix; j++)
				{
					if (i != j)
						result[i, j] = matrix1[i, j];
					else result[i, j] = replacement;
				}
			});
			return result;
		}
		public static double[,] CalculateTOMParallel(this double[,] adjacencyMatrix)
		{
			int rowsMatrix = adjacencyMatrix.GetLength(0);
			int colsMatrix = adjacencyMatrix.GetLength(1); 

			double[,] L = MatrixMultiplicationParallel(adjacencyMatrix, adjacencyMatrix);
			double[] K = ColSum(adjacencyMatrix);
			double[,] _TOM = new double[rowsMatrix, colsMatrix];

			Parallel.For(0, rowsMatrix, i =>
			{
				double[] connectivityMinimums = new double[2]; 
				for(int j = 0; j < colsMatrix; j++)
				{
					connectivityMinimums[0] = K[i];
					connectivityMinimums[1] = K[j];
					double numerator = L[i, j] + adjacencyMatrix[i, j];
					double denominator = connectivityMinimums.Min() + 1 - adjacencyMatrix[i, j];
					_TOM[i, j] = numerator / denominator; 
				}
			});
			return _TOM.ReplaceDiagonalParallel(1); 
		}         

	} 
	public class HClust
	{
		public string ClusteringAlgorithm { get; set; }
		public double[,] Merge { get; set; }
		public double[] Height { get; set; }
		public double[] Order { get; set; }
		public string[] Labels { get; set; }
		public string Method { get; set; }
		public string DistanceMetric { get; set; }
	}
}