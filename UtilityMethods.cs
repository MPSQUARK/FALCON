using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;


namespace MachineLearningSpectralFittingCode
{
    public class UtilityMethods
    {

        public static float[] ReadData(string Path)
        {

            string readText = File.ReadAllText(Path);

            readText = readText.Replace("\n", " ");
            float[] Data = Array.ConvertAll(readText.Split(' '), float.Parse);

            //for (int i = 0; i < Data.Length/2; i++)
            //{
            //    Console.WriteLine(Data[i]);
            //}

            return Data;
        } 

        public static float Degree2Radians(float Deg)
        {
            return Deg * Constants.Deg2RadFactor;
        }

        public static float Mpc2cm(float Mpc)
        {
            return 3.08567758128e+24f * Mpc;
        }

        public static float Median(float[] array)
        {
            if (array.Length % 2 == 0)
            {
                float[] arr = array.OrderBy(x => x).ToArray();
                return (arr[(int)(array.Length * 0.5f - 1)] + arr[(int)(array.Length * 0.5f)]) * 0.5f;
            }

            return array.OrderBy(x => x).ToArray()[(int)(array.Length * 0.5f)];
        }
    }

    public sealed class Tridiagonal
    {
        public float[] Solve(float[,] matrix, float[] d)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            int len = d.Length;
            // this must be a square matrix
            if (rows == cols && rows == len)
            {

                float[] b = new float[rows];
                float[] a = new float[rows];
                float[] c = new float[rows];

                // decompose the matrix into its tri-diagonal components
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        if (i == j)
                            b[i] = matrix[i, j];
                        else if (i == (j - 1))
                            c[i] = matrix[i, j];
                        else if (i == (j + 1))
                            a[i] = matrix[i, j];
                    }
                }
                try
                {
                    c[0] = c[0] / b[0];
                    d[0] = d[0] / b[0];

                    for (int i = 1; i < len - 1; i++)
                    {
                        c[i] = c[i] / (b[i] - a[i] * c[i - 1]);
                        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
                    }
                    d[len - 1] = (d[len - 1] - a[len - 1] * d[len - 2]) / (b[len - 1] - a[len - 1] * c[len - 2]);

                    // back-substitution step
                    for (int i = (len - 1); i-- > 0;)
                    {
                        d[i] = d[i] - c[i] * d[i + 1];
                    }

                    return d;
                }
                catch (DivideByZeroException)
                {
                    Console.WriteLine("Division by zero was attempted. Most likely reason is that ");
                    Console.WriteLine("the tridiagonal matrix condition ||(b*i)||>||(a*i)+(c*i)|| is not satisified.");
                    return null;
                }
            }
            else
            {
                Console.WriteLine("Error: the matrix must be square. The vector must be the same size as the matrix length.");
                return null;
            }
        }
    }

    internal sealed class CubicSplineInterpolation 
    {
        private bool baseset = false;
        private int len;
        private List<float> lX;
        private float[] X;
        private float[] Y;

        public CubicSplineInterpolation(float[] _x, float[] _y)
        {
            X = _x;
            Y = _y;

            len = X.Length;
            if (len > 1)
            {

                Console.WriteLine("Successfully set abscissa and ordinate.");
                baseset = true;
                // make a copy of X as a list for later use
                lX = _x.ToList();
            }
        }

        //public override float? Interpolate(float p)
        //{
        //    if (baseset)
        //    {
        //        float? result = 0;
        //        int N = len - 1;

        //        float[] h = Diff(X);
        //        float[] D = Scale(Diff(Y),h, false);
        //        float[] s = Enumerable.Repeat(3f, D.Length).ToArray();
        //        float[] dD3 = Scale(Diff(D),s);
        //        float[] a = Y;

        //        // generate tridiagonal system
        //        float[,] H = new float[N - 1, N - 1];
        //        float[] diagVals = new float[N - 1];
        //        for (int i = 1; i < N; i++)
        //        {
        //            diagVals[i - 1] = 2 * (h[i - 1] + h[i]);
        //        }

        //        H = Diag(H, diagVals);

        //        // H can be null if non-square matrix is passed
        //        if (H != null)
        //        {
        //            for (int i = 0; i < N - 2; i++)
        //            {
        //                H[i, i + 1] = h[i + 1];
        //                H[i + 1, i] = h[i + 1];
        //            }

        //            float[] c = new float[N + 2];
        //            c = Enumerable.Repeat(0f, N + 1).ToArray();

        //            // solve tridiagonal matrix
        //            Tridiagonal le = new Tridiagonal();
        //            float[] solution = le.Solve(H, dD3);

        //            for (int i = 1; i < N; i++)
        //            {
        //                c[i] = solution[i - 1];
        //            }

        //            float[] b = new float[N];
        //            float[] d = new float[N];
        //            for (int i = 0; i < N; i++)
        //            {
        //                b[i] = D[i] - (h[i] * (c[i + 1] + 2f * c[i])) / 3f;
        //                d[i] = (c[i + 1] - c[i]) / (3f * h[i]);
        //            }

        //            float Rx;

        //            try
        //            {
        //                // point p may be outside abscissa's range
        //                // if it is, we return null
        //                Rx = X.First(m => m >= p);
        //            }
        //            catch
        //            {
        //                return null;
        //            }

        //            // at this stage we know that Rx contains a valid value
        //            // find the index of the value close to the point required to be interpolated for
        //            int iRx = lX.IndexOf(Rx);

        //            if (iRx == -1)
        //                return null;

        //            if (iRx == len - 1 && X[iRx] == p)
        //                return Y[len - 1];

        //            if (iRx == 0)
        //                return Y[0];

        //            iRx = lX.IndexOf(Rx) - 1;
        //            Rx = p - X[iRx];
        //            result = a[iRx] + Rx * (b[iRx] + Rx * (c[iRx] + Rx * d[iRx]));

        //            return result;
        //        }
        //        else
        //        {
        //            return null;
        //        }


        //    }
        //    else
        //        return null;
        //}

        public float? Interpolate(float p)
        {
            if (baseset)
            {
                float? result = 0;
                int N = len - 1;

                float[] h = Diff(X);
                float[] D = Scale(Diff(Y), h, false);
                float[] s = Enumerable.Repeat(3f, D.Length).ToArray();
                float[] dD3 = Scale(Diff(D), s);
                float[] a = Y;

                // generate tridiagonal system
                float[,] H = new float[N - 1, N - 1];
                float[] diagVals = new float[N - 1];
                for (int i = 1; i < N; i++)
                {
                    diagVals[i - 1] = 2 * (h[i - 1] + h[i]);
                }

                H = Diag(H, diagVals);

                // H can be null if non-square matrix is passed
                if (H != null)
                {
                    for (int i = 0; i < N - 2; i++)
                    {
                        H[i, i + 1] = h[i + 1];
                        H[i + 1, i] = h[i + 1];
                    }

                    float[] c = new float[N + 2];
                    c = Enumerable.Repeat(0f, N + 1).ToArray();

                    // solve tridiagonal matrix
                    Tridiagonal le = new Tridiagonal();
                    float[] solution = le.Solve(H, dD3);

                    for (int i = 1; i < N; i++)
                    {
                        c[i] = solution[i - 1];
                    }

                    float[] b = new float[N];
                    float[] d = new float[N];
                    for (int i = 0; i < N; i++)
                    {
                        b[i] = D[i] - (h[i] * (c[i + 1] + 2f * c[i])) / 3f;
                        d[i] = (c[i + 1] - c[i]) / (3f * h[i]);
                    }

                    float Rx;

                    try
                    {
                        // point p may be outside abscissa's range
                        // if it is, we return null
                        Rx = X.First(m => m >= p);
                    }
                    catch
                    {
                        return null;
                    }

                    // at this stage we know that Rx contains a valid value
                    // find the index of the value close to the point required to be interpolated for
                    int iRx = lX.IndexOf(Rx);

                    if (iRx == -1)
                        return 0;

                    if (iRx == len - 1 && X[iRx] == p)
                        return Y[len - 1];

                    if (iRx == 0)
                        return Y[0];

                    iRx = lX.IndexOf(Rx) - 1;
                    Rx = p - X[iRx];
                    result = a[iRx] + Rx * (b[iRx] + Rx * (c[iRx] + Rx * d[iRx]));

                    return result;
                }
                else
                {
                    return null;
                }


            }
            else
                return null;
        }

        private float[,] Diag(float[,] matrix, float[] diagVals)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            // the matrix has to be scare
            if (rows == cols)
            {
                float[,] diagMatrix = new float[rows, cols];
                int k = 0;
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                    {
                        if (i == j)
                        {
                            diagMatrix[i, j] = diagVals[k];
                            k++;
                        }
                        else
                        {
                            diagMatrix[i, j] = 0;
                        }
                    }
                return diagMatrix;
            }
            else
            {
                Console.WriteLine("Diag should be used on square matrix only.");
                return null;
            }


        }

        private float[] Diff(float[] array)
        {
            int len = array.Length - 1;
            float[] diffsArray = new float[len];
            for (int i = 1; i <= len; i++)
            {
                diffsArray[i - 1] = array[i] - array[i - 1];
            }
            return diffsArray;
        }

        public float[] Scale(float[] array, float[] scalor, bool mult = true)
        {
            int len = array.Length;
            float[] scaledArray = new float[len];

            if (mult)
            {
                for (int i = 0; i < len; i++)
                {
                    scaledArray[i] = array[i] * scalor[i];
                }
            }
            else
            {
                for (int i = 0; i < len; i++)
                {
                    if (scalor[i] != 0)
                    {
                        scaledArray[i] = array[i] / scalor[i];
                    }
                    else
                    {
                        // basic fix to prevent division by zero
                        scalor[i] = 0.00001f;
                        scaledArray[i] = array[i] / scalor[i];

                    }
                }
            }

            return scaledArray;
        }

    }

    //internal interface IInterpolate
    //{
    //    float? Interpolate(float p);
    //}

    //internal abstract class Interpolation : IInterpolate
    //{

    //    public Interpolation(float[] _x, float[] _y)
    //    {
    //        int xLength = _x.Length;
    //        if (xLength == _y.Length && xLength > 1 && _x.Distinct().Count() == xLength)
    //        {
    //            x = _x;
    //            y = _y;
    //        }
    //    }

    //    // cubic spline relies on the abscissa values to be sorted
    //    public Interpolation(float[] _x, float[] _y, bool checkSorted = true)
    //    {
    //        int xLength = _x.Length;
    //        if (checkSorted)
    //        {
    //            if (xLength == _y.Length && xLength > 1 && _x.Distinct().Count() == xLength && Enumerable.SequenceEqual(SortedList(_x), _x.ToList()))
    //            {
    //                x = _x;
    //                y = _y;
    //            }
    //        }
    //        else
    //        {
    //            if (xLength == _y.Length && xLength > 1 && _x.Distinct().Count() == xLength)
    //            {
    //                x = _x;
    //                y = _y;
    //            }
    //        }
    //    }

    //    public float[] X
    //    {
    //        get
    //        {
    //            return x;
    //        }
    //    }

    //    public float[] Y
    //    {
    //        get
    //        {
    //            return y;
    //        }
    //    }

    //    public abstract float? Interpolate(float p);

    //    private float[] x;
    //    private float[] y;

    //    private List<T> SortedList<T>(this T[] array)
    //    {
    //        List<T> l = array.ToList();
    //        l.Sort();
    //        return l;

    //    }


    //}

    

}
