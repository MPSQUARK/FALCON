using ILGPU;
using ILGPU.Algorithms;
using ILGPU.Runtime;
using System;
using System.Linq;
using System.Threading.Tasks;

namespace FALCON.vector
{

    /* LOG :
     *      - Access Slice                          : WORKING
     *      - Access Value                          : WORKING
     *      - Consecutive Operation                 : WORKING
     *      - Consecutive Operation 2D              : NOT TESTED
     *      - Consecutive Compound Operation 2D     : NOT TESTED
     *      - Dot Product                           : WORKING
     *      - Fill                                  : WORKING
     *      - Scalar Operation                      : WORKING
     *      - Scalar Compound Operation             : WORKING
     *      - Normalise                             : NOT TESTED
     *      
    */

    public class Vector
    {
        public Vector(float[] value, int columns = 1)
        {
            this.Value = value;
            this.Columns = columns;
        }
        
        public float[] Value { get; set; } // Determines the values in the Vector
        public int Columns { get; set; } // Defines the number of Columns in a Vector STARTS AT 1


        // Creates a Uniform Vector where all values are = Value
        public static Vector Fill(float Value, int Size, int Columns = 1)
        {
            return new Vector(Enumerable.Repeat(Value, Size).ToArray(), Columns);
        }


        // Access 1 Value from 2D Vector
        public static float AccessVal(Vector vector, int row, int col)
        {
            return vector.Value[row*vector.Columns + col];
        }
        // Access a ROW or COLUMN from a 2D Vector
        public static Vector AccessSlice(Accelerator gpu, Vector vector, int row_col_index , char row_col)
        {

            if (vector.Columns == 1)
            {
                throw new Exception("Input Vector must be a 2D Vector");
            }

            int[] ChangeSelectLength;
            int OutPutVectorLength;

            switch (row_col)
            {
                case 'r':
                    ChangeSelectLength = new int[5] { 0, 1, row_col_index, 0, vector.Columns };
                    OutPutVectorLength = vector.Columns;
                    break;
                case 'c':
                    ChangeSelectLength = new int[5] { 1, 0, 0, row_col_index, vector.Columns };
                    OutPutVectorLength = vector.Value.Length / vector.Columns;
                    break;    
                default:
                    throw new Exception("Invalid slice char selector, choose 'r' for row or 'c' for column");
            }

            AcceleratorStream Stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<int>>(AccessSliceKernal);

            var buffer = gpu.Allocate1D(new float[OutPutVectorLength]);
            var buffer2 = gpu.Allocate1D(vector.Value);
            var buffer3 = gpu.Allocate1D(ChangeSelectLength);

            kernelWithStream(Stream, OutPutVectorLength, buffer.View, buffer2.View, buffer3.View);

            Stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(Stream);

            buffer.Dispose();
            buffer2.Dispose();
            buffer3.Dispose();

            Stream.Dispose();

            return new Vector(Output);
        }
        public static Vector AccessSliceAsync(Accelerator gpu, Vector vector, int row_col_index, char row_col)
        {            

            if (vector.Columns == 1)
            {
                throw new Exception("Input Vector must be a 2D Vector");
            }

            int[] ChangeSelectLength;
            int OutPutVectorLength;

            switch (row_col)
            {
                case 'r':
                    ChangeSelectLength = new int[5] { 0, 1, row_col_index, 0, vector.Columns };
                    OutPutVectorLength = vector.Columns;
                    break;
                case 'c':
                    ChangeSelectLength = new int[5] { 1, 0, 0, row_col_index, vector.Columns };
                    OutPutVectorLength = vector.Value.Length / vector.Columns;
                    break;
                default:
                    throw new Exception("Invalid slice char selector, choose 'r' for row or 'c' for column");
            }


            AcceleratorStream Stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<int>>(AccessSliceKernal);

            var buffer = gpu.Allocate1D(new float[OutPutVectorLength]);
            var buffer2 = gpu.Allocate1D(vector.Value);
            var buffer3 = gpu.Allocate1D(ChangeSelectLength);

            kernelWithStream(Stream, OutPutVectorLength, buffer.View, buffer2.View, buffer3.View);

            Stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(Stream);

            buffer.Dispose();
            buffer2.Dispose();
            buffer3.Dispose();

            Stream.Dispose();

            return new Vector(Output);
        }

        // KERNEL
        static void AccessSliceKernal(Index1D index, ArrayView<float> OutPut, ArrayView<float> Input, ArrayView<int> ChangeSelectLength)
        {
            OutPut[index] = Input[
                index * ChangeSelectLength[0] * ChangeSelectLength[4] + // iRcL
                index * ChangeSelectLength[1] +                         // iCc
                ChangeSelectLength[2] * ChangeSelectLength[4] +         // RsL
                ChangeSelectLength[3]];                                 // Cs
        }



        // SCALAR OPERATIONS : Vector * Scalar, Vector / Scalar, Vector +|- Scalar
        public static Vector ScalarOperation(Accelerator gpu, Vector vector, float scalar, string operation = "*")
        {

            AcceleratorStream Stream = gpu.CreateStream();

            var buffer = gpu.Allocate1D(new float[vector.Value.Length]);
            var buffer2 = gpu.Allocate1D(vector.Value);

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float>(ScalarProductKernal);

            switch (operation)
            {
                case "*":
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float>(ScalarProductKernal);
                    break;
                case "/":
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float>(ScalarDivideKernal);
                    break;
                case "+":
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float>(ScalarSumKernal);
                    break;
                case "^*":  // flip the Vector e.g. 1/Vector then multiply by Scalar
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float>(ScalarProductInvVecKernal);
                    break;
            }

            kernelWithStream(Stream, buffer.IntExtent, buffer.View, buffer2.View, scalar);

            Stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(Stream);

            buffer.Dispose();
            buffer2.Dispose();

            Stream.Dispose();

            return new Vector(Output);
        }
        // KERNELS
        static void ScalarProductKernal(Index1D index, ArrayView<float> OutPut, ArrayView<float> Input, float Scalar )
        {
            OutPut[index] = Input[index] * Scalar;
        }
        static void ScalarDivideKernal(Index1D index, ArrayView<float> OutPut, ArrayView<float> Input, float Scalar)
        {
            OutPut[index] = Input[index] / Scalar;
        }
        static void ScalarSumKernal(Index1D index, ArrayView<float> OutPut, ArrayView<float> Input, float Scalar)
        {
            OutPut[index] = Input[index] + Scalar;
        }
        static void ScalarProductInvVecKernal(Index1D index, ArrayView<float> OutPut, ArrayView<float> Input, float Scalar)
        {
            OutPut[index] = Scalar / Input[index];
        }



        // SCALAR OPERATIONS : Vector * Scalar, Vector / Scalar, Vector +|- Scalar
        public static double[] ScalarOperation(Accelerator gpu, Vector vector, double scalar, char operation = '*')
        {

            AcceleratorStream Stream = gpu.CreateStream();

            var buffer = gpu.Allocate1D(new double[vector.Value.Length]);
            var buffer2 = gpu.Allocate1D(vector.Value);

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<double>, ArrayView<float>, double>(ScalarProduct_DKernal);

            switch (operation)
            {
                case '*':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<double>, ArrayView<float>, double>(ScalarProduct_DKernal);
                    break;
                case '/':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<double>, ArrayView<float>, double>(ScalarDivide_DKernal);
                    break;
                case '+':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<double>, ArrayView<float>, double>(ScalarSum_DKernal);
                    break;
            }

            kernelWithStream(Stream, buffer.IntExtent, buffer.View, buffer2.View, scalar);

            Stream.Synchronize();

            double[] Output = buffer.GetAsArray1D(Stream);

            buffer.Dispose();
            buffer2.Dispose();

            Stream.Dispose();

            return Output;
        }
        // KERNELS
        static void ScalarProduct_DKernal(Index1D index, ArrayView<double> OutPut, ArrayView<float> Input, double Scalar)
        {
            OutPut[index] = Input[index] * Scalar;
        }
        static void ScalarDivide_DKernal(Index1D index, ArrayView<double> OutPut, ArrayView<float> Input, double Scalar)
        {
            OutPut[index] = Input[index] / Scalar;
        }
        static void ScalarSum_DKernal(Index1D index, ArrayView<double> OutPut, ArrayView<float> Input, double Scalar)
        {
            OutPut[index] = Input[index] + Scalar;
        }
        


        // COMPOUND SCALAR OPERATIONS : Vector * Scalar1 +|- Scaler2, Vector / Scalar1 +|- Scalar2
        /// <summary>
        /// 
        /// </summary>
        /// <param name="gpu"></param>
        /// <param name="vector"></param>
        /// <param name="Multiple"></param>
        /// <param name="Adder"></param>
        /// <param name="operation"></param>
        /// <returns></returns>
        public static Vector ScalarCompoundOperation(Accelerator gpu, Vector vector, float Multiple, float Adder, string operation = "*+")
        {

            AcceleratorStream Stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float, float>(ScalarProductSumKernal);

            var buffer = gpu.Allocate1D(new float[vector.Value.Length]);
            var buffer2 = gpu.Allocate1D(vector.Value);

            switch (operation)
            {
                case "*+":
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float, float>(ScalarProductSumKernal);
                    break;
                case "/+":
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float, float>(ScalarDivideSumKernal);
                    break;
                case "+*":
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float, float>(ScalarSumProductKernal);
                    break;
                case "+/":
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float, float>(ScalarSumDivideKernal);
                    break;
            }

            kernelWithStream(Stream, buffer.IntExtent, buffer.View, buffer2.View, Multiple, Adder);

            Stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(Stream);

            buffer.Dispose();
            buffer2.Dispose();

            Stream.Dispose();

            return new Vector(Output);
        }
        // KERNELS
        static void ScalarProductSumKernal(Index1D index, ArrayView<float> OutPut, ArrayView<float> Input, float Multiplicator, float Adder)
        {
            OutPut[index] = Input[index] * Multiplicator + Adder;
        }
        static void ScalarDivideSumKernal(Index1D index, ArrayView<float> OutPut, ArrayView<float> Input, float Divisor, float Adder)
        {
            OutPut[index] = Input[index] / Divisor + Adder;
        }
        static void ScalarSumProductKernal(Index1D index, ArrayView<float> OutPut, ArrayView<float> Input, float Multiplicator, float Adder)
        {
            OutPut[index] = (Input[index] + Adder) * Multiplicator;
        }
        static void ScalarSumDivideKernal(Index1D index, ArrayView<float> OutPut, ArrayView<float> Input, float Divisor, float Adder)
        {
            OutPut[index] = (Input[index] + Adder) / Divisor ;
        }


        // Multiplies 2 Vectors Element by Element
        public static Vector ConsecutiveOperation(Accelerator gpu, Vector vectorA, Vector vectorB, char operation = '*')
        {
            if (vectorA.Value.Length != vectorB.Value.Length)
            {
                throw new IndexOutOfRangeException("Vector A and Vector B provided MUST be of EQUAL length" );
            }

            AcceleratorStream Stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>>(ConsecutiveProductKernal);

            var buffer = gpu.Allocate1D(vectorA.Value); // Input
            var buffer2 = gpu.Allocate1D(vectorB.Value); // Input
            var buffer3 = gpu.Allocate1D(new float[vectorA.Value.Length]); // Output


            switch (operation)
            {
                case '*':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>>(ConsecutiveProductKernal);
                    break;
                case '+':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>>(ConsecutiveAdditionKernal);
                    break;
                case '-':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>>(ConsecutiveSubtractKernal);
                    break;
                case '/':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>>(ConsecutiveDivisionKernal);
                    break;
            }

            kernelWithStream(Stream, buffer.IntExtent, buffer.View, buffer2.View, buffer3.View);

            Stream.Synchronize();

            float[] Output = buffer3.GetAsArray1D(Stream);

            buffer.Dispose();
            buffer2.Dispose();
            buffer3.Dispose();

            Stream.Dispose();

            return new Vector(Output);
        }
        // KERNEL
        static void ConsecutiveProductKernal(Index1D index, ArrayView<float> InputA, ArrayView<float> InputB, ArrayView<float> OutPut )
        {

            OutPut[index] = InputA[index] * InputB[index];

        }
        static void ConsecutiveAdditionKernal(Index1D index, ArrayView<float> InputA, ArrayView<float> InputB, ArrayView<float> OutPut)
        {

            OutPut[index] = InputA[index] + InputB[index];

        }
        static void ConsecutiveSubtractKernal(Index1D index, ArrayView<float> InputA, ArrayView<float> InputB, ArrayView<float> OutPut)
        {

            OutPut[index] = InputA[index] - InputB[index];

        }
        static void ConsecutiveDivisionKernal(Index1D index, ArrayView<float> InputA, ArrayView<float> InputB, ArrayView<float> OutPut)
        {

            OutPut[index] = InputA[index] / InputB[index];

        }


        // Multiplies 2 Vectors Element by Element, where Vector A is 2D and Vector B is 1D
        public static Vector ConsecutiveOperation2D(Accelerator gpu, Vector vectorA, Vector vectorB, char operation = '*')
        {
            if (vectorA.Columns == 1)
            {
                throw new Exception("vectorA should be a 2D flattened vector of columns > 1");
            }
            if (vectorB.Columns > 1)
            {
                throw new Exception("vectorB should be a 1D vector of columns = 1");
            }
            if (vectorA.Columns != vectorB.Value.Length)
            {
                throw new Exception($"Length of VectorB : {vectorB.Value.Length} does NOT match number of columns in VectorA : {vectorA.Columns}");
            }

            AcceleratorStream Stream = gpu.CreateStream();
            
            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, float>(ConsecutiveProduct2DKernel);
            
            var buffer = gpu.Allocate1D(new float[vectorA.Value.Length]); // Output
            var buffer2 = gpu.Allocate1D(vectorA.Value); // Input
            var buffer3 = gpu.Allocate1D(vectorB.Value); // Input

            switch (operation)
            {
                case '*':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, float>(ConsecutiveProduct2DKernel);
                    break;
                case '+':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, float>(ConsecutiveSum2DKernel);
                    break;
                case '-':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, float>(ConsecutiveSubtract2DKernel);
                    break;
                case '/':
                    kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, float>(ConsecutiveDivide2DKernel);
                    break;
            }

            kernelWithStream(Stream, buffer.IntExtent, buffer.View, buffer2.View, buffer3.View, vectorA.Columns);

            Stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(Stream);

            buffer.Dispose();
            buffer2.Dispose();
            buffer3.Dispose();

            Stream.Dispose();

            return new Vector(Output, vectorA.Columns);
        }
        // KERNEL
        static void ConsecutiveProduct2DKernel(Index1D index, ArrayView<float> Output, ArrayView<float> InputA, ArrayView<float> InputB, float Cols)
        {
            Output[index] = InputA[index] * InputB[(int)XMath.RoundAwayFromZero(index % Cols)]; 
        }
        static void ConsecutiveSum2DKernel(Index1D index, ArrayView<float> Output, ArrayView<float> InputA, ArrayView<float> InputB, float Cols)
        {
            Output[index] = InputA[index] + InputB[(int)XMath.RoundAwayFromZero(index % Cols)];
        }
        static void ConsecutiveSubtract2DKernel(Index1D index, ArrayView<float> Output, ArrayView<float> InputA, ArrayView<float> InputB, float Cols)
        {
            Output[index] = InputA[index] - InputB[(int)XMath.RoundAwayFromZero(index % Cols)];
        }
        static void ConsecutiveDivide2DKernel(Index1D index, ArrayView<float> Output, ArrayView<float> InputA, ArrayView<float> InputB, float Cols)
        {
            Output[index] = InputA[index] / InputB[(int)XMath.RoundAwayFromZero(index % Cols)];
        }

        // Multiplies 2 Vectors Element by Element, where Vector A is 2D and Vector B is 1D
        public static Vector ConsecutiveCompoundOperation2D(Accelerator gpu, Vector vectorA, Vector vectorB, float scalar, string operation = "**")
        {
            if (vectorA.Columns == 1)
            {
                throw new Exception("vectorA should be a 2D flattened vector of columns > 1");
            }
            if (vectorB.Columns > 1)
            {
                throw new Exception("vectorB should be a 1D vector of columns = 1");
            }
            if ((vectorA.Value.Length / vectorA.Columns) != vectorB.Value.Length)
            {
                throw new Exception($"Length of VectorB : {vectorB.Value.Length} does NOT match number of columns in VectorA : {(vectorA.Value.Length / vectorA.Columns)}");
            }

            AcceleratorStream Stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, float, float>(ConsecutiveC_VDSP_2DKernel);

            var buffer = gpu.Allocate1D(new float[vectorA.Value.Length]); // Output
            var buffer2 = gpu.Allocate1D(vectorA.Value); // Input
            var buffer3 = gpu.Allocate1D(vectorB.Value); // Input


            kernelWithStream(Stream, buffer.IntExtent, buffer.View, buffer2.View, buffer3.View, vectorA.Columns, scalar);

            Stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(Stream);

            buffer.Dispose();
            buffer2.Dispose();
            buffer3.Dispose();

            Stream.Dispose();

            return new Vector(Output, vectorA.Columns);
        }
        // KERNEL
        static void ConsecutiveC_VDSP_2DKernel(Index1D index, ArrayView<float> Output, ArrayView<float> InputA, ArrayView<float> InputB, float Cols, float Scalar)
        {
            Output[index] = (InputA[index] / InputB[(int)(index / Cols)]) * Scalar;
        }


        public static Vector Diff(Accelerator gpu, Vector vectorA)
        {
            AcceleratorStream stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>>(DiffKernel);

            MemoryBuffer1D<float, Stride1D.Dense> buffer = gpu.Allocate1D(new float[vectorA.Value.Length - 1]); // Output
            MemoryBuffer1D<float, Stride1D.Dense> buffer2 = gpu.Allocate1D(vectorA.Value); //  Input

            kernelWithStream(stream, buffer.IntExtent, buffer.View, buffer2.View);

            stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(stream);

            buffer.Dispose();
            buffer2.Dispose();

            stream.Dispose();

            return new Vector(Output);
        }
        static void DiffKernel(Index1D index, ArrayView<float> Output, ArrayView<float> Input)
        {
            Output[index] = Input[index + 1] - Input[index];
        }

        public static Vector Diff_LogMult(Accelerator gpu, Vector vectorA, float scalar)
        {
            AcceleratorStream stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float>(DiffLogMultKernel);

            MemoryBuffer1D<float, Stride1D.Dense> buffer = gpu.Allocate1D(new float[vectorA.Value.Length - 1]); // Output
            MemoryBuffer1D<float, Stride1D.Dense> buffer2 = gpu.Allocate1D(vectorA.Value); //  Input

            kernelWithStream(stream, buffer.IntExtent, buffer.View, buffer2.View, scalar);

            stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(stream);

            buffer.Dispose();
            buffer2.Dispose();

            stream.Dispose();

            return new Vector(Output);
        }
        static void DiffLogMultKernel(Index1D index, ArrayView<float> Output, ArrayView<float> Input, float scalar)
        {
            Output[index] = (XMath.Log(Input[index + 1] - Input[index], XMath.E)) * scalar;
        }

        public static Vector Diff_Log(Accelerator gpu, Vector vectorA)
        {
            AcceleratorStream stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>>(DiffLogKernel);

            MemoryBuffer1D<float, Stride1D.Dense> buffer = gpu.Allocate1D(new float[vectorA.Value.Length - 1]); // Output
            MemoryBuffer1D<float, Stride1D.Dense> buffer2 = gpu.Allocate1D(vectorA.Value); //  Input

            kernelWithStream(stream, buffer.IntExtent, buffer.View, buffer2.View);

            stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(stream);

            buffer.Dispose();
            buffer2.Dispose();

            stream.Dispose();

            return new Vector(Output);
        }
        static void DiffLogKernel(Index1D index, ArrayView<float> Output, ArrayView<float> Input)
        {
            Output[index] = XMath.Log(Input[index + 1], XMath.E) - XMath.Log(Input[index], XMath.E);
        }

        public static Vector Power(Accelerator gpu, Vector vectorA, float pow)
        {
            if (pow == 2f)
            {
                return Vector.ConsecutiveOperation(gpu, vectorA, vectorA, '*');
            }

            AcceleratorStream stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float>(PowerKernel);

            MemoryBuffer1D<float, Stride1D.Dense> buffer = gpu.Allocate1D(new float[vectorA.Value.Length]); // Output
            MemoryBuffer1D<float, Stride1D.Dense> buffer2 = gpu.Allocate1D(vectorA.Value); //  Input

            kernelWithStream(stream, buffer.IntExtent, buffer.View, buffer2.View, pow);

            stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(stream);

            buffer.Dispose();
            buffer2.Dispose();

            stream.Dispose();

            return new Vector(Output);

        }
        static void PowerKernel(Index1D index, ArrayView<float> Output, ArrayView<float> Input, float pow)
        {
            Output[index] = XMath.Pow(Input[index],pow);
        }

        public static Vector Power(Accelerator gpu, Vector vectorA, bool AbsSecPow)
        {
            AcceleratorStream stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>>(PowerAbsKernel);

            MemoryBuffer1D<float, Stride1D.Dense> buffer = gpu.Allocate1D(new float[vectorA.Value.Length]); // Output
            MemoryBuffer1D<float, Stride1D.Dense> buffer2 = gpu.Allocate1D(vectorA.Value); //  Input

            kernelWithStream(stream, buffer.IntExtent, buffer.View, buffer2.View);

            stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(stream);

            buffer.Dispose();
            buffer2.Dispose();

            stream.Dispose();

            return new Vector(Output);
        }
        static void PowerAbsKernel(Index1D index, ArrayView<float> Output, ArrayView<float> Input)
        {
            Output[index] = Input[index] * XMath.Abs(Input[index]);
        }

        public static Vector MultiplySumAxZero(Accelerator gpu, Vector vectorA, Vector vectorB)
        {
            AcceleratorStream stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, int, int>(MultiplySumAxZeroKernel);

            MemoryBuffer1D<float, Stride1D.Dense>
                buffer = gpu.Allocate1D(new float[vectorA.Columns]), // Output
                buffer2 = gpu.Allocate1D(vectorA.Value), //  Input
                buffer3 = gpu.Allocate1D(vectorB.Value); //  Input


            kernelWithStream(stream, vectorA.Columns, buffer.View, buffer2.View, buffer3.View, vectorA.Columns, vectorA.Value.Length / vectorA.Columns );

            stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(stream);

            buffer.Dispose();
            buffer2.Dispose();

            stream.Dispose();

            return new Vector(Output, vectorA.Columns);
        }
        static void MultiplySumAxZeroKernel(Index1D index, ArrayView<float> Output, ArrayView<float> InputA, ArrayView<float> InputB, int columns, int rows)
        {
            for (int i = 0; i < rows; i++)
            {
                Output[index] += InputA[i * columns + index] * InputB[i * columns + index];
            }

        }

        public static Vector TenToPowerVector(Accelerator gpu, float[] vector)
        {
            AcceleratorStream stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>>(TenToPowerVectorKernel);

            MemoryBuffer1D<float, Stride1D.Dense> buffer = gpu.Allocate1D(vector); // IO

            kernelWithStream(stream, vector.Length, buffer.View);

            stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(stream);

            buffer.Dispose();

            stream.Dispose();

            return new Vector(Output, 1);

        }
        static void TenToPowerVectorKernel(Index1D index, ArrayView<float> IO)
        {
            IO[index] = XMath.Pow(10f, IO[index]);
        }

        public static Vector InvSqrt(Accelerator gpu, float[] vector)
        {
            AcceleratorStream stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>>(InvSqrtKernel);

            MemoryBuffer1D<float, Stride1D.Dense> buffer = gpu.Allocate1D(vector); // IO

            kernelWithStream(stream, buffer.IntExtent, buffer.View);

            stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(stream);

            buffer.Dispose();

            stream.Dispose();

            return new Vector(Output, 1);

        }
        static void InvSqrtKernel(Index1D index, ArrayView<float> IO)
        {
            IO[index] = XMath.Rsqrt(IO[index]);
        }



        // DOT PRODUCT : Vector dot Scalar, Vector dot Vector
        /// <summary>
        /// Calculates the Dot product between either
        /// | 1) Vector and Vector | 
        /// | 2) Vector and Scalar | 
        /// </summary>
        /// <param name="gpu"></param>
        /// <param name="vectorA"></param>
        /// <param name="param"></param>
        /// <returns>float</returns>
        public static float DotProduct(Accelerator gpu, Vector vectorA, Vector vectorB )
        {
            return ConsecutiveOperation(gpu, vectorA, vectorB, '*').Value.Sum();
        }
        public static float DotProduct(Accelerator gpu, Vector vectorA, float scalar)
        {
            return ScalarOperation(gpu, vectorA, scalar).Value.Sum();
        }


        // NORMALISE VECTOR
        /// <summary>
        /// Normalises an input Vector Class between 0 and 1
        /// </summary>
        /// <param name="gpu"></param>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Vector Normalise(Accelerator gpu, Vector vector, float Offset = 0f)
        {
            float Min = vector.Value.Min();
            float Max = vector.Value.Max();

            AcceleratorStream Stream = gpu.CreateStream();

            var kernelWithStream = gpu.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, float, float, float>(NormaliseKernel);

            var buffer = gpu.Allocate1D(new float[vector.Value.Length]);
            var buffer2 = gpu.Allocate1D(vector.Value);

            kernelWithStream(Stream, buffer.IntExtent, buffer.View, buffer2.View, Min, Max, Offset);

            Stream.Synchronize();

            float[] Output = buffer.GetAsArray1D(Stream);

            buffer.Dispose();
            buffer2.Dispose();

            Stream.Dispose();

            return new Vector(Output);
        }
        // KERNEL
        public static void NormaliseKernel(Index1D index, ArrayView<float> Output, ArrayView<float> Input, float Min, float Max, float Offset)
        {
            Output[index] = ((Input[index] - Min) / (Max - Min)) + Offset;
        }





    }




}
