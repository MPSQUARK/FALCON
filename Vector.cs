using ILGPU;
using ILGPU.Runtime;
using System;
using System.Linq;

namespace MachineLearningSpectralFittingCode
{

    /* LOG :
     *      - Access Slice                          : WORKING
     *      - Access Value                          : WORKING
     *      - Consecutive Product                   : WORKING
     *      - Dot Product                           : WORKING
     *      - Fill                                  : WORKING
     *      - Scalar Compound Operation             : WORKING
     *      - Scalar Operation                      : WORKING
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
        public static Vector Fill(AcceleratorId acceleratorId, float Value, int Size, int Columns = 1)
        {
            using var context = new Context();
            using var accelerator = Accelerator.Create(context, acceleratorId);

            var kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, float>(FillKernel);

            var buffer = accelerator.Allocate<float>(Size); // Input
            buffer.MemSetToZero();

            kernel(buffer.Length, buffer.View, Value);
            accelerator.Synchronize();

            float[] Output = buffer.GetAsArray();
            buffer.Dispose();

            return new Vector(Output, Columns);
        }
        // KERNEL
        static void FillKernel(Index1 index, ArrayView<float> OutPut, float Value)
        {
            OutPut[index] = Value;
        }



        // Access 1 Value from 2D Vector
        public static float AccessVal(Vector vector, int row, int col)
        {
            return vector.Value[row*vector.Columns + col];
        }
        // Access a ROW or COLUMN from a 2D Vector
        public static Vector AccessSlice(AcceleratorId acceleratorId, Vector vector, int row_col_index , char row_col)
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

            using var context = new Context();
            using var accelerator = Accelerator.Create(context, acceleratorId);

            var kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<int>>(AccessSliceKernal);

            var buffer = accelerator.Allocate<float>(OutPutVectorLength);
            var buffer2 = accelerator.Allocate<float>(vector.Value.Length);
            var buffer3 = accelerator.Allocate<int>(5);

            buffer.MemSetToZero();
            buffer2.MemSetToZero();
            buffer3.MemSetToZero();

            buffer2.CopyFrom(vector.Value, 0, 0, vector.Value.Length);
            buffer3.CopyFrom(ChangeSelectLength, 0, 0, ChangeSelectLength.Length);

            kernel(OutPutVectorLength, buffer.View, buffer2.View, buffer3.View);

            accelerator.Synchronize();

            float[] Output = buffer.GetAsArray();

            buffer.Dispose();
            buffer2.Dispose();
            buffer3.Dispose();

            return new Vector(Output);
        }
        // KERNEL
        static void AccessSliceKernal(Index1 index, ArrayView<float> OutPut, ArrayView<float> Input, ArrayView<int> ChangeSelectLength)
        {
            OutPut[index] = Input[
                index * ChangeSelectLength[0] * ChangeSelectLength[4] + // iRcL
                index * ChangeSelectLength[1] +                         // iCc
                ChangeSelectLength[2] * ChangeSelectLength[4] +         // RsL
                ChangeSelectLength[3]];                                 // Cs
        }



        // SCALAR OPERATIONS : Vector * Scalar, Vector / Scalar, Vector +|- Scalar
        public static Vector ScalarOperation(AcceleratorId acceleratorId, Vector vector, float scalar, char operation = '*')
        {

            using var context = new Context();
            using var accelerator = Accelerator.Create(context, acceleratorId);
            
            Action<ILGPU.Index1, ILGPU.ArrayView<float>, ILGPU.ArrayView<float>, float> kernel;
            var buffer = accelerator.Allocate<float>(vector.Value.Length);
            var buffer2 = accelerator.Allocate<float>(vector.Value.Length);
            buffer.MemSetToZero();
            buffer2.MemSetToZero();


            switch (operation)
            {
                case '*':
                    kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, float>(ScalarProductKernal);

                    buffer2.CopyFrom(vector.Value, 0, 0, vector.Value.Length);
                    kernel(buffer.Length, buffer.View, buffer2.View, scalar);

                    break;
                case '/':
                    kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, float>(ScalarDivideKernal);

                    buffer2.CopyFrom(vector.Value, 0, 0, vector.Value.Length);
                    kernel(buffer.Length, buffer.View, buffer2.View, scalar);

                    break;
                case '+':
                    kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, float>(ScalarSumKernal);

                    buffer2.CopyFrom(vector.Value, 0, 0, vector.Value.Length);
                    kernel(buffer.Length, buffer.View, buffer2.View, scalar);

                    break;
            }

            accelerator.Synchronize();

            float[] Output = buffer.GetAsArray();

            buffer.Dispose();
            buffer2.Dispose();

            return new Vector(Output);
        }
        // KERNELS
        static void ScalarProductKernal(Index1 index, ArrayView<float> OutPut, ArrayView<float> Input, float Scalar )
        {
            OutPut[index] = Input[index] * Scalar;
        }
        static void ScalarDivideKernal(Index1 index, ArrayView<float> OutPut, ArrayView<float> Input, float Scalar)
        {
            OutPut[index] = Input[index] / Scalar;
        }
        static void ScalarSumKernal(Index1 index, ArrayView<float> OutPut, ArrayView<float> Input, float Scalar)
        {
            OutPut[index] = Input[index] + Scalar;
        }



        // COMPOUND SCALAR OPERATIONS : Vector * Scalar1 +|- Scaler2, Vector / Scalar1 +|- Scalar2
        public static Vector ScalarCompoundOperation(AcceleratorId acceleratorId, Vector vector, float Multiple, float Adder, string operation = "*+")
        {

            using var context = new Context();
            using var accelerator = Accelerator.Create(context, acceleratorId);
            
            Action<ILGPU.Index1, ILGPU.ArrayView<float>, ILGPU.ArrayView<float>, float, float> kernel;
            var buffer = accelerator.Allocate<float>(vector.Value.Length);
            var buffer2 = accelerator.Allocate<float>(vector.Value.Length);
            buffer.MemSetToZero();
            buffer2.MemSetToZero();


            switch (operation)
            {
                case "*+":
                    kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, float, float>(ScalarProductSumKernal);

                    buffer2.CopyFrom(vector.Value, 0, 0, vector.Value.Length);
                    kernel(buffer.Length, buffer.View, buffer2.View, Multiple, Adder);

                    break;
                case "/+":
                    kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, float, float>(ScalarDivideSumKernal);

                    buffer2.CopyFrom(vector.Value, 0, 0, vector.Value.Length);
                    kernel(buffer.Length, buffer.View, buffer2.View, Multiple, Adder);

                    break;
            }

            accelerator.Synchronize();

            float[] Output = buffer.GetAsArray();

            buffer.Dispose();
            buffer2.Dispose();

            return new Vector(Output);
        }
        // KERNELS
        static void ScalarProductSumKernal(Index1 index, ArrayView<float> OutPut, ArrayView<float> Input, float Multiplicator, float Adder)
        {
            OutPut[index] = Input[index] * Multiplicator + Adder;
        }
        static void ScalarDivideSumKernal(Index1 index, ArrayView<float> OutPut, ArrayView<float> Input, float Divisor, float Adder)
        {
            OutPut[index] = Input[index] / Divisor + Adder;
        }



        // Multiplies 2 Vectors Element by Element
        public static Vector ConsecutiveProduct(AcceleratorId acceleratorId, Vector vectorA, Vector vectorB)
        {
            if (vectorA.Value.Length != vectorB.Value.Length)
            {
                throw new IndexOutOfRangeException("Vector A and Vector B provided MUST be of EQUAL length" );
            }

            using var context = new Context();
            using var accelerator = Accelerator.Create(context, acceleratorId);
            
            var kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>>(ConsecutiveProductKernal);

            var buffer = accelerator.Allocate<float>(vectorA.Value.Length); // Input
            var buffer2 = accelerator.Allocate<float>(vectorA.Value.Length); // Input
            var buffer3 = accelerator.Allocate<float>(vectorA.Value.Length); // Output
            buffer.MemSetToZero();
            buffer2.MemSetToZero();
            buffer3.MemSetToZero();

            buffer.CopyFrom(vectorA.Value, 0, 0, vectorA.Value.Length);
            buffer2.CopyFrom(vectorB.Value, 0, 0, vectorB.Value.Length);

            kernel(buffer.Length, buffer.View, buffer2.View, buffer3.View);

            accelerator.Synchronize();

            float[] Output = buffer3.GetAsArray();

            buffer.Dispose();
            buffer2.Dispose();
            buffer3.Dispose();

            return new Vector(Output);
        }
        // KERNEL
        static void ConsecutiveProductKernal(Index1 index, ArrayView<float> InputA, ArrayView<float> InputB, ArrayView<float> OutPut )
        {

            OutPut[index] = InputA[index] * InputB[index];

        }



        // DOT PRODUCT : Vector dot Scalar, Vector dot Vector
        public static float DotProduct(AcceleratorId acceleratorId, Vector vectorA, object param )
        {
            if (param.GetType() == typeof(Vector))
            {
                return ConsecutiveProduct(acceleratorId, vectorA, (Vector)param).Value.Sum();
            }
            else if (param.GetType() == typeof(float))
            {
                return ScalarOperation(acceleratorId, vectorA, (float)param).Value.Sum();
            }
            else 
            {
                throw new Exception("Wrong input param object, accepts only T<Vector> and T<float>");
            }
        }





    }




}
