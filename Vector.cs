using ILGPU;
using ILGPU.Runtime;
using System;
using System.Linq;

namespace MachineLearningSpectralFittingCode
{
    public class Vector
    {
        public Vector(float[] value, int columns = 1)
        {
            this.Value = value;
            this.Columns = columns;
        }

        public float[] Value { get; set; } // Determines the values in the Vector
        public int Columns { get; set; } // Defines the number of Columns in a Vector STARTS AT 1

        public static float AccessVal(Vector vector, int row, int col)
        {
            return vector.Value[row*vector.Value.Length + col];
        }

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

            using (var context = new Context())
            {
                using (var accelerator = Accelerator.Create(context, acceleratorId))
                {
                    var kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<int>>(AccessSliceKernal);

                    var buffer = accelerator.Allocate<float>(OutPutVectorLength);
                    var buffer2 = accelerator.Allocate<float>(vector.Value.Length);
                    var buffer3 = accelerator.Allocate<int>(5);

                    buffer.MemSetToZero();
                    buffer2.MemSetToZero();
                    buffer3.MemSetToZero();

                    buffer2.CopyFrom(vector.Value, 0, 0, vector.Value.Length);
                    buffer3.CopyFrom(ChangeSelectLength , 0, 0, ChangeSelectLength.Length);

                    kernel(OutPutVectorLength, buffer.View, buffer2.View, buffer3.View);

                    accelerator.Synchronize();

                    float[] Output = buffer.GetAsArray();

                    buffer.Dispose();
                    buffer2.Dispose();
                    buffer3.Dispose();

                    return new Vector(Output);
                }
            }
        }

        static void AccessSliceKernal(Index1 index, ArrayView<float> OutPut, ArrayView<float> Input, ArrayView<int> ChangeSelectLength)
        {
            OutPut[index] = Input[
                index * ChangeSelectLength[0] * ChangeSelectLength[4] + // iRcL
                index * ChangeSelectLength[1] +                         // iCc
                ChangeSelectLength[2] * ChangeSelectLength[4] +         // RsL
                ChangeSelectLength[3]];                                 // Cs
        }

        public static Vector ScalarProduct(AcceleratorId acceleratorId, Vector vector, float scalar)
        {
            using(var context = new Context())
            {
                using (var accelerator = Accelerator.Create(context, acceleratorId))
                {
                    var kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, float>(ScalarProductKernal);

                    var buffer  = accelerator.Allocate<float>(vector.Value.Length);
                    var buffer2 = accelerator.Allocate<float>(vector.Value.Length);
                    buffer.MemSetToZero();
                    buffer2.MemSetToZero();

                    buffer2.CopyFrom(vector.Value, 0, 0, vector.Value.Length);

                    kernel(buffer.Length, buffer.View, buffer2.View , scalar);

                    accelerator.Synchronize();

                    float[] Output = buffer.GetAsArray();

                    buffer.Dispose();
                    buffer2.Dispose();

                    return new Vector(Output);
                }
            }
        }

        static void ScalarProductKernal(Index1 index, ArrayView<float> OutPut, ArrayView<float> Input, float Scalar )
        {

            OutPut[index] = Input[index] * Scalar;

        }

        public static Vector ConsecutiveProduct(AcceleratorId acceleratorId, Vector vectorA, Vector vectorB)
        {
            if (vectorA.Value.Length != vectorB.Value.Length)
            {
                throw new IndexOutOfRangeException("Vector A and Vector B provided MUST be of EQUAL length" );
            }


            using (var context = new Context())
            {
                using (var accelerator = Accelerator.Create(context, acceleratorId))
                {
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
            }
        }

        static void ConsecutiveProductKernal(Index1 index, ArrayView<float> InputA, ArrayView<float> InputB, ArrayView<float> OutPut )
        {

            OutPut[index] = InputA[index] * InputB[index];

        }

        public static float DotProduct(AcceleratorId acceleratorId, Vector vectorA, object param )
        {
            if (param.GetType() == typeof(Vector))
            {
                return ConsecutiveProduct(acceleratorId, vectorA, (Vector)param).Value.Sum();
            }
            else if (param.GetType() == typeof(float))
            {
                return ScalarProduct(acceleratorId, vectorA, (float)param).Value.Sum();
            }
            else 
            {
                throw new Exception("Wrong input param object, accepts only T<Vector> and T<float>");
            }
        }





    }




}
