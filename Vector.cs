using ILGPU;
using ILGPU.Runtime;
using System;
using System.Linq;

namespace MachineLearningSpectralFittingCode
{
    public class Vector
    {
        public Vector(float[] value)
        {
            this.value = value;
        }

        public float[] value { get; set; }

        public static Vector ScalarProduct(AcceleratorId acceleratorId, Vector vector, float scalar)
        {
            using(var context = new Context())
            {
                using (var accelerator = Accelerator.Create(context, acceleratorId))
                {
                    var kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, float>(ScalarProductKernal);

                    var buffer  = accelerator.Allocate<float>(vector.value.Length);
                    var buffer2 = accelerator.Allocate<float>(vector.value.Length);
                    buffer.MemSetToZero();
                    buffer2.MemSetToZero();

                    // NEED TO ASSIGN vector.value to buffer2
                    buffer2.CopyFrom(vector.value, 0, 0, vector.value.Length);

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
            if (vectorA.value.Length != vectorB.value.Length)
            {
                throw new IndexOutOfRangeException("Vector A and Vector B provided MUST be of EQUAL length" );
            }


            using (var context = new Context())
            {
                using (var accelerator = Accelerator.Create(context, acceleratorId))
                {
                    var kernel = accelerator.LoadAutoGroupedStreamKernel<Index1, ArrayView<float>, ArrayView<float>, ArrayView<float>>(ConsecutiveProductKernal);

                    var buffer = accelerator.Allocate<float>(vectorA.value.Length); // Input
                    var buffer2 = accelerator.Allocate<float>(vectorA.value.Length); // Input
                    var buffer3 = accelerator.Allocate<float>(vectorA.value.Length); // Output
                    buffer.MemSetToZero();
                    buffer2.MemSetToZero();
                    buffer3.MemSetToZero();

                    // NEED TO ASSIGN vector.value to buffer2
                    buffer.CopyFrom(vectorA.value, 0, 0, vectorA.value.Length);
                    buffer2.CopyFrom(vectorB.value, 0, 0, vectorB.value.Length);

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
                return ConsecutiveProduct(acceleratorId, vectorA, (Vector)param).value.Sum();
            }
            else if (param.GetType() == typeof(float))
            {
                return ScalarProduct(acceleratorId, vectorA, (float)param).value.Sum();
            }
            else 
            {
                throw new Exception("Wrong input param object, accepts only T<Vector> and T<float>");
            }
        }



    }




}
