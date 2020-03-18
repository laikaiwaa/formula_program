using Formulism.BaseMathMethod;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Formulism.MathematicalMethod
{
    /// <summary>
    /// 傅裡葉變換
    /// </summary>
    public static class DiscreteFourierTransform
    {
        /// <summary>
        /// 從時域變換到頻域的正向變換
        /// </summary>
        /// <param name="x">輸入時域的採樣點值</param>
        /// <returns></returns>
        public static List<double[,]> Transform(double[,] x)
        {
            List<double[,]> y = new List<double[,]>();
            int NLength = x.Length;
            double[,] YReal = new double[NLength, NLength];
            double[,] YImage = new double[NLength, NLength];
            for (int k = 0; k < NLength; k++)
            {
                for (int n = 0; n < NLength; n++)
                {
                    YReal[n, k] = Math.Cos(-n * 2 * k * Math.PI / NLength);
                    YImage[n, k] = Math.Sin(-n * 2 * k * Math.PI / NLength);
                }
            }
            YReal = Matrix.Product(x, YReal);
            YImage = Matrix.Product(x, YImage);
            y.Add(YReal); y.Add(YImage);
            return y;
        }
        /// <summary>
        /// 從頻域變換到時域的逆向變換
        /// </summary>
        /// <param name="y">輸入頻域的值</param>
        /// <returns></returns>
        public static List<double[,]> InverseTransform(List<double[,]> y)
        {
            List<double[,]> x = new List<double[,]>();
            if (y.Count == 1)
            {
                y.Add(new double[1, y[0].Length]);
            }
            int NLength = y[0].Length;
            //List<double[,]> x = new List<double[1, NLength]>;
            double[,] XReal = new double[NLength, NLength];
            double[,] XImage = new double[NLength, NLength];
            for (int k = 0; k < NLength; k++)
            {
                for (int n = 0; n < NLength; n++)
                {
                    XReal[k, n] = Math.Cos(n * 2 * k * Math.PI / NLength) / NLength;
                    XImage[k, n] = Math.Sin(n * 2 * k * Math.PI / NLength) / NLength;
                }
            }
            XReal = Matrix.Subtract(Matrix.Product(y[0], XReal), Matrix.Product(y[1], XImage));
            XImage = Matrix.Plus(Matrix.Product(y[1], XReal), Matrix.Product(y[0], XImage));
            x.Add(XReal); x.Add(XImage);
            return x;
        }
    }
}
