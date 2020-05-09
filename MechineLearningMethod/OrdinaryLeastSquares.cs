using Formulism.BaseMathMethod;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Formulism.MechineLearningMethod
{
    public static class OrdinaryLeastSquares
    {
        /// <summary>
        /// 計算線性回歸的參數值，返回Dim*1（參數維度數*1）的矩陣
        /// </summary>
        /// <param name="x">輸入Point*Dim（觀測點數*參數維度數）觀測值矩陣</param>
        /// <param name="y">真實值，Point*1（觀測點數*1）矩陣</param>
        /// <returns></returns>
        public static double[,] Beta(double[,] x, double[,] y)
        {
            double[,] beta  = Matrix.Product(
                Matrix.Product(
                    Matrix.Inverse(
                        Matrix.Product(
                            Matrix.Trans(x), x
                            )
                            ), Matrix.Trans(x)
                            ), y);
            return beta;
        }

    }
}
