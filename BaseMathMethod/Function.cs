using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Formulism.BaseMathMethod
{
    /// <summary>
    /// 通用分佈類型
    /// </summary>
    public static class Function
    {
        /// <summary>
        /// 單項式
        /// </summary>
        /// <param name="x">輸入值</param>
        /// <param name="n">次數</param>
        /// <param name="a">係數</param>
        /// <returns></returns>
        public static double SinglePolynomial(double x, int n, double a)
        {
            double pow = 1;
            for(int i = 0; i < n; i++)
            {
                pow = pow * x;
            }
            return pow;
        }
        /// <summary>
        /// 多元線性函數
        /// </summary>
        /// <param name="x">Npoint*Dimension 的輸入矩陣</param>
        /// <param name="beta">1*Dimension 的參數矩陣</param>
        /// <returns></returns>
        public static double[,] Linear(double[,] x, double[,] beta)
        {
            double[,] y = Matrix.Product(beta, Matrix.Trans(x));
            return y;
        }
        /// <summary>
        /// 正態分佈
        /// </summary>
        /// <param name="x">x,輸入值</param>
        /// <param name="e">e,平均值</param>
        /// <param name="sigema">sigema,方差</param>
        /// <returns></returns>
        public static double Normal(double x, double e, double sigema)
        {
            double fx = 1 / Math.Sqrt(2 * Math.PI) / sigema * Math.Exp(-(x - e) * (x - e) / 2 / (sigema * sigema)) ;
            return fx;
        }
        public static double Uniform(int Length)
        {
            double fx = 1 / Length;
            return fx;
        }
    }
}
