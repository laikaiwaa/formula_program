using Formulism.MechineLearningMethod;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Formulism.BaseMathMethod
{
    /// <summary>
    /// 常用統計方法
    /// </summary>
    public static class Statistics
    {
        /// <summary>
        /// 離散數據的期望值及方差
        /// </summary>
        /// <param name="x">離散數組</param>
        /// <returns>
        /// <para>期望值</para>
        /// <para>樣本方差</para>
        /// <para>無偏估計樣本方差</para>
        /// </returns>
        public static List<double> DiscreteStatistic(double[,] x)
        {
            List<double> basestatistic = new List<double>();
            double ex=0;
            double var = 0;
            double var2 = 0;
            for(int i = 0; i < x.Length; i++)
            {
                ex = ex + x[0, i] / x.Length;
            }
            for (int i = 0; i < x.Length; i++)
            {
                var = var + (x[0, i]-ex)* (x[0, i] - ex) / x.Length;
                var2 = var2 + (x[0, i] - ex) * (x[0, i] - ex) / (x.Length-1);
            }
            basestatistic.Add(ex); basestatistic.Add(var); basestatistic.Add(var2);
            return basestatistic;
        } 
        /// <summary>
        /// 計算各階的自相關係數
        /// </summary>
        /// <param name="x">輸入序列，1*Point矩陣</param>
        /// <returns></returns>
        public static double[,] DiscreteAutocorrelation(double[,] x)
        {
            int length = x.Length;List<double> basestatic = DiscreteStatistic(x);
            double ex = basestatic[0];
            double var = basestatic[1];
            double[,] r = new double[1, length];
            for(int k = 0; k <= length; k++)
            {
                for(int i=0;i< length - k; i++)
                {
                    r[0, k] = r[0, k] + (x[0, i] - ex) * (x[0, i+k] - ex) / (var * length);
                } 
            }
            return r;
        }
        public static double[,] DiscretePartialAutocorrelation(double[,] x)
        {
            int length = x.Length; 
            double[,] r = new double[1, length];
            r[0, 0] = 1;r[0, 1] = DiscreteAutocorrelation(x)[0,1];
            for(int k = 2; k < length; k++)
            {
                double[,] Px = new double[length - k,k - 1];
                double[,] PyL = new double[length - k, 1];
                double[,] PyR = new double[length - k, 1];
                for (int i = 0; i < length - k; i++)
                {
                    for(int j = 0; j < k - 1; j++)
                    {
                        Px[i, j] = x[0,j +1+i];
                    }
                    PyL[i,0] = x[0, i];
                    PyR[i,0] = x[0,  i + k];

                }
                r[0,k]=
                Matrix.Product(
                Matrix.Trans(Matrix.Subtract(PyL, Matrix.Product(Px,OrdinaryLeastSquares.Beta(Px, PyL)))),
                Matrix.Subtract(PyR, Matrix.Product(Px, OrdinaryLeastSquares.Beta(Px, PyR)))
                )[0, 0] / Math.Sqrt(
                Matrix.Product(
                    Matrix.Product(
                    Matrix.Trans(Matrix.Subtract(PyL, Matrix.Product(Px, OrdinaryLeastSquares.Beta(Px, PyL)))),
                    Matrix.Subtract(PyL, Matrix.Product(Px, OrdinaryLeastSquares.Beta(Px, PyL)))
                    ),
                    Matrix.Product(
                    Matrix.Trans(Matrix.Subtract(PyR, Matrix.Product(Px, OrdinaryLeastSquares.Beta(Px, PyR)))),
                    Matrix.Subtract(PyR, Matrix.Product(Px, OrdinaryLeastSquares.Beta(Px, PyR)))
                    )
                )[0, 0]);
            }
            return r;
        }
        /// <summary>
        /// 計算序列的YM方法的偏自相關係數
        /// </summary>
        /// <param name="x">輸入序列，1*Point矩陣</param>
        /// <returns></returns>
        public static List<double[,]> DiscretePAbyYuleWaliker(double[,] x)
        {
            double[,] r = DiscreteAutocorrelation(x);
            int leng = r.Length;
            List<double[,]> AllPlevalPACSerise = new List<double[,]>();
            double[,] PlevalPartialAutoCorre = new double[1, leng];
            AllPlevalPACSerise.Add(PlevalPartialAutoCorre);
            for (int p = 1; p < leng; p++)
            {
                double[,] SingR = new double[1, p];
                for (int i = 0; i < p; i++)
                {
                    SingR[0, i] = r[0, i+1];
                }
                double[,] MatrixR = new double[p, p];
                for (int j = 0; j < p; j++)
                {
                    for(int k = 0; k < p; k++)
                    {
                        MatrixR[j, k] = r[0, Math.Abs(j-k)];
                    }
                }
                double[,] PPACTemp = Matrix.Product(Matrix.Inverse(MatrixR), Matrix.Trans(SingR));
                
                PlevalPartialAutoCorre[0,p]=PPACTemp[ p -1,0];
                
                AllPlevalPACSerise.Add(PPACTemp);
            }
            AllPlevalPACSerise[0] = PlevalPartialAutoCorre;


            return AllPlevalPACSerise;
        }
        }
}
