using Formulism.BaseMathMethod;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Formulism.MechineLearningMethod
{
    /// <summary>
    /// 邏輯回歸
    /// </summary>
    public static class Logistic
    {

        /// <summary>
        /// 使用Pi函數
        /// </summary>
        /// <param name="x">輸入值</param> 
        public static double Pi(double x)
        {
            return Math.Exp(x) / (1 + Math.Exp(x));
        }
        /// <summary>
        /// 使用極大似然函數的一階導數
        /// </summary>
        /// <param name="x">X，維度n*k</param> 
        ///  <param name="y">Y,維度n*1,取值1或0</param> 
        ///   <param name="c">B,初始权值</param>
        public static double[,] DerivativeOne(double[,] x, double[,] y, double[,] c)
        {
            int am = x.GetLength(0), an = x.GetLength(1);
            double[,] B = c;
            double[,] Ja = Matrix.EXP(new double[am, 1]);
            double[,] Jb = Matrix.EXP(new double[1, an]);
            double[,] N = Matrix.EXP(new double[1, am]);

            double[,] FirstDerivate =
            Matrix.Subtract(
                Matrix.Product(
                    Matrix.Trans(y), x
                ),
                Matrix.Product(
                    N,
                    Matrix.HadamardProduct(
                            Matrix.Product(
                                Matrix.HadamardProduct(
                                    Matrix.Reciprocal(
                                        Matrix.Plus(
                                            Matrix.EXP(Matrix.Product(x, B)), Ja
                                        )
                                    ),
                                    Matrix.EXP(Matrix.Product(x, B))
                                ),
                                Jb
                            ),
                            x
                     )
                )
            )
            ;

            return FirstDerivate;
        }
        /// <summary>
        /// 使用極大似然函數的二階導數
        /// </summary>
        /// <param name="x">X，維度n*k</param> 
        ///  <param name="y">Y,維度n*1,取值1或0</param> 
        ///  <param name="B">B,初始權值</param>
        public static double[,] DerivativeTwo(double[,] x, double[,] y, double[,] b)
        {
            int am = x.GetLength(0), an = x.GetLength(1);
            double[,] B = b;
            double[,] Ja = Matrix.EXP(new double[am, 1]);
            double[,] Jb = Matrix.EXP(new double[1, an]);
            double[,] N = Matrix.EXP(new double[1, am]);
            double[,] Zero = new double[an, an];

            double[,] SecondDerivate =
            Matrix.Subtract(
                Zero,
                Matrix.Product(
                    Matrix.Trans(x),
                    Matrix.HadamardProduct(
                        Matrix.Product(
                            Matrix.HadamardProduct(
                                Matrix.Reciprocal(
                                    Matrix.HadamardProduct(
                                        Matrix.Plus(
                                            Matrix.EXP(Matrix.Product(x, B)), Ja
                                        ),
                                        Matrix.Plus(
                                            Matrix.EXP(Matrix.Product(x, B)), Ja
                                        )
                                    )
                                ),
                                Matrix.EXP(Matrix.Product(x, B))
                            ),
                            Jb
                        ),
                        x
                    )
                )
            )
            ;

            return SecondDerivate;
        }
        /// <summary>
        /// B的变动值
        /// </summary>
        /// <param name="inputx">inputx，維度n*k</param> 
        ///  <param name="inputy">inputy,維度n*1,取值1或0</param> 
        ///  <param name="b">B,初始權值</param>
        public static double[,] Db(double[,] inputx, double[,] inputy, double[,] b, out double[,] dfirst, out double[,] dsecond)
        {
            dfirst = DerivativeOne(inputx, inputy, b);
            dsecond = DerivativeTwo(inputx, inputy, b);

            double[,] Deltab =
Matrix.Subtract(
                    b
                    ,

                        Matrix.Product(

                                    Matrix.Inverse(
                                        dsecond
                                    )

                                ,
                                Matrix.Trans(dfirst)
                            )
                )
                ;
            return Deltab;
        }
        /// <summary>
        /// ObservedFisherInformation
        /// </summary>
        /// <param name="x">X，維度n*k</param> 
        ///  <param name="y">Y,維度n*1,取值1或0</param>  
        ///  <param name="b">B,維度n*1的權值</param>  
        public static double[,] ObservedFisherInformation(double[,] x, double[,] y, double[,] b)
        {
            int am = x.GetLength(0), an = x.GetLength(1);
            double[,] Result = new double[an, an];
            double[,] Zero = new double[an, an];
            Result = Matrix.Subtract(Zero, DerivativeTwo(x, y, b));
            return Result;
        }
        /// <summary>
        /// B參數的標準誤差
        /// </summary>
        /// <param name="x">X，維度n*k</param> 
        ///  <param name="y">Y,維度n*1,取值1或0</param> 
        ///  <param name="b">B,維度n*1的權值</param>
        public static double[,] StandardErrorMatrix(double[,] x, double[,] y, double[,] b)
        {
            int am = x.GetLength(0), an = x.GetLength(1);
            double[,] Result = new double[an, an];
            double[,] Zero = new double[an, an];
            Result = Matrix.Sqrt(Matrix.Inverse(ObservedFisherInformation(x, y, b)));
            return Result;
        }
        /// <summary>
        /// P-值
        /// </summary>
        /// <param name="x">自变量，維度n*k</param>
        /// <param name="y">因变量，維度n*1,取值1或0</param>
        /// <param name="b">B，維度n*1的估计權值</param>
        /// <returns></returns>
        public static double[,] PValue(double[,] x, double[,] y, double[,] b)
        {
            int leng = b.GetLength(0);
            double[,] P = new double[leng, 1];
            for (int i = 0; i < leng; i++)
            {
                P[i, 0] = Calculas.Integral(-100, b[i, 0] / StandardErrorMatrix(x, y, b)[i, i], "Formulism.BaseMathMethod.Function", "Normal",0.00001, new object[]{
                0,1});
                if (P[i, 0] > 0.5)
                {
                    P[i, 0] = (1 - P[i, 0]) * 2;
                }
                else
                {
                    P[i, 0] = P[i, 0] * 2;
                }
            }
            return P;
        }
    }
}
