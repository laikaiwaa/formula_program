using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Formulism
{
    /// <summary>
    /// 通用矩陣運算
    /// </summary>
    public static class Matrix
    {
        /// <summary>
        /// 矩陣轉置
        /// </summary>
        /// <param name="a">輸入矩陣</param>
        /// <returns></returns>
        public static double[,] Trans(double[,] a)
        {
            int am = a.GetLength(0), an = a.GetLength(1);
            double[,] c = new double[an, am];
            for (int i = 0; i < am; i++)
            {
                for (int j = 0; j < an; j++)
                {
                    c[j, i] = a[i, j];
                }
            }
            return c;
        }
        /// <summary>
        /// 方陣的逆
        /// </summary>
        /// <param name="s">輸入矩陣</param>
        /// <returns></returns>
        public static double[,] Inverse(double[,] s)
        {
            double r = 0;
            double r2 = 0;

            int am = s.GetLength(0), an = s.GetLength(1);
            double[,] c = new double[am, an];

            double[,] a = new double[am, an];

            for (int i = 0; i < am; i++)
            {
                for (int j = 0; j < an; j++)
                {
                    a[j, i] = s[i, j];
                }
            }

            if (am == an)
            {
                //單位矩陣
                for (int i = 0; i < am; i++)
                {
                    for (int j = 0; j < an; j++)
                    {
                        if (i == j)
                        {
                            c[i, j] = 1;
                        }
                        else
                        {
                            c[i, j] = 0;
                        }
                    }
                }

                //向下運算

                for (int i = 0; i < am; i++)
                {
                    r = a[i, i];
                    for (int n = 0; n < an; n++)
                    {
                        a[i, n] = a[i, n] / r;

                        c[i, n] = c[i, n] / r;
                    }
                    for (int m = i + 1; m < am; m++)
                    {
                        r2 = a[m, i];
                        for (int n = 0; n < an; n++)
                        {
                            a[m, n] = a[m, n] - r2 * a[i, n];
                            c[m, n] = c[m, n] - r2 * c[i, n];
                        }
                    }
                }
                //向上運算
                for (int i = am - 1; i >= 0; i = i - 1)
                {
                    r = a[i, i];
                    for (int n = an - 1; n >= 0; n--)
                    {
                        a[i, n] = a[i, n] / r;

                        c[i, n] = c[i, n] / r;
                    }
                    for (int m = i - 1; m >= 0; m--)
                    {
                        r2 = a[m, i];
                        for (int n = an - 1; n >= 0; n--)
                        {
                            a[m, n] = a[m, n] - r2 * a[i, n];
                            c[m, n] = c[m, n] - r2 * c[i, n];
                        }
                    }
                }
                return c;
            }
            else
            {
                return c;
            }
        }
        /// <summary>
        /// 矩陣乘法
        /// </summary>
        /// <param name="a">輸入左矩陣</param>
        /// <param name="b">輸入右矩陣</param>
        /// <returns></returns>
        public static double[,] Product(double[,] a, double[,] b)
        {

            int am = a.GetLength(0), an = a.GetLength(1), bm = b.GetLength(0), bn = b.GetLength(1);
            double[,] c = new double[am, bn];
            if (an == bm)
            {
                for (int i = 0; i < am; i++)
                {
                    for (int j = 0; j < bn; j++)
                    {
                        c[i, j] = 0;
                        for (int k = 0; k < an; k++)
                        {
                            c[i, j] = a[i, k] * b[k, j] + c[i, j];
                        }
                    }
                }
                return c;
            }
            else
            {
                return c;
            }
        }
        /// <summary>
        /// 矩陣的Hadamard乘法
        /// </summary>
        /// <param name="a">輸入左矩陣</param>
        /// <param name="b">輸入右矩陣</param>
        /// <returns></returns>
        public static double[,] HadamardProduct(double[,] a, double[,] b)
        {
            int am = a.GetLength(0), an = a.GetLength(1), bm = b.GetLength(0), bn = b.GetLength(1);
            double[,] c = new double[am, bn];
            if (am == bm && an == bn)
            {
                for (int i = 0; i < am; i++)
                {
                    for (int j = 0; j < bn; j++)
                    {
                        c[i, j] = a[i, j] * b[i, j];
                    }
                }
                return c;
            }
            else
            {
                return c;
            }
        }
        /// <summary>
        /// 矩陣加法
        /// </summary>
        /// <param name="a">輸入左矩陣</param>
        /// <param name="b">輸入右矩陣</param>
        /// <returns></returns>
        public static double[,] Plus(double[,] a, double[,] b)
        {
            int am = a.GetLength(0), an = a.GetLength(1);
            double[,] c = new double[am, an];

            for (int i = 0; i < am; i++)
            {
                for (int j = 0; j < an; j++)
                {
                    c[i, j] = a[i, j] + b[i, j];
                }
            }
            return c;
        }
        /// <summary>
        /// 矩陣減法
        /// </summary>
        /// <param name="a">輸入左矩陣</param>
        /// <param name="b">輸入右矩陣</param>
        /// <returns></returns>
        public static double[,] Subtract(double[,] a, double[,] b)
        {
            int am = a.GetLength(0), an = a.GetLength(1), bm = b.GetLength(0), bn = b.GetLength(1);
            double[,] c = new double[am, bn];
            if (am == bm && an == bn)
            {
                for (int i = 0; i < am; i++)
                {
                    for (int j = 0; j < bn; j++)
                    {
                        c[i, j] = a[i, j] - b[i, j];
                    }
                }
                return c;
            }
            else
            {
                return c;
            }
        }
        /// <summary>
        /// 指數矩陣
        /// </summary>
        /// <param name="a">輸入矩陣</param>
        /// <returns></returns>
        public static double[,] EXP(double[,] a)
        {
            int am = a.GetLength(0), an = a.GetLength(1);
            double[,] c = new double[am, an];

            for (int i = 0; i < am; i++)
            {
                for (int j = 0; j < an; j++)
                {
                    c[i, j] = Math.Exp(a[i, j]);
                }
            }
            return c;
        }
        /// <summary>
        /// 矩陣元素倒數
        /// </summary>
        /// <param name="a">輸入矩陣</param>
        /// <returns></returns>
        public static double[,] Reciprocal(double[,] a)
        {
            int am = a.GetLength(0), an = a.GetLength(1);
            double[,] c = new double[am, an];

            for (int i = 0; i < am; i++)
            {
                for (int j = 0; j < an; j++)
                {
                    if (a[i, j] == 0)
                    {
                        c[i, j] = 0;
                    }
                    else
                    {
                        c[i, j] = 1.0 / (a[i, j]);

                    }

                }
            }
            return c;
        }
        /// <summary>
        /// 矩陣元素平方根
        /// </summary>
        /// <param name="a">輸入矩陣</param>
        /// <returns></returns>
        public static double[,] Sqrt(double[,] a)
        {
            int am = a.GetLength(0), an = a.GetLength(1);
            double[,] c = new double[am, an];
            for (int i = 0; i < am; i++)
            {
                for (int j = 0; j < an; j++)
                {
                    c[i, j] = Math.Sqrt(a[i, j]);
                }
            }
            return c;
        } 
    }
}
