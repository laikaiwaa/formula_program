using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Formulism.BaseMathMethod
{
    /// <summary>
    /// 通用矩陣運算
    /// </summary>
    public static class Matrix
    {
        //public static double[,] InverseQR(double[,] x)
        //{
        //    return new double[,];
        //}
        public static double[,] InverseSVD(double[,] x)
        {
            int m = x.GetLength(0);int n = x.GetLength(1);
            List<double[,]> usv = SVD(x);
            double[,] sigmain = new double[m,n];
            for(int i = 0; i < m; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    if (usv[1][i, j] > 0.000000001)
                    {
                        sigmain[i, j] = 1 / usv[1][i, j];
                    }
                    else { sigmain[i, j] = 0; }
                }
            }


            double[,] ainverfack = Product(Product(usv[2], Trans(sigmain)), Trans(usv[0]));
      
            return ainverfack;
        }
        /// <summary>
        /// 奇異值分解singular value decomposition 
        /// </summary>
        /// <param name="x">輸入矩陣</param>
        /// <returns></returns>
        public static List<double[,]> SVD(double[,] x)
        {
            int m = x.GetLength(0); int n = x.GetLength(1); 

            List<double[,]> svd = new List<double[,]>();
            List<double[,]> veigen = Eigen(Product(Trans(x), x), 10000);
            double[,] v = veigen[1];  ;
            List<double[,]> ueigen = Eigen(Product(x, Trans(x)), 10000);
            double[,] u = ueigen[1];  ;

            double[,] sigma = Product(Product(Trans(u), x), v);

              
            svd.Add(u);svd.Add(sigma);svd.Add(v);
            return svd;
        }
        /// <summary>
        /// 計算矩陣行列式
        /// </summary>
        /// <param name="x">輸入矩陣</param>
        /// <returns></returns>
        public static double Det(double[,] x)
        {
            int xm = x.GetLength(0);
            double detsum = 0;double detP = 0; double detN = 0;
            for (int m = 0; m < xm; m++)
            {
                double detPm = x[0,m]; double detNm = x[0, m];
                for (int k = 1; k < xm; k++)
                {
                    int newx = k;int newPy = m + k; int newNy = m - k;
                    if (newPy >= xm) { newPy = newPy - xm; }
                    if (newNy <0) { newNy = newNy + xm; }
                    detPm = detPm * x[newx, newPy]; detNm = detNm * x[newx, newNy];
                }
                detP = detP + detPm; detN = detN + detNm;
            }
            detsum = detsum + detP - detN;

            return detsum;
        } 
        /// <summary>
        /// 計算矩陣的特征向量及特征值（特征值分解）
        /// </summary>
        /// <param name="x">數據矩陣</param>
        /// <param name="ite">迭代次數</param>
        /// <returns></returns>
        public static List<double[,]> Eigen(double[,] x,int ite)
        {
            int dimesion = x.GetLength(0); int Npoint = x.GetLength(1);
            if (Det(x) == 0)
            {
                x = Plus(x, Lambda(0.000000001, Diag(EXP(new double[dimesion, Npoint]))));
            }
            double[,] Q = SchmidtOrthogonalization(x);
            double[,] bigQ = Q;
            double[,] xi = x;
            for(int i = 0; i < ite; i++)
            {
                xi = Product(Product(Trans(Q), xi), Q);
                if (Det(xi) ==0)
                {
                    xi = Plus(xi, Lambda(0.000000001, Diag(EXP(new double[dimesion, Npoint]))));
                }
                Q = SchmidtOrthogonalization(xi);
                bigQ = Product(bigQ, Q);
                //精確值
                double[,]err=Subtract(x, Product(bigQ, Product(xi,Trans(bigQ))) );
                double errr = 0;
                for(int m=0;m< dimesion; m++)
                {
                    for(int n = 0; n < dimesion; n++)
                    {
                        errr = errr +Math.Abs(err[m, n]);
                    }
                }
                if (errr < 0.000000001) { break; };
            }
            List<double[,]> eig = new List<double[,]> { xi, bigQ };
            return eig;
        }
        /// <summary>
        /// 施密特正交化
        /// </summary>
        /// <param name="x">輸入方陣</param>
        /// <returns></returns>
        public static double[,] SchmidtOrthogonalization(double[,] x)
        {
            int Dimesion = x.GetLength(0);
            List<double[,]> Vectorsp = new List<double[,]>();
            for (int i = 0; i < Dimesion; i++)
            {
                double[,] columnvec = new double[Dimesion, 1];
                for (int j = 0; j < Dimesion; j++)
                {
                    columnvec[j, 0] = x[j, i];
                }
                Vectorsp.Add(columnvec);
            }
            List<double[,]> Qlist = new List<double[,]>();
            for (int i = 0; i < Dimesion; i++)
            {
                double[,] Newcolum = Vectorsp[i];
                for (int j = 0; j < i; j++)
                {
                    Newcolum = Subtract(Newcolum, Projection(Vectorsp[i], Qlist[j]));
                }
                double templamda = Math.Sqrt(Product(Trans(Newcolum), Newcolum)[0, 0]);
                double[,] NewBasic = Lambda(1 / templamda, Newcolum);
                Qlist.Add(NewBasic);
            }
            double[,] Q = new double[Dimesion, Dimesion];
            for (int i = 0; i < Dimesion; i++)
            {
                for (int j = 0; j < Dimesion; j++)
                {
                    Q[j, i] = Qlist[i][j, 0];
                }
            }
            return Q;
        }
        /// <summary>
        /// 向量间的投影算法
        /// </summary>
        /// <param name="x">被投影向量</param>
        /// <param name="lay">地面向量</param>
        /// <returns></returns>
        public static double[,] Projection(double[,] x, double[,] lay)
        {
            int m = x.GetLength(0); int n = x.GetLength(1);
            double[,] shadow = new double[m,n];
            if (m==1)
            {
                shadow = Lambda(Product(x, Trans(lay))[0, 0] / Product(lay, Trans(lay))[0, 0], lay);
            }
            else
            {
                shadow = Lambda(Product(Trans(x), lay)[0, 0] / Product(Trans(lay), lay)[0, 0], lay);
            }
            return shadow;
        }
        /// <summary>
        /// 數值乘矩陣所有元素
        /// </summary>
        /// <param name="lambda">數值</param>
        /// <param name="a">輸入矩陣</param>
        /// <returns></returns>
        public static double[,] Lambda(double lambda,double[,] a)
        {
            int am = a.GetLength(0); int an = a.GetLength(1);
            double[,] c = new double[am, an];
            for (int i = 0; i < am; i++)
            {
                for (int j = 0; j < an; j++)
                {
                    c[i, j] = a[i, j]* lambda;
                }
            }
            return c;
        }
        /// <summary>
        /// 对角因素
        /// </summary>
        /// <param name="a">輸入矩陣</param>
        /// <returns></returns>
        public static double[,] Diag(double[,] a)
        {
            int am = a.GetLength(0);
            double[,] c = new double[am, am];
            for (int i = 0; i < am; i++)
            {
                for (int j = 0; j < am; j++)
                {
                    if (i == j)
                    {
                        c[i, j] = a[i, j];
                    }else
                    {
                        c[i, j] = 0;
                    }
                    
                }
            }
            return c;
        }
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
        /// <summary>
        /// 控制台輸出矩陣
        /// </summary>
        /// <param name="x">矩陣</param>
        public static void Print(double[,] x)
        {
            int m = x.GetLength(0); int n = x.GetLength(1);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(Math.Round(x[i, j], 4) + " ");
                }
                Console.WriteLine("");
            }
            Console.WriteLine("");
        }
    }
}
