using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;


namespace Formulism.BaseMathMethod
{
    //test
    /// <summary>
    /// 通用微積分運算
    /// </summary>
    public static class Calculas
    {
        public delegate double Function(double x, double e, double sigema);
        /// <summary>
        /// 積分
        /// </summary>
        /// <param name="a">初始x</param>
        /// <param name="b">結束x</param>
        /// <param name="classname">函數所屬類的名稱</param>
        /// <param name="methodname">函數名稱</param>
        /// <param name="acc">積分精確度</param>
        /// <param name="param">函數基本參數</param>
        /// <returns></returns>
        public static double Integral(double a, double b, string classname, string methodname,double acc, object[] param)
        {

            long length = Convert.ToInt64(Math.Floor((b - a) / acc));

            Type t = Type.GetType(classname);
            MethodInfo minfo = t.GetMethod(methodname);

            double x = a;
            double sum = 0;

            for (int i = 1; i <= length; i++)
            {
                //sum =sum+ (double)minfo.Invoke(t, ObjJoin(x, param)) * acc;
                
                x = x + acc;
            }
            return sum;
        }
        public static double[,] JacobiMatrix(double[,] InitialthetaPosition,double[,] x, double Detalx, string Classname, string Methodname, object[] Param)
        {
            Type t = Type.GetType(Classname);
            MethodInfo minfo = t.GetMethod(Methodname);
            int mdimension = x.GetLength(0);
            int ndimension = InitialthetaPosition.GetLength(1);
            double[,] Jacobi = new double[mdimension, ndimension];
            for (int i = 0; i < ndimension; i++)
            {
                double[,] FuctionXi = (double[,])minfo.Invoke(t, new object[] { x, InitialthetaPosition });
                InitialthetaPosition[0, i] = InitialthetaPosition[0, i] + Detalx;
                double[,] FuctionXideltay = (double[,])minfo.Invoke(t, new object[] { x, InitialthetaPosition });
                double[,] FuctionXidy = Matrix.Subtract(FuctionXideltay, FuctionXi);
                InitialthetaPosition[0, i] = InitialthetaPosition[0, i] - Detalx;
                for (int j = 0; j < mdimension; j++)
                {
                    Jacobi[j, i] = FuctionXidy[0,j] / Detalx;
                }
            }
            return Jacobi;
        }
        public static double[,]  Derivative(double Point, double Detalx,int Rank, string Classname, string Methodname, object[] Param)
        {
            Type t = Type.GetType(Classname);
            MethodInfo minfo = t.GetMethod(Methodname);
            double[,] y = new double[1, Rank + 1];
            double[,] derativ = new double[Rank+1, Rank+1 ];
            for (int i = 0; i < Rank+1; i++)
            {
                y[0, i] = (double)minfo.Invoke(t, ObjJoin(Point+ Detalx*i, Param));
                derativ[0, i] = y[0, i];
            } 
            for(int j = 1; j < Rank+1; j++)
            { 
                for (int k = 0; k < Rank+1-j; k++)
                {
                    derativ[j, k] =( derativ[j - 1, k + 1]- derativ[j - 1, k ])/ Detalx;
                    if (derativ[j, k] < Detalx)
                    {
                        //derativ[j, k] = Detalx;
                    }
                }
            }
            
            
            return derativ;
        }
        /// <summary>
        /// 階乘運算
        /// </summary>
        /// <param name="n">輸入階數</param>
        /// <returns></returns>
        public static double Factorial(int n)
        {
            double Fact = 1;
            for(int i = 1; i < n+1;i++)
            {
                Fact = Fact * i;
            }
            return Fact;
        }

        /// <summary>
        /// 添加對象到object[]前邊
        /// </summary>
        /// <param name="a">需要添加的值</param>
        /// <param name="b">被添加的object[]</param>
        /// <returns></returns>
        public static object[] ObjJoin(double a, object[] b)
        {
            int l = b.Length;
            object[] n = new object[l + 1];
            for (int i = 0; i < l; i++)
            {
                n[i + 1] = b[i];
            }
            n[0] = a;
            return n;
        }

    }
}
