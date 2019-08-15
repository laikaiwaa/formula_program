using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace Formulism
{
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
        /// <param name="param">函數基本參數</param>
        /// <returns></returns>
        public static double Integral(double a, double b, string classname, string methodname, object[] param)
        {
            double acc = 0.00001;
            long length = Convert.ToInt64(Math.Floor((b - a) / acc));

            Type t = Type.GetType(classname);
            MethodInfo minfo = t.GetMethod(methodname);

            double x = a;
            double sum = 0;

            for (int i = 1; i <= length; i++)
            {
                //sum =sum+ (double)minfo.Invoke(t, ObjJoin(x, param)) * acc;
                sum = sum + Distribution.Normal(x, 0, 1) * acc;
                x = x + acc;
            }
            return sum;
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
