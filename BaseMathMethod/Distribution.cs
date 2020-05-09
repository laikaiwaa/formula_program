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
    public static class Distribution
    {
        /// <summary>
        /// 正態分佈
        /// </summary>
        /// <param name="x">x,輸入值</param>
        /// <param name="e">e,平均值</param>
        /// <param name="sigema">sigema,方差</param>
        /// <returns></returns>
        public static double Normal(double x, double e, double sigema)
        {
            double fx =
            1 / Math.Sqrt(2 * Math.PI) / sigema * Math.Exp(-(x - e) * (x - e) / 2 / (sigema * sigema));
            return fx;
        }
        public static double Uniform(int Length)
        {
            double fx = 1 / Length;
            return fx;
        }
    }
}
