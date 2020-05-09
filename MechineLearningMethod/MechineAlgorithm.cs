using Formulism.BaseMathMethod;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using Formulism.MechineLearningMethod;

namespace Formulism.MechineLearningMethod
{
    public static class MechineAlgorithm
    {
        /// <summary>
        /// adaboost
        /// </summary>
        /// <param name="x">Dimension*Npoint。輸入矩陣</param>
        /// <param name="y">1*Npoint。標籤矩陣</param>
        /// <returns></returns>
        public static double[,] Adaboost(double[,] x, double[,] y)
        {
            int Dimesion = x.GetLength(0);int Npoint = x.GetLength(1);
            double[,] initialweight =Matrix.Lambda(1.0/Npoint, Matrix.EXP( new double[1, Npoint]));
            //分类器数目
            double[,] weeknum = new double[1, Dimesion];
            int modusize = 0;
            for (int i = 0; i < Dimesion; i++)
            {
                weeknum[0,i] = 5;
                modusize = modusize + Convert.ToInt32(weeknum[0, i]);
            }
            double[,] weekmodel= new double[1, modusize];
            //错误率
            double[,] e = new double[1, modusize]; double[,] modelweigh = new double[1, modusize];
            for (int k = 0; k < Dimesion; k++)
            {
                for (int i = 0; i < Convert.ToInt32(weeknum[0, k]); i++)
                {
                    int disize =Convert.ToInt32(weeknum[0, k]);
                    e[0, i+ k * disize] = 0; weekmodel[0, i + k*disize] = 0.5+i;
                    for (int j = 0; j < Npoint; j++)
                    {
                        double morigh = (x[k, j] - weekmodel[0, i + k * disize]);
                        if(morigh>=0) { morigh = 1; } else { morigh = -1; }
                        if (morigh * y[0, j] < 0) { e[0, i + k * disize] = e[0, i + k * disize] + initialweight[0, j]; }
                    }
                    //模型權值
                    modelweigh[0, i + k * disize] = 0.5 * Math.Log((1 - e[0, i + k * disize]) / e[0, i + k * disize]);
                    //更新數據點權值
                    for (int j = 0; j < Npoint; j++)
                    {
                        double morigh = (x[k, j] - weekmodel[0, i + k * disize]);
                        if (morigh >= 0) { morigh = 1; } else { morigh = 0; }
                        initialweight[0, j] = initialweight[0, j]*Math.Exp(-1* morigh * modelweigh[0, i + k * disize]*y[0, j]);
                    }
                    //數據點權值正規化
                    initialweight = Matrix.Lambda(1 / (Matrix.Product(initialweight, Matrix.Trans(Matrix.EXP( new double[1, Npoint])))[0, 0]), initialweight);
                }
            }
            //強分類器
            double[,] hx = new double[1, Npoint];
            for (int j = 0; j < Npoint; j++)
            {
               
                for (int m = 0; m < Dimesion; m++)
                {
                    int disize2 = Convert.ToInt32(weeknum[0, m]);
                    
                    for (int n=0;n< disize2; n++)
                    {
                        double morigh = (x[m, j] - weekmodel[0,n + m * disize2]);
                        if (morigh >= 0) { morigh = 1; } else { morigh = -1; }
                        hx[0, j] = hx[0, j] + morigh * modelweigh[0, n + m * disize2];
                    }
                }
                if (hx[0, j] >= 0) { hx[0, j] = 1; } else { hx[0, j] = -1; }
            }

             return e;
        }
        /// <summary>
        /// 支持向量機
        /// </summary>
        /// <param name="x">Dimension*Npoint，輸入矩陣</param>
        /// <param name="y">1*Npoint，分類矩陣</param>
        /// <returns></returns>
        public static List<double[,]> SupportVectorMechines(double[,] x, double[,] y)
        {
            int Dimension = x.GetLength(0); int Npoint = x.GetLength(1);
            double[,] Lambda = new double[1, Npoint];
            double[,] a = new double[Npoint, 1];

           
            
            for(int jj = 0; jj < 100; jj++)
            {

                double[,] wt = Matrix.Lambda(-1, Matrix.Product(Matrix.HadamardProduct(x, Matrix.Product(Matrix.EXP(new double[Dimension, 1]), y)), a));
                double[,] bt = Matrix.Subtract(Matrix.HadamardProduct(Matrix.EXP(new double[Npoint, 1]), Matrix.Reciprocal(Matrix.Trans(y))), Matrix.Product(Matrix.Trans(x), wt));

                double[,] ywxb = Matrix.HadamardProduct(y, Matrix.Plus(Matrix.Product(Matrix.Trans(wt), x), Matrix.Trans(bt)));
                int bih = 0;
                
                for (int k = 0; k < Npoint; k++)
                {
                    if((ywxb[0,k]>1 & a[k, 0] != 0)||(ywxb[0, k] ==1 &( a[k, 0] >0 ))) { bih = k; }
                }



                for (int asi=bih;asi< Npoint; asi++)
                {
                    double[,] xx = Matrix.Product(Matrix.Trans(x), x);

                    double[,] tempayx = Matrix.Product(
                    Matrix.Product(
                    Matrix.HadamardProduct(Matrix.Trans(a), y),
                    Matrix.Trans(x)
                    ), x);
                    int asib = asi + 1; if (asib >= Npoint) {
                        double dob = Matrix.Product(y, a)[0, 0];
                        a[asi, 0] = (a[asi, 0] * y[ 0, asi] - dob) / y[0, asi];
                        //break;
                        asib = asib- Npoint;
                    }
                    double constri = a[asi, 0] * y[0, asi] + a[asib, 0] * y[0, asib];
                    double temp = tempayx[0, asi] - tempayx[0, asib];
                    double k = 2 * xx[asi, asib] - xx[asib, asib] - xx[asi, asi];

                    


                    a[asi, 0] = (1 - y[0, asib] / y[0, asi] - y[0, asib] * temp + a[asib, 0] * k) / k;



                    a[asib, 0] = (constri - a[asi, 0] * y[0, asi]) / y[0, asib];
                    //asi = asi + 1;
                    
        




                }

                double min=0;
                for (int asi2 = 0; asi2 < Npoint   ; asi2++)
                {
                    if (a[asi2, 0] < -1) {  a[asi2, 0]=-1; }
                    
                    if (a[asi2, 0] > 0) { a[asi2, 0] = 0;} 
                 }
                //for (int asi3 = 0; asi3 < Npoint; asi3++)
                //{
                //    if (a[asi3, 0] !=min) { a[asi3, 0] = 0; }
                //}
            }
            double[,] w =
                Matrix.Lambda(
                    -1,
                    Matrix.Product(
                        Matrix.HadamardProduct(
                            x,
                            Matrix.Product(
                                Matrix.EXP(new double[Dimension, 1]),
                                y
                            )
                        ),
                        a
                    )
                );
            w = Matrix.Lambda(1 / Math.Sqrt(Matrix.Product(Matrix.Trans(w), w)[0, 0]),w);
            double[,]  b =
                Matrix.Subtract(
                    Matrix.HadamardProduct(
                    Matrix.EXP(new double[Npoint, 1]),
                    Matrix.Reciprocal(Matrix.Trans(y))
                    ),
                    Matrix.Product(Matrix.Trans(x), w)
                    );
            double bav = Matrix.Product(Matrix.Trans( a), b)[0,0]/Matrix.Product(Matrix.Trans(a), Matrix.EXP(new double[Npoint, 1]))[0,0];
            double[,] outbav = { { bav } };
            List<double[,]> module=new List<double[,]>();module.Add(w);module.Add(outbav);
            return module;
        }
       
        /// <summary>
        /// 主成分分析法
        /// </summary>
        /// <param name="x">維度*觀測數，輸入矩陣</param>
        /// <param name="iterator">迭代次數</param>
        /// <returns></returns>
        public static double[,] PrincipalComponentAnalysis(double[,] x,int iterator)
        {
            double Npoint = x.GetLength(1);int Dimension = x.GetLength(0);
            double[,] CovX= Matrix.Lambda(1/Npoint ,Matrix.Product(x, Matrix.Trans(x)) );
            //if (Matrix.Det(CovX) == 0)
            //{
            //    CovX=Matrix.Plus(
            //        CovX, Matrix.Diag(
            //           Matrix.Lambda(0.0000001, Matrix.EXP(new double[Dimension, Dimension]))
            //            )
            //        )
            //        ;
            //}
            double[,] eigenvector = Matrix.Eigen(CovX, iterator)[1];
            double[,] NewCom = Matrix.Product(Matrix.Trans(x), eigenvector);
            return NewCom;
        }
        /// <summary>
        /// 高斯混合模型
        /// </summary>
        /// <param name="x">输入矩阵</param>
        /// <param name="ModelNum">模型数量</param>
        /// <param name="Iterator">迭代次数</param>
        /// <returns></returns>
        public static List<double[,]> GMM(double[,] x, int ModelNum,int Iterator)
        {
            //初始化參數
            List<double> av = Statistics.DiscreteStatistic(x);
            double[] GMMtheta = { av[0], av[1] };
            int Npoint = x.GetLength(1);
            List<double[,]> NewParam = new List<double[,]>();
            List<object[]> ParamList = new List<object[]>(Npoint);
            for (int md = 0; md < ModelNum; md++)
            {
                object[] TempParam = new object[2];
                for (int nd = 0; nd < 2; nd++)
                {
                    TempParam[nd] =GMMtheta[nd] + 0.1 * md;
                }
                ParamList.Add(TempParam);
            }
            
            //参数更新
            double[,] mu = new double[1, ModelNum];
            double[,] sigma = new double[1, ModelNum];
            double[,] SumKmodelweight = new double[1, ModelNum];
            double[,] InitalModelWeight = new double[1, ModelNum];
            InitalModelWeight = Matrix.Lambda(1 / (double)ModelNum, Matrix.EXP(InitalModelWeight));
            double[,] ProbyXTheta =Matrix.Product(Matrix.EXP(new double[Npoint,1]) ,InitalModelWeight);

            for (int iter=0; iter< Iterator; iter++)
            {
                ProbyXTheta = ExpectationMaximization(x, ModelNum, ProbyXTheta, "Formulism.BaseMathMethod.Function", "Normal", ParamList);

                for (int i = 0; i < ModelNum; i++)
                {
                    SumKmodelweight[0, i] = 0;
                    for (int j = 0; j < Npoint; j++)
                    {
                        SumKmodelweight[0, i] = SumKmodelweight[0, i] + ProbyXTheta[j, i];
                    }
                    if(SumKmodelweight[0, i] != 0) {
                        mu[0, i] = Matrix.Product(x, ProbyXTheta)[0, i] / SumKmodelweight[0, i];
                        sigma[0, i] =(double) ParamList[i][1];
                        double sigmacalting = 0;
                        for (int j2 = 0; j2 < Npoint; j2++)
                        {
                            sigmacalting = sigmacalting + ProbyXTheta[j2, i] * (x[0, j2] - mu[0, i]) * (x[0, j2] - mu[0, i]);
                        }
                        if (sigmacalting != 0)
                        {
                            sigma[0, i] = sigmacalting / SumKmodelweight[0, i];
                        }
                    }
                    ParamList[i] =new object[] {  mu[0, i], sigma[0, i] } ;
                }
            }
            foreach(object[] p in ParamList)
            {
                NewParam.Add(new double[,] { { (double)p[0], (double)p[1] } });
            }
            //参数输出
            return NewParam;
        }
        /// <summary>
        /// 最大期望法
        /// </summary>
        /// <param name="x">輸入矩陣</param>
        /// <param name="ModelNum">模型數量</param>
        /// <param name="ModelWeight">N*模型數量，各點在各模型的權值</param>
        /// <param name="Classname">調用函數的類名</param>
        /// <param name="Methodname">調用函數的方法名</param>
        /// <param name="ParamList">模型數量*參數個數，調用方法的參數</param>
        /// <returns></returns>
        public static double[,] ExpectationMaximization(double[,] x,int ModelNum,double[,] ModelWeight, string Classname, string Methodname, List<object[]> ParamList)
        {
            //初始化參數
            int Npoint = x.GetLength(1);
            int NofParam = ParamList[0].Length;
            
            
            //調用結果
            Type t = Type.GetType(Classname);
            MethodInfo minfo = t.GetMethod(Methodname);

            double[,] ProXbyTheta = new double[Npoint, ModelNum];
            for(int i = 0; i < Npoint; i++)
            {
                for(int j = 0; j < ModelNum; j++)
                {
                    object[] ijParam = new object[NofParam+1];
                    ijParam[0] = x[0, i];
                    for (int k=0;k< NofParam;k++)
                    {
                        ijParam[k + 1] = ParamList[j][k];
                    }
                    ProXbyTheta[i, j] = (double)minfo.Invoke(t, ijParam);
                }
            }
            //計算新的各個模型的權重
            double[,] SumModelPro = Matrix.Product(
                Matrix.HadamardProduct(ProXbyTheta, ModelWeight)
                , Matrix.Trans(Matrix.EXP(new double[1, ModelNum]))
                );
            double[,] ProbyXTheta = new double[Npoint, ModelNum];
            for (int i = 0; i < Npoint; i++)
            {
                for(int j = 0; j < ModelNum; j++)
                {
                    if (SumModelPro[i, 0]==0)
                    {
                        ProbyXTheta[i, j] = 0;
                    }
                    else
                    {
                        ProbyXTheta[i, j] = ProXbyTheta[i, j] * ModelWeight[i, j] / SumModelPro[i, 0];
                    }
                }
            }
            return ProbyXTheta;
        }
        /// <summary>
        /// levenbergMarquardt方法
        /// </summary>
        /// <param name="InitialthetaPosition">參數初始值</param>
        /// <param name="x">輸入值矩陣</param>
        /// <param name="y">輸出值矩陣</param>
        /// <param name="Lambda">LM因子</param>
        /// <param name="Detalx">精確度</param>
        /// <param name="Classname">調用函數的方法名</param>
        /// <param name="Methodname">調用函數名</param>
        /// <param name="Param">調用函數的參數</param>
        /// <returns></returns>
        public static double[,] LevenbergMarquardt(double[,] InitialthetaPosition, double[,] x, double[,] y, double Lambda, double Detalx, string Classname, string Methodname, object[] Param)
        {
            Type t = Type.GetType(Classname);
            MethodInfo minfo = t.GetMethod(Methodname);
            double[,] Yestimate = (double[,])minfo.Invoke(t, new object[] { x, InitialthetaPosition });
            double[,] Jacobi = Calculas.JacobiMatrix(InitialthetaPosition, x, Detalx, Classname, Methodname, Param);
            double[,] DeltaTheata =
                Matrix.Product(
                Matrix.Inverse(
                    Matrix.Plus(Matrix.Product(Matrix.Trans(Jacobi), Jacobi),
                        Matrix.Lambda(Lambda,Matrix.Diag(Matrix.Product(Matrix.Trans(Jacobi), Jacobi)))
                        )
                    ),
                Matrix.Product(Matrix.Trans(Jacobi), Matrix.Subtract(y, Matrix.Trans(Yestimate))));
            double[,] NewThetaPosition = InitialthetaPosition;
            for (int i = 0; i < InitialthetaPosition.GetLength(1); i++)
            {
                NewThetaPosition[0, i] = InitialthetaPosition[0, i] + DeltaTheata[i, 0] * 0.1;
            }
            return NewThetaPosition;
        }
        /// <summary>
            /// 梯度求解
            /// </summary>
            /// <param name="InitialPosition">1*維度，初始點</param>
            /// <param name="x">輸入值矩陣</param>
            /// <param name="y">輸出值矩陣</param>
            /// <param name="ClassName">調用函數所屬的類名</param>
            /// <param name="MethodName">調用函數所屬的方法名</param>
            /// <param name="GradsStepsize">梯度測算步長</param>
            /// <param name="Stepsize">步進長度</param>
            /// <returns></returns>
        public static List<double[,]> GradientDescent(double[,] InitialPosition, double[,] x, double[,] y, string ClassName, string MethodName,double GradsStepsize, double Stepsize)
        {
            object[] c = { x, InitialPosition };
            Type t = Type.GetType(ClassName);
            MethodInfo minfo = t.GetMethod(MethodName);
            double[,] erroryx0 = Matrix.Subtract((double[,])minfo.Invoke(t, c),y);
            double[,] fxo = Matrix.Product(erroryx0, Matrix.Trans(erroryx0));
            int Dimensionofgrad = InitialPosition.Length;

            double[,] grad = new double[1, Dimensionofgrad];
            for (int i = 0; i < Dimensionofgrad; i++)
            {
                InitialPosition[0, i] = InitialPosition[0, i] + GradsStepsize;
                object[] cgrad = { x, InitialPosition };
                double[,] errorypargrad = Matrix.Subtract((double[,])minfo.Invoke(t, c), y);
                double[,] errorypargradsum = Matrix.Product(erroryx0, Matrix.Trans(Matrix.Subtract(errorypargrad, erroryx0)));
                grad[0, i] = errorypargradsum[0, 0] / GradsStepsize;
                InitialPosition[0, i] = InitialPosition[0, i] - GradsStepsize;
            }
            double sumsqrtgrad = 0;

            for (int j = 0; j < Dimensionofgrad; j++)
            {
                sumsqrtgrad = sumsqrtgrad + grad[0, j] * grad[0, j];
            }
            for (int j = 0; j < Dimensionofgrad; j++)
            {
                grad[0, j] = grad[0, j] / Math.Sqrt(sumsqrtgrad) * Stepsize;
            }
            //sepcial  grad method
            //object[] Partc = { x, InitialPosition };
            //double[,] errorydx = Matrix.Subtract((double[,])minfo.Invoke(t, Partc), y);
            //double[,] grad = Matrix.Product(errorydx,x);

            //double[,] stepmatrix = new double[1, Dimensionofgrad];
            //for (int i = 0; i < Dimensionofgrad; i++)
            //{
            //    stepmatrix[0, i] = Stepsize;
            //}

            //double[,] newposition = Matrix.Subtract(InitialPosition,Matrix.HadamardProduct( grad, stepmatrix));
            double[,] newposition = Matrix.Subtract(InitialPosition, grad);
            List<double[,]> result = new List<double[,]>();
            result.Add(fxo); result.Add(grad); result.Add(newposition);
            return result;
        }
        /// <summary>
        /// 步長自適應梯度下降法
        /// </summary>
        /// <param name="InitialPosition">1*維度，初始點</param>
        /// <param name="x">輸入值矩陣</param>
        /// <param name="y">輸出值矩陣</param>
        /// <param name="ClassName">調用函數所屬的類名</param>
        /// <param name="MethodName">調用函數所屬的方法名</param>
        /// <param name="GradsStepsize">梯度測算步長</param>
        /// <param name="Stepsize">步進長度</param>
        /// <param name="Iterator">迭代次數</param>
        /// <param name="Accurace">要求準確度</param>
        /// <param name="Printmodel">各步輸出模式</param>
        /// <returns></returns>
        public static List<double[,]> FastGradientDescent(double[,] InitialPosition, double[,] x, double[,] y, string ClassName, string MethodName, double GradsStepsize, double Stepsize,int Iterator,double Accurace,int Printmodel)
        {
            double retainerrora = 0; double stepcutrate = 0; double stepcutrateretain = 0;double stepsize = Stepsize;
            List<double[,]> result = new List<double[,]>();
            int lastiterator = Iterator - 1;
            for (int i = 0; i < Iterator; i++)
            {

                if (stepcutrate > 0 || (stepcutrate < 0 && stepcutrate < stepcutrateretain)) { stepsize = stepsize * (1 - 0.8); }
                if ((stepcutrate < 0 && stepcutrate > stepcutrateretain) || stepcutrate < 10) { stepsize = stepsize * (1 + 0.2); }
                if (stepcutrate * stepcutrateretain < 0) { stepsize = stepsize / 2; }
                List<double[,]> re = MechineAlgorithm.GradientDescent(InitialPosition, x, y, ClassName, MethodName, GradsStepsize, stepsize);

                InitialPosition = re[2]; stepcutrateretain = stepcutrate; stepcutrate = (re[0][0, 0] - retainerrora) / stepsize;
                retainerrora = re[0][0, 0];
                if (re[0][0, 0] < Accurace || i== lastiterator) {
                    double[,] NewStepsize = { { stepsize } }; double[,] IteratorNum = { { Convert.ToDouble(i) } };
                    result.Add(re[0]); result.Add(re[1]); result.Add(re[2]); result.Add(NewStepsize); result.Add(IteratorNum);
                    if (re[0][0, 0] < Accurace || i == lastiterator ) { break; }
                }
            }
            return result;
        }
        /// <summary>
        /// 高斯牛頓法
        /// </summary>
        /// <param name="InitialthetaPosition">參數初始矩陣</param>
        /// <param name="x">輸入矩陣</param>
        /// <param name="y">輸出矩陣</param>
        /// <param name="Detalx">精確度</param>
        /// <param name="Classname">調用函數的類名</param>
        /// <param name="Methodname">調用函數的方法名</param>
        /// <param name="Param">調用參數</param>
        /// <returns></returns>
        public static double[,] GaussNewton(double[,] InitialthetaPosition, double[,] x, double[,] y, double Detalx, string Classname, string Methodname, object[] Param)
        {
            Type t = Type.GetType(Classname);
            MethodInfo minfo = t.GetMethod(Methodname);
            double[,] Yestimate = (double[,])minfo.Invoke(t, new object[] { x, InitialthetaPosition });
            double[,] Jacobi =Calculas.JacobiMatrix(InitialthetaPosition, x, Detalx, Classname, Methodname, Param);
            double[,] DeltaTheata =
                Matrix.Product(
                Matrix.Inverse(Matrix.Product(Matrix.Trans(Jacobi), Jacobi)),
                Matrix.Product(Matrix.Trans(Jacobi), Matrix.Subtract(y, Matrix.Trans(Yestimate))));
            double[,] NewThetaPosition= InitialthetaPosition;
            for(int i=0;i< InitialthetaPosition.GetLength(1); i++)
            {
                NewThetaPosition[0, i] = InitialthetaPosition[0, i] + DeltaTheata[i, 0]*0.1;
            }
            return NewThetaPosition;
        }
    }
}
