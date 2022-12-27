using MathNet.Numerics.LinearAlgebra;
using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

class NetworkSyntas
{
    class Program
    {
        #region Данные и загрузка и вывод
        //Данные для ФунРекАлгРКР
        static public int NumberStream { get; set; }//коичество потоков
        static public Matrix<double> RouteMatrinput { get; set; }//МаршрутМатр

        static public Vector<double> Vinput { get; set; }//число ЕКР доступных в линии

        static public Vector<double> Ainput { get; set; }//Интенсивность поступления заявок

        static public Vector<double> Binput { get; set; }//Интенсивность поступления заявок в канальном ресурсе

        static public double OptVelEr { get; set; }//ДопВелОшиб

        static public Vector<double> VESinput { get; set; }//Вектор весов

        static public Vector<double> RKRinput { get; set; }//Вектор RKR


        /// <summary>
        /// Загрузка данных из файла
        /// </summary>
        private static void LoadData()
        {
            try
            {
                Console.Write("Введите количество потоков:");
                NumberStream = Int32.Parse(Console.ReadLine());
                double[] potoki = new double[NumberStream];
                int k = 0;
                using (var streamreader = new StreamReader(@"Datainput.txt"))
                {
                    while (!streamreader.EndOfStream)
                    {
                        potoki = streamreader.ReadLine().Split('/')[0].Split(' ', (char)NumberStream).Select(n => double.Parse(n, CultureInfo.InvariantCulture)).ToArray();
                        VESinput = CreateVector.DenseOfArray(potoki);
                        potoki = streamreader.ReadLine().Split('/')[0].Split(' ', (char)NumberStream).Select(n => double.Parse(n)).ToArray();
                        RKRinput = CreateVector.DenseOfArray(potoki);
                        double[,] routmatr = new double[RKRinput.Count, RKRinput.Count];
                        potoki = streamreader.ReadLine().Split('/')[0].Split(' ', (char)((char)NumberStream * (char)NumberStream)).Select(n => double.Parse(n)).ToArray();
                        Parallel.For(0, RKRinput.Count, (i) =>
                        {
                            for (int j = 0; j < RKRinput.Count; j++)
                            {
                                routmatr[i, j] = potoki[k];
                                k++;
                            }
                        });
                        RouteMatrinput = CreateMatrix.DenseOfArray(routmatr);
                        potoki = streamreader.ReadLine().Split('/')[0].Split(' ', (char)NumberStream).Select(n => double.Parse(n)).ToArray();
                        Vinput = CreateVector.DenseOfArray(potoki);
                        potoki = streamreader.ReadLine().Split('/')[0].Split(' ', (char)NumberStream).Select(n => double.Parse(n)).ToArray();
                        Ainput = CreateVector.DenseOfArray(potoki);
                        potoki = streamreader.ReadLine().Split('/')[0].Split(' ', (char)NumberStream).Select(n => double.Parse(n)).ToArray();
                        Binput = CreateVector.DenseOfArray(potoki);
                        potoki = streamreader.ReadLine().Split('/')[0].Split(' ', (char)NumberStream).Select(n => double.Parse(n)).ToArray();
                        OptVelEr = Math.Pow(potoki[0], potoki[1]);
                    }
                }
            }
            catch (Exception)
            {
                Console.WriteLine("Неправильные данные");
            }

        }
        /// <summary>
        /// Вывод вектора на экран
        /// </summary>
        /// <param name="vector"></param>
        private static void PrintVector(Vector<double> vector, string NameVector)
        {
            Console.Write("Вектор:" + NameVector + "={");
            for (int i = 0; i < vector.Count; i++)
            {
                Console.Write(vector[i] + " ");
            }
            Console.WriteLine("}");
        }
        #endregion
        #region Main
        static void Main(string[] args)
        {
            try
            {
                object[] res = new object[4];
                LoadData();
                //Matrix<double> k = CreateMatrix.Dense(3, 3, (double)(0));
                //k[0, 1] = 2;
                //k[1, 0] = 3;
                //Console.WriteLine(k);
                Algorithm();
                //Task task = new Task(async () => res = await CalcMultNetworkMetricProsMMZ(Ainput, Binput, Vinput, VESinput, RKRinput, RouteMatrinput, OptVelEr));
                //task.Start();
                //Console.WriteLine(res[0]);
                //Console.WriteLine(res[1]);
                //Console.WriteLine(res[2]);
                //Console.WriteLine(res[3]);

            }
            catch (Exception)
            {
                Console.WriteLine("Неправильные данные");
            }

            Console.ReadLine();
        }
        async static void Algorithm()
        {
            Console.WriteLine(await CalcMultNetworkMetricProsMMZ(Ainput, Binput, Vinput, VESinput, RKRinput, RouteMatrinput, OptVelEr));
        }
        #endregion
        //
        #region ФунРекАлгРКР
        /// <summary>
        /// ФунРекАлгРКР
        /// </summary>
        private static async Task<Vector<double>> RecursionAlghoritm(Vector<double> A, Vector<double> B, double V, Vector<double> WeightVector, Vector<double> VarRKRVector)
        {
            double N = 0;
            Vector<double> VectorP = CreateVector.Dense((int)(V + 1), (double)0);
            Vector<double> VectorRO = CreateVector.Dense((int)(V + 1), (double)0);
            Vector<double> VectorPi = CreateVector.Dense(A.Count, (double)0);
            Vector<double> VectorM = CreateVector.Dense(A.Count, (double)0);
            VectorP[0] = 1;
            int iter_P = 1;
            int iter_B = 0;
            while (iter_P <= V)
            {
                VectorP[iter_P] = 0;
                iter_B = 0;
                while (iter_B <= A.Count - 1)
                {
                    if (iter_P - B[iter_B] >= 0)
                    {
                        VectorP[iter_P] += (((double)1 / iter_P)) * A[iter_B] * B[iter_B] * VectorP[iter_P - (int)B[iter_B]] * (1 - (await FunctionLoss(iter_P - B[iter_B], B[iter_B], V, WeightVector[iter_B], VarRKRVector[iter_B])));
                    }
                    iter_B++;
                }
                iter_P++;
            }
            N = VectorP.Sum();
            VectorP.Divide(N, VectorRO);
            iter_B = 0;
            while (iter_B <= A.Count - 1)
            {
                for (int i = 0; i <= V; i++)
                {
                    VectorPi[iter_B] += VectorRO[i] * (await FunctionLoss(i, B[iter_B], V, WeightVector[iter_B], VarRKRVector[iter_B]));
                }
                VectorM[iter_B] = A.At(iter_B) * B.At(iter_B) * (1 - VectorPi.At(iter_B));
                iter_B++;
            }
            return VectorPi;
        }
        #endregion
        #region ФунПотРКР

        /// <summary>
        /// ФунПотРКР
        /// </summary>
        /// <param name="varRKR"></param>
        /// <returns></returns>
        async private static Task<double> FunctionLoss(double iter, double B, double V, double Weight, double varRKR)
        {
            double probability = 0;
            switch ((int)varRKR)
            {
                case 0:
                    {
                        probability = await RKR0(iter, B, V, Weight, probability);
                        break;
                    }
                case 1:
                    {
                        probability = await RKR1(iter, B, V, Weight, probability);
                        break;
                    }
                case 2:
                    {
                        probability = await RKR2(iter, B, V, Weight, probability);
                        break;
                    }
                default:
                    break;
            }
            return probability;
        }

        /// <summary>
        /// Возвращает массив ФунПотРКР если ВарРКР=2
        /// </summary>
        /// <returns></returns>
        async private static Task<double> RKR2(double iter, double B, double V, double Weight, double probability)
        {

            if (iter <= V - B)
            {
                if (Weight == 1)
                    probability = 0;
                else
                    probability = Math.Pow((iter / V), -Math.Log(1 - Weight));
            }
            else
                probability = 1;
            await Task.Delay(0);
            return probability;
        }
        /// <summary>
        /// Возвращает массив ФунПотРКР если ВарРКР=1
        /// </summary>
        /// <returns></returns>
        async private static Task<double> RKR1(double iter, double B, double V, double Weight, double probability)
        {

            if ((iter / V) < Weight & (iter <= V - B))
                probability = 0;
            else
                probability = 1;
            await Task.Delay(0);
            return probability;
        }
        /// <summary>
        /// Возвращает массив ФунПотРКР если ВарРКР=0
        /// </summary>
        /// <returns></returns>
        async private static Task<double> RKR0(double iter, double B, double V, double Weight, double probability)
        {

            if (iter <= V - B)
                probability = 0;
            else
                probability = 1;

            if (Weight >= 1)
                Weight = 1;
            else if (Weight < 0)
                Weight = 0;

            await Task.Delay(0);
            return probability;
        }
        #endregion

        #region МЕТОД ПРОСЕЯННОЙ НАГРУЗКИ МСС С РКР
        #region ФунУмноженВектНаПроизвЭлемСтолбМатр
        /// <summary>
        /// Умножение элемента вектора на столбец матрицы
        /// </summary>
        /// <param name="intensityvector"></param>
        /// <param name="probabilityMatr"></param>
        /// <param name="j"></param>
        /// <param name="i"></param>
        /// <returns></returns>
        async static private Task<Vector<double>> MultiplColumnMatrOnElemVector(Vector<double> intensityvector, Matrix<double> probabilityMatr)
        {
            Vector<double> result = CreateVector.Dense(intensityvector.Count, (double)1);
            for (int i = 0; i < result.Count; i++)
            {
                for (int j = 0; j < probabilityMatr.ToColumnArrays().Length; j++)
                {
                    result[i] *= probabilityMatr[j, i];
                }
                result[i] *= intensityvector[i];
            }

            await Task.Delay(0);
            return result;
        }
        #endregion
        #region ФунОбновлСтрокиМатрB
        /// <summary>
        /// Функция обновления строки матрицы B
        /// </summary>
        /// <param name="B"></param>
        /// <param name="u"></param>
        /// <param name="newrowmatrB"></param>
        /// <param name="activeVector"></param>
        /// <param name="i"></param>
        /// <returns></returns>
        async private static Task<Matrix<double>> UpdateColumnMatr(Matrix<double> BF, int u, Vector<double> newrowmatrB, Vector<double> activeVector)
        {
            for (int j = 0; j < activeVector.Count; j++)
            {
                if (activeVector[j] != 0)
                {
                    BF[u, j] = newrowmatrB[j];
                }
            }
            await Task.Delay(0);
            return BF;
        }
        #endregion
        #region ФунВычОтнОшибВычисл
        /// <summary>
        /// Функция вычисления матрицы Относительной Ошибки Вычисления
        /// </summary>
        /// <param name="B"></param>
        /// <param name="newB"></param>
        /// <returns></returns>
        async private static Task<Matrix<double>> RelativityErrorCalc(Matrix<double> B, Matrix<double> newB)
        {
            Matrix<double> relatercalc = CreateMatrix.Dense(B.RowCount, B.ColumnCount, (double)0);
            Parallel.For(0, newB.RowCount, (i) =>
            {
                int j = 0;
                for (j = 0; j < newB.ColumnCount; j++)
                {
                    if (newB[i, j] != 0)
                        relatercalc[i, j] = Math.Abs(newB[i, j] - B[i, j]) / newB[i, j];
                    else
                        relatercalc[i, j] = 0;
                }
            });
            await Task.Delay(0);
            return relatercalc;
        }
        #endregion
        #region ФунФунОбновлМатрB
        /// <summary>
        /// Функция обновления матрицы B
        /// </summary>
        /// <param name="V"></param>
        /// <param name="RouteMatr"></param>
        /// <param name="B"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="VES"></param>
        /// <param name="RKR"></param>
        /// <returns></returns>
        async static private Task<Matrix<double>> FunctionUpdateMatrB(Vector<double> V, Matrix<double> RouteMatr, Matrix<double> BF, Vector<double> a, Vector<double> b, Vector<double> VES, Vector<double> RKR)
        {

            #region Вспомогатльные величины
            Matrix<double> A = CreateMatrix.DenseOfDiagonalVector(a);//Диагональная матрица
            Matrix<double> B = CreateMatrix.DenseOfMatrix(BF);
            Matrix<double> adiag = CreateMatrix.DenseOfDiagonalVector(a);
            int iter = 0;//u
            #endregion
            #region Основные данные
            Vector<double> VectorPokasActLink;//ВекторПоказАктЗвено
            Vector<double> VectorIntensActPotok;//ВектИнтАктПотоков
            Matrix<double> MatrVerNeProceiv;//МатрВерНеПросеив
            Vector<double> VectorIntProcPotok;//ВектИнтПросПотоков
            Vector<double> NewRowMatrB;//НовСтрокаМатрB
            #endregion
            while (iter <= V.Count - 1)
            {
                VectorPokasActLink = CreateVector.Dense<double>(V.Count, 1);
                VectorPokasActLink[iter] = 0;
                A = CreateMatrix.DenseOfDiagonalVector(VectorPokasActLink);
                VectorIntensActPotok = CreateVector.DenseOfVector((adiag.Multiply(RouteMatr.Transpose().Column(iter))));
                MatrVerNeProceiv = CreateMatrix.DenseOfMatrix(A.Multiply(B).SubtractFrom(1));
                VectorIntProcPotok = CreateVector.DenseOfVector(await MultiplColumnMatrOnElemVector(VectorIntensActPotok, MatrVerNeProceiv));
                NewRowMatrB = CreateVector.DenseOfVector(await RecursionAlghoritm(VectorIntProcPotok, b, V[iter], VES, RKR));
                B = await UpdateColumnMatr(B, iter, NewRowMatrB, RouteMatr.Transpose().Column(iter));
                iter++;
            }
            await Task.Delay(0);
            return B;
        }
        #endregion
        #region ФунРасчВерПотДляПотока
        /// <summary>
        /// Функция для расчета Вероятности потери для потока
        /// </summary>
        /// <param name="ColumnFromB"></param>
        /// <param name="ColumnFromRouteMatr"></param>
        /// <returns></returns>
        async static private Task<double> CalcProbabilForStream(Vector<double> ColumnFromB, Vector<double> ColumnFromRouteMatr)
        {
            double VerNoPotForStream = 1;
            Parallel.For(0, ColumnFromB.Count, (iter) =>
            {
                if (ColumnFromRouteMatr[iter] == 1)
                    VerNoPotForStream *= ColumnFromB[iter];
            });
            await Task.Delay(0);
            return 1 - VerNoPotForStream;
        }
        #endregion
        #region ФунРасчМультСетиССетПросММЗ
        /// <summary>
        /// Функция расчета мультисервисной сети МетПрос ММЗ
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="V"></param>
        /// <param name="VES"></param>
        /// <param name="RKR"></param>
        /// <param name="RouterMatr"></param>
        /// <param name="OptionalSizeEr"></param>
        /// <returns></returns>
        async static private Task<Matrix<double>> CalcMultNetworkMetricProsMMZ(Vector<double> a, Vector<double> b, Vector<double> V, Vector<double> VES, Vector<double> RKR, Matrix<double> RouterMatr, double OptionalSizeEr)
        {
            Matrix<double> MatrixB = CreateMatrix.Dense(RouterMatr.RowCount, RouterMatr.ColumnCount, (double)0);
            Matrix<double>[] Bcopy = new Matrix<double>[30];
            Parallel.For(0, Bcopy.Length, (i) =>
            {
                Bcopy[i] = CreateMatrix.Dense(MatrixB.RowCount, MatrixB.ColumnCount, (double)0);
            });
            MatrixB.CopyTo(Bcopy[0]);
            int iter = 1;//номер шага
            Matrix<double> newB = CreateMatrix.DenseOfMatrix(await FunctionUpdateMatrB(V, RouterMatr, MatrixB, a, b, VES, RKR));
            newB.CopyTo(Bcopy[iter]);
            Matrix<double> MatrRelatErCalc = CreateMatrix.DenseOfMatrix(await RelativityErrorCalc(MatrixB, newB));
            double SumCalcEr = MatrRelatErCalc.RowSums().Sum();
            newB.CopyTo(MatrixB);
            while (SumCalcEr > OptionalSizeEr)
            {
                iter++;
                newB = await FunctionUpdateMatrB(V, RouterMatr, MatrixB, a, b, VES, RKR);
                MatrRelatErCalc = await RelativityErrorCalc(MatrixB, newB);
                newB.CopyTo(Bcopy[iter]);
                SumCalcEr = MatrRelatErCalc.RowSums().Sum();
                newB.CopyTo(MatrixB);
            }
            Vector<double> VectorPi = CreateVector.Dense(a.Count, (double)0);
            Parallel.For(0, a.Count, async (i) =>
            {
                VectorPi[i] = await CalcProbabilForStream(MatrixB.Column(i).SubtractFrom(1), RouterMatr.Column(i));
            });
            Console.WriteLine(VectorPi);
            Console.WriteLine(iter);
            return MatrixB;
        }
        #endregion
        #endregion
    }
}

