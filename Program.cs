using System;

namespace newone
{
    class Program
    {
        static void Main(string[] args)
        {
            var algo = new Algorithm();
            algo.Calc();
            Console.ReadLine();
        }
    }

    public class Algorithm
    {
        private double[,] u2 = new double[2, 2];
        private double[,] u3 = new double[3, 3];
        private double[,] u4 = new double[4, 4];
        private double[,] u5 = new double[5, 5];

        private double Det3()  // todo wrong indexation
        {
            var det3 = u3[0, 0] * u3[1, 1] * u3[2, 2] + u3[2, 0] * u3[0, 1] * u3[1, 2]
                + u3[1, 0] * u3[2, 1] * u3[0, 2] - u3[2, 0] * u3[1, 1] * u3[0, 2]
                - u3[1, 0] * u3[0, 1] * u3[2, 2] - u3[0, 0] * u3[2, 1] * u3[1, 2];
            return det3;
        }

        private double Det4() // todo wrong indexation
        {
            double aa = 0;

            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    u3[k, l] = u4[k, l];
                    aa = aa + u4[3, 3] * Det3();
                }
            }

            for (int k = 0; k < 3; k++)
            {
                u3[k, 2] = u4[k, 3];
                aa = aa - u4[3, 2] * Det3();
            }

            for (int k = 0; k < 3; k++)
            {
                u3[k, 1] = u4[k, 2];
                aa = aa + u4[3, 1] * Det3();
            }

            for (int k = 0; k < 3; k++)
            {
                u3[k, 0] = u4[k, 1];
                aa = aa - u4[3, 0] * Det3();
            }

            return aa;
        }

        private double Det5()  // todo wrong indexation
        {
            double aa = 0;

            for (int k = 0; k < 4; k++)
            {
                for (int l = 0; l < 4; l++)
                {
                    u4[k, l] = u5[k, l];
                    aa = aa + u5[4, 4] * Det4();
                }
            }

            for (int k = 0; k < 4; k++)
            {
                u4[k, 3] = u5[k, 4];
                aa = aa - u5[4, 3] * Det4();
            }

            for (int k = 0; k < 4; k++)
            {
                u4[k, 2] = u5[k, 3];
                aa = aa + u5[4, 2] * Det4();
            }

            for (int k = 0; k < 4; k++)
            {
                u4[k, 1] = u5[k, 2];
                aa = aa - u5[4, 1] * Det4();
            }

            for (int k = 0; k < 4; k++)
            {
                u4[k, 0] = u5[k, 1];
                aa = aa + u5[4, 0] * Det4();
            }

            return aa;
        }


        private int Sign(double x)
        {
            if (x > 0)
                return 1;

            if (x < 0)
                return -1;

            return 0;
        }

        // функцiя, яка апроксимується
        private double Func(double x)
        {
            return Math.Cos(x);
        }

        // апроксимацiйна вага w
        private double W(double x)
        {
           return 1; // найкраще абсолютне наближення
           // return Func(x); //найкраще вiдносне наближення
        }

        // апроксимуючий полiном P
        private double P(double x)
        {
            double PP = 0;
            double xx = 1;
            for (int i = 0; i <= m; i++)
            {
                PP = PP + aa[i] * xx;
                xx = xx * x;
            }

            return PP;
        }

        double a;
        double b;
        double eps; // задана точнiсть
        double ro;
        double r1;
        double d0;
        double d1;
        double xx;
        double dx; // крок для визначення максимального вiдхилення на промiжку
        double xm; // {точка, в якiй максимальне вiдхилення} 
        double d;
        double tt;
        private double[] t0 = new double[5]; //  {точки альтернансу: попередня i наступна iтерацiї}
        private double[] t1 = new double[5]; //  {точки альтернансу: попередня i наступна iтерацiї}
        double nm; //{кiлькiсть точок на пормiжку для визначення максимального вiдхилення}
        int m, m1, v; // {m - порядок полiнома; m=3} 
        private double[] aa = new double[5]; // { aa[0]-aa[3] - коефiцiєнти шуканого апроксимуючого полiнома, aa[4] - miu
        private double[,] cc = new double[5, 5]; // коефiцiєнти системи рiвнянь для aa // TODO: indexes 1 and 0!!!!!!
        private double[] bb = new double[5]; //{вiльнi члени цiєї системи рiвнянь} 
        private double[,] u50 = new double[5, 5];
        int j; // номер iтерацiї 


        // Основна програма
        public void Calc()
        {
            m = 3;
            m1 = m + 1;
            a = 0;
            b = 1;
            eps = 0.01;
            d = (b - a) / m1;
            nm = 200;
            dx = (b - a) / nm;

            for (int k = 0; k <= m1; k++)
            {
                t0[k] = a + d * k;
            }

            ro = 1.0;
            aa[4] = 0; // початковi значення для запуску iтерацiйного циклу
            j = 0;

            while (ro - Math.Abs(aa[4]) > eps * Math.Abs(aa[4]))
            {
                for (int i = 0; i < m + 2; i++)
                {
                    tt = t0[i];
                    cc[i, 0] = 1;

                    for (int k = 1; k <= m; k++)
                        cc[i, k] = cc[i, k] * tt;

                    cc[i, 4] = W(tt) * (1.0 - 2.0 * ((i + 1) % 2.0));
                    bb[i] = Func(tt);
                }

                for (int i = 0; i < 5; i++)
                {
                    for (int k = 0; k < 5; k++)
                        u5[i, k] = cc[i, k]; // TODO: check!!!!!!!!
                }

                for (int i = 0; i < 5; i++)
                {
                    for (int k = 0; k < 5; k++)
                        u50[i, k] = u5[i, k];  // TODO: check!!!!!!!!
                }

                d0 = Det5(); // основний визначник

                if (d0 == 0)
                    Console.WriteLine($"Система рiвнянь вироджена, iтерацiя: {j}");

                for (int k = 0; k < 5; k++)
                {
                    for (int i = 0; i < 5; i++)
                    {
                        for (int l = 0; l < 5; l++)
                        {
                            u5[i, l] = u50[i, l];
                        }
                    }

                    for (int i = 0; i < 5; i++)
                    {
                        u5[i, k] = bb[i];
                    }

                    d1 = Det5(); //  { додатковi визначники}

                    aa[k] = d1 / d0; 
                }

                ro = 0; xm = 0; xx = a;

                while (xx <= b)
                {
                    r1 = Math.Abs(Func(xx) - P(xx));

                    if (r1 > ro)
                    {
                        ro = r1;
                        xm = xx;
                    }

                    xx = xx + dx;
                }

                if (ro - Math.Abs(aa[4]) > eps * Math.Abs(aa[4]))
                {
                    if (xm < t0[0])
                    {
                        // Випадок 2
                        if (Sign(Func(xm) - P(xm)) == Sign(Func(t0[0]) - P(t0[0])))
                        {
                            t0[0] = xm;
                        }
                        else
                        {
                            for (int k = m1; k >= 0; k--)
                            {
                                t0[k] = t0[k - 1];
                            }
                            t0[0] = xm;
                        }
                    }
                    else if (xm > t0[m1])
                    {
                        // Випадок 3
                        if (Sign(Func(xm) - P(xm)) == Sign(Func(t0[m1]) - P(t0[m1])))
                        {
                            t0[m1] = xm;
                        }
                        else
                        {
                            for (int k = 0; k < m; k++)
                            {
                                t0[k] = t0[k + 1];
                            }
                            t0[m1] = xm;
                        }
                    }
                    else
                    {
                        // Випадок 1
                        v = 0;


                        while (v < m + 1 && xm >= t0[v + 1])
                        {
                            v = v + 1;
                            // WUT
                        }
                        if (Sign(Func(xm) - P(xm)) == Sign(Func(t0[v]) - P(t0[v])))
                        {
                            t0[v] = xm;
                        }
                        else
                        {
                            t0[v + 1] = xm;
                        }
                    }
                }


                j = j + 1;
            }

            Console.WriteLine($"Кiлькiсть iтерацiй:{j}, miu: {aa[4]}. Результат: ");
            for (int i = 0; i < m; i++)
            {
                Console.WriteLine($"a({i})={aa[i]}");
            }
        }

    }
}
