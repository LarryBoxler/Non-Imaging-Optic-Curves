using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Non_Imaging_Optics
{
    //secant method for solving algebraic equations
    public class Secant
    {
        public delegate double Function(double x); //declare a delegate that takes double and returns double

        public static void secant(int step_number, double point1, double point2, Function f)
        {
            double p2, p1, p0, prec = .0001f; //set precision to .0001
            int i;
            p0 = f(point1);
            p1 = f(point2);
            p2 = p1 - f(p1) * (p1 - p0) / (f(p1) - f(p0)); //secant formula

            for (i = 0; System.Math.Abs(p2) > prec && i < step_number; i++) //iterate till precision goal is not met or the maximum //number of steps is reached
            {
                p0 = p1;
                p1 = p2;
                p2 = p1 - f(p1) * (p1 - p0) / (f(p1) - f(p0));
            }
            if (i < step_number)
                Console.WriteLine(p2); //method converges
            else
                Console.WriteLine("{0}.The method did not converge", p2);//method does not converge
        }
    }

    //Simpson integration algorithm
    //calculate the integral of f(x) between x=a and x=b by spliting the interval in step_number steps
    class Integral
    {
        public delegate double Function(double x); //declare a delegate that takes and returns double 
        public static double integral(Function f, double a, double b, int step_number)
        {
            double sum = 0;
            double step_size = (b - a) / step_number;
            for (int i = 0; i < step_number; i = i + 2) //Simpson algorithm samples the integrand in several point which significantly improves //precision.
                sum = sum + (f(a + i * step_size) + 4 * f(a + (i + 1) * step_size) + f(a + (i + 2) * step_size)) * step_size / 3; //divide the area under f(x)     //into step_number rectangles and sum their areas 
            return sum;
        }
    }

    //fourth order Runge Kutte method for y'=f(t,y);
    //solve first order ode in the interval (a,b) with a given initial condition at x=a and fixed step h.
    class Runge
    {
        public delegate double Function(double t, double y); //declare a delegate that takes a double and returns double
        public static void runge(double a, double b, double value, double step, Function f)
        {
            double t, w, k1, k2, k3, k4;
            t = a;
            w = value;
            for (int i = 0; i < (b - a) / step; i++)
            {
                k1 = step * f(t, w);
                k2 = step * f(t + step / 2, w + k1 / 2);
                k3 = step * f(t + step / 2, w + k2 / 2);
                k4 = step * f(t + step, w + k3);
                w = w + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                t = a + i * step;
                Console.WriteLine("{0} {1} ", t, w);
            }
        }
    }
}

