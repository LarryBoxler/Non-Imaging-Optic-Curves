using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Non_Imaging_Optics;

namespace ConsoleApplication1
{
    class Program
    {
        static void Main(string[] args)
        {

            NI2dVector i = new NI2dVector(new NI2dPoint(9.984,10.093), new NI2dPoint(7.445,3.667));
            NI2dVector n = new NI2dVector(-1,0);

            //List<NI2dPoint> winmpPoints = NIUtilities.convergingOvalptpt(F, n1, v1,G,n2, S);

            //foreach (NI2dPoint point in winmpPoints)
            //{
            //    Console.WriteLine("{0},{1}", point.x1, point.x2);
            //}

            NI2dVector r = NIUtilities.refractedRay(i, n, 1.489, 1.000);
            Console.WriteLine("{0},{1}", r.x1, r.x2);

            Console.ReadKey(); 

        }
    }
}
