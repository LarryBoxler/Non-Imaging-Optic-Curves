using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Non_Imaging_Optics
{
    public class NI3dVector
    {
        public double x1;
        public double x2;
        public double x3;

        public NI3dVector(double t, double u, double v)
        {
            x1 = t;
            x2 = u;
            x3 = v;
        }

        public double Magnitude
        {
            get
            {
                return Math.Sqrt((x1 * x1) + (x2 * x2) + (x3 * x3));
            }

            set
            {
                if (value < 0)
                { throw new ArgumentOutOfRangeException("value", value, NEGATIVE_MAGNITUDE); }

                //if (this == origin)
                //{ throw new ArgumentException(ORIGIN_VECTOR_MAGNITUDE, "this"); }
                double magnitude = Magnitude;
                x1 = x1 * (value / magnitude);
                x2 = x2 * (value / magnitude);
                x3 = x3 * (value / magnitude);
            }
        }

        private const string NEGATIVE_MAGNITUDE =
        "The magnitude of a Vector must be a positive value, (i.e. greater than 0)";

        private const string ORIGIN_VECTOR_MAGNITUDE =
        "Cannot change the magnitude of Vector(0,0,0)";

        public static NI3dVector operator *(NI3dVector v1, double s1)
        {
            return
            (
                new NI3dVector
                (
                    v1.x1 * s1,
                    v1.x2 * s1,
                    v1.x3 * s1
                )
            );
        }

        public void Normalize()
        {
           if ( this.Magnitude == 0 )
           {
              throw new DivideByZeroException( NORMALIZE_0 );
           }
           else
           {
                double inverse = 1 / this.Magnitude;
                this.x1 = this.x1 * inverse;
                this.x2 = this.x2 * inverse;
                this.x3 = this.x3 * inverse;
           }
        }

        private const string NORMALIZE_0 = "Can not normalize a vector when" +
            "it's magnitude is zero";
       
        public static double DotProduct(NI3dVector v1, NI3dVector v2)
        {
            return
            (
                v1.x1 * v2.x1 +
                v1.x2 * v2.x2 +
                v1.x3 * v2.x3
            );

        }

        public double DotProduct(NI3dVector v2)
        {
            return DotProduct(this, v2);
        }

        public static double AngleBetween(NI3dVector v1, NI3dVector v2)
        {
            return Math.Acos(NI3dVector.DotProduct(v1, v2) / (v1.Magnitude * v2.Magnitude));
        }

        public double AngleBetween(NI3dVector v2)
        {
            return AngleBetween(this, v2);
        }




    }

    public class NI2dVector
    {
        public double x1;
        public double x2;

        public NI2dVector(double t, double u)
        {
            x1 = t;
            x2 = u;
        }

        public NI2dVector(NI2dPoint pt2, NI2dPoint pt1)
        {
            x1 = pt2.x1 - pt1.x1;
            x2 = pt2.x2 - pt1.x2;
        }

        public NI2dVector(NI2dPoint pt1)
        {
            x1 = pt1.x1;
            x2 = pt1.x2;   
        }

        public static NI2dVector operator *(NI2dVector v1, double s1)
        {
            return
            (
                new NI2dVector
                (
                    v1.x1 * s1,
                    v1.x2 * s1
                )
            );
        }

        public static NI2dVector operator *(double s1, NI2dVector v1)
        {
            return v1 * s1;
        }

        public static NI2dVector operator +(NI2dVector v1, NI2dVector v2)
        {
            return
            (
                new NI2dVector
                (
                    v1.x1 + v2.x1,
                    v1.x2 + v2.x2
                )
             );  
        }

        public static NI2dVector operator -(NI2dVector v1, NI2dVector v2)
        {
            return
            (
                new NI2dVector
                (
                    v1.x1 - v2.x1,
                    v1.x2 - v2.x2
                )
             );
        }

               public static NI2dVector operator /(NI2dVector v1, double c1)
        {
            return
            (
                new NI2dVector
                (
                    v1.x1/c1,
                    v1.x2/c1
                )
             );
        }

        public double Magnitude
        {
            get
            {
                return Math.Sqrt((x1 * x1) + (x2 * x2));
            }

            set
            {
                if (value < 0)
                { throw new ArgumentOutOfRangeException("value", value, NEGATIVE_MAGNITUDE); }

                //if (this == origin)
                //{ throw new ArgumentException(ORIGIN_VECTOR_MAGNITUDE, "this"); }
                double magnitude = Magnitude;
                x1 = x1 * (value / magnitude);
                x2 = x2 * (value / magnitude);
            }
        }

        private const string NEGATIVE_MAGNITUDE =
        "The magnitude of a Vector must be a positive value, (i.e. greater than 0)";

        private const string ORIGIN_VECTOR_MAGNITUDE =
        "Cannot change the magnitude of Vector(0,0,0)";

        public void Normalize()
        {
            if (this.Magnitude == 0)
            {
                throw new DivideByZeroException(NORMALIZE_0);
            }
            else
            {
                double inverse = 1 / this.Magnitude;
                this.x1 = this.x1 * inverse;
                this.x2 = this.x2 * inverse;
            }
        }

        private const string NORMALIZE_0 = "Can not normalize a vector when" +
            "it's magnitude is zero";


        public static double DotProduct(NI2dVector v1, NI2dVector v2)
        {
            return
            (
                v1.x1 * v2.x1 +
                v1.x2 * v2.x2 
            );

        }

        public double DotProduct(NI2dVector v2)
        {
            return DotProduct(this, v2);
        }

        public static double AngleBetween(NI2dVector v1, NI2dVector v2)
        {
            return Math.Acos(NI2dVector.DotProduct(v1, v2) / (v1.Magnitude * v2.Magnitude));
        }

        public double AngleBetween(NI2dVector v2)
        {
            return AngleBetween(this, v2);
        }

        public static double PlaneAngleBetween(NI2dVector v1, NI2dVector v2)
        {
            if (((v2.x1 * v1.x2) - (v2.x2 * v1.x1)) >= 0)
            {
                return NI2dVector.AngleBetween(v1, v2);
            }

            else
                return ((2 * Math.PI) - NI2dVector.AngleBetween(v1, v2));
        }

        public double PlaneAngleBetween(NI2dVector v2)
        {
            return PlaneAngleBetween(this, v2);
        }

        public double AngleH()
        {
            NI2dVector v2 = new NI2dVector(1,0);
            return PlaneAngleBetween(this, v2);
        }

        public void RotateVector(double alpha)
        {
            NIRotMatrix rotation = new NIRotMatrix(alpha);
            double X = this.x1;
            double Y = this.x2;
            this.x1 = rotation.r11 * X + rotation.r12 * Y;
            this.x2 = rotation.r21 * X + rotation.r22 * Y;
        }

        public NI2dPoint ConverttoPoint()
        {
            return new NI2dPoint(this);
        }
        
    }

    public class NI3dPoint
    {
        public double x1;
        public double x2;
        public double x3;

        public NI3dPoint(double x, double y, double z)
        {
            x1 = x;
            x2 = y;
            x3 = z;
        }

        public static double DistanceBetweenPoints(NI3dPoint point1, NI3dPoint point2)
        {
            return Math.Sqrt((Math.Abs(point1.x1 - point2.x1) * Math.Abs(point1.x1 - point2.x1)) + 
                             (Math.Abs(point1.x2 - point2.x2) * Math.Abs(point1.x2 - point2.x2)) +
                             (Math.Abs(point1.x3 - point2.x3) * Math.Abs(point1.x3 - point2.x3)));
        }

    }

    public class NI2dPoint
    {
        public double x1;
        public double x2;

        public NI2dPoint(double x, double y)
        {
            x1 = x;
            x2 = y;
        }

        public NI2dPoint(NI2dVector v1)
        {
            x1 = v1.x1;
            x2 = v1.x2;
        }

        public static NI2dPoint operator +(NI2dPoint pt1, NI2dPoint pt2)
        {
            return
            (
                new NI2dPoint
                (
                    pt1.x1 + pt2.x1,
                    pt1.x2 + pt2.x2
                )
             );
        }

        public static NI2dPoint operator -(NI2dPoint pt1, NI2dPoint pt2)
        {
            return
            (
                new NI2dPoint
                (
                    pt1.x1 - pt2.x1,
                    pt1.x2 - pt2.x2
                )
             );
        }


        public static NI2dPoint operator *(NI2dPoint pt1, double s1)
        {
            return
            (
                new NI2dPoint
                (
                    s1*pt1.x1,
                    s1*pt1.x2 
                )
             );
        }

        public static NI2dPoint operator *(double s1, NI2dPoint pt1)
        {
            return pt1 * s1;                       
        }




        public static double DistanceBetweenPoints(NI2dPoint point1, NI2dPoint point2)
        {
            return Math.Sqrt((Math.Abs(point1.x1 - point2.x1) * Math.Abs(point1.x1 - point2.x1)) +
                             (Math.Abs(point1.x2 - point2.x2) * Math.Abs(point1.x2 - point2.x2)));
        }

        public NI2dVector ConverttoVector()
        {
            return new NI2dVector(this);
        }

    }

    public class NIRotMatrix
    {
        public double r11;
        public double r12;
        public double r21;
        public double r22;

        public NIRotMatrix(double alpha)
        {
            r11 = Math.Cos(alpha);
            r12 = -1 * Math.Sin(alpha);
            r21 = Math.Sin(alpha);
            r22 = Math.Cos(alpha);
        }

    }

    public static class NIUtilities
    {
        
        public static NI2dPoint line2DIntersection(NI2dPoint P, NI2dVector v, NI2dPoint Q, NI2dVector u)
        {
            NIRotMatrix rotation = new NIRotMatrix(Math.PI / 2);
            u.RotateVector(Math.PI / 2);
            NI2dVector up = u;
            NI2dVector qp = new NI2dVector(Q - P);
            NI2dPoint isl = P + (((qp.DotProduct(up)) / v.DotProduct(up)) * v).ConverttoPoint();
            return isl;
        }

        public static List<NI2dPoint> tiltedParabola(double alpha, NI2dPoint F, NI2dPoint P, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> parabolapts = new List<NI2dPoint>();

            double distancePF = NI2dPoint.DistanceBetweenPoints(P, F);
            NI2dPoint PFminus = P - F;

            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dVector v1 = new NI2dVector(Math.Cos(alpha), Math.Sin(alpha));
                NI2dPoint pt1 = new NI2dPoint(Math.Cos(phi+alpha), Math.Sin(phi+alpha));
                NI2dPoint parabPoint = (((distancePF - PFminus.ConverttoVector().DotProduct(v1)) / (1 - Math.Cos(phi))) * pt1)+F;
                parabolapts.Add(parabPoint);
            }

            return parabolapts;
        }

        public static List<NI2dPoint> Ellipse( NI2dPoint F, NI2dPoint G, NI2dPoint P, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> ellipsepts = new List<NI2dPoint>();
            double distanceFP = NI2dPoint.DistanceBetweenPoints(F, P);
            double distancePG = NI2dPoint.DistanceBetweenPoints(P, G);
            double distanceFG = NI2dPoint.DistanceBetweenPoints(F, G);
            double alpha = new NI2dVector(G, F).AngleH();
   

            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dPoint pt1 = new NI2dPoint(Math.Cos(phi+alpha), Math.Sin(phi+alpha));
                NI2dPoint ellipsePoint = (((distanceFP + distancePG) * (distanceFP + distancePG) - (distanceFG * distanceFG)) / (2 * (distanceFP + distancePG) - (2 * distanceFG * Math.Cos(phi)))) * pt1 + F;
                ellipsepts.Add(ellipsePoint);                        
            }

            return ellipsepts;
        }

        public static List<NI2dPoint> Hyperbola(NI2dPoint F, NI2dPoint G, NI2dPoint P, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> hyperbolapoints = new List<NI2dPoint>();
            double distanceFP = NI2dPoint.DistanceBetweenPoints(F, P);
            double distancePG = NI2dPoint.DistanceBetweenPoints(P, G);
            double distanceFG = NI2dPoint.DistanceBetweenPoints(F, G);
            double alpha = new NI2dVector(G, F).AngleH();

            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dPoint pt1 = new NI2dPoint(Math.Cos(phi + alpha), Math.Sin(phi + alpha));
                NI2dPoint hyperbolaPoint = (((distanceFP - distancePG) * (distanceFP - distancePG) - (distanceFG * distanceFG)) / (2 * (distanceFP - distancePG) - (2 * distanceFG * Math.Cos(phi)))) * pt1 + F;
                hyperbolapoints.Add(hyperbolaPoint);
            }

            return hyperbolapoints;

        }

        public static List<NI2dPoint> windingInvolute(NI2dPoint P, NI2dPoint F, double r, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> involutePts = new List<NI2dPoint>();
            double phiP = new NI2dVector(P, F).AngleH() + Math.Asin(r / NI2dPoint.DistanceBetweenPoints(P, F));
            double K = Math.Sqrt((NI2dPoint.DistanceBetweenPoints(P, F) * NI2dPoint.DistanceBetweenPoints(P, F)) - (r * r)) + r*phiP;
            
            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dPoint pt1 = new NI2dPoint(Math.Sin(phi), -1.0*Math.Cos(phi));
                NI2dPoint pt2 = new NI2dPoint(Math.Cos(phi), Math.Sin(phi));
                NI2dPoint involutePoint = (r * pt1) + ((K - (r * phi)) * pt2) + F;
                involutePts.Add(involutePoint);
            }

            return involutePts;
        }

        public static List<NI2dPoint> unwindingInvolute(NI2dPoint P, NI2dPoint F, double r, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> involutePts = new List<NI2dPoint>();
            double phiP = new NI2dVector(P, F).AngleH() - Math.Asin(r / NI2dPoint.DistanceBetweenPoints(P, F));
           
            if (phiP < 0)
                phiP = phiP + 2 * Math.PI;

            double K = Math.Sqrt((NI2dPoint.DistanceBetweenPoints(P, F) * NI2dPoint.DistanceBetweenPoints(P, F)) - (r * r)) - r * phiP;

            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dPoint pt1 = new NI2dPoint(-1*Math.Sin(phi), Math.Cos(phi));
                NI2dPoint pt2 = new NI2dPoint(Math.Cos(phi), Math.Sin(phi));
                NI2dPoint involutePoint = (r * pt1) + ((K + (r * phi)) * pt2) + F;
                involutePts.Add(involutePoint);
            }

            return involutePts;
        }

        public static List<NI2dPoint> windingmacroParabola(double alpha, NI2dPoint F, double r, NI2dPoint P, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> parabolaPts = new List<NI2dPoint>();
            NI2dVector v1 = new NI2dVector(P, F);
            NI2dVector v2 = new NI2dVector(Math.Cos(alpha), Math.Sin(alpha));
            double phiP = v1.AngleBetween(v2) + Math.Asin(r/NI2dPoint.DistanceBetweenPoints(P,F));

            double t1 = Math.Sqrt((NI2dPoint.DistanceBetweenPoints(P, F) * NI2dPoint.DistanceBetweenPoints(P, F)) - (r * r));
            double t2 = (1 - Math.Cos(phiP));
            double t3 = (1 + phiP - Math.Sin(phiP));
            double K = (t1 * t2) + (r * t3);

            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dPoint pt1 = new NI2dPoint(Math.Sin(phi + alpha), -1*Math.Cos(phi+alpha));
                NI2dPoint pt2 = new NI2dPoint(Math.Cos(phi + alpha), Math.Sin(phi + alpha));
                double c1 = (K + r * (Math.Sin(phi) - 1 - phi)) / (1 - Math.Cos(phi));
                NI2dPoint parabpoint = (r * pt1) + (c1 * pt2) + F;
                parabolaPts.Add(parabpoint);                  
            }

            return parabolaPts;
        }

        public static List<NI2dPoint> unwindingmacroParabola(double alpha, NI2dPoint F, double r, NI2dPoint P, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> parabolaPts = new List<NI2dPoint>();
            NI2dVector v1 = new NI2dVector(P, F);
            NI2dVector v2 = new NI2dVector(Math.Cos(alpha), Math.Sin(alpha));
            double phiP = v1.AngleBetween(v2) - Math.Asin(r / NI2dPoint.DistanceBetweenPoints(P, F));
            if (phiP < 0.0)
                phiP = 2 * Math.PI + phiP;

            double t1 = Math.Sqrt((NI2dPoint.DistanceBetweenPoints(P, F) * NI2dPoint.DistanceBetweenPoints(P, F)) - (r * r));
            double t2 = (1 - Math.Cos(phiP));
            double t3 = (1 - phiP + Math.Sin(phiP));
            double K = (t1 * t2) + (r * t3);

            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dPoint pt1 = new NI2dPoint(-1 * Math.Sin(phi + alpha), Math.Cos(phi + alpha));
                NI2dPoint pt2 = new NI2dPoint(Math.Cos(phi + alpha), Math.Sin(phi + alpha));
                double c1 = (K + r * (phi - 1 - Math.Sin(phi))) / (1 - Math.Cos(phi));
                NI2dPoint parabpoint = (r * pt1) + (c1 * pt2) + F;
                parabolaPts.Add(parabpoint);
            }

            return parabolaPts;
        }

        public static List<NI2dPoint> windingmacroEllipse(NI2dPoint F, double r,  NI2dPoint G, NI2dPoint P, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> ellipsepts = new List<NI2dPoint>();
            double distanceFP = NI2dPoint.DistanceBetweenPoints(F, P);
            double f = NI2dPoint.DistanceBetweenPoints(F, G);
            double alpha = new NI2dVector(G, F).AngleH();
            NI2dVector v1 = new NI2dVector(Math.Cos(alpha), Math.Sin(alpha));
            NI2dVector v2 = new NI2dVector(P, F);
            double phiP = v2.AngleBetween(v1) + Math.Asin(r / distanceFP);
            double tp = Math.Sqrt((distanceFP * distanceFP) - (r * r));
            double K = tp + (r * phiP) + Math.Sqrt((f * f) + (r * r) + (tp * tp) - (2 * f * ((tp * Math.Cos(phiP)) + (r * Math.Sin(phiP)))));



            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dPoint pt1 = new NI2dPoint(Math.Sin(phi + alpha), -1 * Math.Cos(phi + alpha));
                NI2dPoint pt2 = new NI2dPoint(Math.Cos(phi + alpha), Math.Sin(phi + alpha));
                double c1numerator = ((K - r * phi) * (K - r * phi)) + 2 * f * r * Math.Sin(phi) - (f * f) - (r * r);
                double c1denominator = 2 * (K - (r*phi) - (f*Math.Cos(phi)));
                double c1 = c1numerator / c1denominator;
                NI2dPoint ellipsePoint = (r * pt1) + (c1 * pt2) + F;               
                ellipsepts.Add(ellipsePoint);
            }

            return ellipsepts;
        }

        public static List<NI2dPoint> unwindingmacroEllipse(NI2dPoint F, double r, NI2dPoint G, NI2dPoint P, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> ellipsepts = new List<NI2dPoint>();
            double distanceFP = NI2dPoint.DistanceBetweenPoints(F, P);
            double f = NI2dPoint.DistanceBetweenPoints(F, G);
            double alpha = new NI2dVector(G, F).AngleH();
            NI2dVector v1 = new NI2dVector(Math.Cos(alpha), Math.Sin(alpha));
            NI2dVector v2 = new NI2dVector(P, F);
            double phiP = v2.AngleBetween(v1) - Math.Asin(r / distanceFP);
            
            if (phiP < 0)
                phiP = phiP + 2 * Math.PI;

            double tp = Math.Sqrt((distanceFP * distanceFP) - (r * r));
            double K = tp - (r * phiP) + Math.Sqrt((f * f) + (r * r) + (tp * tp) - (2 * f * ((tp * Math.Cos(phiP)) - (r * Math.Sin(phiP)))));



            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dPoint pt1 = new NI2dPoint(-1*Math.Sin(phi + alpha), Math.Cos(phi + alpha));
                NI2dPoint pt2 = new NI2dPoint(Math.Cos(phi + alpha), Math.Sin(phi + alpha));
                double c1numerator = ((K + r * phi) * (K + r * phi)) - 2 * f * r * Math.Sin(phi) - (f * f) - (r * r);
                double c1denominator = 2 * (K + (r * phi) - (f * Math.Cos(phi)));
                double c1 = c1numerator / c1denominator;
                NI2dPoint ellipsePoint = (r * pt1) + (c1 * pt2) + F;
                ellipsepts.Add(ellipsePoint);
            }

            return ellipsepts;
        }

        public static List<NI2dPoint> cartesianOvalParallel(NI2dPoint F, double n1, double n2, NI2dPoint P, double alpha, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> ovalpts = new List<NI2dPoint>();
            double distanceFP = NI2dPoint.DistanceBetweenPoints(F, P);
            NI2dVector v1 = new NI2dVector(P, F);
            NI2dVector v2 = new NI2dVector(Math.Cos(alpha), Math.Sin(alpha));
            double phiP = v1.AngleBetween(v2);

            for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
            {
                NI2dPoint pt1 = new NI2dPoint(Math.Cos(phi+alpha), Math.Sin(phi + alpha));
                double c1 = (distanceFP * (n1 - (n2 * Math.Cos(phiP)))) / (n1 - (n2 * Math.Cos(phi)));
                NI2dPoint ovalPoint = (c1 * pt1) + F;
                ovalpts.Add(ovalPoint);
            }

            return ovalpts;
        }


        public static List<NI2dPoint> cartesianOvalConverge(NI2dPoint F, double n1, NI2dPoint G, double n2, NI2dPoint P, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> ovalpts = new List<NI2dPoint>();
            double distanceFP = NI2dPoint.DistanceBetweenPoints(F, P);
            double distancePG = NI2dPoint.DistanceBetweenPoints(P, G);
            double K = (n1 * distanceFP) + (n2 * distancePG);
            double f = NI2dPoint.DistanceBetweenPoints(F, G);
            double phiC=0;

            
            if(Math.Abs(K) <= n1*f)
                phiC = Math.Acos(K/(f*n1));
            else if(Math.Abs(K) > n1*f)
                phiC = 0;
            
            NI2dVector v1 = new NI2dVector(P, F);
            NI2dVector v2 = new NI2dVector(G, F);
            if (n2 * f < K && K < n1 * f && v1.AngleBetween(v2) <= phiC)
            {
                for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
                {
                    double D = (f * n1 - K * Math.Cos(phi)) * (f * n1 - K * Math.Cos(phi)) + (K * K - f * f * n2 * n2) * Math.Sin(phi) * Math.Sin(phi);
                    double alpha = new NI2dVector(G, F).AngleH();
                    double c1numerator = K * n1 - f * n2 * n2 * Math.Cos(phi) - n2 * Math.Sqrt(D);
                    double c1denominator = n1 * n1 - n2 * n2;
                    double c1 = c1numerator / c1denominator;
                    NI2dPoint pt1 = new NI2dPoint(Math.Cos(phi + alpha), Math.Sin(phi + alpha));
                    NI2dPoint ovalPoint = (c1 * pt1) + F;
                    ovalpts.Add(ovalPoint);
                }
            }         

            return ovalpts;
        }

        public static List<NI2dPoint> cartesianOvalDiverge(NI2dPoint F, double n1, NI2dPoint G, double n2, NI2dPoint P, double phi_start, double phi_end, double phi_step)
        {
            List<NI2dPoint> ovalpts = new List<NI2dPoint>();
            double distanceFP = NI2dPoint.DistanceBetweenPoints(F, P);
            double distancePG = NI2dPoint.DistanceBetweenPoints(P, G);
            double K = (n1 * distanceFP) - (n2 * distancePG);
            double f = NI2dPoint.DistanceBetweenPoints(F, G);
            double phiC = 0;


            if (Math.Abs(K) <= n1 * f)
                phiC = Math.Acos(K / (f * n1));
            else if (Math.Abs(K) > n1 * f)
                phiC = 0;

            NI2dVector v1 = new NI2dVector(P, F);
            NI2dVector v2 = new NI2dVector(G, F);
            if (K > -1*n2*f && v1.AngleBetween(v2) >= phiC)
            {
                for (double phi = phi_start; phi < phi_end + phi_step; phi += phi_step)
                {
                    double D = (f * n1 - K * Math.Cos(phi)) * (f * n1 - K * Math.Cos(phi)) + (K * K - f * f * n2 * n2) * Math.Sin(phi) * Math.Sin(phi);
                    double alpha = new NI2dVector(G, F).AngleH();
                    double c1numerator = K * n1 - f * n2 * n2 * Math.Cos(phi) + n2 * Math.Sqrt(D);
                    double c1denominator = n1 * n1 - n2 * n2;
                    double c1 = c1numerator / c1denominator;
                    NI2dPoint pt1 = new NI2dPoint(Math.Cos(phi + alpha), Math.Sin(phi + alpha));
                    NI2dPoint ovalPoint = (c1 * pt1) + F;
                    ovalpts.Add(ovalPoint);
                }
            }

            return ovalpts;
        }

        public static NI2dPoint convergingOvalptpt(NI2dPoint F, double n1, NI2dVector v, NI2dPoint G, double n2, double S)
        {
            NI2dPoint ovalpt = new NI2dPoint(0, 0);
            v.Normalize();
            NI2dVector v1 = new NI2dVector(F,G);
            double C1 = n1 * S + (n2 * n2) * v1.DotProduct(v);
            double C2 = S * S - (n2 * n2) * v1.DotProduct(v1);
            double gamma = 1.0;
            double sigma = 0.0;
            
            if (n1 > n2)
                sigma = -1 * gamma;
            else if (n1 < n2)
                sigma = 1 * gamma;

            double c1numerator = C1 + sigma * Math.Sqrt(C2 * (n2 * n2 - n1 * n1) + C1 * C1);
            double c1denominator = n1 * n1 - n2 * n2;
            double c1 = c1numerator / c1denominator;

            ovalpt = F + c1 * v.ConverttoPoint();

            return ovalpt;

        }

        public static NI2dPoint divergingOvalptpt(NI2dPoint F, double n1, NI2dVector v, NI2dPoint G, double n2, double S)
        {
            NI2dPoint ovalpt = new NI2dPoint(0, 0);
            v.Normalize();
            NI2dVector v1 = new NI2dVector(F, G);
            double C1 = n1 * S + (n2 * n2) * v1.DotProduct(v);
            double C2 = S * S - (n2 * n2) * v1.DotProduct(v1);
            double gamma = -1.0;
            double sigma = 0.0;

            if (n1 > n2)
                sigma = -1 * gamma;
            else if (n1 < n2)
                sigma = 1 * gamma;

            double c1numerator = C1 + sigma * Math.Sqrt(C2 * (n2 * n2 - n1 * n1) + C1 * C1);
            double c1denominator = n1 * n1 - n2 * n2;
            double c1 = c1numerator / c1denominator;

            ovalpt = F + c1 * v.ConverttoPoint();

            return ovalpt;

        }

        public static NI2dPoint reflectingOvalptpt(NI2dPoint F, NI2dVector v, NI2dPoint G, double n, double S)
        {
            NI2dPoint ovalpt = new NI2dPoint(0, 0);
            v.Normalize();
            NI2dVector v1 = new NI2dVector(F, G);
            double c1numerator = (S / n) * (S / n) - v1.DotProduct(v1);
            double c1denominator = 2*((S/n) + v1.DotProduct(v));
            double c1 = c1numerator/c1denominator;

            ovalpt = F + c1 * v.ConverttoPoint();

            return ovalpt;

        }

        public static NI2dPoint planewaveOvalpt(NI2dPoint F, double n1, NI2dVector v, NI2dPoint Q, double n2, NI2dVector n, double S)
        {
            NI2dPoint ovalpt = new NI2dPoint(0, 0);
            NI2dVector v1 = new NI2dVector(Q, F);
            double c1numerator = S - n2 * (v1.DotProduct(n));
            double c1denominator = n1 - n2 * (v.DotProduct(n));
            double c1 = c1numerator/c1denominator;

            ovalpt = F + c1 * v.ConverttoPoint();

            return ovalpt;

        }

        public static NI2dVector reflectedRay(NI2dVector i, NI2dVector n)
        {
            NI2dVector r = new NI2dVector(0, 0);
            i.Normalize();
            n.Normalize();
            r = i - 2*(i.DotProduct(n))*n;

            return r;
        }

        public static NI2dVector refractedRay(NI2dVector i, NI2dVector ns, double n1, double n2)
        {
            NI2dVector r = new NI2dVector(0, 0);
            NI2dVector n = new NI2dVector(0, 0);
            i.Normalize();
            ns.Normalize();

            double ncheck = i.DotProduct(ns);
            if (ncheck >= 0.0)
                n = ns;
            if (ncheck < 0.0)
                n = -1 * ns;

            double delta = 1 - ((n1 / n2) * (n1 / n2) * (1 - (i.DotProduct(n) * i.DotProduct(n))));

            if (delta > 0)
                r = (n1 / n2) * i + ((-1 * (i.DotProduct(n)) * (n1 / n2) + Math.Sqrt(delta)) * n);
            if (delta <= 0)
                r = reflectedRay(i, ns);
            
            return r;
        }

        public static NI2dVector refractedNormal(NI2dVector i, NI2dVector r, double n1, double n2)
        {
            NI2dVector n = new NI2dVector(0, 0);
            i.Normalize();
            r.Normalize();

            n = ((n1 * i) - (n2 * r)) / ((i - r).Magnitude);
            return n;
        }

        public static NI2dVector reflectedNormal(NI2dVector i, NI2dVector r)
        {
            return (i-r)/((i-r).Magnitude);
        }
        
    }

}



