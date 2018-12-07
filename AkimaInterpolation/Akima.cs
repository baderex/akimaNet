using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AkimaInterpolation
{
    public class Akima
    {
        List<double> xo = new List<double>();
        List<double> yo = new List<double>();
        List<double> zo = new List<double>();
        List<double> x = new List<double>();
        List<double> y = new List<double>();
        List<double> z = new List<double>();
        List<double> iwk = new List<double>();
        int ncp = 25;
        int ndp = 0;
        int nx = 40;
        int ny = 40; 


        public List<double> interpp(List<double> _x, List<double> _y, List<double> _z, String duplicate = "error")
        {
            x = _x;
            y = _y;
            z = _z;
            ndp = x.Count;
            int iwkLen = (31 + ncp) * ndp + nx * ny;
            int wk = 5 * ndp;

            xo = seq(x.Min(),x.Max(),nx);
            yo = seq(y.Min(), y.Max(), ny);

            zo = idsfft();

            return zo; 
        }

        public List<double> idsfft()
        {
            int ncp0 = ncp;
            int ndp0 = ndp;
            int nxi0 = nx;
            int nyi0 = ny;

            if(ncp > ndp)
            {
                Console.WriteLine("NCP must be less than ndp");
                return null;
            }

            if (ndp < 4)
            {
                Console.WriteLine("NDP must be > 4");
                return null;
            }

            if (nyi0 < 1 || nxi0 < 1)
            {
                Console.WriteLine("NX or Ny must be > 0");
                return null;
            }

            iwk[0] = ncp0;
            iwk[1] = ndp0;
            iwk[2] = nxi0;
            iwk[3] = nyi0;

            // allocation of storage areas in the iwk array.  (for md = 1, 2, 3)
            int jwipt = 16;
            int jwiwl = 6 * ndp0 + 1;
            double jwngp0 = jwiwl - 1;
            int jwipl = 24 * ndp0 + 1;
            int jwiwp = 30 * ndp0 + 1;
            int jwipc = 27 * ndp0 + 1;
            double jwigp0 = Math.Max(31, 27 + ncp0) * ndp0;

            idtang(ndp0,0, iwk[jwipt], 0, iwk[jwipl], iwk[jwiwl], iwk[jwiwp]);
            idcldp(ndp0, ncp0, iwk[jwipc]);

            return null;
        }

        void idtang(int ndp0, double nt, double jwipt, double nl, double jwipl, double jwiwl,double jwiwp){

            iwk[4] = nt;
            iwk[5] = nl;
        }

        void idcldp(int ndp0, int ncp0, double jwipc)
        {

        }

        void idptip(List<double> XD, List<double> YD, List<double> ZD, int NDP, int NT, List<int> IPT, int NL, List<int> IPL, List<double> PDD, int ITI, double XII, double  YII, out double ZII, out bool MISSII)
        {
            // THIS SUBROUTINE PERFORMS PUNCTUAL INTERPOLATION OR EXTRAPOLA-
            // TION, I.E., DETERMINES THE Z VALUE AT A POINT.
            //
            // THE INPUT PARAMETERS ARE
            //     XD,YD,ZD = ARRAYS OF DIMENSION NDP CONTAINING THE X,
            //           Y, AND Z COORDINATES OF THE DATA POINTS, WHERE
            //           NDP IS THE NUMBER OF THE DATA POINTS,
            //     NT  = NUMBER OF TRIANGLES,
            //     IPT = INTEGER ARRAY OF DIMENSION 3*NT CONTAINING THE
            //           POINT NUMBERS OF THE VERTEXES OF THE TRIANGLES,
            //     NL  = NUMBER OF BORDER LINE SEGMENTS,
            //     IPL = INTEGER ARRAY OF DIMENSION 3*NL CONTAINING THE
            //           POINT NUMBERS OF THE END POINTS OF THE BORDER
            //           LINE SEGMENTS AND THEIR RESPECTIVE TRIANGLE
            //           NUMBERS,
            //     PDD = ARRAY OF DIMENSION 5*NDP CONTAINING THE PARTIAL
            //           DERIVATIVES AT THE DATA POINTS,
            //     ITI = TRIANGLE NUMBER OF THE TRIANGLE IN WHICH LIES
            //           THE POINT FOR WHICH INTERPOLATION IS TO BE
            //           PERFORMED,
            //     XII,YII = X AND Y COORDINATES OF THE POINT FOR WHICH
            //           INTERPOLATION IS TO BE PERFORMED.
            //
            // THE OUTPUT PARAMETERs are
            //     ZII     = INTERPOLATED Z VALUE.
            //     MISSII = LOCIGAL INDICATING MISSING VALUE
            

            // DECLARATION STATEMENTS
            List<double> X = new List<double>(3);
            List<double> Y = new List<double>(3);
            List<double> Z = new List<double>(3);
            List<double> PD = new List<double>(15);
            List<double> ZU = new List<double>(3);
            List<double> ZV = new List<double>(3);
            List<double> ZUU = new List<double>(3);
            List<double> ZUV = new List<double>(3);
            List<double> ZVV = new List<double>(3);

            List<double> LU = new List<double>();
            List<double> LV = new List<double>();

            // ?? initializing output
            MISSII = false;
            ZII = double.NaN;

            // ??: value by reference globaly
            //COMMON / IDPI / ITPV
            int ITPV = new int();

            // ??: value by reference locally
            //EQUIVALENCE(P5, P50)
            double P5 = new double();
            double P50 = new double();

            // PRELIMINARY PROCESSING
            int IT0 = ITI;
            int NTL = NT + NL;

            if (IT0 <= NTL)
            {
                // CALCULATION OF ZII BY INTERPOLATION.
                // CHECKS IF THE NECESSARY COEFFICIENTS HAVE BEEN CALCULATED.
                if (IT0 == ITPV)
                {
                    // CONVERTS XII AND YII TO U - V SYSTEM.
                    double DX = XII - X0;
                    double DY = YII - Y0;
                    double U = AP * DX + BP * DY;
                    double V = CP * DX + DP * DY;
                    // EVALUATES THE POLYNOMIAL.
                    double P0 = P00 + V * (P01 + V * (P02 + V * (P03 + V * (P04 + V * P05))));
                    double P1 = P10 + V * (P11 + V * (P12 + V * (P13 + V * P14)));
                    double P2 = P20 + V * (P21 + V * (P22 + V * P23));
                    double P3 = P30 + V * (P31 + V * P32);
                    double P4 = P40 + V * P41;
                    ZII = P0 + U * (P1 + U * (P2 + U * (P3 + U * (P4 + U * P5))));
                    MISSII = false;
                    return;
                }
                else
                { };

            } else {
                // EXTRAPOLATION OR MISSING VALUE WANTED?
                // ?? if MISSII get initialized as F, then this part is useless. It also might be used as an input...
                if (MISSII)
                {
                    ZII = 0;
                    return; // ?? correct usage?
                }

                int IL1 = IT0 / NTL;
                int IL2 = IT0 - IL1 * NTL;

                if (IL1 == IL2) {
                    // CALCULATION OF ZII BY EXTRAPOLATION IN THE RECTANGLE.
                    // CHECKS IF THE NECESSARY COEFFICIENTS HAVE BEEN CALCULATED.
                    if (IT0 == ITPV)
                    {
                        //GO TO 50
                    }
                    else
                    {
                        // LOADS COORDINATE AND PARTIAL DERIVATIVE VALUES AT THE END
                        // POINTS OF THE BORDER LINE SEGMENT.
                        int JIPL = 3 * (IL1 - 1);
                        int JPD = 0;

                        //initialie intermediate variables
                        int IDP;
                        int JPDD;

                        //DO 43  I = 1,2
                        JIPL = JIPL + 1;
                        IDP = IPL[JIPL];
                        X[I] = XD[IDP];
                        Y[I] = YD[IDP];
                        Z[I] = ZD[IDP];
                        JPDD = 5 * (IDP - 1);
                        //DO 42  KPD = 1,5
                        JPD = JPD + 1;
                        JPDD = JPDD + 1;
                        PD[JPD] = PDD[JPDD];
                        //  CONTINUE
                        // CONTINUE

                        // DETERMINES THE COEFFICIENTS FOR THE COORDINATE SYSTEM
                        // TRANSFORMATION FROM THE X-Y SYSTEM TO THE U-V SYSTEM
                        // AND VICE VERSA.
                        X0 = X(1);
      Y0 = Y(1)
      A = Y(2) - Y(1)
      B = X(2) - X(1)
      C = -B
      D = A
      AD = A * D
      BC = B * C
      DLT = AD - BC
      AP = D / DLT
      BP = -B / DLT
      CP = -BP
      DP = AP
C CONVERTS THE PARTIAL DERIVATIVES AT THE END POINTS OF THE
C BORDER LINE SEGMENT FOR THE U-V COORDINATE SYSTEM.
   45 AA = A * A
      ACT2 = 2.0 * A * C
      CC = C * C
      AB = A * B
      ADBC = AD + BC
      CD = C * D
      BB = B * B
      BDT2 = 2.0 * B * D
      DD = D * D
      DO 46  I = 1,2
        JPD = 5 * I
        ZU(I) = A * PD(JPD - 4) + C * PD(JPD - 3)
        ZV(I) = B * PD(JPD - 4) + D * PD(JPD - 3)
        ZUU(I) = AA * PD(JPD - 2) + ACT2 * PD(JPD - 1) + CC * PD(JPD)
        ZUV(I) = AB * PD(JPD - 2) + ADBC * PD(JPD - 1) + CD * PD(JPD)
        ZVV(I) = BB * PD(JPD - 2) + BDT2 * PD(JPD - 1) + DD * PD(JPD)
   46 CONTINUE
C CALCULATES THE COEFFICIENTS OF THE POLYNOMIAL.
   47 P00 = Z(1)
      P10 = ZU(1)
      P01 = ZV(1)
      P20 = 0.5 * ZUU(1)
      P11 = ZUV(1)
      P02 = 0.5 * ZVV(1)
      H1 = Z(2) - P00 - P01 - P02
      H2 = ZV(2) - P01 - ZVV(1)
      H3 = ZVV(2) - ZVV(1)
      P03 = 10.0 * H1 - 4.0 * H2 + 0.5 * H3
      P04 = -15.0 * H1 + 7.0 * H2 - H3
      P05 = 6.0 * H1 - 3.0 * H2 + 0.5 * H3
      H1 = ZU(2) - P10 - P11
      H2 = ZUV(2) - P11
      P12 = 3.0 * H1 - H2
      P13 = -2.0 * H1 + H2
      P21 = 0.0
      P23 = -ZUU(2) + ZUU(1)
      P22 = -1.5 * P23
      ITPV = IT0
C CONVERTS XII AND YII TO U - V SYSTEM.
   50 DX = XII - X0
      DY = YII - Y0
      U = AP * DX + BP * DY
      V = CP * DX + DP * DY
C EVALUATES THE POLYNOMIAL.
   51 P0 = P00 + V * (P01 + V * (P02 + V * (P03 + V * (P04 + V * P05))))
      P1 = P10 + V * (P11 + V * (P12 + V * P13))
      P2 = P20 + V * (P21 + V * (P22 + V * P23))
      ZII = P0 + U * (P1 + U * P2)
      RETURN
                    }
                } else {
                    //GO TO 60
                }
            };


        }

        public List<double> seq(double min, double max, int length)
        {
            List<double> result = new List<double>();

            double by = (max - min) / (length -1 );
            int n = Convert.ToInt32((max-min) / by);
            for (int i=0; i <= n; i++)
            {
                double val = min + i * by;
                result.Add(val);
            }

            return result;
        }
    }
}
