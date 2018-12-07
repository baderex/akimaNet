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

            double LU;
            double LV;

            // ?? initializing output
            MISSII = false;
            ZII = double.NaN;

            // ??: value by reference globaly
            //COMMON / IDPI / ITPV
            int ITPV = 0; // ?? create a public property 

            // ??: value by reference locally
            //EQUIVALENCE(P5, P50)
            double P5;
            double P50;

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
                {
                    // LOADS COORDINATE AND PARTIAL DERIVATIVE VALUES AT THE
                    // VERTEXES.
                    int JIPT = 3 * (IT0 - 1);
                    int JPD = 0;

                    //DO I = 1,3
                    for (int I = 0; I < 3; I++)
                    {
                        JIPT = JIPT + 1;
                        int IDP = IPT[JIPT];
                        X[I] = XD[IDP];
                        Y[I] = YD[IDP];
                        Z[I] = ZD[IDP];
                        int JPDD = 5 * (IDP - 1);
                        //DO  KPD = 1,5
                        for (int KPD = 0; KPD < 3; KPD++)
                        {
                            JPD = JPD + 1;
                            JPDD = JPDD + 1;
                            PD[JPD] = PDD[JPDD];
                        }
                    }

                    // DETERMINES THE COEFFICIENTS FOR THE COORDINATE SYSTEM
                    // TRANSFORMATION FROM THE X-Y SYSTEM TO THE U-V SYSTEM
                    // AND VICE VERSA.
                    double X0 = X[1];
                    double Y0 = Y[1];
                    double A = X[2] - X0;
                    double B = X[3] - X0;
                    double C = Y[2] - Y0;
                    double D = Y[3] - Y0;
                    double AD = A * D;
                    double BC = B * C;
                    double DLT = AD - BC;
                    double AP = D / DLT;
                    double BP = -B / DLT;
                    double CP = -C / DLT;
                    double DP = A / DLT;

                    // CONVERTS THE PARTIAL DERIVATIVES AT THE VERTEXES OF THE
                    // TRIANGLE FOR THE U - V COORDINATE SYSTEM.
                    double AA = A * A;
                    double ACT2 = 2.0 * A * C;
                    double CC = C * C;
                    double AB = A * B;
                    double ADBC = AD + BC;
                    double CD = C * D;
                    double BB = B * B;
                    double BDT2 = 2.0 * B * D;
                    double DD = D * D;

                    // ?? 1-index issue
                    //DO 26  I = 1,3
                    for (int I = 1; I < 3; I++)
                    {
                        JPD = 5 * I;
                        ZU[I - 1] = A * PD[JPD - 4 - 1] + C * PD[JPD - 3 - 1];
                        ZV[I - 1] = B * PD[JPD - 4 - 1] + D * PD[JPD - 3 - 1];
                        ZUU[I - 1] = AA * PD[JPD - 2 - 1] + ACT2 * PD[JPD - 1 - 1] + CC * PD[JPD - 1];
                        ZUV[I - 1] = AB * PD[JPD - 2 - 1] + ADBC * PD[JPD - 1 - 1] + CD * PD[JPD - 1];
                        ZVV[I - 1] = BB * PD[JPD - 2 - 1] + BDT2 * PD[JPD - 1 - 1] + DD * PD[JPD - 1];
                    }

                    // CALCULATES THE COEFFICIENTS OF THE POLYNOMIAL.
                    double P00 = Z[1];
                    double P10 = ZU[1];
                    double P01 = ZV[1];
                    double P20 = 0.5 * ZUU[1];
                    double P11 = ZUV[1];
                    double P02 = 0.5 * ZVV[1];
                    double H1 = Z[2] - P00 - P10 - P20;
                    double H2 = ZU[2] - P10 - ZUU[1];
                    double H3 = ZUU[2] - ZUU[1];
                    double P30 = 10.0 * H1 - 4.0 * H2 + 0.5 * H3;
                    double P40 = -15.0 * H1 + 7.0 * H2 - H3;
                    P50 = 6.0 * H1 - 3.0 * H2 + 0.5 * H3;
                    H1 = Z[3] - P00 - P01 - P02;
                    H2 = ZV[3] - P01 - ZVV[1];
                    H3 = ZVV[3] - ZVV[1];
                    double P03 = 10.0 * H1 - 4.0 * H2 + 0.5 * H3;
                    double P04 = -15.0 * H1 + 7.0 * H2 - H3;
                    double P05 = 6.0 * H1 - 3.0 * H2 + 0.5 * H3;
                    LU = Math.Sqrt(AA + CC);
                    LV = Math.Sqrt(BB + DD);
                    double THXU = Math.Atan2(C, A);
                    double THUV = Math.Atan2(D, B) - THXU;
                    double CSUV = Math.Cos(THUV);
                    double P41 = 5.0 * LV * CSUV / LU * P50;
                    double P14 = 5.0 * LU * CSUV / LV * P05;
                    H1 = ZV[2] - P01 - P11 - P41;
                    H2 = ZUV[2] - P11 - 4.0 * P41;
                    double P21 = 3.0 * H1 - H2;
                    double P31 = -2.0 * H1 + H2;
                    H1 = ZU[3] - P10 - P11 - P14;
                    H2 = ZUV[3] - P11 - 4.0 * P14;
                    double P12 = 3.0 * H1 - H2;
                    double P13 = -2.0 * H1 + H2;
                    double THUS = Math.Atan2(D - C, B - A) - THXU;
                    double THSV = THUV - THUS;
                    AA = Math.Sin(THSV) / LU;
                    BB = -Math.Cos(THSV) / LU;
                    CC = Math.Sin(THUS) / LV;
                    DD = Math.Cos(THUS) / LV;
                    double AC = AA * CC;
                    AD = AA * DD;
                    BC = BB * CC;
                    double G1 = AA * AC * (3.0 * BC + 2.0 * AD);
                    double G2 = CC * AC * (3.0 * AD + 2.0 * BC);
                    H1 = -AA * AA * AA * (5.0 * AA * BB * P50 + [4.0 * BC + AD] * P41) - CC * CC * CC * (5.0 * CC * DD * P05 + (4.0 * AD + BC) * P14);
                    H2 = 0.5 * ZVV[2] - P02 - P12;
                    H3 = 0.5 * ZUU[3] - P20 - P21;
                    double P22 = (G1 * H2 + G2 * H3 - H1) / (G1 + G2);
                    double P32 = H2 - P22;
                    double P23 = H3 - P22;
                    ITPV = IT0;

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

                };

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
                        // CONVERTS XII AND YII TO U - V SYSTEM.
                        double DX = XII - X0;
                        double DY = YII - Y0;
                        double U = AP * DX + BP * DY;
                        double V = CP * DX + DP * DY;

                        // EVALUATES THE POLYNOMIAL.
                        double P0 = P00 + V * (P01 + V * (P02 + V * (P03 + V * (P04 + V * P05))));
                        double P1 = P10 + V * (P11 + V * (P12 + V * P13));
                        double P2 = P20 + V * (P21 + V * (P22 + V * P23));
                        ZII = P0 + U * (P1 + U * P2);
                        return;
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

                        // ?? 1-index issue
                        //DO  I = 1,2
                        for (int I = 0; I < 2; I++)
                        {
                            JIPL = JIPL + 1;
                            IDP = IPL[JIPL];
                            X[I] = XD[IDP];
                            Y[I] = YD[IDP];
                            Z[I] = ZD[IDP];
                            JPDD = 5 * (IDP - 1);
                            //DO  KPD = 1,5
                            for (int KPD = 0; KPD < 5; KPD++) {
                                JPD = JPD + 1;
                                JPDD = JPDD + 1;
                                PD[JPD] = PDD[JPDD];
                            }
                        }

                        // DETERMINES THE COEFFICIENTS FOR THE COORDINATE SYSTEM
                        // TRANSFORMATION FROM THE X-Y SYSTEM TO THE U-V SYSTEM
                        // AND VICE VERSA.
                        double X0 = X[1];
                        double Y0 = Y[1];
                        double A = Y[2] - Y[1];
                        double B = X[2] - X[1];
                        double C = -B;
                        double D = A;
                        double AD = A * D;
                        double BC = B * C;
                        double DLT = AD - BC;
                        double AP = D / DLT;
                        double BP = -B / DLT;
                        double CP = -BP;
                        double DP = AP;

                        // CONVERTS THE PARTIAL DERIVATIVES AT THE END POINTS OF THE
                        // BORDER LINE SEGMENT FOR THE U-V COORDINATE SYSTEM.
                        double AA = A * A;
                        double ACT2 = 2.0 * A * C;
                        double CC = C * C;
                        double AB = A * B;
                        double ADBC = AD + BC;
                        double CD = C * D;
                        double BB = B * B;
                        double BDT2 = 2.0 * B * D;
                        double DD = D * D;

                        // ?? chnage from 1-index to 0-index
                        //DO  I=1,3
                        for (int I = 1; I < 4; I++)
                        {
                            JPD = 5 * I;
                            ZU[I - 1] = A * PD[JPD - 4 - 1] + C * PD[JPD - 3 - 1];
                            ZV[I - 1] = B * PD[JPD - 4 - 1] + D * PD[JPD - 3 - 1];
                            ZUU[I - 1] = AA * PD[JPD - 2 - 1] + ACT2 * PD[JPD - 1 - 1] + CC * PD[JPD - 1];
                            ZUV[I - 1] = AB * PD[JPD - 2 - 1] + ADBC * PD[JPD - 1 - 1] + CD * PD[JPD - 1];
                            ZVV[I - 1] = BB * PD[JPD - 2 - 1] + BDT2 * PD[JPD - 1 - 1] + DD * PD[JPD - 1];
                        }

                        // CALCULATES THE COEFFICIENTS OF THE POLYNOMIAL.
                        double P00 = Z[1];
                        double P10 = ZU[1];
                        double P01 = ZV[1];
                        double P20 = 0.5 * ZUU[1];
                        double P11 = ZUV[1];
                        double P02 = 0.5 * ZVV[1];
                        double H1 = Z[2] - P00 - P01 - P02;
                        double H2 = ZV[2] - P01 - ZVV[1];
                        double H3 = ZVV[2] - ZVV[1];
                        double P03 = 10.0 * H1 - 4.0 * H2 + 0.5 * H3;
                        double P04 = -15.0 * H1 + 7.0 * H2 - H3;
                        double P05 = 6.0 * H1 - 3.0 * H2 + 0.5 * H3;
                        H1 = ZU[2] - P10 - P11;
                        H2 = ZUV[2] - P11;
                        double P12 = 3.0 * H1 - H2;
                        double P13 = -2.0 * H1 + H2;
                        double P21 = 0.0;
                        double P23 = -ZUU[2] + ZUU[1];
                        double P22 = -1.5 * P23;
                        ITPV = IT0;

                        // CONVERTS XII AND YII TO U - V SYSTEM.
                        double DX = XII - X0;
                        double DY = YII - Y0;
                        double U = AP * DX + BP * DY;
                        double V = CP * DX + DP * DY;

                        // EVALUATES THE POLYNOMIAL.
                        double P0 = P00 + V * (P01 + V * (P02 + V * (P03 + V * (P04 + V * P05))));
                        double P1 = P10 + V * (P11 + V * (P12 + V * P13));
                        double P2 = P20 + V * (P21 + V * (P22 + V * P23));
                        ZII = P0 + U * (P1 + U * P2);
                        return;
                    }
                } else {
                    // CALCULATION OF ZII BY EXTRAPOLATION IN THE TRIANGLE.
                    // CHECKS IF THE NECESSARY COEFFICIENTS HAVE BEEN CALCULATED.
                    if (IT0 == ITPV)
                    {
                        // GO TO 70 !!
                        // CONVERTS XII AND YII TO U - V SYSTEM.
                        double U = XII - X[1];
                        double V = YII - Y[1];

                        // EVALUATES THE POLYNOMIAL.
                        double P0 = P00 + V * (P01 + V * P02);
                        double P1 = P10 + V * P11;
                        ZII = P0 + U * (P1 + U * P20);
                        return;
                    }
                    else
                    {
                        // LOADS COORDINATE AND PARTIAL DERIVATIVE VALUES AT THE VERTEX
                        // OF THE TRIANGLE.
                        int JIPL = 3 * IL2 - 2;
                        int IDP = IPL[JIPL];
                        X[1] = XD[IDP];
                        Y[1] = YD[IDP];
                        Z[1] = ZD[IDP];
                        int JPDD = 5 * (IDP - 1);

                        //DO KPD = 1,5
                        for (int KPD = 0; KPD < 5; KPD++)
                        {
                            JPDD = JPDD + 1;
                            PD[KPD] = PDD[JPDD];
                        }

                        // CALCULATES THE COEFFICIENTS OF THE POLYNOMIAL.
                        double P00 = Z[1];
                        double P10 = PD[1];
                        double P01 = PD[2];
                        double P20 = 0.5 * PD[3];
                        double P11 = PD[4];
                        double P02 = 0.5 * PD[5];
                        ITPV = IT0;

                        // CONVERTS XII AND YII TO U - V SYSTEM.
                        double U = XII - X[1];
                        double V = YII - Y[1];

                        // EVALUATES THE POLYNOMIAL.
                        double P0 = P00 + V * (P01 + V * P02);
                        double P1 = P10 + V * P11;
                        ZII = P0 + U * (P1 + U * P20);
                        return;
                    }
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
