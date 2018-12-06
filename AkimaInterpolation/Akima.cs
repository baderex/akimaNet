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


        public List<double> interpp(List<double> _x, List<double> _y, List<double> _z,String duplicate = "error")
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
