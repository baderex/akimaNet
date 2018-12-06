using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AkimaInterpolation
{
    class Program
    {
        static void Main(string[] args)
        {
            Akima a = new Akima();

            List<double> z = a.seq(5, 100, 4);

            Console.WriteLine("Z Count = " + z.Count);
            Console.WriteLine(String.Join(", ", z));

            Console.ReadLine();
        }
    }
}
