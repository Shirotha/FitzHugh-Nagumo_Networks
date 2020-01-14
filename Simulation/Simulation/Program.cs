using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Simulation
{
    class Program
    {
        static void Main(string[] args)
        {
            var sim = new Simulator();
            var cli = new CLI(sim);

            cli.Execute(args);

            Console.Write(sim.OutputDirectory);
        }
    }
}
