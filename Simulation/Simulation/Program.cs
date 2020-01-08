using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Simulation
{
    class Program
    {
        // TODO: set parameters via command line
        static void Main(string[] args)
        {
            var sim = new Simulator();
            Console.Write(sim.SaveTrajectory().ToString());
        }
    }
}
