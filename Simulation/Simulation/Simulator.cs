using System;
using System.IO;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using System.Xml.Serialization;

namespace Simulation
{
    [Serializable]
    public class Simulator
    {
        [Serializable]
        [StructLayout(LayoutKind.Sequential)]
        public struct Node
        {
            public double U;
            public double V;

            public override string ToString() => $"({U}, {V})";
        }
        [CLIParameter("N")]
        public int NodeCount = 100;
        [CLIParameter("P")]
        public int CouplingRadius = 1;

        [CLIParameter("sigma")]
        public double CouplingStrength = 0.1;
        [CLIParameter("epsilon")]
        public double TimeScale = 0.01;
        [CLIParameter("D")]
        public double DiffusionConstant = 0.5;

        [CLIParameter("dt")]
        public double DeltaTime = 0.001;

        [CLIParameter("t0")]
        public double WarmupTime = 2000.0;
        int WarmupCount;
        [CLIParameter("dtm")]
        public double MeasureInterval = 0.1;
        int IntervalSize;
        [CLIParameter("t1")]
        public double MaximumTime = 2060.0;
        int Measurements;

        [CLIParameter("tau")]
        public double DelayTime = 0.0;
        int DelaySteps;

        [CLIParameter("s")]
        public int Seed = (int)DateTime.Now.Ticks;
        Random[] RNG;

        [CLIParameter("a0")]
        public double BifurcationParameterMean = 1.05;
        [CLIParameter("da")]
        public double BifurcationParameterVariance = 0;
        double[] bifurcationParameter;

        double CouplingConstant;
        double NoiseAmplitude;

        Node[] SwapBuffer;
        Node[][] DelayBuffer;
        Node[][] Data;

        bool hasRun = false;
        bool createdOutput = false;
        Guid guid = Guid.NewGuid();

        public string OutputDirectory => createdOutput ? guid.ToString() : null;

        int Cores;

        double EnsurePhase(double x) => x <= 1.0 ? 2.0 - x : x;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        int Mod(int x) => x < 0 ? x + NodeCount : x >= NodeCount ? x - NodeCount : x;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        double Cube(double x) => x * x * x;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        double RandomNormal(int core) => Math.Cos(2.0 * Math.PI * RNG[core].NextDouble()) * Math.Sqrt(-2.0 * Math.Log(1.0 - RNG[core].NextDouble()));

        void Init2DArray(ref Node[][] nodes, int length)
        {
            nodes = new Node[length][];
            for (int i = 0; i < length; ++i)
                nodes[i] = new Node[NodeCount];
        }

        void Copy2DArray(ref Node[][] source, ref Node[][] destination)
        {
            for (int i = 0; i < Measurements; ++i)
                Array.Copy(source[i], destination[i], NodeCount);
        }

        void Setup()
        {
            Cores = Environment.ProcessorCount;

            WarmupCount = Convert.ToInt32(WarmupTime / DeltaTime);
            IntervalSize = Convert.ToInt32(MeasureInterval / DeltaTime);
            Measurements = Convert.ToInt32((MaximumTime - WarmupTime) / MeasureInterval);

            if (Measurements <= 0)
                throw new ArgumentOutOfRangeException("MaximumTime should be bigger then WarmupTime");

            DelaySteps = Convert.ToInt32(Math.Ceiling(DelayTime / DeltaTime));

            RNG = new Random[Cores];
            for (int i = 0; i < Cores; ++i)
                RNG[i] = new Random(Seed + i);

            bifurcationParameter = new double[NodeCount];
            for (int i = 0; i < NodeCount; ++i)
                bifurcationParameter[i] = EnsurePhase(BifurcationParameterMean + RandomNormal(i % Cores) * BifurcationParameterVariance);

            CouplingConstant = CouplingStrength / (2.0 * CouplingRadius);
            NoiseAmplitude = Math.Sqrt(2.0 * DiffusionConstant);

            SwapBuffer = new Node[NodeCount];
            Init2DArray(ref DelayBuffer, DelaySteps + 1);
            Init2DArray(ref Data, Measurements + 1);

        }

        void Run()
        {
            for (int i = 0; i < WarmupCount; ++i)
                Step();

            Array.Copy(DelayBuffer[DelaySteps], Data[0], NodeCount);
            for (int i = 1; i <= Measurements; ++i)
            {
                for (int j = 0; j < IntervalSize; ++j)
                    Step();

                Array.Copy(DelayBuffer[DelaySteps], Data[i], NodeCount);
            }

            hasRun = true;
        }

        void Step()
        {
            Parallel.For(0, Cores, delegate (int core)
            {
                for (int i = core; i < NodeCount; i += Cores)
                {
                    double coupling = 0.0;
                    int maxJ = i + CouplingRadius;
                    for (int j = i - CouplingRadius; j <= maxJ; ++j)
                        coupling += DelayBuffer[0][Mod(j)].U - DelayBuffer[DelaySteps][i].U;

                    SwapBuffer[i].U = DelayBuffer[DelaySteps][i].U + DeltaTime * (
                        DelayBuffer[DelaySteps][i].U - Cube(DelayBuffer[DelaySteps][i].U) / 3.0 - DelayBuffer[DelaySteps][i].V + coupling * CouplingConstant) / TimeScale;
                    SwapBuffer[i].V = DelayBuffer[DelaySteps][i].V + DeltaTime * (
                        DelayBuffer[DelaySteps][i].U + bifurcationParameter[i] + RandomNormal(core) * NoiseAmplitude);
                }
            });

            for (int i = 0; i < DelaySteps; ++i)
                Array.Copy(DelayBuffer[i + 1], DelayBuffer[i], NodeCount);
            Array.Copy(SwapBuffer, DelayBuffer[DelaySteps], NodeCount);
        }

        public void GenerateTrajectory(ref Node[][] data)
        {
            Setup();
            Run();

            Copy2DArray(ref Data, ref data);
        }

        string MakeDirectory()
        {
            var dir = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Data", guid.ToString());
            Directory.CreateDirectory(dir);

            createdOutput = true;

            return dir;
        }

        void EnsureRun()
        {
            if (!hasRun)
            {
                Setup();
                Run();
            }
        }

        [CLIAction("cfg")]
        public void Config()
        {
            EnsureRun();

            var dir = MakeDirectory();

            var xml = new XmlSerializer(typeof(Simulator));
            using (var config = File.Create(Path.Combine(dir, "Config.xml")))
                xml.Serialize(config, this);
        }

        [CLIAction("traj")]
        public void Trajectory()
        {
            EnsureRun();

            var dir = MakeDirectory();

            using (var data = File.Create(Path.Combine(dir, "Data.bin")))
            using (var writer = new BinaryWriter(data))
                for (int i = 0; i <= Measurements; ++i)
                    for (int j = 0; j < NodeCount; ++j)
                    {
                        writer.Write(Data[i][j].U);
                        // writer.Write(Data[i][j].V);
                    }
        }
    }
}
