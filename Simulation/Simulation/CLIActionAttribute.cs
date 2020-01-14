using System;

namespace Simulation
{
    [AttributeUsage(AttributeTargets.Method)]
    class CLIActionAttribute : Attribute
    {
        public string Shorthand { get; set; }

        public CLIActionAttribute(string shorthand)
        {
            Shorthand = shorthand;
        }
    }
}
