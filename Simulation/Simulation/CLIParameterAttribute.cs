using System;

namespace Simulation
{
    [AttributeUsage(AttributeTargets.Field)]
    class CLIParameterAttribute : Attribute
    {
        public string Shorthand { get; set; }

        public CLIParameterAttribute(string shorthand)
        {
            Shorthand = shorthand;
        }
    }
}
