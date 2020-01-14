using System;
using System.Collections.Generic;
using System.Reflection;

namespace Simulation
{
    public class CLI
    {
        class Parameter
        {
            public FieldInfo Field;

            public string Name;
            public string Shorthand;

            public bool IsMatch(string input) => input.Equals($"--{Name}", StringComparison.OrdinalIgnoreCase) || input.Equals($"-{Shorthand}", StringComparison.OrdinalIgnoreCase);

            public void Set(object target, object value) => Field.SetValue(target, value);

            public void Set(object target, string value) => Field.SetValue(target, Convert.ChangeType(value, Field.FieldType));
        }

        class Action
        {
            public MethodInfo Method;

            public string Name;
            public string Shorthand;

            public bool IsMatch(string input) => input.Equals(Name, StringComparison.OrdinalIgnoreCase) || input.Equals(Shorthand, StringComparison.OrdinalIgnoreCase);

            public void Invoke(object target) => Method.Invoke(target, null);
        }

        object target;

        List<Parameter> parameters = new List<Parameter>();
        List<Action> actions = new List<Action>();

        public CLI(object target)
        {
            this.target = target;

            var type = target.GetType();
            foreach (var field in type.GetFields())
                if (field.GetCustomAttribute(typeof(CLIParameterAttribute)) is CLIParameterAttribute attr)
                    parameters.Add(new Parameter
                    {
                        Field = field,
                        Name = field.Name,
                        Shorthand = attr.Shorthand
                    });

            foreach (var method in type.GetMethods())
                if (method.GetCustomAttribute(typeof(CLIActionAttribute)) is CLIActionAttribute attr)
                    actions.Add(new Action
                    {
                        Method = method,
                        Name = method.Name,
                        Shorthand = attr.Shorthand
                    });
        }

        public void Execute(string[] args)
        {
            bool handled;
            for (int i = 0; i < args.Length; ++i)
            {
                handled = false;
                foreach (var p in parameters)
                    if (p.IsMatch(args[i]))
                    {
                        p.Set(target, args[++i]);
                        handled = true;
                        break;
                    }
                if (handled)
                    continue;

                foreach (var a in actions)
                    if (a.IsMatch(args[i]))
                    {
                        a.Invoke(target);
                        handled = true;
                        break;
                    }
            }
        }
    }
}
