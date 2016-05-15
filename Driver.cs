﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    class Driver
    {
        static void Main()
        {
            int N = 32;
            double dt = 0.05;
            double nu = 1;

            double tf = 1;
            double t = 0;
        
            double[] u0 = new double[(int)Math.Pow(N, 3)];
            double[] v0 = new double[(int)Math.Pow(N, 3)];
            double[] w0 = new double[(int)Math.Pow(N, 3)];

            CavityDomain omega = new CavityDomain(N);
            FluidSolver ffd = new FluidSolver(omega, dt, nu, u0, v0, w0);

            while (t < tf)
            {                
                t += dt;
                ffd.time_step();

                ffd.export_vtk(String.Concat("lid_driven_cavity_", t, ".vtk"));
                Console.WriteLine("Time t = {0}", t);
            }
        }
    }
}