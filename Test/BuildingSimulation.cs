using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FastFluidSolverMT;

/*
 * BuildingSimulation.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 * modified 2018 by Christoph Waibel <chwaibel@student.ethz.ch>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */
 
namespace Test
{
    /// <summary>
    /// Driver to simulate wind flow around 3 buildings.
    /// </summary>
    internal class BuildingSimulation
    {
        internal static void Run()
        {
            // Set mesh parameters, here we ask for Nx cells in the x direction
            // Ny cells in the y direction and Nz cells in the z direction (ignoring ghost cells)
            int Nx = 100;
            int Ny = 45;
            int Nz = 30;

            // Set time step and viscosity
            double dt = 0.5;
            double nu = 1e-3;

            // Set simulation start and finish times
            double tf = 100;
            double t = 0;

            // Set initial conditions
            double[, ,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];

            // Create empty arrays for body forces
            double[, ,] f_x = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] f_y = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] f_z = new double[Nx + 2, Ny + 2, Nz + 1];

            // Create structure for solver parameters
            FluidSolver.solver_struct solver_params = new FluidSolver.solver_struct();

            solver_params.tol = 1e-4;
            solver_params.min_iter = 1;
            solver_params.max_iter = 100;
            solver_params.verbose = false;
            solver_params.backtrace_order = 2;

            // Create domain
            WindInflow omega = new WindInflow(Nx + 2, Ny + 2, Nz + 2, 100, 45, 30);

            // Add buildings
            omega.add_obstacle(15, 20, 15, 20, 0, 5);
            omega.add_obstacle(15, 18, 18, 20, 5, 12);

            omega.add_obstacle(18, 25, 25, 30, 0, 8);
            omega.add_obstacle(30, 40, 15, 25, 0, 10);

            // Create FFD solver
            FluidSolver ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_params);

            // Create post processor and export initial conditions and geometry information
            PostProcessor pp = new PostProcessor(ffd, omega);
            
            int tstep = 0;
            //pp.export_data_vtk(String.Concat("city_test_", tstep, ".vtk"), 0, false);
            //pp.export_geometry_vtk("city_test_geometry.vtk", 0);

            // Run time loop
            while (t < tf)
            {
                t += dt;
                tstep++;

                Console.WriteLine("Time t = {0}", t);

                double[,,] p_t2 = new double[ffd.p.GetLength(0), ffd.p.GetLength(1), ffd.p.GetLength(2)];
                Array.Copy(ffd.p, 0, p_t2, 0, ffd.p.Length);
                double[,,] u_t2 = new double[ffd.u.GetLength(0), ffd.u.GetLength(1), ffd.u.GetLength(2)];
                Array.Copy(ffd.u, 0, u_t2, 0, ffd.u.Length);
                double[,,] v_t2 = new double[ffd.v.GetLength(0), ffd.v.GetLength(1), ffd.v.GetLength(2)];
                Array.Copy(ffd.v, 0, v_t2, 0, ffd.v.Length);
                double[,,] w_t2 = new double[ffd.w.GetLength(0), ffd.w.GetLength(1), ffd.w.GetLength(2)];
                Array.Copy(ffd.w, 0, w_t2, 0, ffd.w.Length);

                // Solve single time step and export results
                ffd.time_step(f_x, f_y, f_z);
                //pp.export_data_vtk(String.Concat("city_test_", tstep, ".vtk"), t, false);

                double[] p_residuals;
                double[,,] p_t1 = ffd.p;
                Utilities.calculate_residuals(p_t1, p_t2, out p_residuals);
                Console.WriteLine("p residuals: {0};{1};{2}", p_residuals[0], p_residuals[1], p_residuals[2]);
                double[] u_residuals;
                double[,,] u_t1 = ffd.u;
                Utilities.calculate_residuals(u_t1, u_t2, out u_residuals);
                Console.WriteLine("u residuals: {0};{1};{2}", u_residuals[0], u_residuals[1], u_residuals[2]);
                double[] v_residuals;
                double[,,] v_t1 = ffd.v;
                Utilities.calculate_residuals(v_t1, v_t2, out v_residuals);
                Console.WriteLine("v residuals: {0};{1};{2}", v_residuals[0], v_residuals[1], v_residuals[2]);
                double[] w_residuals;
                double[,,] w_t1 = ffd.w;
                Utilities.calculate_residuals(w_t1, w_t2, out w_residuals);
                Console.WriteLine("w residuals: {0};{1};{2}", w_residuals[0], w_residuals[1], w_residuals[2]);



                if (t % 1 == 1) Console.ReadKey();
            }
        }
    }
}
