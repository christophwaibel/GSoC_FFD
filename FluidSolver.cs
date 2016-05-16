﻿using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace FastFluidSolver
{
    class FluidSolver
    {
        const int MAX_ITER = 50; //maximum number of iterations for Gauss-Seidel solver
        const double TOL = 1e-5; //maximum relative error for Gauss-Seidel solver

        private double[] u; // x component of velocity
        private double[] v; // y component of velocity
        private double[] w; // z component of velocity
        private double[] p; // pressure

        private double[] u_sratch; //scratch arrays for velocities
        private double[] v_sratch;
        private double[] w_sratch;

        private double dt;  //time step
        private static int N; //number of points in each coordiate direction
        private double h;   //spacing in each corrdiate direction
        private double nu;  //fluid viscosity

        private bool verbose;

        Domain omega;

        void initialize() { }
        void add_force() { }

        /****************************************************************************
         * Constructor
         ****************************************************************************/
        public FluidSolver(Domain omega, double dt, double nu, double[] u0, double[] v0, double[] w0, bool verbose)
        {
            N = omega.N;
            h = omega.h;
            this.dt = dt;
            this.nu = nu;

            this.omega = omega;
            this.verbose = verbose;

            u = u0;
            v = v0;
            w = w0;

            p = new double[(int) Math.Pow(N, 3)];
        }
        /****************************************************************************
         * Diffusion step. Diffuse solve diffusion equation x_t = L(x) using second 
         * order finite difference in space and backwards Euler in time
         ****************************************************************************/
        void diffuse(ref double[] x)
        {
            double[] x_old = new double[x.Length];
            x.CopyTo(x_old, 0);

            gs_solve(1 + 6 * nu * dt / Math.Pow(h, 2), -dt * nu / Math.Pow(h, 2), x_old, ref x, 0);
        }

        /*****************************************************************************
         * Projection step. Solves Poisson equation L(p) = div(u_old) using second order finte
         * difference and updates the velocities, u = u_old - grad(p)
         ****************************************************************************/
        void project()
        {
            double[] div = new double[(int) Math.Pow(N,3)];

            // calculate div(w) using second order finite differences
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        if (omega.obstacle[cell_index(i, j, k, N)] == 0) //node not inside an obstacle
                        {
                            if (omega.boundary_nodes[cell_index(i, j, k, N)] == 0) //node not on boundary, use second order finite differnce
                            {
                                div[cell_index(i, j, k, N)] = (u[cell_index(i + 1, j, k, N)] - u[cell_index(i - 1, j, k, N)] +
                                    v[cell_index(i, j + 1, k, N)] - v[cell_index(i, j - 1, k, N)] + w[cell_index(i, j, k + 1, N)] +
                                    w[cell_index(i, j, k - 1, N)]) / (2 * h);
                            }
                            else //use first order finite difference
                            {
                                int nx = omega.boundary_normal_x[cell_index(i, j, k, N)];
                                int ny = omega.boundary_normal_x[cell_index(i, j, k, N)];
                                int nz = omega.boundary_normal_x[cell_index(i, j, k, N)];

                                //calculate each partial derivative individually
                                double ux = 0;
                                double vy = 0;
                                double wz = 0;

                                switch (nx)
                                {
                                    case 0:
                                        ux = (u[cell_index(i + 1, j, k, N)] - u[cell_index(i - 1, j, k, N)]) / (2 * h);
                                        break;

                                    case -1:
                                        ux = (u[cell_index(i + 1, j, k, N)] - u[cell_index(i, j, k, N)]) /h;
                                        break;

                                    case 1:
                                        ux = (u[cell_index(i, j, k, N)] - u[cell_index(i - 1, j, k, N)]) / h;
                                        break;
                                }

                                switch (ny)
                                {
                                    case 0:
                                        vy = (v[cell_index(i, j + 1, k, N)] - v[cell_index(i, j - 1, k, N)]) / (2 * h);
                                        break;

                                    case -1:
                                        vy = (v[cell_index(i, j + 1, k, N)] - v[cell_index(i, j, k, N)]) / h;
                                        break;

                                    case 1:
                                        vy = (v[cell_index(i, j, k, N)] - v[cell_index(i, j - 1, k, N)]) / h;
                                        break;
                                }

                                switch (nz)
                                {
                                    case 0:
                                        wz = (w[cell_index(i, j + 1, k, N)] - w[cell_index(i, j - 1, k, N)]) / (2 * h);
                                        break;

                                    case -1:
                                        wz = (w[cell_index(i, j, k + 1, N)] - w[cell_index(i, j, k, N)]) / h;
                                        break;

                                    case 1:
                                        wz = (w[cell_index(i, j, k, N)] - w[cell_index(i, j, k - 1, N)]) / h;
                                        break;
                                }

                                div[cell_index(i, j, k, N)] = ux + vy + wz;
                            }
                        }
                    }
                }
            }

            gs_solve(-6 / Math.Pow(h, 2), 1 / Math.Pow(h, 2), div, ref p, 1);

            //update velocity by adding calculate grad(p), calculated using second order finite difference
            //only need to add to interior points, as the velocity at boundary points has already been fixed
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        if (omega.obstacle[cell_index(i, j, k, N)] == 0) //node not inside an obstacle
                        {
                            if (omega.boundary_nodes[cell_index(i, j, k, N)] == 0) //node not on boundary, use second order finite differnce
                            {
                                u[cell_index(i, j, k, N)] -= (p[cell_index(i + 1, j, k, N)] - p[cell_index(i - 1, j, k, N)]) / (2 * h);
                                v[cell_index(i, j, k, N)] -= (p[cell_index(i, j + 1, k, N)] - p[cell_index(i, j - 1, k, N)]) / (2 * h);
                                w[cell_index(i, j, k, N)] -= (p[cell_index(i, j, k + 1, N)] - p[cell_index(i, j, k - 1, N)]) / (2 * h);
                            }
                        }
                    }
                }
            }
        }

        /***************************************************************************
         * Advection step. Uses a first order backtrace to update x.
         * 
         * @inputs
         * double[] x - reference to quantity to advect, can be a velocity, a
         * concentration or temperature
         * double[] x0 - initial state of x before advection
         * double[] velx - x velocity (either before convection if x is a velocity.
         * or after convection otherwise)
         * double[] vely - y velocity (either before convection if x is a velocity.
         * or after convection otherwise)
         * double[] velz - z velocity (either before convection if x is a velocity.
         * or after convection otherwise)
         **************************************************************************/
        void advect(ref double[] x, double[] x0, double[] velx, double[] vely,  double[] velz)
        {
        }

        /*****************************************************************************
         * Solves the banded system given by the finite difference method applied
         * to the Poisson or diffusion equation using the iterative Gauss-Seidel method.
         * @inputs
         * double a - coefficient along diagonal entry
         * double c - coefficient of all other nonzero entries
         * double[] b - right hand side
         * double[] x - reference to array in which to store solution
         * int boundary type - 0 corresponds to Dirichlet, 1 to homogeneous Neumann 
         ****************************************************************************/
        void gs_solve(double a, double c, double[] b, ref double[] x, int boundary_type)
        {
            int iter = 0;
            double res = 2 * TOL;
            while (iter < MAX_ITER && res > TOL)
            {
                res = 0;

                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            if (omega.obstacle[cell_index(i, j, k, N)] == 0) //if node not inside obstacle
                            {
                                double x_old = x[cell_index(i, j, k, N)];

                                if (omega.boundary_nodes[cell_index(i, j, k, N)] == 0) //if not on boundary, second order finite difference
                                {
                                    x[cell_index(i, j, k, N)] = (b[cell_index(i, j, k, N)] - c * (x[cell_index(i - 1, j, k, N)] +
                                            x[cell_index(i + 1, j, k, N)] + x[cell_index(i, j - 1, k, N)] + x[cell_index(i, j + 1, k, N)] +
                                            x[cell_index(i, j, k - 1, N)] + x[cell_index(i, j, k + 1, N)])) / a;
                                }
                                else if (boundary_type == 1)//if on boundary and homogeneous Neumann boundary conditions
                                {
                                    int nx = omega.boundary_normal_x[cell_index(i, j, k, N)];
                                    int ny = omega.boundary_normal_x[cell_index(i, j, k, N)];
                                    int nz = omega.boundary_normal_x[cell_index(i, j, k, N)];

                                    x[cell_index(i, j, k, N)] = (b[cell_index(i, j, k, N)] - c * (Math.Abs(1 + nx) * x[cell_index(i - 1, j, k, N)] +
                                            Math.Abs(1 - nx) * x[cell_index(i + 1, j, k, N)] + Math.Abs(1 + ny) * x[cell_index(i, j - 1, k, N)] +
                                            Math.Abs(1 - ny) * x[cell_index(i, j + 1, k, N)] + Math.Abs(1 + nz) * x[cell_index(i, j, k - 1, N)] +
                                            Math.Abs(1 - nz) * x[cell_index(i, j, k + 1, N)])) / a;

                                }

                                res += Math.Pow((x_old - x[cell_index(i, j, k, N)]), 2);
                            }
                        }
                    }
                }

                res = Math.Sqrt(res) / Math.Pow(N, 3);
                iter++;
            }

            if (verbose)
            {
                Console.WriteLine("Gauss-Seidel solver completed with residual of {0} in {1} iterations", res, iter);
            }
        }

        /*********************************************************************************
         * Applies the boundary conditions from the domain omega to the velocities
         ********************************************************************************/
        void apply_boundary_conditions()
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        if (omega.boundary_nodes[cell_index(i, j, k, N)] == 1)
                        {
                            u[cell_index(i, j, k, N)] = omega.boundary_u[cell_index(i, j, k, N)];
                            v[cell_index(i, j, k, N)] = omega.boundary_v[cell_index(i, j, k, N)];
                            w[cell_index(i, j, k, N)] = omega.boundary_w[cell_index(i, j, k, N)];
                        }
                    }
                }
            }
        }

        /***************************************************************************
         * Takes the x, y, z indices of a cell and returns the global coordinate
         **************************************************************************/
        public static int cell_index(int x, int y, int z, int N)
        {
            return (x >= 0 && y >= 0 && z >= 0 && x < N && y < N && z < N) ? x + y * N + z * (int)Math.Pow(N, 2) : 0;
        }

        /*****************************************************************************
         * Perform a single time step
         *****************************************************************************/
        public void time_step()
        {
            apply_boundary_conditions();

            diffuse(ref u);
            diffuse(ref v);
            diffuse(ref w);

            project();
        }

        /*****************************************************************************
         * Export data to a VTK file for visualization, based on file format guide here:
         * http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
         ****************************************************************************/
        public void export_vtk(String fname) 
        {
            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics data\n");
                sw.WriteLine("ASCII");
                /*sw.WriteLine("DATASET STRUCTURED_POINTS");
                sw.WriteLine("DIMENSIONS {0} {1} {2}", N, N, N);//TO DO change to accept different domain sizes
                sw.WriteLine("ORIGIN {0} {1} {2}", 0, 0, 0);
                sw.WriteLine("SPACING {0} {1} {2}", h, h, h);*/
                sw.WriteLine("DATASET STRUCTURED_GRID");
                sw.WriteLine("DIMENSIONS {0} {1} {2}", N, N, N);
                sw.WriteLine("POINTS {0} double", Math.Pow(N, 3));
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0} {1} {2}", h * i, h * j, h * k);
                        }
                    }
                }

                sw.WriteLine("POINT_DATA {0}", Math.Pow(N, 3));
                sw.WriteLine("VECTORS velocity double");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0} {1} {2}", u[cell_index(i, j, k, N)], v[cell_index(i, j, k, N)], w[cell_index(i, j, k, N)]);
                        }
                    }
                }

                sw.WriteLine("SCALARS pressure double {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0}", p[cell_index(i, j, k, N)]);
                        }
                    }
                }

                sw.WriteLine("SCALARS nx int {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0}", omega.boundary_normal_x[cell_index(i, j, k, N)]);
                        }
                    }
                }

                sw.WriteLine("SCALARS ny int {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0}", omega.boundary_normal_y[cell_index(i, j, k, N)]);
                        }
                    }
                }

                sw.WriteLine("SCALARS nz int {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0}", omega.boundary_normal_z[cell_index(i, j, k, N)]);
                        }
                    }
                }
            }
        }
    }
}
