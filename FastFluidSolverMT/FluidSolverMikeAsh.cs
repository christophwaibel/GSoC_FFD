using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolverMT
{
    public class FluidSolverMikeAsh
    {
        public int FFD_btnVal7;


        int G1, G2, G3;     //grid size. not x,y,z, coz evrth being put into 1-dim array

        int nCells;         //number to know the total amount of cells in the grid
        int[] solid;        // solids for Finite Solid Method. the index of the cell, either 0 or 1, 0 being no solid, 1 being solid
        int[] solidx;
        int[] solidy;
        int[] solidz;
        int iter = 5;

        float dt;                   //time step
        float diffusion;            //diffusion value
        float viscosity;            //viscosity value
        float densityIn;            //inflow density
        public float velocityIn;            //inflow velocity
        float[] density;            //density [] value, per cell
        float[] densityScratch; //density scratch []. scratch, so things can be changes. sth like a copy

        public float[] Vx;					//velocity [] in direction X
        public float[] Vy;					//velocity [] in direction Y
        public float[] Vz;                  //velocity [] in direction Z
        float[] Vx0;                //velocity  scratch [] in X
        float[] Vy0;                //velocity  scratch [] in Y
        float[] Vz0;                //velocity  scratch [] in Z


        float[] loads;				//loads per cell. only important, if it is a solid


        public FluidSolverMikeAsh(int GridSizeX, int GridSizeY, int GridSizeZ)
        {
            G1 = GridSizeX;
            G2 = GridSizeY;
            G3 = GridSizeZ;

            int nCellsX, nCellsY, nCellsZ;
            nCellsX = G1;
            nCellsY = G2;
            nCellsZ = G3;
            nCells = nCellsX * nCellsY * nCellsZ;

            density = new float[nCells];
            densityScratch = new float[nCells];
            solid = new int[nCells];
            solidx = new int[nCells];
            solidy = new int[nCells];
            solidz = new int[nCells];
            Vx = new float[nCells];
            Vy = new float[nCells];
            Vz = new float[nCells];
            Vx0 = new float[nCells];
            Vy0 = new float[nCells];
            Vz0 = new float[nCells];

            loads = new float[nCells];

            for (int i = 0; i < nCells; i++)
            {
                Vx[i] = 0f;
                Vy[i] = 0f;
                Vz[i] = 0f;
                Vx0[i] = 0f;
                Vy0[i] = 0f;
                Vz0[i] = 0f;
                density[i] = 0.0f;
                densityScratch[i] = 0.0f;
                solid[i] = 0;
                solidx[i] = 0;
                solidy[i] = 0;
                solidz[i] = 0;

                loads[i] = 0f;
            }

            for (int i = 0; i < G1; i++)
            {
                for (int j = 0; j < G2; j++)
                {
                    for (int k = 0; k < G3; k++)
                    {
                        solidx[IndexCell(i, j, k)] = i;
                        solidy[IndexCell(i, j, k)] = j;
                        solidz[IndexCell(i, j, k)] = k;
                    }
                }
            }

        }



        /*
		void drawAllIntBnd(){
			int index;
			float col;
			PVector vec;
			pushStyle();
			for(int k=0; k<G3; k++){
				for (int j=0; j<G2; j++){
					for (int i=0; i<G1; i++){
						index = IndexCell(i,j,k);
						if (solid[index]==1){

							
							pushMatrix();
							translate(i*FFD_ScaleXYZ,j*FFD_ScaleXYZ,k*FFD_ScaleXYZ);
							vec = new PVector(Vx[index], Vy[index], Vz[index]);
							col =vec.mag();
							loads[index] += vec.mag();	//no wrong... coz here,you lose the info, in which direction the vec is pointing.
							//println(Vz[index]);
							fill(abs(col)*1000);
							//println(abs(col)*1000);
							box(FFD_ScaleXYZ-1);
							
							//draw the vertex! to check the direction of the box. if it will move potentially
							//consider mass?
							popMatrix();							
						}

					}
				}
			}
			popStyle();
		}
			*/

        /*
		void drawSection(int xyzSection, int XYZcase){
			strokeWeight(2);			
			int index;
			
			switch (XYZcase){
			case 0:
				for(int k=0; k<G3; k++){
					for (int j=0; j<G2; j++){
						index = IndexCell(xyzSection,j,k);	
						if (FFD_tglVal1 == true){
							beginShape();
							vertex(xyzSection*FFD_ScaleXYZ,j*FFD_ScaleXYZ,k*FFD_ScaleXYZ);
							vertex(xyzSection*FFD_ScaleXYZ+(Vx[index]*FFD_ScaleXYZ*FFD_ScaleArrows), 
									j*FFD_ScaleXYZ+(Vy[index]*FFD_ScaleXYZ*FFD_ScaleArrows), 
									k*FFD_ScaleXYZ+(Vz[index]*FFD_ScaleXYZ*FFD_ScaleArrows));
							endShape();
						}
						if (FFD_tglVal2 == true){
							noFill();
							pushMatrix();
							translate(xyzSection*FFD_ScaleXYZ,j*FFD_ScaleXYZ,k*FFD_ScaleXYZ);
							box(density[index]);
							popMatrix();
						}
						if (FFD_tglVal3 == true){
							if (solid[index] == 1){
								fill(0,80);
								pushMatrix();
								translate(xyzSection*FFD_ScaleXYZ,j*FFD_ScaleXYZ,k*FFD_ScaleXYZ);
								box(FFD_ScaleXYZ-1);
								popMatrix();
							}						
						}
					}	
				}
				break;
			case 1:
				for(int k=0; k<G3; k++){
					for (int i=0; i<G1; i++){
						index = IndexCell(i, xyzSection, k);	
						if (FFD_tglVal1 == true){
							beginShape();
							vertex(i*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ, k*FFD_ScaleXYZ);
							vertex(i*FFD_ScaleXYZ+(Vx[index]*FFD_ScaleXYZ*FFD_ScaleArrows),
									xyzSection*FFD_ScaleXYZ+(Vy[index]*FFD_ScaleXYZ*FFD_ScaleArrows),  
									k*FFD_ScaleXYZ+(Vz[index]*FFD_ScaleXYZ*FFD_ScaleArrows));
							endShape();
						}
						if (FFD_tglVal2 == true){
							noFill();
							pushMatrix();
							translate(i*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ, k*FFD_ScaleXYZ);
							box(density[index]);
							popMatrix();
						}
						if (FFD_tglVal3 == true){
							if (solid[index] == 1){
								fill(0,80);
								pushMatrix();
								translate(i*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ, k*FFD_ScaleXYZ);
								box(FFD_ScaleXYZ-1);
								popMatrix();
							}						
						}
					}	
				}
				break;
			case 2:
				for(int j=0; j<G2; j++){
					for (int i=0; i<G1; i++){
						index = IndexCell(i, j, xyzSection);	
						if (FFD_tglVal1 == true){
							beginShape();
							vertex(i*FFD_ScaleXYZ, j*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ);
							vertex(i*FFD_ScaleXYZ+(Vx[index]*FFD_ScaleXYZ*FFD_ScaleArrows),
									j*FFD_ScaleXYZ+(Vy[index]*FFD_ScaleXYZ*FFD_ScaleArrows),
									xyzSection*FFD_ScaleXYZ+(Vz[index]*FFD_ScaleXYZ*FFD_ScaleArrows));
							endShape();
						}
						if (FFD_tglVal2 == true){
							noFill();
							pushMatrix();
							translate(i*FFD_ScaleXYZ, j*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ);
							box(density[index]);
							popMatrix();
						}
						if (FFD_tglVal3 == true){
							if (solid[index] == 1){
								fill(0,80);
								pushMatrix();
								translate(i*FFD_ScaleXYZ, j*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ);
								box(FFD_ScaleXYZ-1);
								popMatrix();
							}						
						}
					}	
				}
				break;
			}
			
		}
        */
        /*
		void drawSectionColors(int xyzSection, int XYZcase){
			pushStyle();
			strokeWeight(1);			
			int index;
			colorMode(RGB,255);
			float RR, GG, BB;
			float Top, Low, Third;
			Top = 0.3f;
			Low = 0f;
			Third = (Top-Low)/5f;
			PVector vec;
			
			switch (XYZcase){
			case 0:
				for(int k=0; k<G3; k++){
					for (int j=0; j<G2; j++){
						index = IndexCell(xyzSection,j,k);	
						if (FFD_tglVal1 == true){
							vec = new PVector(Vx[index], Vy[index], Vz[index]);
							if (vec.mag() > Third && vec.mag() <= (2f*Third)){
								RR = (vec.mag()-Third) * (255f/Third);
								GG = 0f;
								BB = 255-((vec.mag()-Third) * (255f/Third));
							}else if(vec.mag() > (2f*Third)){
								RR = 255f;
								GG = ((vec.mag()-(2f*Third)) * (255f/Third));
								BB = 0f;
							}else {
								RR = 0f;
								GG = 0f;
								BB = 255f;
							}
							stroke(RR, GG, BB);
							pushMatrix();
							translate(xyzSection*FFD_ScaleXYZ,j*FFD_ScaleXYZ,k*FFD_ScaleXYZ);
							beginShape(LINES);
							vertex(0,0,0);
							vec.normalize();
							vec.scaleTo(FFD_ScaleXYZ);
							vertex(vec.x, vec.y, vec.z);
							
							vertex(vec.x, vec.y, vec.z);
							vec.div(1.5f);
							vertex(vec.x, vec.y, vec.z+FFD_ScaleXYZ/5);
							
							vec.normalize();
							vec.scaleTo(FFD_ScaleXYZ);
							vertex(vec.x, vec.y, vec.z);
							vec.div(1.5f);
							vertex(vec.x, vec.y, vec.z-FFD_ScaleXYZ/5);
							endShape();
							popMatrix();	
						}
						if (FFD_tglVal2 == true){
							noFill();
							pushMatrix();
							translate(xyzSection*FFD_ScaleXYZ,j*FFD_ScaleXYZ,k*FFD_ScaleXYZ);
							box(density[index]);
							popMatrix();
						}
						if (FFD_tglVal3 == true){
							if (solid[index] == 1){
								fill(0,80);
								pushMatrix();
								translate(xyzSection*FFD_ScaleXYZ,j*FFD_ScaleXYZ,k*FFD_ScaleXYZ);
								box(FFD_ScaleXYZ-1);
								popMatrix();
							}						
						}
					}	
				}
				break;
			case 1:
				for(int k=0; k<G3; k++){
					for (int i=0; i<G1; i++){
						index = IndexCell(i, xyzSection, k);	
						if (FFD_tglVal1 == true){
							vec = new PVector(Vx[index], Vy[index], Vz[index]);
							if (vec.mag() > Third && vec.mag() <= (2f*Third)){
								RR = (vec.mag()-Third) * (255f/Third);
								GG = 0f;
								BB = 255-((vec.mag()-Third) * (255f/Third));
							}else if(vec.mag() > (2f*Third)){
								RR = 255f;
								GG = ((vec.mag()-(2f*Third)) * (255f/Third));
								BB = 0f;
							}else {
								RR = 0f;
								GG = 0f;
								BB = 255f;
							}
							stroke(RR, GG, BB);
							pushMatrix();
							translate(i*FFD_ScaleXYZ,xyzSection*FFD_ScaleXYZ,k*FFD_ScaleXYZ);
							beginShape(LINES);
							vertex(0,0,0);
							vec.normalize();
							vec.scaleTo(FFD_ScaleXYZ);
							vertex(vec.x, vec.y, vec.z);
							
							vertex(vec.x, vec.y, vec.z);
							vec.div(1.5f);
							vertex(vec.x+FFD_ScaleXYZ/5, vec.y, vec.z);
							
							vec.normalize();
							vec.scaleTo(FFD_ScaleXYZ);
							vertex(vec.x, vec.y, vec.z);
							vec.div(1.5f);
							vertex(vec.x-FFD_ScaleXYZ/5, vec.y, vec.z);
							endShape();
							popMatrix();
						}
						if (FFD_tglVal2 == true){
							noFill();
							pushMatrix();
							translate(i*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ, k*FFD_ScaleXYZ);
							box(density[index]);
							popMatrix();
						}
						if (FFD_tglVal3 == true){
							if (solid[index] == 1){
								fill(0,80);
								pushMatrix();
								translate(i*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ, k*FFD_ScaleXYZ);
								box(FFD_ScaleXYZ-1);
								popMatrix();
							}						
						}
					}	
				}
				break;
			case 2:
				for(int j=0; j<G2; j++){
					for (int i=0; i<G1; i++){
						index = IndexCell(i, j, xyzSection);	
						if (FFD_tglVal1 == true){
							vec = new PVector(Vx[index], Vy[index], Vz[index]);
							if (vec.mag() > Third && vec.mag() <= (2f*Third)){
								RR = (vec.mag()-Third) * (255f/Third);
								GG = 0f;
								BB = 255-((vec.mag()-Third) * (255f/Third));
							}else if(vec.mag() > (2f*Third)){
								RR = 255f;
								GG = ((vec.mag()-(2f*Third)) * (255f/Third));
								BB = 0f;
							}else {
								RR = 0f;
								GG = 0f;
								BB = 255f;
							}
							stroke(RR, GG, BB);
							pushMatrix();
							translate(i*FFD_ScaleXYZ,j*FFD_ScaleXYZ,xyzSection*FFD_ScaleXYZ);
							beginShape(LINES);
							vertex(0,0,0);
							vec.normalize();
							vec.scaleTo(FFD_ScaleXYZ);
							vertex(vec.x, vec.y, vec.z);
							
							vertex(vec.x, vec.y, vec.z);
							vec.div(1.5f);
							vertex(vec.x, vec.y+FFD_ScaleXYZ/5, vec.z);
							
							vec.normalize();
							vec.scaleTo(FFD_ScaleXYZ);
							vertex(vec.x, vec.y, vec.z);
							vec.div(1.5f);
							vertex(vec.x, vec.y-FFD_ScaleXYZ/5, vec.z);
							endShape();
							popMatrix();
						}
						if (FFD_tglVal2 == true){
							noFill();
							pushMatrix();
							translate(i*FFD_ScaleXYZ, j*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ);
							box(density[index]);
							popMatrix();
						}
						if (FFD_tglVal3 == true){
							if (solid[index] == 1){
								fill(0,80);
								pushMatrix();
								translate(i*FFD_ScaleXYZ, j*FFD_ScaleXYZ, xyzSection*FFD_ScaleXYZ);
								box(FFD_ScaleXYZ-1);
								popMatrix();
							}						
						}
					}	
				}
				break;
			}
			popStyle();
		}
	*/



        public void FluidRun()
        {

            addInternalBoundary();
            addWind();
            OpDiffuse(1, Vx0, Vx, viscosity, dt, iter);
            OpDiffuse(2, Vy0, Vy, viscosity, dt, iter);
            OpDiffuse(3, Vz0, Vz, viscosity, dt, iter);

            OpProject(Vx0, Vy0, Vz0, Vx, Vy, iter);

            OpAdvect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt);
            OpAdvect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt);
            OpAdvect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt);

            OpProject(Vx, Vy, Vz, Vx0, Vy0, iter);

            OpDiffuse(0, densityScratch, density, diffusion, dt, iter);
            OpAdvect(0, density, densityScratch, Vx, Vy, Vz, dt);

            //OpLoads();
            OpForce();
        }

        void reset()
        {
            for (int i = 0; i < nCells; i++)
            {
                Vx[i] = 0f;
                Vy[i] = 0f;
                Vz[i] = 0f;
                Vx0[i] = 0f;
                Vy0[i] = 0f;
                Vz0[i] = 0f;
                density[i] = 0.0f;
                densityScratch[i] = 0.0f;
                solid[i] = 0;

                loads[i] = 0f;
            }
        }

        void addDye(int x, int y, int z, float amount)
        {
            density[IndexCell(x, y, z)] += amount;
        }

        void addVelocity(int x, int y, int z, float VxAdd, float VyAdd, float VzAdd)
        {
            int index = IndexCell(x, y, z);
            Vx[index] += VxAdd;
            Vy[index] += VyAdd;
            Vz[index] += VzAdd;
        }

        void addWind()
        {
            //			for(int i=((G1/2)-2); i<((G1/2)+2); i++){
            //				for (int k=0; k<G3-1; k++){
            //					
            //					addDye(i,1,k,densityIn);
            //					addVelocity(i,1,k,0.0f, velocityIn, 0.0f);
            //				}
            //			}

            //			for(int i=((G1/2)-2); i<((G1/2)+2); i++){
            //				for (int j=((G2/2)-2); j<((G2/2)+2); j++){
            //					
            //					addDye(i,j,1,densityIn);
            //					addVelocity(i,j,1, 0.0f, 0.0f,velocityIn);
            //				}
            //			}

            //			for(int i=((G1/2)-2); i<((G1/2)+2); i++){
            //				for (int j=((G2/2)-2); j<((G2/2)+2); j++){
            //					
            //					addDye(i,j,1,densityIn);
            //					addVelocity(i,j,1, 0.0f, 0.0f,velocityIn);
            //				}
            //			}
            for (int i = 0; i < G1; i++)
            {
                for (int j = 0; j < G2; j++)
                {

                    addDye(i, j, 1, densityIn);
                    addVelocity(i, j, 1, 0.0f, 0.0f, velocityIn);
                }
            }

        }

        void addInternalBoundary()
        {

            for (int i = 0; i < G1; i++)
            {
                for (int j = 0; j < G2; j++)
                {
                    solid[IndexCell(i, j, G3 / 2)] = 0;
                }
            }

            for (int i = ((G1 / 2) - (FFD_btnVal7 / 2)); i < ((G1 / 2) + (FFD_btnVal7 / 2)); i++)
            {
                for (int j = ((G2 / 2) - (FFD_btnVal7 / 2)); j < ((G2 / 2) + (FFD_btnVal7 / 2)); j++)
                {
                    solid[IndexCell(i, j, G3 / 2)] = 1;
                }
            }


        }

        void OpForce()
        {


            PVector normal;
            PVector vectorAdd;
            float[] VxBuff = new float[Vx.Length];
            float[] VyBuff = new float[Vy.Length];
            float[] VzBuff = new float[Vz.Length];
            int a = 1;

            for (int m = 0; m < solid.Length; m++)
            {
                if (solid[m] == 1)
                {
                    VxBuff[m] = 0;
                    VyBuff[m] = 0;
                    VzBuff[m] = 0;
                    if (solidx[m] > 1)
                    {
                        normal = new PVector(-1, 0, 0);
                        vectorAdd = weightedVector(solidx[m] - a, solidy[m], solidz[m], normal);
                        VxBuff[m] += vectorAdd.x;
                        VyBuff[m] += vectorAdd.y;
                        VzBuff[m] += vectorAdd.z;
                    }
                    if (solidx[m] < G1 - 1)
                    {
                        normal = new PVector(1, 0, 0);
                        vectorAdd = weightedVector(solidx[m] + a, solidy[m], solidz[m], normal);
                        VxBuff[m] += vectorAdd.x;
                        VyBuff[m] += vectorAdd.y;
                        VzBuff[m] += vectorAdd.z;
                    }
                    if (solidy[m] > 1)
                    {
                        normal = new PVector(0, -1, 0);
                        vectorAdd = weightedVector(solidx[m], solidy[m] - a, solidz[m], normal);
                        VxBuff[m] += vectorAdd.x;
                        VyBuff[m] += vectorAdd.y;
                        VzBuff[m] += vectorAdd.z;
                    }
                    if (solidy[m] < G2 - 1)
                    {
                        normal = new PVector(0, 1, 0);
                        vectorAdd = weightedVector(solidx[m], solidy[m] + a, solidz[m], normal);
                        VxBuff[m] += vectorAdd.x;
                        VyBuff[m] += vectorAdd.y;
                        VzBuff[m] += vectorAdd.z;
                    }
                    if (solidz[m] > 1)
                    {
                        normal = new PVector(0, 0, -1);
                        vectorAdd = weightedVector(solidx[m], solidy[m], solidz[m] - a, normal);
                        VxBuff[m] += vectorAdd.x;
                        VyBuff[m] += vectorAdd.y;
                        VzBuff[m] += vectorAdd.z;
                    }
                    if (solidy[m] < G3 - 1)
                    {
                        normal = new PVector(0, 0, 1);
                        vectorAdd = weightedVector(solidx[m], solidy[m], solidz[m] + a, normal);
                        VxBuff[m] += vectorAdd.x;
                        VyBuff[m] += vectorAdd.y;
                        VzBuff[m] += vectorAdd.z;
                    }
                }
            }

            for (int m = 0; m < solid.Length; m++)
            {
                if (solid[m] == 1)
                {
                    Vx[m] = VxBuff[m];
                    Vy[m] = VyBuff[m];
                    Vz[m] = VzBuff[m];
                }
            }



        }

        PVector weightedVector(int x, int y, int z, PVector normal)
        {
            int ix;
            float angle;
            float vectorWeight;
            PVector fluid, vectorAdd;

            ix = IndexCell(x, y, z);

            fluid = new PVector(Vx[ix], Vy[ix], Vz[ix]);
            angle = PVector.angleBetween(normal, fluid);
            angle = PVector.RadianToDegree(angle);

            vectorWeight = angle <= 90f ? 1f - ((1f / 90f) * angle) : (1f / 90f) * (angle - 90f);

            //println(angle + " " + fluid.mag() + " " + vectorWeight);
            vectorAdd = PVector.multiply(fluid, vectorWeight);
            if (PVector.magnitude(vectorAdd) > 0)
            {
                //println("hi");
            }
            else
            {
                vectorAdd = new PVector(0, 0, 0);
            }
            return vectorAdd;
        }




        void OpDiffuse(int b, float[] cell, float[] cell0, float diff, float dt, int iter)
        {
            float a;
            if (b == 1) a = dt * diff * (G2 - 2) * (G3 - 2);
            else if (b == 2) a = dt * diff * (G1 - 2) * (G3 - 2);
            else if (b == 3) a = dt * diff * (G1 - 2) * (G2 - 2);
            else a = dt * diff * (G1 - 2) * (G2 - 2);
            LinearSolver(b, cell, cell0, a, 1 + 6 * a, iter);
        }

        void OpProject(float[] velocX, float[] velocY, float[] velocZ, float[] p, float[] div, int iter)
        {

            int[] ic = new int[6];

            for (int k = 1; k < G3 - 1; k++)
            {
                for (int j = 1; j < G2 - 1; j++)
                {
                    for (int i = 1; i < G1 - 1; i++)
                    {
                        if (solid[IndexCell(i, j, k)] == 0)
                        {
                            ic[0] = IndexCell(i + 1, j, k);
                            ic[1] = IndexCell(i - 1, j, k);
                            ic[2] = IndexCell(i, j + 1, k);
                            ic[3] = IndexCell(i, j - 1, k);
                            ic[4] = IndexCell(i, j, k + 1);
                            ic[5] = IndexCell(i, j, k - 1);
                            for (int cc = 0; cc < 5; cc++)
                            {
                                if (solid[ic[cc]] == 1) ic[cc] = IndexCell(i, j, k);
                            }
                            div[IndexCell(i, j, k)] = -0.5f * (
                                                velocX[ic[0]] -
                                                velocX[ic[1]] +
                                                velocY[ic[2]] -
                                                velocY[ic[3]] +
                                                velocZ[ic[4]] -
                                                velocZ[ic[5]]
                                                ) / G1;
                            p[IndexCell(i, j, k)] = 0;
                        }
                    }
                }
            }
            Boundary(0, div);
            Boundary(0, p);
            LinearSolver(0, p, div, 1, 6, iter);

            for (int k = 1; k < G3 - 1; k++)
            {
                for (int j = 1; j < G2 - 1; j++)
                {
                    for (int i = 1; i < G1 - 1; i++)
                    {
                        if (solid[IndexCell(i, j, k)] == 0)
                        {
                            ic[0] = IndexCell(i + 1, j, k);
                            ic[1] = IndexCell(i - 1, j, k);
                            ic[2] = IndexCell(i, j + 1, k);
                            ic[3] = IndexCell(i, j - 1, k);
                            ic[4] = IndexCell(i, j, k + 1);
                            ic[5] = IndexCell(i, j, k - 1);
                            for (int cc = 0; cc < 5; cc++)
                            {
                                if (solid[ic[cc]] == 1) ic[cc] = IndexCell(i, j, k);
                            }
                            velocX[IndexCell(i, j, k)] -= 0.5f * (p[ic[0]] -
                                                                p[ic[1]]) * G1;
                            velocY[IndexCell(i, j, k)] -= 0.5f * (p[ic[2]] -
                                                                p[ic[3]]) * G2;
                            velocZ[IndexCell(i, j, k)] -= 0.5f * (p[ic[4]] -
                                                                p[ic[5]]) * G3;
                        }
                    }
                }
            }
            Boundary(1, velocX);
            Boundary(2, velocY);
            Boundary(3, velocZ);

        }

        void OpAdvect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float[] velocZ, float dt)
        {
            float i0, i1, j0, j1, k0, k1;
            float dtx = dt * (G1 - 2);
            float dty = dt * (G2 - 2);
            float dtz = dt * (G3 - 2);
            float s0, s1, t0, t1, u0, u1;
            float tmp1, tmp2, tmp3, x, y, z;
            float G1float = G1;
            float G2float = G2;
            float G3float = G3;
            float ifloat, jfloat, kfloat;
            int i, j, k;

            for (k = 1, kfloat = 1f; k < G3 - 1; k++, kfloat++)
            {
                for (j = 1, jfloat = 1f; j < G2 - 1; j++, jfloat++)
                {
                    for (i = 1, ifloat = 1f; i < G1 - 1; i++, ifloat++)
                    {
                        if (solid[IndexCell(i, j, k)] == 0)
                        {
                            tmp1 = dtx * velocX[IndexCell(i, j, k)];
                            tmp2 = dty * velocY[IndexCell(i, j, k)];
                            tmp3 = dtz * velocZ[IndexCell(i, j, k)];
                            x = ifloat - tmp1;
                            y = jfloat - tmp2;
                            z = kfloat - tmp3;

                            if (x < 0.5f) x = 0.5f;
                            if (x > G1float - 1.5f) x = G1float - 1.5f;
                            i0 = (float)Math.Floor(x);
                            i1 = i0 + 1.0f;
                            if (y < 0.5f) y = 0.5f;
                            if (y > G2float - 1.5f) y = G2float - 1.5f;
                            j0 = (float)Math.Floor(y);
                            j1 = j0 + 1.0f;
                            if (z < 0.5f) z = 0.5f;
                            if (z > G3float - 1.5f) z = G3float - 1.5f;
                            k0 = (float)Math.Floor(z);
                            k1 = k0 + 1.0f;

                            s1 = x - i0;
                            s0 = 1.0f - s1;
                            t1 = y - j0;
                            t0 = 1.0f - t1;
                            u1 = z - k0;
                            u0 = 1.0f - u1;

                            int i0i = (int)i0;
                            int i1i = (int)i1;
                            int j0i = (int)j0;
                            int j1i = (int)j1;
                            int k0i = (int)k0;
                            int k1i = (int)k1;

                            d[IndexCell(i, j, k)] =
                                    s0 * (t0 * (u0 * d0[IndexCell(i0i, j0i, k0i)]
                                                + u1 * d0[IndexCell(i0i, j0i, k1i)])
                                        + (t1 * (u0 * d0[IndexCell(i0i, j1i, k0i)]
                                                + u1 * d0[IndexCell(i0i, j1i, k1i)])))
                                + s1 * (t0 * (u0 * d0[IndexCell(i1i, j0i, k0i)]
                                                + u1 * d0[IndexCell(i1i, j0i, k1i)])
                                        + (t1 * (u0 * d0[IndexCell(i1i, j1i, k0i)]
                                                + u1 * d0[IndexCell(i1i, j1i, k1i)])));
                        }
                    }
                }
            }
            Boundary(b, d);
        }

        void Boundary(int b, float[] cell)
        {
            for (int j = 1; j < G2 - 1; j++)
            {
                for (int i = 1; i < G1 - 1; i++)
                {
                    cell[IndexCell(i, j, 0)] = b == 3 ? -cell[IndexCell(i, j, 1)] : cell[IndexCell(i, j, 1)];
                    cell[IndexCell(i, j, G3 - 1)] = b == 3 ? -cell[IndexCell(i, j, G3 - 2)] : cell[IndexCell(i, j, G3 - 2)];
                }
            }

            for (int k = 1; k < G3 - 1; k++)
            {
                for (int i = 1; i < G1 - 1; i++)
                {
                    cell[IndexCell(i, 0, k)] = b == 2 ? -cell[IndexCell(i, 1, k)] : cell[IndexCell(i, 1, k)];
                    cell[IndexCell(i, G2 - 1, k)] = b == 2 ? -cell[IndexCell(i, G2 - 2, k)] : cell[IndexCell(i, G2 - 2, k)];
                }
            }

            for (int k = 1; k < G3 - 1; k++)
            {
                for (int j = 1; j < G2 - 1; j++)
                {
                    cell[IndexCell(0, j, k)] = b == 1 ? -cell[IndexCell(1, j, k)] : cell[IndexCell(1, j, k)];
                    cell[IndexCell(G1 - 1, j, k)] = b == 1 ? -cell[IndexCell(G1 - 2, j, k)] : cell[IndexCell(G1 - 2, j, k)];
                }
            }


            //first bottom row, first corner
            cell[IndexCell(0, 0, 0)] = 0.33f * (cell[IndexCell(1, 0, 0)] +
                                                        cell[IndexCell(0, 1, 0)] +
                                                        cell[IndexCell(0, 0, 1)]);
            //first bottom row, second corner
            cell[IndexCell(0, G2 - 1, 0)] = 0.33f * (cell[IndexCell(1, G2 - 1, 0)] +
                                                        cell[IndexCell(0, G2 - 2, 0)] +
                                                        cell[IndexCell(0, G2 - 1, 1)]);
            //first top row, first corner
            cell[IndexCell(0, 0, G3 - 1)] = 0.33f * (cell[IndexCell(1, 0, G3 - 1)] +
                                                        cell[IndexCell(0, 1, G3 - 1)] +
                                                        cell[IndexCell(0, 0, G3 - 2)]);
            //first top row, second corner
            cell[IndexCell(0, G2 - 1, G3 - 1)] = 0.33f * (cell[IndexCell(1, G2 - 1, G3 - 1)] +
                                                        cell[IndexCell(0, G2 - 2, G3 - 1)] +
                                                        cell[IndexCell(0, G2 - 1, G3 - 2)]);
            //last bottom row, first corner
            cell[IndexCell(G1 - 1, 0, 0)] = 0.33f * (cell[IndexCell(G1 - 2, 0, 0)] +
                                                        cell[IndexCell(G1 - 1, 1, 0)] +
                                                        cell[IndexCell(G1 - 1, 0, 1)]);
            //last bottom row, second corner
            cell[IndexCell(G1 - 1, G2 - 1, 0)] = 0.33f * (cell[IndexCell(G1 - 2, G1 - 1, 0)] +
                                                        cell[IndexCell(G1 - 1, G2 - 2, 0)] +
                                                        cell[IndexCell(G1 - 1, G2 - 1, 1)]);
            //last top row, first corner
            cell[IndexCell(G1 - 1, 0, G3 - 1)] = 0.33f * (cell[IndexCell(G1 - 2, 0, G3 - 1)] +
                                                        cell[IndexCell(G1 - 1, 1, G3 - 1)] +
                                                        cell[IndexCell(G1 - 1, 0, G3 - 2)]);
            //last top row, second corner
            cell[IndexCell(G1 - 1, G2 - 1, G3 - 1)] = 0.33f * (cell[IndexCell(G1 - 2, G2 - 1, G3 - 1)] +
                                                        cell[IndexCell(G1 - 1, G2 - 2, G3 - 1)] +
                                                        cell[IndexCell(G1 - 1, G2 - 1, G3 - 2)]);

            BoundaryInternal(b, cell);
        }

        void BoundaryInternal(int b, float[] cell)
        {
            int index;
            for (int m = 0; m < solid.Length; m++)
            {
                if (solid[m] == 1)
                {
                    switch (b)
                    {
                        case 1:
                            if (solidx[m] > 0)
                            {
                                index = IndexCell(solidx[m] - 1, solidy[m], solidz[m]);
                                if (solid[index] == 0) cell[index] = -cell[index];
                            }
                            if (solidx[m] < G1)
                            {
                                index = IndexCell(solidx[m] + 1, solidy[m], solidz[m]);
                                if (solid[index] == 0) cell[index] = -cell[index];
                            }
                            break;
                        case 2:
                            if (solidy[m] > 0)
                            {
                                index = IndexCell(solidx[m], solidy[m] - 1, solidz[m]);
                                if (solid[index] == 0) cell[index] = -cell[index];
                            }
                            if (solidy[m] < G2)
                            {
                                index = IndexCell(solidx[m], solidy[m] + 1, solidz[m]);
                                if (solid[index] == 0) cell[index] = -cell[index];
                            }
                            break;
                        case 3:
                            if (solidz[m] > 0)
                            {
                                index = IndexCell(solidx[m], solidy[m], solidz[m] - 1);
                                if (solid[index] == 0) cell[index] = -cell[index];
                            }
                            if (solidz[m] < G3)
                            {
                                index = IndexCell(solidx[m], solidy[m], solidz[m] + 1);
                                if (solid[index] == 0) cell[index] = -cell[index];
                            }
                            break;
                    }

                    cell[IndexCell(solidx[m], solidy[m], solidz[m])] = 0f;
                }
            }
        }

        void LinearSolver(int b, float[] cell, float[] cell0, float a, float c, int iter)
        {
            int[] ic = new int[6];

            float cRecip = 1.0f / c;
            for (int m = 0; m < iter; m++)
            {
                for (int k = 1; k < G3 - 1; k++)
                {
                    for (int j = 1; j < G2 - 1; j++)
                    {
                        for (int i = 1; i < G1 - 1; i++)
                        {
                            if (solid[IndexCell(i, j, k)] == 0)
                            {
                                ic[0] = IndexCell(i + 1, j, k);
                                ic[1] = IndexCell(i - 1, j, k);
                                ic[2] = IndexCell(i, j + 1, k);
                                ic[3] = IndexCell(i, j - 1, k);
                                ic[4] = IndexCell(i, j, k + 1);
                                ic[5] = IndexCell(i, j, k - 1);
                                for (int cc = 0; cc < 5; cc++)
                                {
                                    if (solid[ic[cc]] == 1) ic[cc] = IndexCell(i, j, k);
                                }
                                cell[IndexCell(i, j, k)] =
                                        (cell0[IndexCell(i, j, k)] + a *
                                                                    (cell[ic[0]] +
                                                                    cell[ic[1]] +
                                                                    cell[ic[2]] +
                                                                    cell[ic[3]] +
                                                                    cell[ic[4]] +
                                                                    cell[ic[5]]
                                                                    )) * cRecip;
                            }
                        }
                    }
                }
                Boundary(b, cell);
            }
        }

        int IndexCell(int x, int y, int z)
        {
            return (x) + (y * (G1)) + (z * (G2) * (G1));
        }






        /// <summary>
        /// Vector Operations.
        /// </summary>
        private protected class PVector
        {
            public float x, y, z;

            public PVector(float x, float y, float z)
            {
                this.x = x;
                this.y = y;
                this.z = z;
            }


            /// <summary>
            /// Calculates and returns the angle (in radians) between two vectors.
            /// </summary>
            /// <param name="v1">PVector: the x, y, and z components of a PVector</param>
            /// <param name="v2">PVector: the x, y, and z components of a PVector</param>
            /// <returns>float</returns>
            public static float angleBetween(PVector v1, PVector v2)
            {
                float dP = dotProduct(v1, v2);
                float v1mag = magnitude(v1);
                float v2mag = magnitude(v2);

                return (dP / (v1mag * v2mag));
            }

            public static float RadianToDegree(float angle)
            {
                return (float)(angle * (180.0 / Math.PI));
            }

            public static PVector multiply(PVector v, float scalar)
            {
                PVector v_ = new PVector(v.x, v.y, v.z);
                v_.x *= scalar;
                v_.y *= scalar;
                v_.z *= scalar;

                return v_;
            }


            public static float dotProduct(PVector v1, PVector v2)
            {
                return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
            }

            public static float magnitude(PVector v)
            {
                return ((float)Math.Sqrt(Math.Pow(v.x, 2) + Math.Pow(v.y, 2) + Math.Pow(v.z, 2)));
            }


            public static double Norm(double[] x)
            {
                double sum = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    sum += Math.Pow(x[i], 2);
                }

                return Math.Sqrt(sum);

            }


            public static double[] Centroid(double[][] X)
            {
                double[] centroid = new double[X[0].Length];
                for (int i = 0; i < X[0].Length; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < X.Length; j++)
                    {
                        sum += X[j][i];
                    }
                    centroid[i] = sum / X[0].Length;
                }
                return centroid;
            }

        }
    }



}
