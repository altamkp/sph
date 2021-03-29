/* Imperial College London
 * AERO-96021 HPC Coursework Assignment 2021
 * Kin Pui Tam - 01338789
 * Date of submission - 24/03/2021
 *
 * SPH class implementation
 */

#include "sph.h"
#include <cmath>
#include <cstdlib>
#include <math.h>

#define _USE_MATH_DEFINES

// Custom constructor
SPH::SPH(unsigned int N, double dt, double T, double h, int rank, int size, bool showx, bool showv, bool showe)
    : N(N), dt(dt), T(T), h(h), showx(showx), showv(showv), showe(showe)
{
    // Initialise matrix x and uniformly place out particles
    InitialiseX(rank, size);

    // Continue with simulation if the input N matches one of the switch cases in InitialiseX()
    if((matchN == 1) & (x311ok == 1))
	{
	    // Begin simulation
	    if(rank == 0)
		{
		    cout << "Begin simulating case:" << endl;
		}
	    Physics(rank, size);
	}
}

// Destructor
SPH::~SPH()
{
    // nothing to do
}

// Precalculate densities and pressure at the start
void SPH::PreCalc(int rank, int size)
{
    // Loop [i] - particle [i] out of n in the current rank
    for(unsigned int i = 0; i < maxN; i++)
	{
	    // Reset rho
	    if(i < n)
		{
		    rho[i] = 0;
		}

	    // Loop [d] - loop through ranks
	    for(int d = 0; d < size; d++)
		{
		    // Get size n2 for the loop j that is based on rank d
		    unsigned int n2 = CalcN(d, size);

		    // Reset x2 for the next rank to compare to
		    delete[] x2;
		    x2 = new double[n2 * 2];

		    // Set up x2 for the rank that we are communicating with, x == x2 for that rank
		    for(unsigned int c = 0; c < n2 * 2; c++)
			{
			    if(rank == d)
				{
				    x2[c] = x[c];
				}
			    else
				{
				    x2[c] = 0;
				}
			}

		    // Send x from rank [d] to x2 in the current rank
		    MPI_Bcast(x2, n2 * 2, MPI_DOUBLE, d, MPI_COMM_WORLD);

		    if(i < n)
			{
			    // Loop [j] - particle [j] out of n2 in rank [d]
			    for(unsigned int j = 0; j < n2; j++)
				{
				    // Calculate r_ij = x_i - x_j
				    for(unsigned int k = 0; k < 2; k++)
					{
					    r[k] = x[2 * i + k] - x2[2 * j + k];
					}

				    // Calculate q = ||r_ij|| / h
				    q = cblas_dnrm2(2, r, 1) / h;

				    // Calculate phi_d based on the value of q
				    if(q < 1)
					{
					    phi_d = 4 / (M_PI * pow(h, 2)) * pow((1 - pow(q, 2)), 3);
					}
				    else
					{
					    phi_d = 0;
					}

				    // Update rho
				    rho[i] += m * phi_d;
				}
			}
		}
	}

    for(unsigned int i = 0; i < n; i++)
	{
	    // Set rhostore to store the value of rho since rho will change in each time step
	    rhostore[i] = rho[i];

	    // Calculate pressure
	    p[i] = k * (rho[i] - rho0);
	}
}

// Time step loop for calculating physics
void SPH::Physics(int rank, int size)
{
    // ------------------------ //
    // -----INITIALISATION----- //
    // ------------------------ //

    /* Initiate simulation time. nt and it are used for the condition of the while loop
     * since t < T for some reason returns 1 when t = 2.0000 and T = 2.0000
     * @param t		Simulation time
     * @param nt	Total number of time steps
     * @param it	Index of current time step
     * */
    double t = 0;
    unsigned int nt = T / dt;
    unsigned int it = 0;

    // Add initial noise to test cases
    InitialiseV(rank, size);

    // Update global velocity matrix
    if(rank == 0)
	{
	    vall = new double[N * 2];
	}
    UpdateRank0(rank, size, v, vall, 2);

    /* File output setup - NEED TO CREATE FOLDER OF NAME "OutputN<N>" FIRST FOR FILES TO SHOW
     * E.g. create a folder "OutputN4" for --ic-four-particles.
     * */
    string* filename;
    ofstream* xOut = nullptr;
    ofstream* eOut = nullptr;
    if(rank == 0)
	{
	    filename = new string[N];
	    xOut = new ofstream[N];
	    eOut = new ofstream[1];
	    for(unsigned int i = 0; i < N; i++)
		{
		    filename[i] = to_string(i) + "_N" + to_string(N) + ".txt";
		    xOut[i].open(("./OutputN" + to_string(N) + "/x" + filename[i]).c_str());
		}
	    eOut[0].open(("./OutputN" + to_string(N) + "/energy.txt").c_str());
	}

    // --------------------------- //
    // -----INITIAL CONDITION----- //
    // --------------------------- //

    // Initialise parameters before starting
    ResetParams();

    // Initialise density arrays
    rhostore = new double[n];
    rho = new double[n];
    p = new double[n];

    // Calculate density, assuming m = 1 at the start
    PreCalc(rank, size);

    // Rescale mass
    if(rank == 0)
	{
	    // Calculate the sum of density of all particles
	    double rhoSum = GetRhoSum(size);

	    // Rescale mass
	    m = N * rho0 / rhoSum;

	    // Send the rescaled mass to all other ranks
	    for(int i = 1; i < size; i++)
		{
		    MPI_Send(&m, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
    else
	{
	    // Receive the rescaled mass from rank 0
	    MPI_Recv(&m, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    PreCalc(rank, size);

    /* Begin loop through time until t reaches T set by user.
     * The condition (t < T) was originally used but for some reason it returns true when
     * t = 2.0000 and T = 2.0000, meaning that an extra step would be calculated than required.
     * Using (it < nt) solves the problem since they are integer values
     * */
    while(it < nt)
	{
	    // Announce new time step reached in the console
	    if(rank == 0)
		{
		    cout << "--------------------------------------------------------------------------------" << endl;
		    cout << "Progress:" << endl;
		    cout << "Current time step: " << setw(6) << setprecision(4) << t + dt << " / " << T << endl;
		    ShowProgress(t);
		    cout << left << setw(20) << caseName << "N = " << setw(5) << N << "; h = " << setw(5) << h << endl;
		}

	    // Loop [i] - particle [i] out of n in the current rank
	    for(unsigned int i = 0; i < maxN; i++)
		{
		    // Reset parameters from the previous time step
		    if(i < n)
			{
			    rho[i] = 0;
			    for(unsigned int k = 0; k < 2; k++)
				{
				    Fp[2 * i + k] = 0;
				    Fv[2 * i + k] = 0;
				}
			}

		    // Loop [d] - loop through ranks
		    for(int d = 0; d < size; d++)
			{
			    // Get size n2 for the loop [j] that base on rank [d]
			    unsigned int n2 = CalcN(d, size);

			    // Reset x2 and v2 for the next rank to compare to
			    delete[] x2;
			    delete[] v2;
			    delete[] rho2;
			    delete[] p2;
			    x2 = new double[n2 * 2];
			    v2 = new double[n2 * 2];
			    rho2 = new double[n2];
			    p2 = new double[n2];

			    // ------------------------------------- //
			    // -----CROSS-PROCESS COMMUNICATION----- //
			    // ------------------------------------- //

			    // Set up x2 for the rank that we are communicating with, x == x2 for that rank
			    for(unsigned int c = 0; c < n2; c++)
				{
				    if(rank == d)
					{
					    for(unsigned int k = 0; k < 2; k++)
						{
						    x2[2 * c + k] = x[2 * c + k];
						    v2[2 * c + k] = v[2 * c + k];
						}
					    rho2[c] = rhostore[c];
					    p2[c] = p[c];
					}
				    else
					{
					    for(unsigned int k = 0; k < 2; k++)
						{
						    x2[2 * c + k] = 0;
						    v2[2 * c + k] = 0;
						}
					    rho2[c] = 0;
					    p2[c] = 0;
					}
				}

			    // Send x and v from rank [d] to x2 and v2 in the current rank
			    MPI_Bcast(x2, n2 * 2, MPI_DOUBLE, d, MPI_COMM_WORLD);
			    MPI_Bcast(v2, n2 * 2, MPI_DOUBLE, d, MPI_COMM_WORLD);
			    MPI_Bcast(rho2, n2, MPI_DOUBLE, d, MPI_COMM_WORLD);
			    MPI_Bcast(p2, n2, MPI_DOUBLE, d, MPI_COMM_WORLD);

			    /* Only calculate physics if i < n, ie array index within
			     * the number of particles in this rank
			     * */
			    if(i < n)
				{
				    // Loop [j] - particle [j] out of n2 in rank [d]
				    for(unsigned int j = 0; j < n2; j++)
					{
					    // -------------------------------------------- //
					    // -----CALCULATE CROSS-PROCESS PARAMETERS----- //
					    // -------------------------------------------- //

					    // Calculate r_ij = x_i - x_j and v_ij = v_i - v_j
					    for(unsigned int k = 0; k < 2; k++)
						{
						    r[k] = x[2 * i + k] - x2[2 * j + k];
						    vij[k] = v[2 * i + k] - v2[2 * j + k];
						}

					    // Calculate q = ||r_ij|| / h
					    q = cblas_dnrm2(2, r, 1) / h;

					    // --------------------------- //
					    // -----CALCULATE DENSITY----- //
					    // --------------------------- //

					    // Calculate phi_d based on the value of q
					    if(q < 1)
						{
						    phi_d = 4 / (M_PI * pow(h, 2)) * pow((1 - pow(q, 2)), 3);
						}
					    else
						{
						    phi_d = 0;
						}

					    // Update rho
					    rho[i] += m * phi_d;

					    // check if particle [rank,i] is the same particle as [d,j]
					    bool ijequal = false;
					    if((rank == d) & (i == j))
						{
						    ijequal = true;
						}

					    // ---------------------------------- //
					    // -----CALCULATE PRESSURE FORCE----- //
					    // ---------------------------------- //

					    // Calculate div_phi_p based on the value of q
					    if((q < 1) & !ijequal)
						{
						    double alpha = -30 / (M_PI * pow(h, 3)) * pow((1 - q), 2) / q;
						    cblas_dscal(2, alpha, r, 1);
						    for(unsigned int k = 0; k < 2; k++)
							{
							    div_phi_p[k] = r[k];
							}
						}
					    else
						{
						    for(unsigned int k = 0; k < 2; k++)
							{
							    div_phi_p[k] = 0;
							}
						}

					    // Calculate coefficients
					    double alpha = -m / rho2[j] * (p[i] + p2[j]) / 2;
					    for(unsigned int k = 0; k < 2; k++)
						{
						    Fp_i[k] = 0;
						}

					    // Update pressure force
					    cblas_daxpy(2, alpha, div_phi_p, 1, Fp_i, 1);
					    for(unsigned int k = 0; k < 2; k++)
						{
						    Fp[2 * i + k] += Fp_i[k];
						}

					    // --------------------------------- //
					    // -----CALCULATE VISCOUS FORCE----- //
					    // --------------------------------- //

					    // Calculate div2_phi_v based on the value of q
					    if((q < 1) & !ijequal)
						{
						    div2_phi_v = 40 / (M_PI * pow(h, 4)) * (1 - q);
						}
					    else
						{
						    div2_phi_v = 0;
						}

					    // Update viscous force
					    double theta = -mu * m / rho2[j] * div2_phi_v;
					    cblas_dscal(2, theta, vij, 1);
					    for(unsigned int k = 0; k < 2; k++)
						{
						    Fv[2 * i + k] += vij[k];
						}

					} // End of for loop [j]
				} // End of if-statement (i < n)
			} // End of for loop [d]
		} // End of for loop [i]

	    for(unsigned int i = 0; i < n; i++)
		{
		    // Update rhostore to hold rho for this time step
		    rhostore[i] = rho[i];

		    // ---------------------------- //
		    // -----CALCULATE PRESSURE----- //
		    // ---------------------------- //

		    // Calculate pressure
		    p[i] = k * (rho[i] - rho0);

		    // --------------------------------- //
		    // -----CALCULATE GRAVITY FORCE----- //
		    // --------------------------------- //

		    // Calculate gravity force
		    Fg[2 * i + 1] = -rho[i] * g;
		}

	    // -------------------------- //
	    // -----TIME INTEGRATION----- //
	    // -------------------------- //
	    for(unsigned int i = 0; i < n; i++)
		{
		    for(unsigned int k = 0; k < 2; k++)
			{
			    // Calculate acceleration
			    a[2 * i + k] = (Fp[2 * i + k] + Fv[2 * i + k] + Fg[2 * i + k]) / rho[i];

			    // Calculate velocity
			    v[2 * i + k] += a[2 * i + k] * dt;
			    if(t == 0) // Times velocity by 0.5 for the first time step
				{
				    v[2 * i + k] *= 0.5;
				}

			    /* Boundary check - checks if particle is within distance h from any of the
			     * four boundaries and check if its velocity is point towards that boundary.
			     * Reverse its velocity and multiply e if this is true.
			     * */
			    if(k == 1) // Check for top and bottom boundaries
				{
				    if(((x[2 * i + 1] < h) & (v[2 * i + 1] < 0)) |
				        ((x[2 * i + 1] > 1 - h) & (v[2 * i + 1] > 0)))
					{
					    v[2 * i + 1] *= -e;
					}
				}
			    else // Check for side boundaries
				{
				    if(((x[2 * i] < h) & (v[2 * i] < 0)) | ((x[2 * i] > 1 - h) & (v[2 * i] > 0)))
					{
					    v[2 * i] *= -e;
					}
				}

			    // Calculate position
			    x[2 * i + k] += v[2 * i + k] * dt;
			}
		}

	    // Compile the full position matrix xall of all particles at rank 0 for printing
	    UpdateRank0(rank, size, x, xall, 2);
	    // Update velocities of all particles for printing purpose as requested by user
	    if(showv)
		{
		    UpdateRank0(rank, size, v, vall, 2);
		}

	    // -------------------------- //
	    // -----CALCULATE ENERGY----- //
	    // -------------------------- //

	    /* Initialise parameters for calculating energy
	     * @param vnorm2	Norm of the velocity squared
	     * @param ytotal	Total height of all particles
	     * @param vtemp		Temporary vector to store each particle's velocity
	     * */
	    double vnorm2 = 0;
	    double ytotal = 0;
	    double* vtemp = new double[2];

	    // Calculate vnorm2 and ytotal for the current rank
	    for(unsigned int i = 0; i < n; i++)
		{
		    // Assign particle's velocity to the temporary array
		    for(unsigned int k = 0; k < 2; k++)
			{
			    vtemp[k] = v[2 * i + k];
			}

		    // Calculate norm of the velocity squared
		    vnorm2 += cblas_ddot(2, vtemp, 1, vtemp, 1);

		    // Calculate total height of all particles
		    ytotal += x[2 * i + 1];
		}

	    /* Calculate the contributions of KE and PE from the current rank
	     * @param EkLoc	Local contribution of kinetic energy
	     * @param EpLoc Local contribution of potential energy
	     * */
	    double EkLoc = 0.5 * m * vnorm2;
	    double EpLoc = m * g * ytotal;

	    /* Initialise array to store contributions of KE and PE from all ranks at rank 0.
	     * Initialises at all ranks to prevent uninitialise error in Makefile
	     * @param EkLocs	Array to store contributions of KE from all ranks
	     * @parma EpLocs	Array to store contirbutions of PE from all ranks
	     * */
	    double* EkLocs = nullptr;
	    double* EpLocs = nullptr;
	    if(rank == 0)
		{
		    EkLocs = new double[size];
		    EpLocs = new double[size];
		}

	    // Collate local contributions
	    MPI_Gather(&EkLoc, 1, MPI_DOUBLE, EkLocs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	    MPI_Gather(&EpLoc, 1, MPI_DOUBLE, EpLocs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	    // Computes results at the first process
	    double Ek, Ep;
	    if(rank == 0)
		{
		    Ek = cblas_dasum(size, EkLocs, 1);
		    Ep = cblas_dasum(size, EpLocs, 1);
		    if(showe)
			{
			    cout << "Ek = " << Ek << " ; Ep = " << Ep << endl;
			}
		}

	    // --------------------- //
	    // -----FILE OUTPUT----- //
	    // --------------------- //

	    // File output and console display
	    if(rank == 0)
		{
		    if(N < 5)
			{
			    if(showx)
				{
				    DispXall(xall);
				}
			    if(showv)
				{
				    DispVall(vall);
				}
			}

		    // Output to x_i.txt
		    for(unsigned int i = 0; i < N; i++)
			{
			    if(xOut[i].good())
				{
				    for(unsigned int k = 0; k < 2; k++)
					{
					    xOut[i] << setw(16) << xall[2 * i + k];
					}
				    xOut[i] << endl;
				}
			}

		    // Output to output.txt at every second of simulation time
		    if(fabs(fmod(t, 1) - 1.0) < dt)
			{
			    // Initialise file output for output.txt
			    ofstream xOutAll(
			        ("./OutputN" + to_string(N) + "/output.txt").c_str(), ios::out | ios::trunc);

			    // Output to output.txt if the file can be opened
			    if(xOutAll.good())
				{
				    for(unsigned int i = 0; i < N; i++)
					{
					    for(unsigned int k = 0; k < 2; k++)
						{
						    xOutAll << setw(16) << setprecision(6) << fixed << xall[2 * i + k];
						}
					    xOutAll << endl;
					}
				    xOutAll.close();
				}
			    else
				{
				    xOutAllok = 0;
				}
			}

		    // Output to energy.txt
		    if(eOut[0].good())
			{
			    eOut[0].precision(4);
			    eOut[0] << setw(16) << fixed << t << setw(16) << Ek << setw(16) << Ep << setw(16) << Ep + Ek
			            << endl;
			}
		}

	    // ------------------------ //
	    // -----NEXT TIME STEP----- //
	    // ------------------------ //

	    // Prepare for next time step
	    t += dt;
	    it++;
	    if(rank == 0)
		{
		    cout << endl;
		}
	}

    // ----------------------------------- //
    // -----FILE OUTPUT ERROR MESSAGE----- //
    // ----------------------------------- //
    if(rank == 0)
	{
	    cout << "--------------------------------------------------------------------------------" << endl;
		cout << "Progress:" << endl;
		cout << "Simulation completed." << endl;
		cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 100%" << endl;
	    if(!xOut[0].good() | !eOut[0].good() | (xOutAllok == 0))
		{
		    cout << "Simulation completed but the following files were not created:" << endl;
		    if(!xOut[0].good())
			{
			    cout << "File x<N>.txt output error." << endl;
			}
		    if(!eOut[0].good())
			{
			    cout << "File energy.txt output error." << endl;
			}
		    if(xOutAllok == 0)
			{
			    cout << "File output.txt output error." << endl;
			}
		    cout << endl;
		    cout << "Check if a folder name 'OutputN<N>' is created in the directory where this "
		            "program is located."
		         << endl;
		    cout << "E.g.: folder name for ic-one-particle should be 'OutputN1'\n" << endl;
		}
	    else
		{
		    cout << "All files were successfully written.\n" << endl;
		}
	}

    // ----------------------------- //
    // -----SIMULATION COMPLETE----- //
    // ----------------------------- //

    // Notify the end of the simulation
    if(rank == 0)
	{
	    cout << R"(
   __ _                          _        _
  / ___| ___   _ __ ___   _ __  | |  ___ | |_  ___
 | |    / _ \ | '_ ` _ \ | '_ \ | | / _ \| __|/ _ \
 | |___| (_) || | | | | || |_) || ||  __/| |_|  __/
  \____|\___/ |_| |_| |_|| .__/ |_| \___| \__|\___|
                         |_|
	        )"
	         << endl;
	}

    // De-allocate memory
    if(rank == 0)
	{
	    // delete[] xall;
	    // delete[] vall;
	}
    // delete[] div_phi_p;

    // Close opened files
    if(rank == 0)
	{
	    if(eOut[0].good())
		{
		    eOut[0].close();
		}

	    for(unsigned int i = 0; i < N; i++)
		{
		    if(xOut[i].good())
			{
			    xOut[i].close();
			}
		}
	}
}
