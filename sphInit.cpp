/* Imperial College London
 * AERO-96021 HPC Coursework Assignment 2021
 * Kin Pui Tam - 01338789
 * Date of submission - 24/03/2021
 *
 * SPH class initial particle placement
 */

#include "sph.h"
#include <cmath>

// Initialise matrix x and uniformly place out particles
void SPH::InitialiseX(int rank, int size)
{
    /* Allocate local array size, where n is the number of particles
     * @param n		Number of local particles
     * @param N 	Numer of global particles
     * @param rmd	Remainder of N / size
     * @param maxN	Maximum number local particles n a rank could be assigned
     * */
    rmd = N % size;
    n = CalcN(rank, size);
    if(rmd > 0)
	{
	    maxN = N / size + 1;
	}
    else
	{
	    maxN = N / size;
	}

    // Initialise x to store particle position
    x = new double[n * 2];

    if(rank == 0)
	{
	    xall = new double[N * 2];

	    switch(N)
		{
		case 1: // One particle validation case - [0.5, 0.5]
		    {
			caseName = "One particle - ";
			xall[0] = 0.5;
			xall[1] = 0.5;
			break;
		    }
		case 2: // Two particles validation case - [0.5, 0.5], [0.5, h]
		    {
			caseName = "Two particles - ";
			xall[0] = 0.5;
			xall[1] = 0.5;
			xall[2] = 0.5;
			xall[3] = h;
			break;
		    }
		case 3: // Three particles validation case - [0.5, 0.5], [0.495, h], [0.505, h]
		    {
			caseName = "Three particles - ";
			xall[0] = 0.5;
			xall[1] = 0.5;
			xall[2] = 0.495;
			xall[3] = h;
			xall[4] = 0.505;
			xall[5] = h;
			break;
		    }
		case 4: // Four particles validation case - [0.505, 0.5], [0.515, 0.5], [0.51, 0.45], [0.5, 0.45]
		    {
			caseName = "Four particles - ";
			xall[0] = 0.505;
			xall[1] = 0.5;
			xall[2] = 0.515;
			xall[3] = 0.5;
			xall[4] = 0.51;
			xall[5] = 0.45;
			xall[6] = 0.5;
			xall[7] = 0.45;
			break;
		    }
		case 5: // Extra validation case - N = 5
		    {
			caseName = "Extra test case - ";
			for(unsigned int i = 0; i < N; i++)
			    {
				xall[2 * i] = i * 0.005 + 0.25;
				xall[2 * i + 1] = i * 0.01 + 0.5;
			    }
			break;
		    }
		case 9: // Extra validation case - N = 9
		    {
			caseName = "Extra test case - ";
			for(unsigned int i = 0; i < N; i++)
			    {
				xall[2 * i] = i * 0.01 + 0.25;
				xall[2 * i + 1] = h + i * 0.01;
			    }
			break;
		    }
		case 20: // Extra validation case - N = 20
		    {
			caseName = "Extra test case - ";
			for(unsigned int i = 0; i < N; i++)
			    {
				xall[2 * i] = i * 0.01 + 0.25;
				xall[2 * i + 1] = h + i * 0.01;
			    }
			break;
		    }
		case 99: // Extra validation case - N = 99
		    {
			caseName = "Extra test case - ";
			for(unsigned int i = 0; i < N; i++)
			    {
				xall[2 * i] = h + h * i;
				xall[2 * i + 1] = 0.5;
			    }
			break;
		    }
		case 121: // Extra validation case - N = 121
		    {
			caseName = "Extra test case - ";
			for(unsigned int i = 0; i < sqrt(N); i++)
			    {
				for(unsigned int j = 0; j < sqrt(N); j++)
				    {
					xall[int(sqrt(N)) * 2 * i + 2 * j] = i * h + h;
					xall[int(sqrt(N)) * 2 * i + 2 * j + 1] = j * h + h;
				    }
			    }
			break;
		    }
		case 198: // Extra validation case - N = 198
		    {
			caseName = "Extra test case - ";
			for(unsigned int i = 0; i < N / 2; i++)
			    {
				xall[2 * i] = h + h * i;
				xall[2 * i + 1] = 0.5;
				xall[2 * 99 + 2 * i] = h + h * i;
				xall[2 * 99 + 2 * i + 1] = 0.5 + 0.1;
			    }
			break;
		    }
		case 400: // Dam break test case - grid of particles occupying the region [0, 0.2]^2
		    {
			caseName = "Dam break - ";
			for(unsigned int i = 0; i < sqrt(N); i++)
			    {
				for(unsigned int j = 0; j < sqrt(N); j++)
				    {
					xall[int(sqrt(N)) * 2 * i + 2 * j] = i * h + h;
					xall[int(sqrt(N)) * 2 * i + 2 * j + 1] = j * h + h;
				    }
			    }
			break;
		    }
		case 21 * 31: // Block drop test case - grid of particles occupying the region [0.1, 0.3] x [0.3, 0.6]
		    {
			caseName = "Block drop - ";
			for(unsigned int i = 0; i < 21; i++)
			    {
				for(unsigned int j = 0; j < 31; j++)
				    {
					xall[31 * 2 * i + 2 * j] = i * 0.01 + 0.1;
					xall[31 * 2 * i + 2 * j + 1] = j * 0.01 + 0.3;
				    }
			    }
			break;
		    }
		case 311: // droplet test case - particles occupying a circle of radius 0.1, centred at point [0.5, 0.7]
		    {
			caseName = "Droplet - ";
			ifstream x_N311("N311.txt");
			string line;
			unsigned int i = 0;
			/*
			 * Read N311.txt containing the initial placement of each particle
			 * generated using MATLAB, code given below. This approach is adopted as 
			 * it is most efficient in creating a grid of particles that are evenly distributed
			 * within the circle. 
			 *
			 * 	xpre = zeros(N, 2);
			 * 	xpre(1,:) = [0.5, 0.7];
			 * 	temp = 1;
			 *
			 * 	for i = 1 : 10
			 * 		temp = temp + (2*i-1)^2;
			 *
			 * 		for j = 1 : 2*i+1
			 * 			xtemp = xpre(1,1) + (j-i-1)*0.01;
			 *
			 * 			for k = 1 : 2*i+1
			 * 				ytemp = xpre(1,2) + (k-i-1)*0.01;
			 * 				rSpawn = sqrt((xtemp - xpre(1,1))^2 + (ytemp - xpre(1,2))^2);
			 *
			 * 				if rSpawn < 0.1
			 *   				if (round(abs(xtemp - xpre(1,1)),4) == i*0.01 
			 * 						|| round(abs(ytemp * - xpre(1,2)),4) == i*0.01) 
			 * 						xpre(temp +(2*i+1)*(j-1)+(k-1),1) = xtemp; 
			 * 						xpre(temp +(2*i+1)*(j-1)+(k-1),2) = ytemp;
			 * 					end
			 *				end
			 * 			end
			 *		end
			 *	end
			 *
			 * x(:,1) = nonzeros(xpre(:,1));
			 * x(:,2) = nonzeros(xpre(:,2));
			 */
			if(x_N311.good())
			{
				while(true)
				{
					x_N311 >> line;
					xall[i] = stod(line);
					if(x_N311.eof())
					{
						break;
					}
					i++;
				}
				x_N311.close();
			}
			else
			{
				x311ok = 0;
				cout << "-------------------------------------ERROR--------------------------------------" << endl;
				cout << "File N311.txt is not available.\n" << endl;
			}
			break;
		    }

		default:
		    {
			cout << "-------------------------------------ERROR--------------------------------------" << endl;
			cout << "Particle placement for N = " << N << " is undefined.\n" << endl;
			matchN = 0;
			break;
		    }
		}
	}
    // Broadcast matchN and x311ok from rank 0 to other ranks
    MPI_Bcast(&matchN, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x311ok, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Continue with simulation if the input N matches one of the switch cases above
    if((matchN == 1) & (x311ok == 1))
	{
	    // Assign local position matrix for each process
	    AssignLocalX(rank, size);

	    // Console output
	    if(rank == 0)
		{
		    // Display name of test/validation case
		    cout << caseName;
		    cout << "N = " << N << "\n\n";

		    // Display matrix x (position of particles) in the console
		    DispXall(xall);
		}
	}
}

// Initialise local matrix v and introduce noise for test cases (e.g. dam break)
void SPH::InitialiseV(int rank, int size)
{
    v = new double[n * 2];
    for(unsigned int i = 0; i < n; i++)
	{
	    for(unsigned int k = 0; k < 2; k++)
		{
		    if(N > 4) // For test cases
			{
			    // Set initial velocities to within [-0.005, 0.005]
			    v[2 * i + k] = double(rand()) / RAND_MAX * 0.01 - 0.005;
			}
		    else // For validation cases
			{
			    // Set initiaial velocities to 0
			    v[2 * i + k] = 0;
			}
		}
	}
}