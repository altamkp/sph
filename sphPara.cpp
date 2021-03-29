/* Imperial College London
 * AERO-96021 HPC Coursework Assignment 2021
 * Kin Pui Tam - 01338789
 * Date of submission - 24/03/2021
 *
 * SPH member functions for implementing parallelisation
 */

#include "sph.h"

/* @breif	The following input parameters have the same definition throughtout the whole class
 * 			and are therefore documented only once here.
 * @param	rank	The index of the current rank
 * @param	size	The total number of assigned processes
 * */

// Calculate the number of particles assigned to the current rank
int SPH::CalcN(int rank, int size)
{
    if(rmd > 0 && rank < rmd)
	{
	    return N / size + 1;
	}
    else
	{
	    return N / size;
	}
}

// Assign local particle locations for each process
void SPH::AssignLocalX(int rank, int size)
{
    for(int d = 0; d < size; d++)
	{
		// Get number of particles n2 in the rank [d]
	    unsigned int n2 = CalcN(d, size);

	    if(rank == 0)
		{
		    // Create temporary arrray to hold the subarray of xall to send to rank [i]
		    double* xtemp = new double[n2 * 2];

		    // Assign x from xall for rank 0
			if(d == 0)
			{
			    for(unsigned int j = 0; j < n2; j++)
				{
				    for(unsigned int k = 0; k < 2; k++)
					{
					    x[2 * j + k] = xall[2 * j + k];
					}
				}
			}
			// Send the correct part of xall to x from rank 0 to rank [d]
		    else
			{
			    // Calculate the array index of xall to refer to
			    int loc;

			    for(unsigned int j = 0; j < n2; j++)
				{
				    if((rmd == 0) | (d < rmd))
					{
					    loc = 2 * (n2 * d + j);
					}
				    else
					{
					    loc = 2 * (n2 * d + rmd + j);
					}

				    for(unsigned int k = 0; k < 2; k++)
					{
					    xtemp[2 * j + k] = xall[loc + k];
					}
				}
			    MPI_Send(xtemp, n2 * 2, MPI_DOUBLE, d, 0, MPI_COMM_WORLD);
			}
		}

		// For each rank [d] that is not 0, receive part of xallfrom rank 0 
	    if((rank == d) & (d > 0))
		{
		    MPI_Recv(x, n2 * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
}

/* @brief Get array index for div_phi_p
 * @param i 	The the i'th particle out of n in the current rank
 * @param d		The d'th rank that this rank is communicating with
 * @param n2	The number of particles n2 assigned to rank [d]
 * @param j		The j'th particle out of n2 in rank [d]
 * */
int SPH::GetLoc_N(int i, int d, int n2, int j)
{
    if((rmd == 0) | ((rmd != 0) & (d < rmd)))
	{
	    return N * i + n2 * d + j;
	}
    else
	{
	    return N * i + n2 * d + rmd + j;
	}
}

/* @brief	Get rank-independent global array index for reference
 * 			Definition of d, n2 and j is the same as above
 * */
int SPH::GetLoc_n(int d, int n2, int j)
{
    if((rmd == 0) | ((rmd != 0) & (d < rmd)))
	{
	    return n2 * d + j;
	}
    else
	{
	    return n2 * d + rmd + j;
	}
}

// Calculate the sum of density of all particles
double SPH::GetRhoSum(int size)
{
    double sum = 0;

    for(int i = 0; i < size; i++)
	{
	    int ntemp = CalcN(i, size);
	    sum += cblas_dasum(ntemp, rho, 1);
	}

    return sum;
}

// Reset parameters before the start of the simulation
void SPH::ResetParams()
{
    rho = new double[n];
    rho2 = new double[n];
    rhostore = new double[n];
    p = new double[n];
    p2 = new double[n];
	pstore = new double[n];
    a = new double[n * 2];
	div_phi_p = new double[2];
    //div2_phi_v = new double[N * n];
    Fp = new double[n * 2];
    Fv = new double[n * 2];
    Fg = new double[n * 2];
    for(unsigned int i = 0; i < n; i++)
	{
	    Fg[2 * i] = 0;
	}
}

/* @brief		Update the global matrix _all with column size col
 * @param vloc	The local array to send to rank 0
 * @param vglo	The global array to update
 * @param col	The number of columns in that array
 * 				col = 1 for scalars (e.g. rho, p)
 * 				col = 2 for vectors (e.g. x, Fp)
 * */
void SPH::UpdateGlobalArray(int rank, int size, double* vloc, double* vglo, int col)
{
    for(int i = 0; i < size; i++)
	{
	    if((rank == i) & (rank != 0))
		{
		    // Send array to rank 0 for all ranks apart from 0
		    MPI_Send(vloc, n * col, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

		    // Receive globlal array from rank 0
		    MPI_Recv(vglo, N * col, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

	    // Rank 0 receives local array from all other ranks
	    else if(rank == 0)
		{
		    // Reset n for rank 0 such that the correct number of entries are received
		    int n1 = CalcN(i, size);

		    // Store the velocities from rank i in temporary arrar temp
		    double* temp = new double[n1 * col];
		    if(i != 0)
			{
			    MPI_Recv(temp, n1 * col, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		    else
			{
			    for(int j = 0; j < n1; j++)
				{
				    for(int k = 0; k < col; k++)
					{
					    temp[col * j + k] = vloc[col * j + k];
					    // cout << temp[col * j + k] << endl;
					}
				}
			}

		    // Populate the global array vglo
		    for(int j = 0; j < n1; j++)
			{
			    int loc;
			    if(rmd == 0)
				{
				    loc = col * (n1 * i + j);
				    // cout << loc << endl;
				    for(int k = 0; k < col; k++)
					{
					    vglo[loc + k] = temp[col * j + k];
					    // cout << vglo[loc+k] << endl;
					}
				}
			    else
				{
				    if(i < rmd)
					{
					    loc = col * (n1 * i + j);
					    for(int k = 0; k < col; k++)
						{
						    vglo[loc + k] = temp[col * j + k];
						}
					}
				    else
					{
					    loc = col * (n1 * i + j + rmd);
					    for(int k = 0; k < col; k++)
						{
						    vglo[loc + k] = temp[col * j + k];
						}
					}
				}
			}

		    // Deallocate memory for xtemp
		    // delete[] temp;
		}
	}

    // Send global velocity matrix to all other ranks
    if(rank == 0)
	{
	    for(int i = 1; i < size; i++)
		{
		    MPI_Send(vglo, N * col, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
}

/* @brief		An updated version of the function UpdateGlobalArray() above,
 * 				only formulates global array at rank 0
 * @param	send	The array to send to rank 0
 * @param 	recv	The array to receive the array sent from other ranks
 * */

void SPH::UpdateRank0(int rank, int size, double* send, double* recv, int col)
{
    for(int i = 0; i < size; i++)
	{
	    if((rank == i) & (rank != 0))
		{
		    // Send array to rank 0 for all ranks apart from 0
		    MPI_Send(send, n * col, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}

	    // Rank 0 receives local array from all other ranks
	    else if(rank == 0)
		{
		    // Reset n for rank 0 such that the correct number of entries are received
		    int n1 = CalcN(i, size);

		    // Store the velocities from rank i in temporary arrar temp
		    double* temp = new double[n1 * col];
		    if(i != 0)
			{
			    MPI_Recv(temp, n1 * col, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		    else
			{
			    for(int j = 0; j < n1; j++)
				{
				    for(int k = 0; k < col; k++)
					{
					    temp[col * j + k] = send[col * j + k];
					}
				}
			}

		    // Populate the global array vglo
		    for(int j = 0; j < n1; j++)
			{
			    int loc;
			    if(rmd == 0)
				{
				    loc = col * (n1 * i + j);
				    // cout << loc << endl;
				    for(int k = 0; k < col; k++)
					{
					    recv[loc + k] = temp[col * j + k];
					    // cout << vglo[loc+k] << endl;
					}
				}
			    else
				{
				    if(i < rmd)
					{
					    loc = col * (n1 * i + j);
					    for(int k = 0; k < col; k++)
						{
						    recv[loc + k] = temp[col * j + k];
						}
					}
				    else
					{
					    loc = col * (n1 * i + j + rmd);
					    for(int k = 0; k < col; k++)
						{
						    recv[loc + k] = temp[col * j + k];
						}
					}
				}
			}

		    // Deallocate memory for xtemp
		    // delete[] temp;
		}
	}
}
