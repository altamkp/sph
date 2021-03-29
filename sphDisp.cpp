/* Imperial College London
 * AERO-96021 HPC Coursework Assignment 2021
 * Kin Pui Tam - 01338789
 * Date of submission - 24/03/2021
 *
 * SPH member functions for displaying matrices
 */

#include "sph.h"

// Display global positional matrix
void SPH::DispXall(double* xall)
{
    cout << "Position of particles x =" << endl;
    cout << setw(16) << right << "x" << setw(13) << right << "y" << endl;
    for(unsigned int i = 0; i < N; i++)
	{
	    cout << setw(6) << i << ": [";
	    for(unsigned int j = 0; j < 2; j++)
		{
		    cout << setw(12) << fixed << setprecision(4) << xall[2 * i + j];
		}
	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display global velocity matrix
void SPH::DispVall(double* vall)
{
    cout << "Velocity of particles v =" << endl;
    cout << setw(15) << right << "x" << setw(12) << right << "y" << endl;
    for(unsigned int i = 0; i < N; i++)
	{
	    cout << setw(6) << i << ": [";
	    for(unsigned int j = 0; j < 2; j++)
		{
		    cout << setw(12) << fixed << setprecision(5) << vall[2 * i + j];
		}
	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display r_ij = x_i - x_j
void SPH::DispR(double* r, int i, int j)
{
    cout << "r" << i << j << " = [ ";
    for(unsigned int k = 0; k < 2; k++)
	{
	    cout << setw(12) << r[k];
	}
    cout << " ]\n" << endl;
}

// Display q_ij = ||r_ij|| / h
void SPH::DispQ(double q, int i, int j)
{
    cout << "q" << i << j << " = " << setw(4) << q << "\n" << endl;
}

// Display v_ij = v_i - v_j
void SPH::DispVij(double* vij, int i, int j)
{
    cout << "v" << i << j << " = [ ";
    for(unsigned int k = 0; k < 2; k++)
	{
	    cout << setw(8) << vij[k];
	}
    cout << " ]\n" << endl;
}

// Display densities
void SPH::DispRho(double* rho, int rank)
{
    cout << right << "Rank: " << rank << " - rho =" << endl;
    for(unsigned int i = 0; i < n; i++)
	{
	    cout << setw(6) << i << ": [";
	    cout << setw(12) << setprecision(0) << fixed << rho[i];
	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display global density
void SPH::DispRhoAll(double* rhoall)
{
    cout << right << "Global rho = " << endl;
    for(unsigned int i = 0; i < N; i++)
	{
	    cout << setw(6) << i << ": [";
	    cout << setw(12) << setprecision(0) << fixed << rhoall[i];
	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display pressures
void SPH::DispP(double* p, int rank)
{
    cout << right << "Rank: " << rank << " - p =" << endl;
    for(unsigned int i = 0; i < n; i++)
	{
	    cout << setw(6) << right << i << ": [";
	    cout << setw(12) << setprecision(0) << fixed << p[i];
	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display global pressure
void SPH::DispPAll(double* pall)
{
    cout << right << "Global p = " << endl;
    for(unsigned int i = 0; i < N; i++)
	{
	    cout << setw(6) << i << ": [";
	    cout << setw(12) << setprecision(0) << fixed << pall[i];
	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display pressure forces
void SPH::DispFp(double* Fp, int rank)
{
    cout << "Rank: " << rank << " - Fp =" << endl;
    cout << setw(15) << right << "x" << setw(12) << right << "y" << endl;
    for(unsigned int i = 0; i < n; i++)
	{
	    cout << setw(6) << i << ": [";
	    for(unsigned int j = 0; j < 2; j++)
		{
		    cout << setw(12) << setprecision(0) << fixed << Fp[2 * i + j];
		}

	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display viscous forces
void SPH::DispFv(double* Fv, int rank)
{
    cout << left << "Rank: " << rank << " - Fv =" << endl;
    cout << setw(15) << right << "x" << setw(12) << right << "y" << endl;
    for(unsigned int i = 0; i < n; i++)
	{
	    cout << setw(6) << i << ": [";
	    for(unsigned int j = 0; j < 2; j++)
		{
		    cout << setw(12) << setprecision(0) << fixed << Fv[2 * i + j];
		}

	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display gravity forces
void SPH::DispFg(double* Fg, int rank)
{
    cout << "Fg =" << endl;
    cout << setw(15) << right << "x" << setw(12) << right << "y" << endl;
    for(unsigned int i = 0; i < n; i++)
	{
	    cout << setw(6) << i << ": [";
	    for(unsigned int j = 0; j < 2; j++)
		{
		    cout << setw(12) << setprecision(0) << fixed << Fg[2 * i + j];
		}

	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display accelerations
void SPH::DispA(double* a, int rank)
{
    cout << right << "Rank: " << rank << " - a =" << endl;
    cout << setw(15) << right << "x" << setw(6) << right << "y" << endl;
    for(unsigned int i = 0; i < n; i++)
	{
	    cout << setw(6) << i << ": [";
	    for(unsigned int j = 0; j < 2; j++)
		{
		    cout << setw(6) << setprecision(2) << fixed << a[2 * i + j];
		}

	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display velocities
void SPH::DispV(double* v, int rank)
{
    cout << right << "Rank: " << rank << " - v =" << endl;
    cout << setw(15) << right << "x" << setw(12) << right << "y" << endl;
    for(unsigned int i = 0; i < n; i++)
	{
	    cout << setw(6) << i << ": [";
	    for(unsigned int j = 0; j < 2; j++)
		{
		    cout << setw(12) << setprecision(5) << fixed << v[2 * i + j];
		}

	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display positions
void SPH::DispX(double* x, int rank)
{
    cout << right << "Rank: " << rank << " - x =" << endl;
    cout << setw(19) << right << "x" << setw(10) << right << "y" << endl;
    for(unsigned int i = 0; i < n; i++)
	{
	    cout << setw(6) << i << ": [";
	    for(unsigned int j = 0; j < 2; j++)
		{
		    cout << setw(10) << setprecision(4) << fixed << x[2 * i + j];
		}

	    cout << " ]" << endl;
	}
    cout << endl;
}

// Display div_phi_p
void SPH::DispDivPhiP(double* div_phi_p, int rank)
{
    cout << right << "Rank: " << rank << " - div_phi_p =" << endl;
    cout << setw(19) << right << "x" << setw(10) << right << "y" << endl;
    for(unsigned int i = 0; i < N * n; i++)
	{
	    cout << setw(6) << i << ": [";
	    for(unsigned int j = 0; j < 2; j++)
		{
		    cout << setw(10) << setprecision(4) << fixed << div_phi_p[2 * i + j];
		}

	    cout << " ]" << endl;
	}
    cout << endl;
}

// Show progress bar
void SPH::ShowProgress(double t)
{
    string progress = "";

    for(int i = 0; i < int(t / T * 50); i++)
	{
	    progress += ">";
	}
	progress += " " + to_string(int(t / T * 100)) + "%";

    cout << progress << endl;
}