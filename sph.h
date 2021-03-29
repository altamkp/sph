/* Imperial College London
 * AERO-96021 HPC Coursework Assignment 2021
 * Kin Pui Tam - 01338789
 * Date of submission - 24/03/2021
 *
 * SPH class header file
 */

#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include "mpi.h"
#include <cblas.h>
using namespace std;

/*
 * @class 	SPH
 * 
 * @brief 	A class that holds all parameters and member function for simulating
 * 			a 2D smoothed particle hydrodynamic formulation of the Navier-Stokes.
 * 			This class is constructed in cwMain.cpp after the user has defined the
 * 			validation/test case to simulate.
 * 			
 * @file	sph.cpp
 * 			Holds the core of this class for both implementation and simulation.
 * 
 * @file	sphInit.cpp
 * 			Holds functions for initialising particles positions and velocities
 * 			before the simulation.
 * 
 * @file	sphPara.cpp
 * 			Holds functions for implementing parallelisation of the programme.
 * 
 * @file	sphDisp.cpp
 * 			Holds functions for displaying vectors or matrices onto the console
 * 			for prompting and debugging purposes.
 * 
 * @param 	rank	The index of the current rank
 * @param	size	The total number of assigned processes
 */
class SPH
{
public:
	// Initalise public constants
	unsigned int N = 0;			// Number of particles in total
	unsigned int n = 0;			// Number of particles assigned to the current rank
	double dt = 1e-4;			// Time step size
    double T = 2.0000;			// Total integration time
	const double h = 0.01; 		// Radius of influence
	
	// Parameters for parallelisation
	int rmd = 0;				// = N % size, remainder of total number of particles divided by available processes
	unsigned int maxN = 0;		// Maximum number of particles that could be assigned to a rank
	
	// Parameters for checking file I/O availability
	int matchN = 1;				// Checks if the input N by the user matches one of the switch cases in InitiliaseX()
	int x311ok = 1;				// Checks if the file N311.txt can be read
	int xOutAllok = 1;			// Checks if the file output.txt can be written to
	
	// Boolean checks for consle display
	bool showx = false;			// Boolean check for displaying positions in the console
	bool showv = false;			// Boolean check for displaying velocities in the console
	bool showe = false;			// Boolean check for displaying energies in the console
 	
	// Class constructors and destructor
	
	// Default constructor - one particle
    SPH() = default;
	// Custom constructor - based on number of particles
    SPH(unsigned int N, double dt, double T, double h, int rank, int size, bool showx, bool showv, bool showe); 
	// Destructor
    ~SPH();

    // Public functions for simulation
    void Physics(int rank, int size); 				// Time step loop for calculating the physics
	void PreCalc(int rank, int size);				// Pre-calculate densities and pressures at the start
	
	// Public functions for implementing parallelisation
	int CalcN(int rank, int size);					// Calculate the number of particles assigned to current rank
	void AssignLocalX(int rank, int size);			// Assign local particle locations for each process
	int GetLoc_N(int i, int d, int n2, int j);		// Get current rank-dependent global array index for reference
	int GetLoc_n(int d, int n2, int j);				// Get current rank-independent global array index for reference
	double GetRhoSum(int size);						// Get sum of density of all particles
	void ResetParams();								// Reset parameters before simulation
	void UpdateGlobalArray(int rank, int size, double* vloc, double* vglo, int col);	// [obsolete] Update global matrix
	void UpdateRank0(int rank, int size, double* send, double* recv, int col);			// Update global matrix in rank 0 only
	
	// Public functions for displaying vectors/matrices in the console
	void DispXall(double* xall);					// Display global x 
	void DispVall(double* vall);					// Display global v 
    void DispR(double* r, int i, int j);			// Display vector r
    void DispQ(double q, int i, int j);				// Display double q
    void DispVij(double* vij, int i, int j);		// Display vector v_ij
    void DispRho(double* rho, int rank);			// Display vector rho
	void DispRhoAll(double* rhoall);				// [obsolete] Display global rho
    void DispP(double* p, int rank);				// Display vector p
	void DispPAll(double* pall);					// [obsolete] Display global p
    void DispFp(double* Fp, int rank);				// Display matrix Fp
    void DispFv(double* Fv, int rank);				// Display matrix Fv
    void DispFg(double* Fg, int rank);				// Display matrix Fg
    void DispA(double* a, int rank);				// Display matrix a
    void DispV(double* v, int rank);				// Display matrix v
	void DispX(double* x, int rank);				// Display matrix x
	void DispDivPhiP(double* div_phi_p, int rank);	// [obsolete] Display matrix div_phi_p
	void ShowProgress(double t);					// Display progress bar
	
private:
	// Initialise constants
    const double k = 2000.0;    // Gas constant
    const double rho0 = 1000.0; // Resting density
    const double mu = 1.0;      // Viscosity
    const double g = 9.81;      // Acceleration due to gravity
    const double e = 0.5;       // Coefficient of resdtitution

    // Initialise initial conditions (ie particle locations at t = 0)
	string caseName = "";						// Name of validation/test case
    double* x = new double[n * 2];				// Local particle position
	double* x2 = new double[n * 2];				// Particle position of the rank that the current rank is communicating to 
	double* xall;								// Global particle position - initiated in rank 0 only
    void InitialiseX(int rank, int size); 		// Spread particles uniformly in the grid
	void InitialiseV(int rank, int size);		// Introduce noise for test cases

	// ------------------------------------------------- //
    // -----Define matrices for calculating physics----- //
	// ------------------------------------------------- //
	
	/*
     * @brief Initiate parameters for density calculation
	 * @param m				Original assumed mass
     * @param rho 			Density
	 * @param rho2			Density of the rank that the current rank is communicating with
	 * @param rhoold		Stores the densities from the last time step
     * @param r				Distance between two particles
     * @param q 			Relative distance between two particles
     * @param phi_d 		Kernel density function for density
     */
    double m = 1;
	double* rho = new double[n];
	double* rho2;
	double* rhostore = new double[n];
    double* r = new double[2];
    double q, phi_d;
	
	/*
     * @brief Initiate parameters for pressure force calculation
	 * @param p				Pressure from ideal gas law
	 * @param p2			Pressure of the rank that the current rank is communicating with
	 * @param pold			Stores the pressures from the last time step
     * @param Fp			Pressure force of all particles
     * @param div_phi_p		1D Arrary to store div_phi_p
     * @param divtemp		2D array to store the immediate div_phi_p
     * @param Fp_i			Pressure force on particle i
     */
	double* p = new double[n];  
	double* p2;
	double* pstore = new double[n];
    double* Fp = new double[n * 2]; 
    double* div_phi_p = new double[2];
    double* divtemp = new double[2];
    double* Fp_i = new double[2];
	
	/*
     * @brief Initiate parameters for viscous force calculation
     * @param Fv			Viscous force of all particles
     * @param vij			Difference in velocities between two particles
     * @param div2_phi_v	Laplacian operator
     */
	double* Fv = new double[n * 2];
    double* vij = new double[2];
    double div2_phi_v;// = new double[N * n];
	
	/*
	 * @brief Initiate parameters for gravity force calculation
	 * @param Fg	Gravity force of all particles
	 */ 
	 double* Fg = new double[n * 2];  // Gravity force
	 
	 /*
     * @brief Initiate parameters for position calculation
     * @param a		Acceleration
     * @param v		Velocity, =0 @ t=0 for validation cases, !=0 @ t=0 for test cases
	 * @param v2	Velocity of the rank that the current rank is communicating with
	 * @param vall	Veloctiy of all particles, initialised only at rank 0
     */
    double* a = new double[n * 2];
    double* v = new double[n * 2];
	double* v2 = new double[n * 2];
	double* vall;
};