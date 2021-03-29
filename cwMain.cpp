/* Imperial College London
 * AERO-96021 HPC Coursework Assignment 2021
 * Kin Pui Tam - 01338789
 * Date of submission - 24/03/2021
 *
 * Coursework main script
 *
 * This program allows the user to simulate smoothed particle hydrodynamics (SPH) via the terminal.
 * The implementation of the simulation is achieved by the SPH class, which consists of essential parameters
 * and several member functions.
 */

// Boost external library for user inputs in command line
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Include SPH class header
#include "sph.h"

int main(int argc, char* argv[])
{
    // ---------------------- //
    // -----MPI Settings----- //
    // ---------------------- //

    // Test MPI return values
    int retVal;
    retVal = MPI_Init(&argc, &argv);
    if(retVal != MPI_SUCCESS)
	{
	    cout << "An error occured initialising MPI" << endl;
	    return -1;
	}

    // Retrieve communicator's rank and size
    int rank, size, retRank, retSize;
    retRank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    retSize = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(retRank == MPI_ERR_COMM || retSize == MPI_ERR_COMM)
	{
	    cout << "Invalid communicator" << endl;
	    MPI_Finalize();
	    return 1;
	}

    // -------------------------------------------------- //
    // -----SIMULATION AND TIME INTEGRATION SETTINGS----- //
    // -------------------------------------------------- //

    // Simulation and time integration settings
    unsigned int N = 0; 	// Initialise number of particles
    double dt = 0.0001; 	// Default time step size
    double T = 2; 			// Default total integration time
    double h = 0.01; 		// Default radius of interaction
	
	// Console display settings
	bool showx = false;		// Display positions at each time step
	bool showv = false;		// Display velocities at each time step
	bool showe = false;		// Display energies at each time step
	string trash;			// Trash the value the user entered after confirming change of settings 

    // ------------------------ //
    // -----BOOST SETTINGS----- //
    // ------------------------ //

    // Use boost_program_options external libraries to allow command line input
    try
	{
	    po::options_description desc("Allowed options");
	    desc.add_options()
			("help", "Produce help message.")
			("ic-one-particle", "Use one particle validation case initial condition.")
			("ic-two-particles", "Use two particles validation case initial condition.")
			("ic-three-particles", "Use three particles validation case initial condition.")
			("ic-four-particles", "Use four particles validation case initial condition.")
			("ic-dam-break", "Use dam-break initial condition.")
			("ic-block-drop", "Use block-drop initial condition.")
			("ic-droplet", "Use droplet initial condition.")
			("dt", po::value<double>(), "Change time-step to use.")
			("T", po::value<double>(), "Change total integration time.")
			("h", po::value<double>(), "Change radius of influence of each particle.")
			("showx", "Show position of each particle at all time steps, validation cases only")
			("showv", "Show velocity of each particle at all time steps, validation cases only")
			("showe", "Show kinectic energy and potential energy at all time steps")
			("N", po::value<double>(), "[Dev] Change number of particles N directly");

	    po::variables_map vm;
	    po::store(po::parse_command_line(argc, argv, desc), vm);
	    po::notify(vm);

	    if(vm.count("help"))
		{
		    if(rank == 0)
			cout << desc << endl;
		    MPI_Finalize();
		    return 0;
		}
		
	    if(vm.count("ic-one-particle"))
		{
		    N = 1;
		}
		
	    if(vm.count("ic-two-particles"))
		{
		    N = 2;
		}
		
	    if(vm.count("ic-three-particles"))
		{
		    N = 3;
		}
		
	    if(vm.count("ic-four-particles"))
		{
		    N = 4;
		}
		
	    if(vm.count("ic-dam-break"))
		{
		    N = 400;
		}
		
	    if(vm.count("ic-block-drop"))
		{
		    N = 651;
		}
		
	    if(vm.count("ic-droplet"))
		{
		    N = 311;
		}
		
	    if(vm.count("dt"))
		{
		    if(rank == 0)
			{
			    cout << "Time step size dt [s] is set from" << endl;
				cout.precision(4);
			    cout << setw(6) << fixed << right << dt << setw(6) << "to ";
			}
		    dt = vm["dt"].as<double>();
		    if(rank == 0)
			{
			    cout << setw(6) << fixed << dt << endl;
			}
		}
		
	    if(vm.count("T"))
		{
		    if(rank == 0)
			{
			    cout << "Total integration time T [s] is set from" << endl;
				cout.precision(4);
			    cout << setw(6) << fixed << right << T << setw(6) << "to ";
			}
		    T = vm["T"].as<double>();
		    if(rank == 0)
			{
			    cout << setw(6) << fixed << T << endl;
			}
		}
		
	    if(vm.count("h"))
		{
		    if(rank == 0)
			{
			    cout << "Radius of interaction h [m] is set from" << endl;
				cout.precision(4);
			    cout << setw(6) << fixed << right << h << setw(6) << "to ";
			}
		    h = vm["h"].as<double>();
		    if(rank == 0)
			{
			    cout << setw(6) << fixed << h << endl;
			}
		}
		
	    if(vm.count("N"))
		{
		    if(rank == 0)
			{
			    cout << "Number of particles N is set from" << endl;
			    cout << setw(6) << right << N << setw(6) << "to ";
			}
		    N = vm["N"].as<double>();
		    if(rank == 0)
			{
			    cout << setw(6) << N << endl;
			}
		}
		
		if((vm.count("T")) | (vm.count("dt")) | (vm.count("h")) | (vm.count("N")))
		{
			if(rank == 0)
			{
				cout << "Press ENTER to continue..." << endl;
				getline(cin, trash);
			}
		}
		
		if(vm.count("showx"))
		{
			showx = true;
		}
		if(vm.count("showv"))
		{
			showv = true;
		}
		if(vm.count("showe"))
		{
			showe = true;
		}
	}
    // Output messages if user input is invalid
    catch(exception& e)
	{
	    if(rank == 0)
		cerr << "Error: " << e.what() << endl;
	}
    catch(...)
	{
	    if(rank == 0)
		cerr << "Excpetion of unknown type!\n";
	}
    if(rank == 0)
	cout << endl;

    // -------------------------- //
    // -----INPUT VALIDATION----- //
    // -------------------------- //

    // Output error message if user didnt not enter validation/test case
    if(N == 0)
	{
	    if(rank == 0)
		{
		    cout << "-------------------------------------ERROR--------------------------------------" << endl;
		    cout << "Did not enter a validation or test case for simuation." << endl;
		    cout << endl;
		}
	}
    // Output error message if number of particles is smaller than number of processes, N < size
    else if(int(N) < size)
	{
	    if(rank == 0)
		{
		    cout << "-------------------------------------ERROR--------------------------------------" << endl;
		    cout << "Number of particles cannot be smaller than the number of allocated processors." << endl;
		    cout << "Current number of particles N = " << N << endl;
		    cout << "Current number of processes = " << size << endl;
		    cout << endl;
		}
	}
    else
	{
	    SPH sph = SPH(N, dt, T, h, rank, size, showx, showv, showe);
	}

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
