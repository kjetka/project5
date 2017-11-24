#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(int numberOfArguments, char **argumentList){


    // Initial values setting up system
    int nrUnitCellsEachDirection =3;
    double initialTemperature = UnitConverter::temperatureFromSI(300.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms
    //double sigma = UnitConverter::lengthFromAngstroms(3.405)
    int timeLimit = 1e5;
    //IF we are using the command line for input variables:
 /*
    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2)
        initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3)
        latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));
*/

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.
    /*
    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
*/
    // setting up system
    System system(nrUnitCellsEachDirection);

    system.createFCCLattice(latticeConstant, initialTemperature);
    //system.potential().setEpsilon(1.0);

    //system.potential().setSigma(1.0);
    system.potential().setEpsilon(UnitConverter::temperatureFromSI(119.8));
    system.potential().setSigma(UnitConverter::lengthFromAngstroms(3.405));
    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;
    IO movie("../results/movie_c.xyz"); // To write the state to file. here: ofstream "../results/movie.xyz"

    cout << setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" << endl;

    for(int timestep=0; timestep<timeLimit; timestep++) {
        system.step(dt);
        statisticsSampler.sample(system); // system - same as *this within a object.
        if( timestep % 100 == 0 ) {
            // Print the timestep every 100 timesteps
            propertiesFile << setw(20) << system.steps() <<
                    setw(20) << system.time() <<
                    setw(20) << statisticsSampler.temperature() <<
                    setw(20) << statisticsSampler.kineticEnergy() <<
                    setw(20) << statisticsSampler.potentialEnergy() <<
                    setw(20) << statisticsSampler.totalEnergy() << endl;
        }
        movie.saveState(system);
    }
    cout << "check if applyPeriodicBoundaryConditions works for diffusion"<< endl;
    movie.close();

    return 0;
}
