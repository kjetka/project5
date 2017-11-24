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

int main(){


    // Initial values setting up system
    int nrUnitCellsEachDirection =5;

    vector<double> Temperatures_si = {50.0,85.0,300.0,500.0};
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms
    //double sigma = UnitConverter::lengthFromAngstroms(3.405)
    int timeLimit = 1e4;
    //IF we are using the command line for input variables:

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.
    /*
    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
*/


    for(int temperature_current:Temperatures_si){

        double initialTemperature = UnitConverter::temperatureFromSI(temperature_current); // measured in Kelvin

        // setting up system
        System system(nrUnitCellsEachDirection);

        system.createFCCLattice(latticeConstant, initialTemperature);
        //system.potential().setEpsilon(1.0);

        //system.potential().setSigma(1.0);
        system.potential().setEpsilon(UnitConverter::temperatureFromSI(119.8));
        system.potential().setSigma(UnitConverter::lengthFromAngstroms(3.405));
        system.removeTotalMomentum();

        StatisticsSampler statisticsSampler;
        string movietitle = "../results/movie_T_"+to_string(temperature_current)+".xyz";
        IO movie(movietitle.c_str());

        /*cout << setw(20) << "Timestep" <<
                setw(20) << "Time" <<
                setw(20) << "Temperature" <<
                setw(20) << "KineticEnergy" <<
                setw(20) << "PotentialEnergy" <<
                setw(20) << "TotalEnergy" << endl;
        */

        for(int timestep=0; timestep<timeLimit; timestep++) {
            system.step(dt);
            statisticsSampler.sample(system); // system - same as *this within a object.
            if( timestep % 50 == 0 ) {
                // Print the timestep every 100 timesteps
                /*cout << setw(20) << system.steps() <<
                        setw(20) << system.time() <<
                        setw(20) << statisticsSampler.temperature() <<
                        setw(20) << statisticsSampler.kineticEnergy() <<
                        setw(20) << statisticsSampler.potentialEnergy() <<
                        setw(20) << statisticsSampler.totalEnergy() << endl;
                */

            movie.saveState(system);
            }
        }
        //cout << "check if applyPeriodicBoundaryConditions works for diffusion"<< endl;
        movie.close();
    }
    return 0;
}
