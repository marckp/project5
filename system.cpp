#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include <iostream>

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions()
{
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    double  x_size  = m_systemSize.x();
    double  y_size  = m_systemSize.y();
    double  z_size  = m_systemSize.z();
    for(Atom *atom : m_atoms)
    {
        // Update x dimension
        if (atom->position.x() <  0) atom->position.setX(atom->position.x() + x_size);
        if (atom->position.x() >= x_size) atom->position.setX(atom->position.x() - x_size);

        if (atom->position.y() <  0) atom->position.setY(atom->position.y() + y_size);
        if (atom->position.y() >= y_size) atom->position.setX(atom->position.y() - y_size);

        if (atom->position.z() <  0) atom->position.setZ(atom->position.z() + z_size);
        if (atom->position.z() >= z_size) atom->position.setZ(atom->position.z() - z_size);
    }

}

void System::removeTotalMomentum()
{
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    m_momentum.zeros();
    for(Atom *atom : m_atoms)
    {
        m_momentum += atom->mass() * atom->velocity;
    }
    vec3 scaling_factor = m_momentum/m_atoms.size();
    for(Atom *atom : m_atoms)
    {
        atom->velocity -= scaling_factor/atom->mass();
    }
    m_momentum.zeros();
    for(Atom *atom : m_atoms)
    {
        m_momentum += atom->mass() * atom->velocity;
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature)
{
    setSystemSize(vec3(0, 0, 0) + numberOfUnitCellsEachDimension * latticeConstant); // Remember to set the correct system size!
    vec3 latticePositions[4];
    for (int i=0; i < numberOfUnitCellsEachDimension; ++i)
    {
        for (int j=0; j < numberOfUnitCellsEachDimension; ++j)
        {
            for (int k=0; k < numberOfUnitCellsEachDimension; ++k)
            {
                Atom *atom1
            }
        }

    }
}

// You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10)
//for(int i=0; i<100; i++)
//{
//    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); // mu of Argon
//    // random number in the interval [0,10]
//    double x = Random::nextDouble(0, 10);
//    double y = Random::nextDouble(0, 10);
//    double z = Random::nextDouble(0, 10);
//    atom->position.set(x,y,z);
//    atom->resetVelocityMaxwellian(temperature);
//    m_atoms.push_back(atom);
//}


void System::calculateForces()
{
    for(Atom *atom : m_atoms)
    {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt)
{
    // find new partical configuration by applying forces
    m_integrator.integrate(*this, dt);

    // Apply periodic boundary conditions
    applyPeriodicBoundaryConditions();

    m_steps++;
    m_time += dt;
}
