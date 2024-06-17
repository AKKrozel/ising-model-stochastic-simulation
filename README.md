# ising-model-stochastic-simulation

This project is a stochastic simulation of the ising model of ferromagnetism. A spin lattice is randomly updataed millions of times according to the Boltzmann Distribution at a specific temperature. In the simulation, a external magnetic field is applied in one direction and then gradually in the other directiona and back again in order to produce an effect known as hysteresis. An animation is created displaying how the spin lattice state evolves between each value of the external magnetic field along with a hysteresis curve that shows how the spin lattice's net magnetic field tends to lag behind the external magnetic field due to interaction between spins.

# Usage

The file hysteresis.cpp is used to produce output spin data that is then used by highres-fluid-animation.ipynb to create images and animations.

All of these files used should be fairly easy to run using the command line and a Jupyter Notebook. Be sure to allow the files eachother by placing the .cpp and .h file in the same directory. It will also be necessary to provide appropriate file paths in Ising_Animation.ipyb.

# Animation

To see the end result of this project's pipeline, check out Ising_Anim.mp4.
