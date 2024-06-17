#include "random.h"

#include <cmath>
#include <fstream>
#include <iostream>

std::vector<std::vector<double>> simulation_step(std::vector<std::vector<double>> inputState, double J, double kBT, double& M, double H, random_number_generator& rng) {
    int Nx = inputState.size() - 2;
    int Ny = inputState[0].size() - 2;
    std::vector<std::vector<double>> outputState = inputState;

    // Choose a random spin
    int row = rng.uniform_int(1, Ny + 1);
    int col = rng.uniform_int(1, Nx + 1);

    // Calculate the change in energy if the spin is flipped
    double spin = outputState[row][col];
    double delta_E = 2.0 * J * spin * (
        outputState[row + 1][col] + outputState[row - 1][col] +
        outputState[row][col + 1] + outputState[row][col - 1]) + 2.0 * H * spin;

    // flip probabalistically based on change in energy
    if (delta_E <= 0.0) {
        outputState[row][col] *= -1.0;
        M += (2.0 / (Nx * Ny)) * outputState[row][col];
    }
    else {
        double P_a = std::exp(-delta_E / kBT);
        if (P_a > rng.uniform()) {
            outputState[row][col] *= -1.0;
            M += (2.0 / (Nx * Ny)) * outputState[row][col];
        }
    }


    // Apply periodic boundary conditions
    if (row == 1) {
        outputState[Ny + 1][col] = outputState[1][col];
    }
    if (row == Ny) {
        outputState[0][col] = outputState[Ny][col];
    }
    if (col == 1) {
        outputState[row][Nx + 1] = outputState[row][1];
    }
    if (col == Nx) {
        outputState[row][0] = outputState[row][Nx];
    }

    return outputState;
}

void capture_state(std::vector<std::vector<double>> inputState, std::ofstream& output_file) {
    int Nx = inputState.size() - 2;
    int Ny = inputState[0].size() - 2;
    for (int row = 1; row <= Ny; row++) {
        for (int col = 1; col <= Nx; col++) {
            output_file << inputState[row][col];
            if (col != Nx) {
                output_file << ", ";
            }
        }
        output_file << std::endl;
    }
    output_file.close();
}

void run_simulation(double kBT, double H_low, double H_high, double H_step, int Nx, int Ny) {
    
    random_number_generator rng;

    std::ofstream output_file0("hysteresis.csv");

    double J = 1.0; //coupling constant
    int equilibrationSteps = 1000000;
    int measurementSteps = 50000;

    // Initialize state vector with all spins down
    std::vector<double> initRow(Ny + 2, -1.0);
    std::vector<std::vector<double>> state(Nx + 2, initRow);

    // Initialize magnetization
    double M = 0.0;
    for (int row = 1; row <= Ny; row++) {
        for (int col = 1; col <= Nx; col++) {
            M += state[row][col];
        }
    }
    M *= 1.0 / (Nx * Ny);

    // Sweep magnetic field from H_start to H_end and back and record average magnetization along the way
    for (double H = H_low; H <= H_high; H += H_step) {
        for (int step = 0; step < equilibrationSteps; ++step) {
            state = simulation_step(state, J, kBT, M, H, rng);
        }
        double M_avg = 0.0;
        for (int step = 0; step < measurementSteps; ++step) {
            state = simulation_step(state, J, kBT, M, H, rng);
            M_avg += M;
        }
        M_avg /= measurementSteps;
        output_file0 << H << ", " << M_avg << std::endl;

        std::ofstream output_file1("Ising_output_grid_up" + std::to_string(H) + ".csv");
        capture_state(state, output_file1);
        std::cout << std::to_string(H) << std::endl;

    }
    for (double H = H_high; H >= H_low; H -= H_step) {
        for (int step = 0; step < equilibrationSteps; ++step) {
            state = simulation_step(state, J, kBT, M, H, rng);
        }
        double M_avg = 0.0;
        for (int step = 0; step < measurementSteps; ++step) {
            state = simulation_step(state, J, kBT, M, H, rng);
            M_avg += M;
        }
        M_avg /= measurementSteps;
        output_file0 << std::to_string(H) << ", " << std::to_string(M_avg) << std::endl;

        std::ofstream output_file1("Ising_output_grid_down" + std::to_string(H) + ".csv");
        capture_state(state, output_file1);
        std::cout << std::to_string(H) << std::endl;

    }
    output_file0.close();
}

int main() {

    double kBT = 1.5; // Product of Boltmann constant with temperature

    // spin grid dimensions
    int Nx = 200;
    int Ny = 200;

    // behavior of the external magnetic field
    double H_low = -0.35;  
    double H_high = 0.35;     
    double H_step = 0.001;

    run_simulation(kBT, H_low, H_high, H_step, Nx, Ny);

    return 0;
}

