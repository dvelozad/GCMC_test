# GCMC Example

This repository contains a simple Grand Canonical Monte Carlo (GCMC) simulation for Lennard-Jones particles written in C++. The code uses reduced Lennard-Jones units and reads all simulation parameters from a plain text input file.

## Building

Compile using a C++11 compiler. For example:

```bash
g++ -std=c++11 gcmc.cpp -o gcmc
```

## Running

Provide an input file describing the simulation parameters. An example is provided in `example_input.txt`.

```bash
./gcmc example_input.txt
```

During the run the program prints progress information every 1000 Monte Carlo steps showing the current number of particles and system energy.
