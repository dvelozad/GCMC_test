#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <unordered_map>

struct Vec3 {
    double x, y, z;
};

// Apply periodic boundary conditions
void pbc(Vec3 &v, double box) {
    if (v.x < 0) v.x += box; else if (v.x >= box) v.x -= box;
    if (v.y < 0) v.y += box; else if (v.y >= box) v.y -= box;
    if (v.z < 0) v.z += box; else if (v.z >= box) v.z -= box;
}

// Minimum image distance squared
double dist2(const Vec3 &a, const Vec3 &b, double box) {
    double dx = a.x - b.x; dx -= box * std::round(dx / box);
    double dy = a.y - b.y; dy -= box * std::round(dy / box);
    double dz = a.z - b.z; dz -= box * std::round(dz / box);
    return dx*dx + dy*dy + dz*dz;
}

// Lennard-Jones potential in reduced units with cutoff and shift
double lj_potential(double r2, double cutoff2) {
    if (r2 >= cutoff2) return 0.0;
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    double inv_r12 = inv_r6 * inv_r6;
    double vc = 4.0 * (inv_r12 - inv_r6);
    return vc;
}

class GCMC {
public:
    // Simulation parameters
    double box_length = 10.0; // reduced units
    double temperature = 1.0; // epsilon/k_B
    double mu = -3.0;         // reduced chemical potential
    int steps = 1000;
    double max_disp = 0.1;    // max displacement for translations
    double cutoff = 2.5;      // LJ cutoff
    unsigned int seed = 1234;

    std::vector<Vec3> positions;
    std::mt19937 rng;
    std::uniform_real_distribution<double> uniform01{0.0,1.0};

    void load(const std::string &filename) {
        std::ifstream f(filename);
        if (!f) {
            throw std::runtime_error("Cannot open input file");
        }
        std::string key; double val;
        while (f >> key >> val) {
            if (key == "box_length") box_length = val;
            else if (key == "temperature") temperature = val;
            else if (key == "mu") mu = val;
            else if (key == "steps") steps = static_cast<int>(val);
            else if (key == "max_disp") max_disp = val;
            else if (key == "cutoff") cutoff = val;
            else if (key == "seed") seed = static_cast<unsigned int>(val);
        }
    }

    double beta() const { return 1.0 / temperature; }

    double total_energy() const {
        double e = 0.0;
        double cutoff2 = cutoff * cutoff;
        for (size_t i=0; i<positions.size(); ++i) {
            for (size_t j=i+1; j<positions.size(); ++j) {
                double r2 = dist2(positions[i], positions[j], box_length);
                e += lj_potential(r2, cutoff2);
            }
        }
        return e;
    }

    double particle_energy(const Vec3 &pos) const {
        double e = 0.0;
        double cutoff2 = cutoff * cutoff;
        for (const auto &p : positions) {
            double r2 = dist2(pos, p, box_length);
            e += lj_potential(r2, cutoff2);
        }
        return e;
    }

    void initialize() {
        rng.seed(seed);
    }

    void run() {
        initialize();
        double beta_mu = beta() * mu;
        std::uniform_real_distribution<double> posdist(0.0, box_length);

        for (int step = 0; step < steps; ++step) {
            double move_type = uniform01(rng);

            if (move_type < 0.33) { // insertion
                Vec3 new_pos{posdist(rng), posdist(rng), posdist(rng)};
                double dE = particle_energy(new_pos);
                double prob = std::exp(beta_mu - beta() * dE) * box_length * box_length * box_length / static_cast<double>(positions.size() + 1);
                if (uniform01(rng) < prob) {
                    positions.push_back(new_pos);
                }
            } else if (move_type < 0.66) { // deletion
                if (!positions.empty()) {
                    std::uniform_int_distribution<int> pick(0, positions.size()-1);
                    int idx = pick(rng);
                    Vec3 pos = positions[idx];
                    double dE = particle_energy(pos);
                    double prob = std::exp(-beta_mu + beta() * dE) * static_cast<double>(positions.size()) / (box_length * box_length * box_length);
                    if (uniform01(rng) < prob) {
                        positions.erase(positions.begin() + idx);
                    }
                }
            } else { // translation
                if (!positions.empty()) {
                    std::uniform_int_distribution<int> pick(0, positions.size()-1);
                    int idx = pick(rng);
                    Vec3 old_pos = positions[idx];
                    Vec3 new_pos = old_pos;
                    new_pos.x += (uniform01(rng) - 0.5) * 2.0 * max_disp;
                    new_pos.y += (uniform01(rng) - 0.5) * 2.0 * max_disp;
                    new_pos.z += (uniform01(rng) - 0.5) * 2.0 * max_disp;
                    pbc(new_pos, box_length);
                    double dE = particle_energy(new_pos) - particle_energy(old_pos);
                    if (uniform01(rng) < std::exp(-beta() * dE)) {
                        positions[idx] = new_pos;
                    }
                }
            }
            // Optionally output progress every 1000 steps
            if ((step+1) % 1000 == 0) {
                std::cout << "Step " << step+1 << " N=" << positions.size() << " E=" << total_energy() << "\n";
            }
        }
    }
};

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input_file" << std::endl;
        return 1;
    }
    GCMC sim;
    try {
        sim.load(argv[1]);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    sim.run();
    return 0;
}

