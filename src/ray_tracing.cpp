#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <random>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <chrono>

using namespace std;

struct Point {
    // Represents a point in 3D space (x, y, z)
    double x;
    double y;
    double z;
};


double** create_window(int const& N) {
    double** matrix = new double*[N];

    for (int i = 0; i < N; i++) {
        matrix[i] = new double[N];

        for (int j = 0; j < N; j++) {
            matrix[i][j] = 0.0;
        }
    }

    return matrix;
}

Point direction_sample() {
    Point V;
    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr_psi(0.0, 2.0*M_PI);
    uniform_real_distribution<double> distr_cos_theta(-1.0, 1.0);

    double psi = distr_psi(eng);
    double cos_theta = distr_cos_theta(eng);
    double sin_theta = sqrt(1 - pow(cos_theta, 2));
    double cos_psi = cos (psi); 
    double sin_psi = sin (psi);
    V.x = sin_theta * cos_psi;
    V.y = sin_theta * sin_psi;
    V.z = cos_theta;

    return V;
}

Point arithmetic_vectors(Point& vector_one, Point& vector_two, bool const& flag) {
    // flag == true => addition
    Point output;
    if (flag) {
        output.x = vector_one.x + vector_two.x;
        output.y = vector_one.y + vector_two.y;
        output.z = vector_one.z + vector_two.z;
    } else {
        output.x = vector_one.x - vector_two.x;
        output.y = vector_one.y - vector_two.y;
        output.z = vector_one.z - vector_two.z;
    }

    return output;
}

void mult_vec_scal(Point& output, Point& vector, double const& scalar) {
    output.x = vector.x * scalar;
    output.y = vector.y * scalar;
    output.z = vector.z * scalar;
}

void divide_vec_scal(Point& output, Point& vector, double const& scalar) {
    output.x = vector.x / scalar;
    output.y = vector.y / scalar;
    output.z = vector.z / scalar;
}

double dotp(Point& vector_one, Point& vector_two) {
    return (vector_one.x * vector_two.x) + (vector_one.y * vector_two.y) + (vector_one.z * vector_two.z);
}

double magnitude(Point& vector) {
    return sqrt (pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2));
}

void write_to_file(double** const& G, int const& window_size) {
    ofstream myfile;

    myfile.open("../out/sphere.txt");

    myfile << "[" << window_size << "], ";

    myfile << "[";
    for (int i = 0; i < window_size; i++) {
        if (i == 0) {
            myfile << "[";
        } else {
            myfile << ", [";
        }
        for (int j = 0; j < window_size; j++) {
            if (j == window_size - 1) {
                myfile << G[i][j];
            } else {
                myfile << G[i][j] << ", ";
            }
        }
        myfile << "]";
    }
    myfile << "]";
    myfile.close();
}

int main() {
    // allocate G[1...n][1...n]
    double window_size = 1000;
    double num_rays = 1e8;
    // radius of sphere
    double R = 6;

    double psi;
    double cos_theta;
    double sin_theta;

    double** G = create_window(window_size);

    Point V;
    Point W;
    Point N;
    Point I;
    Point S;
    
    // GIVEN - CHANGE TO INPUT ARGS
    double W_max = 10.0;
    // distance from window to object
    double W_y = 10.0;
    Point C = {0, 12, 0};
    Point L = {4, 4, -1};
    for (int i = 1; i < num_rays + 1; i++) {
        auto start = std::chrono::steady_clock::now();
        while (true) {
            // Sample random V from unit sphere
            V = direction_sample();
            mult_vec_scal(W, V, (W_y/V.y));

            if (abs (W.x) < W_max && abs (W.z) < W_max && pow(dotp(V, C), 2) + pow(R, 2) - dotp(C, C) > 0) {
                break;
            }
        };

        double t = dotp(V, C) - sqrt (pow(dotp(V, C), 2) + pow(R, 2) - dotp(C, C));
        mult_vec_scal(I, V, t);
        
        Point I_minus_C = arithmetic_vectors(I, C, false);
        divide_vec_scal(N, I_minus_C, magnitude(I_minus_C));
        
        Point L_minus_I = arithmetic_vectors(L, I, false);
        divide_vec_scal(S, L_minus_I, magnitude(L_minus_I));

        double b = max(0.0, dotp(S, N));

        int local_i = ((W.z + W_max) / (2*W_max)) * window_size;
        int local_j = ((W.x + W_max) / (2*W_max)) * window_size;

        G[local_i][local_j] += b;

        auto end = std::chrono::steady_clock::now();
	    std::chrono::duration<double> elapsed_time = end - start;
	    double duration = elapsed_time.count();
	    std::cout << "ray " << i << " -> " << 1/duration << " rays/s" << std::endl;
    }

    write_to_file(G, window_size);
}