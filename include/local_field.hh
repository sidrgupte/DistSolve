#pragma once

#include <vector>

struct LocalField {
    int nx;
    int ny;
    std::vector<double> data;

    LocalField(int nx_, int ny_)
        : nx(nx_),
          ny(ny_),
          data(static_cast<std::size_t>(nx_ + 2) * (ny_ + 2), 0.0)
    {}

    double& at(int i, int j) {
        return data[static_cast<std::size_t>(i) * (ny + 2) + j];
    }

    const double& at(int i, int j) const {
        return data[static_cast<std::size_t>(i) * (ny + 2) + j];
    }
};