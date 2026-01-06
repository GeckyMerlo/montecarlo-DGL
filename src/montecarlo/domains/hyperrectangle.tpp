#include <cmath>
#include "../geometry.hpp"
#include "integration_domain.hpp"
#include <array>
#include <utility>

using namespace geom;
using namespace std;

template <size_t dim>
HyperRectangle<dim>::HyperRectangle(array<double, dim> &dims):
    dimensions(dims)
{};


template <size_t dim>
auto HyperRectangle<dim>::getBounds() const -> Bounds<dim> {
    Bounds<dim> bounds;
    for(size_t i = 0; i < dim; ++i)
        bounds[i] = make_pair(-dimensions[i]/2, dimensions[i]/2);
    return bounds;
}

template <size_t dim>
double HyperRectangle<dim>::getBoxVolume() const{
    double volume = 1;
    for(size_t i = 0; i < dim; ++i) {
        volume = volume * dimensions[i];
    }
    return volume;
}

template <size_t dim>
bool HyperRectangle<dim>::isInside(const Point<dim> &point) const{
    bool inside = true;
        for(size_t i = 0; i < dim; ++i) {
            if (point[i] < -dimensions[i]/2 || point[i] > dimensions[i]/2) {
                inside = false;
            }
        }
    return inside;
}