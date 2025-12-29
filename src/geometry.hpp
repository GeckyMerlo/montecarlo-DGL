#ifndef MONTECARLO_1_GEOMETRY_HPP
#define MONTECARLO_1_GEOMETRY_HPP


#include <array>
#include <cstddef>

using namespace std;
namespace geom {

template <int dim>
class Point
{
public:
    Point() { coords.fill(0.0); }

    double& operator[](std::size_t i)       { return coords[i]; }
    double  operator[](std::size_t i) const { return coords[i]; }

    std::size_t dimension() const { return dim; }

private:
    std::array<double, dim> coords;
};

template <int dim>
class Bounds
{
public:
    Bounds()
    {
        for (std::size_t i = 0; i < dim; ++i)
        {
            bounds[i] = std::make_pair(0.0, 0.0);
        }
    }

    const std::pair<double, double> &operator[](std::size_t i) const 
    {
        return bounds[i];
    }

    std::pair<double, double> &operator[](std::size_t i)
    {
        return bounds[i];
    }

private:
    array<pair<double, double>, dim> bounds;
};
                                                                                        

} // namespace geom

#endif // MONTECARLO_1_GEOMETRY_HPP