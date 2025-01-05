/// @date 2025-01-02
/// @file main.cpp
/// @author Ma Pengfei (code@pengfeima.cn)
/// @version 0.1
/// @copyright Copyright (c) 2025 Ma Pengfei
/// 
/// @brief 
/// 
///

#include <dolfin.h>
#include "Poisson2D.h"
#include "FacetIntegration.h"
using namespace dolfin;

class DirichletBoundary : public SubDomain
{
    bool inside(const Array<double> &x, bool on_boundary) const
    {
        return on_boundary;
    }
};

class Force : public Expression
{
public:
    Force() : Expression(2) {}

    void eval(Array<double> &values, const Array<double> &x) const
    {
        values[0] = x[0];
        values[1] = x[1];
    }
};
class Force2 : public Expression
{
public:
    Force2() : Expression(2) {}

    void eval(Array<double> &values, const Array<double> &x) const
    {
        values[0] = x[0]*x[0];
        values[1] = x[1]*x[1];
    }
};

int main_2()
{
    // auto mesh = std::make_shared<Mesh>(
    //     UnitSquareMesh::create({{32, 32}}, CellType::Type::triangle));

    auto mesh = std::make_shared<dolfin::Mesh>();
    {
        dolfin::XDMFFile mesh_file_1("/home/fenics/mesh/benchmark/cylinder/cylinder_40.xdmf");
        mesh_file_1.read(*mesh);
        mesh_file_1.close();
    }

    auto bdry = std::make_shared<MeshFunction<std::size_t>>(
        mesh, mesh->topology().dim() - 1, 0);

    auto V = std::make_shared<Poisson2D::FunctionSpace>(mesh);

    auto disp = std::make_shared<Function>(V);
    auto force = std::make_shared<Function>(V);
    disp->interpolate(Force());
    force->interpolate(Force2());

    DirichletBoundary dirichlet_boundary;

    dirichlet_boundary.mark(*bdry, 1);

    FacetIntegration facet_integration(mesh, bdry, 1);

    File file("force.pvd");
    file << *force;

    File file_bdry("bdry.pvd");
    file_bdry << *bdry;

    facet_integration.fun4(disp, force);

    return 0;
}
