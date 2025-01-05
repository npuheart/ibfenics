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
#include "main.h"
#include "Poisson2D.h"
#include "VTUIO.h"
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

class FacetIntegration
{
private:
    std::vector<std::size_t> _facets;
    const size_t num_points = 1;
    const std::vector<std::vector<double>> some_points_all = {{0.5, 0.5}, {0.0, 0.5}, {0.5, 0.0}};

public:
    // Find all id of facets on the given sub_domain
    FacetIntegration(const Mesh &mesh, const MeshFunction<std::size_t> &bdry, size_t sub_domain)
    {
        // Make sure we have the facet - cell connectivity
        const std::size_t D = mesh.topology().dim();
        mesh.init(D - 1); // Is it neccessary ?
        mesh.init(D);
        mesh.init(D - 1, D);

        // Build set of boundary facets
        // TODO: find facets on the boundary
        dolfin_assert(_facets.empty());
        for (FacetIterator facet(mesh); !facet.end(); ++facet)
        {
            if (bdry[*facet] == sub_domain)
                _facets.push_back(facet->index());
        }
    }

    auto fun3(const std::shared_ptr<Function> &disp, const std::shared_ptr<Function> &force)
    {
        // TODO: check if disp and force are on the same function space
        auto _function_space = disp->function_space();
        const Mesh &mesh = *_function_space->mesh();
        const std::size_t D = mesh.topology().dim();
        const GenericDofMap &dofmap = *_function_space->dofmap();
        const FiniteElement &element = *_function_space->element();
        const ufc::finite_element &ufc_element = *element.ufc_element();

        ufc::cell ufc_cell;
        auto value_size = disp->value_size();
        auto space_dimension = element.space_dimension();
        std::vector<double> basis_values(value_size * space_dimension);
        std::vector<double> coefficients_disp(element.space_dimension());
        std::vector<double> coefficients_force(element.space_dimension());

        // to hold the results
        std::vector<double> values_disp_all(_facets.size() * num_points * value_size);
        std::vector<double> values_force_all(_facets.size() * num_points * value_size);

        for (std::size_t f = 0; f < _facets.size(); ++f)
        {
            // Create facet
            const Facet facet(mesh, _facets[f]);

            // Get cell to which facet belongs.
            dolfin_assert(facet.num_entities(D) > 0);
            const std::size_t cell_index = facet.entities(D)[0];
            const Cell cell(mesh, cell_index);

            // Get local index of facet with respect to the cell
            const std::size_t local_facet = cell.index(facet);
            // TODO: Get Gauss Quadrature Rule for local facet.
            const auto some_points = some_points_all[local_facet];
            printf("local_facet: %ld\n", local_facet);

            auto dofs = dofmap.cell_dofs(cell_index);
            disp->vector()->get_local(coefficients_disp.data(), dofs.size(), dofs.data());
            force->vector()->get_local(coefficients_force.data(), dofs.size(), dofs.data());

            std::vector<double> vertex_coordinates(6);
            cell.get_vertex_coordinates(vertex_coordinates);
            printf("vertex_coordinates: %f %f, %f %f, %f %f\n",
                   vertex_coordinates[0], vertex_coordinates[1], vertex_coordinates[2],
                   vertex_coordinates[3], vertex_coordinates[4], vertex_coordinates[5]);

            std::vector<double> some_points_local(some_points.size());
            linearInterpolation<2>(some_points_local.data(), some_points.data(), vertex_coordinates.data(), num_points);
            std::vector<double> ref_vertex_basis_values(2* space_dimension * num_points);
            ufc_element.evaluate_reference_basis(ref_vertex_basis_values.data(), num_points, some_points.data());

            std::vector<double> values_disp(value_size * num_points);
            std::vector<double> values_force(value_size * num_points);

            for (size_t k = 0; k < num_points; k++)
            {
                for (std::size_t i = 0; i < space_dimension; ++i)
                {
                    for (std::size_t j = 0; j < value_size; ++j)
                    {
                        values_disp[2 * k + j] += coefficients_disp[i] * ref_vertex_basis_values[2* space_dimension* k + 2 * i + j];
                        values_force[2 * k + j] += coefficients_force[i] * ref_vertex_basis_values[2* space_dimension* k + 2 * i + j];
                    }
                }
            }

            for (size_t k = 0; k < num_points; k++)
            {
                for (std::size_t j = 0; j < value_size; ++j)
                {
                    values_disp_all[f * num_points * value_size + value_size * k + j] = values_disp[value_size * k + j];
                    values_force_all[f * num_points * value_size + value_size * k + j] = values_force[value_size * k + j];
                }
            }
            // for (size_t i = 0; i < values.size(); i++)
            // {
            //     printf("%f ", values[i]);
            // }
            // printf("\n");

            for (size_t i = 0; i < 2 * num_points; i++)
            {
                printf("%.e ", some_points_local[i]);
            }
            printf("\n");
                        for (size_t i = 0; i < 2 * num_points; i++)
            {
                printf("%.e ", values_disp[i]);
            }
            printf("\n");
                        for (size_t i = 0; i < 2 * num_points; i++)
            {
                printf("%.e ", values_force[i] - values_disp[i]*values_disp[i]);
            }
            printf("\n");
            // for (size_t i = 0; i < 2 * num_points; i++)
            // {
            //     printf("%f ", some_points[i]);
            // }
            // printf("\n");
        }
        return std::make_pair(values_disp_all, values_force_all);
    }

    void fun4(const std::shared_ptr<Function> &disp, const std::shared_ptr<Function> &force)
    {
        auto ret = fun3(disp, force);
        std::vector<std::array<double, 3>> points0;
        std::vector<std::array<double, 3>> points;
        std::vector<std::array<double, 3>> velocities;

        size_t num_points = ret.first.size() / disp->value_size();
        for (size_t i = 0; i < num_points; i++)
        {
            points0.push_back({ret.first[2 * i], ret.first[2 * i + 1], 0.0});
            points.push_back({ret.first[2 * i], ret.first[2 * i + 1], 0.0});
            velocities.push_back({ret.first[2 * i], ret.first[2 * i + 1], 0.0});
        }
        std::vector<std::array<double, 3>> data1 = points;
        std::vector<std::array<double, 3>> data2 = velocities;

        std::vector<std::string> tags{"position"};
        IO::write_particles_to_vtu("a.vtu", points, tags, data2);
    }
    ~FacetIntegration()
    {
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

    FacetIntegration facet_integration(*mesh, *bdry, 1);

    File file("force.pvd");
    file << *force;

    File file_bdry("bdry.pvd");
    file_bdry << *bdry;

    facet_integration.fun4(disp, force);

    std::vector<double> some_points = {0.5, 0.5, 0.1, 0.2};
    auto result = vectors_to_strings(some_points);
    std::cout << result[0] << std::endl;

    return 0;
}
