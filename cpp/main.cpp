#include <dolfin.h>
#include "main.h"
#include "Poisson2D.h"
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

class FacetIntegration
{
private:
    std::vector<std::size_t> _facets;

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

    void fun1(std::shared_ptr<Function> disp)
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
        std::vector<double> coefficients(element.space_dimension());

        // NOTE:For a given cell
        size_t cell_index = 9;
        const Cell cell(mesh, cell_index);
        auto dofs = dofmap.cell_dofs(cell_index);
        disp->vector()->get_local(coefficients.data(), dofs.size(), dofs.data());

        std::vector<double> vertex_coordinates(6);
        cell.get_vertex_coordinates(vertex_coordinates);
        printf("vertex_coordinates: %f %f %f %f %f %f\n",
               vertex_coordinates[0], vertex_coordinates[1], vertex_coordinates[2],
               vertex_coordinates[3], vertex_coordinates[4], vertex_coordinates[5]);

        // NOTE:For a given point
        const size_t num_points = 2;
        const std::vector<double> some_points = {0.5, 0.5, 0.1, 0.2};
        std::vector<double> some_points_local(some_points.size());
        linearInterpolation<2>(some_points_local.data(), some_points.data(), vertex_coordinates.data(), num_points);
        std::vector<double> ref_vertex_basis_values(12 * num_points);
        ufc_element.evaluate_reference_basis(ref_vertex_basis_values.data(), num_points, some_points.data());

        std::vector<double> values(value_size * num_points);

        for (size_t k = 0; k < num_points; k++)
        {
            for (std::size_t i = 0; i < space_dimension; ++i)
            {
                for (std::size_t j = 0; j < value_size; ++j)
                {
                    values[2 * k + j] += coefficients[i] * ref_vertex_basis_values[12 * k + 2 * i + j];
                }
            }
        }

        for (size_t i = 0; i < values.size(); i++)
        {
            printf("%f ", values[i]);
        }
        printf("\n");

        for (size_t i = 0; i < 2 * num_points; i++)
        {
            printf("%f ", some_points_local[i]);
        }
        printf("\n");
        for (size_t i = 0; i < 2 * num_points; i++)
        {
            printf("%f ", some_points[i]);
        }
        printf("\n");
    }

    void fun(std::shared_ptr<Function> disp, std::shared_ptr<Function> force)
    {
        // TODO: check if disp and force are on the same function space
        auto _function_space = disp->function_space();
        const Mesh &mesh = *_function_space->mesh();
        const GenericDofMap &dofmap = *_function_space->dofmap();
        const FiniteElement &element = *_function_space->element();

        // TODO: find facets on the boundary
        // init_facets(mesh.mpi_comm());

        const std::size_t D = mesh.topology().dim();
        mesh.init(D - 1); // Is it neccessary ?
        mesh.init(D);
        mesh.init(D - 1, D);

        ufc::cell ufc_cell;
        std::vector<double> coordinate_dofs;

        auto value_size = disp->value_size();
        auto space_dimension = element.space_dimension();
        std::vector<double> basis_values(value_size * space_dimension);

        for (std::size_t f = 0; f < _facets.size(); ++f)
        {
            const Facet facet(mesh, _facets[f]);
            const std::size_t cell_index = facet.entities(D)[0];
            const Cell cell(mesh, cell_index);
            const size_t facet_local_index = cell.index(facet);
            // TODO: get the Gauss quadrature points on local facet

            // TODO: evaluate disp and force at the quadrature points
            cell.get_coordinate_dofs(coordinate_dofs);
            cell.get_cell_data(ufc_cell, facet_local_index);

            std::vector<double> w;
            disp->restrict(w.data(), element, cell, coordinate_dofs.data(), ufc_cell);
            auto cell_dofs = dofmap.cell_dofs(cell.index());

            // TODO: calculate values on gauss points
            // element->evaluate_basis_all(
            //     basis_values.data(),
            //     point.coordinates(),
            //     coordinate_dofs.data(),
            //     ufc_cell.orientation);
            // 计算出函数值

            // TODO: disp->vector()->getitem(0), basis_values, cell_dofs
            // 遍历层数:
            //   space_dimension
            //   num_gauss_points
            //   value_size
            std::vector<double> disp_values;
            size_t num_gauss_points = 10;
            for (size_t i = 0; i < space_dimension; i++)
            {
                for (size_t j = 0; j < num_gauss_points; j++)
                {
                    for (size_t k = 0; k < value_size; k++)
                    {
                        disp_values[D * j + k] = basis_values[i] * disp->vector()->getitem(D * cell_dofs[i] + k);
                    }
                }
            }
        }
    }

    std::pair<std::vector<double>, std::vector<std::size_t>> basis_values_gauss(
        const Point &point,
        const Cell &cell,
        const Function &f)
    {
        auto element = f.function_space()->element();

        auto value_size = f.value_size();
        auto space_dimension = element->space_dimension();
        std::vector<double> basis_values(value_size * space_dimension);

        ufc::cell ufc_cell;
        cell.get_cell_data(ufc_cell);

        std::vector<double> coordinate_dofs;
        cell.get_coordinate_dofs(coordinate_dofs);

        element->evaluate_basis_all(
            basis_values.data(),
            point.coordinates(),
            coordinate_dofs.data(),
            ufc_cell.orientation);
        /// 局部的dofmap就行了
        auto cell_dofmap = f.function_space()->dofmap()->cell_dofs(cell.index());
        std::vector<size_t> cell_dofmap_vector;
        for (size_t i = 0; i < cell_dofmap.size(); i++)
        {
            cell_dofmap_vector.push_back(cell_dofmap[i]);
        }

        return std::make_pair(basis_values, cell_dofmap_vector);
    }

    ~FacetIntegration()
    {
    }
};

int main()
{
    auto mesh = std::make_shared<Mesh>(
        UnitSquareMesh::create({{32, 32}}, CellType::Type::triangle));

    auto bdry = std::make_shared<MeshFunction<std::size_t>>(
        mesh, mesh->topology().dim() - 1, 0);

    auto V = std::make_shared<Poisson2D::FunctionSpace>(mesh);

    auto force = std::make_shared<Function>(V);
    force->interpolate(Force());

    DirichletBoundary dirichlet_boundary;

    dirichlet_boundary.mark(*bdry, 1);

    FacetIntegration facet_integration(*mesh, *bdry, 1);

    File file("force.pvd");
    file << *force;

    File file_bdry("bdry.pvd");
    file_bdry << *bdry;

    facet_integration.fun1(force);

    return 0;
}
