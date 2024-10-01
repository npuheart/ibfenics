/// @date 2024-07-12
/// @file IB.cpp
/// @author Ma Pengfei (code@pengfeima.cn)
/// @version 0.1
/// @copyright Copyright (c) 2024 Ma Pengfei
/// 
/// @brief 封装IBInterpolation和IBMesh类，提供python接口
/// 
///
#include "IBInterpolation.h"
#include "IBMesh.h"

#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


using namespace dolfin;

namespace py = pybind11;
PYBIND11_MODULE(cpp, m)
{
	py::class_<double2>(m, "double2")
		.def_readonly("x", &double2::x)
		.def_readonly("y", &double2::y);

    py::class_<IBMesh, std::shared_ptr<IBMesh>>(m, "IBMesh")
        .def(py::init<std::array<dolfin::Point, 2>, std::array<size_t, 2>, size_t>())
		.def("mesh",&IBMesh::mesh)
		.def("get_hash",&IBMesh::get_hash)
 		.def("build_map", &IBMesh::build_map)
		.def("assign_value_on_dofs", &IBMesh::assign_value_on_dofs)
		.def("extract_dofs", &IBMesh::extract_dofs)
		.def("assign_dofs", &IBMesh::assign_dofs)
		.def("evaluate", &IBMesh::evaluate)
		;

    py::class_<IBInterpolation>(m, "IBInterpolation")
        .def(py::init<std::shared_ptr<IBMesh>, std::shared_ptr<Mesh>>())
		.def("solid_to_fluid_derivative", &IBInterpolation::solid_to_fluid_derivative)
		.def("solid_to_fluid", &IBInterpolation::solid_to_fluid)
		.def("fluid_to_solid", &IBInterpolation::fluid_to_solid)
		.def("evaluate_current_points", &IBInterpolation::evaluate_current_points)
		.def("evaluate_weights", &IBInterpolation::evaluate_weights)
		;
	
	#ifdef VERSION_INFO
		m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
	#else
		m.attr("__version__") = "dev";
	#endif
}