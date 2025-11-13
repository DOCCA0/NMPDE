cmake_minimum_required(VERSION 3.12.0)
project(elliptic LANGUAGES CXX C)

# Set the path to the Boost library
set(Boost_DIR "/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/boost/1.76.0/lib/cmake/Boost-1.76.0")

# Set the path to the deal.II library
set(DEAL_II_DIR "/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/dealii/9.5.1/lib/cmake/deal.II")

include(common/cmake-common.cmake)

add_executable(elliptic src_elliptic/elliptic_main.cpp src_elliptic/Elliptic.cpp)
deal_ii_setup_target(elliptic)

add_executable(elliptic_1D src_elliptic_1D/elliptic_main.cpp src_elliptic_1D/Elliptic.cpp)
deal_ii_setup_target(elliptic_1D)

add_executable(parabolic src_parabolic/parabolic_main.cpp src_parabolic/Parabolic.cpp)
deal_ii_setup_target(parabolic)

add_executable(stokes src_stokes/stokes_main.cpp src_stokes/Stokes.cpp)
deal_ii_setup_target(stokes)

add_executable(domain_decomposition src_domain_decomposition/domain_decomposition_main.cpp src_domain_decomposition/Poisson2D.cpp)
deal_ii_setup_target(domain_decomposition)

add_executable(linear_elasticity src_linear_elasticity/linear_elasticity_main.cpp src_linear_elasticity/LinearElasticity.cpp)
deal_ii_setup_target(linear_elasticity)

add_executable(FK src_fk/main.cpp src_fk/FK_solver.cpp)
deal_ii_setup_target(FK)

add_executable(nonlinear_td src_nonlinear_td/NonLinear_main.cpp src_nonlinear_td/HeatNonLinear.cpp)
deal_ii_setup_target(nonlinear_td)

add_executable(1D_Mesh src_1D_Mesh/mesh_main.cpp)
deal_ii_setup_target(1D_Mesh)

add_executable(2D_Mesh src_2D_Mesh/mesh_main.cpp)
deal_ii_setup_target(2D_Mesh)

add_executable(nonlinear_steady src_nonlinear_steady/NonLinear_main.cpp src_nonlinear_steady/NonLinearDiffusion.cpp)
deal_ii_setup_target(nonlinear_steady)

add_executable(domain_decomposition_1D src_domain_decomposition_1D/domain_decomposition_main.cpp src_domain_decomposition_1D/Poisson2D.cpp)
deal_ii_setup_target(domain_decomposition_1D)


add_executable(stokes_td src_stokes_td/stokes_td_main.cpp src_stokes_td/stokes-time.cpp)
deal_ii_setup_target(stokes_td)

add_executable(stokes_si_td src_stokes_si_td/stokes_td_main.cpp src_stokes_si_td/stokes-time.cpp)
deal_ii_setup_target(stokes_si_td)  


add_executable(parabolic_1D src_parabolic_1D/parabolic_main.cpp src_parabolic_1D/Parabolic.cpp)
deal_ii_setup_target(parabolic_1D)






