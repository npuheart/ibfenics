# . ~/spack/share/spack/setup-env.sh
# spack env deactivate
spack load gcc@11.4.0
spack env activate -p kokkos
spack load cmake
spack load py-pip