include ../../../Makefile.${ARCH}

# Names of files that should be compiled and 
# inserted to libgridbuilder.a:

OBJS = sep_init.o\
	sep_base_class_shock.o sep_base_class_particle_splitter.o\
	sep_distrib_energy_container.o sep_distrib_energy_mono.o sep_distrib_energy_kappa.o sep_distrib_energy_powerlaw.o\
	sep_distrib_pitch_container.o sep_distrib_pitch_mono.o sep_distrib_pitch_isotropic.o\
	sep_distrib_wave_energy_base_class.o sep_distrib_wave_energy_powerlaw.o\
	sep_fields_cartesian_double_shock.o sep_fields_cartesian_homogeneous.o\
	sep_fields_cylindrical_linecurrent.o\
	sep_fields_spherical_mann.o sep_fields_spherical_parker.o\
	sep_fields_container.o sep_injectors_common.o\
	sep_lagr_species.o\
	sep_mesh_logical.o\
	sep_object_wrapper.o\
	sep_operator_accumulation_base.o sep_operator_anisotropy.o\
	sep_operator_density.o sep_operator_differential_flux.o sep_operator_Dmumu.o sep_operator_energy_channels.o sep_operator_mfp.o sep_operator_particle_energy.o\
	sep_operator_plasma_state.o sep_operator_reduced_scalar.o\
	sep_operator_spatial_lineout.o sep_operator_spatial_lineout_energy_spectrum.o\
	sep_operator_shock_electrons.o sep_operator_shock_mesh.o sep_operator_shock_variables.o sep_operator_spatial_slice.o sep_operator_spherical_surface_mesh.o\
	sep_operator_timeseries_energy.o sep_operator_volume.o\
	sep_operator_wave_cross_helicity.o sep_operator_wave_energy_total.o\
	sep_operator_wavelength_mesh.o sep_operator_wavelength_mesh_density.o\
	sep_operator_wavelength_mesh_particles.o sep_operator_wavelength_mesh_scalar.o sep_operator_wavelength_mesh_volume.o\
	sep_particle_scatterer.o sep_particle_species.o\
	sep_propagate.o\
	sep_register_bconds.o sep_register_dataoperators.o sep_register_distributions.o sep_register_lagr_injectors.o sep_register_lagr_propagators.o\
	sep_register_particle_injectors.o sep_register_particle_propagators.o\
	sep_simcontrol.o sep_shock_accelerator.o sep_shock_drift_acceleration.o sep_shock_paraboloid.o sep_shock_planar.o sep_shock_spherical.o\
	sep_user.o\
	sep_wavelength_mesh_builder.o

OBJS_VLSV_EXTRACT=vlsv_extract.o

# Make targets

default:
	${MAKE} lib${SIM}.a
	${MAKE} vlsv_extract_bin

clean:
	rm -rf *.o *.a *~
	rm -f ../lib/lib${SIM}.a

lib${SIM}.a: ${OBJS}
	${AR} r lib${SIM}.a ${OBJS}
	ln -f -s ${CURDIR}/lib${SIM}.a -t ../../lib

vlsv_extract_bin: ${OBJS_VLSV_EXTRACT}

# Dependencies

override INCS+=${INC} ${INC_BOOST} ${INC_PARGRID} ${INC_VLSV} ${INC_ZOLTAN}
override INCS+=-I../../include

override FLAGS+=-std=c++0x

INCS_SEP_INIT=-I../../gridbuilder -I../../dataoperator

DEPS_SEP_BASE_SHOCK=sep_base_class_shock.h sep_base_class_shock.cpp
DEPS_SEP_BASE_PARTICLE_SPLITTER=sep_base_class_particle_splitter.h sep_base_class_particle_splitter.cpp
DEPS_SEP_DISTRIB_ENERGY_CONTAINER=sep_distrib_energy_container.h sep_distrib_energy_container.cpp
DEPS_SEP_DISTRIB_ENERGY_KAPPA=${DEPS_SEP_DISTRIB_ENERGY_CONTAINER} sep_distrib_energy_kappa.h sep_distrib_energy_kappa.cpp
DEPS_SEP_DISTRIB_ENERGY_MONO=${DEPS_SEP_DISTRIB_ENERGY_CONTAINER} sep_distrib_energy_mono.h sep_distrib_energy_mono.cpp
DEPS_SEP_DISTRIB_ENERGY_POWERLAW=${DEPS_SEP_DISTRIB_ENERGY_CONTAINER} sep_distrib_energy_powerlaw.h sep_distrib_energy_powerlaw.cpp
DEPS_SEP_DISTRIB_PITCH_CONTAINER=sep_distrib_pitch_container.h sep_distrib_pitch_container.cpp
DEPS_SEP_DISTRIB_PITCH_ISOTROPIC=${DEPS_SEP_DISTRIB_PITCH_CONTAINER} sep_distrib_pitch_isotropic.h sep_distrib_pitch_isotropic.cpp
DEPS_SEP_DISTRIB_PITCH_MONO=sep_distrib_pitch_mono.h sep_distrib_pitch_mono.cpp
DEPS_SEP_DISTRIB_WAVE_BASE=sep_distrib_wave_energy_base_class.h sep_distrib_wave_energy_base_class.cpp
DEPS_SEP_DISTRIB_WAVE_POWERLAW=${DEPS_SEP_DISTRIB_WAVE_BASE} sep_distrib_wave_energy_powerlaw.h sep_distrib_wave_energy_powerlaw.cpp
DEPS_SEP_FIELDS_CARTESIAN_DOUBLE_SHOCK=${DEPS_SEP_FIELDS_CONTAINER} sep_distrib_wave_energy_base_class.h sep_fields_cartesian_double_shock.h sep_fields_cartesian_double_shock.cpp
DEPS_SEP_FIELDS_CARTESIAN_HOMOG=${DEPS_SEP_FIELDS_CONTAINER} sep_fields_cartesian_homogeneous.h sep_fields_cartesian_homogeneous.cpp
DEPS_SEP_FIELDS_CYLINDRICAL_LINECURRENT=${DEPS_SEP_FIELDS_CONTAINER} sep_fields_cylindrical_linecurrent.h sep_fields_cylindrical_linecurrent.cpp
DEPS_SEP_FIELDS_SPHERICAL_PARKER=${DEPS_SEP_FIELDS_CONTAINER} sep_fields_spherical_parker.h sep_fields_spherical_parker.cpp
DEPS_SEP_FIELDS_SPHERICAL_MANN=${DEPS_SEP_FIELDS_CONTAINER} sep_fields_spherical_mann.h sep_fields_spherical_mann.cpp
DEPS_SEP_FIELDS_CONTAINER=sep_fields_container.h sep_fields_container.cpp
DEPS_SEP_INIT=../../gridbuilder/logically_cartesian_builder.h sep_mesh_logical.h sep_particle_list.h sep_init.h sep_init.cpp
DEPS_SEP_INJ_COMMON=sep_namespaces.h sep_injectors_common.h sep_injectors_common.cpp
DEPS_SEP_LAGR_ACCUMULATOR=sep_accumulation_stretched.h sep_lagr_accumulator.h sep_wavelength_mesh_builder.h
DEPS_SEP_LAGR_BCOND=../../include/particle_list_util.h sep_lagr_boundary_cond.h
DEPS_SEP_LAGR_DEFINITION=sep_lagr_definition.h
DEPS_SEP_LAGR_INJ_DEFAULT=sep_distrib_wave_energy_base_class.h sep_lagr_injector_default.h
DEPS_SEP_LAGR_INJ_NG1994=sep_lagr_injector_ng_reames_1994.h sep_propagate.h
DEPS_SEP_LAGR_PROPAG_RK2=sep_lagr_accelerator.h sep_lagr_propagator_rk2.h
DEPS_SEP_LAGR_PROPAG_SHOCK_RK2=sep_injection_buffer.h sep_lagr_accelerator.h sep_lagr_propagator_shock_rk2.h
DEPS_SEP_LAGR_SPECIES=sep_lagr_species.h sep_lagr_species.cpp
DEPS_SEP_MESH_LOGICAL=sep_mesh_logical.h sep_mesh_logical.cpp
DEPS_SEP_OWRAPPER=sep_object_wrapper.h sep_object_wrapper.cpp
DEPS_SEP_OP_ACCUM_BASE=sep_operator_accumulation_base.h sep_operator_accumulation_base.cpp
DEPS_SEP_OP_ANISOTROPY=sep_fields_container.h sep_accumulation_stretched.h sep_operator_anisotropy.h sep_operator_anisotropy.cpp
DEPS_SEP_OP_DENSITY=sep_accumulation_stretched.h sep_operator_density.h sep_operator_density.cpp
DEPS_SEP_OP_DIFF_FLUX=${DEPS_SEP_OP_ACCUM_BASE} sep_operator_differential_flux.h sep_operator_differential_flux.cpp
DEPS_SEP_OP_DMUMU=sep_operator_differential_flux.h sep_particle_scatterer.h sep_operator_spatial_slice.h sep_operator_Dmumu.h sep_operator_Dmumu.cpp
DEPS_SEP_OP_ECHANNELS=${DEPS_SEP_OP_ACCUM_BASE} sep_fields_container.h sep_guiding_center_theory.h sep_operator_energy_channels.h sep_operator_energy_channels.cpp
DEPS_SEP_OP_MFP=${DEPS_SEP_OP_DIFF_FLUX} sep_operator_mfp.h sep_operator_mfp.cpp
DEPS_SEP_OP_PARTICLE_ENERGY=${DEPS_SEP_OP_ACCUM_BASE} sep_fields_container.h sep_guiding_center_theory.h sep_operator_particle_energy.h sep_operator_particle_energy.cpp
DEPS_SEP_OP_PLASMA_STATE=sep_coordinate_transform.h sep_operator_plasma_state.h sep_operator_plasma_state.cpp
DEPS_SEP_OP_REDUCED_SCALAR=sep_operator_reduced_scalar.h sep_operator_reduced_scalar.cpp
DEPS_SEP_OP_SHOCK_ELECTRONS=sep_distrib_energy_container.h sep_coordinate_transform.h sep_operator_shock_electrons.h sep_operator_shock_electrons.cpp
DEPS_SEP_OP_SHOCK_MESH=sep_coordinate_transform.h sep_operator_shock_mesh.h sep_operator_shock_mesh.cpp
DEPS_SEP_OP_SHOCK_VARS=${DEPS_SEP_OP_SHOCK_MESH} sep_operator_shock_variables.h sep_operator_shock_variables.cpp
DEPS_SEP_OP_SL_ENERGY_SPECTRUM=${DEPS_SEP_OP_SPATIAL_LINEOUT} sep_operator_spatial_lineout_energy_spectrum.h sep_operator_spatial_lineout_energy_spectrum.cpp
DEPS_SEP_OP_SPATIAL_LINEOUT=sep_mesh_logical.h sep_operator_spatial_lineout.h sep_operator_spatial_lineout.cpp
DEPS_SEP_OP_SPATIAL_SLICE=sep_operator_spatial_slice.h sep_operator_spatial_slice.cpp
DEPS_SEP_OP_SPHERICAL_SURF_MESH=sep_operator_spherical_surface_mesh.h sep_operator_spherical_surface_mesh.cpp
DEPS_SEP_OP_TIMESERIES_ENERGY=sep_operator_timeseries_energy.h sep_operator_timeseries_energy.cpp
DEPS_SEP_OP_WAVE_HELICITY=sep_operator_wave_cross_helicity.h sep_operator_wave_cross_helicity.cpp
DEPS_SEP_OP_WAVE_ENERGY=sep_propagate.h sep_operator_wave_energy_total.h sep_operator_wave_energy_total.cpp
DEPS_SEP_OP_WLM=sep_operator_wavelength_mesh.h sep_operator_wavelength_mesh.cpp
DEPS_SEP_OP_WLM_DENSITY=sep_operator_wavelength_mesh_density.h sep_operator_wavelength_mesh_density.cpp
DEPS_SEP_OP_WLM_PARTICLES=${DEPS_SEP_LAGR_SPECIES} ${DEPS_SEP_LAGR_DEFINITION} ${DEPS_SEP_PARTICLE_SPECIES} ${DEPS_SEP_PARTICLE_DEFINITION}
DEPS_SEP_OP_WLM_PARTICLES+=sep_fields_container.h sep_operator_wavelength_mesh_particles.h sep_operator_wavelength_mesh_particles.cpp
DEPS_SEP_OP_WLM_SCALAR=sep_operator_wavelength_mesh_scalar.h sep_operator_wavelength_mesh_scalar.cpp
DEPS_SEP_OP_WLM_VOLUME=sep_propagate.h sep_operator_wavelength_mesh_volume.h sep_operator_wavelength_mesh_volume.cpp
DEPS_SEP_OP_VOLUME=sep_propagate.h sep_operator_volume.h sep_operator_volume.cpp
DEPS_SEP_PARTICLE_ACCUMULATOR=sep_particle_accumulator.h sep_accumulation_stretched.h
DEPS_SEP_PARTICLE_BCOND=../../include/particle_list_util.h sep_particle_boundary_cond.h sep_particle_boundary_cond_none.h
DEPS_SEP_PARTICLE_BCOND_DOUBLE_SHOCK=sep_particle_boundary_cond_double_shock.h
DEPS_SEP_PARTICLE_BCOND_SHOCK=../../include/particle_list_util.h sep_particle_boundary_cond_shock.h
DEPS_SEP_PARTICLE_DEFINITION=sep_particle_definition.h
DEPS_SEP_PARTICLE_INJ_CELL=sep_particle_injector_cell.h sep_injectors_common.h
DEPS_SEP_PARTICLE_INJ_HOMOG=sep_particle_injector_homogeneous.h sep_injectors_common.h
DEPS_SEP_PARTICLE_INJ_NG1994=sep_particle_injector_ng_reames_1994.h sep_propagate.h sep_injectors_common.h
DEPS_SEP_PARTICLE_INJ_SHOCK=sep_particle_injector_shock.h sep_propagate.h sep_injectors_common.h
DEPS_SEP_PARTICLE_PROPAG_RK2=${DEPS_SEP_SHOCK_DSA} sep_particle_propagator_rk2.h sep_particle_propagator_shock_rk2.h sep_particle_propagator_coronal_rk2.h sep_particle_accelerator.h
DEPS_SEP_PARTICLE_SCATTERER=sep_fields_container.h sep_particle_scatterer.h sep_particle_scatterer.cpp sep_accumulation_stretched.h
DEPS_SEP_PARTICLE_SPECIES=sep_particle_species.h sep_particle_species.cpp
DEPS_SEP_PARTICLE_SPLITTER=${DEPS_SEP_BASE_PARTICLE_SPLITTER} sep_particle_splitter.h
DEPS_SEP_PROPAGATE=sep_injection_buffer.h sep_debugging.h sep_particle_list.h sep_propagate.h sep_propagate.cpp
DEPS_SEP_REG_BCONDS=sep_register_bconds.h sep_register_bconds.cpp
DEPS_SEP_REG_DATAOPS=sep_register_dataoperators.h sep_register_dataoperators.cpp
DEPS_SEP_REG_DISTRIBS=sep_register_distributions.h sep_register_distributions.cpp
DEPS_SEP_REG_LAGR_INJ=sep_register_lagr_injectors.h sep_register_lagr_injectors.cpp
DEPS_SEP_REG_LAGR_PROPAG=sep_register_lagr_propagators.h sep_register_lagr_propagators.cpp
DEPS_SEP_REG_PINJ=sep_register_particle_injectors.h sep_register_particle_injectors.cpp
DEPS_SEP_REG_PPROPAG=sep_register_particle_propagators.h sep_register_particle_propagators.cpp
DEPS_SEP_SHOCK_ACC=sep_shock_accelerator.h sep_shock_accelerator.cpp
DEPS_SEP_SHOCK_DSA=sep_shock_drift_acceleration.h sep_shock_drift_acceleration.cpp
DEPS_SEP_SHOCK_PARABOLOID=sep_shock_paraboloid.h sep_shock_paraboloid.cpp
DEPS_SEP_SHOCK_PLANAR=${DEPS_SEP_BASE_SHOCK} ../../include/triangulated_plane.h sep_shock_planar.h sep_shock_planar.cpp
DEPS_SEP_SHOCK_SPHERICAL=${DEPS_SEP_BASE_SHOCK} ../../include/triangulated_sphere.h sep_shock_spherical.h sep_shock_spherical.cpp
DEPS_SEP_SIMCONTROL=../../include/main.h sep_fields_container.h sep_simcontrol.h sep_simcontrol.cpp
DEPS_SEP_USER=../../include/user.h sep_init.h sep_propagate.h sep_user.cpp
DEPS_SEP_WLM_BUILDER=sep_wavelength_mesh_builder.h sep_wavelength_mesh_builder.cpp

DEPS_ALL=${DEPS_SEP_SIMCONTROL} ${DEPS_SEP_OWRAPPER}

DEPS_SEP_DATA_OPS_COMMON=../../include/dataoperator.h ../../dataoperator/dataoperator.cpp

DEPS_SEP_DATA_OPS=${DEPS_SEP_DATA_OPS_COMMON}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_ACCUM_BASE}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_ANISOTROPY}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_DENSITY}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_DIFF_FLUX}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_DMUMU}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_ECHANNELS}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_MFP}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_PARTICLE_ENERGY}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_PLASMA_STATE}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_REDUCED_SCALAR}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_SHOCK_ELECTRONS}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_SHOCK_MESH}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_SHOCK_VARS}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_SPATIAL_LINEOUT}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_SPATIAL_SLICE}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_SL_ENERGY_SPECTRUM}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_SPHERICAL_SURF_MESH}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_TIMESERIES_ENERGY}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_WAVE_HELICITY}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_WAVE_ENERGY}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_WLM}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_WLM_DENSITY}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_WLM_PARTICLES}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_WLM_SCALAR}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_WLM_VOLUME}
DEPS_SEP_DATA_OPS+=${DEPS_SEP_OP_VOLUME}

DEPS_SEP_DISTRIBS_ENERGY=
DEPS_SEP_DISTRIBS_ENERGY+=${DEPS_SEP_DISTRIB_ENERGY_CONTAINER}
DEPS_SEP_DISTRIBS_ENERGY+=${DEPS_SEP_DISTRIB_ENERGY_MONO}
DEPS_SEP_DISTRIBS_ENERGY+=${DEPS_SEP_DISTRIB_ENERGY_POWERLAW}
DEPS_SEP_DISTRIBS_ENERGY+=${DEPS_SEP_DISTRIB_ENERGY_KAPPA}

DEPS_SEP_DISTRIBS_PITCH=
DEPS_SEP_DISTRIBS_PITCH+=${DEPS_SEP_DISTRIB_PITCH_CONTAINER}
DEPS_SEP_DISTRIBS_PITCH+=${DEPS_SEP_DISTRIB_PITCH_ISOTROPIC}
DEPS_SEP_DISTRIBS_PITCH+=${DEPS_SEP_DISTRIB_PITCH_MONO}

DEPS_SEP_DISTRIBS_WAVES=
DEPS_SEP_DISTRIBS_WAVES+=${DEPS_SEP_DISTRIB_WAVE_BASE}
DEPS_SEP_DISTRIBS_WAVES+=${DEPS_SEP_DISTRIB_WAVE_POWERLAW}

DEPS_SEP_DISTRIBS=
DEPS_SEP_DISTRIBS+=${DEPS_SEP_DISTRIBS_ENERGY}
DEPS_SEP_DISTRIBS+=${DEPS_SEP_DISTRIBS_PITCH}
DEPS_SEP_DISTRIBS+=${DEPS_SEP_DISTRIBS_WAVES}

DEPS_SEP_FIELDS=
DEPS_SEP_FIELDS+=${DEPS_SEP_FIELDS_CONTAINER}
DEPS_SEP_FIELDS+=${DEPS_SEP_FIELDS_CARTESIAN_DOUBLE_SHOCK}
DEPS_SEP_FIELDS+=${DEPS_SEP_FIELDS_CARTESIAN_HOMOGENEOUS}
DEPS_SEP_FIELDS+=${DEPS_SEP_FIELDS_CYLINDRICAL_LINECURRENT}
DEPS_SEP_FIELDS+=${DEPS_SEP_FIELDS_SPHERICAL_PARKER}
DEPS_SEP_FIELDS+=${DEPS_SEP_FIELDS_SPHERICAL_MANN}

DEPS_SEP_PINJ+=${DEPS_SEP_PARTICLE_INJ_CELL}
DEPS_SEP_PINJ+=${DEPS_SEP_PARTICLE_INJ_HOMOG}
DEPS_SEP_PINJ+=${DEPS_SEP_PARTICLE_INJ_NG1994}
DEPS_SEP_PINJ+=${DEPS_SEP_PARTICLE_INJ_SHOCK}

DEPS_SEP_LAGR_INJ+=${DEPS_SEP_LAGR_INJ_DEFAULT}
DEPS_SEP_LAGR_INJ+=${DEPS_SEP_LAGR_INJ_NG1994}

DEPS_SEP_BCONDS=
DEPS_SEP_BCONDS+=${DEPS_SEP_PARTICLE_BCOND}
DEPS_SEP_BCONDS+=${DEPS_SEP_PARTICLE_BCOND_DOUBLE_SHOCK}
DEPS_SEP_BCONDS+=${DEPS_SEP_PARTICLE_BCOND_SHOCK}
DEPS_SEP_BCONDS+=${DEPS_SEP_LAGR_BCOND}

DEPS_SEP_PARTICLES=
DEPS_SEP_PARTICLES+=${DEPS_SEP_PARTICLE_ACCUMULATOR}
DEPS_SEP_PARTICLES+=${DEPS_SEP_PARTICLE_DEFINITION}
DEPS_SEP_PARTICLES+=${DEPS_SEP_PARTICLE_SPECIES}
DEPS_SEP_PARTICLES+=${DEPS_SEP_LAGR_ACCUMULATOR}
DEPS_SEP_PARTICLES+=${DEPS_SEP_LAGR_DEFINITION}
DEPS_SEP_PARTICLES+=${DEPS_SEP_LAGR_SPECIES}

DEPS_SEP_REG_LAGR_PROPAG+=${DEPS_SEP_LAGR_PROPAG_RK2}
DEPS_SEP_REG_LAGR_PROPAG+=${DEPS_SEP_LAGR_PROPAG_SHOCK_RK2}

DEPS_SEP_REG_PPROPAG+=${DEPS_SEP_PARTICLE_PROPAG_RK2}

DEPS_SEP_INIT+=${DEPS_SEP_FIELDS}
DEPS_SEP_INIT+=${DEPS_SEP_PARTICLES}
DEPS_SEP_INIT+=${DEPS_SEP_PROPAGATORS}
DEPS_SEP_INIT+=${DEPS_SEP_REG_BCONDS}
DEPS_SEP_INIT+=${DEPS_SEP_REG_DATAOPS}
DEPS_SEP_INIT+=${DEPS_SEP_REG_DISTRIBS}
DEPS_SEP_INIT+=${DEPS_SEP_REG_PINJ}
DEPS_SEP_INIT+=${DEPS_SEP_REG_PPROPAG}

DEPS_SEP_OWRAPPER+=${DEPS_SEP_BASE_SHOCK}
DEPS_SEP_OWRAPPER+=${DEPS_SEP_DISTRIB_WAVE_BASE}
DEPS_SEP_OWRAPPER+=${DEPS_SEP_FIELDS_CONTAINER}
DEPS_SEP_OWRAPPER+=${DEPS_SEP_DISTRIB_ENERGY_CONTAINER}
DEPS_SEP_OWRAPPER+=${DEPS_SEP_DISTRIB_PITCH_CONTAINER}
DEPS_SEP_OWRAPPER+=${DEPS_SEP_BASE_PARTICLE_SPLITTER}

# Compilation rules

.cpp.o:
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INCS}

sep_base_class_shock.o: ${DEPS_SEP_BASE_SHOCK} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_base_class_shock.cpp ${INCS}

sep_base_class_particle_splitter.o: ${DEPS_SEP_BASE_PARTICLE_SPLITTER} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_base_class_particle_splitter.cpp ${INCS}

sep_distrib_energy_container.o: ${DEPS_SEP_DISTRIB_ENERGY_CONTAINER} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${FLAGS} -c sep_distrib_energy_container.cpp ${INCS}

sep_distrib_energy_kappa.o: ${DEPS_SEP_DISTRIB_ENERGY_KAPPA} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -std=c++0x -c sep_distrib_energy_kappa.cpp ${INCS}

sep_distrib_energy_mono.o: ${DEPS_SEP_DISTRIB_ENERGY_MONO} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_distrib_energy_mono.cpp ${INCS}

sep_distrib_energy_powerlaw.o: ${DEPS_SEP_DISTRIB_ENERGY_POWERLAW} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_distrib_energy_powerlaw.cpp ${INCS}

sep_distrib_pitch_container.o: ${DEPS_SEP_DISTRIB_PITCH_CONTAINER} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${FLAGS} -c sep_distrib_pitch_container.cpp ${INCS}

sep_distrib_pitch_mono.o: ${DEPS_SEP_DISTRIB_PITCH_MONO} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_distrib_pitch_mono.cpp ${INCS}

sep_distrib_pitch_isotropic.o: ${DEPS_SEP_DISTRIB_PITCH_ISOTROPIC} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_distrib_pitch_isotropic.cpp ${INCS}

sep_distrib_wave_energy_base.o: ${DEPS_SEP_DISTRIB_WAVE_BASE} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c sep_distrib_wave_energy_base.cpp ${INCS}

sep_distrib_wave_energy_powerlaw.o: ${DEPS_SEP_DISTRIB_WAVE_POWERLAW} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c sep_distrib_wave_energy_powerlaw.cpp ${INCS}

sep_fields_cartesian_double_shock.o: ${DEPS_SEP_FIELDS_CARTESIAN_DOUBLE_SHOCK} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_fields_cartesian_double_shock.cpp ${INCS}

sep_fields_cartesian_homogeneous.o: ${DEPS_SEP_FIELDS_CARTESIAN_HOMOG} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_fields_cartesian_homogeneous.cpp ${INCS}

sep_fields_cylindrical_linecurrent.o: ${DEPS_SEP_FIELDS_CYLINDRICAL_LINECURRENT} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_fields_cylindrical_linecurrent.cpp ${INCS}

sep_fields_spherical_parker.o: ${DEPS_SEP_FIELDS_SPHERICAL_PARKER} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_fields_spherical_parker.cpp ${INCS}

sep_fields_spherical_mann.o: ${DEPS_SEP_FIELDS_SPHERICAL_MANN} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_fields_spherical_mann.cpp ${INCS}

sep_fields_container.o: ${DEPS_SEP_FIELDS_CONTAINER} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${FLAGS} -c sep_fields_container.cpp ${INCS}

sep_init.o: ${DEPS_SEP_INIT} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_init.cpp ${INCS} ${INCS_SEP_INIT}

sep_injectors_common.o: ${DEPS_SEP_INJ_COMMON} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_injectors_common.cpp ${INCS}

sep_lagr_species.o: ${DEPS_LAGR_SPECIES} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_lagr_species.cpp ${INCS}

sep_mesh_logical.o: ${DEPS_SEP_MESH_LOGICAL} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_mesh_logical.cpp ${INCS}

sep_object_wrapper.o: ${DEPS_SEP_OWRAPPER} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_object_wrapper.cpp ${INCS}

sep_operator_accumulation_base.o: ${DEPS_SEP_OP_ACCUM_BASE} ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_accumulation_base.cpp ${INCS}

sep_operator_anisotropy.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_ACCUM_BASE} ${DEPS_SEP_OP_ANISOTROPY} ${DEPS_ALL} 
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_anisotropy.cpp ${INCS}

sep_operator_density.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_ACCUM_BASE} ${DEPS_SEP_OP_DENSITY} ${DEPS_ALL} 
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_density.cpp ${INCS}

sep_operator_differential_flux.o:${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_DIFF_FLUX} ${DEPS_ALL}
	 ${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_differential_flux.cpp ${INCS}

sep_operator_Dmumu.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_DMUMU} ${DEPS_ALL}
	 ${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_Dmumu.cpp ${INCS}	

sep_operator_energy_channels.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_ACCUM_BASE} ${DEPS_SEP_OP_ECHANNELS} ${DEPS_ALL} 
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_energy_channels.cpp ${INCS} 

sep_operator_mfp.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_MFP} ${DEPS_ALL} 
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_mfp.cpp ${INCS}

sep_operator_particle_energy.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_PARTICLE_ENERGY} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_particle_energy.cpp ${INCS}

sep_operator_plasma_state.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_PLASMA_STATE} ${DEPS_ALL} 
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_plasma_state.cpp ${INCS}

sep_operator_reduced_scalar.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_REDUCED_SCALAR} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_reduced_scalar.cpp ${INCS}

sep_operator_spatial_lineout.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SPATIAL_LINEOUT} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_spatial_lineout.cpp ${INCS}

sep_operator_spatial_slice.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SPATIAL_SLICE} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_spatial_slice.cpp ${INCS}

sep_operator_spatial_lineout_energy_spectrum.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SL_ENERGY_SPECTRUM} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_spatial_lineout_energy_spectrum.cpp ${INCS}

sep_operator_shock_electrons.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SHOCK_MESH} ${DEPS_SEP_OP_SHOCK_ELECTRONS} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_shock_electrons.cpp ${INCS}	

sep_operator_shock_mesh.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SHOCK_MESH} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_shock_mesh.cpp ${INCS}

sep_operator_shock_variables.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SHOCK_VARS} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_shock_variables.cpp ${INCS}

sep_operator_spherical_surface_mesh.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SPHERICAL_SURF_MESH} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_spherical_surface_mesh.cpp ${INCS}

sep_operator_timeseries_energy.o: ${DEPS_SEP_OP_TIMESERIES_ENERGY}  ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_timeseries_energy.cpp ${INCS}

sep_operator_volume.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_VOLUME} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_volume.cpp ${INCS}

sep_operator_wave_cross_helicity.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_WAVE_HELICITY} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_wave_cross_helicity.cpp ${INCS}

sep_operator_wave_energy_total.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_WAVE_ENERGY} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_wave_energy_total.cpp ${INCS}

sep_operator_wavelength_mesh.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SPATIAL_SLICE} ${DEPS_SEP_OP_WLM} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_wavelength_mesh.cpp ${INCS}

sep_operator_wavelength_mesh_density.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SPATIAL_SLICE} ${DEPS_SEP_OP_WLM_DENSITY} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_wavelength_mesh_density.cpp ${INCS}

sep_operator_wavelength_mesh_particles.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SPATIAL_SLICE} ${DEPS_SEP_OP_WLM_PARTICLES} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_wavelength_mesh_particles.cpp ${INCS}

sep_operator_wavelength_mesh_scalar.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SPATIAL_SLICE} ${DEPS_SEP_OP_WLM_SCALAR} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_wavelength_mesh_scalar.cpp ${INCS}

sep_operator_wavelength_mesh_volume.o: ${DEPS_SEP_DATA_OPS_COMMON} ${DEPS_SEP_OP_SPATIAL_SLICE} ${DEPS_SEP_OP_WLM_VOLUME} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_operator_wavelength_mesh_volume.cpp ${INCS}

sep_particle_scatterer.o: ${DEPS_SEP_PARTICLE_SCATTERER} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_particle_scatterer.cpp ${INCS}

sep_particle_species.o: ${DEPS_SEP_PARTICLE_SPECIES} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_particle_species.cpp ${INCS}

sep_propagate.o: ${DEPS_SEP_PROPAGATE} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_propagate.cpp ${INCS}

sep_register_bconds.o: ${DEPS_SEP_REG_BCONDS} ${DEPS_SEP_BCONDS} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_register_bconds.cpp ${INCS}

sep_register_dataoperators.o: ${DEPS_SEP_DATA_OPS} ${DEPS_SEP_REG_DATAOPS} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_register_dataoperators.cpp ${INCS} -I../../dataoperator

sep_register_distributions.o: ${DEPS_SEP_DISTRIBS} ${DEPS_SEP_REG_DISTRIBS} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_register_distributions.cpp ${INCS}

sep_register_lagr_injectors.o: ${DEPS_SEP_REG_LAGR_INJ} ${DEPS_SEP_LAGR_INJ} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_register_lagr_injectors.cpp ${INCS}

sep_register_lagr_propagators.o: ${DEPS_SEP_REG_LAGR_PROPAG} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_register_lagr_propagators.cpp ${INCS}

sep_register_particle_injectors.o: ${DEPS_SEP_REG_PINJ} ${DEPS_SEP_PINJ} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_register_particle_injectors.cpp ${INCS}

sep_register_particle_propagators.o: ${DEPS_SEP_REG_PPROPAG} ${DEPS_ALL} ${DEPS_SEP_PARTICLE_SPLITTER}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_register_particle_propagators.cpp ${INCS}

sep_shock_accelerator.o: ${DEPS_SEP_SHOCK_ACC} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_shock_accelerator.cpp ${INCS}

sep_shock_drift_acceleration.o: ${DEPS_SEP_SHOCK_DSA} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_shock_drift_acceleration.cpp ${INCS}

sep_shock_paraboloid.o: ${DEPS_SEP_SHOCK_PARABOLOID} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_shock_paraboloid.cpp ${INCS}

sep_shock_planar.o: ${DEPS_SEP_SHOCK_PLANAR} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_shock_planar.cpp ${INCS}

sep_shock_spherical.o: ${DEPS_SEP_SHOCK_SPHERICAL} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_shock_spherical.cpp ${INCS}

sep_simcontrol.o: ${DEPS_SEP_SIMCONTROL} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${FLAGS} -c sep_simcontrol.cpp ${INCS}

sep_user.o: ${DEPS_SEP_USER} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c sep_user.cpp ${INCS}

sep_wavelength_mesh_builder.o: ${DEPS_SEP_WLM_BUILDER} ${DEPS_ALL}
	${CMP} ${CXXFLAGS} ${FLAGS} -c sep_wavelength_mesh_builder.cpp ${INCS}

vlsv_extract.o: vlsv_extract.cpp
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c vlsv_extract.cpp ${INCS}
	${CMP} -o vlsv_extract vlsv_extract.o ${LIB_VLSV}
