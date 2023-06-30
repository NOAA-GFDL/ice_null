!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Ice Null Model Component.
!*
!* Ice Null is free software: you can redistribute it and/or modify it
!* under the terms of the GNU Lesser General Public License as published
!* by the Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* Ice Null is distributed in the hope that it will be useful, but
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with Ice Null.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module ice_model_mod

use   ice_albedo_mod, only:  ice_albedo_init, ice_albedo

use ocean_albedo_mod, only:  compute_ocean_albedo_new

use  ocean_rough_mod, only:  compute_ocean_roughness, fixed_ocean_roughness

use fms2_io_mod, only : file_exists, open_file, close_file, write_restart, read_data, &
                        variable_exists, FmsNetcdfDomainFile_t,FmsNetcdfFile_t, read_new_restart 

use          fms_mod, only: mpp_pe, mpp_root_pe, mpp_npes, write_version_number, stdlog,   &
                            error_mesg, FATAL, check_nml_error, &
                            NOTE, WARNING, stdout, &
                            clock_flag_default


use mpp_mod,          only: mpp_chksum, mpp_clock_id, CLOCK_COMPONENT, &
                            CLOCK_LOOP, CLOCK_ROUTINE, mpp_clock_begin, mpp_clock_end

#ifdef INTERNAL_FILE_NML
use          mpp_mod, only: input_nml_file
#else
use          fms_mod, only: open_namelist_file
#endif

use    constants_mod, only: hlv, hlf, tfreeze, pi, radius

use  mpp_domains_mod, only: domain1d, domain2d, mpp_define_domains, mpp_get_compute_domain, &
                            mpp_get_compute_domains, mpp_get_domain_components, mpp_get_pelist, &
                            mpp_define_layout, mpp_define_io_domain, mpp_global_field, &
                            XUPDATE, YUPDATE, CORNER, mpp_set_domain_symmetry, BGRID_NE

use diag_manager_mod, only: diag_axis_init, register_diag_field, send_data

use time_manager_mod, only: time_type, operator(+)

use mosaic2_mod,       only: get_mosaic_ntiles, get_mosaic_grid_sizes
use grid2_mod,         only: get_grid_comp_area, get_grid_size, get_grid_cell_vertices, &
                            get_grid_cell_centers, get_grid_cell_area


!use  amip_interp_mod, only: amip_interp_type, amip_interp_new
use coupler_types_mod,only: coupler_1d_bc_type, coupler_2d_bc_type, coupler_3d_bc_type
implicit none
private

public :: ice_data_type, ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
public :: ice_model_init, ice_model_end, update_ice_model_fast, ice_stock_pe, cell_area
public :: ice_model_restart
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum, ice_data_type_chksum
public :: update_ice_atm_deposition_flux, share_ice_domains
public :: unpack_ocean_ice_boundary, exchange_slow_to_fast_ice, set_ice_surface_fields
public :: ice_model_fast_cleanup, unpack_land_ice_boundary
public :: exchange_fast_to_slow_ice, update_ice_model_slow

type ice_data_type
  type(domain2d)                        :: Domain
  type(domain2d)                        :: slow_Domain_NH

   logical                              :: pe
   logical                              :: slow_ice_pe
   logical                              :: fast_ice_pe
   logical                              :: shared_slow_fast_PEs = .true.

   integer, pointer, dimension(:)       :: pelist =>NULL()
   integer, pointer, dimension(:)       :: slow_pelist =>NULL()
   integer, pointer, dimension(:)       :: fast_pelist =>NULL()

   real,    pointer, dimension(:)       :: glon_bnd =>NULL(), &
                                           glat_bnd =>NULL(), &
                                           lon_bnd =>NULL() , &
                                           lat_bnd =>NULL()

   real,    pointer, dimension(:,:)     :: lon =>NULL(), &
                                           lat =>NULL(), &
                                           SST_C =>NULL()

   logical, pointer, dimension(:,:)     :: mask =>NULL()

   logical, pointer, dimension(:,:,:)   :: ice_mask =>NULL()

   real,    pointer, dimension(:,:,:,:) :: temp =>NULL()

   real,    pointer, dimension(:,:,:)   :: part_size =>NULL(), &
                                           t_surf =>NULL(), &
                                           albedo =>NULL(), &
                                           albedo_vis_dir =>NULL(), &
                                           albedo_nir_dir =>NULL(), &
                                           albedo_vis_dif =>NULL(), &
                                           albedo_nir_dif =>NULL(), &
                                           rough_mom =>NULL(),&
                                           rough_heat =>NULL(), &
                                           rough_moist =>NULL(),  &
                                           frazil =>NULL(),  &
                                           u_surf =>NULL(),  &
                                           v_surf =>NULL()

   real,    pointer, dimension(:,:,:)   :: flux_u_bot =>NULL(), &
                                           flux_v_bot =>NULL(), &
                                           flux_t_bot =>NULL(),   &
                                           flux_q_bot =>NULL(), &
                                           flux_lh_bot =>NULL(), &
                                           flux_sw_bot =>NULL(), &
                                           flux_sw_vis_bot =>NULL(), &
                                           flux_sw_dir_bot =>NULL(), &
                                           flux_sw_dif_bot =>NULL(), &
                                           flux_sw_vis_dir_bot =>NULL(), &
                                           flux_sw_vis_dif_bot =>NULL(), &
                                           flux_sw_nir_dir_bot =>NULL(), &
                                           flux_sw_nir_dif_bot =>NULL(), &
                                           flux_lw_bot =>NULL(), &
                                           lprec_bot =>NULL(), &
                                           fprec_bot =>NULL(), &
                                           runoff_bot =>NULL()

   real,    pointer, dimension(:,:  )   :: flux_u =>NULL(), &
                                           flux_v =>NULL(), &
                                           flux_t =>NULL(), &
                                           flux_q =>NULL(), &
                                           flux_lh =>NULL(), &
                                           flux_sw =>NULL(), &
                                           flux_sw_vis =>NULL(), &
                                           flux_sw_dir =>NULL(), &
                                           flux_sw_dif =>NULL(), &
                                           flux_sw_vis_dir =>NULL(), &
                                           flux_sw_vis_dif =>NULL(), &
                                           flux_sw_nir_dir =>NULL(), &
                                           flux_sw_nir_dif =>NULL(), &
                                           flux_lw =>NULL(), &
                                           lprec =>NULL(), &
                                           fprec =>NULL(), &
                                           p_surf =>NULL(), &
                                           runoff =>NULL(), &
                                           calving =>NULL(), &
                                           stress_mag =>NULL(), &
                                           ustar_berg =>NULL(), &
                                           area_berg =>NULL(), &
                                           mass_berg =>NULL(), &
                                           runoff_hflx =>NULL(), &
                                           calving_hflx =>NULL(), &
                                           area =>NULL(), &
                                           flux_salt =>NULL()
  logical, pointer, dimension(:,:) :: maskmap =>NULL()   ! A pointer to an array indicating which
                                                         ! logical processors are actually used for
                                                         ! the ocean code. The other logical
                                                         ! processors would be all land points and
                                                         ! are not assigned to actual processors.
                                                         ! This need not be assigned if all logical
                                                         ! processors are used. This variable is dummy and need
                                                         ! not to be set, but it is needed to pass compilation.

   integer                              :: avg_kount

   real,    pointer, dimension(:,:,:)   :: thickness =>NULL()

   type (time_type)                     :: Time_Init, Time,  &
                                           Time_step_fast,   &
                                           Time_step_slow
   integer, dimension(3)              :: axes
   type(coupler_3d_bc_type)           :: ocean_fields       ! array of fields used for additional tracers
   type(coupler_2d_bc_type)           :: ocean_fluxes       ! array of fluxes used for additional tracers
   type(coupler_3d_bc_type)           :: ocean_fluxes_top   ! array of fluxes for averaging
   real,    pointer, dimension(:,:)   :: mi                  =>NULL() ! This is needed for the wave model. It is introduced here,
   integer :: flux_uv_stagger = -999

end type ice_data_type

type :: ocean_ice_boundary_type
  real, dimension(:,:),   pointer :: u =>NULL(), &
                                     v =>NULL(), &
                                     t =>NULL(), &
                                     s =>NULL(), &
                                     frazil =>NULL(), &
                                     sea_level =>NULL()
  real, dimension(:,:,:), pointer :: data =>NULL()
  integer                         :: xtype
  type(coupler_2d_bc_type)        :: fields     ! array of fields used for additional tracers
  integer                         :: stagger = BGRID_NE
end type ocean_ice_boundary_type

type :: atmos_ice_boundary_type
  real, dimension(:,:,:), pointer :: u_flux =>NULL(), &
                                     v_flux =>NULL(), &
                                     u_star =>NULL(), &
                                     t_flux =>NULL(), &
                                     q_flux =>NULL(), &
                                     lw_flux =>NULL(), &
                                     sw_flux_vis_dir =>NULL(), &
                                     sw_flux_vis_dif =>NULL(), &
                                     sw_flux_nir_dir =>NULL(), &
                                     sw_flux_nir_dif =>NULL(), &
                                     sw_down_vis_dir =>NULL(), &
                                     sw_down_vis_dif =>NULL(), &
                                     sw_down_nir_dir =>NULL(), &
                                     sw_down_nir_dif =>NULL(), &
                                     lprec =>NULL(), &
                                     fprec =>NULL()
  real, dimension(:,:,:), pointer :: dhdt =>NULL(), &
                                     dedt =>NULL(), &
                                     drdt =>NULL(), &
                                     coszen =>NULL(), &
                                     p =>NULL(), &
                                     data =>NULL()
  integer                         :: xtype
  type(coupler_3d_bc_type)        :: fluxes     ! array of fluxes used for additional tracers
end type atmos_ice_boundary_type

type :: land_ice_boundary_type
  real, dimension(:,:),   pointer :: runoff =>NULL(), &
                                     calving =>NULL(), &
                                     runoff_hflx =>NULL(), &
                                     calving_hflx =>NULL()
  real, dimension(:,:,:), pointer :: data =>NULL()
  integer :: xtype
end type land_ice_boundary_type

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

character(len=80) :: restart_format = 'amip ice model restart format 02'
logical :: module_is_initialized = .false.
logical :: stock_warning_issued  = .false.


integer                           :: isc, iec, jsc, jec ! compute domain

! iceClock variables from ice_type.F90 for compatibility with ice_sis interfaces
integer :: iceClock, iceClock1, iceCLock3

! id's for diagnostics
integer :: id_sh, id_lh, id_sw, id_lw, id_snofl, id_rain, &
           id_runoff, id_calving, id_evap, id_fax, id_fay, &
           id_sw_vis, id_sw_dir, id_sw_dif, id_sw_vis_dir, id_sw_vis_dif, &
           id_sw_nir_dir, id_sw_nir_dif
logical :: sent

!-----------------------------------------------------------------------
!
!  use_climo_ice           = use monthly climatological amip ice mask
!  use_annual_ice          = use annual climatology amip ice mask
!  temp_ice_freeze         = temperature at which sea ice melts
!  no_ice                  = run with no ice (only open water)
!  use_leads               = use fraction ice coverage (i.e., leads) if it exists
!  roughness_ice           = constant roughness for all ice
!  specified_ice_thickness = constant thickness for specified ice

integer, parameter :: num_lev  = 1
integer, parameter :: num_part = 2

real    :: diff                     = 2.092
real    :: thickness_min            = 0.10
real    :: specified_ice_thickness  = 2.0
real    :: temp_ice_freeze          = -1.66    ! was 271.5
real    :: roughness_ice            = 1.e-4
logical :: no_ice                   = .false.
logical :: use_leads                = .false.
logical :: use_climo_ice            = .false.
logical :: use_annual_ice           = .false.
integer, dimension(2) :: layout     = (/ 0, 0 /)
integer, dimension(2) :: io_layout  = (/ 1, 1 /)
character(len=64) :: interp_method  = "conservative" ! default conservative scheme
logical :: do_netcdf_restart        = .true.
character(len=128) :: axisname_x    = 'xt'  ! x-axis name of temperature grid
character(len=128) :: axisname_y    = 'yt'  ! y-axis name of temperature grid
character(len=128) :: axisname_xb   = 'xb'  ! x-axis bounds name of temperature grid
character(len=128) :: axisname_yb   = 'yb'  ! y-axis bounds name of temperature grid

namelist /ice_model_nml/ do_netcdf_restart, diff, thickness_min, &
                         specified_ice_thickness,                &
                         temp_ice_freeze, roughness_ice,         &
                         use_climo_ice, use_annual_ice,          &
                         no_ice, use_leads, layout, interp_method,  &
                         axisname_x, axisname_y, axisname_xb, axisname_yb, &
                         io_layout

real, parameter :: latent = hlv + hlf
!type(amip_interp_type), save :: Amip
real, allocatable, dimension(:,:) ::  cell_area  ! grid cell area; sphere frac.

integer :: id_restart_albedo
integer :: mlon, mlat, npart ! global grid size
type(FmsNetcdfDomainFile_t), save :: Ice_restart
character(len=8) :: dim_names(4)

! interface for fast ice for compatibility with SIS2
interface update_ice_model_fast ! overload to support old interface
     module procedure update_ice_model_fast_new
end interface

contains

!=============================================================================================
  subroutine ice_model_init( Ice, Time_Init, Time, Time_step_fast, Time_step_slow, Verona_coupler, &
                             concurrent_ice, gas_fluxes, gas_fields_ocn)
    type(ice_data_type), intent(inout) :: Ice
    type(time_type)    , intent(in)    :: Time_Init, Time, Time_step_fast, Time_step_slow
    logical,   optional, intent(in)    :: Verona_coupler
    logical,    optional, intent(in)    :: Concurrent_ice ! for compatibility with SIS2. For SIS1
                                                          ! there is no difference between fast and
                                                          ! slow ice PEs.
    type(coupler_1d_bc_type), &
                optional, intent(in)    :: gas_fluxes     ! If present, this type describes the
                                              ! additional gas or other tracer fluxes between the
                                              ! ocean, ice, and atmosphere, and can be used to
                                              ! spawn related internal variables in the ice model.
    type(coupler_1d_bc_type), &
                 optional, intent(in)   :: gas_fields_ocn  ! If present, this type describes the
                                              ! ocean and surface-ice fields that will participate
                                              ! in the calculation of additional gas or other
                                              ! tracer fluxes, and can be used to spawn related
                                              ! internal variables in the ice model.

    real, allocatable, dimension(:,:)   :: lonv, latv, rmask
    real, allocatable, dimension(:)     :: glon, glat
    real, allocatable, dimension(:)     :: xb, yb ! 1d global grid for diag_mgr
    real, allocatable, dimension(:,:)   :: tmpx, tmpy, tmp_2d
    integer                             :: io, ierr, unit, siz(4), logunit
    integer                             :: nlon, nlat, i, j, k
    character(len=80)                   :: control
    character(len=80)                   :: domainname
    character(len=256)                  :: err_mesg
    character(len=64)                   :: fname = 'INPUT/ice_model.res.nc'
    character(len=256)                  :: grid_filename='INPUT/grid_spec.nc'
    type(FmsNetcdfFile_t)               :: grid_fileobj, ocn_mosaic_fileobj
    character(len=256)                  :: ocn_mosaic_filename, tile_file
    character(len=256)                  :: axo_file      ! atmosXocean exchange grid file
    integer                             :: nx(1), ny(1)
    integer                             :: ntiles, n, m
    integer                             :: grid_version
    integer, parameter                  :: VERSION_0 = 0  ! grid file with field geolon_t
    integer, parameter                  :: VERSION_1 = 1  ! grid file with field x_T
    integer, parameter                  :: VERSION_2 = 2  ! mosaic file

    if(module_is_initialized) return

    read(input_nml_file, ice_model_nml, iostat=io)
    ierr = check_nml_error(io, 'ice_model_nml')

    call write_version_number (version, tagname)
    if ( mpp_pe() == mpp_root_pe() ) then
       logunit = stdlog()
       write (logunit, nml=ice_model_nml)
    endif

   !if (num_part /= 2) call error_mesg ('ice_model_init','this version must have num_part = 2', FATAL)
   !if (num_lev  /= 1) call error_mesg ('ice_model_init','this version must have num_lev = 1', FATAL)

    !--- get the grid size
    call get_grid_size('OCN',1,nlon,nlat)

    !-------------------- domain decomposition -----------------------------
    if( layout(1).EQ.0 .AND. layout(2).EQ.0 ) &
         call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes(), layout )
    if( layout(1).NE.0 .AND. layout(2).EQ.0 )layout(2) = mpp_npes()/layout(1)
    if( layout(1).EQ.0 .AND. layout(2).NE.0 )layout(1) = mpp_npes()/layout(2)
    domainname = 'AMIP Ice'
    call mpp_define_domains( (/1,nlon,1,nlat/), layout, Ice%Domain, name=domainname )
    call mpp_define_io_domain (Ice%Domain, io_layout)
    Ice%slow_Domain_NH = Ice%domain
    call mpp_get_compute_domain( Ice%Domain, isc, iec, jsc, jec )

    allocate ( glon     (nlon  ), glat     (nlat  )  )
    allocate ( lonv (isc:iec+1, jsc:jec+1) )
    allocate ( latv (isc:iec+1, jsc:jec+1) )
    allocate(  rmask(isc:iec,jsc:jec) )
    allocate(  Ice%glon_bnd(nlon+1),    Ice%glat_bnd(nlat+1)    )
    allocate ( Ice%lon_bnd (isc:iec+1),  &
               Ice%lat_bnd (jsc:jec+1),  &
               Ice%lon (isc:iec, jsc:jec), &
               Ice%lat (isc:iec, jsc:jec), &
               Ice%SST_C(isc:iec, jsc:jec) )

!  ---- set up local grid -----
   call get_grid_cell_vertices('OCN', 1, Ice%glon_bnd, Ice%glat_bnd)
   call get_grid_cell_centers ('OCN', 1, glon, glat)
   call get_grid_cell_vertices('OCN', 1, lonv, latv, Ice%domain)
   call get_grid_cell_centers ('OCN', 1, Ice%lon, Ice%lat, Ice%domain)

   !--- for conservation interpolation, the grid should be rectangular ----
   if(trim(interp_method) == "conservative" ) then
      err_mesg = 'Bilinear interpolation must be used for a tripolar grid'
      do i=isc, iec
         if(any(glon(i) /= Ice%lon(i,:)))  &
              call error_mesg ('ice_model_init',err_mesg,FATAL)
      enddo
      do j=jsc,jec
         if(any(glat(j) /= Ice%lat(:,j)))  &
              call error_mesg ('ice_model_init',err_mesg,FATAL)
      enddo
   endif

    !---------------- read ice cell areas from grid_spec.nc or ----------------
    !---------------- calculate the area for mosaic grid file  ----------------
    allocate (cell_area(isc:iec, jsc:jec))
    cell_area = 0.0
    call get_grid_cell_area('OCN', 1, cell_area, Ice%domain)

   !   ------ define Data masks ------
   if(open_file(grid_fileobj, grid_filename, 'read')) then

     if(variable_exists(grid_fileobj, 'wet')) then
       call read_data(grid_fileobj, "wet",      rmask)
     else if(variable_exists(grid_fileobj, 'ocn_mosaic_file') ) then ! read from mosaic file
       call read_data(grid_fileobj, "ocn_mosaic_file", ocn_mosaic_filename)
       ocn_mosaic_filename = "INPUT/"//trim(ocn_mosaic_filename)

       if(open_file(ocn_mosaic_fileobj, ocn_mosaic_filename, 'read')) then
         ntiles = get_mosaic_ntiles(ocn_mosaic_fileobj)
         if(ntiles .NE. 1) call error_mesg('ice_model_init', ' ntiles should be 1 for ocean mosaic, contact developer', FATAL)
         rmask = 0.0
         call get_grid_comp_area('OCN', 1, rmask, domain=Ice%Domain)
         rmask = rmask/cell_area
         do j = jsc, jec
            do i = isc, iec
               if(rmask(i,j) == 0.0) cell_area(i,j) = 0.0
            end do
         end do
       else
         call error_mesg("ice_model_init", "ocn_mosaic_file:"//trim(ocn_mosaic_filename)//"not found", FATAL)
       endif

     else
       call error_mesg('ice_model_init','wet and ocn_mosaic_file does not exist in file '//trim(grid_filename), FATAL )
     end if

   else
     call error_mesg("ice_model_init", "grid file: "//trim(grid_filename)//" not found", FATAL)
   end if

    !--- xb and yb is for diagnostics --------------------------------------
   allocate(xb(nlon+1), yb(nlat+1) )
   allocate ( tmpx(isc:iec+1, nlat+1) )
   call mpp_set_domain_symmetry(Ice%Domain, .TRUE.)
   call mpp_global_field(Ice%Domain, lonv, tmpx, flags=YUPDATE, position=CORNER)
   allocate ( tmp_2d(isc:iec+1, jsc:jec+1) )
   tmp_2d = 0
   tmp_2d(isc:iec+1,jsc) = sum(tmpx,2)/(nlat+1);
   deallocate(tmpx)
   allocate ( tmpx(nlon+1, jsc:jec+1) )

   call mpp_global_field(Ice%Domain, tmp_2d, tmpx, flags=XUPDATE, position=CORNER)
   xb = tmpx(:,jsc)
   deallocate(tmpx, tmp_2d)

   allocate ( tmpy(nlon+1, jsc:jec+1) )
   call mpp_global_field(Ice%Domain, latv, tmpy, flags=XUPDATE, position=CORNER)
   allocate ( tmp_2d(isc:iec+1, jsc:jec+1) )
   tmp_2d = 0
   tmp_2d(isc,jsc:jec+1) = sum(tmpy,1)/(nlon+1);
   deallocate(tmpy)
   allocate ( tmpy(isc:iec+1, nlat+1) )
   call mpp_global_field(Ice%Domain, tmp_2d, tmpy, flags=YUPDATE, position=CORNER)
   yb = tmpy(isc,:)
   deallocate(tmpy, tmp_2d)
   call mpp_set_domain_symmetry(Ice%Domain, .FALSE.)

    !--- define the ice data -----------------------------------------------
   Ice%lon = Ice%lon*pi/180.
   Ice%lat = Ice%lat*pi/180.
   Ice%glon_bnd = Ice%glon_bnd*pi/180.
   Ice%glat_bnd = Ice%glat_bnd*pi/180.
   Ice%lon_bnd    = Ice%glon_bnd(isc:iec+1)
   Ice%lat_bnd    = Ice%glat_bnd(jsc:jec+1)
    !-----------------------------------------------------------------------

    allocate ( Ice%ice_mask    (isc:iec, jsc:jec, num_part)         , &
               Ice%temp        (isc:iec, jsc:jec, num_part, num_lev), &
               Ice%part_size   (isc:iec, jsc:jec, num_part)         , &
               Ice%albedo      (isc:iec, jsc:jec, num_part)         , &
            Ice%albedo_vis_dir (isc:iec, jsc:jec, num_part)         , &
            Ice%albedo_nir_dir (isc:iec, jsc:jec, num_part)         , &
            Ice%albedo_vis_dif (isc:iec, jsc:jec, num_part)         , &
            Ice%albedo_nir_dif (isc:iec, jsc:jec, num_part)         , &
               Ice%rough_mom   (isc:iec, jsc:jec, num_part)         , &
               Ice%rough_heat  (isc:iec, jsc:jec, num_part)         , &
               Ice%rough_moist (isc:iec, jsc:jec, num_part)         , &
               Ice%u_surf      (isc:iec, jsc:jec, num_part)         , &
               Ice%v_surf      (isc:iec, jsc:jec, num_part)         , &
               Ice%thickness   (isc:iec, jsc:jec, num_part)         , &
               Ice%mask        (isc:iec, jsc:jec)                   , &
               Ice%SST_C(isc:iec, jsc:jec) )

    Ice%t_surf => Ice%temp (:,:,:,1)

    Ice%Time           = Time
    Ice%Time_init      = Time_init
    Ice%Time_step_fast = Time_step_fast
    Ice%Time_step_slow = Time_step_slow
    Ice%avg_kount = 0
    Ice%SST_C = 0.0
    Ice%mask = .false.
    where ( rmask > 0 ) Ice%mask = .true.


    allocate ( Ice%flux_u_bot  (isc:iec, jsc:jec, num_part) , &
               Ice%flux_v_bot  (isc:iec, jsc:jec, num_part) , &
               Ice%flux_t_bot  (isc:iec, jsc:jec, num_part) , &
               Ice%flux_q_bot  (isc:iec, jsc:jec, num_part) , &
               Ice%flux_lh_bot (isc:iec, jsc:jec, num_part) , &
               Ice%flux_sw_bot (isc:iec, jsc:jec, num_part) , &
        Ice%flux_sw_vis_bot    (isc:iec, jsc:jec, num_part) , &
        Ice%flux_sw_dir_bot    (isc:iec, jsc:jec, num_part) , &
        Ice%flux_sw_dif_bot    (isc:iec, jsc:jec, num_part) , &
        Ice%flux_sw_vis_dir_bot(isc:iec, jsc:jec, num_part) , &
        Ice%flux_sw_vis_dif_bot(isc:iec, jsc:jec, num_part) , &
        Ice%flux_sw_nir_dir_bot(isc:iec, jsc:jec, num_part) , &
        Ice%flux_sw_nir_dif_bot(isc:iec, jsc:jec, num_part) , &
               Ice%flux_lw_bot (isc:iec, jsc:jec, num_part) , &
               Ice%lprec_bot   (isc:iec, jsc:jec, num_part) , &
               Ice%fprec_bot   (isc:iec, jsc:jec, num_part) , &
               Ice%runoff_bot  (isc:iec, jsc:jec, num_part) , &
               Ice%frazil      (isc:iec, jsc:jec, num_part)   )

    allocate ( Ice%flux_u    (isc:iec, jsc:jec) , &
               Ice%flux_v    (isc:iec, jsc:jec) , &
               Ice%flux_t    (isc:iec, jsc:jec) , &
               Ice%flux_q    (isc:iec, jsc:jec) , &
               Ice%flux_lh   (isc:iec, jsc:jec) , &
               Ice%flux_sw   (isc:iec, jsc:jec) , &
         Ice%flux_sw_vis     (isc:iec, jsc:jec) , &
         Ice%flux_sw_dir     (isc:iec, jsc:jec) , &
         Ice%flux_sw_dif     (isc:iec, jsc:jec) , &
         Ice%flux_sw_vis_dir (isc:iec, jsc:jec) , &
         Ice%flux_sw_vis_dif (isc:iec, jsc:jec) , &
         Ice%flux_sw_nir_dir (isc:iec, jsc:jec) , &
         Ice%flux_sw_nir_dif (isc:iec, jsc:jec) , &
               Ice%flux_lw   (isc:iec, jsc:jec) , &
               Ice%lprec     (isc:iec, jsc:jec) , &
               Ice%fprec     (isc:iec, jsc:jec) , &
               Ice%p_surf    (isc:iec, jsc:jec) , &
               Ice%runoff    (isc:iec, jsc:jec) , &
               Ice%calving   (isc:iec, jsc:jec) , &
             Ice%runoff_hflx (isc:iec, jsc:jec) , &
             Ice%calving_hflx(isc:iec, jsc:jec) , &
             Ice%area        (isc:iec, jsc:jec) , &
             Ice%mi          (isc:iec, jsc:jec) , &
               Ice%flux_salt (isc:iec, jsc:jec)   )
Ice%flux_u = 0.0
Ice%flux_v = 0.0
Ice%flux_t = 0.0
Ice%flux_q = 0.0
Ice%flux_lh = 0.0
Ice%flux_sw = 0.0
Ice%flux_sw_vis = 0.0
Ice%flux_sw_dir = 0.0
Ice%flux_sw_dif = 0.0
Ice%flux_sw_vis_dir = 0.0
Ice%flux_sw_vis_dif = 0.0
Ice%flux_sw_nir_dir = 0.0
Ice%flux_sw_nir_dif = 0.0
Ice%flux_lw = 0.0
Ice%lprec = 0.0
Ice%fprec = 0.0
Ice%p_surf = 0.0
Ice%runoff = 0.0
Ice%calving = 0.0
Ice%runoff_hflx = 0.0
Ice%calving_hflx = 0.0
Ice%area = 0.0
Ice%mi = 0.0
Ice%flux_salt = 0.0

Ice%area = cell_area  * 4*PI*RADIUS*RADIUS
Ice%mi   = 0.0


    !-----------------------------------------------------------------------
    !  -------- read restart --------
!if (do_netcdf_restart) call ice_register_restart(Ice, 'ice_model.res.nc')

if (open_file(Ice_restart, "INPUT/ice_model.res.nc", 'read', Ice%domain)) then
  !if (mpp_pe() == mpp_root_pe()) call error_mesg ('ice_model_mod', &
  !         'Reading NetCDF formatted restart file: INPUT/ice_model.res.nc', NOTE)

   call read_new_restart(Ice_restart)

   ! have to double check variable name, before id_restart_albedo was passed in uninitialized
   if (variable_exists(Ice_restart, "albedo")) then
      if (mpp_pe() == mpp_root_pe()) call error_mesg ('ice_model_mod', &
                'Initializing diffuse and direct streams to albedo', NOTE)
    ! initialize ocean values - ice values initialized below
      Ice%albedo_vis_dir(:,:,1) = Ice%albedo(:,:,1)
      Ice%albedo_nir_dir(:,:,1) = Ice%albedo(:,:,1)
      Ice%albedo_vis_dif(:,:,1) = Ice%albedo(:,:,1)
      Ice%albedo_nir_dif(:,:,1) = Ice%albedo(:,:,1)
   endif

else
    if (file_exists('INPUT/ice_model.res')) then
       if (mpp_pe() == mpp_root_pe()) call error_mesg ('ice_model_mod', &
            'Reading native formatted restart file.', NOTE)

       open(newunit=unit, file='INPUT/ice_model.res')

       read (unit) control

       ! must use correct restart version with native format
       if (trim(control) /= trim(restart_format)) call error_mesg &
            ('ice_model_init', 'invalid restart format', FATAL)

       read  (unit) mlon, mlat, npart

       !     ---- restart resolution must be consistent with input args ---
       if (mlon /= nlon .or. mlat /= nlat .or. npart /= 2)  &
            call error_mesg ('ice_model_init',           &
            'incorrect resolution on restart', FATAL)

       read(unit) Ice%part_size
       read(unit) Ice%temp
       read(unit) Ice%thickness
       read(unit) Ice%albedo

       read(unit) Ice%albedo_vis_dir
       read(unit) Ice%albedo_nir_dir
       read(unit) Ice%albedo_vis_dif
       read(unit) Ice%albedo_nir_dif

       read(unit) Ice%rough_mom
       read(unit) Ice%rough_heat
       read(unit) Ice%rough_moist
       read(unit) Ice%u_surf
       read(unit) Ice%v_surf
       read(unit) Ice%frazil
       read(unit) Ice%flux_u_bot
       read(unit) Ice%flux_v_bot

       close(unit)

    else

       !     ----- no restart - no ice -----

       Ice%temp      = tfreeze + temp_ice_freeze
       Ice%thickness = 0.0
       Ice%part_size         = 0.0
       Ice%part_size (:,:,1) = 1.0
       Ice%albedo     = 0.14
     ! initialize ocean values - ice values initialized below
       Ice%albedo_vis_dir(:,:,1) = Ice%albedo(:,:,1)
       Ice%albedo_nir_dir(:,:,1) = Ice%albedo(:,:,1)
       Ice%albedo_vis_dif(:,:,1) = Ice%albedo(:,:,1)
       Ice%albedo_nir_dif(:,:,1) = Ice%albedo(:,:,1)
       Ice%rough_mom  = 0.0004
       Ice%rough_heat = 0.0004
       Ice%rough_moist= 0.0004
       Ice%u_surf     = 0.0
       Ice%v_surf     = 0.0
       Ice%frazil     = 0.0

       !     --- open water roughness (initially fixed) ---

       call fixed_ocean_roughness ( Ice%mask, Ice%rough_mom  (:,:,1), &
            Ice%rough_heat (:,:,1), &
            Ice%rough_moist(:,:,1)  )

    endif
endif

!! set ice partiton values to that of Ice%albedo.
   Ice%albedo_vis_dir(:,:,2) = Ice%albedo(:,:,2)
   Ice%albedo_nir_dir(:,:,2) = Ice%albedo(:,:,2)
   Ice%albedo_vis_dif(:,:,2) = Ice%albedo(:,:,2)
   Ice%albedo_nir_dif(:,:,2) = Ice%albedo(:,:,2)

    !---- initialization of ice mask (actually where ice exists) -----
    Ice%ice_mask = .false.
    where (Ice%mask     (:,:)           .and.     &
         Ice%part_size(:,:,2) > 0.0   .and.     &
         Ice%thickness(:,:,2) >= thickness_min) &
         Ice%ice_mask(:,:,2) = .true.


    call ice_albedo_init (tfreeze)

!   if(trim(interp_method) == "conservative") then
!      Amip = amip_interp_new ( Ice%lon_bnd,     Ice%lat_bnd,  &
!           Ice%mask(:,:),                 &
!           interp_method = interp_method, &
!           use_climo=use_climo_ice,       &
!           use_annual=use_annual_ice     )
!   else if(trim(interp_method) == "bilinear") then
!      Amip = amip_interp_new ( Ice%lon,     Ice%lat,          &
!           Ice%mask(:,:),                 &
!           interp_method = interp_method, &
!           use_climo=use_climo_ice,       &
!           use_annual=use_annual_ice     )
!   else
!      call error_mesg('ice_model_init', 'interp_method should be conservative or bilinear', &
!           FATAL)
!   endif

    ! --- diagnostics ---

!   call ice_diag_init (Ice, xb, yb)

    !Balaji
    ! iceClock computation added for sis2 interface compatibility
    iceClock = mpp_clock_id( 'Ice', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    iceClock1 = mpp_clock_id( 'Ice: bot to top', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    iceClock3 = mpp_clock_id( 'Ice: update fast', flags=clock_flag_default, grain=CLOCK_ROUTINE )


    !--- release the memory ------------------------------------------------
    deallocate(lonv, latv, glon, glat, rmask, xb, yb )

    module_is_initialized = .true.

  end subroutine ice_model_init
!=============================================================================================
subroutine update_ice_model_fast_old( Ice, Atmos_boundary)
type (ice_data_type), intent(inout) :: Ice
type(atmos_ice_boundary_type), intent(in) :: Atmos_boundary

return
end subroutine update_ice_model_fast_old
!=============================================================================================

subroutine unpack_ocean_ice_boundary(Ocean_boundary, Ice)
  type(ocean_ice_boundary_type),  intent(inout) :: Ocean_boundary
  type(ice_data_type),            intent(inout) :: Ice

 call ice_bottom_to_ice_top (Ice, Ocean_boundary%t, Ocean_boundary%u, Ocean_boundary%v,        &
                             Ocean_boundary%frazil, Ocean_boundary, Ocean_boundary%s, Ocean_boundary%sea_level)
  return
end subroutine unpack_ocean_ice_boundary
!=============================================================================================
subroutine update_ice_model_slow( Ice, runoff, calving, runoff_hflx, calving_hflx, p_surf)
type(ice_data_type),           intent(inout) :: Ice
real, dimension(:,:), intent(in), optional :: runoff, calving
real, dimension(:,:), intent(in), optional :: runoff_hflx, calving_hflx
real, dimension(:,:),           intent(in), optional :: p_surf

return
end subroutine update_ice_model_slow
!=============================================================================================
subroutine ice_model_end(Ice)
type(ice_data_type), intent(inout) :: Ice
character(len=64) :: fname='RESTART/ice_model.res.nc'
integer :: unit

if(.not.module_is_initialized) return

module_is_initialized = .false.

end subroutine ice_model_end

!#######################################################################
subroutine ice_model_restart(Ice, time_stamp)
  type (ice_data_type),     intent(inout), optional :: Ice
  character(len=*),         intent(in), optional :: time_stamp

  return
end subroutine ice_model_restart

!=============================================================================================

 subroutine ice_stock_pe(Ice, index, value)
 type(ice_data_type), intent(in) :: Ice
 integer, intent(in) :: index
 real, intent(out)   :: value

 value = 0.0
 return
 end subroutine ice_stock_pe
!=============================================================================================

subroutine ice_data_type_chksum(id, timestep, data_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_data_type), intent(in) :: data_type

    return
end subroutine ice_data_type_chksum

!=============================================================================================

subroutine ocn_ice_bnd_type_chksum(id, timestep, bnd_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ocean_ice_boundary_type), intent(in) :: bnd_type

    return
end subroutine ocn_ice_bnd_type_chksum

!=============================================================================================
subroutine atm_ice_bnd_type_chksum(id, timestep, bnd_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(atmos_ice_boundary_type), intent(in) :: bnd_type

    return
end subroutine atm_ice_bnd_type_chksum
!=============================================================================================
subroutine lnd_ice_bnd_type_chksum(id, timestep, bnd_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_ice_boundary_type), intent(in) :: bnd_type

    return
end subroutine lnd_ice_bnd_type_chksum
!=============================================================================================
subroutine update_ice_atm_deposition_flux(Ice_boundary, Ice )
  type(atmos_ice_boundary_type), intent(inout):: Ice_boundary
  type(ice_data_type),           intent(inout):: Ice
  !Do nothing! This interface is only active for newer coupled models that use SIS2
end subroutine update_ice_atm_deposition_flux
!=============================================================================================
subroutine share_ice_domains(Ice)
  type(ice_data_type), intent(inout) :: Ice
  ! This is a null interface
end subroutine share_ice_domains
!=============================================================================================
subroutine exchange_slow_to_fast_ice(Ice)
  type(ice_data_type), intent(in) :: Ice
end subroutine exchange_slow_to_fast_ice
!=============================================================================================
subroutine exchange_fast_to_slow_ice(Ice)
  type(ice_data_type), intent(in) :: Ice
end subroutine exchange_fast_to_slow_ice
!=============================================================================================
subroutine set_ice_surface_fields(Ice)
  type(ice_data_type), intent(inout) :: Ice
  ! Null routines to provide similar interfaces to those used by SIS2.
  ! In SIS this functionality is combined into ice_bottom_to_ice_top.
end subroutine set_ice_surface_fields
!=============================================================================================
subroutine ice_model_fast_cleanup(Ice)
  type(ice_data_type), intent(in) :: Ice
end subroutine ice_model_fast_cleanup
!=============================================================================================
subroutine unpack_land_ice_boundary(Ice, LIB)
  type(ice_data_type),          intent(inout) :: Ice ! The publicly visible ice data type.
  type(land_ice_boundary_type), intent(in)    :: LIB ! The land ice boundary type that is being unpacked.
! logical :: sent
! integer :: i, j, i_off, j_off

! if (associated(LIB%runoff)) then ! save liquid runoff for ocean
!    i_off = LBOUND(LIB%runoff,1) - isc
!    j_off = LBOUND(LIB%runoff,2) - jsc
!    do j = jsc, jec
!       do i = isc, iec
!          Ice % runoff(i,j)  = LIB%runoff(i+i_off, j+j_off)
!       enddo
!    enddo
! else
!    do j = jsc, jec
!       do i = isc, iec
!          Ice % runoff(i,j) = 0.0
!       enddo
!    enddo
! endif

! if (associated(LIB%calving)) then ! save frozen runoff for ocean
!    i_off = LBOUND(LIB%calving,1) - isc
!    j_off = LBOUND(LIB%calving,2) - jsc
!    do j = jsc, jec
!       do i = isc, iec
!          Ice % calving(i,j) = LIB%calving(i+i_off, j+j_off)
!       enddo
!    enddo
! else
!    do j = jsc, jec
!       do i = isc, iec
!          Ice % calving(i,j) = 0.0
!       enddo
!    enddo
! endif

! if (associated(LIB%runoff_hflx)) then ! save liquid runoff hflx for ocean
!    do j = jsc, jec
!       do i = isc, iec
!          Ice % runoff_hflx(i,j)  = LIB%runoff_hflx(i+i_off, j+j_off)
!       enddo
!    enddo
! else
!    do j = jsc, jec
!       do i = isc, iec
!          Ice % runoff_hflx(i,j) = 0.0
!       enddo
!    enddo
! endif

! if (associated(LIB%calving_hflx)) then ! save frozen runoff hflx for ocean
!    do j = jsc, jec
!       do i = isc, iec
!          Ice % calving_hflx(i,j) = LIB%calving_hflx(i+i_off, j+j_off)
!       enddo
!    enddo
! else
!    do j = jsc, jec
!       do i = isc, iec
!          Ice % calving_hflx(i,j) = 0.0
!       enddo
!    enddo
! endif
end subroutine unpack_land_ice_boundary
!=============================================================================================
!\brief: dummy version of ice_sis routine that is called by the new interfaces,
!! but is not used by ice_null.
 subroutine ice_bottom_to_ice_top(Ice, t_surf_ice_bot, u_surf_ice_bot, v_surf_ice_bot, &
                                    frazil_ice_bot, Ocean_ice_boundary, s_surf_ice_bot, sea_lev_ice_bot )
    type (ice_data_type),                    intent(inout) :: Ice
    real, dimension(isc:iec,jsc:jec),           intent(in) :: t_surf_ice_bot, u_surf_ice_bot
    real, dimension(isc:iec,jsc:jec),           intent(in) :: v_surf_ice_bot, frazil_ice_bot
    type(ocean_ice_boundary_type),           intent(inout) :: Ocean_ice_boundary
    real, dimension(isc:iec,jsc:jec), intent(in), optional :: s_surf_ice_bot, sea_lev_ice_bot
end subroutine ice_bottom_to_ice_top
!=============================================================================================
!\brief new coupler interface to provide ocean surface data to atmosphere
! calls a dummy version of ice_bottom_to_ice_top
subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_boundary
    type (ice_data_type),          intent(inout) :: Ice

    call mpp_clock_begin(iceClock)
    call mpp_clock_begin(iceClock1)
    call ice_bottom_to_ice_top(Ice, Ocean_boundary%t, Ocean_boundary%u, Ocean_boundary%v,        &
                                Ocean_boundary%frazil, Ocean_boundary, Ocean_boundary%s, Ocean_boundary%sea_level  )
    call mpp_clock_end(iceClock1)
    call mpp_clock_end(iceClock)

end subroutine update_ice_model_slow_up
!=============================================================================================
!\brief new coupler interface to do slow ice processes:  dynamics, transport, mass
!
subroutine update_ice_model_slow_dn ( Atmos_boundary, Land_boundary, Ice )
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
    type(land_ice_boundary_type),  intent(inout) :: Land_boundary
    type (ice_data_type),          intent(inout) :: Ice

    call unpack_land_ice_boundary( Ice, Land_boundary )

    call update_ice_model_slow (Ice, Ice%runoff, Ice%calving, Ice%runoff_hflx, Ice%calving_hflx, Ice%p_surf)

end subroutine update_ice_model_slow_dn
!=============================================================================================
!\brief new coupler interface to do fast ice processes
!\note the update_ice_model_fast_old routine is different from the one in ice_sis
! However, the interface is the same in both versions of ice_model.F90.
subroutine update_ice_model_fast_new ( Atmos_boundary, Ice )
   type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
   type (ice_data_type),          intent(inout) :: Ice
   call mpp_clock_begin(iceClock)
   call mpp_clock_begin(iceClock3)

   call update_ice_model_fast_old(Ice, Atmos_boundary)

   call mpp_clock_end(iceClock3)
   call mpp_clock_end(iceClock)

end subroutine update_ice_model_fast_new
!=============================================================================================
end module ice_model_mod
