module ice_model_mod

  use mpp_domains_mod,  only: mpp_update_domains, domain2D
  use fms_mod,          only: error_mesg, open_namelist_file, stdlog, check_nml_error, close_file
  use time_manager_mod, only: time_type, get_time
  use coupler_types_mod,only: coupler_2d_bc_type, coupler_3d_bc_type
  use constants_mod,    only: PI, RADIUS

  use ice_grid_mod,     only: t_to_uv, set_ice_grid, cell_area, wett
  use ice_grid_mod,     only: isc, iec, jsc, jec, isd, ied, jsd, jed, km

  implicit none
  private

  public :: ice_data_type
  public :: ocean_ice_boundary_type
  public :: atmos_ice_boundary_type
  public :: land_ice_boundary_type
  public :: ice_model_init
  public :: ice_model_end
  public :: ice_bottom_to_ice_top
  public :: update_ice_model_fast
  public :: ice_stock_pe
  public :: cell_area
  public :: update_ice_model_slow
  public :: update_ice_model_slow_up
  public :: update_ice_model_slow_dn
  public :: ice_model_restart
  public :: ocn_ice_bnd_type_chksum
  public :: atm_ice_bnd_type_chksum
  public :: lnd_ice_bnd_type_chksum
  public :: ice_data_type_chksum

  !
  ! the following four types are for data exchange with the new coupler
  ! they are defined here but declared in coupler_main and allocated in flux_init
  !
  type ice_data_type
     type(domain2D)                     :: Domain
     type(time_type)                    :: Time_Init, Time
     type(time_type)                    :: Time_step_fast, Time_step_slow
     integer                            :: avg_kount
     logical                            :: pe
     integer, pointer, dimension(:)     :: pelist              =>NULL() ! Used for flux-exchange.
     logical, pointer, dimension(:,:)   :: mask                =>NULL() ! where ice can be
     logical, pointer, dimension(:,:,:) :: ice_mask            =>NULL() ! where ice actually is
     real,    pointer, dimension(:,:,:) :: part_size           =>NULL()
     real,    pointer, dimension(:,:,:) :: part_size_uv        =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo              =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_vis_dir      =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_nir_dir      =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_vis_dif      =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_nir_dif      =>NULL()
     real,    pointer, dimension(:,:,:) :: rough_mom           =>NULL()
     real,    pointer, dimension(:,:,:) :: rough_heat          =>NULL()
     real,    pointer, dimension(:,:,:) :: rough_moist         =>NULL()
     real,    pointer, dimension(:,:,:) :: t_surf              =>NULL()
     real,    pointer, dimension(:,:,:) :: u_surf              =>NULL()
     real,    pointer, dimension(:,:,:) :: v_surf              =>NULL()
     real,    pointer, dimension(:,:)   :: sea_lev             =>NULL()
     real,    pointer, dimension(:,:)   :: s_surf              =>NULL()
     real,    pointer, dimension(:,:)   :: u_ocn               =>NULL()
     real,    pointer, dimension(:,:)   :: v_ocn               =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_u_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_v_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_t_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_q_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_lw_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_vis_dir_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_vis_dif_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_nir_dir_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_nir_dif_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_lh_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: lprec_top           =>NULL()
     real,    pointer, dimension(:,:,:) :: fprec_top           =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_u              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_v              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_t              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_q              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_lw             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_vis_dir     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_vis_dif     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_nir_dir     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_nir_dif     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_lh             =>NULL()
     real,    pointer, dimension(:,:  ) :: lprec               =>NULL()
     real,    pointer, dimension(:,:  ) :: fprec               =>NULL()
     real,    pointer, dimension(:,:  ) :: p_surf              =>NULL()
     real,    pointer, dimension(:,:  ) :: runoff              =>NULL()
     real,    pointer, dimension(:,:  ) :: calving             =>NULL()
     real,    pointer, dimension(:,:  ) :: runoff_hflx         =>NULL()
     real,    pointer, dimension(:,:  ) :: calving_hflx        =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_salt           =>NULL()
     real,    pointer, dimension(:,:)   :: lwdn                =>NULL()
     real,    pointer, dimension(:,:  ) :: swdn                =>NULL() ! downward long/shortwave
     real,    pointer, dimension(:,:,:) :: pen                 =>NULL()
     real,    pointer, dimension(:,:,:) :: trn                 =>NULL() ! ice optical parameters
     real,    pointer, dimension(:,:,:) :: tmelt               =>NULL()
     real,    pointer, dimension(:,:,:) :: bmelt               =>NULL()
     real,    pointer, dimension(:,:,:) :: h_snow              =>NULL()
     real,    pointer, dimension(:,:,:) :: h_ice               =>NULL()
     real,    pointer, dimension(:,:,:) :: t_ice1              =>NULL()
     real,    pointer, dimension(:,:,:) :: t_ice2              =>NULL()
     real,    pointer, dimension(:,:)   :: u_ice               =>NULL()
     real,    pointer, dimension(:,:)   :: v_ice               =>NULL()
     real,    pointer, dimension(:,:)   :: sig11               =>NULL()
     real,    pointer, dimension(:,:)   :: sig22               =>NULL()
     real,    pointer, dimension(:,:)   :: sig12               =>NULL()
     real,    pointer, dimension(:,:)   :: frazil              =>NULL()
     real,    pointer, dimension(:,:)   :: bheat               =>NULL()
     real,    pointer, dimension(:,:)   :: qflx_lim_ice        =>NULL()
     real,    pointer, dimension(:,:)   :: qflx_res_ice        =>NULL()
     real,    pointer, dimension(:,:)   :: area                =>NULL()
     real,    pointer, dimension(:,:)   :: mi                  =>NULL() ! This is needed for the wave model. It is introduced here,
                                                                        ! because flux_ice_to_ocean cannot handle 3D fields. This may be
                                                                        ! removed, if the information on ice thickness can be derived from
                                                                        ! eventually from h_ice outside the ice module.
     logical, pointer, dimension(:,:)   :: maskmap             =>NULL() ! A pointer to an array indicating which
                                                                        ! logical processors are actually used for
                                                                        ! the ocean code. The other logical
                                                                        ! processors would be all land points and
                                                                        ! are not assigned to actual processors.
                                                                        ! This need not be assigned if all logical
                                                                        ! processors are used
     integer, dimension(3)              :: axes
     type(coupler_3d_bc_type)           :: ocean_fields       ! array of fields used for additional tracers
     type(coupler_2d_bc_type)           :: ocean_fluxes       ! array of fluxes used for additional tracers
     type(coupler_3d_bc_type)           :: ocean_fluxes_top   ! array of fluxes for averaging
  end type ice_data_type
  type :: ocean_ice_boundary_type
     real, dimension(:,:),   pointer :: u         =>NULL()
     real, dimension(:,:),   pointer :: v         =>NULL()
     real, dimension(:,:),   pointer :: t         =>NULL()
     real, dimension(:,:),   pointer :: s         =>NULL()
     real, dimension(:,:),   pointer :: frazil    =>NULL()
     real, dimension(:,:),   pointer :: sea_level =>NULL()
     real, dimension(:,:,:), pointer :: data      =>NULL() ! collective field for "named" fields above
     integer                         :: xtype              ! REGRID, REDIST or DIRECT used by coupler
     type(coupler_2d_bc_type)        :: fields     ! array of fields used for additional tracers
  end type 

  type :: atmos_ice_boundary_type 
     real, dimension(:,:,:), pointer :: u_flux  =>NULL()
     real, dimension(:,:,:), pointer :: v_flux  =>NULL()
     real, dimension(:,:,:), pointer :: u_star  =>NULL()
     real, dimension(:,:,:), pointer :: t_flux  =>NULL()
     real, dimension(:,:,:), pointer :: q_flux  =>NULL()
     real, dimension(:,:,:), pointer :: lw_flux =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_vis_dir =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_vis_dif =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_nir_dir =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_nir_dif =>NULL()
     real, dimension(:,:,:), pointer :: lprec   =>NULL()
     real, dimension(:,:,:), pointer :: fprec   =>NULL()
     real, dimension(:,:,:), pointer :: dhdt    =>NULL()
     real, dimension(:,:,:), pointer :: dedt    =>NULL()
     real, dimension(:,:,:), pointer :: drdt    =>NULL()
     real, dimension(:,:,:), pointer :: coszen  =>NULL()
     real, dimension(:,:,:), pointer :: p       =>NULL()
     real, dimension(:,:,:), pointer :: data    =>NULL()
     integer                         :: xtype
     type(coupler_3d_bc_type)        :: fluxes     ! array of fluxes used for additional tracers
  end type

  type :: land_ice_boundary_type
     real, dimension(:,:),   pointer :: runoff  =>NULL()
     real, dimension(:,:),   pointer :: calving =>NULL()
     real, dimension(:,:),   pointer :: runoff_hflx  =>NULL()
     real, dimension(:,:),   pointer :: calving_hflx =>NULL()
     real, dimension(:,:,:), pointer :: data    =>NULL() ! collective field for "named" fields above
     integer                         :: xtype            ! REGRID, REDIST or DIRECT used by coupler
  end type

  !--- namelist interface --------------
  real    :: mom_rough_ice  = 1.0e-4     ! momentum same, cd10=(von_k/ln(10/z0))^2
  real    :: heat_rough_ice = 1.0e-4     ! heat roughness length
  real    :: kmelt          = 6e-5*4e6   ! ocean/ice heat flux constant
  real    :: ks             = 0.31       ! snow conductivity (W/mK)
  real    :: alb_sno        = 0.85       ! snow albedo (less if melting)
  real    :: alb_ice        = 0.5826     ! ice albedo (less if melting)
  real    :: pen_ice        = 0.3        ! part unreflected solar penetrates ice
  real    :: opt_dep_ice    = 0.67       ! ice optical depth
  real    :: t_range_melt   = 1.0        ! melt albedos scaled in over T range
  real    :: ice_bulk_salin = 0.0        ! ice bulk salinity (for ocean salt flux)
  real    :: p0             = 2.75e4     ! ice strength parameter
  real    :: c0             = 20.0       ! another ice strength parameter
  real    :: cdw            = 3.24e-3    ! water/ice drag coefficient
  real    :: wd_turn        = 25.0       ! water/ice drag turning angle
  real    :: h_lo_lim       = 0.0        ! min ice thickness for temp. calc.
  integer :: nsteps_dyn     = 432        ! dynamics steps per slow timestep
  integer :: nsteps_adv     = 8          ! advection steps per slow timestep
  integer :: num_part       = 6          ! number of ice grid partitions
                                         ! partition 1 is open water
                                         ! partitions 2 to num_part-1 are
                                         !   thickness limited ice categories
                                         ! partition num_part is unlimited ice
  logical :: atmos_winds = .true.        ! wind stress from atmosphere model over t points and has wrong sign
  logical :: slab_ice    = .false.       ! do old-style GFDL slab ice?
  logical :: spec_ice    = .false.       ! old-style GFDL slab ice with SST, ice thickness and conc. from data
  logical :: do_ice_restore  = .false.   ! restore sea-ice toward climatology
  logical :: do_ice_limit    = .false.   ! limit sea ice to max_ice_limit
  real    :: max_ice_limit   = 4.0       ! maximum sea ice height(m),
                                         ! if do_ice_limit is true
                                         ! TK: default chosen based on observed
                                         !     ice thickness data used by climate
                                         !     group, which range up to 7.31 m
  real    :: ice_restore_timescale = 5.0 ! time scale for restoring ice (days)
  logical :: conservation_check = .true. ! check for heat and h2o conservation
  logical :: slp2ocean          = .false.! apply sea level pressure to ocean surface
  logical :: cm2_bugs           = .false.! keep cm2 bugs for reproducibility        
  logical :: verbose            = .false.! control printing message, will slow model down when turn true
  logical :: do_icebergs        = .false.! call iceberg code to modify calving field
  logical :: add_diurnal_sw     = .false.! apply an additional diurnal cycle to shortwave radiation
  integer :: layout(2)          = (/0, 0/)
  integer :: io_layout(2)       = (/0, 0/)

  real    :: channel_viscosity  = 0.     ! viscosity used in one-cell wide channels to parameterize transport (m^2/s)
  real    :: smag_ocn           = 0.15   ! Smagorinksy coefficient for viscosity (dimensionless)
  real    :: ssh_gravity        = 9.81   ! Gravity parameter used in channel viscosity parameterization (m/s^2)
  real    :: chan_cfl_limit     = 0.25   ! CFL limit for channel viscosity parameterization (dimensionless)
  logical :: do_sun_angle_for_alb = .false.! find the sun angle for ocean albed in the frame of the ice model
  ! mask_table contains information for masking domain ( n_mask, layout and mask_list).
  !   A text file to specify n_mask, layout and mask_list to reduce number of processor
  !   usage by masking out some domain regions which contain all land points.
  !   The default file name of mask_table is "INPUT/ice_mask_table". Please note that
  !   the file name must begin with "INPUT/". The first
  !   line of mask_table will be number of region to be masked out. The second line
  !   of the mask_table will be the layout of the model. User need to set ice_model_nml
  !   variable layout to be the same as the second line of the mask table.
  !   The following n_mask line will be the position of the processor to be masked out.
  !   The mask_table could be created by tools check_mask.
  !   For example the mask_table will be as following if n_mask=2, layout=4,6 and
  !   the processor (1,2) and (3,6) will be masked out.
  !     2
  !     4,6
  !     1,2
  !     3,6

  character(len=128) :: mask_table = "INPUT/ice_mask_table"


  namelist /ice_model_nml/ mom_rough_ice, heat_rough_ice, p0, c0, cdw, wd_turn,  &
                           kmelt, alb_sno, alb_ice, pen_ice, opt_dep_ice,        &
                           nsteps_dyn, nsteps_adv, num_part, atmos_winds,        &
                           slab_ice, spec_ice, ice_bulk_salin, layout,           & 
                           do_ice_restore, do_ice_limit, max_ice_limit,          &
                           ice_restore_timescale, slp2ocean, conservation_check, &
                           t_range_melt, cm2_bugs, ks, h_lo_lim, verbose,        &
                           do_icebergs, add_diurnal_sw, io_layout, channel_viscosity,&
                           smag_ocn, ssh_gravity, chan_cfl_limit, do_sun_angle_for_alb, &
                           mask_table


contains

  subroutine ice_model_init (Ice, Time_Init, Time, Time_step_fast, Time_step_slow )

    type (ice_data_type), intent(inout) :: Ice
    type (time_type)    , intent(in)    :: Time_Init      ! starting time of model integration
    type (time_type)    , intent(in)    :: Time           ! current time
    type (time_type)    , intent(in)    :: Time_step_fast ! time step for the ice_model_fast
    type (time_type)    , intent(in)    :: Time_step_slow ! time step for the ice_model_slow

    integer           :: io, ierr, nlon, nlat, npart, unit, log_unit, k
    integer           :: sc, dy, i, j, lay_out(2)
    integer           :: id_restart, id_restart_albedo, id_restart_flux_sw
    real              :: dt_slow
    character(len=64) :: restart_file
    !
    ! read namelist and write to logfile
    !
    unit = open_namelist_file()
    read  (unit, ice_model_nml,iostat=io)
    write (stdlog(), ice_model_nml)
    ierr = check_nml_error(io, 'ice_model_nml')
    call close_file(unit)

    if (spec_ice) then
       slab_ice = .true.
       nsteps_dyn = 0
       nsteps_adv = 0
    end if
    if (slab_ice) num_part = 2 ! open water and ice ... but never in same place

    call get_time(Time_step_slow, sc, dy); dt_slow=864e2*dy+sc

    if( ASSOCIATED(Ice%maskmap) ) then
       lay_out(1) = size(Ice%maskmap,1); lay_out(2) = size(Ice%maskmap,2)
       call set_ice_grid(Ice%domain, dt_slow, nsteps_dyn, nsteps_adv, num_part, lay_out, io_layout, Ice%maskmap  )
    else
       call set_ice_grid(Ice%domain, dt_slow, nsteps_dyn, nsteps_adv, num_part, layout, io_layout )
    end if

    allocate ( Ice % mask     (isc:iec, jsc:jec)       , &
         Ice % ice_mask       (isc:iec, jsc:jec, km)   , &
         Ice % t_surf         (isc:iec, jsc:jec, km)   , &
         Ice % s_surf         (isc:iec, jsc:jec)       , &
         Ice % sea_lev        (isd:ied, jsd:jed)       , &
         Ice % part_size      (isd:ied, jsd:jed, km)   , &
         Ice % part_size_uv   (isc:iec, jsc:jec, km)   , &
         Ice % u_surf         (isc:iec, jsc:jec, km)   , &
         Ice % v_surf         (isc:iec, jsc:jec, km)   , &
         Ice % u_ocn          (isd:ied, jsd:jed)       , &
         Ice % v_ocn          (isd:ied, jsd:jed)       , &
         Ice % rough_mom      (isc:iec, jsc:jec, km)   , &
         Ice % rough_heat     (isc:iec, jsc:jec, km)   , &
         Ice % rough_moist    (isc:iec, jsc:jec, km)   , &
         Ice % albedo         (isc:iec, jsc:jec, km)   , &                
         Ice % albedo_vis_dir (isc:iec, jsc:jec, km)   , &
         Ice % albedo_nir_dir (isc:iec, jsc:jec, km)   , &
         Ice % albedo_vis_dif (isc:iec, jsc:jec, km)   , &
         Ice % albedo_nir_dif (isc:iec, jsc:jec, km)   )

    allocate ( Ice % flux_u_top   (isd:ied, jsd:jed, km) ,       &
         Ice % flux_v_top         (isd:ied, jsd:jed, km) ,       &
         Ice % flux_t_top         (isc:iec, jsc:jec, km) ,       &
         Ice % flux_q_top         (isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_vis_dir_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_vis_dif_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_nir_dir_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_nir_dif_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_lw_top        (isc:iec, jsc:jec, km) ,       &
         Ice % flux_lh_top        (isc:iec, jsc:jec, km) ,       &
         Ice % lprec_top          (isc:iec, jsc:jec, km) ,       &
         Ice % fprec_top          (isc:iec, jsc:jec, km)   )

    allocate ( Ice % flux_u    (isc:iec, jsc:jec ) ,       &
         Ice % flux_v          (isc:iec, jsc:jec ) ,       &
         Ice % flux_t          (isc:iec, jsc:jec ) ,       &
         Ice % flux_q          (isc:iec, jsc:jec ) ,       &
         Ice % flux_sw_vis_dir (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_vis_dif (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_nir_dir (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_nir_dif (isc:iec, jsc:jec) ,        &
         Ice % flux_lw         (isc:iec, jsc:jec ) ,       &
         Ice % flux_lh         (isc:iec, jsc:jec ) ,       &
         Ice % lprec           (isc:iec, jsc:jec ) ,       &
         Ice % fprec           (isc:iec, jsc:jec ) ,       &
         Ice % p_surf          (isc:iec, jsc:jec ) ,       &
         Ice % runoff          (isc:iec, jsc:jec ) ,       &
         Ice % calving         (isc:iec, jsc:jec ) ,       &
         Ice % runoff_hflx     (isc:iec, jsc:jec ) ,       &
         Ice % calving_hflx    (isc:iec, jsc:jec ) ,       &
         Ice % flux_salt       (isc:iec, jsc:jec ) ,       &
         Ice % lwdn            (isc:iec, jsc:jec ) ,       &
         Ice % swdn            (isc:iec, jsc:jec )         )
    allocate ( Ice % frazil (isc:iec, jsc:jec), Ice % bheat  (isc:iec, jsc:jec), &
               Ice % u_ice  (isd:ied, jsd:jed), Ice % v_ice  (isd:ied, jsd:jed), &
               Ice % sig11  (isd:ied, jsd:jed), Ice % sig22  (isd:ied, jsd:jed), &
               Ice % sig12  (isd:ied, jsd:jed)                               )
    allocate ( Ice % tmelt  (isc:iec, jsc:jec, 2:km), Ice % bmelt  (isc:iec, jsc:jec, 2:km) , &
               Ice % pen    (isc:iec, jsc:jec, 2:km), Ice % trn    (isc:iec, jsc:jec, 2:km) , &
               Ice % h_snow (isd:ied, jsd:jed, 2:km), Ice % h_ice  (isd:ied, jsd:jed, 2:km) , &
               Ice % t_ice1 (isd:ied, jsd:jed, 2:km), Ice % t_ice2 (isd:ied, jsd:jed, 2:km)   )
    allocate ( Ice % qflx_lim_ice  (isc:iec, jsc:jec) , Ice % qflx_res_ice  (isc:iec, jsc:jec)   )
    allocate ( Ice % area          (isc:iec, jsc:jec) )
    allocate ( Ice % mi            (isc:iec, jsc:jec) )


    Ice % flux_sw_vis_dir = 0.
    Ice % flux_sw_vis_dif = 0.
    Ice % flux_sw_nir_dir = 0.
    Ice % flux_sw_nir_dif = 0.
    Ice % flux_lh         = 0. 
    Ice % lwdn            = 0.
    Ice % swdn            = 0.
    Ice % flux_u_top      = 0. 
    Ice % flux_v_top      = 0.
    Ice % flux_lh_top     = 0.
    Ice % flux_lw_top     = 0.
    Ice % flux_q_top      = 0.
    Ice % flux_sw_nir_dif_top = 0.
    Ice % flux_sw_nir_dir_top = 0.
    Ice % flux_sw_vis_dif_top = 0.
    Ice % flux_sw_vis_dir_top = 0.
    Ice % flux_t_top      = 0.
    Ice % fprec_top       = 0.
    Ice % lprec_top       = 0.
    Ice % pen             = 0.
    Ice % s_surf          = 0.
    Ice % u_surf          = 0.
    Ice % v_surf          = 0.
    Ice % sea_lev         = 0.
    Ice % u_ocn           = 0.
    Ice % v_ocn           = 0.
    Ice % u_ice           = 0.
    Ice % v_ice           = 0.
    Ice % sig11           = 0.
    Ice % sig12           = 0.
    Ice % sig22           = 0.
    Ice % h_snow          = 0.
    Ice % h_ice           = 0.
    Ice % t_ice1          =-5.
    Ice % t_ice2          =-5.
    Ice % mi              = 0.
    Ice % area            = cell_area * 4*PI*RADIUS*RADIUS

    do j = jsc, jec
       do i = isc, iec
          if( wett(i,j) > 0.5 ) then
             Ice % mask(i,j) = .true.
          else
             Ice % mask(i,j) = .false.
          end if
       enddo
    enddo

    Ice % Time              = Time
    Ice % Time_Init         = Time_Init
    Ice % Time_step_fast    = Time_step_fast
    Ice % Time_step_slow    = Time_step_slow
    Ice % avg_kount         = 0
    Ice % part_size         = 0.0
    Ice % part_size (:,:,1) = 1.0
    Ice % albedo            = 0.0
    Ice % albedo_vis_dir    = 0.0         
    Ice % albedo_nir_dir    = 0.0       
    Ice % albedo_vis_dif    = 0.0       
    Ice % albedo_nir_dif    = 0.0       
    Ice % rough_mom         = 1.0e-4
    Ice % rough_heat        = 1.0e-4
    Ice % rough_moist       = 1.0e-4
    Ice % t_surf            = 270.
    Ice % flux_u            = 0.0 
    Ice % flux_v            = 0.0
    Ice % flux_t            = 0.0 
    Ice % flux_q            = 0.0
    Ice % flux_lw           = 0.0
    Ice % flux_salt         = 0.0 
    Ice % lprec             = 0.0
    Ice % fprec             = 0.0
    Ice % p_surf            = 1.0e5
    Ice % runoff            = 0.0
    Ice % calving           = 0.0
    Ice % runoff_hflx       = 0.0
    Ice % calving_hflx      = 0.0
    Ice % frazil            = 0.0
    Ice % tmelt             = 0.0
    Ice % bmelt             = 0.0
    Ice % qflx_lim_ice      = 0.0
    Ice % qflx_res_ice      = 0.0
    Ice % part_size_uv(:,:,1) = 1.0
    do k=2,km
       call t_to_uv(Ice%part_size(:,:,k),Ice%part_size_uv(:,:,k))
       Ice%part_size_uv (:,:,1) = Ice%part_size_uv(:,:,1)-Ice%part_size_uv (:,:,k)
    end do

  end subroutine ice_model_init

  !#######################################################################

  subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_boundary
    type (ice_data_type),          intent(inout) :: Ice

    return

  end subroutine update_ice_model_slow_up

  subroutine update_ice_model_slow_dn ( Atmos_boundary, Land_boundary, Ice )
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
    type(land_ice_boundary_type),  intent(inout) :: Land_boundary
    type (ice_data_type),          intent(inout) :: Ice

    return

  end subroutine update_ice_model_slow_dn

  subroutine update_ice_model_fast ( Atmos_boundary, Ice )
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
    type (ice_data_type),          intent(inout) :: Ice

    return

  end subroutine update_ice_model_fast

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_bottom_to_ice_top - prepare surface state for atmosphere fast physics    !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_bottom_to_ice_top (Ice, t_surf_ice_bot, u_surf_ice_bot, v_surf_ice_bot, &
                                    frazil_ice_bot, Ocean_ice_boundary, s_surf_ice_bot, sea_lev_ice_bot )
    type (ice_data_type),                    intent(inout) :: Ice
    real, dimension(isc:iec,jsc:jec),           intent(in) :: t_surf_ice_bot, u_surf_ice_bot
    real, dimension(isc:iec,jsc:jec),           intent(in) :: v_surf_ice_bot, frazil_ice_bot
    type(ocean_ice_boundary_type),           intent(inout) :: Ocean_ice_boundary
    real, dimension(isc:iec,jsc:jec), intent(in), optional :: s_surf_ice_bot, sea_lev_ice_bot

    return

  end subroutine ice_bottom_to_ice_top

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! update_ice_model_slow - do ice dynamics, transport, and mass changes         !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine update_ice_model_slow (Ice, runoff, calving, &
                                         runoff_hflx, calving_hflx, p_surf)

    type (ice_data_type),             intent(inout)        :: Ice
    real, dimension(isc:iec,jsc:jec), intent(in), optional :: runoff, calving
    real, dimension(isc:iec,jsc:jec), intent(in), optional :: runoff_hflx, calving_hflx
    real, dimension(:,:,:),           intent(in), optional :: p_surf ! obsolete

    return

  end subroutine update_ice_model_slow

  subroutine ice_stock_pe(Ice, index, value)
    type(ice_data_type) :: Ice
    integer, intent(in) :: index
    real, intent(out)   :: value
  
    value = 0.0
    return
  end subroutine ice_stock_pe

  subroutine ice_model_restart(Ice, time_stamp)
    type (ice_data_type),     intent(inout), optional :: Ice
    character(len=*), intent(in), optional :: time_stamp

    return

  end subroutine ice_model_restart

  subroutine ice_model_end (Ice)
    type (ice_data_type), intent(inout) :: Ice

    deallocate(Ice % mask, Ice % ice_mask, Ice % t_surf, Ice % s_surf, Ice % sea_lev )
    deallocate(Ice % part_size, Ice % part_size_uv, Ice % u_surf, Ice % v_surf )
    deallocate(Ice % u_ocn, Ice % v_ocn ,  Ice % rough_mom, Ice % rough_heat )
    deallocate(Ice % rough_moist, Ice % albedo, Ice % flux_u_top, Ice % flux_v_top )
    deallocate(Ice % flux_t_top, Ice % flux_q_top, Ice % flux_lw_top )
    deallocate(Ice % flux_lh_top, Ice % lprec_top, Ice % fprec_top, Ice % flux_u )
    deallocate(Ice % flux_v, Ice % flux_t, Ice % flux_q, Ice % flux_lw )
    deallocate(Ice % flux_lh, Ice % lprec, Ice % fprec, Ice % p_surf, Ice % runoff ) 
    deallocate(Ice % calving, Ice % runoff_hflx, Ice % calving_hflx )
    deallocate(Ice % flux_salt)
    deallocate(Ice % lwdn)
    deallocate(Ice % swdn)
    deallocate(Ice % frazil )
    deallocate(Ice % bheat, Ice % u_ice, Ice % v_ice, Ice % sig11, Ice % sig22 )
    deallocate(Ice % sig12, Ice % tmelt, Ice % bmelt, Ice % pen, Ice % trn )
    deallocate(Ice % h_snow, Ice % h_ice, Ice % t_ice1, Ice % t_ice2  )
    deallocate(Ice % qflx_lim_ice, Ice % qflx_res_ice )
    deallocate(Ice % flux_sw_vis_dir, Ice % flux_sw_vis_dif )
    deallocate(Ice % flux_sw_nir_dir, Ice % flux_sw_nir_dif )

  end subroutine ice_model_end

  subroutine ocn_ice_bnd_type_chksum(id, timestep, Ocean_ice_boundary)
    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ocean_ice_boundary_type), intent(in) :: Ocean_ice_boundary
    return
  end subroutine ocn_ice_bnd_type_chksum

  subroutine atm_ice_bnd_type_chksum(id, timestep, Atmos_ice_boundary)
    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(atmos_ice_boundary_type), intent(in) :: Atmos_ice_boundary
    return
  end subroutine atm_ice_bnd_type_chksum

  subroutine lnd_ice_bnd_type_chksum(id, timestep, Land_ice_boundary)
    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_ice_boundary_type), intent(in) :: Land_ice_boundary
    return
  end subroutine lnd_ice_bnd_type_chksum

  subroutine ice_data_type_chksum(id, timestep, Ice_data)
    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_data_type), intent(in) :: Ice_data
    return
  end subroutine ice_data_type_chksum

end module ice_model_mod
