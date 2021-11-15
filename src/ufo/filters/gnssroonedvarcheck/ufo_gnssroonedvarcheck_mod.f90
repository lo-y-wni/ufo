! (C) Copyright 2017-2020 Met Office
! 
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> The main Fortran module for implementing the GNSS-RO onedvar check

module ufo_gnssroonedvarcheck_mod

use, intrinsic :: iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : fckit_log
use iso_c_binding
use kinds
use missing_values_mod
use obsspace_mod
use oops_variables_mod
use ufo_geovals_mod
use ufo_vars_mod
use ufo_gnssro_bendmetoffice_mod
use ufo_gnssroonedvarcheck_utils_mod, only: &
    Ops_RealSortQuick, deallocate_singleob, &
    allocate_singleob, allocate_singlebg, deallocate_singlebg, find_unique, &
    singlebg_type, singleob_type
use ufo_gnssroonedvarcheck_get_bmatrix_mod, only: bmatrix_type
use ufo_gnssroonedvarcheck_do1dvar_mod, only: Ops_GPSRO_Do1DVar_BA

implicit none
private
public :: ufo_gnssroonedvarcheck_create
public :: ufo_gnssroonedvarcheck_delete
public :: ufo_gnssroonedvarcheck_apply

type, public :: ufo_gnssroonedvarcheck
  character(len=800)        :: qcname            !< name of the filter
  type(c_ptr)               :: obsdb             !< pointer to the observation space
  type(fckit_configuration) :: conf              !< contents of the yaml file
  integer(c_int)            :: onedvarflag       !< flag uased by the qc manager for a 1D-var check
  logical                   :: capsupersat       !< Whether to remove super-saturation (wrt ice?)
  real(kind_real)           :: cost_funct_test   !< Threshold value for the cost function convergence test
  real(kind_real)           :: Delta_ct2         !< Threshold used in calculating convergence
  real(kind_real)           :: Delta_factor      !< Threshold used in calculating convergence
  integer                   :: n_iteration_test  !< Maximum number of iterations in the 1DVar
  real(kind_real)           :: OB_test           !< Threshold for the O-B throughout the profile
  real(kind_real)           :: y_test            !< Threshold on distance between observed and solution bending angles
  character(len=800)        :: bmatrix_filename  !< Location of the B-matrix file
  logical                   :: pseudo_ops        !< Whether to use pseudo levels in forward operator
  logical                   :: vert_interp_ops   !< Whether to use ln(p) or exner in vertical interpolation
  real(kind_real)           :: min_temp_grad     !< The minimum vertical temperature gradient allowed
end type ufo_gnssroonedvarcheck

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup the main GNSS-RO onedvar object in Fortran
!!
!! \details Makes a call to the main setup routine.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_gnssroonedvarcheck_create(self, obsspace, bmatrix_filename, &
                                         capsupersat, cost_funct_test, &
                                         Delta_ct2, Delta_factor, min_temp_grad, &
                                         n_iteration_test, OB_test, pseudo_ops, &
                                         vert_interp_ops, y_test, onedvarflag)

  implicit none
  type(ufo_gnssroonedvarcheck), intent(inout) :: self              !< gnssroonedvarcheck main object
  type(c_ptr), value, intent(in)              :: obsspace          !< observation database pointer
  character(len=*), intent(in)                :: bmatrix_filename  !< Location of the B-matrix file
  logical(c_bool), intent(in)                 :: capsupersat       !< Whether to remove super-saturation (wrt ice?)
  real(c_float), intent(in)                   :: cost_funct_test   !< Threshold value for the cost function convergence test
  real(c_float), intent(in)                   :: Delta_ct2         !< Threshold used in calculating convergence
  real(c_float), intent(in)                   :: Delta_factor      !< Threshold used in calculating convergence
  real(c_float), intent(in)                   :: min_temp_grad     !< The minimum vertical temperature gradient allowed
  integer(c_int), intent(in)                  :: n_iteration_test  !< Maximum number of iterations in the 1DVar
  real(c_float), intent(in)                   :: OB_test           !< Threshold for the O-B throughout the profile
  logical(c_bool), intent(in)                 :: pseudo_ops        !< Whether to use pseudo levels in forward operator
  logical(c_bool), intent(in)                 :: vert_interp_ops   !< Whether to use ln(p) or exner in vertical interpolation
  real(c_float), intent(in)                   :: y_test            !< Threshold on distance between observed and solution bending angles
  integer(c_int), intent(in)                  :: onedvarflag       !< flag for qc manager

  character(len=800) :: message

  self % obsdb = obsspace
  self % onedvarflag = onedvarflag

  self % bmatrix_filename = bmatrix_filename
  self % capsupersat = capsupersat
  self % cost_funct_test = cost_funct_test
  self % Delta_ct2 = Delta_ct2
  self % Delta_factor = Delta_factor
  self % min_temp_grad = min_temp_grad
  self % n_iteration_test = n_iteration_test
  self % OB_test = OB_test
  self % pseudo_ops = pseudo_ops
  self % vert_interp_ops = vert_interp_ops
  self % y_test = y_test

  write(message, '(A)') 'GNSS-RO 1D-Var check: input parameters are:'
  call fckit_log % debug(message)
  write(message, '(2A)') 'bmatrix_filename = ', bmatrix_filename
  call fckit_log % debug(message)
  write(message, '(A,L1)') 'capsupersat = ', capsupersat
  call fckit_log % debug(message)
  write(message, '(A,F16.8)') 'cost_funct_test = ', cost_funct_test
  call fckit_log % debug(message)
  write(message, '(A,F16.8)') 'Delta_ct2 = ', Delta_ct2
  call fckit_log % debug(message)
  write(message, '(A,F16.8)') 'Delta_factor = ', Delta_factor
  call fckit_log % debug(message)
  write(message, '(A,F16.8)') 'min_temp_grad = ', min_temp_grad
  call fckit_log % debug(message)
  write(message, '(A,I7)') 'n_iteration_test = ', n_iteration_test
  call fckit_log % debug(message)
  write(message, '(A,F16.8)') 'OB_test = ', OB_test
  call fckit_log % debug(message)
  write(message, '(A,L1)') 'pseudo_ops = ', pseudo_ops
  call fckit_log % debug(message)
  write(message, '(A,L1)') 'vert_interp_ops = ', vert_interp_ops
  call fckit_log % debug(message)
  write(message, '(A,F16.8)') 'y_test = ', y_test
  call fckit_log % debug(message)

end subroutine ufo_gnssroonedvarcheck_create

! ------------------------------------------------------------------------------
!> Delete the main GNSS-RO onedvar object in Fortran
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_gnssroonedvarcheck_delete(self)

  implicit none
  type(ufo_gnssroonedvarcheck), intent(inout) :: self !< gnssroonedvarcheck main object

end subroutine ufo_gnssroonedvarcheck_delete

! ------------------------------------------------------------------------------
!> The main routine that applys the GNSS-RO onedvar filter
!!
!! \details Heritage : 
!!
!! This routine is called from the c++ apply method.  The filter performs 
!! a 1D-Var minimization
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_gnssroonedvarcheck_apply(self, geovals, apply)

  implicit none

  ! Subroutine arguments
  type(ufo_gnssroonedvarcheck), intent(inout) :: self    !< gnssroonedvarcheck main object
  type(ufo_geovals), intent(in)              :: geovals  !< model values at observation space
  logical, intent(in)                        :: apply(:) !< qc manager flags

  ! Local parameters
  logical, parameter :: verboseOutput = .FALSE.          ! Whether to output extra debugging information

  ! Local variables
  integer :: nobs                                        ! Number of observations to be processed
  type(ufo_geoval), pointer          :: q                ! Model background values of specific humidity
  type(ufo_geoval), pointer          :: prs              ! Model background values of air pressure
  type(ufo_geoval), pointer          :: theta_heights    ! Model heights of levels containing specific humidity
  type(ufo_geoval), pointer          :: rho_heights      ! Model heights of levels containing air pressure
  type(singlebg_type)                :: Back             ! Model background fields
  type(singleob_type)                :: Ob               ! The profile of observations
  type(bmatrix_type)                 :: b_matrix         ! Background-error covariance matrix
  real(kind_real), allocatable       :: obsLat(:)             ! Latitude of the observation
  real(kind_real), allocatable       :: obsLon(:)             ! Longitude of the observation
  real(kind_real), allocatable       :: impact_param(:)       ! Impact parameter of the observation
  real(kind_real), allocatable       :: radius_curv(:)        ! Earth's radius of curvature at the observation tangent point
  real(kind_real), allocatable       :: undulation(:)         ! Undulation - height of the geoid above the ellipsoid
  real(kind_real), allocatable       :: obs_bending_angle(:)  ! Observed bending angle
  real(kind_real), allocatable       :: obs_err(:)            ! Observation error, taken from a previous filter
  integer, allocatable               :: qc_flags(:)           ! QC flags to be updated
  integer, allocatable               :: obsSatid(:)           ! Satellite identifier for each observation
  integer, allocatable               :: obsOrigC(:)           ! Originating centre for each observation
  integer(c_size_t), allocatable     :: record_number(:)      ! Number used to identify unique profiles in the data
  real(kind_real), allocatable       :: sort_key(:)           ! Key for the sorting (based on record number and impact parameter)
  integer, allocatable               :: index_vals(:)         ! Indices of sorted observation
  integer, allocatable               :: unique(:)             ! Set of unique profile numbers
  integer                            :: start_point           ! Starting index of the current profile
  integer                            :: current_point         ! Ending index of the current profile
  integer                            :: iprofile              ! Loop variable, profile number
  integer                            :: nobs_profile          ! Number of observations in the profile
  character(len=800)                 :: Message               ! Message to be output
  logical                            :: BAerr                 ! Has there been an error in the bending angle calculation?
  integer                            :: iband                 ! Selected latitude band of the B-matrix
  integer                            :: iseason               ! Selected season of the B-matrix
  integer                            :: ipoint                ! Loop variable, observation point
  real(kind_real)                    :: dfs                   ! Degrees of freedom for signal in profile
  real(kind_real)                    :: O_Bdiff               ! Average RMS(O-B) for profile
  real(kind_real), allocatable       :: Tb(:)                 ! Calculated background temperature (derived from p,q)
  real(kind_real), allocatable       :: Ts(:)                 ! 1DVar solution temperature

  ! Get the obs-space information
  nobs = obsspace_get_nlocs(self % obsdb)
  allocate(obsLon(nobs))
  allocate(obsLat(nobs))
  allocate(impact_param(nobs))
  allocate(radius_curv(nobs))
  allocate(undulation(nobs))
  allocate(obs_bending_angle(nobs))
  allocate(record_number(nobs))
  allocate(obsSatid(nobs))
  allocate(obsOrigC(nobs))
  allocate(qc_flags(nobs))
  allocate(obs_err(nobs))

  call obsspace_get_db(self % obsdb, "MetaData", "longitude", obsLon)
  call obsspace_get_db(self % obsdb, "MetaData", "latitude", obsLat)
  call obsspace_get_db(self % obsdb, "MetaData", "impact_parameter", impact_param)
  call obsspace_get_db(self % obsdb, "MetaData", "earth_radius_of_curvature", radius_curv)
  call obsspace_get_db(self % obsdb, "MetaData", "geoid_height_above_reference_ellipsoid", undulation)
  call obsspace_get_db(self % obsdb, "ObsValue", "bending_angle", obs_bending_angle)
  call obsspace_get_recnum(self % obsdb, record_number)
  call obsspace_get_db(self % obsdb, "MetaData", "occulting_sat_id", obsSatid)
  call obsspace_get_db(self % obsdb, "MetaData", "originating_center", obsOrigC)
  call obsspace_get_db(self % obsdb, "FortranQC", "bending_angle", qc_flags)
  call obsspace_get_db(self % obsdb, "FortranERR", "bending_angle", obs_err)

  ! get variables from geovals
  call ufo_geovals_get_var(geovals, var_q, q)               ! specific humidity
  call ufo_geovals_get_var(geovals, var_prsi, prs)          ! pressure
  call ufo_geovals_get_var(geovals, var_z, theta_heights)   ! Geopotential height of the normal model levels
  call ufo_geovals_get_var(geovals, var_zi, rho_heights)    ! Geopotential height of the pressure levels

  ! Read in the B-matrix and background profiles
  call b_matrix % get(self % bmatrix_filename, prs % nval, q % nval)
  call allocate_singlebg(Back, prs % nval, q % nval)
  allocate(Tb(q % nval), Ts(q % nval))

  ! Read through the record numbers in order to find a profile of observations
  ! Each profile shares the same record number
  allocate(sort_key(nobs))
  sort_key = record_number + (impact_param / MAXVAL(impact_param))
  call Ops_RealSortQuick(sort_key, index_vals)
  call find_unique(record_number, unique)

  ! For every profile that we have found, perform a 1DVar minimisation
  current_point = 1
  do iprofile = 1, size(unique)
    start_point = current_point
    WRITE (Message, '(A,I0)') 'ObNumber ', iprofile
    call fckit_log % info(Message)
    WRITE (Message, '(A,F12.2)') 'Latitude ', obsLat(index_vals(start_point))
    call fckit_log % info(Message)
    WRITE (Message, '(A,F12.2)') 'Longitude ', obsLon(index_vals(start_point))
    call fckit_log % info(Message)
    WRITE (Message, '(A,I0)') 'Processing centre ', obsOrigC(index_vals(start_point))
    call fckit_log % info(Message)
    WRITE (Message, '(A,I0)') 'Sat ID ', obsSatid(index_vals(start_point))
    call fckit_log % info(Message)

    ! Work out which observations belong to the current profile
    do current_point = start_point, nobs
      if (unique(iprofile) /= record_number(index_vals(current_point))) exit
    end do

    ! Load the geovals into the background structure
    Back % za(:) = rho_heights % vals(:, index_vals(start_point))
    Back % zb(:) = theta_heights % vals(:, index_vals(start_point))
    Back % p(:) = prs % vals(:, index_vals(start_point))
    Back % q(:) = q % vals(:, index_vals(start_point))

    ! Allocate the observations structure
    nobs_profile = current_point - start_point
    call allocate_singleob(Ob, nobs_profile, prs % nval, q % nval)

    ! Load the observations information into the obsevations structure
    Ob % id = record_number(index_vals(start_point))
    Ob % latitude = obsLat(index_vals(start_point))
    Ob % longitude = obsLon(index_vals(start_point))
    Ob % niter = 0  ! We haven't yet run 1DVar
    Ob % jcost = missing_value(Ob % jcost)
    Ob % bendingangle(:) % value = obs_bending_angle(index_vals(start_point:current_point-1))
    Ob % bendingangle(:) % oberr = obs_err(index_vals(start_point:current_point-1))
    Ob % impactparam(:) % value = impact_param(index_vals(start_point:current_point-1))
    Ob % qc_flags(:) = qc_flags(index_vals(start_point:current_point-1))
    Ob % ro_rad_curv % value = radius_curv(index_vals(start_point))
    Ob % ro_geoid_und % value = undulation(index_vals(start_point))

    ! Choose the latitude band and season of the B-matrix information
    iseason = 1    ! Temporary -only one season at present!
    iband = 1
    DO
      IF (b_matrix % band_up_lim(iband) > Ob % latitude .OR. &
          iband == b_matrix % nband) EXIT
      iband = iband + 1
    END DO

!   Call the code to set up the 1D-Var calculation
    call Ops_GPSRO_Do1DVar_BA(prs % nval,              &   ! Number of pressure levels
                              q % nval,                &   ! Number of specific humidity levels
                              b_matrix % inverse(iseason,iband,:,:), &  ! Inverse of the b-matrix
                              b_matrix % sigma(iseason,iband,:), &      ! Standard deviations of the b-matrix
                              Back,                    &   ! Structure containing the model background information
                              Ob,                      &   ! Structure containing the observation information
                              self % pseudo_ops,       &   ! Whether to use pseudo-levels in calculation
                              self % vert_interp_ops,  &   ! Whether to interpolate using ln(p) or exner
                              self % min_temp_grad,    &   ! Minimum vertical temperature gradient allowed
                              self % cost_funct_test,  &   ! Threshold value for the cost function convergence test
                              self % y_test,           &   ! Threshold value for the yobs-ysol tes
                              self % n_iteration_test, &   ! Maximum number of iterations
                              self % Delta_factor,     &   ! Delta
                              self % Delta_ct2,        &   ! Delta observations
                              self % OB_test,          &   ! Threshold value for the O-B test
                              self % capsupersat,      &   ! Whether to remove super-saturation
                              BAerr,                   &   ! Whether there are errors in the bending angle calculation
                              Tb,                      &   ! Calculated background temperature
                              Ts,                      &   ! 1DVar solution temperature
                              O_Bdiff,                 &   ! Difference between observations and background for profile
                              DFS)                         ! Estimated degrees of freedom for signal

    ! Flag bad profiles
    do ipoint = 0, nobs_profile-1
      if (qc_flags(start_point + ipoint) > 0) then
        ! Do nothing, since the data are already flagged
      else if (Ob % bendingangle(ipoint+1) % PGEFinal > 0.5) then
        qc_flags(start_point + ipoint) = self % onedvarflag
        Ob % qc_flags(ipoint+1) = self % onedvarflag
      end if
    end do

    write(Message,'(A,2I5,2F10.3,I5,E16.8)') 'Profile stats: ', obsSatid(index_vals(start_point)), &
        obsOrigC(index_vals(start_point)), Ob % latitude, Ob % longitude, &
        Ob % niter, Ob % jcost
    call fckit_log % debug(Message)
    
    if (verboseOutput) then
      do ipoint = 0, nobs_profile-1, 20
          write(Message,'(20I5)') qc_flags(index_vals(start_point+ipoint: &
                                                      min(start_point+ipoint+19, current_point-1)))
          call fckit_log % debug(Message)
      end do
      do ipoint = 0, nobs_profile-1, 10
          write(Message,'(10E16.5)') obs_bending_angle(index_vals(start_point+ipoint: &
                                                               min(start_point+ipoint+9, current_point-1)))
          call fckit_log % debug(Message)
      end do
    end if

    call deallocate_singleob(Ob)
  end do
  call obsspace_put_db(self % obsdb, "FortranQC", "bending_angle", qc_flags)

end subroutine ufo_gnssroonedvarcheck_apply

end module ufo_gnssroonedvarcheck_mod
