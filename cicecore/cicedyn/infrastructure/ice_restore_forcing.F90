!=======================================================================
!
! Time-varying prescribed sea-ice boundary data for MITgcm-like
! restoring. Supported file modes use Fortran-ordered dimensions
! compatible with the CICE NetCDF read path.
!
! daily_netcdf requires one record per day beginning on Jan 1 of
! fyear_init when restore_ice_cycle_year is false, or a full 365-day
! climatology when restore_ice_cycle_year is true.
!
! monthly_netcdf requires one record per month beginning on Jan of
! fyear_init when restore_ice_cycle_year is false, or a 12-record
! climatology when restore_ice_cycle_year is true.
!
! Required variables for all restore formats:
!   aicen(ni, nj, ncat, time)
!   vicen(ni, nj, ncat, time)
!   vsnon(ni, nj, ncat, time)
! Supported thermo representations:
!   full_trcrn: optional trcrn(ni, nj, ntrcr, ncat, time)
!   bc_fields_v1: Tsfc(ni, nj, ncat, time),
!                 Tinz(ni, nj, nilyr, ncat, time),
!                 Sinz(ni, nj, nilyr, ncat, time)
! When neither thermo representation is present, the restoring only
! targets structural state and leaves tracer fields unchanged.
!

#ifdef ncdf
#define USE_NETCDF
#endif

            module ice_restore_forcing

            use ice_kinds_mod
            use ice_blocks, only: nx_block, ny_block, block, get_block, nblocks_x, nblocks_y
            use ice_broadcast, only: broadcast_scalar
            use ice_communicate, only: my_task, master_task
            use ice_constants, only: c0, c1, field_loc_center, field_type_scalar
            use ice_domain, only: distrb_info, nblocks, blocks_ice, ew_boundary_type, ns_boundary_type
            use ice_domain_size, only: ncat, max_blocks, nx_global, ny_global, nilyr
            use ice_exit, only: abort_ice
            use ice_fileunits, only: nu_diag
            use ice_forcing, only: restore_ice_data_type, restore_ice_data_file, &
                                                        restore_ice_cycle_year, restore_ice_use_west, &
                                                        restore_ice_use_east, restore_ice_use_south, &
                                                        restore_ice_use_north, fyear_init
            use ice_calendar, only: msec, mday, mmonth, myear, yday, days_per_year, &
                                                            compute_days_between, daymo
            use ice_gather_scatter, only: scatter_global
            use ice_grid, only: hm
            use ice_read_write, only: ice_open_nc, ice_close_nc, ice_check_nc
            use ice_state, only: aicen, vicen, vsnon, trcrn
            use icepack_intfc, only: icepack_query_parameters

#ifdef USE_NETCDF
            use netcdf, only: nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension, &
                                                nf90_get_var, NF90_MAX_VAR_DIMS
#endif

            implicit none
            private
            public :: restore_forcing_init, restore_forcing_update, &
                               restore_forcing_finalize, &
                               restore_forcing_is_active, restore_forcing_has_trcrn_data, &
                               restore_forcing_uses_bc_fields

            logical (kind=log_kind) :: &
                  restore_forcing_initialized = .false., &
                  restore_forcing_active      = .false., &
                  restore_forcing_has_trcrn   = .false., &
                  restore_forcing_has_bc_thermo = .false.

            integer (kind=int_kind) :: &
                  restore_forcing_ntrcr    = 0, &
                  restore_forcing_fid      = -1, &
                  restore_forcing_nrecords = 0, &
                  restore_forcing_varid_aicen = -1, &
                  restore_forcing_varid_vicen = -1, &
                  restore_forcing_varid_vsnon = -1, &
                  restore_forcing_varid_Tsfc  = -1, &
                  restore_forcing_varid_Tinz  = -1, &
                  restore_forcing_varid_Sinz  = -1, &
                  restore_forcing_varid_trcrn = -1, &
                  restore_forcing_loaded_record_1 = -1, &
                  restore_forcing_loaded_record_2 = -1

            logical (kind=log_kind) :: &
                  restore_forcing_cycle_year = .true., &
                  restore_forcing_use_west   = .true., &
                  restore_forcing_use_east   = .true., &
                  restore_forcing_use_south  = .true., &
                  restore_forcing_use_north  = .true.

            character (char_len) :: &
                  restore_forcing_data_type = 'legacy'

            character (char_len_long) :: &
                  restore_forcing_data_file = 'unknown_restore_ice_file'

            real (kind=dbl_kind), dimension(:,:), allocatable :: &
                  restore_forcing_work_g

            real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable :: &
                  aicen_cache, vicen_cache, vsnon_cache, Tsfc_cache

            real (kind=dbl_kind), dimension(:,:,:,:,:,:), allocatable :: &
                  trcrn_cache, Tinz_cache, Sinz_cache

            integer (kind=int_kind), dimension(:), allocatable :: &
                  restore_forcing_ilo, restore_forcing_ihi, &
                  restore_forcing_jlo, restore_forcing_jhi, &
                  restore_forcing_east_bc, restore_forcing_north_bc

            logical (kind=log_kind), dimension(:), allocatable :: &
                  restore_forcing_west_block, restore_forcing_east_block, &
                  restore_forcing_south_block, restore_forcing_north_block

!=======================================================================

            contains

!=======================================================================

            subroutine restore_forcing_init(ntrcr)

            integer (kind=int_kind), intent(in) :: &
                  ntrcr

            character(len=*), parameter :: subname = '(restore_forcing_init)'

            restore_forcing_initialized = .true.
            restore_forcing_ntrcr = ntrcr

            restore_forcing_data_type = trim(restore_ice_data_type)
            restore_forcing_data_file = trim(restore_ice_data_file)
            restore_forcing_cycle_year = restore_ice_cycle_year
            restore_forcing_use_west = restore_ice_use_west
            restore_forcing_use_east = restore_ice_use_east
            restore_forcing_use_south = restore_ice_use_south
            restore_forcing_use_north = restore_ice_use_north

            restore_forcing_active = .false.
            restore_forcing_has_trcrn = .false.
            restore_forcing_has_bc_thermo = .false.
            restore_forcing_fid = -1
            restore_forcing_varid_aicen = -1
            restore_forcing_varid_vicen = -1
            restore_forcing_varid_vsnon = -1
            restore_forcing_varid_Tsfc = -1
            restore_forcing_varid_Tinz = -1
            restore_forcing_varid_Sinz = -1
            restore_forcing_varid_trcrn = -1
            restore_forcing_loaded_record_1 = -1
            restore_forcing_loaded_record_2 = -1

            if (trim(restore_forcing_data_type) == 'legacy') return

            if (my_task == master_task) then
                  write(nu_diag,*) 'CICE restore forcing:'
                  write(nu_diag,*) '  data_type  = ', trim(restore_forcing_data_type)
                  write(nu_diag,*) '  data_file  = ', trim(restore_forcing_data_file)
                  write(nu_diag,*) '  cycle_year = ', restore_forcing_cycle_year
                  write(nu_diag,*) '  use_west   = ', restore_forcing_use_west
                  write(nu_diag,*) '  use_east   = ', restore_forcing_use_east
                  write(nu_diag,*) '  use_south  = ', restore_forcing_use_south
                  write(nu_diag,*) '  use_north  = ', restore_forcing_use_north
            endif

                  if (trim(restore_forcing_data_type) /= 'daily_netcdf' .and. &
                        trim(restore_forcing_data_type) /= 'monthly_netcdf') then
                  call abort_ice(error_message=subname//' unsupported restore_ice_data_type='// &
                                                trim(restore_forcing_data_type), file=__FILE__, line=__LINE__)
            endif

            if (trim(restore_forcing_data_file) == 'unknown_restore_ice_file') then
                  call abort_ice(error_message=subname//' restore_ice_data_file must be set', &
                                                file=__FILE__, line=__LINE__)
            endif

#ifndef USE_NETCDF
            call abort_ice(error_message=subname//' daily_netcdf requires USE_NETCDF', &
                                          file=__FILE__, line=__LINE__)
#else
            call restore_forcing_cache_block_metadata()
            call restore_forcing_allocate_cache()
            call ice_open_nc(trim(restore_forcing_data_file), restore_forcing_fid)
            call restore_forcing_validate_file()
            restore_forcing_active = .true.
#endif

            end subroutine restore_forcing_init

!=======================================================================

            subroutine restore_forcing_finalize()

#ifdef USE_NETCDF
            if (my_task == master_task) then
                  if (restore_forcing_fid >= 0) call ice_close_nc(restore_forcing_fid)
            endif
#endif

            if (allocated(restore_forcing_work_g)) deallocate(restore_forcing_work_g)
            if (allocated(aicen_cache)) deallocate(aicen_cache)
            if (allocated(vicen_cache)) deallocate(vicen_cache)
            if (allocated(vsnon_cache)) deallocate(vsnon_cache)
            if (allocated(Tsfc_cache)) deallocate(Tsfc_cache)
            if (allocated(trcrn_cache)) deallocate(trcrn_cache)
            if (allocated(Tinz_cache)) deallocate(Tinz_cache)
            if (allocated(Sinz_cache)) deallocate(Sinz_cache)
            if (allocated(restore_forcing_ilo)) deallocate(restore_forcing_ilo)
            if (allocated(restore_forcing_ihi)) deallocate(restore_forcing_ihi)
            if (allocated(restore_forcing_jlo)) deallocate(restore_forcing_jlo)
            if (allocated(restore_forcing_jhi)) deallocate(restore_forcing_jhi)
            if (allocated(restore_forcing_east_bc)) deallocate(restore_forcing_east_bc)
            if (allocated(restore_forcing_north_bc)) deallocate(restore_forcing_north_bc)
            if (allocated(restore_forcing_west_block)) deallocate(restore_forcing_west_block)
            if (allocated(restore_forcing_east_block)) deallocate(restore_forcing_east_block)
            if (allocated(restore_forcing_south_block)) deallocate(restore_forcing_south_block)
            if (allocated(restore_forcing_north_block)) deallocate(restore_forcing_north_block)

            restore_forcing_initialized = .false.
            restore_forcing_active = .false.
            restore_forcing_has_trcrn = .false.
            restore_forcing_has_bc_thermo = .false.
            restore_forcing_fid = -1
            restore_forcing_nrecords = 0
            restore_forcing_loaded_record_1 = -1
            restore_forcing_loaded_record_2 = -1

            end subroutine restore_forcing_finalize

!=======================================================================

            logical (kind=log_kind) function restore_forcing_is_active()

            restore_forcing_is_active = restore_forcing_active

            end function restore_forcing_is_active

!=======================================================================

            logical (kind=log_kind) function restore_forcing_has_trcrn_data()

            restore_forcing_has_trcrn_data = restore_forcing_has_trcrn

            end function restore_forcing_has_trcrn_data

!=======================================================================

            logical (kind=log_kind) function restore_forcing_uses_bc_fields()

            restore_forcing_uses_bc_fields = restore_forcing_has_bc_thermo

            end function restore_forcing_uses_bc_fields

!=======================================================================

            subroutine restore_forcing_update(aicen_rest, vicen_rest, vsnon_rest, trcrn_rest, &
                                                                            Tsfc_rest, Tinz_rest, Sinz_rest)

            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: &
                  aicen_rest, vicen_rest, vsnon_rest

            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: &
                  trcrn_rest

            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout), optional :: &
                  Tsfc_rest

            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout), optional :: &
                  Tinz_rest, Sinz_rest

            integer (kind=int_kind) :: recnum_1, recnum_2
            real (kind=dbl_kind) :: weight_2

            if (.not. restore_forcing_initialized) return
            if (.not. restore_forcing_active) return

            call restore_forcing_select_records(recnum_1, recnum_2, weight_2)
            call restore_forcing_ensure_records(recnum_1, recnum_2)

            ! The restoring code only consumes the physical boundary lines and
            ! their outward ghost/padding cells.  Preparing interior targets here
            ! caused large 4-D/5-D array copies on every CICE timestep.
            call restore_forcing_interpolate_boundaries(weight_2, aicen_rest, vicen_rest, &
                                                         vsnon_rest, trcrn_rest, &
                                                         Tsfc_rest, Tinz_rest, Sinz_rest)

            call restore_forcing_extend_open_edges(aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
            call restore_forcing_disable_selected_edges(aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
            call restore_forcing_apply_land_mask(aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
            if (restore_forcing_has_bc_thermo .and. present(Tsfc_rest) .and. present(Tinz_rest) .and. present(Sinz_rest)) then
                  call restore_forcing_extend_open_edges_bc(Tsfc_rest, Tinz_rest, Sinz_rest)
                  call restore_forcing_apply_land_mask_bc(Tsfc_rest, Tinz_rest, Sinz_rest)
            endif

            end subroutine restore_forcing_update

!=======================================================================

            subroutine restore_forcing_cache_block_metadata()

            integer (kind=int_kind) :: i, j, iblk, npad
            type (block) :: this_block

            if (.not. allocated(restore_forcing_ilo)) then
                  allocate(restore_forcing_ilo(max_blocks), restore_forcing_ihi(max_blocks), &
                           restore_forcing_jlo(max_blocks), restore_forcing_jhi(max_blocks), &
                           restore_forcing_east_bc(max_blocks), restore_forcing_north_bc(max_blocks), &
                           restore_forcing_west_block(max_blocks), restore_forcing_east_block(max_blocks), &
                           restore_forcing_south_block(max_blocks), restore_forcing_north_block(max_blocks))
            endif

            restore_forcing_ilo = 0
            restore_forcing_ihi = 0
            restore_forcing_jlo = 0
            restore_forcing_jhi = 0
            restore_forcing_east_bc = nx_block
            restore_forcing_north_bc = ny_block
            restore_forcing_west_block = .false.
            restore_forcing_east_block = .false.
            restore_forcing_south_block = .false.
            restore_forcing_north_block = .false.

            do iblk = 1, nblocks
                  this_block = get_block(blocks_ice(iblk), iblk)
                  restore_forcing_ilo(iblk) = this_block%ilo
                  restore_forcing_ihi(iblk) = this_block%ihi
                  restore_forcing_jlo(iblk) = this_block%jlo
                  restore_forcing_jhi(iblk) = this_block%jhi

                  restore_forcing_west_block(iblk) = this_block%iblock == 1 .and. &
                                                     trim(ew_boundary_type) /= 'cyclic'
                  restore_forcing_east_block(iblk) = this_block%iblock == nblocks_x .and. &
                                                     trim(ew_boundary_type) /= 'cyclic'
                  restore_forcing_south_block(iblk) = this_block%jblock == 1 .and. &
                                                      trim(ns_boundary_type) /= 'cyclic'
                  restore_forcing_north_block(iblk) = this_block%jblock == nblocks_y .and. &
                                                      trim(ns_boundary_type) /= 'cyclic' .and. &
                                                      trim(ns_boundary_type) /= 'tripole' .and. &
                                                      trim(ns_boundary_type) /= 'tripoleT'

                  if (restore_forcing_east_block(iblk)) then
                        do i = nx_block, 1, -1
                              npad = 0
                              if (this_block%i_glob(i) == 0) then
                                    do j = 1, ny_block
                                          npad = npad + this_block%j_glob(j)
                                    enddo
                              endif
                              if (npad /= 0) restore_forcing_east_bc(iblk) = &
                                   restore_forcing_east_bc(iblk) - 1
                        enddo
                  endif

                  if (restore_forcing_north_block(iblk)) then
                        do j = ny_block, 1, -1
                              npad = 0
                              if (this_block%j_glob(j) == 0) then
                                    do i = 1, nx_block
                                          npad = npad + this_block%i_glob(i)
                                    enddo
                              endif
                              if (npad /= 0) restore_forcing_north_bc(iblk) = &
                                   restore_forcing_north_bc(iblk) - 1
                        enddo
                  endif
            enddo

            end subroutine restore_forcing_cache_block_metadata

!=======================================================================

            subroutine restore_forcing_interpolate_boundaries(weight_2, aicen_rest, vicen_rest, &
                                                               vsnon_rest, trcrn_rest, &
                                                               Tsfc_rest, Tinz_rest, Sinz_rest)

            real (kind=dbl_kind), intent(in) :: weight_2
            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: &
                  aicen_rest, vicen_rest, vsnon_rest
            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: trcrn_rest
            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout), optional :: Tsfc_rest
            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout), optional :: Tinz_rest, Sinz_rest

            integer (kind=int_kind) :: i, j, iblk

            do iblk = 1, nblocks
                  if (restore_forcing_west_block(iblk)) then
                        i = restore_forcing_ilo(iblk)
                        do j = 1, ny_block
                              call restore_forcing_interpolate_cell(i, j, iblk, weight_2, &
                                   aicen_rest, vicen_rest, vsnon_rest, trcrn_rest, &
                                   Tsfc_rest, Tinz_rest, Sinz_rest)
                        enddo
                  endif

                  if (restore_forcing_east_block(iblk)) then
                        i = restore_forcing_ihi(iblk)
                        do j = 1, ny_block
                              call restore_forcing_interpolate_cell(i, j, iblk, weight_2, &
                                   aicen_rest, vicen_rest, vsnon_rest, trcrn_rest, &
                                   Tsfc_rest, Tinz_rest, Sinz_rest)
                        enddo
                  endif

                  if (restore_forcing_south_block(iblk)) then
                        j = restore_forcing_jlo(iblk)
                        do i = 1, nx_block
                              call restore_forcing_interpolate_cell(i, j, iblk, weight_2, &
                                   aicen_rest, vicen_rest, vsnon_rest, trcrn_rest, &
                                   Tsfc_rest, Tinz_rest, Sinz_rest)
                        enddo
                  endif

                  if (restore_forcing_north_block(iblk)) then
                        j = restore_forcing_jhi(iblk)
                        do i = 1, nx_block
                              call restore_forcing_interpolate_cell(i, j, iblk, weight_2, &
                                   aicen_rest, vicen_rest, vsnon_rest, trcrn_rest, &
                                   Tsfc_rest, Tinz_rest, Sinz_rest)
                        enddo
                  endif
            enddo

            end subroutine restore_forcing_interpolate_boundaries

!=======================================================================

            subroutine restore_forcing_interpolate_cell(i, j, iblk, weight_2, &
                                                        aicen_rest, vicen_rest, vsnon_rest, &
                                                        trcrn_rest, Tsfc_rest, Tinz_rest, Sinz_rest)

            integer (kind=int_kind), intent(in) :: i, j, iblk
            real (kind=dbl_kind), intent(in) :: weight_2
            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: &
                  aicen_rest, vicen_rest, vsnon_rest
            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: trcrn_rest
            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout), optional :: Tsfc_rest
            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout), optional :: Tinz_rest, Sinz_rest

            integer (kind=int_kind) :: k, n, nt
            real (kind=dbl_kind) :: weight_1

            weight_1 = c1 - weight_2

            do n = 1, ncat
                  aicen_rest(i,j,n,iblk) = weight_1 * aicen_cache(i,j,n,iblk,1) + &
                                           weight_2 * aicen_cache(i,j,n,iblk,2)
                  vicen_rest(i,j,n,iblk) = weight_1 * vicen_cache(i,j,n,iblk,1) + &
                                           weight_2 * vicen_cache(i,j,n,iblk,2)
                  vsnon_rest(i,j,n,iblk) = weight_1 * vsnon_cache(i,j,n,iblk,1) + &
                                           weight_2 * vsnon_cache(i,j,n,iblk,2)

                  if (restore_forcing_has_trcrn) then
                        do nt = 1, restore_forcing_ntrcr
                              trcrn_rest(i,j,nt,n,iblk) = weight_1 * trcrn_cache(i,j,nt,n,iblk,1) + &
                                                         weight_2 * trcrn_cache(i,j,nt,n,iblk,2)
                        enddo
                  else
                        do nt = 1, restore_forcing_ntrcr
                              trcrn_rest(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk)
                        enddo
                  endif

                  if (restore_forcing_has_bc_thermo) then
                        if (present(Tsfc_rest)) then
                              Tsfc_rest(i,j,n,iblk) = weight_1 * Tsfc_cache(i,j,n,iblk,1) + &
                                                     weight_2 * Tsfc_cache(i,j,n,iblk,2)
                        endif
                        if (present(Tinz_rest)) then
                              do k = 1, nilyr
                                    Tinz_rest(i,j,k,n,iblk) = weight_1 * Tinz_cache(i,j,k,n,iblk,1) + &
                                                             weight_2 * Tinz_cache(i,j,k,n,iblk,2)
                              enddo
                        endif
                        if (present(Sinz_rest)) then
                              do k = 1, nilyr
                                    Sinz_rest(i,j,k,n,iblk) = weight_1 * Sinz_cache(i,j,k,n,iblk,1) + &
                                                             weight_2 * Sinz_cache(i,j,k,n,iblk,2)
                              enddo
                        endif
                  endif
            enddo

            end subroutine restore_forcing_interpolate_cell

!=======================================================================

            subroutine restore_forcing_allocate_cache()

            if (.not. allocated(aicen_cache)) then
                  allocate(aicen_cache(nx_block,ny_block,ncat,max_blocks,2))
                  allocate(vicen_cache(nx_block,ny_block,ncat,max_blocks,2))
                  allocate(vsnon_cache(nx_block,ny_block,ncat,max_blocks,2))
                  allocate(Tsfc_cache(nx_block,ny_block,ncat,max_blocks,2))
                  allocate(trcrn_cache(nx_block,ny_block,restore_forcing_ntrcr,ncat,max_blocks,2))
                  allocate(Tinz_cache(nx_block,ny_block,nilyr,ncat,max_blocks,2))
                  allocate(Sinz_cache(nx_block,ny_block,nilyr,ncat,max_blocks,2))
            endif

            aicen_cache = c0
            vicen_cache = c0
            vsnon_cache = c0
            Tsfc_cache = c0
            trcrn_cache = c0
            Tinz_cache = c0
            Sinz_cache = c0

            if (.not. allocated(restore_forcing_work_g)) then
                  if (my_task == master_task) then
                        allocate(restore_forcing_work_g(nx_global,ny_global))
                  else
                        allocate(restore_forcing_work_g(1,1))
                  endif
            endif

            end subroutine restore_forcing_allocate_cache

!=======================================================================

#ifdef USE_NETCDF
            subroutine restore_forcing_validate_file()

            integer (kind=int_kind) :: time_len
            logical (kind=log_kind) :: found_Tsfc, found_Tinz, found_Sinz

            time_len = 0

            call restore_forcing_validate_var_4d('aicen', ncat, restore_forcing_varid_aicen, time_len)
            call restore_forcing_validate_var_4d('vicen', ncat, restore_forcing_varid_vicen, time_len)
            call restore_forcing_validate_var_4d('vsnon', ncat, restore_forcing_varid_vsnon, time_len)
            call restore_forcing_try_validate_var_5d('trcrn', restore_forcing_ntrcr, ncat, &
                                                                                     restore_forcing_varid_trcrn, time_len, &
                                                                                     restore_forcing_has_trcrn)

            call restore_forcing_try_validate_var_4d('Tsfc', ncat, restore_forcing_varid_Tsfc, &
                                                                                     time_len, found_Tsfc)
            call restore_forcing_try_validate_var_5d('Tinz', nilyr, ncat, restore_forcing_varid_Tinz, &
                                                                                     time_len, found_Tinz)
            call restore_forcing_try_validate_var_5d('Sinz', nilyr, ncat, restore_forcing_varid_Sinz, &
                                                                                     time_len, found_Sinz)

            if ((found_Tsfc .or. found_Tinz .or. found_Sinz) .and. &
                .not. (found_Tsfc .and. found_Tinz .and. found_Sinz)) then
                  call abort_ice(error_message='(restore_forcing_validate_file) incomplete thermodynamic schema: Tsfc, Tinz, and Sinz must be supplied together', &
                                                file=__FILE__, line=__LINE__)
            endif

            if (restore_forcing_has_trcrn .and. found_Tsfc) then
                  call abort_ice(error_message='(restore_forcing_validate_file) ambiguous schema: supply trcrn or Tsfc/Tinz/Sinz, not both', &
                                                file=__FILE__, line=__LINE__)
            endif

            if (found_Tsfc) then
                  restore_forcing_has_bc_thermo = .true.
            endif

            restore_forcing_nrecords = time_len

            if (restore_forcing_nrecords < 1) then
                  call abort_ice(error_message='(restore_forcing_validate_file) no records found', &
                                                file=__FILE__, line=__LINE__)
            endif

            if (trim(restore_forcing_data_type) == 'daily_netcdf') then
                  if (restore_forcing_cycle_year .and. restore_forcing_nrecords /= days_per_year) then
                        call abort_ice(error_message='(restore_forcing_validate_file) cyclic daily_netcdf file must contain exactly one model year', &
                                                      file=__FILE__, line=__LINE__)
                  endif
            else if (trim(restore_forcing_data_type) == 'monthly_netcdf') then
                  if (restore_forcing_cycle_year .and. restore_forcing_nrecords /= 12) then
                        call abort_ice(error_message='(restore_forcing_validate_file) cyclic monthly_netcdf file must contain exactly 12 records', &
                                                      file=__FILE__, line=__LINE__)
                  endif
            endif

            if (my_task == master_task) then
                  write(nu_diag,*) '  nrecords   = ', restore_forcing_nrecords
                  write(nu_diag,*) '  has_trcrn  = ', restore_forcing_has_trcrn
                  write(nu_diag,*) '  has_bc_thermo = ', restore_forcing_has_bc_thermo
            endif

            end subroutine restore_forcing_validate_file

!=======================================================================

            subroutine restore_forcing_validate_var_4d(varname, nk_expected, varid, time_len)

            character(len=*), intent(in) :: varname
            integer (kind=int_kind), intent(in) :: nk_expected
            integer (kind=int_kind), intent(out) :: varid
            integer (kind=int_kind), intent(inout) :: time_len

            integer (kind=int_kind) :: status, ndims, dimids(NF90_MAX_VAR_DIMS), dimlen
            integer (kind=int_kind) :: local_time_len
            character(len=*), parameter :: subname = '(restore_forcing_validate_var_4d)'

            local_time_len = 0
            if (my_task == master_task) then
                  status = nf90_inq_varid(restore_forcing_fid, trim(varname), varid)
                  call ice_check_nc(status, subname//' missing variable '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)

                  status = nf90_inquire_variable(restore_forcing_fid, varid, ndims=ndims, dimids=dimids)
                  call ice_check_nc(status, subname//' inquire variable '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)

                  if (ndims /= 4) then
                        call abort_ice(error_message=subname//' '//trim(varname)//' must have 4 dimensions (x,y,ncat,time)', &
                                                      file=__FILE__, line=__LINE__)
                  endif

                  status = nf90_inquire_dimension(restore_forcing_fid, dimids(1), len=dimlen)
                  call ice_check_nc(status, subname//' inquire x dim '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
                  if (dimlen /= nx_global) then
                        call abort_ice(error_message=subname//' '//trim(varname)//' x dimension mismatch', &
                                                      file=__FILE__, line=__LINE__)
                  endif

                  status = nf90_inquire_dimension(restore_forcing_fid, dimids(2), len=dimlen)
                  call ice_check_nc(status, subname//' inquire y dim '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
                  if (dimlen /= ny_global) then
                        call abort_ice(error_message=subname//' '//trim(varname)//' y dimension mismatch', &
                                                      file=__FILE__, line=__LINE__)
                  endif

                  status = nf90_inquire_dimension(restore_forcing_fid, dimids(3), len=dimlen)
                  call ice_check_nc(status, subname//' inquire category dim '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
                  if (dimlen /= nk_expected) then
                        call abort_ice(error_message=subname//' '//trim(varname)//' category dimension mismatch', &
                                                      file=__FILE__, line=__LINE__)
                  endif

                  status = nf90_inquire_dimension(restore_forcing_fid, dimids(4), len=local_time_len)
                  call ice_check_nc(status, subname//' inquire time dim '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
            endif

            call restore_forcing_sync_scalar(local_time_len)
            if (time_len == 0) then
                  time_len = local_time_len
            elseif (time_len /= local_time_len) then
                  call abort_ice(error_message=subname//' inconsistent time dimension for '//trim(varname), &
                                                file=__FILE__, line=__LINE__)
            endif

            call restore_forcing_sync_scalar(varid)

            end subroutine restore_forcing_validate_var_4d

!=======================================================================

            subroutine restore_forcing_validate_var_5d(varname, ntrcr_expected, nk_expected, varid, time_len)

            character(len=*), intent(in) :: varname
            integer (kind=int_kind), intent(in) :: ntrcr_expected, nk_expected
            integer (kind=int_kind), intent(out) :: varid
            integer (kind=int_kind), intent(inout) :: time_len

            integer (kind=int_kind) :: status, ndims, dimids(NF90_MAX_VAR_DIMS), dimlen
            integer (kind=int_kind) :: local_time_len
            character(len=*), parameter :: subname = '(restore_forcing_validate_var_5d)'

            local_time_len = 0
            if (my_task == master_task) then
                  status = nf90_inq_varid(restore_forcing_fid, trim(varname), varid)
                  call ice_check_nc(status, subname//' missing variable '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)

                  status = nf90_inquire_variable(restore_forcing_fid, varid, ndims=ndims, dimids=dimids)
                  call ice_check_nc(status, subname//' inquire variable '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)

                  if (ndims /= 5) then
                        call abort_ice(error_message=subname//' '//trim(varname)//' must have 5 dimensions (x,y,ntrcr,ncat,time)', &
                                                      file=__FILE__, line=__LINE__)
                  endif

                  status = nf90_inquire_dimension(restore_forcing_fid, dimids(1), len=dimlen)
                  call ice_check_nc(status, subname//' inquire x dim '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
                  if (dimlen /= nx_global) then
                        call abort_ice(error_message=subname//' '//trim(varname)//' x dimension mismatch', &
                                                      file=__FILE__, line=__LINE__)
                  endif

                  status = nf90_inquire_dimension(restore_forcing_fid, dimids(2), len=dimlen)
                  call ice_check_nc(status, subname//' inquire y dim '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
                  if (dimlen /= ny_global) then
                        call abort_ice(error_message=subname//' '//trim(varname)//' y dimension mismatch', &
                                                      file=__FILE__, line=__LINE__)
                  endif

                  status = nf90_inquire_dimension(restore_forcing_fid, dimids(3), len=dimlen)
                  call ice_check_nc(status, subname//' inquire tracer dim '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
                  if (dimlen /= ntrcr_expected) then
                        call abort_ice(error_message=subname//' '//trim(varname)//' tracer dimension mismatch', &
                                                      file=__FILE__, line=__LINE__)
                  endif

                  status = nf90_inquire_dimension(restore_forcing_fid, dimids(4), len=dimlen)
                  call ice_check_nc(status, subname//' inquire category dim '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
                  if (dimlen /= nk_expected) then
                        call abort_ice(error_message=subname//' '//trim(varname)//' category dimension mismatch', &
                                                      file=__FILE__, line=__LINE__)
                  endif

                  status = nf90_inquire_dimension(restore_forcing_fid, dimids(5), len=local_time_len)
                  call ice_check_nc(status, subname//' inquire time dim '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
            endif

            call restore_forcing_sync_scalar(local_time_len)
            if (time_len == 0) then
                  time_len = local_time_len
            elseif (time_len /= local_time_len) then
                  call abort_ice(error_message=subname//' inconsistent time dimension for '//trim(varname), &
                                                file=__FILE__, line=__LINE__)
            endif

            call restore_forcing_sync_scalar(varid)

            end subroutine restore_forcing_validate_var_5d

!=======================================================================

            subroutine restore_forcing_try_validate_var_4d(varname, nk_expected, varid, time_len, found_var)

            character(len=*), intent(in) :: varname
            integer (kind=int_kind), intent(in) :: nk_expected
            integer (kind=int_kind), intent(out) :: varid
            integer (kind=int_kind), intent(inout) :: time_len
            logical (kind=log_kind), intent(out) :: found_var

            integer (kind=int_kind) :: status

            found_var = .true.
            varid = -1

            if (my_task == master_task) then
                  status = nf90_inq_varid(restore_forcing_fid, trim(varname), varid)
                  if (status /= 0) found_var = .false.
            endif

            call broadcast_scalar(found_var, master_task)
            call restore_forcing_sync_scalar(varid)

            if (.not. found_var) return

            call restore_forcing_validate_var_4d(varname, nk_expected, varid, time_len)

            end subroutine restore_forcing_try_validate_var_4d

!=======================================================================

            subroutine restore_forcing_try_validate_var_5d(varname, ntrcr_expected, nk_expected, varid, time_len, found_var)

            character(len=*), intent(in) :: varname
            integer (kind=int_kind), intent(in) :: ntrcr_expected, nk_expected
            integer (kind=int_kind), intent(out) :: varid
            integer (kind=int_kind), intent(inout) :: time_len
            logical (kind=log_kind), intent(out) :: found_var

            integer (kind=int_kind) :: status

            found_var = .true.
            varid = -1

            if (my_task == master_task) then
                  status = nf90_inq_varid(restore_forcing_fid, trim(varname), varid)
                  if (status /= 0) then
                        found_var = .false.
                  endif
            endif

            call broadcast_scalar(found_var, master_task)
            call restore_forcing_sync_scalar(varid)

            if (.not. found_var) return

            call restore_forcing_validate_var_5d(varname, ntrcr_expected, nk_expected, varid, time_len)

            end subroutine restore_forcing_try_validate_var_5d
#endif

!=======================================================================

            subroutine restore_forcing_select_records(recnum_1, recnum_2, weight_2)

            integer (kind=int_kind), intent(out) :: recnum_1, recnum_2
            real (kind=dbl_kind), intent(out) :: weight_2

            integer (kind=int_kind) :: record_index0, month_len
            real (kind=dbl_kind) :: secday, record_position, month_fraction
            character(len=*), parameter :: subname = '(restore_forcing_select_records)'

            call icepack_query_parameters(secday_out=secday)

            if (trim(restore_forcing_data_type) == 'daily_netcdf') then
                  if (restore_forcing_cycle_year) then
                        record_position = real(yday - 1, kind=dbl_kind) + real(msec, kind=dbl_kind) / secday
                  else
                        record_position = real(compute_days_between(fyear_init, 1, 1, myear, mmonth, mday), kind=dbl_kind) + &
                                                      real(msec, kind=dbl_kind) / secday
                  endif
            else if (trim(restore_forcing_data_type) == 'monthly_netcdf') then
                  month_len = daymo(mmonth)
                  month_fraction = (real(mday - 1, kind=dbl_kind) + real(msec, kind=dbl_kind) / secday) / &
                                   real(month_len, kind=dbl_kind)

                  if (restore_forcing_cycle_year) then
                        record_position = real(mmonth - 1, kind=dbl_kind) + month_fraction
                  else
                        record_position = real(12 * (myear - fyear_init) + (mmonth - 1), kind=dbl_kind) + month_fraction
                  endif
            else
                  call abort_ice(error_message=subname//' unsupported restore_ice_data_type='// &
                                                trim(restore_forcing_data_type), file=__FILE__, line=__LINE__)
            endif

            record_index0 = int(record_position)
            weight_2 = record_position - real(record_index0, kind=dbl_kind)
            recnum_1 = record_index0 + 1
            recnum_2 = recnum_1 + 1

            if (restore_forcing_cycle_year) then
                  recnum_1 = mod(recnum_1 - 1, restore_forcing_nrecords) + 1
                  recnum_2 = mod(recnum_2 - 1, restore_forcing_nrecords) + 1
            else
                  if (recnum_1 > restore_forcing_nrecords) then
                        call abort_ice(error_message=subname//' model time exceeds restore_ice_data_file coverage', &
                                                      file=__FILE__, line=__LINE__)
                  endif
                  if (recnum_2 > restore_forcing_nrecords) then
                        if (weight_2 == c0) then
                              recnum_2 = recnum_1
                        else
                              call abort_ice(error_message=subname//' interpolation requires record beyond restore_ice_data_file coverage', &
                                                            file=__FILE__, line=__LINE__)
                        endif
                  endif
            endif

            end subroutine restore_forcing_select_records

!=======================================================================

            subroutine restore_forcing_ensure_records(recnum_1, recnum_2)

            integer (kind=int_kind), intent(in) :: recnum_1, recnum_2

            if (restore_forcing_loaded_record_1 == recnum_1 .and. &
                   restore_forcing_loaded_record_2 == recnum_2) then
                  return
            endif

            call restore_forcing_load_record(recnum_1, 1)
            restore_forcing_loaded_record_1 = recnum_1

            if (recnum_2 == recnum_1) then
                  aicen_cache(:,:,:,:,2) = aicen_cache(:,:,:,:,1)
                  vicen_cache(:,:,:,:,2) = vicen_cache(:,:,:,:,1)
                  vsnon_cache(:,:,:,:,2) = vsnon_cache(:,:,:,:,1)
                  if (restore_forcing_has_bc_thermo) then
                        Tsfc_cache(:,:,:,:,2) = Tsfc_cache(:,:,:,:,1)
                        Tinz_cache(:,:,:,:,:,2) = Tinz_cache(:,:,:,:,:,1)
                        Sinz_cache(:,:,:,:,:,2) = Sinz_cache(:,:,:,:,:,1)
                  endif
                  if (restore_forcing_has_trcrn) then
                        trcrn_cache(:,:,:,:,:,2) = trcrn_cache(:,:,:,:,:,1)
                  endif
            else
                  call restore_forcing_load_record(recnum_2, 2)
            endif
            restore_forcing_loaded_record_2 = recnum_2

            if (my_task == master_task) then
                  write(nu_diag,*) 'restore_forcing loaded records ', recnum_1, recnum_2
            endif

            end subroutine restore_forcing_ensure_records

!=======================================================================

#ifdef USE_NETCDF
            subroutine restore_forcing_load_record(recnum, slot)

            integer (kind=int_kind), intent(in) :: recnum, slot
            integer (kind=int_kind) :: n, nt

            do n = 1, ncat
                  call restore_forcing_read_4d_slice(restore_forcing_varid_aicen, 'aicen', n, recnum, &
                                                                                      aicen_cache(:,:,n,:,slot))
                  call restore_forcing_read_4d_slice(restore_forcing_varid_vicen, 'vicen', n, recnum, &
                                                                                      vicen_cache(:,:,n,:,slot))
                  call restore_forcing_read_4d_slice(restore_forcing_varid_vsnon, 'vsnon', n, recnum, &
                                                                                      vsnon_cache(:,:,n,:,slot))
                  if (restore_forcing_has_bc_thermo) then
                        call restore_forcing_read_4d_slice(restore_forcing_varid_Tsfc, 'Tsfc', n, recnum, &
                                                                                    Tsfc_cache(:,:,n,:,slot))
                        do nt = 1, nilyr
                              call restore_forcing_read_5d_slice(restore_forcing_varid_Tinz, 'Tinz', nt, n, recnum, &
                                                                                                  Tinz_cache(:,:,nt,n,:,slot))
                              call restore_forcing_read_5d_slice(restore_forcing_varid_Sinz, 'Sinz', nt, n, recnum, &
                                                                                                  Sinz_cache(:,:,nt,n,:,slot))
                        enddo
                  endif
                  if (restore_forcing_has_trcrn) then
                        do nt = 1, restore_forcing_ntrcr
                              call restore_forcing_read_5d_slice(restore_forcing_varid_trcrn, 'trcrn', nt, n, recnum, &
                                                                                                  trcrn_cache(:,:,nt,n,:,slot))
                        enddo
                  endif
            enddo

            end subroutine restore_forcing_load_record

!=======================================================================

            subroutine restore_forcing_read_4d_slice(varid, varname, ncat_idx, recnum, work)

            integer (kind=int_kind), intent(in) :: varid, ncat_idx, recnum
            character(len=*), intent(in) :: varname
            real (kind=dbl_kind), dimension(:,:,:), intent(inout) :: work

            integer (kind=int_kind) :: status
            character(len=*), parameter :: subname = '(restore_forcing_read_4d_slice)'

            restore_forcing_work_g = c0
            if (my_task == master_task) then
                  status = nf90_get_var(restore_forcing_fid, varid, restore_forcing_work_g, &
                                                             start=(/1,1,ncat_idx,recnum/), count=(/nx_global,ny_global,1,1/))
                  call ice_check_nc(status, subname//' cannot read '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
            endif

            call scatter_global(work, restore_forcing_work_g, master_task, distrb_info, &
                                                  field_loc_center, field_type_scalar)

            end subroutine restore_forcing_read_4d_slice

!=======================================================================

            subroutine restore_forcing_read_5d_slice(varid, varname, ntrcr_idx, ncat_idx, recnum, work)

            integer (kind=int_kind), intent(in) :: varid, ntrcr_idx, ncat_idx, recnum
            character(len=*), intent(in) :: varname
            real (kind=dbl_kind), dimension(:,:,:), intent(inout) :: work

            integer (kind=int_kind) :: status
            character(len=*), parameter :: subname = '(restore_forcing_read_5d_slice)'

            restore_forcing_work_g = c0
            if (my_task == master_task) then
                  status = nf90_get_var(restore_forcing_fid, varid, restore_forcing_work_g, &
                                                             start=(/1,1,ntrcr_idx,ncat_idx,recnum/), count=(/nx_global,ny_global,1,1,1/))
                  call ice_check_nc(status, subname//' cannot read '//trim(varname), &
                                                      file=__FILE__, line=__LINE__)
            endif

            call scatter_global(work, restore_forcing_work_g, master_task, distrb_info, &
                                                  field_loc_center, field_type_scalar)

            end subroutine restore_forcing_read_5d_slice
#endif

!=======================================================================

            subroutine restore_forcing_extend_open_edges(aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)

            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: &
                  aicen_rest, vicen_rest, vsnon_rest

            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: &
                  trcrn_rest

            integer (kind=int_kind) :: i, j, iblk, n, nt, ilo, ihi, jlo, jhi, ibc

            do iblk = 1, nblocks
                  ilo = restore_forcing_ilo(iblk)
                  ihi = restore_forcing_ihi(iblk)
                  jlo = restore_forcing_jlo(iblk)
                  jhi = restore_forcing_jhi(iblk)

                  if (restore_forcing_west_block(iblk)) then
                        do n = 1, ncat
                              do j = 1, ny_block
                                    do i = 1, ilo
                                          aicen_rest(i,j,n,iblk) = aicen_rest(ilo,j,n,iblk)
                                          vicen_rest(i,j,n,iblk) = vicen_rest(ilo,j,n,iblk)
                                          vsnon_rest(i,j,n,iblk) = vsnon_rest(ilo,j,n,iblk)
                                          do nt = 1, restore_forcing_ntrcr
                                                trcrn_rest(i,j,nt,n,iblk) = trcrn_rest(ilo,j,nt,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif

                  if (restore_forcing_east_block(iblk)) then
                        ibc = restore_forcing_east_bc(iblk)

                        do n = 1, ncat
                              do j = 1, ny_block
                                    do i = ihi, ibc
                                          aicen_rest(i,j,n,iblk) = aicen_rest(ihi,j,n,iblk)
                                          vicen_rest(i,j,n,iblk) = vicen_rest(ihi,j,n,iblk)
                                          vsnon_rest(i,j,n,iblk) = vsnon_rest(ihi,j,n,iblk)
                                          do nt = 1, restore_forcing_ntrcr
                                                trcrn_rest(i,j,nt,n,iblk) = trcrn_rest(ihi,j,nt,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif

                  if (restore_forcing_south_block(iblk)) then
                        do n = 1, ncat
                              do j = 1, jlo
                                    do i = 1, nx_block
                                          aicen_rest(i,j,n,iblk) = aicen_rest(i,jlo,n,iblk)
                                          vicen_rest(i,j,n,iblk) = vicen_rest(i,jlo,n,iblk)
                                          vsnon_rest(i,j,n,iblk) = vsnon_rest(i,jlo,n,iblk)
                                          do nt = 1, restore_forcing_ntrcr
                                                trcrn_rest(i,j,nt,n,iblk) = trcrn_rest(i,jlo,nt,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif

                  if (restore_forcing_north_block(iblk)) then
                        ibc = restore_forcing_north_bc(iblk)

                        do n = 1, ncat
                              do j = jhi, ibc
                                    do i = 1, nx_block
                                          aicen_rest(i,j,n,iblk) = aicen_rest(i,jhi,n,iblk)
                                          vicen_rest(i,j,n,iblk) = vicen_rest(i,jhi,n,iblk)
                                          vsnon_rest(i,j,n,iblk) = vsnon_rest(i,jhi,n,iblk)
                                          do nt = 1, restore_forcing_ntrcr
                                                trcrn_rest(i,j,nt,n,iblk) = trcrn_rest(i,jhi,nt,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif
            enddo

            end subroutine restore_forcing_extend_open_edges

!=======================================================================

            subroutine restore_forcing_extend_open_edges_bc(Tsfc_rest, Tinz_rest, Sinz_rest)

            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: Tsfc_rest
            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: Tinz_rest, Sinz_rest

            integer (kind=int_kind) :: i, j, iblk, n, k, ilo, ihi, jlo, jhi, ibc

            do iblk = 1, nblocks
                  ilo = restore_forcing_ilo(iblk)
                  ihi = restore_forcing_ihi(iblk)
                  jlo = restore_forcing_jlo(iblk)
                  jhi = restore_forcing_jhi(iblk)

                  if (restore_forcing_west_block(iblk)) then
                        do n = 1, ncat
                              do j = 1, ny_block
                                    do i = 1, ilo
                                          Tsfc_rest(i,j,n,iblk) = Tsfc_rest(ilo,j,n,iblk)
                                          do k = 1, nilyr
                                                Tinz_rest(i,j,k,n,iblk) = Tinz_rest(ilo,j,k,n,iblk)
                                                Sinz_rest(i,j,k,n,iblk) = Sinz_rest(ilo,j,k,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif

                  if (restore_forcing_east_block(iblk)) then
                        ibc = restore_forcing_east_bc(iblk)

                        do n = 1, ncat
                              do j = 1, ny_block
                                    do i = ihi, ibc
                                          Tsfc_rest(i,j,n,iblk) = Tsfc_rest(ihi,j,n,iblk)
                                          do k = 1, nilyr
                                                Tinz_rest(i,j,k,n,iblk) = Tinz_rest(ihi,j,k,n,iblk)
                                                Sinz_rest(i,j,k,n,iblk) = Sinz_rest(ihi,j,k,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif

                  if (restore_forcing_south_block(iblk)) then
                        do n = 1, ncat
                              do j = 1, jlo
                                    do i = 1, nx_block
                                          Tsfc_rest(i,j,n,iblk) = Tsfc_rest(i,jlo,n,iblk)
                                          do k = 1, nilyr
                                                Tinz_rest(i,j,k,n,iblk) = Tinz_rest(i,jlo,k,n,iblk)
                                                Sinz_rest(i,j,k,n,iblk) = Sinz_rest(i,jlo,k,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif

                  if (restore_forcing_north_block(iblk)) then
                        ibc = restore_forcing_north_bc(iblk)

                        do n = 1, ncat
                              do j = jhi, ibc
                                    do i = 1, nx_block
                                          Tsfc_rest(i,j,n,iblk) = Tsfc_rest(i,jhi,n,iblk)
                                          do k = 1, nilyr
                                                Tinz_rest(i,j,k,n,iblk) = Tinz_rest(i,jhi,k,n,iblk)
                                                Sinz_rest(i,j,k,n,iblk) = Sinz_rest(i,jhi,k,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif
            enddo

            end subroutine restore_forcing_extend_open_edges_bc

!=======================================================================

            subroutine restore_forcing_disable_selected_edges(aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)

            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: &
                  aicen_rest, vicen_rest, vsnon_rest

            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: &
                  trcrn_rest

            if (.not. restore_forcing_use_west) then
                  call restore_forcing_set_edge_from_state('west', aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
            endif
            if (.not. restore_forcing_use_east) then
                  call restore_forcing_set_edge_from_state('east', aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
            endif
            if (.not. restore_forcing_use_south) then
                  call restore_forcing_set_edge_from_state('south', aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
            endif
            if (.not. restore_forcing_use_north) then
                  call restore_forcing_set_edge_from_state('north', aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
            endif

            end subroutine restore_forcing_disable_selected_edges

!=======================================================================

            subroutine restore_forcing_set_edge_from_state(edge_name, aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)

            character(len=*), intent(in) :: edge_name

            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: &
                  aicen_rest, vicen_rest, vsnon_rest

            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: &
                  trcrn_rest

            integer (kind=int_kind) :: i, j, iblk, n, nt, ilo, ihi, jlo, jhi, ibc

            do iblk = 1, nblocks
                  ilo = restore_forcing_ilo(iblk)
                  ihi = restore_forcing_ihi(iblk)
                  jlo = restore_forcing_jlo(iblk)
                  jhi = restore_forcing_jhi(iblk)

                  if (edge_name == 'west' .and. restore_forcing_west_block(iblk)) then
                        do n = 1, ncat
                              do j = 1, ny_block
                                    do i = 1, ilo
                                          aicen_rest(i,j,n,iblk) = aicen(ilo,j,n,iblk)
                                          vicen_rest(i,j,n,iblk) = vicen(ilo,j,n,iblk)
                                          vsnon_rest(i,j,n,iblk) = vsnon(ilo,j,n,iblk)
                                          do nt = 1, restore_forcing_ntrcr
                                                trcrn_rest(i,j,nt,n,iblk) = trcrn(ilo,j,nt,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif

                  if (edge_name == 'east' .and. restore_forcing_east_block(iblk)) then
                        ibc = restore_forcing_east_bc(iblk)
                        do n = 1, ncat
                              do j = 1, ny_block
                                    do i = ihi, ibc
                                          aicen_rest(i,j,n,iblk) = aicen(ihi,j,n,iblk)
                                          vicen_rest(i,j,n,iblk) = vicen(ihi,j,n,iblk)
                                          vsnon_rest(i,j,n,iblk) = vsnon(ihi,j,n,iblk)
                                          do nt = 1, restore_forcing_ntrcr
                                                trcrn_rest(i,j,nt,n,iblk) = trcrn(ihi,j,nt,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif

                  if (edge_name == 'south' .and. restore_forcing_south_block(iblk)) then
                        do n = 1, ncat
                              do j = 1, jlo
                                    do i = 1, nx_block
                                          aicen_rest(i,j,n,iblk) = aicen(i,jlo,n,iblk)
                                          vicen_rest(i,j,n,iblk) = vicen(i,jlo,n,iblk)
                                          vsnon_rest(i,j,n,iblk) = vsnon(i,jlo,n,iblk)
                                          do nt = 1, restore_forcing_ntrcr
                                                trcrn_rest(i,j,nt,n,iblk) = trcrn(i,jlo,nt,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif

                  if (edge_name == 'north' .and. restore_forcing_north_block(iblk)) then
                        ibc = restore_forcing_north_bc(iblk)
                        do n = 1, ncat
                              do j = jhi, ibc
                                    do i = 1, nx_block
                                          aicen_rest(i,j,n,iblk) = aicen(i,jhi,n,iblk)
                                          vicen_rest(i,j,n,iblk) = vicen(i,jhi,n,iblk)
                                          vsnon_rest(i,j,n,iblk) = vsnon(i,jhi,n,iblk)
                                          do nt = 1, restore_forcing_ntrcr
                                                trcrn_rest(i,j,nt,n,iblk) = trcrn(i,jhi,nt,n,iblk)
                                          enddo
                                    enddo
                              enddo
                        enddo
                  endif
            enddo

            end subroutine restore_forcing_set_edge_from_state

!=======================================================================

            subroutine restore_forcing_apply_land_mask(aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)

            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: &
                  aicen_rest, vicen_rest, vsnon_rest

            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: &
                  trcrn_rest

            integer (kind=int_kind) :: iblk

            do iblk = 1, nblocks
                  if (restore_forcing_west_block(iblk)) then
                        call restore_forcing_mask_range(1, restore_forcing_ilo(iblk), &
                             1, ny_block, iblk, aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
                  endif
                  if (restore_forcing_east_block(iblk)) then
                        call restore_forcing_mask_range(restore_forcing_ihi(iblk), &
                             restore_forcing_east_bc(iblk), 1, ny_block, iblk, &
                             aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
                  endif
                  if (restore_forcing_south_block(iblk)) then
                        call restore_forcing_mask_range(1, nx_block, 1, &
                             restore_forcing_jlo(iblk), iblk, &
                             aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
                  endif
                  if (restore_forcing_north_block(iblk)) then
                        call restore_forcing_mask_range(1, nx_block, &
                             restore_forcing_jhi(iblk), restore_forcing_north_bc(iblk), iblk, &
                             aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)
                  endif
            enddo

            end subroutine restore_forcing_apply_land_mask

!=======================================================================

            subroutine restore_forcing_mask_range(i1, i2, j1, j2, iblk, &
                                                  aicen_rest, vicen_rest, vsnon_rest, trcrn_rest)

            integer (kind=int_kind), intent(in) :: i1, i2, j1, j2, iblk
            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: &
                  aicen_rest, vicen_rest, vsnon_rest
            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: trcrn_rest

            integer (kind=int_kind) :: i, j, n, nt

            do n = 1, ncat
                  do j = j1, j2
                        do i = i1, i2
                              if (hm(i,j,iblk) <= c0) then
                                    aicen_rest(i,j,n,iblk) = c0
                                    vicen_rest(i,j,n,iblk) = c0
                                    vsnon_rest(i,j,n,iblk) = c0
                                    do nt = 1, restore_forcing_ntrcr
                                          trcrn_rest(i,j,nt,n,iblk) = c0
                                    enddo
                              endif
                        enddo
                  enddo
            enddo

            end subroutine restore_forcing_mask_range

!=======================================================================

            subroutine restore_forcing_apply_land_mask_bc(Tsfc_rest, Tinz_rest, Sinz_rest)

            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: Tsfc_rest
            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: Tinz_rest, Sinz_rest

            integer (kind=int_kind) :: iblk

            do iblk = 1, nblocks
                  if (restore_forcing_west_block(iblk)) then
                        call restore_forcing_mask_range_bc(1, restore_forcing_ilo(iblk), &
                             1, ny_block, iblk, Tsfc_rest, Tinz_rest, Sinz_rest)
                  endif
                  if (restore_forcing_east_block(iblk)) then
                        call restore_forcing_mask_range_bc(restore_forcing_ihi(iblk), &
                             restore_forcing_east_bc(iblk), 1, ny_block, iblk, &
                             Tsfc_rest, Tinz_rest, Sinz_rest)
                  endif
                  if (restore_forcing_south_block(iblk)) then
                        call restore_forcing_mask_range_bc(1, nx_block, 1, &
                             restore_forcing_jlo(iblk), iblk, Tsfc_rest, Tinz_rest, Sinz_rest)
                  endif
                  if (restore_forcing_north_block(iblk)) then
                        call restore_forcing_mask_range_bc(1, nx_block, &
                             restore_forcing_jhi(iblk), restore_forcing_north_bc(iblk), iblk, &
                             Tsfc_rest, Tinz_rest, Sinz_rest)
                  endif
            enddo

            end subroutine restore_forcing_apply_land_mask_bc

!=======================================================================

            subroutine restore_forcing_mask_range_bc(i1, i2, j1, j2, iblk, &
                                                     Tsfc_rest, Tinz_rest, Sinz_rest)

            integer (kind=int_kind), intent(in) :: i1, i2, j1, j2, iblk
            real (kind=dbl_kind), dimension (:,:,:,:), intent(inout) :: Tsfc_rest
            real (kind=dbl_kind), dimension (:,:,:,:,:), intent(inout) :: Tinz_rest, Sinz_rest

            integer (kind=int_kind) :: i, j, n, k

            do n = 1, ncat
                  do j = j1, j2
                        do i = i1, i2
                              if (hm(i,j,iblk) <= c0) then
                                    Tsfc_rest(i,j,n,iblk) = c0
                                    do k = 1, nilyr
                                          Tinz_rest(i,j,k,n,iblk) = c0
                                          Sinz_rest(i,j,k,n,iblk) = c0
                                    enddo
                              endif
                        enddo
                  enddo
            enddo

            end subroutine restore_forcing_mask_range_bc

!=======================================================================

            subroutine restore_forcing_sync_scalar(value)

            integer (kind=int_kind), intent(inout) :: value

            call broadcast_scalar(value, master_task)

            end subroutine restore_forcing_sync_scalar

!=======================================================================

            end module ice_restore_forcing
