! ***********************************************************************
!
!   Copyright (C) 2011  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      use interp_1d_def
      use interp_1d_lib
      use custom_eos
      
      implicit none
      real(dp), pointer, save :: T_old(:), rho_old(:), r_old(:) => null()

      
      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% eos_rq% other_eos_frac => custom_eos_frac
         s% eos_rq% other_eos_component => custom_eos_component
         s% eos_rq% other_eos_results => custom_eos_results
         
         ! Include core heating and irradiated atmosphere
         s% other_energy => core_heating
         s% other_surface_PT => irradiated_atm
         
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         integer :: k
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (s% use_other_energy) call core_heating(id, ierr)
         if (ierr /= 0) return
         
         extras_finish_step = keep_going
      end function extras_finish_step
      
      subroutine core_heating(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         
         real(dp) :: cv_core, M_core, dTdt_core, L_radio
         integer :: ival=300
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !dTdt_core = ( s%T(s%nz) - s%T_start(s%nz) )/ s%dt
         dTdt_core = s%T(s%nz) * s%dxh_lnT(s%nz) / s%dt
         
         cv_core = s% x_ctrl(1)
         M_core = s% x_ctrl(2)
         L_radio = s% x_ctrl(3)
         
         !print *, "dTdt_core ", dTdt_core
         
         do k = 1, s%nz-1
           s% extra_heat(k) = 0.0d0
         end do
         s% extra_heat(s%nz) = -cv_core*M_core*dTdt_core/s%dm(s%nz) + L_radio/s%dm(s%nz)
         !ierr = -1
      end subroutine core_heating
      
      subroutine irradiated_atm(id, &
            skip_partials, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
         use const_def, only: dp
         use star_def
         !use star_lib, only: star_get_surf_PT
         integer, intent(in) :: id
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr
         
         type (star_info), pointer :: s
         !real(dp) :: Tint4, Tirr4, cg, tau, f, T4, sqrt3, g
         type(auto_diff_real_4var_order1) :: L, lnR, lnM, lnkap, M, R, kap, g, lnT, lnP
         type(auto_diff_real_4var_order1) :: Tint4, Tirr4, cg, f, T4
         real(dp) :: sqrt3, tau
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         g = standard_cgrav * s%M(1) / pow2(s%R(1))
         sqrt3 = (3d0)
         
         ! Set up auto_diff base variables.
         ! Note that the ordering of which variable goes in which d1val slot doesn't matter.
         ! We just need to use the same ordering when reading the results out.

         L = s%L(1)
         L%d1val1 = 1d0

         lnR = log(s%R(1))
         lnR%d1val2 = 1d0

         lnM = log(s%M(1))
         lnM%d1val3 = 1d0

         lnkap = log(s%opacity(1))
         lnkap%d1val4 = 1d0

         ! Calculate lnT and lnP

         M = exp(lnM)
         R = exp(lnR)
         kap = exp(lnkap)

         f = s% x_ctrl(4)
         tau = s% x_ctrl(5)
         cg = s% x_ctrl(6)
         !if (s% atm_irradiated_kap_v_div_kap_th == 0) then
         !  cg = 1.0d0
         !end if
         Tint4 = L / (4d0 * pi * boltz_sigma * pow2(R))
         Tirr4 = s% irradiation_flux / boltz_sigma
         
         T4 = 0.75d0*Tint4*(2d0/3d0 + tau) + 0.75d0*Tirr4*f*(2d0/3d0 + 1/cg/sqrt3 + (cg/sqrt3 - 1/sqrt3/cg)*exp(-cg*tau*sqrt3))
         lnT = 0.25d0*log(T4)
         lnP = log(g / kap)
         
         lnT_surf = lnT%val
         dlnT_dL = lnT%d1val1
         dlnT_dlnR = lnT%d1val2
         dlnT_dlnM = lnT%d1val3
         dlnT_dlnkap = lnT%d1val4
         
         lnP_surf = lnP%val
         dlnP_dL = lnP%d1val1
         dlnP_dlnR = lnP%d1val2
         dlnP_dlnM = lnP%d1val3
         dlnP_dlnkap = lnP%d1val4
         !ierr = -1
         
         
      end subroutine irradiated_atm
         
         

      end module run_star_extras
      
