!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Follows notably code from JB.
!
! Note: will debug/clean/note and explain tomorrow - also will prepare few questions following up 
!       on the various erroneously applications I couldn't avoid below.
!
! - Selim
!
! Compile with:
!   gfortran -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -O3 -cpp -o evolve_trapped evolve_trapped_predraft.f90
! Then run with
!   ./evolve_trapped
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program evolve_trapped
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)

  integer :: i,j, ii, jj, ll, kout, ind1, ind2
  integer, parameter :: maxstep_md = 100000
  integer, parameter :: maxstep_bg = 100000
  integer, parameter :: stepsize=2**3
  real(dp), parameter :: Nstep = 1e-3_dp !7.5/640./4.
  real(dp), parameter :: Nstep_bg = 1e-3_dp !7.5/640./4.

  integer, parameter :: kbins = 300
  integer, parameter :: pibins = 120

  real(dp), dimension(kbins) :: kvector, kbinwidths, k_ic, &
                                lnaH_ic, a_init, a_end, &
                                hubble_ic
  real(dp), dimension(pibins) :: piarray

! Pack the phase space variables
  integer, parameter :: dl = kind(1.d0)
  integer, parameter :: nfld=2
  integer, parameter :: nvar = 6

  complex(dp), dimension(nvar, kbins) :: y

  ! real(dp), parameter :: k_pivot = 2e-3_dp
  real(dp), parameter :: kfactor = 1e-60_dp
  real(dp), parameter :: g2 = 1e2_dp
  real(dp), parameter :: m2 = 1e-4_dp
  real(dp), dimension(1) :: phi_init, dphi_init
  real(dp), parameter :: phi_trap = 15.2e0_dp
  real(dp), parameter :: N_pivot = 55e0_dp
  real(dp), parameter :: pibin = 0.00694854


  real(dp), parameter :: Mpl = 1e0_dp

  real(dp) :: k_pivot

  complex(dp) :: sint

  real(dp) :: hubble_init, N_init, N_end, &
              hubble_pivot, a_pivot, N_trap, &
              a_trap, phi_end, pot_end, N_diff, phi_pivot, &
              dphi_pivot, phi_ic, pot_ic, eps_ic, &
              dphi_ic, N_start

  real(dp) :: polerr, omega_chi, omega_inflaton

  real(dp) :: k

  real(dp), dimension(maxstep_bg,kbins) :: hubble_array, N_array, &
									   phi_array, dphi_array, &
                                       pot_array, dpot_array, &
									   eps_array, lnaH_array

  logical :: bg_finish

  open(unit=97,file='inflaton_bg_traj_long_new2.out')
  open(unit=98,file='inflaton_mode_traj_long_new2.out')
  open(unit=99,file='source_traj_long_new2.out')

  kvector = (/1.00000000e-04,   1.08005237e-04,   1.16651313e-04,&
			 1.25989528e-04,   1.36075289e-04,   1.46968439e-04,&
			 1.58733611e-04,   1.71440614e-04,   1.85164842e-04,&
			 1.99987727e-04,   2.15997219e-04,   2.33288310e-04,&
			 2.51963593e-04,   2.72133877e-04,   2.93918840e-04,&
			 3.17447741e-04,   3.42860186e-04,   3.70306958e-04,&
			 3.99950910e-04,   4.31967930e-04,   4.66547988e-04,&
			 5.03896263e-04,   5.44234355e-04,   5.87801607e-04,&
			 6.34856522e-04,   6.85678294e-04,   7.40568469e-04,&
			 7.99852734e-04,   8.63882844e-04,   9.33038717e-04,&
			 1.00773068e-03,   1.08840192e-03,   1.17553107e-03,&
			 1.26963513e-03,   1.37127243e-03,   1.48104605e-03,&
			 1.59960730e-03,   1.72765966e-03,   1.86596292e-03,&
			 2.01533769e-03,   2.17667025e-03,   2.35091788e-03,&
			 2.53911443e-03,   2.74237657e-03,   2.96191033e-03,&
			 3.19901828e-03,   3.45510729e-03,   3.73169684e-03,&
			 4.03042803e-03,   4.35307336e-03,   4.70154722e-03,&
			 5.07791724e-03,   5.48441658e-03,   5.92345715e-03,&
			 6.39764396e-03,   6.90979055e-03,   7.46293569e-03,&
			 8.06036141e-03,   8.70561248e-03,   9.40251743e-03,&
			 1.01552113e-02,   1.09681601e-02,   1.18461873e-02,&
			 1.27945027e-02,   1.38187331e-02,   1.49249555e-02,&
			 1.61197336e-02,   1.74101565e-02,   1.88038809e-02,&
			 2.03091762e-02,   2.19349740e-02,   2.36909207e-02,&
			 2.55874352e-02,   2.76357701e-02,   2.98480792e-02,&
			 3.22374888e-02,   3.48181763e-02,   3.76054540e-02,&
			 4.06158599e-02,   4.38672559e-02,   4.73789339e-02,&
			 5.11717301e-02,   5.52681486e-02,   5.96924951e-02,&
			 6.44710211e-02,   6.96320794e-02,   7.52062927e-02,&
			 8.12267350e-02,   8.77291280e-02,   9.47520530e-02,&
			 1.02337180e-01,   1.10529514e-01,   1.19377664e-01,&
			 1.28934130e-01,   1.39255613e-01,   1.50403355e-01,&
			 1.62443501e-01,   1.75447489e-01,   1.89492477e-01,&
			 2.04661800e-01,   2.21045463e-01,   2.38740677e-01,&
			 2.57852435e-01,   2.78494135e-01,   3.00788252e-01,&
			 3.24867066e-01,   3.50873446e-01,   3.78961698e-01,&
			 4.09298482e-01,   4.42063797e-01,   4.77452054e-01,&
			 5.15673224e-01,   5.56954091e-01,   6.01539588e-01,&
			 6.49694260e-01,   7.01703829e-01,   7.57876886e-01,&
			 8.18546731e-01,   8.84073340e-01,   9.54845510e-01,&
			 1.03128316e+00,   1.11383983e+00,   1.20300535e+00,&
			 1.29930878e+00,   1.40332154e+00,   1.51566076e+00,&
			 1.63699300e+00,   1.76803818e+00,   1.90957383e+00,&
			 2.06243975e+00,   2.22754295e+00,   2.40586305e+00,&
			 2.59845810e+00,   2.80647085e+00,   3.03113550e+00,&
			 3.27378510e+00,   3.53585937e+00,   3.81891330e+00,&
			 4.12462638e+00,   4.45481252e+00,   4.81143084e+00,&
			 5.19659730e+00,   5.61259726e+00,   6.06189899e+00,&
			 6.54716840e+00,   7.07128478e+00,   7.63735792e+00,&
			 8.24874655e+00,   8.90907830e+00,   9.62227117e+00,&
			 1.03925568e+01,   1.12245057e+01,   1.21230540e+01,&
			 1.30935333e+01,   1.41417017e+01,   1.52737785e+01,&
			 1.64964807e+01,   1.78170632e+01,   1.92433614e+01,&
			 2.07838382e+01,   2.24476338e+01,   2.42446202e+01,&
			 2.61854596e+01,   2.82816678e+01,   3.05456825e+01,&
			 3.29909369e+01,   3.56319397e+01,   3.84843611e+01,&
			 4.15651256e+01,   4.48925126e+01,   4.84862648e+01,&
			 5.23677054e+01,   5.65598646e+01,   6.10876161e+01,&
			 6.59778248e+01,   7.12595063e+01,   7.69639990e+01,&
			 8.31251499e+01,   8.97795155e+01,   9.69665789e+01,&
			 1.04728984e+02,   1.13112788e+02,   1.22167735e+02,&
			 1.31947552e+02,   1.42510267e+02,   1.53918552e+02,&
			 1.66240098e+02,   1.79548012e+02,   1.93921257e+02,&
			 2.09445114e+02,   2.26211693e+02,   2.44320476e+02,&
			 2.63878910e+02,   2.85003044e+02,   3.07818214e+02,&
			 3.32459793e+02,   3.59073989e+02,   3.87818715e+02,&
			 4.18864524e+02,   4.52395623e+02,   4.88610967e+02,&
			 5.27725435e+02,   5.69971109e+02,   6.15598650e+02,&
			 6.64878784e+02,   7.18103909e+02,   7.75589832e+02,&
			 8.37677640e+02,   9.04735724e+02,   9.77161967e+02,&
			 1.05538610e+03,   1.13987227e+03,   1.23112175e+03,&
			 1.32967597e+03,   1.43611969e+03,   1.55108448e+03,&
			 1.67525247e+03,   1.80936041e+03,   1.95420401e+03,&
 			 2.11064268e+03,   2.27960464e+03,   2.46209240e+03,&
			 2.65918874e+03,   2.87206312e+03,   3.10197859e+03,&
			 3.35029934e+03,   3.61849876e+03,   3.90816818e+03,&
			 4.22102632e+03,   4.55892950e+03,   4.92388263e+03,&
			 5.31805113e+03,   5.74377375e+03,   6.20357648e+03,&
			 6.70018750e+03,   7.23655342e+03,   7.81585671e+03,&
			 8.44153460e+03,   9.11729948e+03,   9.84716096e+03,&
			 1.06354496e+04,   1.14868426e+04,   1.24063916e+04,&
			 1.33995527e+04,   1.44722187e+04,   1.56307542e+04,&
			 1.68820332e+04,   1.82334800e+04,   1.96931134e+04,&
			 2.12695939e+04,   2.29722754e+04,   2.48112606e+04,&
			 2.67974609e+04,   2.89426612e+04,   3.12595900e+04,&
			 3.37619944e+04,   3.64647222e+04,   3.93838098e+04,&
			 4.25365773e+04,   4.59417313e+04,   4.96194760e+04,&
			 5.35916329e+04,   5.78817704e+04,   6.25153435e+04,&
			 6.75198452e+04,   7.29249692e+04,   7.87627861e+04,&
			 8.50679342e+04,   9.18778243e+04,   9.92328623e+04,&
			 1.07176689e+05,   1.15756437e+05,   1.25023015e+05,&
			 1.35031404e+05,   1.45840988e+05,   1.57515906e+05,&
			 1.70125428e+05,   1.83744372e+05,   1.98453546e+05,&
			 2.14340223e+05,   2.31498667e+05,   2.50030685e+05,&
			 2.70046235e+05,   2.91664078e+05,   3.15012480e+05,&
			 3.40229977e+05,   3.67466194e+05,   3.96882735e+05,&
			 4.28654141e+05,   4.62968923e+05,   5.00030684e+05,&
			 5.40059328e+05,   5.83292359e+05,   6.29986298e+05,&
			 6.80418197e+05,   7.34887289e+05,   7.93716762e+05,&
			 8.57255673e+05,   9.25881025e+05,   1.00000000e+06/)

  kbinwidths = (/7.70286e-6, 8.31949e-6, 8.98548e-6, 9.70479e-6, &
               0.0000104817, 0.0000113208, 0.000012227, 0.0000132058, 0.000014263, &
               0.0000154048, 0.000016638, 0.0000179699, 0.0000194084, 0.0000209621, &
               0.0000226401, 0.0000244525, 0.00002641, 0.0000285242, 0.0000308076, &
               0.0000332739, 0.0000359375, 0.0000388144, 0.0000419216, 0.0000452775, &
               0.0000489021, 0.0000528168, 0.0000570449, 0.0000616115, 0.0000665436, &
               0.0000718706, 0.000077624, 0.000083838, 0.0000905495, 0.0000977982, &
               0.000105627, 0.000114083, 0.000123215, 0.000133079, 0.000143732, &
               0.000155239, 0.000167666, 0.000181088, 0.000195584, 0.000211241, &
               0.000228152, 0.000246416, 0.000266142, 0.000287447, 0.000310458, &
               0.000335311, 0.000362153, 0.000391145, 0.000422457, 0.000456275, &
               0.000492801, 0.000532251, 0.000574859, 0.000620878, 0.000670581, &
               0.000724262, 0.000782241, 0.000844861, 0.000912495, 0.000985542, &
               0.00106444, 0.00114965, 0.00124168, 0.00134108, 0.00144844, &
               0.00156439, 0.00168962, 0.00182488, 0.00197096, 0.00212874, &
               0.00229915, 0.00248321, 0.00268199, 0.00289669, 0.00312858, &
               0.00337903, 0.00364953, 0.00394168, 0.00425722, 0.00459803, &
               0.00496611, 0.00536366, 0.00579303, 0.00625678, 0.00675765, &
               0.00729861, 0.00788288, 0.00851393, 0.00919549, 0.00993161, &
               0.0107267, 0.0115853, 0.0125128, 0.0135145, 0.0145963, 0.0157648, &
               0.0170268, 0.0183898, 0.019862, 0.021452, 0.0231693, 0.025024, &
               0.0270273, 0.0291909, 0.0315277, 0.0340515, 0.0367774, 0.0397215, &
               0.0429014, 0.0463357, 0.050045, 0.0540512, 0.0583781, 0.0630514, &
               0.0680989, 0.0735503, 0.0794382, 0.0857974, 0.0926657, 0.100084, &
               0.108096, 0.116749, 0.126095, 0.136189, 0.147092, 0.158867, 0.171584, &
               0.18532, 0.200155, 0.216178, 0.233484, 0.252175, 0.272362, 0.294165, &
               0.317714, 0.343148, 0.370617, 0.400286, 0.43233, 0.466939, 0.504319, &
               0.544691, 0.588294, 0.635389, 0.686253, 0.741189, 0.800523, 0.864607, &
               0.933821, 1.00858, 1.08931, 1.17652, 1.2707, 1.37242, 1.48229, &
               1.60095, 1.72911, 1.86753, 2.01703, 2.17849, 2.35289, 2.54124, &
               2.74468, 2.96439, 3.2017, 3.458, 3.73482, 4.03381, 4.35672, 4.70549, &
               5.08217, 5.48901, 5.92842, 6.40301, 6.91558, 7.46919, 8.06712, &
               8.71291, 9.4104, 10.1637, 10.9774, 11.8561, 12.8052, 13.8303, &
               14.9375, 16.1332, 17.4247, 18.8196, 20.3262, 21.9534, 23.7108, &
               25.6089, 27.6589, 29.8731, 32.2645, 34.8474, 37.637, 40.6499, 43.904, &
               47.4186, 51.2146, 55.3145, 59.7425, 64.525, 69.6904, 75.2693, &
               81.2948, 87.8026, 94.8315, 102.423, 110.622, 119.478, 129.042, &
               139.372, 150.529, 162.58, 175.595, 189.651, 204.833, 221.231, &
               238.941, 258.068, 278.727, 301.04, 325.139, 351.167, 379.279, &
               409.641, 442.434, 477.852, 516.105, 557.421, 602.044, 650.239, &
               702.292, 758.512, 819.233, 884.814, 955.646, 1032.15, 1114.77, &
               1204.01, 1300.4, 1404.5, 1516.93, 1638.36, 1769.52, 1911.17, 2064.17, &
               2229.41, 2407.88, 2600.64, 2808.82, 3033.67, 3276.53, 3538.82, &
               3822.11, 4128.08, 4458.54, 4815.46, 5200.95, 5617.3, 6066.98, &
               6552.65, 7077.21, 7643.76, 8255.66, 8916.54, 9630.33, 10401.3, &
               11233.9, 12133.2, 13104.5, 14153.5, 15286.6, 16510.3, 17832., &
               19259.5, 20801.2, 22466.4, 24264.9, 26207.4, 28305.4, 30571.3, &
               33018.6, 35661.8, 38516.6, 41599.9, 44930.1, 48526.9, 52411.6, &
               56607.2, 61138.8, 66033.1, 71319.2, 77028.5/)

  piarray = (/ 0.        ,  0.05279988,  0.10559975,  0.15839963,  0.21119951, &
              0.26399938,  0.31679926,  0.36959914,  0.42239901,  0.47519889, &
              0.52799877,  0.58079864,  0.63359852,  0.68639839,  0.73919827, &
              0.79199815,  0.84479802,  0.8975979 ,  0.95039778,  1.00319765, &
              1.05599753,  1.10879741,  1.16159728,  1.21439716,  1.26719704, &
              1.31999691,  1.37279679,  1.42559667,  1.47839654,  1.53119642, &
              1.5839963 ,  1.63679617,  1.68959605,  1.74239593,  1.7951958 , &
              1.84799568,  1.90079556,  1.95359543,  2.00639531,  2.05919518, &
              2.11199506,  2.16479494,  2.21759481,  2.27039469,  2.32319457, &
              2.37599444,  2.42879432,  2.4815942 ,  2.53439407,  2.58719395, &
              2.63999383,  2.6927937 ,  2.74559358,  2.79839346,  2.85119333, &
              2.90399321,  2.95679309,  3.00959296,  3.06239284,  3.11519272, &
              3.16799259,  3.22079247,  3.27359234,  3.32639222,  3.3791921 , &
              3.43199197,  3.48479185,  3.53759173,  3.5903916 ,  3.64319148, &
              3.69599136,  3.74879123,  3.80159111,  3.85439099,  3.90719086, &
              3.95999074,  4.01279062,  4.06559049,  4.11839037,  4.17119025, &
              4.22399012,  4.27679   ,  4.32958988,  4.38238975,  4.43518963, &
              4.48798951,  4.54078938,  4.59358926,  4.64638913,  4.69918901, &
              4.75198889,  4.80478876,  4.85758864,  4.91038852,  4.96318839, &
              5.01598827,  5.06878815,  5.12158802,  5.1743879 ,  5.22718778, &
              5.27998765,  5.33278753,  5.38558741,  5.43838728,  5.49118716, &
              5.54398704,  5.59678691,  5.64958679,  5.70238667,  5.75518654, &
              5.80798642,  5.86078629,  5.91358617,  5.96638605,  6.01918592, &
              6.0719858 ,  6.12478568,  6.17758555,  6.23038543,  6.28318531/)


  phi_init(1) = 18e0_dp

  bg_finish = .false.

  ! IC for the background and background dynamic evolution
  call initialize_background(y)

  do while (.not. bg_finish .and. i .lt. maxstep_bg)
    do ii=1, stepsize
      call gl10_bg( y, Nstep_bg)
    enddo

    phi_array(i,1:1) = real(y(1:1,1))
    dphi_array(i,1:1) = real(y(2:2,1))
    N_array(i,:) = i*stepsize*Nstep_bg
    pot_array(i,:) = pot(real(y(1:1,1)))
    dpot_array(i,:) = dpotdphi(real(y(1:1,1)))
    hubble_array(i,:) = gethubble(real(y(1:1,1)),real(y(2:2,1)))
    eps_array(i,:) = 0.5e0_dp*(Mpl)**2*dot_product(real(y(2:2,1)),real(y(2:2,1)))

    write(97,*) phi_array(i,:), dphi_array(i,:), N_array(i,:), &
                pot_array(i,:), dpot_array(i,:), hubble_array(i,:), &
                eps_array(i,:)

    if (eps_array(i,1).gt.2.0e0_dp) then

      bg_finish = .true.

      kout = i

      ind1 = locate(eps_array(1:kout,1),1.0e0_dp)

      call polint(eps_array(ind1-4:ind1+4,1),N_array(ind1-4:ind1+4,1),1.0e0_dp,N_end,polerr)
      call polint(eps_array(ind1-4:ind1+4,1),phi_array(ind1-4:ind1+4,1),1.0e0_dp,phi_end,polerr)
      call polint(eps_array(ind1-4:ind1+4,1),pot_array(ind1-4:ind1+4,1),1.0e0_dp,pot_end,polerr)

      N_diff = N_end - N_pivot

      ind2 = locate(N_array(1:kout,1),N_diff)

      call polint(N_array(ind2-4:ind2+4,1),hubble_array(ind2-4:ind2+4,1),N_diff,hubble_pivot,polerr)
      call polint(N_array(ind2-4:ind2+4,1),phi_array(ind2-4:ind2+4,1),N_diff,phi_pivot,polerr)
      call polint(N_array(ind2-4:ind2+4,1),dphi_array(ind2-4:ind2+4,1),N_diff,dphi_pivot,polerr)

    end if

    i = i + 1

  end do


  call initialize_background_modes(y,N_start)

  print*, "k_ic = ", k_ic
  print*, " y = ", y


  j = 0
  do while (real(y(1,1)) .gt. phi_pivot .and. j .lt. maxstep_bg)
    do ll=1,kbins
      k = kvector(ll) * kfactor

      if (real(y(1,ll)) .lt. phi_trap) then
        ! print*,real(y(1,ll))
        call calculate_source(y,sint)
      else
        sint = 0.0e0_dp
	  endif

      do jj=1, stepsize
        call gl10_modes( y(:,ll), Nstep )
      enddo
      ! record data near the trap
      if (j*stepsize*Nstep .gt. 4.0e0_dp .and. &
         j*stepsize*Nstep .lt. 6.0e0_dp) then
        write(98,*) k, j*stepsize*Nstep, real(y(:,ll)), aimag(y(:,ll))
      endif
    end do
    j = j + 1
  end do

  contains

    function pot(phi) result(potential)

      real(dp) :: potential
      real(dp), intent(in) :: phi(:)

      potential = 0.5e0_dp*m2*sum(phi*phi)

    end function pot

    function dpotdphi(phi) result(dpotential)

      real(dp), intent(in) :: phi(:)
      real(dp) :: dpotential(size(phi))

      dpotential = m2 * phi

    end function

    function d2potdphi2(phi) result(d2potential)

      real(dp), intent(in) :: phi
      real(dp) :: d2potential

      d2potential = m2

	end function

    function gethubble(phi,dphi) result(gethub)

      real(dp) :: gethub
      real(dp), intent(in) :: phi(:), dphi(:)

      gethub = pot(phi)/3.0e0_dp/Mpl**2 / &
		(1.0e0_dp - dot_product(dphi, dphi)/6.0e0_dp/Mpl**2)

	  gethub = sqrt(gethub)

    end function

	function dhubbledN(phi, dphi) result(getdhubble)

      real(dp) :: getdhubble
      real(dp), intent(in) :: phi(:), dphi(:)

      getdhubble = -dot_product(dphi, dphi) * gethubble(phi,dphi)/2.0e0_dp/Mpl**2

     end function

!
! 10th order Gauss-Legendre integrator
! Used for pieces of Hamiltonian without exact solutions
! To do : implement entire evolution using this
! To do : pack into smaller vector as required (just chop a, so probably not worth it)
!
    subroutine calculate_source(y,totalsource)

      complex(dp), dimension(:,:), intent(in) :: y
      complex(dp), intent(out) :: totalsource

      real(dp) :: kvalue

      complex(dp) :: source

      integer :: a, b, kindex

      source = 0.0e0_dp
      totalsource = 0.0e0_dp
	  ! like to discuss how to write such source terms correctly (will bring up).
      ! Int{2*k'^2\Chi_{\sqrt{k^2-2kk'*cos(theta)+k'^2}}\Chi_{k'}d(k')d(theta)}
      do a=1, kbins
        do b=1, pibins
          kvalue = sqrt(k*k-2.0e0_dp*k*kvector(a)*kfactor*&
                        cos(piarray(b))+kvector(a)**2*kfactor**2)

          kindex = locate(kvector,kvalue)

          source = source + kvector(a)**2*y(5,a)*y(5,kindex)*pibin
        end do
        totalsource = totalsource + source*kbinwidths(a)
      end do

    end subroutine calculate_source

    subroutine initialize_background(y)

      complex(dp), dimension(nvar,kbins) :: y

      real(dp), dimension(size(y(:,1))) :: dVi
      real(dp) :: phi, dphi, hubble, dhubble, dV1

      real(dp) :: phi_bgstart, eps

      integer :: ll

      ! set H_init from slow roll

      dVi = dpotdphi(phi_init)
	  dV1 = dVi(1)

      hubble_init = pot(phi_init)/6e0_dp/Mpl**2 *(1e0_dp+sqrt(1e0_dp+&
                  2e0_dp/3e0_dp*Mpl**2*dV1**2/pot(phi_init)**2))

      hubble_init = sqrt(hubble_init)

      dphi_init = -dpotdphi(phi_init)/(3e0_dp*hubble_init**2)

      eps = 0.5e0_dp*(Mpl)**2 * sum(dphi_init*dphi_init)

      do ll=1,kbins
        y(1:1,ll) = cmplx(phi_init(1:1),kind=dp)
        y(2:2,ll) = cmplx(dphi_init(1:1),kind=dp)
      end do

	end subroutine initialize_background

    subroutine initialize_background_modes(y,N_ic)
      ! check for how deep into horizon (in terms of e-folds/scale-factor) BD is satisfied!
	  ! set those conditions as the IC

      complex(dp), dimension(nvar,kbins) :: y

      real(dp), intent(out) :: N_ic

      real(dp) :: push_horizon, cond
      real(dp) :: check1, check2, check3
      real(dp) :: lnaH, lnaH_pivot

      integer :: ll

      logical :: bd_satisfied

      do ll=1, kbins

        push_horizon = 1e0_dp

        bd_satisfied = .false.

        k = kvector(ll) * kfactor
        a_init(ll) = k / hubble_pivot / exp(N_end-N_pivot)
        a_end(ll) = a_init(ll) * exp(N_end)

        do i=1,kout
          lnaH_array(i,ll)=log(a_init(ll)*exp(N_array(i,1))*hubble_array(i,1))
        end do

        ! start at least as deep as beyond the trap
        ind1 = locate(phi_array(1:kout,1),phi_trap)

        call polint(phi_array(ind1-4:ind1+4,1),lnaH_array(ind1-4:ind1+4,ll), &
                                 phi_trap, lnaH,polerr)

        ind2 = locate(phi_array(1:kout,1),phi_pivot)

        ! locate the pivot value to scale with k_piv
        call polint(phi_array(ind2-4:ind2+4,1),lnaH_array(ind2-4:ind2+4,ll), &
                    phi_pivot, lnaH_pivot,polerr)

        push_horizon = exp(lnaH_pivot) / exp(lnaH)

        ! print*,"lnaH_array=",lnaH_array(1:kout,1)

        do  while (.not. bd_satisfied)

          push_horizon = push_horizon*2.0e0_dp

          lnaH = lnaH_pivot - log(push_horizon)

          ind2 = locate(lnaH_array(1:kout,ll),lnaH)

          call polint(lnaH_array(ind2-4:ind2+4,ll), phi_array(ind2-4:ind2+4,1), &
                                                             lnaH, phi_ic, polerr)
          call polint(lnaH_array(ind2-4:ind2+4,ll), hubble_array(ind2-4:ind2+4,1), &
													lnaH, hubble_ic(ll), polerr)
          call polint(lnaH_array(ind2-4:ind2+4,ll), dphi_array(ind2-4:ind2+4,1),&
                                                    lnaH,dphi_ic,polerr)
          call polint(lnaH_array(ind2-4:ind2+4,ll), eps_array(ind2-4:ind2+4,1), &
 													lnaH, eps_ic, polerr)
          call polint(lnaH_array(ind2-4:ind2+4,ll), pot_array(ind2-4:ind2+4,1),&
													lnaH, pot_ic, polerr)
          call polint(lnaH_array(ind2-4:ind2+4,ll), N_array(ind2-4:ind2+4,1), &
												    lnaH, N_ic, polerr)

          check1 = abs((eps_ic-2.0e0_dp)/(push_horizon**2))
          check2 = abs(d2potdphi2(phi_ic)/(push_horizon**2 * hubble_ic(ll)**2))
          check3 = g2*(phi_ic - phi_trap)**2/(hubble_ic(ll)**2*push_horizon**2)

          cond = 1e-3_dp

	      if (check1 < cond .and. &
			  check2 < cond .and. &
              check3 < cond) then

            bd_satisfied = .true.

            print*,"BD satisfied, push horizon factor = ",push_horizon

          end if

          k_ic(ll) = push_horizon
          lnaH_ic(ll) = lnaH

	    end do

	    ! inflaton background
	    y(1:1,ll) = cmplx(phi_ic,kind=dp)
	    y(2:2,ll) = cmplx(dphi_ic,kind=dp)

	    ! initializing the modes (inflaton)
	    y(3:3,ll) = cmplx(1.0e0_dp/sqrt(2.0e0_dp*k),kind=dp)
	    y(4:4,ll) = cmplx(0.0e0_dp,-k/exp(lnaH_ic(ll))/sqrt(2.0e0_dp*k),kind=dp)

	    ! initializing the modes (matter field)
	    y(5:5,ll)  = cmplx(1.0e0_dp/sqrt(2.0e0_dp*k),kind=dp)
	    y(6:6,ll) = cmplx(0.0e0_dp,-k/exp(lnaH_ic(ll))/sqrt(2.0e0_dp*k),kind=dp)

      end do

    end subroutine initialize_background_modes


    subroutine gl10_bg( y_comp, dt )

      complex(dp), dimension(nvar,kbins) :: y_comp
      real(dp) :: y_bg(2)
      real(dp) :: dt

      integer, parameter :: s = 5
      integer, parameter :: niter = 8
      real*8 :: g(2,s)

      ! Butcher tableau for 10th order Gauss-Legendre method
      real*8, parameter :: a(s,s) = reshape( (/ &
           0.5923172126404727187856601017997934066D-1, -1.9570364359076037492643214050884060018D-2, &
           1.1254400818642955552716244215090748773D-2, -0.5593793660812184876817721964475928216D-2, &
           1.5881129678659985393652424705934162371D-3,  1.2815100567004528349616684832951382219D-1, &
           1.1965716762484161701032287870890954823D-1, -2.4592114619642200389318251686004016630D-2, &
           1.0318280670683357408953945056355839486D-2, -2.7689943987696030442826307588795957613D-3, &
           1.1377628800422460252874127381536557686D-1,  2.6000465168064151859240589518757397939D-1, &
           1.4222222222222222222222222222222222222D-1, -2.0690316430958284571760137769754882933D-2, &
           4.6871545238699412283907465445931044619D-3,  1.2123243692686414680141465111883827708D-1, &
           2.2899605457899987661169181236146325697D-1,  3.0903655906408664483376269613044846112D-1, &
           1.1965716762484161701032287870890954823D-1, -0.9687563141950739739034827969555140871D-2, &
           1.1687532956022854521776677788936526508D-1,  2.4490812891049541889746347938229502468D-1, &
           2.7319004362580148889172820022935369566D-1,  2.5888469960875927151328897146870315648D-1, &
           0.5923172126404727187856601017997934066D-1 /) , [s,s])
      real, parameter :: b(s) = (/ &
           1.1846344252809454375713202035995868132D-1,  2.3931433524968323402064575741781909646D-1, &
           2.8444444444444444444444444444444444444D-1,  2.3931433524968323402064575741781909646D-1, &
           1.1846344252809454375713202035995868132D-1 /)

      integer :: i,kk

      y_bg(1) = real(y_comp(1,1),kind=dp)
      y_bg(2) = real(y_comp(2,1),kind=dp)

      g = 0.
      do kk=1,niter
        g = matmul(g,a)
		do i=1,s
		  call derivs_bg( y_bg+g(:,i)*dt, g(:,i) )
		  enddo
		enddo
	  y_bg = y_bg + matmul(g,b)*dt

      y_comp(1,1) = cmplx(y_bg(1),kind=dp)
	  y_comp(2,1) = cmplx(y_bg(2),kind=dp)

	end subroutine gl10_bg



    subroutine gl10_modes( ycomplex, dt )

      complex(dp), dimension(nvar) :: ycomplex
	  real(dp), dimension(2*nvar) :: yreal
      real*8 :: dt

      integer, parameter :: s = 5
      integer, parameter :: niter = 8
      real*8 :: g(2*nvar,s)

      ! Butcher tableau for 10th order Gauss-Legendre method
      real*8, parameter :: a(s,s) = reshape( (/ &
           0.5923172126404727187856601017997934066D-1, -1.9570364359076037492643214050884060018D-2, &
           1.1254400818642955552716244215090748773D-2, -0.5593793660812184876817721964475928216D-2, &
           1.5881129678659985393652424705934162371D-3,  1.2815100567004528349616684832951382219D-1, &
           1.1965716762484161701032287870890954823D-1, -2.4592114619642200389318251686004016630D-2, &
           1.0318280670683357408953945056355839486D-2, -2.7689943987696030442826307588795957613D-3, &
           1.1377628800422460252874127381536557686D-1,  2.6000465168064151859240589518757397939D-1, &
           1.4222222222222222222222222222222222222D-1, -2.0690316430958284571760137769754882933D-2, &
           4.6871545238699412283907465445931044619D-3,  1.2123243692686414680141465111883827708D-1, &
           2.2899605457899987661169181236146325697D-1,  3.0903655906408664483376269613044846112D-1, &
           1.1965716762484161701032287870890954823D-1, -0.9687563141950739739034827969555140871D-2, &
           1.1687532956022854521776677788936526508D-1,  2.4490812891049541889746347938229502468D-1, &
           2.7319004362580148889172820022935369566D-1,  2.5888469960875927151328897146870315648D-1, &
           0.5923172126404727187856601017997934066D-1 /) , [s,s])
      real, parameter :: b(s) = (/ &
           1.1846344252809454375713202035995868132D-1,  2.3931433524968323402064575741781909646D-1, &
           2.8444444444444444444444444444444444444D-1,  2.3931433524968323402064575741781909646D-1, &
           1.1846344252809454375713202035995868132D-1 /)

      integer :: i,kk

      call complextoreal(yreal,ycomplex)

      g = 0.
      do kk=1,niter
         g = matmul(g,a)
         do i=1,s
            call derivs_modes( yreal+g(:,i)*dt, g(:,i) )
         enddo
      enddo
      yreal = yreal + matmul(g,b)*dt

      call realtocomplex(yreal,ycomplex)

    end subroutine gl10_modes

    subroutine derivs_bg(y, yprime)

      real(kind=dp), dimension(:), intent(in) :: y
      real(kind=dp), dimension(:), intent(out) :: yprime

      real(dp), dimension(1) :: dVi
      real(dp) :: phi, dphi, hubble, dhubble, dV1

      phi = y(1)
      dphi = y(2)
      hubble = gethubble(y(1:1),y(2:2))
      dhubble = dhubbledN(y(1:1),y(2:2))
      dVi = dpotdphi(y(1:1))
      dV1 = dVi(1)

      yprime(1) = dphi

      yprime(2) = -((3e0_dp+dhubble/hubble)*dphi+dV1/hubble/hubble)

    end subroutine derivs_bg

!
! Evolution vector in phase space for GL integration
!
    subroutine derivs_modes(y, yprime)

      real(kind=dp), dimension(:), intent(in) :: y
      real(kind=dp), dimension(:), intent(out) :: yprime


      complex(kind=dp), dimension(size(y)/2) :: ycomplex
      complex(kind=dp), dimension(size(y)/2) :: yprime_complex

      real(dp), dimension(1) :: dVi, phi, dphi
      real(dp) :: hubble, dhubble, dV1, eps, scale_factor

      call realtocomplex(y,ycomplex)

      phi = real(ycomplex(1:1))
      dphi = real(ycomplex(2:2))

      hubble = gethubble(phi,dphi)
      dhubble = dhubbledN(phi,dphi)
      dVi = dpotdphi(phi)
      dV1 = dVi(1)
      eps = 0.5e0_dp*(Mpl)**2*sum(dphi*dphi)

      ! background evolution
      yprime_complex(1:1) = ycomplex(2:2)
      yprime_complex(2:2) = cmplx(-((3e0_dp+dhubble/hubble)*dphi+dV1/hubble/hubble),kind=dp)

	  ! modes evolution (inflaton)
      scale_factor = exp(j*stepsize*Nstep)/(hubble_ic(ll)/exp(lnaH_ic(ll)))

      yprime_complex(3:3) = ycomplex(4:4)
      yprime_complex(4:4) = -(1.0e0_dp - eps)*ycomplex(4:4) &
						- (k/scale_factor/hubble)**2*ycomplex(3:3) &
                        + (2.0e0_dp - eps)*ycomplex(3:3) - &
                           d2potdphi2(phi(1))*ycomplex(3:3)/hubble**2 &
						- (g2*(phi-phi_trap)/hubble**2)*sint


      yprime_complex(5:5) = ycomplex(6:6)
      yprime_complex(6:6) = -(1.0e0_dp - eps)*ycomplex(6:6) &
                            - (k/scale_factor/hubble)**2*ycomplex(5:5) &
				            + (2.0e0_dp - eps)*ycomplex(5:5) - &
						   (g2*(phi-phi_trap)**2/hubble**2)*ycomplex(5:5)

      call complextoreal(yprime,yprime_complex)

    end subroutine derivs_modes

    subroutine realtocomplex(yreal,ycomplex)

      real(kind=dp), dimension(:), intent(in) :: yreal
      complex(kind=dp), dimension(:), intent(out) :: ycomplex

      ycomplex(1) = cmplx(yreal(1),yreal(2),kind=dp)
      ycomplex(2) = cmplx(yreal(3),yreal(4),kind=dp)
      ycomplex(3) = cmplx(yreal(5),yreal(6),kind=dp)
      ycomplex(4) = cmplx(yreal(7),yreal(8),kind=dp)
      ycomplex(5) = cmplx(yreal(9),yreal(10),kind=dp)
      ycomplex(6) = cmplx(yreal(11),yreal(12),kind=dp)

    end subroutine realtocomplex

    subroutine complextoreal(yreal,ycomplex)

      real(kind=dp), dimension(:), intent(out) :: yreal
      complex(kind=dp), dimension(:), intent(in) :: ycomplex

      yreal(1) = real(ycomplex(1),kind=dp)
      yreal(2) = aimag(ycomplex(1))
      yreal(3) = real(ycomplex(2),kind=dp)
      yreal(4) = aimag(ycomplex(2))
      yreal(5) = real(ycomplex(3),kind=dp)
      yreal(6) = aimag(ycomplex(3))
      yreal(7) = real(ycomplex(4),kind=dp)
      yreal(8) = aimag(ycomplex(4))
      yreal(9) = real(ycomplex(5),kind=dp)
      yreal(10) = aimag(ycomplex(5))
      yreal(11) = real(ycomplex(6),kind=dp)
      yreal(12) = aimag(ycomplex(6))

	end subroutine complextoreal

    ! Polynomial interpolation
    ! Given array XA and YA (of same length) and given a value X, returns
    !value Y such that if P(x) is a polynomial st P(XA_i)=YA_i, then Y=P(X) ---
    !and an error estimate DY
    SUBROUTINE polint(xa,ya,x,y,dy)
      IMPLICIT NONE
      real(dp), DIMENSION(:), INTENT(IN) :: xa,ya
      real(dp), INTENT(IN) :: x
      real(dp), INTENT(OUT) :: y,dy
      INTEGER*4 :: m,n,ns
      INTEGER*4, DIMENSION(1) :: imin
      real(dp), DIMENSION(size(xa)) :: c,d,den,ho,absho

      if (size(xa)==size(ya)) then
        n=size(xa)
      else
        print*,"error in polint (wrong array size)"
      !  call raise%fatal_code(&
      !    'Wrong array sizes in polint',&
      !    __FILE__, __LINE__)
	  end if

      c=ya
      d=ya
      ho=xa-x
      absho=abs(ho)
      imin=minloc(absho(:))
      ns=imin(1)
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
        den(1:n-m)=ho(1:n-m)-ho(1+m:n)
        if (any(den(1:n-m) == 0.0)) then
          print*,"calculation failure in polint"
          ! call raise%fatal_code(&
          !   'Calculation failure in polint.',&
          !   __FILE__, __LINE__)
		end if
		den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
		d(1:n-m)=ho(1+m:n)*den(1:n-m)
		c(1:n-m)=ho(1:n-m)*den(1:n-m)
		if (2*ns < n-m) then
		  dy=c(ns+1)
		else
		  dy=d(ns)
		  ns=ns-1
		end if
		y=y+dy
	  end do
    !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
    END SUBROUTINE polint

    pure FUNCTION locate(xx,x)
      IMPLICIT NONE
      real(dp), DIMENSION(:), INTENT(IN) :: xx
      real(dp), INTENT(IN) :: x
      INTEGER*4 :: locate
      INTEGER*4 :: n,jl,jm,ju
      LOGICAL :: ascnd
      n=size(xx)
      ascnd = (xx(n) >= xx(1))
      jl=0
      ju=n+1
      do
        if (ju-jl <= 1) exit
          jm=(ju+jl)/2
        if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
		else
          ju=jm
        end if
      end do
      if (x == xx(1)) then
        locate=1
	  else if (x == xx(n)) then
        locate=n-1
      else
        locate=jl
      end if
      !  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
    END FUNCTION locate


  end program evolve_trapped
