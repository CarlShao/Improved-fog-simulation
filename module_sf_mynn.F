!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_sf_mynn

!-------------------------------------------------------------------
!Modifications implemented by Joseph Olson NOAA/GSD/AMB - CU/CIRES
!for WRFv3.4, v3.4.1, v3.5.1, v3.6, v3.7.1, and v3.9:
!
!   BOTH LAND AND WATER:
!1) Calculation of stability parameter (z/L) taken from Li et al. (2010 BLM)
!   for first iteration of first time step; afterwards, exact calculation.
!2) Fixed isflux=0 option to turn off scalar fluxes, but keep momentum
!   fluxes for idealized studies (credit: Anna Fitch).
!3) Kinematic viscosity now varies with temperature
!4) Uses Monin-Obukhov flux-profile relationships more consistent with
!   those used in the MYNN PBL code.
!5) Allows negative QFX, similar to MYJ scheme
!
!   LAND only:
!1) iz0tlnd option is now available with the following options:
!   (default) =0: Zilitinkevich (1995) 
!             =1: Czil_new (modified according to Chen & Zhang 2008)
!             =2: Modified Yang et al (2002, 2008) - generalized for all landuse
!             =3: constant zt = z0/7.4 (original form; Garratt 1992)
!             =4: Pan et al. (1994) with RUC mods for z_q, zili for z_t
!2) Relaxed u* minimum from 0.1 to 0.01
!
!   WATER only:
!1) isftcflx option is now available with the following options:
!   (default) =0: z0, zt, and zq from the COARE algorithm. Set COARE_OPT (below) to
!                 3.0 (Fairall et al. 2003, default)
!                 3.5 (Edson et al 2013) 
!             =1: z0 from Davis et al (2008), zt & zq from COARE 3.0/3.5
!             =2: z0 from Davis et al (2008), zt & zq from Garratt (1992)
!             =3: z0 from Taylor and Yelland (2004), zt and zq from COARE 3.0/3.5
!             =4: z0 from Zilitinkevich (2001), zt & zq from COARE 3.0/3.5
!
!   SNOW/ICE only:
!1) Added Andreas (2002) snow/ice parameterization for thermal and
!   moisture roughness to help reduce the cool/moist bias in the arctic
!   region. Also added a z0 mod for snow (Andreas et al. 2005, BLM), which
!
! Misc:
!   2) added a more elaborate diagnostic for u10 & V10 for high vertical resolution
!      model configurations.
!
! New for v3.9:
!   - option for stochastic parameter perturbations (SPP) 
!
!NOTE: This code was primarily tested in combination with the RUC LSM.
!      Performance with the Noah (or other) LSM is relatively unknown.
!-------------------------------------------------------------------
!For WRF
  USE module_model_constants, only: &
       &g, p1000mb, cp, xlv, ep_2, r_d, r_v, rcp, cpv 

  USE module_bl_mynn,   only: tv0, b1, b2, p608, ev, rd, & !, mym_condensation
       &esat_blend, xl_blend, qsat_blend
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!For non-WRF
!   REAL    , PARAMETER :: g            = 9.81
!   REAL    , PARAMETER :: r_d          = 287.
!   REAL    , PARAMETER :: cp           = 7.*r_d/2.
!   REAL    , PARAMETER :: r_v          = 461.6
!   REAL    , PARAMETER :: cpv          = 4.*r_v
!   REAL    , PARAMETER :: rcp          = r_d/cp
!   REAL    , PARAMETER :: XLV          = 2.5E6
!   REAL    , PARAMETER :: XLF          = 3.50E5
!   REAL    , PARAMETER :: p1000mb      = 100000.
!   REAL    , PARAMETER :: EP_2         = r_d/r_v

  REAL, PARAMETER :: xlvcp=xlv/cp, ep_3=1.-ep_2 
  REAL, PARAMETER :: wmin=0.1    ! Minimum wind speed
  REAL, PARAMETER :: VCONVC=1.25
  REAL, PARAMETER :: SNOWZ0=0.011
  REAL, PARAMETER :: COARE_OPT=3.0  ! 3.0 or 3.5
  REAL, PARAMETER :: am = 5.0
  REAL, PARAMETER :: bm = am/6.5
  REAL, PARAMETER :: B_m = ((1.0-bm)/bm)**(1.0/3.0)
  REAL, PARAMETER :: ah = 5.0
  REAL, PARAMETER :: bh = 5.0
  REAL, PARAMETER :: B_h = 5.0**0.5
  REAL, PARAMETER :: C_h = 3.0
  REAL, PARAMETER :: cm_a = -3.0*am/bm
  REAL, PARAMETER :: cm_b = am*B_m/bm/2.0
  REAL, PARAMETER :: cm_c = 1.0-B_m+B_m**2.0
  REAL, PARAMETER :: cm_d = (2.0-B_m)/B_m/(3.0**0.5)
  REAL, PARAMETER :: ch_a = 5.0**0.5*0.5
  REAL, PARAMETER :: ch_b = C_h-B_h
  REAL, PARAMETER :: ch_c = C_h+B_h
  REAL, PARAMETER :: ch_d = ch_b/ch_c
  !For debugging purposes:
  LOGICAL, PARAMETER :: debug_code = .false.

CONTAINS

!-------------------------------------------------------------------
  SUBROUTINE mynn_sf_init_driver(allowed_to_read)

    LOGICAL, INTENT(in) :: allowed_to_read

    !Fill the PSIM and PSIH tables. The subroutine "sfclayinit" 
    !can be found in module_sf_sfclay.F. This subroutine returns
    !the forms from Dyer and Hicks (1974).
       
!    CALL sfclayinit(allowed_to_read)

  END SUBROUTINE mynn_sf_init_driver

!-------------------------------------------------------------------
   SUBROUTINE SFCLAY_mynn(                                         &
                     U3D,V3D,T3D,QV3D,P3D,dz8w,                    &
                     CP,G,ROVCP,R,XLV,PSFCPA,CHS,CHS2,CQS2,CPM,    &
                     ZNT,UST,PBLH,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
                     XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
                     U10,V10,TH2,T2,Q2,SNOWH,                      &
                     GZ1OZ0,WSPD,BR,ISFFLX,DX,                     &
                     SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                     KARMAN,itimestep,ch,th3d,pi3d,qc3d,rho3d,     &
                     tsq,qsq,cov,sh3d,el_pbl,qcg,                  &
                     icloud_bl,qc_bl,cldfra_bl,                    &
                     spp_pbl,pattern_spp_pbl,                      &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte,                    &
                     ustm,ck,cka,cd,cda,isftcflx,iz0tlnd,          &
                     bl_mynn_cloudpdf,z0zt,zaz0,za_sf,psim_sf,psih_sf, &
                     dth_sf, dqv_sf, minpsi_m, minpsi_h, li_opt, cb05_opt)
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
!-- U3D         3D u-velocity interpolated to theta points (m/s)
!-- V3D         3D v-velocity interpolated to theta points (m/s)
!-- T3D         3D temperature (K)
!-- QV3D        3D water vapor mixing ratio (Kg/Kg)
!-- P3D         3D pressure (Pa)
!-- RHO3D       3D density (kg/m3) 
!-- dz8w        3D dz between full levels (m)
!-- CP          heat capacity at constant pressure for dry air (J/kg/K)
!-- G           acceleration due to gravity (m/s^2)
!-- ROVCP       R/CP
!-- R           gas constant for dry air (J/kg/K)
!-- XLV         latent heat of vaporization for water (J/kg)
!-- PSFCPA      surface pressure (Pa)
!-- ZNT         roughness length (m)
!-- UST         u* in similarity theory (m/s)
!-- USTM        u* in similarity theory (m/s) w* added to WSPD. This is
!               used to couple with TKE scheme but not in MYNN.
!               (as of now, USTM = UST in this version)
!-- PBLH        PBL height from previous time (m)
!-- MAVAIL      surface moisture availability (between 0 and 1)
!-- ZOL         z/L height over Monin-Obukhov length
!-- MOL         T* (similarity theory) (K)
!-- RMOL        Reciprocal of M-O length (/m)
!-- REGIME      flag indicating PBL regime (stable, unstable, etc.)
!-- PSIM        similarity stability function for momentum
!-- PSIH        similarity stability function for heat
!-- XLAND       land mask (1 for land, 2 for water)
!-- HFX         upward heat flux at the surface (W/m^2)
!-- QFX         upward moisture flux at the surface (kg/m^2/s)
!-- LH          net upward latent heat flux at surface (W/m^2)
!-- TSK         surface temperature (K)
!-- FLHC        exchange coefficient for heat (W/m^2/K)
!-- FLQC        exchange coefficient for moisture (kg/m^2/s)
!-- CHS         heat/moisture exchange coefficient for LSM (m/s)
!-- QGH         lowest-level saturated mixing ratio
!-- QSFC        qv (specific humidity) at the surface
!-- QSFCMR      qv (mixing ratio) at the surface
!-- U10         diagnostic 10m u wind
!-- V10         diagnostic 10m v wind
!-- TH2         diagnostic 2m theta (K)
!-- T2          diagnostic 2m temperature (K)
!-- Q2          diagnostic 2m mixing ratio (kg/kg)
!-- SNOWH       Snow height (m)
!-- GZ1OZ0      log((z1+ZNT)/ZNT) where ZNT is roughness length 
!-- WSPD        wind speed at lowest model level (m/s)
!-- BR          bulk Richardson number in surface layer
!-- ISFFLX      isfflx=1 for surface heat and moisture fluxes
!-- DX          horizontal grid size (m)
!-- SVP1        constant for saturation vapor pressure (=0.6112 kPa)
!-- SVP2        constant for saturation vapor pressure (=17.67 dimensionless)
!-- SVP3        constant for saturation vapor pressure (=29.65 K)
!-- SVPT0       constant for saturation vapor pressure (=273.15 K)
!-- EP1         constant for virtual temperature (Rv/Rd - 1) (dimensionless)
!-- EP2         constant for spec. hum. calc (Rd/Rv = 0.622) (dimensionless)
!-- EP3         constant for spec. hum. calc (1 - Rd/Rv = 0.378 ) (dimensionless)
!-- KARMAN      Von Karman constant
!-- ck          enthalpy exchange coeff at 10 meters
!-- cd          momentum exchange coeff at 10 meters
!-- cka         enthalpy exchange coeff at the lowest model level
!-- cda         momentum exchange coeff at the lowest model level
!-- isftcflx    =0: z0, zt, and zq from COARE3.0/3.5 (Fairall et al 2003/Edson et al 2013)
!   (water      =1: z0 from Davis et al (2008), zt & zq from COARE3.0/3.5
!    only)      =2: z0 from Davis et al (2008), zt & zq from Garratt (1992)
!               =3: z0 from Taylor and Yelland (2004), zt and zq from COARE 3.0/3.5
!               =4: z0 from Zilitinkevich (2001), zt & zq from COARE 3.0/3.5
!-- iz0tlnd     =0: Zilitinkevich (1995) with Czil=0.10, 
!   (land       =1: Czil_new (modified according to Chen & Zhang 2008)
!    only)      =2: Modified Yang et al (2002, 2008) - generalized for all landuse
!               =3: constant zt = z0/7.4 (Garratt 1992)
!               =4: Pan et al (1994) for zq; ZIlitintevich for zt
!-- bl_mynn_cloudpdf =0: Mellor & Yamada
!                    =1: Kuwano et al.
!-- el_pbl      = mixing length from PBL scheme (meters)
!-- Sh3d        = Stability finction for heat (unitless)
!-- cov         = T'q' from PBL scheme
!-- tsq         = T'T' from PBL scheme
!-- qsq         = q'q' from PBL scheme
!-- icloud_bl  = namelist option for subgrid scale cloud/radiation feedback
!-- qc_bl       = subgrid scale (bloundary layer) clouds
!-- cldfra_bl   = subgridscale cloud fraction
!
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!=================================================================
! SCALARS
!===================================
      INTEGER,  INTENT(IN)   ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte
      INTEGER,  INTENT(IN)   ::        itimestep
      REAL,     INTENT(IN)   ::   SVP1,SVP2,SVP3,SVPT0,minpsi_m,minpsi_h
      REAL,     INTENT(IN)   ::        EP1,EP2,KARMAN
      REAL,     INTENT(IN)   ::        CP,G,ROVCP,R,XLV,DX
!NAMELIST OPTIONS:
      INTEGER,  INTENT(IN)   ::        ISFFLX
      INTEGER,  INTENT(IN)   ::        li_opt, cb05_opt
      INTEGER,  OPTIONAL,  INTENT(IN)   ::     ISFTCFLX, IZ0TLND,&
                                                bl_mynn_cloudpdf,&
                                                icloud_bl
      INTEGER,  INTENT(IN),OPTIONAL    ::    spp_pbl

!===================================
! 3D VARIABLES
!===================================
      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                           dz8w, &
                                                             QV3D, &
                                                              P3D, &
                                                              T3D, &
                                                             QC3D, &
                                                          U3D,V3D, &
                          RHO3D,th3d,pi3d,tsq,qsq,cov,sh3d,el_pbl

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme ) ::   qc_bl, &
                                                        cldfra_bl
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN),OPTIONAL ::pattern_spp_pbl
!===================================
! 2D VARIABLES
!===================================
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN   )               ::             MAVAIL, &
                                                             PBLH, &
                                                            XLAND, &
                                                              TSK, &
                                                              QCG, &
                                                           PSFCPA ,&
                                                            SNOWH


      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(OUT  )               ::            U10,V10, &
                                                        TH2,T2,Q2

      REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme )              , &
                INTENT(OUT)     ::              ck,cka,cd,cda,ustm
!
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)               ::             REGIME, &
                                                              HFX, &
                                                              QFX, &
                                                               LH, &
                                                         MOL,RMOL, &
                                                        QSFC, QGH, &
                                                              ZNT, &
                                                              ZOL, &
                                                              UST, &
                                                              CPM, &
                                                             CHS2, &
                                                             CQS2, &
                                                              CHS, &
                                                               CH, &
                                                        FLHC,FLQC, &
                                                   GZ1OZ0,WSPD,BR, &
                                                        PSIM,PSIH, &
                                                z0zt, zaz0, za_sf, &
                                     psim_sf,psih_sf,dth_sf,dqv_sf

!ADDITIONAL OUTPUT
!JOE-begin
      REAL,     DIMENSION( ims:ime, jms:jme )    ::     &
                                 BulkRi,wstar,qstar,resist,logres
!JOE-end 
!===================================
! 1D LOCAL ARRAYS
!===================================
      REAL,     DIMENSION( its:ite ) ::                       U1D, &
                                                              V1D, &
                                                        U1D2,V1D2, & !level2 winds
                                                             QV1D, &
                                                              P1D, &
                                                         T1D,QC1D, &
                                                            RHO1D, &
                                                           dz8w1d, & !level 1 height
                                                           dz2w1d    !level 2 height

      REAL,     DIMENSION( its:ite ) ::                  rstoch1D

      ! VARIABLE FOR PASSING TO MYM_CONDENSATION
      REAL,     DIMENSION(kts:kts+1 ) :: dummy1,dummy2,dummy3,dummy4, &
                                         dummy5,dummy6,dummy7,dummy8, &
                                         dummy9,dummy10,dummy11,      &
                                         dummy12,dummy13,dummy14

      REAL,     DIMENSION( its:ite ) ::  vt1,vq1
      REAL,     DIMENSION(kts:kts+1) ::  thl, qw, vt, vq
      REAL                           ::  ql

      INTEGER ::  I,J,K,itf,jtf,ktf
!-----------------------------------------------------------

      itf=MIN0(ite,ide-1)
      jtf=MIN0(jte,jde-1)
      ktf=MIN0(kte,kde-1)

      DO J=jts,jte
        DO i=its,ite
           dz8w1d(I) = dz8w(i,kts,j)
           dz2w1d(I) = dz8w(i,kts+1,j)
           U1D(i) =U3D(i,kts,j)
           V1D(i) =V3D(i,kts,j)
           !2nd model level winds - for diags with high-res grids
           U1D2(i) =U3D(i,kts+1,j)
           V1D2(i) =V3D(i,kts+1,j)
           QV1D(i)=QV3D(i,kts,j)
           QC1D(i)=QC3D(i,kts,j)
           P1D(i) =P3D(i,kts,j)
           T1D(i) =T3D(i,kts,j)
           RHO1D(i)=RHO3D(i,kts,j)
           if (spp_pbl==1) then
               rstoch1D(i)=pattern_spp_pbl(i,kts,j)
           else
               rstoch1D(i)=0.0
           endif
        ENDDO

        IF (itimestep==1) THEN
           DO i=its,ite
              vt1(i)=0.
              vq1(i)=0.
              UST(i,j)=MAX(0.025*SQRT(U1D(i)*U1D(i) + V1D(i)*V1D(i)),0.001)
              MOL(i,j)=0.     ! Tstar
              QSFC(i,j)=QV3D(i,kts,j)/(1.+QV3D(i,kts,j))
              qstar(i,j)=0.0
           ENDDO
        ELSE
           DO i=its,ite
              DO k = kts,kts+1
                ql = qc3d(i,k,j)/(1.+qc3d(i,k,j))
                qw(k) = qv3d(i,k,j)/(1.+qv3d(i,k,j)) + ql
                thl(k) = th3d(i,k,j)-xlvcp*ql/pi3d(i,k,j)
                dummy1(k)=dz8w(i,k,j)
                dummy2(k)=thl(k)
                dummy3(k)=qw(k)
                dummy4(k)=p3d(i,k,j)
                dummy5(k)=pi3d(i,k,j)
                dummy6(k)=tsq(i,k,j)
                dummy7(k)=qsq(i,k,j)
                dummy8(k)=cov(i,k,j)
                dummy9(k)=Sh3d(i,k,j)
                dummy10(k)=el_pbl(i,k,j)
                dummy14(k)=th3d(i,k,j)
                if(icloud_bl > 0) then
                  dummy11(k)=qc_bl(i,k,j)
                  dummy12(k)=cldfra_bl(i,k,j)
                else
                  dummy11(k)=0.0
                  dummy12(k)=0.0
                endif
                dummy13(k)=0.0     !sgm
              ENDDO

              ! NOTE: The last grid number is kts+1 instead of kte.
              CALL mym_condensation (kts,kts+1, dx,     &
                   &            dummy1,dummy2,dummy3,   &
                   &            dummy4,dummy5,dummy6,   &
                   &            dummy7,dummy8,dummy9,   &
                   &            dummy10,bl_mynn_cloudpdf,&
                   &            dummy11,dummy12,        &
                   &            PBLH(i,j),HFX(i,j),     &
                   &            vt(kts:kts+1), vq(kts:kts+1), &
                   &            dummy14,dummy13)

!              ! NOTE: The last grid number is kts+1 instead of kte.
!              CALL mym_condensation (kts,kts+1, dx,     &
!                   &            dz8w(i,kts:kts+1,j),    &
!                   &            thl(kts:kts+1),         &
!                   &            qw(kts:kts+1),          &
!                   &            p3d(i,kts:kts+1,j),     &
!                   &            pi3d(i,kts:kts+1,j),    &
!                   &            tsq(i,kts:kts+1,j),     &
!                   &            qsq(i,kts:kts+1,j),     &
!                   &            cov(i,kts:kts+1,j),     &
!                   &            Sh3d(i,kts:kts+1,j),    & !JOE - cloud PDF testing
!                   &            el_pbl(i,kts:kts+1,j),  & !JOE - cloud PDF testing
!                   &            bl_mynn_cloudpdf,       & !JOE - cloud PDF testing
!                   &            qc_bl2D(i,kts:kts+1),   & !JOE-subgrid BL clouds
!                   &           cldfra_bl2D(i,kts:kts+1),& !JOE-subgrid BL clouds
!                   &            PBLH(i,j),HFX(i,j),     & !JOE-subgrid BL clouds
!                   &            vt(kts:kts+1), vq(kts:kts+1), &
 !                  &            th,sgm)
              vt1(i) = vt(kts)
              vq1(i) = vq(kts)
           ENDDO
        ENDIF

        CALL SFCLAY1D_mynn(                                        &
                J,U1D,V1D,T1D,QV1D,P1D,dz8w1d,rho1d,               &
                U1D2,V1D2,dz2w1d,                                  &
                CP,G,ROVCP,R,XLV,PSFCPA(ims,j),CHS(ims,j),CHS2(ims,j),&
                CQS2(ims,j),CPM(ims,j),PBLH(ims,j), RMOL(ims,j),   &
                ZNT(ims,j),UST(ims,j),MAVAIL(ims,j),ZOL(ims,j),    &
                MOL(ims,j),REGIME(ims,j),PSIM(ims,j),PSIH(ims,j),  &
                XLAND(ims,j),HFX(ims,j),QFX(ims,j),TSK(ims,j),     &
                U10(ims,j),V10(ims,j),TH2(ims,j),T2(ims,j),        &
                Q2(ims,j),FLHC(ims,j),FLQC(ims,j),SNOWH(ims,j),    &
                QGH(ims,j),QSFC(ims,j),LH(ims,j),                  &
                GZ1OZ0(ims,j),WSPD(ims,j),BR(ims,j),ISFFLX,DX,     &
                SVP1,SVP2,SVP3,SVPT0,EP1,EP2,KARMAN,               &
                ch(ims,j),vt1,vq1,qc1d,qcg(ims,j),                 &
                itimestep,                                         &
!JOE-begin additional output
                z0zt(ims,j), zaz0(ims,j),za_sf(ims,j),wstar(ims,j), &
                qstar(ims,j),resist(ims,j),logres(ims,j),          &
!JOE-end
                spp_pbl,rstoch1D,                                  &
                ids,ide, jds,jde, kds,kde,                         &
                ims,ime, jms,jme, kms,kme,                         &
                its,ite, jts,jte, kts,kte                          &
                ,isftcflx,iz0tlnd,                                 &
                USTM(ims,j),CK(ims,j),CKA(ims,j),                  &
                CD(ims,j),CDA(ims,j),psim_sf(ims,j),psih_sf(ims,j), &
      dth_sf(ims,j),dqv_sf(ims,j),minpsi_m,minpsi_h,li_opt,cb05_opt &
                                                                   )
!        psim_sf(ims,j) = PSIM(ims,j)
!        psih_sf(ims,j) = PSIH(ims,j)
      ENDDO

    END SUBROUTINE SFCLAY_MYNN

!-------------------------------------------------------------------
   SUBROUTINE SFCLAY1D_mynn(                                       &
                     J,U1D,V1D,T1D,QV1D,P1D,dz8w1d,rho1d,          &
                     U1D2,V1D2,dz2w1d,                             &
                     CP,G,ROVCP,R,XLV,PSFCPA,CHS,CHS2,CQS2,CPM,    &
                     PBLH,RMOL,ZNT,UST,MAVAIL,ZOL,MOL,REGIME,      &
                     PSIM,PSIH,XLAND,HFX,QFX,TSK,                  &
                     U10,V10,TH2,T2,Q2,FLHC,FLQC,SNOWH,QGH,        &
                     QSFC,LH,GZ1OZ0,WSPD,BR,ISFFLX,DX,             &
                     SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                     KARMAN,ch,vt1,vq1,qc1d,qcg,                   &
                     itimestep,                                    &
!JOE-additional output
                     zratio, zzm1d, za1d, wstar,qstar,resist,logres, &
!JOE-end
                     spp_pbl,rstoch1D,                             &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte                     &
                     ,isftcflx, iz0tlnd,                           &
                     ustm,ck,cka,cd,cda,psim1d,psih1d,   &       
             dth_sf1d,dqv_sf1d,minpsi_m,minpsi_h,li_opt,cb05_opt   &
                     )

!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
! SCALARS
!-----------------------------
      INTEGER,  INTENT(IN) ::        ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte, &
                                     J, itimestep

      REAL,     PARAMETER  :: XKA=2.4E-5   !molecular diffusivity
      REAL,     PARAMETER  :: PRT=1.       !prandlt number
      REAL,     INTENT(IN) :: SVP1,SVP2,SVP3,SVPT0,EP1,EP2
      REAL,     INTENT(IN) :: KARMAN,CP,G,ROVCP,R,XLV,DX

!-----------------------------
! NAMELIST OPTIONS
!-----------------------------
      INTEGER,  INTENT(IN) :: ISFFLX
      REAL, INTENT(IN) :: minpsi_m,minpsi_h
      INTEGER,  INTENT(IN) :: li_opt,cb05_opt
      INTEGER,  OPTIONAL,  INTENT(IN )   ::     ISFTCFLX, IZ0TLND
      INTEGER,    INTENT(IN)             ::     spp_pbl

!-----------------------------
! 1D ARRAYS
!-----------------------------
      REAL,     DIMENSION( ims:ime ), INTENT(IN)    ::     MAVAIL, &
                                                             PBLH, &
                                                            XLAND, &
                                                              TSK, &
                                                           PSFCPA, &
                                                              QCG, &
                                                            SNOWH

      REAL,     DIMENSION( its:ite ), INTENT(IN)   ::     U1D,V1D, &
                                                        U1D2,V1D2, &
                                                         QV1D,P1D, &
                                                         T1D,QC1d, &
                                                    dz8w1d,dz2w1d, &
                                                            RHO1D, &
                                                          vt1,vq1

      REAL,     DIMENSION( ims:ime ), INTENT(INOUT) ::     REGIME, &
                                                       HFX,QFX,LH, &
                                                         MOL,RMOL, &
                                                         QGH,QSFC, &
                                                              ZNT, &
                                                              ZOL, &
                                                              UST, &
                                                              CPM, &
                                                        CHS2,CQS2, &
                                                           CHS,CH, &
                                                        FLHC,FLQC, &
                                                           GZ1OZ0, &
                                                             WSPD, &
                                                               BR, &
                                                        PSIM,PSIH
      REAL,     DIMENSION( its:ite ), INTENT(IN)   ::     rstoch1D

      ! DIAGNOSTIC OUTPUT
      REAL,     DIMENSION( ims:ime ), INTENT(OUT)   ::    U10,V10, &
                                                        TH2,T2,Q2

      REAL, OPTIONAL, DIMENSION( ims:ime )                       , &
                INTENT(OUT)     ::              ck,cka,cd,cda,ustm
!--------------------------------------------
!JOE-additinal output
      REAL,     DIMENSION( ims:ime ) ::        zratio,wstar,qstar, &
                                         resist,logres, zzm1d, za1d
      REAL,     DIMENSION( ims:ime ) ::  psim1d, psih1d, dth_sf1d, dqv_sf1d
!JOE-end
!----------------------------------------------------------------
! LOCAL VARS
!----------------------------------------------------------------
      REAL :: thl1,sqv1,sqc1,exner1,sqvg,sqcg,vv,ww

      REAL, DIMENSION(its:ite) :: &
                 ZA, &    !Height of lowest 1/2 sigma level(m)
                ZA2, &    !Height of 2nd lowest 1/2 sigma level(m)
              THV1D, &    !Theta-v at lowest 1/2 sigma (K)
               TH1D, &    !Theta at lowest 1/2 sigma (K)
               TC1D, &    !T at lowest 1/2 sigma (Celsius)
               TV1D, &    !Tv at lowest 1/2 sigma (K)
               QVSH, &    !qv at lowest 1/2 sigma (spec humidity)
        PSIH2,PSIM2, &    !M-O stability functions at z=2 m
      PSIH10,PSIM10, &    !M-O stability functions at z=10 m
              WSPDI, & 
            z_t,z_q, &    !thermal & moisture roughness lengths
           ZNTstoch, &
             GOVRTH, &    !g/theta
               THGB, &    !theta at ground
              THVGB, &    !theta-v at ground
               PSFC, &    !press at surface (Pa/1000)
             QSFCMR, &    !qv at surface (mixing ratio, kg/kg)
             GZ2OZ0, &    !LOG((2.0+ZNT(I))/ZNT(I))
            GZ10OZ0, &    !LOG((10.+ZNT(I))/ZNT(I))
             GZ2OZt, &    !LOG((2.0+z_t(i))/z_t(i))
            GZ10OZt, &    !LOG((10.+z_t(i))/z_t(i))
             GZ1OZt       !LOG((ZA(I)+z_t(i))/z_t(i))

      INTEGER ::  N,I,K,L,NZOL,NK,NZOL2,NZOL10, ITER, yesno
      INTEGER, PARAMETER :: ITMAX=1
      REAL    ::  PL,THCON,TVCON,E1
      REAL    ::  DTHVDZ,DTHVM,VCONV,RZOL,RZOL2,RZOL10,ZOL2,ZOL10
      REAL    ::  DTG,PSIX,DTTHX,DTHDZ,PSIX10,PSIT,PSIT2,PSIT10, &
                  PSIQ,PSIQ2,PSIQ10
      REAL    ::  FLUXC,VSGD
      REAL    ::  restar,VISC,DQG,OLDUST,OLDTST
!      REAL, PARAMETER :: psilim = minpsi_sf
      REAL, PARAMETER :: psilim = -10.  ! ONLY AFFECTS z/L > 2.0
!      psilim = minpsi_sf
!-------------------------------------------------------------------

      DO I=its,ite
         ! CONVERT GROUND & LOWEST LAYER TEMPERATURE TO POTENTIAL TEMPERATURE:
         ! PSFC cmb
         PSFC(I)=PSFCPA(I)/1000.
         THGB(I)=TSK(I)*(100./PSFC(I))**ROVCP   !(K)              
         ! PL cmb
         PL=P1D(I)/1000.                                                   
         THCON=(100./PL)**ROVCP                                                 
         TH1D(I)=T1D(I)*THCON                   !(Theta, K)
         TC1D(I)=T1D(I)-273.15                  !(T, Celsius)    

         ! CONVERT TO VIRTUAL TEMPERATURE
         QVSH(I)=QV1D(I)/(1.+QV1D(I))        !CONVERT TO SPEC HUM (kg/kg)
         TVCON=(1.+EP1*QVSH(I))
         THV1D(I)=TH1D(I)*TVCON                 !(K)
         TV1D(I)=T1D(I)*TVCON   !(K)

         !RHO1D(I)=PSFCPA(I)/(R*TV1D(I)) !now using value calculated in sfc driver
         ZA(I)=0.5*dz8w1d(I)             !height of first half-sigma level 
         ZA2(I)=dz8w1d(I) + 0.5*dz2w1d(I)    !height of 2nd half-sigma level
         GOVRTH(I)=G/TH1D(I)
      ENDDO

      DO I=its,ite
         IF (TSK(I) .LT. 273.15) THEN
            !SATURATION VAPOR PRESSURE WRT ICE (SVP1=.6112; 10*mb)
            E1=SVP1*EXP(4648*(1./273.15 - 1./TSK(I)) - &
            & 11.64*LOG(273.15/TSK(I)) + 0.02265*(273.15 - TSK(I)))
         ELSE
            !SATURATION VAPOR PRESSURE WRT WATER (Bolton 1980)
            E1=SVP1*EXP(SVP2*(TSK(I)-SVPT0)/(TSK(I)-SVP3))
         ENDIF
         !FOR LAND POINTS, QSFC can come from LSM, ONLY RECOMPUTE OVER WATER
         IF (xland(i).gt.1.5 .or. QSFC(i).le.0.0) THEN   !WATER
            QSFC(I)=EP2*E1/(PSFC(I)-ep_3*E1)             !specific humidity    
            QSFCMR(I)=EP2*E1/(PSFC(I)-E1)                !mixing ratio 
         ELSE                                            !LAND 
            QSFCMR(I)=QSFC(I)/(1.-QSFC(I))
         ENDIF

         ! QGH CHANGED TO USE LOWEST-LEVEL AIR TEMP CONSISTENT WITH MYJSFC CHANGE
         ! Q2SAT = QGH IN LSM
         IF (TSK(I) .LT. 273.15) THEN
            !SATURATION VAPOR PRESSURE WRT ICE
            E1=SVP1*EXP(4648*(1./273.15 - 1./T1D(I)) - &
            &  11.64*LOG(273.15/T1D(I)) + 0.02265*(273.15 - T1D(I)))
         ELSE
            !SATURATION VAPOR PRESSURE WRT WATER (Bolton 1980)
            E1=SVP1*EXP(SVP2*(T1D(I)-SVPT0)/(T1D(I)-SVP3))
         ENDIF
         PL=P1D(I)/1000.
         !QGH(I)=EP2*E1/(PL-ep_3*E1)    !specific humidity
         QGH(I)=EP2*E1/(PL-E1)          !mixing ratio
         CPM(I)=CP*(1.+0.84*QV1D(I))
      ENDDO

      DO I=its,ite
         WSPD(I)=SQRT(U1D(I)*U1D(I)+V1D(I)*V1D(I))     

         !account for partial condensation
         exner1=(p1d(I)/p1000mb)**ROVCP
         sqc1=qc1d(I)/(1.+qc1d(I))         !lowest mod level cloud water spec hum
         sqv1=QVSH(I)                      !lowest mod level water vapor spec hum
         thl1=TH1D(I)-xlvcp/exner1*sqc1
         sqvg=qsfc(I)                      !sfc water vapor spec hum
         sqcg=qcg(I)/(1.+qcg(I))           !sfc cloud water spec hum

         vv = thl1-THGB(I)
         !TGS:ww = mavail(I)*(sqv1-sqvg) + (sqc1-sqcg)
         ww = (sqv1-sqvg) + (sqc1-sqcg)

         !TGS:THVGB(I)=THGB(I)*(1.+EP1*QSFC(I)*MAVAIL(I)) 
         THVGB(I)=THGB(I)*(1.+EP1*QSFC(I))

         DTHDZ=(TH1D(I)-THGB(I))
         DTHVDZ=(THV1D(I)-THVGB(I))
         !DTHVDZ= (vt1(i) + 1.0)*vv + (vq1(i) + tv0)*ww

         !--------------------------------------------------------
         !  Calculate the convective velocity scale (WSTAR) and 
         !  subgrid-scale velocity (VSGD) following Beljaars (1995, QJRMS) 
         !  and Mahrt and Sun (1995, MWR), respectively
         !-------------------------------------------------------
         !  Use Beljaars over land and water
         fluxc = max(hfx(i)/RHO1D(i)/cp                    &
         &    + ep1*THVGB(I)*qfx(i)/RHO1D(i),0.)
         WSTAR(I) = vconvc*(g/TSK(i)*pblh(i)*fluxc)**.33

         !--------------------------------------------------------
         ! Mahrt and Sun low-res correction
         ! (for 13 km ~ 0.37 m/s; for 3 km == 0 m/s)
         !--------------------------------------------------------
         VSGD = 0.32 * (max(dx/5000.-1.,0.))**.33
         WSPD(I)=SQRT(WSPD(I)*WSPD(I)+WSTAR(I)*WSTAR(I)+vsgd*vsgd)
         WSPD(I)=MAX(WSPD(I),wmin)

         !--------------------------------------------------------
         ! CALCULATE THE BULK RICHARDSON NUMBER OF SURFACE LAYER, 
         ! ACCORDING TO AKB(1976), EQ(12). 
         !--------------------------------------------------------
         BR(I)=GOVRTH(I)*ZA(I)*DTHVDZ/(WSPD(I)*WSPD(I))
         !SET LIMITS ACCORDING TO Li et al. (2010) Boundary-Layer Meteorol (p.158)
         BR(I)=MAX(BR(I),-20.0)
         BR(I)=MIN(BR(I),2.0)
                
         ! IF PREVIOUSLY UNSTABLE, DO NOT LET INTO REGIMES 1 AND 2 (STABLE)
         !if (itimestep .GT. 1) THEN
         !    IF(MOL(I).LT.0.)BR(I)=MIN(BR(I),0.0)
         !ENDIF
     
         !IF(I .eq. 2)THEN
         !  write(*,1006)"BR:",BR(I)," fluxc:",fluxc," vt1:",vt1(i)," vq1:",vq1(i)
         !  write(*,1007)"XLAND:",XLAND(I)," WSPD:",WSPD(I)," DTHVDZ:",DTHVDZ," WSTAR:",WSTAR(I)
         !ENDIF

      ENDDO

 1006   format(A,F7.3,A,f9.4,A,f9.5,A,f9.4)
 1007   format(A,F2.0,A,f6.2,A,f7.3,A,f7.2)

!--------------------------------------------------------------------      
!--------------------------------------------------------------------      
!--- BEGIN ITERATION LOOP (ITMAX=5); USUALLY CONVERGES IN TWO PASSES
!--------------------------------------------------------------------
!--------------------------------------------------------------------

 DO I=its,ite

   ITER = 1
   DO WHILE (ITER .LE. ITMAX)

      !COMPUTE KINEMATIC VISCOSITY (m2/s) Andreas (1989) CRREL Rep. 89-11
      !valid between -173 and 277 degrees C.
      VISC=1.326e-5*(1. + 6.542e-3*TC1D(I) + 8.301e-6*TC1D(I)*TC1D(I) &
                        - 4.84e-9*TC1D(I)*TC1D(I)*TC1D(I))

      IF((XLAND(I)-1.5).GE.0)THEN
          !--------------------------------------
          ! WATER
          !--------------------------------------
          ! CALCULATE z0 (znt)
          !--------------------------------------
          IF ( PRESENT(ISFTCFLX) ) THEN
             IF ( ISFTCFLX .EQ. 0 ) THEN
                IF (COARE_OPT .EQ. 3.0) THEN
                  !COARE 3.0 (MISLEADING SUBROUTINE NAME)
                  CALL charnock_1955(ZNT(i),UST(i),WSPD(i),visc,ZA(I))
                ELSE
                  !COARE 3.5
                  CALL edson_etal_2013(ZNT(i),UST(i),WSPD(i),visc,ZA(I))
                ENDIF
             ELSEIF ( ISFTCFLX .EQ. 1 .OR. ISFTCFLX .EQ. 2 ) THEN
                CALL davis_etal_2008(ZNT(i),UST(i))
             ELSEIF ( ISFTCFLX .EQ. 3 ) THEN
                CALL Taylor_Yelland_2001(ZNT(i),UST(i),WSPD(i))
             ELSEIF ( ISFTCFLX .EQ. 4 ) THEN
                IF (COARE_OPT .EQ. 3.0) THEN
                  !COARE 3.0 (MISLEADING SUBROUTINE NAME)
                  CALL charnock_1955(ZNT(i),UST(i),WSPD(i),visc,ZA(I))
                ELSE
                  !COARE 3.5
                  CALL edson_etal_2013(ZNT(i),UST(i),WSPD(i),visc,ZA(I))
                ENDIF
             ENDIF
          ELSE
             !DEFAULT TO COARE 3.0/3.5
             IF (COARE_OPT .EQ. 3.0) THEN
                !COARE 3.0
                CALL charnock_1955(ZNT(i),UST(i),WSPD(i),visc,ZA(I))
             ELSE
                !COARE 3.5
                CALL edson_etal_2013(ZNT(i),UST(i),WSPD(i),visc,ZA(I))
             ENDIF
          ENDIF

          ! add stochastic perturbaction of ZNT
          if (spp_pbl==1) then
             ZNTstoch(I)  = MAX(ZNT(I) + 1.5 * ZNT(I) * rstoch1D(i), 1e-6)
          else
             ZNTstoch(I)  = ZNT(I)
          endif

          !COMPUTE ROUGHNESS REYNOLDS NUMBER (restar) USING NEW ZNT
          ! AHW: Garrattt formula: Calculate roughness Reynolds number
          !      Kinematic viscosity of air (linear approx to
          !      temp dependence at sea level)
          restar=MAX(ust(i)*ZNTstoch(i)/visc, 0.1)

          !--------------------------------------
          !CALCULATE z_t and z_q
          !--------------------------------------
          IF ( PRESENT(ISFTCFLX) ) THEN
             IF ( ISFTCFLX .EQ. 0 ) THEN
                IF (COARE_OPT .EQ. 3.0) THEN
                   CALL fairall_etal_2003(z_t(i),z_q(i),restar,UST(i),visc)
                ELSE
                   !presumably, this will be published soon, but hasn't yet
                   CALL fairall_etal_2014(z_t(i),z_q(i),restar,UST(i),visc,rstoch1D(i),spp_pbl)
                ENDIF
             ELSEIF ( ISFTCFLX .EQ. 1 ) THEN
                IF (COARE_OPT .EQ. 3.0) THEN
                   CALL fairall_etal_2003(z_t(i),z_q(i),restar,UST(i),visc)
                ELSE
                   CALL fairall_etal_2014(z_t(i),z_q(i),restar,UST(i),visc,rstoch1D(i),spp_pbl)
                ENDIF
             ELSEIF ( ISFTCFLX .EQ. 2 ) THEN
                CALL garratt_1992(z_t(i),z_q(i),ZNTstoch(i),restar,XLAND(I))
             ELSEIF ( ISFTCFLX .EQ. 3 ) THEN
                IF (COARE_OPT .EQ. 3.0) THEN
                   CALL fairall_etal_2003(z_t(i),z_q(i),restar,UST(i),visc)
                ELSE
                   CALL fairall_etal_2014(z_t(i),z_q(i),restar,UST(i),visc,rstoch1D(i),spp_pbl)
                ENDIF
             ELSEIF ( ISFTCFLX .EQ. 4 ) THEN
                CALL zilitinkevich_1995(ZNTstoch(i),z_t(i),z_q(i),restar,&
                                   UST(I),KARMAN,XLAND(I),IZ0TLND,spp_pbl,rstoch1D(i))
             ENDIF
          ELSE
             !DEFAULT TO COARE 3.0/3.5
             IF (COARE_OPT .EQ. 3.0) THEN
                CALL fairall_etal_2003(z_t(i),z_q(i),restar,UST(i),visc)
             ELSE
                CALL fairall_etal_2014(z_t(i),z_q(i),restar,UST(i),visc,rstoch1D(i),spp_pbl)
             ENDIF
          ENDIF
 
       ELSE

          ! add stochastic perturbaction of ZNT
          if (spp_pbl==1) then
             ZNTstoch(I)  = MAX(ZNT(I) + 1.5 * ZNT(I) * rstoch1D(i), 1e-6)
          else
             ZNTstoch(I)  = ZNT(I)
          endif

          !--------------------------------------
          ! LAND
          !--------------------------------------
          !COMPUTE ROUGHNESS REYNOLDS NUMBER (restar) USING DEFAULT ZNT
          restar=MAX(ust(i)*ZNTstoch(i)/visc, 0.1)

          !--------------------------------------
          !GET z_t and z_q
          !--------------------------------------
          !CHECK FOR SNOW/ICE POINTS OVER LAND
          !IF ( ZNTSTOCH(i) .LE. SNOWZ0  .AND.  TSK(I) .LE. 273.15 ) THEN
          IF ( SNOWH(i) .GE. 0.1) THEN
             CALL Andreas_2002(ZNTSTOCH(i),visc,ust(i),z_t(i),z_q(i))
          ELSE
             IF ( PRESENT(IZ0TLND) ) THEN
                IF ( IZ0TLND .LE. 1 .OR. IZ0TLND .EQ. 4) THEN
                   !IF IZ0TLND==4, THEN PSIQ WILL BE RECALCULATED USING
                   !PAN ET AL (1994), but PSIT FROM ZILI WILL BE USED.
                   CALL zilitinkevich_1995(ZNTSTOCH(i),z_t(i),z_q(i),restar,&
                                  UST(I),KARMAN,XLAND(I),IZ0TLND,spp_pbl,rstoch1D(i))
                ELSEIF ( IZ0TLND .EQ. 2 ) THEN
                   CALL Yang_2008(ZNTSTOCH(i),z_t(i),z_q(i),UST(i),MOL(I),&
                                  qstar(I),restar,visc,XLAND(I))
                ELSEIF ( IZ0TLND .EQ. 3 ) THEN
                   !Original MYNN in WRF-ARW used this form:
                   CALL garratt_1992(z_t(i),z_q(i),ZNTSTOCH(i),restar,XLAND(I))
                ENDIF
             ELSE
                !DEFAULT TO ZILITINKEVICH
                CALL zilitinkevich_1995(ZNTSTOCH(i),z_t(i),z_q(i),restar,&
                                        UST(I),KARMAN,XLAND(I),0,spp_pbl,rstoch1D(i))
             ENDIF
          ENDIF

       ENDIF
       zratio(i)=zntstoch(i)/z_t(i)
       zzm1d(i)=ZA(I)/ZNTstoch(I)
       za1d(i) = ZA(I)
       !ADD RESISTANCE (SOMEWHAT FOLLOWING JIMENEZ ET AL. (2012)) TO PROTECT AGAINST
       !EXCESSIVE FLUXES WHEN USING A LOW FIRST MODEL LEVEL (ZA < 10 m).        
       !Formerly: GZ1OZ0(I)= LOG(ZA(I)/ZNTstoch(I))
       GZ1OZ0(I)= LOG((ZA(I)+ZNTstoch(I))/ZNTstoch(I))
       GZ1OZt(I)= LOG((ZA(I)+z_t(i))/z_t(i))           
       GZ2OZ0(I)= LOG((2.0+ZNTstoch(I))/ZNTstoch(I))                                        
       GZ2OZt(I)= LOG((2.0+z_t(i))/z_t(i))                                        
       GZ10OZ0(I)=LOG((10.+ZNTstoch(I))/ZNTstoch(I)) 
       GZ10OZt(I)=LOG((10.+z_t(i))/z_t(i)) 

     !--------------------------------------------------------------------      
     !--- DIAGNOSE BASIC PARAMETERS FOR THE APPROPRIATE STABILITY CLASS:
     !                                                                                
     !    THE STABILITY CLASSES ARE DETERMINED BY BR (BULK RICHARDSON NO.).
     !                                                                                
     !    CRITERIA FOR THE CLASSES ARE AS FOLLOWS:                                   
     !                                                                                
     !        1. BR .GE. 0.2;                                                         
     !               REPRESENTS NIGHTTIME STABLE CONDITIONS (REGIME=1),               
     !                                                                                
     !        2. BR .LT. 0.2 .AND. BR .GT. 0.0;                                       
     !               REPRESENTS DAMPED MECHANICAL TURBULENT CONDITIONS                
     !               (REGIME=2),                                                      
     !                                                                                
     !        3. BR .EQ. 0.0                                                          
     !               REPRESENTS FORCED CONVECTION CONDITIONS (REGIME=3),              
     !                                                                                
     !        4. BR .LT. 0.0                                                          
     !               REPRESENTS FREE CONVECTION CONDITIONS (REGIME=4).                
     !                                                                                
     !--------------------------------------------------------------------
     IF (BR(I) .GT. 0.0) THEN
        IF (BR(I) .GT. 0.2) THEN        
            !---CLASS 1; STABLE (NIGHTTIME) CONDITIONS:                                    
            REGIME(I)=1.
        ELSE
            !---CLASS 2; DAMPED MECHANICAL TURBULENCE:
            REGIME(I)=2.
        ENDIF

        !COMPUTE z/L
        !CALL Li_etal_2010(ZOL(I),BR(I),ZA(I)/ZNTstoch(I),zratio(I))
!        IF (ITER .EQ. 1 .AND. itimestep .LE. 1) THEN
         if (li_opt.eq.0) then
            CALL Li_etal_2010(ZOL(I),BR(I),ZA(I)/ZNTstoch(I),zratio(I))
         end if
         if (li_opt.eq.1) then
            CALL Li_etal_2014(ZOL(I),BR(I),ZA(I),ZNTstoch(I),z_t(I))
         end if
!        ELSE
!           ZOL(I)=ZA(I)*KARMAN*G*MOL(I)/(TH1D(I)*MAX(UST(I)*UST(I),0.0001))
!           ZOL(I)=MAX(ZOL(I),0.0)
!           ZOL(I)=MIN(ZOL(I),2.)
!        ENDIF
 
        !COMPUTE PSIM and PSIH
        IF((XLAND(I)-1.5).GE.0)THEN                                            
           ! WATER
           !CALL PSI_Suselj_Sood_2010(PSIM(I),PSIH(I),ZOL(I))
           !CALL PSI_Beljaars_Holtslag_1991(PSIM(I),PSIH(I),ZOL(I))
           !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
           if (ZOL(I).le.1.0) then
              CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
           end if
           if (ZOL(I).gt.1.0) then
              if (cb05_opt.eq.0) then
                 CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
              end if
              if (cb05_opt.eq.1) then
                 CALL PSI_CB2005(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
              end if
              if (cb05_opt.eq.2) then
                 CALL PSI_G2007(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
              end if
           end if
        ELSE
           ! LAND  
           !CALL PSI_Beljaars_Holtslag_1991(PSIM(I),PSIH(I),ZOL(I))
           !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
           !CALL PSI_Zilitinkevich_Esau_2007(PSIM(I),PSIH(I),ZOL(I))
           if (ZOL(I).le.1.0) then
              CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
           end if
           if (ZOL(I).gt.1.0) then
              if (cb05_opt.eq.0) then
                 CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
              end if
              if (cb05_opt.eq.1) then
                 CALL PSI_CB2005(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
              end if
              if (cb05_opt.eq.2) then
                 CALL PSI_G2007(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
              end if
           end if
        ENDIF              

        ! LOWER LIMIT ON PSI IN STABLE CONDITIONS
!        PSIM10(I)=MAX(10./ZA(I)*PSIM(I), psilim)
!        PSIH10(I)=MAX(10./ZA(I)*PSIH(I), psilim)
!        PSIM2(I)=MAX(2./ZA(I)*PSIM(I), psilim)
!        PSIH2(I)=MAX(2./ZA(I)*PSIH(I), psilim)
!        PSIM(I)=MAX(PSIM(I),psilim)
!        PSIH(I)=MAX(PSIH(I),psilim)
                                     
        PSIM10(I)=MAX(10./ZA(I)*PSIM(I), minpsi_m)
        PSIH10(I)=MAX(10./ZA(I)*PSIH(I), minpsi_h)
        PSIM2(I)=MAX(2./ZA(I)*PSIM(I), minpsi_m)
        PSIH2(I)=MAX(2./ZA(I)*PSIH(I), minpsi_h)
        PSIM(I)=MAX(PSIM(I),minpsi_m)
        PSIH(I)=MAX(PSIH(I),minpsi_h)
        ! 1.0 over Monin-Obukhov length
        RMOL(I)= ZOL(I)/ZA(I)

     ELSEIF(BR(I) .EQ. 0.) THEN                  
        !=========================================================  
        !-----CLASS 3; FORCED CONVECTION/NEUTRAL:                                                
        !=========================================================
        REGIME(I)=3.

        PSIM(I)=0.0
        PSIH(I)=PSIM(I)
        PSIM10(I)=0.
        PSIH10(I)=PSIM10(I)
        PSIM2(I)=0.
        PSIH2(I)=PSIM2(I)

        !ZOL(I)=0.
        IF(UST(I) .LT. 0.01)THEN
          ZOL(I)=BR(I)*GZ1OZ0(I)
        ELSE
          ZOL(I)=KARMAN*GOVRTH(I)*ZA(I)*MOL(I)/(MAX(UST(I)*UST(I),0.001))
        ENDIF
        RMOL(I) = ZOL(I)/ZA(I)

     ELSEIF(BR(I) .LT. 0.)THEN
        !==========================================================
        !-----CLASS 4; FREE CONVECTION:                                                  
        !==========================================================
        REGIME(I)=4.

        !COMPUTE z/L
        !CALL Li_etal_2010(ZOL(I),BR(I),ZA(I)/ZNTstoch(I),zratio(I))
        !IF (ITER .EQ. 1 .AND. itimestep .LE. 1) THEN
        if (li_opt.eq.0) then
           CALL Li_etal_2010(ZOL(I),BR(I),ZA(I)/ZNTstoch(I),zratio(I))
        end if
        if (li_opt.eq.1) then
           CALL Li_etal_2014(ZOL(I),BR(I),ZA(I),ZNTstoch(I),z_t(I))
        end if
        !ELSE
        !   ZOL(I)=ZA(I)*KARMAN*G*MOL(I)/(TH1D(I)*MAX(UST(I)*UST(I),0.001))
        !   ZOL(I)=MAX(ZOL(I),-19.999)
        !   ZOL(I)=MIN(ZOL(I),0.0)
        !ENDIF

        ZOL10=10./ZA(I)*ZOL(I)
        ZOL2=2./ZA(I)*ZOL(I)
        ZOL(I)=MIN(ZOL(I),0.)
        ZOL(I)=MAX(ZOL(I),-19.9999)
        ZOL10=MIN(ZOL10,0.)
        ZOL10=MAX(ZOL10,-19.9999)
        ZOL2=MIN(ZOL2,0.)
        ZOL2=MAX(ZOL2,-19.9999)
        NZOL=INT(-ZOL(I)*100.)
        RZOL=-ZOL(I)*100.-NZOL
        NZOL10=INT(-ZOL10*100.)
        RZOL10=-ZOL10*100.-NZOL10
        NZOL2=INT(-ZOL2*100.)
        RZOL2=-ZOL2*100.-NZOL2

        !COMPUTE PSIM and PSIH
        IF((XLAND(I)-1.5).GE.0)THEN                                            
           ! WATER
           !CALL PSI_Suselj_Sood_2010(PSIM(I),PSIH(I),ZOL(I))
           !CALL PSI_Hogstrom_1996(PSIM(I),PSIH(I),ZOL(I), z_t(I), ZNTstoch(I), ZA(I))
           !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
           if (ZOL(I).le.1.0) then
              CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
              CALL PSI_DyerHicks(PSIM2(I),PSIH2(I),ZOL2,z_t(I),ZNTstoch(I),2.)
           end if
           if (ZOL(I).gt.1.0) then
              if (cb05_opt.eq.0) then
                 CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
                 CALL PSI_DyerHicks(PSIM2(I),PSIH2(I),ZOL2,z_t(I),ZNTstoch(I),2.)
              end if
              if (cb05_opt.eq.1) then
                 CALL PSI_CB2005(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
                 CALL PSI_CB2005(PSIM2(I),PSIH2(I),ZOL2,z_t(I),ZNTstoch(I),2.)
              end if
              if (cb05_opt.eq.2) then
                 CALL PSI_G2007(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
                 CALL PSI_G2007(PSIM2(I),PSIH2(I),ZOL2,z_t(I),ZNTstoch(I),2.)
              end if
           end if
        ELSE           
           ! LAND  
           !CALL PSI_Hogstrom_1996(PSIM(I),PSIH(I),ZOL(I), z_t(I), ZNTstoch(I), ZA(I))
           !CALL PSI_Businger_1971(PSIM(I),PSIH(I),ZOL(I))
           if (ZOL(I).le.1.0) then
              CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
              CALL PSI_DyerHicks(PSIM2(I),PSIH2(I),ZOL2,z_t(I),ZNTstoch(I),2.)
           end if
           if (ZOL(I).gt.1.0) then
              if (cb05_opt.eq.0) then
                 CALL PSI_DyerHicks(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
                 CALL PSI_DyerHicks(PSIM2(I),PSIH2(I),ZOL2,z_t(I),ZNTstoch(I),2.)
              end if
              if (cb05_opt.eq.1) then
                 CALL PSI_CB2005(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
                 CALL PSI_CB2005(PSIM2(I),PSIH2(I),ZOL2,z_t(I),ZNTstoch(I),2.)
              end if
              if (cb05_opt.eq.2) then
                 CALL PSI_G2007(PSIM(I),PSIH(I),ZOL(I),z_t(I),ZNTstoch(I),ZA(I))
                 CALL PSI_G2007(PSIM2(I),PSIH2(I),ZOL2,z_t(I),ZNTstoch(I),2.)
              end if
           end if
        ENDIF              

        PSIM10(I)=10./ZA(I)*PSIM(I)
        PSIH10(I)=10./ZA(I)*PSIH(I)

        !---LIMIT PSIH AND PSIM IN THE CASE OF THIN LAYERS AND
        !---HIGH ROUGHNESS.  THIS PREVENTS DENOMINATOR IN FLUXES
        !---FROM GETTING TOO SMALL
        !PSIH(I)=MIN(PSIH(I),0.9*GZ1OZt(I))    !JOE: less restricitive over forest/urban.
        PSIH(I)=MIN(PSIH(I),0.9*GZ1OZ0(I))
        PSIM(I)=MIN(PSIM(I),0.9*GZ1OZ0(I))
        !PSIH2(I)=MIN(PSIH2(I),0.9*GZ2OZt(I))  !JOE: less restricitive over forest/urban.
        PSIH2(I)=MIN(PSIH2(I),0.9*GZ2OZ0(I))
        PSIM2(I)=MIN(PSIM2(I),0.9*GZ2OZ0(I))
        PSIM10(I)=MIN(PSIM10(I),0.9*GZ10OZ0(I))
        PSIH10(I)=MIN(PSIH10(I),0.9*GZ10OZ0(I))

        RMOL(I) = ZOL(I)/ZA(I)  

     ENDIF

     !------------------------------------------------------------
     !-----COMPUTE THE FRICTIONAL VELOCITY:                                           
     !------------------------------------------------------------
     !     ZA(1982) EQS(2.60),(2.61).
      PSIX=GZ1OZ0(I)-PSIM(I)
      PSIX10=GZ10OZ0(I)-PSIM10(I)
      ! TO PREVENT OSCILLATIONS AVERAGE WITH OLD VALUE 
      OLDUST = UST(I)
      UST(I)=0.5*UST(I)+0.5*KARMAN*WSPD(I)/PSIX 
      !NON-AVERAGED: UST(I)=KARMAN*WSPD(I)/PSIX
     
      ! Compute u* without vconv for use in HFX calc when isftcflx > 0           
      WSPDI(I)=MAX(SQRT(U1D(I)*U1D(I)+V1D(I)*V1D(I)), wmin)
      IF ( PRESENT(USTM) ) THEN
         USTM(I)=0.5*USTM(I)+0.5*KARMAN*WSPDI(I)/PSIX
      ENDIF

      IF ((XLAND(I)-1.5).LT.0.) THEN        !LAND
         UST(I)=MAX(UST(I),0.005)  !Further relaxing this limit - no need to go lower
         !Keep ustm = ust over land.
         IF ( PRESENT(USTM) ) USTM(I)=UST(I)
      ENDIF

     !------------------------------------------------------------
     !-----COMPUTE THE THERMAL AND MOISTURE RESISTANCE (PSIQ AND PSIT):                                           
     !------------------------------------------------------------
      ! LOWER LIMIT ADDED TO PREVENT LARGE FLHC IN SOIL MODEL
      ! ACTIVATES IN UNSTABLE CONDITIONS WITH THIN LAYERS OR HIGH Z0
      GZ1OZt(I)= LOG((ZA(I)+z_t(i))/z_t(i))           
      GZ2OZt(I)= LOG((2.0+z_t(i))/z_t(i))                                        

      PSIT =MAX(GZ1OZt(I)-PSIH(I) ,1.)
      PSIT2=MAX(GZ2OZt(I)-PSIH2(I),1.)
      resist(I)=PSIT
      logres(I)=GZ1OZt(I)

      PSIQ=MAX(LOG((ZA(I)+z_q(i))/z_q(I))-PSIH(I) ,1.0)
      PSIQ2=MAX(LOG((2.0+z_q(i))/z_q(I))-PSIH2(I) ,1.0)

      IF((XLAND(I)-1.5).LT.0)THEN    !Land only
         IF ( IZ0TLND .EQ. 4 ) THEN
            CALL Pan_etal_1994(PSIQ,PSIQ2,UST(I),PSIH(I),PSIH2(I),&
                       & KARMAN,ZA(I))
         ENDIF
      ENDIF

      !----------------------------------------------------
      !COMPUTE THE TEMPERATURE SCALE (or FRICTION TEMPERATURE, T*)
      !----------------------------------------------------
      !DTG=TH1D(I)-THGB(I) !SWITCH TO THETA-V                                                   
      DTG=THV1D(I)-THVGB(I)
      OLDTST=MOL(I)
      MOL(I)=KARMAN*DTG/PSIT/PRT
      !t_star(I) = -HFX(I)/(UST(I)*CPM(I)*RHO1D(I))
      !t_star(I) = MOL(I)
      !----------------------------------------------------
      !COMPUTE THE MOISTURE SCALE (or q*)
      DQG=(QVSH(i)-qsfc(i))*1000.   !(kg/kg -> g/kg)
      qstar(I)=KARMAN*DQG/PSIQ/PRT

      !CHECK FOR CONVERGENCE
      IF (ITER .GE. 2) THEN
         !IF (ABS(OLDUST-UST(I)) .lt. 0.01) THEN
         IF (ABS(OLDTST-MOL(I)) .lt. 0.01) THEN
            ITER = ITER+ITMAX
         ENDIF

         !IF () THEN
         !  print*,"ITER:",ITER
         !  write(*,1001)"REGIME:",REGIME(I)," z/L:",ZOL(I)," U*:",UST(I)," Tstar:",MOL(I)
         !  write(*,1002)"PSIM:",PSIM(I)," PSIH:",PSIH(I)," W*:",WSTAR(I)," DTHV:",THV1D(I)-THVGB(I)
         !  write(*,1003)"CPM:",CPM(I)," RHO1D:",RHO1D(I)," L:",ZOL(I)/ZA(I)," DTH:",TH1D(I)-THGB(I)
         !  write(*,1004)"Z0/Zt:",zratio(I)," Z0:",ZNTstoch(I)," Zt:",z_t(I)," za:",za(I)
         !  write(*,1005)"Re:",restar," MAVAIL:",MAVAIL(I)," QSFC(I):",QSFC(I)," QVSH(I):",QVSH(I)
         !  print*,"VISC=",VISC," Z0:",ZNTstoch(I)," T1D(i):",T1D(i)
         !  write(*,*)"============================================="
         !ENDIF
      ENDIF

      ITER = ITER + 1

   ENDDO  ! end ITERATION-loop

 ENDDO     ! end i-loop

 1000   format(A,F6.1, A,f6.1, A,f5.1, A,f7.1)
 1001   format(A,F2.0, A,f10.4,A,f5.3, A,f11.5)
 1002   format(A,f7.2, A,f7.2, A,f7.2, A,f10.3)
 1003   format(A,f7.2, A,f7.2, A,f10.3,A,f10.3)
 1004   format(A,f11.3,A,f9.7, A,f9.7, A,f6.2, A,f10.3)
 1005   format(A,f9.2,A,f6.4,A,f7.4,A,f7.4)

      !----------------------------------------------------------
      !  COMPUTE SURFACE HEAT AND MOISTURE FLUXES
      !----------------------------------------------------------
 DO I=its,ite

    !For computing the diagnostics and fluxes (below), whether the fluxes
    !are turned off or on, we need the following:
    PSIX=GZ1OZ0(I)-PSIM(I)
    PSIX10=GZ10OZ0(I)-PSIM10(I)

    PSIT =MAX(GZ1OZt(I)-PSIH(I), 1.0)
    PSIT2=MAX(GZ2OZt(I)-PSIH2(I),1.0)
    PSIT10=MAX(GZ10OZ0(I)-PSIH10(I), 1.0)

    psih1d(I) = PSIH(I)
    psim1d(I) = PSIM(I)
  
    PSIQ=MAX(LOG((ZA(I)+z_q(i))/z_q(I))-PSIH(I) ,1.0)
    PSIQ2=MAX(LOG((2.0+z_q(i))/z_q(I))-PSIH2(I) ,1.0)
    PSIQ10=MAX(GZ10OZ0(I)-PSIH10(I),1.0)

    IF (ISFFLX .LT. 1) THEN                                                

       QFX(i)  = 0.                                                              
       HFX(i)  = 0.    
       FLHC(I) = 0.                                                             
       FLQC(I) = 0.                                                             
       LH(I)   = 0.                                                             
       CHS(I)  = 0.                                                             
       CH(I)   = 0.                                                             
       CHS2(i) = 0.                                                              
       CQS2(i) = 0.                                                              
       IF(PRESENT(ck)  .and. PRESENT(cd) .and. &
         &PRESENT(cka) .and. PRESENT(cda)) THEN
           Ck(I) = 0.
           Cd(I) = 0.
           Cka(I)= 0.
           Cda(I)= 0.
       ENDIF
    ELSE

      IF((XLAND(I)-1.5).LT.0)THEN !LAND Only  
         IF ( IZ0TLND .EQ. 4 ) THEN
            CALL Pan_etal_1994(PSIQ,PSIQ2,UST(I),PSIH(I),PSIH2(I),&
                       & KARMAN,ZA(I))
         ENDIF
      ENDIF

      !------------------------------------------
      ! CALCULATE THE EXCHANGE COEFFICIENTS FOR HEAT (FLHC)
      ! AND MOISTURE (FLQC)
      !------------------------------------------
      FLQC(I)=RHO1D(I)*MAVAIL(I)*UST(I)*KARMAN/PSIQ
      FLHC(I)=RHO1D(I)*CPM(I)*UST(I)*KARMAN/PSIT
      !OLD WAY:
      !DTTHX=ABS(TH1D(I)-THGB(I))                                            
      !IF(DTTHX.GT.1.E-5)THEN                                                   
      !   FLHC(I)=CPM(I)*RHO1D(I)*UST(I)*MOL(I)/(TH1D(I)-THGB(I))          
      !ELSE                                                                     
      !   FLHC(I)=0.                                                             
      !ENDIF   

      !----------------------------------
      ! COMPUTE SURFACE MOISTURE FLUX:
      !----------------------------------
      dqv_sf1d(I) = QSFCMR(I)-QV1D(I)
      QFX(I)=FLQC(I)*(QSFCMR(I)-QV1D(I))          
      !JOE: QFX(I)=MAX(QFX(I),0.)   !originally did not allow neg QFX 
      IF (QFX(I).lt.-0.02) then          
         QFX(I)=MAX(QFX(I),-0.02)   !allows small neg QFX, like MYJ 
         dqv_sf1d(I) = -0.02/FLQC(I)
      END IF
      LH(I)=XLV*QFX(I)

      !----------------------------------
      ! COMPUTE SURFACE HEAT FLUX:
      !----------------------------------
      dth_sf1d(I)=THGB(I)-TH1D(I)
      IF(XLAND(I)-1.5.GT.0.)THEN      !WATER                                           
         HFX(I)=FLHC(I)*(THGB(I)-TH1D(I))                                
         IF ( PRESENT(ISFTCFLX) ) THEN
            IF ( ISFTCFLX.NE.0 ) THEN
               ! AHW: add dissipative heating term
               HFX(I)=HFX(I)+RHO1D(I)*USTM(I)*USTM(I)*WSPDI(I)
            ENDIF
         ENDIF
      ELSEIF(XLAND(I)-1.5.LT.0.)THEN  !LAND                               
         HFX(I)=FLHC(I)*(THGB(I)-TH1D(I))
         IF (HFX(I).lt.-250.) THEN                                
            HFX(I)=MAX(HFX(I),-250.) 
            dth_sf1d(I)=-250.0/FLHC(I)
         END IF                                     
      ENDIF

      !CHS(I)=UST(I)*KARMAN/(ALOG(KARMAN*UST(I)*ZA(I) &
      !       /XKA+ZA(I)/ZL)-PSIH(I))

      CHS(I)=UST(I)*KARMAN/PSIT

      ! The exchange coefficient for cloud water is assumed to be the 
      ! same as that for heat. CH is multiplied by WSPD.

      !ch(i)=chs(i)
      ch(i)=flhc(i)/( cpm(i)*RHO1D(i) )

      !THESE ARE USED FOR 2-M DIAGNOSTICS ONLY
      CQS2(I)=UST(I)*KARMAN/PSIQ2
      CHS2(I)=UST(I)*KARMAN/PSIT2

      IF(PRESENT(ck)  .and. PRESENT(cd) .and. &
        &PRESENT(cka) .and. PRESENT(cda)) THEN
         Ck(I)=(karman/psix10)*(karman/psiq10)
         Cd(I)=(karman/psix10)*(karman/psix10)
         Cka(I)=(karman/psix)*(karman/psiq)
         Cda(I)=(karman/psix)*(karman/psix)
      ENDIF

   ENDIF !end ISFFLX option

   !-----------------------------------------------------
   !COMPUTE DIAGNOSTICS
   !-----------------------------------------------------
   !COMPUTE 10 M WNDS
   !-----------------------------------------------------
   ! If the lowest model level is close to 10-m, use it
   ! instead of the flux-based diagnostic formula.
   if (ZA(i) .le. 7.0) then
      ! high vertical resolution
      if(ZA2(i) .gt. 7.0 .and. ZA2(i) .lt. 13.0) then
         !use 2nd model level
         U10(I)=U1D2(I)
         V10(I)=V1D2(I)
      else
         U10(I)=U1D(I)*log(10./ZNTstoch(I))/log(ZA(I)/ZNTstoch(I))
         V10(I)=V1D(I)*log(10./ZNTstoch(I))/log(ZA(I)/ZNTstoch(I))
      endif
   elseif(ZA(i) .gt. 7.0 .and. ZA(i) .lt. 13.0) then
      !moderate vertical resolution
      !U10(I)=U1D(I)*PSIX10/PSIX
      !V10(I)=V1D(I)*PSIX10/PSIX
      !use neutral-log:
      U10(I)=U1D(I)*log(10./ZNTstoch(I))/log(ZA(I)/ZNTstoch(I))
      V10(I)=V1D(I)*log(10./ZNTstoch(I))/log(ZA(I)/ZNTstoch(I))
   else
      ! very coarse vertical resolution
      U10(I)=U1D(I)*PSIX10/PSIX
      V10(I)=V1D(I)*PSIX10/PSIX
   endif

   !-----------------------------------------------------
   !COMPUTE 2m T, TH, AND Q
   !THESE WILL BE OVERWRITTEN FOR LAND POINTS IN THE LSM
   !-----------------------------------------------------
   DTG=TH1D(I)-THGB(I) 
   TH2(I)=THGB(I)+DTG*PSIT2/PSIT
   !***  BE CERTAIN THAT THE 2-M THETA IS BRACKETED BY
   !***  THE VALUES AT THE SURFACE AND LOWEST MODEL LEVEL.
   IF ((TH1D(I)>THGB(I) .AND. (TH2(I)<THGB(I) .OR. TH2(I)>TH1D(I))) .OR. &
       (TH1D(I)<THGB(I) .AND. (TH2(I)>THGB(I) .OR. TH2(I)<TH1D(I)))) THEN
       TH2(I)=THGB(I) + 2.*(TH1D(I)-THGB(I))/ZA(I)
   ENDIF
   T2(I)=TH2(I)*(PSFC(I)/100.)**ROVCP

   Q2(I)=QSFCMR(I)+(QV1D(I)-QSFCMR(I))*PSIQ2/PSIQ
   !***  BE CERTAIN THAT THE 2-M MIXING RATIO IS BRACKETED BY
   !***  THE VALUES AT THE SURFACE AND LOWEST MODEL LEVEL.
   IF ((QV1D(I)>QSFCMR(I) .AND. (Q2(I)<QSFCMR(I) .OR. Q2(I)>QV1D(I))) .OR. &
       (QV1D(I)<QSFCMR(I) .AND. (Q2(I)>QSFCMR(I) .OR. Q2(I)<QV1D(I)))) THEN
       Q2(I)=QSFCMR(I) + 2.*(QV1D(I)-QSFCMR(I))/ZA(I)
   ENDIF

   IF ( debug_code ) THEN
      yesno = 0
      IF (HFX(I) > 1200. .OR. HFX(I) < -700.)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            ITER-ITMAX," ITERATIONS",I,J, "HFX: ",HFX(I)
            yesno = 1
      ENDIF
      IF (LH(I)  > 1200. .OR. LH(I)  < -700.)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            ITER-ITMAX," ITERATIONS",I,J, "LH: ",LH(I)
            yesno = 1
      ENDIF
      IF (UST(I) < 0.0 .OR. UST(I) > 4.0 )THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&                                                                     
            ITER-ITMAX," ITERATIONS",I,J, "UST: ",UST(I)
            yesno = 1
      ENDIF
      IF (WSTAR(I)<0.0 .OR. WSTAR(I) > 6.0)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            ITER-ITMAX," ITERATIONS",I,J, "WSTAR: ",WSTAR(I)
            yesno = 1
      ENDIF
      IF (RHO1D(I)<0.0 .OR. RHO1D(I) > 1.6 )THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            ITER-ITMAX," ITERATIONS",I,J, "rho: ",RHO1D(I)
            yesno = 1
      ENDIF
      IF (QSFC(I)*1000. <0.0 .OR. QSFC(I)*1000. >40.)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&
            ITER-ITMAX," ITERATIONS",I,J, "QSFC: ",QSFC(I)
            yesno = 1
      ENDIF
      IF (PBLH(I)<0. .OR. PBLH(I)>6000.)THEN
            print*,"SUSPICIOUS VALUES IN MYNN SFCLAYER",&                                                              
            ITER-ITMAX," ITERATIONS",I,J, "PBLH: ",PBLH(I)
            yesno = 1
      ENDIF

      IF (yesno == 1) THEN
        print*," OTHER INFO:"
        write(*,1001)"REGIME:",REGIME(I)," z/L:",ZOL(I)," U*:",UST(I),&
              " Tstar:",MOL(I)
        write(*,1002)"PSIM:",PSIM(I)," PSIH:",PSIH(I)," W*:",WSTAR(I),&
              " DTHV:",THV1D(I)-THVGB(I)
        write(*,1003)"CPM:",CPM(I)," RHO1D:",RHO1D(I)," L:",&
              ZOL(I)/ZA(I)," DTH:",TH1D(I)-THGB(I)
        write(*,1004)"Z0/Zt:",zratio(I)," Z0:",ZNTstoch(I)," Zt:",z_t(I),&
              " za:",za(I)
        write(*,1005)"Re:",restar," MAVAIL:",MAVAIL(I)," QSFC(I):",&
              QSFC(I)," QVSH(I):",QVSH(I)
        print*,"PSIX=",PSIX," Z0:",ZNTstoch(I)," T1D(i):",T1D(i)
        write(*,*)"============================================="
      ENDIF
   ENDIF

 ENDDO !end i-loop

END SUBROUTINE SFCLAY1D_mynn
!-------------------------------------------------------------------          
   SUBROUTINE zilitinkevich_1995(Z_0,Zt,Zq,restar,ustar,KARMAN,&
       & landsea,IZ0TLND2,spp_pbl,rstoch)

       ! This subroutine returns the thermal and moisture roughness lengths
       ! from Zilitinkevich (1995) and Zilitinkevich et al. (2001) over
       ! land and water, respectively. 
       !
       ! MODS:
       ! 20120705 : added IZ0TLND option. Note: This option was designed
       !            to work with the Noah LSM and may be specific for that
       !            LSM only. Tests with RUC LSM showed no improvements. 

       IMPLICIT NONE
       REAL, INTENT(IN) :: Z_0,restar,ustar,KARMAN,landsea
       INTEGER, OPTIONAL, INTENT(IN)::  IZ0TLND2
       REAL, INTENT(OUT) :: Zt,Zq
       REAL :: CZIL  !=0.100 in Chen et al. (1997)
                     !=0.075 in Zilitinkevich (1995)
                     !=0.500 in Lemone et al. (2008)
       INTEGER,  INTENT(IN)  ::    spp_pbl
       REAL,     INTENT(IN)  ::    rstoch


       IF (landsea-1.5 .GT. 0) THEN    !WATER

          !THIS IS BASED ON Zilitinkevich, Grachev, and Fairall (2001;
          !Their equations 15 and 16).
          IF (restar .LT. 0.1) THEN
             Zt = Z_0*EXP(KARMAN*2.0)
             Zt = MIN( Zt, 6.0e-5)
             Zt = MAX( Zt, 2.0e-9)
             Zq = Z_0*EXP(KARMAN*3.0)
             Zq = MIN( Zq, 6.0e-5)
             Zq = MAX( Zq, 2.0e-9)
          ELSE
             Zt = Z_0*EXP(-KARMAN*(4.0*SQRT(restar)-3.2))
             Zt = MIN( Zt, 6.0e-5)
             Zt = MAX( Zt, 2.0e-9)
             Zq = Z_0*EXP(-KARMAN*(4.0*SQRT(restar)-4.2))
             Zq = MIN( Zt, 6.0e-5)
             Zq = MAX( Zt, 2.0e-9)
          ENDIF

       ELSE                             !LAND

          !Option to modify CZIL according to Chen & Zhang, 2009
          IF ( IZ0TLND2 .EQ. 1 ) THEN
             CZIL = 10.0 ** ( -0.40 * ( Z_0 / 0.07 ) )
          ELSE
             CZIL = 0.075 !0.10
          END IF

          Zt = Z_0*EXP(-KARMAN*CZIL*SQRT(restar))
          Zt = MIN( Zt, Z_0/2.)

          Zq = Z_0*EXP(-KARMAN*CZIL*SQRT(restar))
          Zq = MIN( Zq, Z_0/2.)

! perturb thermal and moisture roughness lenth by +/-50%
! uses same perturbation pattern for perturbing cloud fraction 
! and turbulent mixing length (module_sf_mynn.F), but 
! twice the amplitude; 
! multiplication with -1.0 anticorrelates patterns
          if (spp_pbl==1) then
             Zt = Zt + Zt * 2.0 * rstoch
             Zt = MAX(Zt, 0.001)
             Zq = Zt
          endif

       ENDIF
                   
       return

   END SUBROUTINE zilitinkevich_1995
!--------------------------------------------------------------------
   SUBROUTINE Pan_etal_1994(PSIQ,PSIQ2,ustar,psih,psih2,KARMAN,Z1)

       ! This subroutine returns the resistance (PSIQ) for moisture
       ! exchange. This is a modified form originating from Pan et al. 
       ! (1994) but modified according to tests in both the RUC model 
       ! and WRF-ARW. Note that it is very similar to Carlson and
       ! Boland (1978) model (include below in comments) but has an
       ! extra molecular layer (a third layer) instead of two layers. 

       IMPLICIT NONE
       REAL, INTENT(IN) :: Z1,ustar,KARMAN,psih,psih2
       REAL, INTENT(OUT) :: psiq,psiq2
       REAL, PARAMETER :: Cpan=1.0 !was 20.8 in Pan et al 1994
       REAL, PARAMETER :: ZL=0.01  
       REAL, PARAMETER :: ZMUs=0.2E-3
       REAL, PARAMETER :: XKA = 2.4E-5

         !PAN et al. (1994): 3-layer model, as in paper:
         !ZMU = Cpan*XKA/(KARMAN*UST(I)) 
         !PSIQ =MAX(KARMAN*ustar*ZMU/XKA + LOG((KARMAN*ustar*ZL + XKA)/XKA + &
         !     & Z1/ZL) - PSIH,2.0)
         !PSIQ2=MAX(KARMAN*ustar*ZMU/XKA + LOG((KARMAN*ustar*ZL + XKA)/XKA + &
         !     & 2./ZL) - PSIH2,2.0)                                       
         !MODIFIED FORM:
         PSIQ =MAX(KARMAN*ustar*ZMUs/XKA + LOG((KARMAN*ustar*Z1)/XKA + &
              & Z1/ZL) - PSIH,2.0)
         PSIQ2=MAX(KARMAN*ustar*ZMUs/XKA + LOG((KARMAN*ustar*2.0)/XKA + &
              & 2./ZL) - PSIH2,2.0)

         !CARLSON AND BOLAND (1978): 2-layer model
         !PSIQ =MAX(LOG(KARMAN*ustar*Z1/XKA + Z1/ZL)-PSIH  ,2.0)
         !PSIQ2=MAX(LOG(KARMAN*ustar*2./XKA + 2./ZL)-PSIH2 ,2.0)

    END SUBROUTINE Pan_etal_1994
!--------------------------------------------------------------
   SUBROUTINE davis_etal_2008(Z_0,ustar)

    !a.k.a. : Donelan et al. (2004)
    !This formulation for roughness length was designed to match 
    !the labratory experiments of Donelan et al. (2004).
    !This is an update version from Davis et al. 2008, which
    !corrects a small-bias in Z_0 (AHW real-time 2012).

       IMPLICIT NONE
       REAL, INTENT(IN)  :: ustar
       REAL, INTENT(OUT)  :: Z_0
       REAL :: ZW, ZN1, ZN2
       REAL, PARAMETER :: G=9.81, OZO=1.59E-5

       !OLD FORM: Z_0 = 10.*EXP(-10./(ustar**(1./3.)))
       !NEW FORM:

       ZW  = MIN((ustar/1.06)**(0.3),1.0)
       ZN1 = 0.011*ustar*ustar/G + OZO
       ZN2 = 10.*exp(-9.5*ustar**(-.3333)) + &
             0.11*1.5E-5/AMAX1(ustar,0.01)
       Z_0 = (1.0-ZW) * ZN1 + ZW * ZN2

       Z_0 = MAX( Z_0, 1.27e-7)  !These max/mins were suggested by
       Z_0 = MIN( Z_0, 2.85e-3)  !Davis et al. (2008)
                   
       return

   END SUBROUTINE davis_etal_2008
!--------------------------------------------------------------------
   SUBROUTINE Taylor_Yelland_2001(Z_0,ustar,wsp10)

    !This formulation for roughness length was designed account for 
    !wave steepness.

       IMPLICIT NONE
       REAL, INTENT(IN)  :: ustar,wsp10
       REAL, INTENT(OUT) :: Z_0
       REAL, parameter  :: g=9.81, pi=3.14159265
       REAL :: hs, Tp, Lp

       !hs is the significant wave height
        hs = 0.0248*(wsp10**2.)
       !Tp dominant wave period
        Tp = 0.729*MAX(wsp10,0.1)
       !Lp is the wavelength of the dominant wave
        Lp = g*Tp**2/(2*pi)

       Z_0 = 1200.*hs*(hs/Lp)**4.5
       Z_0 = MAX( Z_0, 1.27e-7)  !These max/mins were suggested by
       Z_0 = MIN( Z_0, 2.85e-3)  !Davis et al. (2008)
                   
       return

   END SUBROUTINE Taylor_Yelland_2001
!--------------------------------------------------------------------
   SUBROUTINE charnock_1955(Z_0,ustar,wsp10,visc,zu)
 
    !This version of Charnock's relation employs a varying
    !Charnock parameter, similar to COARE3.0 [Fairall et al. (2003)].
    !The Charnock parameter CZC is varied from .011 to .018 
    !between 10-m wsp = 10 and 18. 

       IMPLICIT NONE
       REAL, INTENT(IN)  :: ustar, visc, wsp10, zu
       REAL, INTENT(OUT) :: Z_0
       REAL, PARAMETER   :: G=9.81, CZO2=0.011
       REAL              :: CZC    !variable charnock "constant"   
       REAL              :: wsp10m ! logarithmically calculated 10 m

       wsp10m = wsp10*log(10./1e-4)/log(zu/1e-4)
       CZC = CZO2 + 0.007*MIN(MAX((wsp10m-10.)/8., 0.), 1.0)

       Z_0 = CZC*ustar*ustar/G + (0.11*visc/MAX(ustar,0.05))
       Z_0 = MAX( Z_0, 1.27e-7)  !These max/mins were suggested by
       Z_0 = MIN( Z_0, 2.85e-3)  !Davis et al. (2008)

       return

   END SUBROUTINE charnock_1955
!--------------------------------------------------------------------
   SUBROUTINE edson_etal_2013(Z_0,ustar,wsp10,visc,zu)

     !This version of Charnock's relation employs a varying
     !Charnock parameter, taken from COARE 3.5 [Edson et al. (2001, JPO)].
     !The Charnock parameter CZC is varied from about .005 to .028
     !between 10-m wind speeds of 6 and 19 m/s.

       IMPLICIT NONE
       REAL, INTENT(IN)  :: ustar, visc, wsp10, zu
       REAL, INTENT(OUT) :: Z_0
       REAL, PARAMETER   :: G=9.81
       REAL, PARAMETER   :: m=0.017, b=-0.005
       REAL              :: CZC    ! variable charnock "constant"
       REAL              :: wsp10m ! logarithmically calculated 10 m

       wsp10m = wsp10*log(10/1e-4)/log(zu/1e-4)
       wsp10m = MIN(19., wsp10m)
       CZC    = m*wsp10m + b
       CZC    = MAX(CZC, 0.0)

       Z_0 = CZC*ustar*ustar/G + (0.11*visc/MAX(ustar,0.07))
       Z_0 = MAX( Z_0, 1.27e-7)  !These max/mins were suggested by
       Z_0 = MIN( Z_0, 2.85e-3)  !Davis et al. (2008)

       return

   END SUBROUTINE edson_etal_2013
!--------------------------------------------------------------------
   SUBROUTINE garratt_1992(Zt,Zq,Z_0,Ren,landsea)

    !This formulation for the thermal and moisture roughness lengths
    !(Zt and Zq) relates them to Z0 via the roughness Reynolds number (Ren).
    !This formula comes from Fairall et al. (2003). It is modified from
    !the original Garratt-Brutsaert model to better fit the COARE/HEXMAX
    !data. The formula for land uses a constant ratio (Z_0/7.4) taken
    !from Garratt (1992).

       IMPLICIT NONE
       REAL, INTENT(IN)  :: Ren, Z_0,landsea
       REAL, INTENT(OUT) :: Zt,Zq
       REAL :: Rq
       REAL, PARAMETER  :: e=2.71828183

       IF (landsea-1.5 .GT. 0) THEN    !WATER

          Zt = Z_0*EXP(2.0 - (2.48*(Ren**0.25)))
          Zq = Z_0*EXP(2.0 - (2.28*(Ren**0.25)))

          Zq = MIN( Zq, 5.5e-5)
          Zq = MAX( Zq, 2.0e-9)
          Zt = MIN( Zt, 5.5e-5)
          Zt = MAX( Zt, 2.0e-9) !same lower limit as ECMWF
       ELSE                            !LAND
          Zq = Z_0/(e**2.)      !taken from Garratt (1980,1992)
          Zt = Zq
       ENDIF
                   
       return

    END SUBROUTINE garratt_1992
!--------------------------------------------------------------------
    SUBROUTINE fairall_etal_2003(Zt,Zq,Ren,ustar,visc)

    !This formulation for thermal and moisture roughness length (Zt and Zq)
    !as a function of the roughness Reynolds number (Ren) comes from the
    !COARE3.0 formulation, empirically derived from COARE and HEXMAX data
    ![Fairall et al. (2003)]. Edson et al. (2004; JGR) suspected that this
    !relationship overestimated the scalar roughness lengths for low Reynolds
    !number flows, so an optional smooth flow relationship, taken from Garratt
    !(1992, p. 102), is available for flows with Ren < 2.
    !
    !This is for use over water only.

       IMPLICIT NONE
       REAL, INTENT(IN)  :: Ren,ustar,visc
       REAL, INTENT(OUT) :: Zt,Zq

       IF (Ren .le. 2.) then

          Zt = (5.5e-5)*(Ren**(-0.60))
          Zq = Zt
          !FOR SMOOTH SEAS, CAN USE GARRATT
          !Zq = 0.2*visc/MAX(ustar,0.1)
          !Zq = 0.3*visc/MAX(ustar,0.1)

       ELSE

          !FOR ROUGH SEAS, USE COARE
          Zt = (5.5e-5)*(Ren**(-0.60))
          Zq = Zt

       ENDIF

       Zt = MIN(Zt,1.0e-4)
       Zt = MAX(Zt,2.0e-9)

       Zq = MIN(Zt,1.0e-4)
       Zq = MAX(Zt,2.0e-9)

       return

    END SUBROUTINE fairall_etal_2003
!--------------------------------------------------------------------
    SUBROUTINE fairall_etal_2014(Zt,Zq,Ren,ustar,visc,rstoch,spp_pbl)

    !This formulation for thermal and moisture roughness length (Zt and Zq)
    !as a function of the roughness Reynolds number (Ren) comes from the
    !COARE 3.5/4.0 formulation, empirically derived from COARE and HEXMAX data
    ![Fairall et al. (2014? coming soon, not yet published as of July 2014)].
    !This is for use over water only.

       IMPLICIT NONE
       REAL, INTENT(IN)  :: Ren,ustar,visc,rstoch
       INTEGER, INTENT(IN):: spp_pbl
       REAL, INTENT(OUT) :: Zt,Zq

       !Zt = (5.5e-5)*(Ren**(-0.60))
       Zt = MIN(1.6E-4, 5.8E-5/(Ren**0.72))
       Zq = Zt

       IF (spp_pbl ==1) THEN
          Zt = MAX(Zt + Zt*2.0*rstoch,2.0e-9)
          Zq = MAX(Zt + Zt*2.0*rstoch,2.0e-9)
       ELSE
          Zt = MAX(Zt,2.0e-9)
          Zq = MAX(Zt,2.0e-9)
       ENDIF

       return

    END SUBROUTINE fairall_etal_2014
!--------------------------------------------------------------------
    SUBROUTINE Yang_2008(Z_0,Zt,Zq,ustar,tstar,qst,Ren,visc,landsea)

    !This is a modified version of Yang et al (2002 QJRMS, 2008 JAMC) 
    !and Chen et al (2010, J of Hydromet). Although it was originally 
    !designed for arid regions with bare soil, it is modified 
    !here to perform over a broader spectrum of vegetation.
    !
    !The original formulation relates the thermal roughness length (Zt) 
    !to u* and T*:
    !  
    ! Zt = ht * EXP(-beta*(ustar**0.5)*(ABS(tstar)**0.25))
    !
    !where ht = Renc*visc/ustar and the critical Reynolds number 
    !(Renc) = 70. Beta was originally = 10 (2002 paper) but was revised 
    !to 7.2 (in 2008 paper). Their form typically varies the
    !ratio Z0/Zt by a few orders of magnitude (1-1E4).
    !
    !This modified form uses beta = 1.5 and a variable Renc (function of Z_0),
    !so zt generally varies similarly to the Zilitinkevich form (with Czil = 0.1)
    !for very small or negative surface heat fluxes but can become close to the
    !Zilitinkevich with Czil = 0.2 for very large HFX (large negative T*).
    !Also, the exponent (0.25) on tstar was changed to 1.0, since we found
    !Zt was reduced too much for low-moderate positive heat fluxes.
    !
    !This should only be used over land!

       IMPLICIT NONE
       REAL, INTENT(IN)  :: Z_0, Ren, ustar, tstar, qst, visc, landsea
       REAL              :: ht,     &! roughness height at critical Reynolds number
                            tstar2, &! bounded T*, forced to be non-positive
                            qstar2, &! bounded q*, forced to be non-positive
                            Z_02,   &! bounded Z_0 for variable Renc2 calc
                            Renc2    ! variable Renc, function of Z_0
       REAL, INTENT(OUT) :: Zt,Zq
       REAL, PARAMETER  :: Renc=300., & !old constant Renc
                           beta=1.5,  & !important for diurnal variation
                           m=170.,    & !slope for Renc2 function
                           b=691.       !y-intercept for Renc2 function

       Z_02 = MIN(Z_0,0.5)
       Z_02 = MAX(Z_02,0.04)
       Renc2= b + m*log(Z_02)
       ht     = Renc2*visc/MAX(ustar,0.01)
       tstar2 = MIN(tstar, 0.0)
       qstar2 = MIN(qst,0.0)

       Zt     = ht * EXP(-beta*(ustar**0.5)*(ABS(tstar2)**1.0))
       Zq     = ht * EXP(-beta*(ustar**0.5)*(ABS(qstar2)**1.0))
       !Zq     = Zt

       Zt = MIN(Zt, Z_0/2.0)
       Zq = MIN(Zq, Z_0/2.0)

       return

    END SUBROUTINE Yang_2008
!--------------------------------------------------------------------
    SUBROUTINE Andreas_2002(Z_0,bvisc,ustar,Zt,Zq)

    ! This is taken from Andreas (2002; J. of Hydromet) and 
    ! Andreas et al. (2005; BLM).
    !
    ! This should only be used over snow/ice!

       IMPLICIT NONE
       REAL, INTENT(IN)  :: Z_0, bvisc, ustar
       REAL, INTENT(OUT) :: Zt, Zq
       REAL :: Ren2, zntsno

       REAL, PARAMETER  :: bt0_s=1.25,  bt0_t=0.149,  bt0_r=0.317,  &
                           bt1_s=0.0,   bt1_t=-0.55,  bt1_r=-0.565, &
                           bt2_s=0.0,   bt2_t=0.0,    bt2_r=-0.183

       REAL, PARAMETER  :: bq0_s=1.61,  bq0_t=0.351,  bq0_r=0.396,  &
                           bq1_s=0.0,   bq1_t=-0.628, bq1_r=-0.512, &
                           bq2_s=0.0,   bq2_t=0.0,    bq2_r=-0.180

      !Calculate zo for snow (Andreas et al. 2005, BLM)                                                                     
       zntsno = 0.135*bvisc/ustar + &
               (0.035*(ustar*ustar)/9.8) * &
               (5.*exp(-1.*(((ustar - 0.18)/0.1)*((ustar - 0.18)/0.1))) + 1.)                                                
       Ren2 = ustar*zntsno/bvisc

       ! Make sure that Re is not outside of the range of validity
       ! for using their equations
       IF (Ren2 .gt. 1000.) Ren2 = 1000. 

       IF (Ren2 .le. 0.135) then

          Zt = zntsno*EXP(bt0_s + bt1_s*LOG(Ren2) + bt2_s*LOG(Ren2)**2)
          Zq = zntsno*EXP(bq0_s + bq1_s*LOG(Ren2) + bq2_s*LOG(Ren2)**2)

       ELSE IF (Ren2 .gt. 0.135 .AND. Ren2 .lt. 2.5) then

          Zt = zntsno*EXP(bt0_t + bt1_t*LOG(Ren2) + bt2_t*LOG(Ren2)**2)
          Zq = zntsno*EXP(bq0_t + bq1_t*LOG(Ren2) + bq2_t*LOG(Ren2)**2)

       ELSE

          Zt = zntsno*EXP(bt0_r + bt1_r*LOG(Ren2) + bt2_r*LOG(Ren2)**2)
          Zq = zntsno*EXP(bq0_r + bq1_r*LOG(Ren2) + bq2_r*LOG(Ren2)**2)

       ENDIF

       return

    END SUBROUTINE Andreas_2002
!--------------------------------------------------------------------
    SUBROUTINE PSI_Hogstrom_1996(psi_m, psi_h, zL, Zt, Z_0, Za)

    ! This subroutine returns the stability functions based off
    ! of Hogstrom (1996).

       IMPLICIT NONE
       REAL, INTENT(IN)  :: zL, Zt, Z_0, Za
       REAL, INTENT(OUT) :: psi_m, psi_h
       REAL  :: x, x0, y, y0, zmL, zhL

       zmL = Z_0*zL/Za  
       zhL = Zt*zL/Za

       IF (zL .gt. 0.) THEN  !STABLE (not well tested - seem large)

          psi_m = -5.3*(zL - zmL)
          psi_h = -8.0*(zL - zhL)
 
       ELSE                 !UNSTABLE

          x = (1.-19.0*zL)**0.25
          x0= (1.-19.0*zmL)**0.25
          y = (1.-11.6*zL)**0.5
          y0= (1.-11.6*zhL)**0.5

          psi_m = 2.*LOG((1.+x)/(1.+x0)) + &
                    &LOG((1.+x**2.)/(1.+x0**2.)) - &
                    &2.0*ATAN(x) + 2.0*ATAN(x0)
          psi_h = 2.*LOG((1.+y)/(1.+y0))

       ENDIF
                   
       return

    END SUBROUTINE PSI_Hogstrom_1996
!--------------------------------------------------------------------
    SUBROUTINE PSI_DyerHicks(psi_m, psi_h, zL, Zt, Z_0, Za)

    ! This subroutine returns the stability functions based off
    ! of Hogstrom (1996), but with different constants compatible
    ! with Dyer and Hicks (1970/74?). This formulation is used for
    ! testing/development by Nakanishi (personal communication).

       IMPLICIT NONE
       REAL, INTENT(IN)  :: zL, Zt, Z_0, Za
       REAL, INTENT(OUT) :: psi_m, psi_h
       REAL  :: x, x0, y, y0, zmL, zhL

       zmL = Z_0*zL/Za  !Zo/L
       zhL = Zt*zL/Za   !Zt/L

       IF (zL .gt. 0.) THEN  !STABLE

          psi_m = -5.0*(zL - zmL)
          psi_h = -5.0*(zL - zhL)
 
       ELSE                 !UNSTABLE

          x = (1.-16.*zL)**0.25
          x0= (1.-16.*zmL)**0.25

          y = (1.-16.*zL)**0.5
          y0= (1.-16.*zhL)**0.5

          psi_m = 2.*LOG((1.+x)/(1.+x0)) + &
                    &LOG((1.+x**2.)/(1.+x0**2.)) - & 
                    &2.0*ATAN(x) + 2.0*ATAN(x0)
          psi_h = 2.*LOG((1.+y)/(1.+y0))

       ENDIF
                   
       return

    END SUBROUTINE PSI_DyerHicks
!--------------------------------------------------------------------
    SUBROUTINE PSI_CB2005(psi_m, psi_h, zL, Zt, Z_0, Za)

    ! This subroutine returns the stability functions based off
    ! of Hogstrom (1996), but with different constants compatible
    ! with Dyer and Hicks (1970/74?). This formulation is used for
    ! testing/development by Nakanishi (personal communication).

       IMPLICIT NONE
       REAL, INTENT(IN)  :: zL, Zt, Z_0, Za
       REAL, INTENT(OUT) :: psi_m, psi_h
       REAL  :: x, x0, y, y0, zmL, zhL

       zmL = Z_0*zL/Za  !Zo/L
       zhL = Zt*zL/Za   !Zt/L

       IF (zL .gt. 0.) THEN  !STABLE

!          psi_m = -5.0*(zL - zmL)
!          psi_h = -5.0*(zL - zhL)
           x = zL+(1.0+zL**2.5)**(1.0/2.5)
           x0 = zmL+(1.0+zmL**2.5)**(1.0/2.5)

           y = zL+(1.0+zL**1.1)**(1.0/1.1)
           y0 = zhL+(1.0+zhL**1.1)**(1.0/1.1)

           psi_m = -6.1*LOG(x/x0)
           psi_h = -5.3*LOG(y/y0)
       ELSE                 !UNSTABLE

          x = (1.-16.*zL)**0.25
          x0= (1.-16.*zmL)**0.25

          y = (1.-16.*zL)**0.5
          y0= (1.-16.*zhL)**0.5

          psi_m = 2.*LOG((1.+x)/(1.+x0)) + &
                    &LOG((1.+x**2.)/(1.+x0**2.)) - &
                    &2.0*ATAN(x) + 2.0*ATAN(x0)
          psi_h = 2.*LOG((1.+y)/(1.+y0))

       ENDIF

       return

    END SUBROUTINE PSI_CB2005
!--------------------------------------------------------------------

    SUBROUTINE PSI_G2007(psi_m, psi_h, zL, Zt, Z_0, Za)

    ! This subroutine returns the stability functions based off
    ! of Hogstrom (1996), but with different constants compatible
    ! with Dyer and Hicks (1970/74?). This formulation is used for
    ! testing/development by Nakanishi (personal communication).

       IMPLICIT NONE
       REAL, INTENT(IN)  :: zL, Zt, Z_0, Za
       REAL, INTENT(OUT) :: psi_m, psi_h
       REAL  :: x, x0, y, y0, zmL, zhL
       REAL  ::  psim_zol, psim_0, psih_zol, psih_0
       REAL  ::  partx_log, partx0_log, partx_atan, partx0_atan

       zmL = Z_0*zL/Za  !Zo/L
       zhL = Zt*zL/Za   !Zt/L

       IF (zL .gt. 0.) THEN  !STABLE
!          psi_m = -5.0*(zL - zmL)
!          psi_h = -5.0*(zL - zhL)
           x = (1.0+zL)**(1.0/3.0)
           x0 = (1.0+zmL)**(1.0/3.0)

           partx_log = 2.0*LOG((x+B_m)/(1.0+B_m))-LOG((x**2.0-B_m*x+B_m**2.0)/cm_c)
           partx0_log = 2.0*LOG((x0+B_m)/(1.0+B_m))-LOG((x0**2.0-B_m*x0+B_m**2.0)/cm_c)
           partx_atan = 2.0*3.0**0.5*(ATAN((2.0*x-B_m)/B_m/3.0**0.5)-ATAN(cm_d))
           partx0_atan = 2.0*3.0**0.5*(ATAN((2.0*x0-B_m)/B_m/3.0**0.5)-ATAN(cm_d))
           psim_zol = cm_a*(x-1.0)+cm_b*(partx_log+partx_atan)
           psim_0 = cm_a*(x0-1.0)+cm_b*(partx0_log+partx0_atan)
           psi_m = psim_zol-psim_0

           psih_zol = -2.5*LOG(zL**2.0+C_h*zL+1.0)+ch_a*LOG(((2.0*zL+C_h-B_h)/(2.0*zL+C_h+B_h))/ch_d)
           psih_0 = -2.5*LOG(zhL**2.0+C_h*zhL+1.0)+ch_a*LOG(((2.0*zhL+C_h-B_h)/(2.0*zhL+C_h+B_h))/ch_d)
           psi_h = psih_zol-psih_0

       ELSE                 !UNSTABLE

          x = (1.-16.*zL)**0.25
          x0= (1.-16.*zmL)**0.25

          y = (1.-16.*zL)**0.5
          y0= (1.-16.*zhL)**0.5

          psi_m = 2.*LOG((1.+x)/(1.+x0)) + &
                    &LOG((1.+x**2.)/(1.+x0**2.)) - &
                    &2.0*ATAN(x) + 2.0*ATAN(x0)
          psi_h = 2.*LOG((1.+y)/(1.+y0))

       ENDIF

       return

    END SUBROUTINE PSI_G2007

!--------------------------------------------------------------------
    SUBROUTINE PSI_Beljaars_Holtslag_1991(psi_m, psi_h, zL)

    ! This subroutine returns the stability functions based off
    ! of Beljaar and Holtslag 1991, which is an extension of Holtslag
    ! and Debruin 1989.

       IMPLICIT NONE
       REAL, INTENT(IN)  :: zL
       REAL, INTENT(OUT) :: psi_m, psi_h
       REAL, PARAMETER  :: a=1., b=0.666, c=5., d=0.35

       IF (zL .lt. 0.) THEN  !UNSTABLE

          WRITE(*,*)"WARNING: Universal stability functions from"
          WRITE(*,*)"        Beljaars and Holtslag (1991) should only"
          WRITE(*,*)"        be used in the stable regime!"
          psi_m = 0.
          psi_h = 0.
 
       ELSE                 !STABLE

          psi_m = -(a*zL + b*(zL -(c/d))*exp(-d*zL) + (b*c/d))
          psi_h = -((1.+.666*a*zL)**1.5 + &
                  b*(zL - (c/d))*exp(-d*zL) + (b*c/d) -1.)

       ENDIF
                   
       return

    END SUBROUTINE PSI_Beljaars_Holtslag_1991
!--------------------------------------------------------------------
    SUBROUTINE PSI_Zilitinkevich_Esau_2007(psi_m, psi_h, zL)

    ! This subroutine returns the stability functions come from
    ! Zilitinkevich and Esau (2007, BM), which are formulatioed from the
    ! "generalized similarity theory" and tuned to the LES DATABASE64
    ! to determine their dependence on z/L.

       IMPLICIT NONE
       REAL, INTENT(IN)  :: zL
       REAL, INTENT(OUT) :: psi_m, psi_h
       REAL, PARAMETER  :: Cm=3.0, Ct=2.5

       IF (zL .lt. 0.) THEN  !UNSTABLE

          WRITE(*,*)"WARNING: Universal stability function from"
          WRITE(*,*)"        Zilitinkevich and Esau (2007) should only"
          WRITE(*,*)"        be used in the stable regime!"
          psi_m = 0.
          psi_h = 0.
 
       ELSE                 !STABLE

          psi_m = -Cm*(zL**(5./6.))
          psi_h = -Ct*(zL**(4./5.))

       ENDIF
                   
       return

    END SUBROUTINE PSI_Zilitinkevich_Esau_2007
!--------------------------------------------------------------------
    SUBROUTINE PSI_Businger_1971(psi_m, psi_h, zL)

    ! This subroutine returns the flux-profile relationships
    ! of Businger el al. 1971.

       IMPLICIT NONE
       REAL, INTENT(IN)  :: zL
       REAL, INTENT(OUT) :: psi_m, psi_h
       REAL  :: x, y
       REAL, PARAMETER  ::  Pi180 = 3.14159265/180.

       IF (zL .lt. 0.) THEN  !UNSTABLE

          x = (1. - 15.0*zL)**0.25
          y = (1. - 9.0*zL)**0.5

          psi_m = LOG(((1.+x)/2.)**2.) + &
                 &LOG((1.+x**2.)/2.) - &
                 &2.0*ATAN(x) + Pi180*90.
          psi_h = 2.*LOG((1.+y)/2.)

       ELSE                 !STABLE

          psi_m = -4.7*zL
          psi_h = -(4.7/0.74)*zL

       ENDIF
                   
       return

    END SUBROUTINE PSI_Businger_1971
!--------------------------------------------------------------------
    SUBROUTINE PSI_Suselj_Sood_2010(psi_m, psi_h, zL)

    !This subroutine returns flux-profile relatioships based off
    !of Lobocki (1993), which is derived from the MY-level 2 model.
    !Suselj and Sood (2010) applied the surface layer length scales
    !from Nakanishi (2001) to get this new relationship. These functions
    !are more agressive (larger magnitude) than most formulations. They
    !showed improvement over water, but untested over land.

       IMPLICIT NONE
       REAL, INTENT(IN)  :: zL
       REAL, INTENT(OUT) :: psi_m, psi_h
       REAL, PARAMETER  :: Rfc=0.19, Ric=0.183, PHIT=0.8

       IF (zL .gt. 0.) THEN  !STABLE

          psi_m = -(zL/Rfc + 1.1223*EXP(1.-1.6666/zL))
          !psi_h = -zL*Ric/((Rfc**2.)*PHIT) + 8.209*(zL**1.1091)
          !THEIR EQ FOR PSI_H CRASHES THE MODEL AND DOES NOT MATCH
          !THEIR FIG 1. THIS EQ (BELOW) MATCHES THEIR FIG 1 BETTER:
          psi_h = -(zL*Ric/((Rfc**2.)*5.) + 7.09*(zL**1.1091))
 
       ELSE                 !UNSTABLE

          psi_m = 0.9904*LOG(1. - 14.264*zL)
          psi_h = 1.0103*LOG(1. - 16.3066*zL)

       ENDIF
                   
       return

    END SUBROUTINE PSI_Suselj_Sood_2010
!--------------------------------------------------------------------
    SUBROUTINE Li_etal_2010(zL, Rib, zaz0, z0zt)

    !This subroutine returns a more robust z/L that best matches
    !the z/L from Hogstrom (1996) for unstable conditions and Beljaars
    !and Holtslag (1991) for stable conditions.

       IMPLICIT NONE
       REAL, INTENT(OUT)  :: zL
       REAL, INTENT(IN) :: Rib, zaz0, z0zt
       REAL :: alfa, beta, zaz02, z0zt2
       REAL, PARAMETER  :: au11=0.045, bu11=0.003, bu12=0.0059, &
                          &bu21=-0.0828, bu22=0.8845, bu31=0.1739, &
                          &bu32=-0.9213, bu33=-0.1057
       REAL, PARAMETER  :: aw11=0.5738, aw12=-0.4399, aw21=-4.901,&
                          &aw22=52.50, bw11=-0.0539, bw12=1.540, &
                          &bw21=-0.669, bw22=-3.282
       REAL, PARAMETER  :: as11=0.7529, as21=14.94, bs11=0.1569,&
                          &bs21=-0.3091, bs22=-1.303
          
       !set limits according to Li et al (2010), p 157.
       zaz02=zaz0
       IF (zaz0 .lt. 100.0) zaz02=100.
       IF (zaz0 .gt. 100000.0) zaz02=100000.

       !set more limits according to Li et al (2010)
       z0zt2=z0zt
       IF (z0zt .lt. 0.5) z0zt2=0.5
       IF (z0zt .gt. 100.0) z0zt2=100.

       alfa = LOG(zaz02)
       beta = LOG(z0zt2)

       IF (Rib .le. 0.0) THEN
          zL = au11*alfa*Rib**2 + (                   &
               &  (bu11*beta + bu12)*alfa**2 +        &
               &  (bu21*beta + bu22)*alfa    +        &
               &  (bu31*beta**2 + bu32*beta + bu33))*Rib
          !if(zL .LT. -15 .OR. zl .GT. 0.)print*,"VIOLATION Rib<0:",zL
          zL = MAX(zL,-15.) !LIMITS SET ACCORDING TO Li et al (2010)
          zL = MIN(zL,0.)   !Figure 1.
       ELSEIF (Rib .gt. 0.0 .AND. Rib .le. 0.2) THEN
          zL = ((aw11*beta + aw12)*alfa +             &
             &  (aw21*beta + aw22))*Rib**2 +          &
             & ((bw11*beta + bw12)*alfa +             &
             &  (bw21*beta + bw22))*Rib
          !if(zL .LT. 0 .OR. zl .GT. 4)print*,"VIOLATION 0<Rib<0.2:",zL
          zL = MIN(zL,4.) !LIMITS APPROX SET ACCORDING TO Li et al (2010)
          zL = MAX(zL,0.) !THEIR FIGURE 1B.
       ELSE
          zL = (as11*alfa + as21)*Rib + bs11*alfa +   &
             &  bs21*beta + bs22
          !if(zL .LE. 1 .OR. zl .GT. 23)print*,"VIOLATION Rib>0.2:",zL
          zL = MIN(zL,20.) !LIMITS ACCORDING TO Li et al (2010), THIER
                           !FIGUE 1C.
          zL = MAX(zL,1.)
       ENDIF

       return

    END SUBROUTINE Li_etal_2010
SUBROUTINE Li_etal_2014(zol,Rib,zz,z0m,z0h)

 IMPLICIT NONE

 REAL, INTENT(OUT)  :: zol
 REAL, INTENT(IN) :: Rib, zz, z0m, z0h
! REAL :: zol_stable, zol_unstable
 if (Rib.GE.0) then
    zol=zol_stable(Rib,zz,z0m,z0h)
 else
    zol=zol_unstable(Rib,zz,z0m,z0h)
 end if
return

 END SUBROUTINE Li_etal_2014

REAL FUNCTION zol_stable(Rib,zz,z0m,z0h) 
IMPLICIT NONE
real :: Rib,zz,z0m,z0h,ZL, &
        x11,x12,x21,x22,s11,s12,s13,s21,s22,s23,s31,s32,s33, &
        mz,mzn

real, dimension(8) :: XX1
real, dimension(32) :: SS1
! a:1 b:2 c:3 d:4 e:5 f:6 g:7 h:8
real, dimension(6,8) :: coef_ribc_a , coef_ribc_g
real, dimension(3,8) :: coef_ribc_b , coef_ribc_d , coef_ribc_f , coef_ribc_h
real, dimension(4,8) :: coef_ribc_c
real, dimension(5,8) :: coef_ribc_e

real, dimension(7,32) :: coef_zl_a , coef_zl_g
real, dimension(4,32) :: coef_zl_b , coef_zl_d , coef_zl_f , coef_zl_h
real, dimension(5,32) :: coef_zl_c
real, dimension(6,32) :: coef_zl_e

real, dimension(6) :: Ribc_a , Ribc_g
real, dimension(3) :: Ribc_b , Ribc_d , Ribc_f , Ribc_h
real, dimension(4) :: Ribc_c
real, dimension(5) :: Ribc_e

INTEGER :: i
mz=zz/z0m
mzn=z0m/z0h               
x11 = log(log(mz))
x12 = x11**2
x21 = log(mzn)
x22 = x21**2

s21 = log(mz)
s22 = s21**2
s23 = s21**3
s31 = log(mzn)
s32 = s31**2
s33 = s31**3

s11 = Rib
s12 = s11**2
s13 = s11**3

XX1 = (/1.0, x11, x12, x21, x11*x21, x12*x21, x22, x11*x22/)
SS1 = (/1.0, s11, s12, s13, s31, s11*s31, s12*s31, s13*s31, &
        s32, s11*s32, s12*s32, s33, s11*s33, s21, s11*s21, &
        s12*s21, s13*s21, s21*s31, s11*s21*s31,          &
        s12*s21*s31, s21*s32, s11*s21*s32, s21*s33,     &
        s22, s11*s22, s12*s22, s22*s31, s11*s22*s31,     &
        s22*s32, s23, s11*s23, s23*s31/)

if (mz.LT.160 .and. mz.GE.10 .and. mzn.GE.0.607 .and.mzn.LT.100) then !region1***************
    
coef_ribc_a=reshape((/0.3095 , 0.3219 , 0.3545 , 0.4390 , 0.6887 , 1.706 , &
                   -0.2852 , -0.2613 , -0.2569 , -0.3133 , -0.5375 , -1.620 , &
                   0.07955 , 0.06753 , 0.06609 , 0.08619 , 0.1616 , 0.5231 , &
                   0.03388 , 0.04838 , 0.05837 , 0.08930 , 0.1754 , 0.5124 , &
                   -0.01605 , -0.03101 , -0.03934 , -0.07112 , -0.1564 , -0.5026 , &
                   0.0 , 0.003908 , 0.005643 , 0.01403 , 0.03489 , 0.1239 , &
                   0.0 , -0.001780 , -0.003381 , -0.005965 , -0.01277 , -0.03577 , &
                   -0.0001079 , 0.001165 , 0.002194 , 0.003806 , 0.008101 , 0.02238/), shape(coef_ribc_a))
    do i=1,6
        Ribc_a(i) = sum(coef_ribc_a(i,:)*XX1)
    end do
    coef_zl_a=reshape((/-1.134   , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , &
                      31.1 ,   86.35   , -280.4   , 0.0 , 0.0 , -17.32   ,-6.343 , &
                      -71.16   , 0.0 , 3235.0 , 0.0 , 0.0 , 8.773   , 7.66 , &
                      227.4   , 0.0 , -6165.0 , 0.0 , 0.0 , 0.0 , -0.7661  , &
                      -0.2094   , -11.53   , -10.64   , 0.0 , 0.0 , 0.0 , 0.0125 ,  &
                      3.293   , 194.9   , 193.8   , 0.0 , 1.113   , 0.0 , -2.203 , & 
                      -20.11   , -975.4   , -1194.0 , -12.37   , -97.56 , 0.0 , 0.8896  , &
                      14.42   , 1472.0 , 2161.0 , 0.0 , 159.4   , 0.0 , -0.1273 , &
                      0.1476   , -2.535   , -4.603   , 0.0 , 0.0 , 1.919   , -0.008274  , &
                      -0.07325   , 28.24   , 52.02   , 11.99   , 16.33 , 0.0 , 0.3327  , &
                      0.5627   , -61.13   , -110.7   , -15.63   , -25.67 , 0.2679   , -0.04613  , &
                      -0.01178   , -0.2378   , -0.5367   , -0.3157   , -0.6447 , -0.2892   , 0.0 , &
                      0.0218 ,  0.7405   , 1.503   , 0.2948   , 0.9718 , 0.0 , -0.04968 , & 
                      1.405   , 13.6 ,   30.26   , 0.0 , 6.821   , 10.27 , 7.513 , & 
                      -32.47   , -316.2   , -314.9   , 0.0 , -57.13   , 0.0 , 0.0 , &
                      46.59   , 1067.0 , 186.0 , -108.1   , 227.3   , 0.0 , -4.799  , &
                      -38.25   , -1494.0 , 0.0 , 317.8   , -244.0 , 0.0 , 0.5598 , &
                      -0.2286   , 8.023   , 9.038   , 0.0 , 0.9287   , -3.457 , -1.612  , &
                      -1.097   , -91.31   , -87.06   , -12.52   , -17.88 , -1.617   , 0.0 , &
                      -0.3394   , 213.7   , 198.6   , 0.0 , 34.41   , 0.0 , 0.0 , &
                      0.0 , 1.035   , 1.529   , 0.0 , 0.319 ,   -0.07536 , 0.4666 , &  
                      0.0 , -5.072   , -7.439   , -1.025   , -2.452   , 0.0 , 0.0605 ,  &
                      0.0 , 0.03622   , 0.07369   , 0.04669   , 0.08583 , 0.05146   , -0.01808  , & 
                      0.0 , -4.699   , -10.71   , -1.896   , -2.195   , -3.108 , 0.0 , &
                      10.71   , 97.46   , 122.1   , 28.39   , 22.21   , 7.948 , 2.442  , & 
                      0.0 , -152.4   , -76.91   , -14.19   , -31.44   , -2.985 , 0.1584  , & 
                      0.0 , -1.704   , -2.035   , 0.0 , -0.1355   , 0.8751 , 0.0 , &
                      0.0 , 9.069   , 8.248   , 2.214   , 1.976   , 0.3139 , -0.04377 , & 
                      0.0 , -0.09576   , -0.1263   , -0.01472   , -0.04636 , -0.05131   , -0.0694 ,  &
                      -0.007485   , 0.4446   , 1.015   , 0.3069   , 0.1708   , 0.2598   , -0.1675  , &
                      -0.9671   , -7.991   , -10.96   , -3.635   , -1.623   , -0.8513   , -0.2181  , &
                      0.003402   , 0.1138   , 0.1426   , -0.008769   , 0.0 , -0.05427   , 0.05052/), shape(coef_ZL_a))
    if (Rib.LE.Ribc_a(1)) then
        ZL = sum(coef_zl_a(1,:)*SS1)
    else if (Rib.LE.Ribc_a(2)) then
        ZL = sum(coef_zl_a(2,:)*SS1)
    else if (Rib.LE.Ribc_a(3)) then
        ZL = sum(coef_zl_a(3,:)*SS1)
    else if (Rib.LE.Ribc_a(4)) then
        ZL = sum(coef_zl_a(4,:)*SS1)
    else if (Rib.LE.Ribc_a(5)) then
        ZL = sum(coef_zl_a(5,:)*SS1)
    else if (Rib.LE.Ribc_a(6)) then
        ZL = sum(coef_zl_a(6,:)*SS1)
    else if (Rib.GT.Ribc_a(6)) then
        ZL = sum(coef_zl_a(7,:)*SS1)
    end if
else if (mz.LE.10**5 .and. mz.GE.160 .and. mzn.GE.0.607 .and.mzn.LT.100) then  !region 2**************************
       coef_ribc_b=reshape((/0.0 , 0.2002 , 0.4499 ,&
                             0.08606  , 0.0 , 0.0 ,&
                            -0.03048  , -0.01589  , -0.02397  ,&
                             0.09019  , 0.0 , 0.0388  ,&
                            -0.07682  , 0.00367  , -0.01145  ,&
                             0.01693  , 0.0 , 0.0 ,&
                                    0.0 , 0.005057  , 0.0 ,&
                                    0.0 , -0.002399  , 0.0 /),shape(coef_ribc_b))
    do i=1,3
        Ribc_b(i) = sum(coef_ribc_b(i,:)*XX1)
    end do
    coef_zl_b=reshape((/ 0.0 , 0.0 , 0.0 , 0.0 ,&
                       0.0 , 0.0 , 41.53 , 0.0 ,&
                       0.0 , 0.0 , 0.0 , 0.0 ,&
                       0.0 , 0.0 , 0.0 , 0.0 ,&
                       0.0 , 0.0 , -1.616 , -2.570 ,&
                       0.0 , -12.35 , 0.0 , -2.910 ,&
                       0.0 , 0.0 , 0.0 , 0.0 ,&
                       0.0 , 0.0 , 0.0 , 0.0 ,&
                       0.0 , 0.0 , 0.0 , 0.8740 ,&
                       0.0 , 0.5183 , 0.0 , 0.3377 ,&
                       0.0 , 0.0 , 0.0 , 0.0 ,&
                       0.0 , 0.0 , 0.0 , -0.002092 ,&
                       0.0 , 0.0 , 0.0 , -0.01343 ,&
                  0.9996 , 0.8247 , 0.0 , 7.453 ,&
                       0.0 , 0.0 , 15.82 , 5.400 ,&
                   56.57 , 112.5 , -27.37 , -1.623 ,&
                       0.0 , 0.0 , 0.0 , 0.1999 ,&
                 -0.1456 , -0.09054 , 0.0 , 0.0 ,&
                       0.0 , 0.0 , 0.0 , 0.4753 ,&
                  -12.10 , -2.249 , 0.0 , 0.0 ,&
                       0.0 , 0.01653 , 0.0 , -0.2047 ,&
                  0.1303 , 0.0 , 0.02288 , -0.02581 ,&
                       0.0 , 0.0 , 0.0 , 0.0 ,&
                       0.0 , 0.0 , 0.1062 , -0.9043 ,&
                  0.2950 , 0.8326 , -0.9992 , -0.3386 ,&
                       0.0 , -9.554 , 1.560 , 0.04556 ,&
                0.005508 , 0.0 , 0.0 , 0.04682 ,&
                -0.03590 , 0.07022 , 0.0 , -0.01924 ,&
               0.0004067 , -0.001333 , 0.0 , 0.01217 ,&
                       0.0 , 0.0 , 0.0 , 0.03944 ,&
                       0.0 , 0.0 , 0.0 , 0.006516 ,&
                       0.0 , 0.0 , 0.0 , -0.003571/), shape(coef_zl_b))
    if (Rib.LE.Ribc_b(1)) then
        ZL = sum(coef_zl_b(1,:)*SS1)
    else if (Rib.LE.Ribc_b(2)) then
         ZL = sum(coef_zl_b(2,:)*SS1)
    else if (Rib.LE.Ribc_b(3)) then
        ZL = sum(coef_zl_b(3,:)*SS1)
    else if (Rib.GT.Ribc_b(3)) then
        ZL = sum(coef_zl_b(4,:)*SS1)
    end if
    
else if (mz.LT.80 .and. mz.GE.10 .and. mzn.GE.100 .and.mzn.LT.10**7) then !region 3**************************
    coef_ribc_c=reshape((/ 0.3063 , 0.3555 , 0.5064 , 1.638 , &
                          -0.2849 , -0.3002 , -0.4282 , -1.743 , &
                          0.07886 , 0.07855 , 0.1229 , 0.5813 , &
                          0.03104 , 0.02617 , 0.02138 , 0.04471 , &
                         -0.01423 , -0.004769 , 0.0 , -0.01874 , &
                       -0.0005632 , -0.004012 , -0.004410 , 0.0 , &
                        3.684e-06 , -1.298e-05 , 0.0 , 0.0 , &
                       -2.926e-06 , 9.907e-06 , 0.0 , 0.0/),shape(coef_ribc_c))
    do i=1,4
        Ribc_c(i) = sum(coef_ribc_c(i,:)*XX1)
    end do
    coef_zl_c=reshape((/2.001 , 0.0 , -68.85 , -1.514 , 0.0 ,&
                      -0.7876 , 0.0 , 756.9 , 0.0 , 0.0 ,&
                            0.0 , 0.0 , -1100.0 , 0.0 , 0.0 ,&
                        60.42 , 368.9 , 0.0 , 19.63 , 0.0 ,&
                      -0.1401 , 3.514 , 0.0 , 0.5590 , 0.0 ,&
                      -0.1085 , -8.524 , -30.13 , 0.0 , 0.0 ,&
                       -2.065 , -18.05 , 86.99 , 0.0 , 0.0 ,&
                       -2.980 , -4.852 , 5.710 , -2.424 , 0.0 ,&
                      0.01334 , 0.08174 , 0.7274 , -0.002248 , 0.0,&
                      0.02130 , 0.5791 , -2.554 , 0.0 , 0.0 ,&
                       0.1963 , 0.1207 , -0.2169 , 0.1259 , 0.0 ,&
                   -0.0003704 , -0.007021 , 0.01587 , 0.0008267, 0.0002413 ,&
                    -0.002957 , 0.0 , 0.003912 , -0.004141 , 7.107e-05 ,&
                       -1.442 , 1.207 , 76.25 , -8.751 , 0.0 ,&
                        1.047 , -31.68 , -874.1 , 51.96 , 1.905 ,&
                            0.0 , 32.78 , 1636.0 , -76.51 , -1.761 ,&
                            0.0 , -25.65 , -1040.0 , 27.69 , 0.3658 ,&
                            0.0 , -2.096 , 4.942 , -1.349 , -0.05227 ,&
                            0.0 , 2.222 , -17.32 , 1.297 , 0.0 ,&
                       -1.121 , 0.3871 , 14.97 , -0.09621 , 0.0 ,&
                            0.0 , -0.004486 , -0.09096 , 0.0 , 0.0 ,&
                      0.02730 , -0.06669 , 0.2281 , 0.0 , 0.0 ,&
                            0.0 , 0.001086 , -0.002971 , 0.0002192 , 0.0 ,&
                       0.6868 , -0.07632 , -21.66 , 3.734 , 2.165,&
                            0.0 , 14.32 , 232.4 , -6.438 , 0.6139 ,&
                        3.820 , 2.353 , -224.1 , 6.284 , -0.1166 ,&
                     -0.01898 , 0.3396 , -1.724 , 0.2422 ,-0.07307 ,&
                      -0.1228 , -0.3281 , 3.144 , -0.2272 ,0.005656 ,&
                    0.0002845 , -0.0003600 , -0.0004477 , 0.0 , 0.0 ,&
                     -0.06543 , 0.0 , 1.875 , -0.4111 , -0.3134 ,&
                       0.1469 , -1.505 , -18.02 , 0.2556 , 0.0 ,&
                     0.001790 , -0.01529 , 0.1523 , -0.009961 ,0.008105/),shape(coef_zl_c))
    if (Rib.LE.Ribc_c(1)) then
        ZL = sum(coef_zl_c(1,:)*SS1)
    else if (Rib.LE.Ribc_c(2)) then
        ZL = sum(coef_zl_c(2,:)*SS1)
    else if (Rib.LE.Ribc_c(3)) then
        ZL = sum(coef_zl_c(3,:)*SS1)
    else if (Rib.LE.Ribc_c(4)) then
        ZL = sum(coef_zl_c(4,:)*SS1)
    else if (Rib.GT.Ribc_c(4)) then
        ZL = sum(coef_zl_c(5,:)*SS1)
    end if
else if (mz.LE.10**5 .and. mz.GE.80 .and. mzn.GE.100 .and.mzn.LT.10**7) then !region 4**************************
  coef_ribc_d=reshape((/0.09742 , 0.1768 , 0.3636  ,&
                              0.0 , 0.0 , 0.0  ,&
                       -0.01096 , -0.01434 , -0.0224  ,&
                        0.04544 , 0.03558 , 0.04607   ,&
                       -0.03299 , -0.02059 , -0.02506   ,&
                       0.006383 , 0.003327 , 0.004152   ,& 
                              0.0 , 0.0 , 0.0  ,&
                              0.0 , 0.0 , 0.0/), shape(coef_ribc_d))
    do i=1,3
        Ribc_d(i) = sum(coef_ribc_d(i,:)*XX1)
    end do  
  coef_zl_d=reshape((/0.0 , -3.528 , 0.0 , 0.0 ,&
                      0.0 , 0.0 , 0.0 , 0.0 ,&
                      0.0 , 0.0 , 0.0 , -8.306 ,&
                      0.0 , 0.0 , 0.0 , 1.212 ,&
                      0.0 , -0.2511 , -1.018 , 0.0 ,&
                      0.0 , 0.0 , 0.0 , 0.0 ,&
                 -6.267 , -10.06 , 0.0 , 0.0 ,&
                      0.0 , 0.0 , 0.0 , 0.0 ,&
                      0.0 , 0.0 , 0.0 , 0.0 ,&
                0.09808 , 0.1809 , 0.0 , 0.02790 ,&
                      0.0 , 0.0 , 0.0 , 0.0 ,&
                      0.0 , 0.0 , 6.740e-05 , 0.0006853 ,&
                      0.0 , 0.0 , 0.001341 , -0.0009314 ,&
                 0.5961 , 1.375 , -2.404 , 5.253 ,&
                      0.0 , 2.951 , 41.12 , 7.626 ,&
                  18.49 , 68.09 , -48.05 , -0.2889 ,&
                  34.53 , 0.0 , 24.94 , 0.06073 ,&
               -0.08450 , 0.0 , -0.06671 , -0.3959 ,&
                -0.5106 , -1.361 , 0.0 , -0.07098 ,&
                -0.3543 , 0.0 , -0.1319 , 0.003821 ,&
               0.004555 , 0.003711 , 0.006818 , 0.0 ,&
                      0.0 , 0.0 , 0.0 , 0.0 ,&
             -9.402e-05 , 0.0 , -0.0001788 , 0.0 ,&
                0.05628 , -0.02359 , 0.5172 , -0.5006 ,&
                 0.8075 , 0.3050 , -4.023 , -0.7376 ,&
                      0.0 , -3.765 , 2.074 , 0.0 ,&
                      0.0 , -0.001535 , 0.0 , 0.04853 ,&
                0.01631 , 0.07098 , 0.0 , 0.002956 ,&
             -3.800e-05 , -0.0002577 , 0.0 , 0.0 ,&
              -0.001890 , 0.0 , -0.01920 , 0.01968 ,&
               -0.03755 , 0.0 , 0.1250 , 0.02500 ,&
              5.177e-05 , 0.0 , 0.0 , -0.001897/),shape(coef_zl_d))
    if (Rib.LE.Ribc_d(1)) then
        ZL = sum(coef_zl_d(1,:)*SS1)
    else if (Rib.LE.Ribc_d(2)) then
        ZL = sum(coef_zl_d(2,:)*SS1)
    else if (Rib.LE.Ribc_d(3)) then
        ZL = sum(coef_zl_d(3,:)*SS1)
    else if (Rib.GT.Ribc_d(3)) then
        ZL = sum(coef_zl_d(4,:)*SS1)
    end if
else if (mz.LT.40 .and. mz.GE.10 .and. mzn.GE.10**7 .and.mzn.LT.10**11) then !region 5**************************
   coef_ribc_e=reshape((/0.0 , 0.0 , 0.0 , 0.0 , 0.4437 ,&
                         0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,&
                         0.0 , 0.08807 , 0.1219 , 0.1609 , 0.0 ,&
                   0.04825 , 0.05219 , 0.05830 , 0.07789 , 0.1349 ,&
                  -0.01677 , -0.01822 , -0.02373 , -0.04617 , -0.1388 ,&
                 -0.004762 , -0.01245 , -0.01224 , -0.007360 , 0.03347 ,&
                -0.0005212 , -0.0008500 , -0.001081 , -0.001399 , -0.001190 ,&
                 0.0002768 , 0.0007516 , 0.0009539 , 0.001238 , 0.001095/),shape(coef_ribc_e))
    do i=1,5
        Ribc_e(i) = sum(coef_ribc_e(i,:)*XX1)
    end do
  coef_zl_e=reshape((/0.0 ,  0.0 ,  -207.7  ,  -587.1  ,  0.0 ,  0.0 ,&
                      0.0 ,  77.11  ,  880.0 ,  2726.0 ,  7.886  ,  0.0 ,&
                -2.541  ,  -201.2  ,  -1550.0 ,  -3759.0 ,  -0.5889  ,  0.0 ,&
                 25.22  ,  386.1  ,  2201.0 ,  1605.0 ,  0.0 ,  0.0 ,&
              -0.03201  ,  -0.6831  ,  0.0 ,  -9.376  ,  -0.4057  ,  0.0 ,&
                0.1159  ,  0.0 ,  11.61  ,  -4.513  ,  0.0 ,  0.0 ,&
               -0.5745  ,  -7.571  ,  -96.51  ,  70.55  ,  -0.5218  ,  0.0  ,&
               -0.8502  ,  -8.978  ,  0.0 ,  -58.16  ,  0.0 ,  0.0 ,&
                0.00208 ,  0.07136  ,  0.5093  ,  0.1711  ,  0.01745  ,  0.0,&
             -0.001668  ,  0.0 ,  0.8873  ,  -0.9373  ,  -0.01349  ,  0.0 ,&
               0.03737  ,  0.3442  ,  0.2868  ,  1.132  ,  0.01468  ,  0.0,&
            -1.828e-05 ,  0.0 ,  -0.001909  ,  -0.006865  ,  0.0 ,  0.0 ,&
            -0.0003967  ,  -0.003421  ,  -0.004313  ,  -0.001126  ,  0.0 ,0.0 ,&
                0.4298  ,  0.0 ,  189.4  ,  286.9  ,  0.0 ,  0.0 ,&
              -0.03339  ,  -31.72  ,  -543.8  ,  -903.7  ,  0.0 ,  0.0 ,&
               0.05692  ,  2.558  ,  324.0 ,  407.6  ,  0.0 ,  0.0 ,&
                      0.0 ,  0.0 ,  -80.25  ,  260.2  ,  0.0 ,  0.08919  ,&
                -0.0233 ,  0.0 ,  -5.403  ,  0.0 ,  0.0 ,  0.0 ,&
                      0.0 ,  2.695  ,  14.95  ,  14.82  ,  0.2908  ,  0.0 ,&
               -0.3158  ,  -2.449  ,  -1.706  ,  -26.07  ,  0.1992  ,  0.0,&
                      0.0 ,  -0.05044  ,  -0.4221  ,  0.01062  , -0.003177  , 0.0 ,&
              0.007595  ,  0.05465  ,  0.164  ,  0.2099  ,  -0.00933 ,0.0 ,&
                      0.0 ,  -6.869e-05 ,  -0.00111  ,  0.0009863  ,  0.0 ,0.0 ,&
                      0.0 ,  0.3612  ,  -53.83  ,  -44.24  ,  0.7321  , 2.053 ,&
                      0.0 ,  0.0 ,  89.42  ,  98.98  ,  2.304  ,  0.2534  ,&
                 1.793  ,  18.63  ,  34.6 ,  22.67  ,  -2.456  , -0.2585  ,&
                0.00249 ,  0.1236  ,  2.704  ,  -0.01096  ,  -0.09448  , -0.0338 ,&
              -0.05666  ,  -0.837 ,  -4.573  ,  -1.67 ,  0.007636  , 0.004269 ,& 
                      0.0 ,  0.008316  ,  0.0718 ,  -0.01056  ,  0.002124 ,  0.0 ,&
                      0.0 ,  -0.06987  ,  4.95 ,  2.138  ,  0.0 ,  -0.3116 ,&
                   0.129,  0.8756  ,  -3.112  ,  -4.604  ,  0.0 ,  0.1241  ,&
                      0.0 ,  -0.01959  ,  -0.3287  ,  0.054  ,  0.0 , 0.0/),shape(coef_zl_e))
    if (Rib.LE.Ribc_e(1)) then
        ZL = sum(coef_zl_e(1,:)*SS1)
    else if (Rib.LE.Ribc_e(2)) then
        ZL = sum(coef_zl_e(2,:)*SS1)
    else if (Rib.LE.Ribc_e(3)) then
        ZL = sum(coef_zl_e(3,:)*SS1)
    else if (Rib.LE.Ribc_e(4)) then
        ZL = sum(coef_zl_e(4,:)*SS1)
    else if (Rib.LE.Ribc_e(5)) then
        ZL = sum(coef_zl_e(5,:)*SS1)
    else if (Rib.GT.Ribc_e(5)) then
        ZL = sum(coef_zl_e(6,:)*SS1)
    end if
else if (mz.LE.10**5 .and. mz.GE.40 .and. mzn.GE.10**7 .and.mzn.LT.10**11) then !region 6**************************
    coef_ribc_f=reshape((/0.0  , 0.1945   , 0.4288  ,&
                          0.0  , 0.0  , -0.1436   ,&
                          0.0  , 0.0  , 0.01635   ,&
                    0.05594  , 0.03347   , 0.03207   ,&
                   -0.03245  , -0.02116   , -0.01382   ,&
                   0.005037  , 0.002301   , 0.001571   ,&
                 -0.0003654  , 0.0  , 1.326e-05  ,&
                 0.0001135   , 8.92e-05  , -6.424e-06/),shape(coef_ribc_f))

    do i=1,3
        Ribc_f(i) = sum(coef_ribc_f(i,:)*XX1)
    end do
   coef_zl_f=reshape((/0.0 , 0.4383 , 0.0 , -6.744  ,&
                  -7.864 , 0.0 , -41.74 , 8.800  ,&
                       0.0 , 0.0 , 177.0 , -13.03  ,&
                       0.0 , 0.0 , -118.2 , 2.203  ,&
                -0.02699 , 0.0 , 0.0 , -0.1139  ,&
                  0.7414 , -4.810 , -4.006 , -0.06103  ,&
                  -1.114 , 5.094 , -0.5102 , 0.2406  ,&
                       0.0 , -1.159 , 0.0 , -0.04635  ,&
                       0.0 , 0.04547 , 0.0 , 0.01341  ,&
                       0.0 , 0.0 , 0.05670 , -0.002749  ,&
                       0.0 , -0.1233 , 0.1868 , 5.316e-06  ,&
                       0.0 , -0.0005595 , 0.002457 , -0.0001434 ,&
               0.0001281 , 0.002459 , -0.006455 , 0.0  ,&
                  0.2440 , 0.0 , 0.0 , 6.511  ,&
                   1.743 , 0.0 , 27.45 , 6.369  ,&
                   4.749 , 44.44 , -17.37 , -0.1750  ,&
                   11.28 , 0.0 , -7.740 , 0.03419  ,&
                       0.0 , 0.0 , 0.0 , -0.3147  ,&
                 -0.3093 , 0.0 , 0.0 , -0.06781  ,&
                 -0.2208 , -0.6068 , 0.01170 , -0.0002026,&
                       0.0 , -0.005459 , -0.01576 , 0.002444  ,&
                0.003674 , 0.0 , 0.02102 , 0.0002616  ,&
                       0.0 , 0.0 , -1.975e-05 , -5.149e-06  ,&
                 0.04168 , 0.0 , -0.1563 , -0.6219  ,&
                  0.4341 , 0.9983 , -2.085 , -0.5980  ,&
                  0.6518 , -2.874 , 0.3443 , 0.002868  ,&
               -0.002080 , -0.001520 , 0.03278 , 0.03359,&
                       0.0 , 0.01501 , -0.03250 , 0.003178  ,&
               2.895e-05 , 0.0003541 , 0.0005167 ,-0.0001423  ,&
                       0.0 , 0.006587 , 0.008163 , 0.02407  ,&
                -0.01307 , -0.04253 , 0.08540 , 0.01880,&
               1.425e-05 , -0.0003659 , -0.001602 ,-0.001167/),shape(coef_zl_f))
    if (Rib.LE.Ribc_f(1)) then
        ZL = sum(coef_zl_f(1,:)*SS1)
    else if (Rib.LE.Ribc_f(2)) then
        ZL = sum(coef_zl_f(2,:)*SS1)
    else if (Rib.LE.Ribc_f(3)) then
        ZL = sum(coef_zl_f(3,:)*SS1)
    else if (Rib.GT.Ribc_f(3)) then
        ZL = sum(coef_zl_f(4,:)*SS1)
    end if
else if (mz.LT.40 .and. mz.GE.10 .and. mzn.GE.10**11 .and.mzn.LE.1.07*10**13) then  !!region 7**************************
    coef_ribc_g=reshape((/0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0  ,&
                          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,&
                          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,&
                    0.03681 , 0.03655 , 0.03822 , 0.03840 , 0.05616 , 0.1472 ,&
                  -0.007664 , 0.0 , 0.0 , 0.0 , -0.02275 , -0.1144,&
                  -0.005619 , -0.009977 , -0.01036 ,-0.009243 , 0.0 , 0.02796 ,&
                 -0.0001211 , -0.0002691 , -0.0003658 ,-0.0003629 , -0.0005172 , -0.001218 ,&
                          0.0 , 0.0001057 , 0.0001769 , 0.0001471 , 0.0002261 , 0.0005835/), shape(coef_ribc_g))
    do i=1,6
        Ribc_g(i) = sum(coef_ribc_g(i,:)*XX1)
    end do
    coef_zl_g=reshape((/-1.412 , -4.502 , -104.2 , 542.4 , 178.4 , 0.0 , 0.0,&
                         6.658 , 40.44 , 136.3 , -1845.0 , 158.8 , 0.0, 0.0 ,&
                        -5.680 , 37.42 , 233.3 , 2157.0 , -480.9 , 0.0, 0.0 ,&
                         11.90 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,&
                        0.1285 , 0.3067 , 13.80 , -3.691 , -31.49, 0.0 , 0.0 ,&
                       -0.1110 , -5.444 , -37.21 , 3.330 , 47.56, 0.0 , 0.0 ,&
                       -0.2095 , 2.053 , 10.33 , -45.62 , -4.153, 0.0 , 0.0 ,&
                       -0.3181 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,&
                     -0.004693 , 0.05302 , -0.1157 , 0.1434 ,0.3998 , 0.0 , 0.0 ,&
                      0.004467 , 0.0 , 0.5542 , 0.4557 , -0.8692 ,0.0 , 0.0 ,&
                       0.01324 , -0.01586 , -0.2568 , 0.08936 ,0.2504 , 0.0 , 0.0 ,&
                     6.640e-05 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,&
                    -0.0002023 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,&
                        0.7122 , 1.663 , 16.56 , -263.7 , -37.94 ,0.0 , 0.0 ,&
                        -4.599 , -28.10 , 0.0 , 677.7 , -147.8 ,20.56 , 0.0 ,&
                         2.705 , -11.02 , -114.4 , -644.2 , 144.7 ,-13.42 , 0.0 ,&
                             0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 3.002 , 0.0 ,&
                      -0.04962 , 0.1172 , -3.238 , 4.440 , 9.904, -0.5254 , 0.06758 ,&
                       0.01147 , 1.979 , 7.578 , -3.037 , -7.914, 0.0 , 0.0 ,&
                       -0.1621 , -0.7285 , 0.0 , 10.93 , -2.224 , 0.0, 0.003671 ,&
                      0.001459 , -0.02930 , 0.0 , -0.08875 ,-0.1235 , 0.0 , 0.0 ,&
                      0.003514 , 0.01334 , -0.06568 , -0.1436 ,0.1631 , 0.0 , -0.0006967 ,&
                    -2.010e-05 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0002282 , 0.0,&
                      0.003692 , -0.4475 , 0.0 , 32.93 , 0.0 , 0.0 , 0.0,&
                         1.299 , 5.193 , -0.1495 , -56.58 , 23.51 , -2.349 , 0.6983 ,&
                        0.6516 , 5.593 , 18.12 , 53.14 , -2.645 ,0.6280 , -0.1455 ,&
                             0.0 , -0.009728 , 0.1670 , -0.9510 , -0.7278 , 0.2176 , 0.0 ,&
                      -0.03414 , -0.3375 , -0.6387 , 0.0 , -0.1801, 0.02067 , 0.0 ,&
                     2.840e-05 , 0.003470 , 0.004280 , 0.02119, 0.008599 , -0.005396 , -0.0004282 ,&
                     0.0006293 , 0.0 , 0.0 , 0.0 , 0.0 , -0.4148 , 0.0 ,&
                      -0.02559 , 0.0 , 0.0 , 0.0 , 0.0 , 0.02245 , 0.0 ,&
                             0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.01163 , 0.0/),shape(coef_zl_g))
    if (Rib.LE.Ribc_g(1)) then
        ZL = sum(coef_zl_g(1,:)*SS1)
    else if (Rib.LE.Ribc_g(2)) then
        ZL = sum(coef_zl_g(2,:)*SS1)
    else if (Rib.LE.Ribc_g(3)) then
        ZL = sum(coef_zl_g(3,:)*SS1)
    else if (Rib.LE.Ribc_g(4)) then
        ZL = sum(coef_zl_g(4,:)*SS1)
    else if (Rib.LE.Ribc_g(5)) then
        ZL = sum(coef_zl_g(5,:)*SS1)
    else if (Rib.LE.Ribc_g(6)) then
        ZL = sum(coef_zl_g(6,:)*SS1)
    else if (Rib.GT.Ribc_g(6)) then
        ZL = sum(coef_zl_g(7,:)*SS1)
    end if
    
else if (mz.LE.10**5 .and. mz.GE.40 .and. mzn.GE.10**11 .and.mzn.LE.1.07*10**13) then  !!region 8**************************
    coef_ribc_h=reshape((/0.0 , 0.0 , 0.5775  ,&
                          0.0 , 0.0 , -0.2236 ,&
                          0.0 , 0.0 , 0.03477 ,&
                    0.05139 , 0.04919 , 0.03805 ,&
                   -0.02991 , -0.01970 , -0.01617 ,&
                   0.004664 , 0.002011 , 0.001770 ,&
                 -0.0002135 , -0.0003325 , -2.191e-05 ,&
                  6.535e-05 , 7.974e-05 , 1.067e-05/),shape(coef_ribc_h))
    do i=1,3
        Ribc_h(i) = sum(coef_ribc_h(i,:)*XX1)
    end do
    coef_zl_h=reshape((/-3.13 , -49.55 , 0.0 , 0.0  ,&
                         5.26 , 97.14 , 0.0 , 0.0  ,&
                       -29.85 , 352.5 , 10.72 , 0.0  ,&
                        57.04 , -573.4 , 0.0 , 0.0  ,&
                       0.2176 , 2.052 , 0.0 , 0.0  ,&
                     -0.00898 , -21.41 , 0.0 , 0.0  ,&
                       -1.756 , 13.12 , 0.0 , 0.0  ,&
                       -1.663 , 20.82 , -1.354 , 0.0  ,&
                    -0.007271 , 0.1357 , -0.06227 , 0.0  ,&
                       0.0304 , 0.238 , 0.0 , -0.01477  ,&
                      0.05349 , -0.7316 , 0.08799 , -0.001292,&
                    8.978e-05 , -0.003367 , 0.002359 , 0.0  ,&
                   -0.0006252 , 0.006023 , -0.002387 ,0.0003921  ,&
                       0.9846 , 14.57 , -0.2492 , 0.0  ,&
                       -1.011 , 0.0 , 19.79 , 0.0  ,&
                        14.45 , 0.0 , -18.86 , -0.8522  ,&
                        4.433 , -54.39 , 9.463 , 0.1065  ,&
                     -0.05083 , -0.8911 , 0.0 , 0.0  ,&
                      -0.2604 , 1.478 , 0.0 , 0.374  ,&
                      -0.2977 , 2.13 , -0.3291 , 0.004036  ,&
                     0.001361 , -0.000936 , 0.0 , 0.002528  ,&
                      0.00375 , -0.04272 , 0.01369 , -0.006853,&
                   -1.464e-05 , 0.0001939 , -0.000241 ,-8.747e-05  ,&
                    -0.004659 , -1.165 , 0.0 , 0.0  ,&
                       0.6393 , 0.0 , -1.689 , -0.4307  ,&
                            0.0 , -3.616 , 1.036 , 0.01469  ,&
                            0.0 , 0.06747 , 0.00194 , 0.001642  ,&
                            0.0 , 0.01581 , -0.02897 , 0.0  ,&
                            0.0 , -0.0003126 , 0.0008316 , 0.0  ,&
                    0.0008014 , 0.03485 , 0.01694 , 0.0  ,&
                     -0.01934 , 0.0 , 0.06734 , 0.01348  ,&
                            0.0 , -0.001713 , -0.001447 ,0.0/),shape(coef_zl_h))
    if (Rib.LE.Ribc_h(1)) then
        ZL = sum(coef_zl_f(1,:)*SS1)
    else if (Rib.LE.Ribc_h(2)) then
        ZL = sum(coef_zl_f(2,:)*SS1)
    else if (Rib.LE.Ribc_h(3)) then
        ZL = sum(coef_zl_f(3,:)*SS1)
    else if (Rib.GT.Ribc_h(3)) then
        ZL = sum(coef_zl_f(4,:)*SS1)
    end if
end if
zol_stable= ZL*s11

END FUNCTION zol_stable

REAL FUNCTION zol_unstable(Rib,zz,z0m,z0h)

real :: Rib, zz, z0m, z0h, L0m, L0h, &
        kvon, Rprt, zs_z0m, zs, YFIT, &
        x0, x11, x21, x22, x23, x31, x32, x33
    
    kvon=0.4
    Rprt=1.0
    zs_z0m = 16.7
    zs=zs_z0m*z0m
    L0m = log(zz/z0m)
    L0h = log(zz/z0h)  
    x0  = 1.0
    x11 = (1-1/(1.0-Rib))
    x21 = 1.0/(log(zz/z0m))
    x22 = x21**2
    x23 = x21**3
    x31 =1/(log(zz/z0h))
    x32 = x31**2
    x33 = x31**3
        
    if ( (Rib.GT.-2 .and. Rib.LE.-0) .and. zz/z0m.LE.80 .and. z0m/z0h.LE.10 )then !region 1
       YFIT = -116.93   *       x23*x31 + &
               126.9    *       x22*x32 + &
             -115.09    *       x21*x32 + &
              111.64    *       x22*x31 + &
              -3.572    *       x11*x32 + &  
              3.8304    *       x11*x22 + &
              -4.5297   *       x11*x21*x31 + &
              25.776    *       x32 + &
              0.45351   *       x22 + &
              -26.473   *       x21*x31 + & 
              6.6361    *       x11*x31 + & 
              -3.3234   *       x11*x21 + & 
              -0.65868  *       x11 + & 
              1.0393    *       x0

    else if ( (Rib.GT.-2 .and. Rib.LE.-0) .and. zz/z0m.LE.80 .and. z0m/z0h.GT.10 ) then !region 2
       YFIT = 14.904    *       x11*x32 + &  
              6.8022    *       x32 + &
              2.6097    *       x22 + &
              -6.1598   *       x21*x31 + & 
              -2.5799   *       x11*x31 + & 
              -0.98398  *       x11*x21 + & 
              -0.19337  *       x11 + & 
              0.87712   *       x0
    
    else if ( (Rib.GT.-2 .and. Rib.LE.-0) .and. zz/z0m.GT.80 .and.z0m/z0h.LE.10 ) then !region 3
       YFIT = 60.525    *       x21*x32 + &
             -128.39    *       x22*x31 + &
              33.947    *       x22 + &
              5.4086    *       x11*x31 + & 
               -4.91    *       x11*x21 + & 
               1.516    *       x31 + &
              -6.877    *       x21 + &
            -0.21653    *       x11 + & 
              1.2529    *       x0

   else if ( (Rib.GT.-2 .and. Rib.LE.-0) .and. zz/z0m.GT.80 .and. z0m/z0h.GT.10 ) then !region 4
        YFIT= 4.5577    *       x32 + &
              5.4274    *       x22 + &
              -1.5934   *       x11*x21 + & 
              -0.46985  *       x31 + &
              -1.4933   *       x21 + &
              -0.184    *       x11 + & 
              1.0526    *       x0

    else if ( Rib.LE.-2 .and. zz/z0m.LE.80 .and. z0m/z0h.LE.10 ) then !region 5
        YFIT= 133.3       *       x22*x32 + &
             12.548       *       x11*x22*x31 + &
             -34.644      *   x23 + &
             -121.31      *   x21*x32 + &
             -4.1189      *   x11*x32 + &  
               26.6       *       x32 + &
             33.905       *       x22 + &
             3.2532       *   x11*x31 + & 
             -3.4678      *   x11*x21 + & 
             -8.6705      *   x21 + &
             1.3343       *     x0

    else if ( Rib.LE.-2 .and. zz/z0m.LE.80 .and. z0m/z0h.GT.10 ) then !region 6
        YFIT= -21.221   *       x21*x32 + &
              9.7024    *       x22*x31 + &
              5.4758    *       x11*x21*x31 + &
              23.793    *       x32 + &
              2.8282    *       x22 + &
             -12.674    *       x21*x31 + & 
             -1.3558    *       x11*x21 + & 
             -1.7794    *       x31 + &
            -0.42907    *       x11 + & 
              1.0925    *       x0

    else if ( Rib.LE.-2 .and. zz/z0m.GT.80 .and. z0m/z0h.LE.10 ) then !region 7
        YFIT=-343.68    *       x23 + &
             -637.21    *       x21*x32 + &
              935.47    *       x22*x31 + &
              152.13    *       x32 + &
              85.433    *       x22 + &
              -212.9    *       x21*x31 + & 
              1.6941    *       x11*x31 + & 
             -5.0047    *       x21 + &
            -0.50742    *       x11 + & 
              1.4164    *       x0

   else if ( Rib.LE.-2 .and. zz/z0m.GT.80 .and. z0m/z0h.GT.10 ) then !region 8
        YFIT=-35.397    *       x21*x32 + &
               32.746   *       x32 + &
               10.744   *       x22 + &
               -19.461  *       x21*x31 + & 
               -3.3425  *       x21 + &
               -0.53463 *       x11 + & 
               1.3388   *       x0
   end if
       
    zol_unstable = YFIT*L0m*(L0m/L0h)*Rib

END FUNCTION zol_unstable

!-------------------------------------------------------------------
!---- add pbl modules so they can be optimized in pbl code
!-------------------------------------------------------------------
  SUBROUTINE  mym_condensation (kts,kte,  &
    &            dx, dz,                  &
    &            thl, qw,                 &
    &            p,exner,                 &
    &            tsq, qsq, cov,           &
    &            Sh, el, bl_mynn_cloudpdf,&
    &            qc_bl1D, cldfra_bl1D,    &
    &            PBLH1,HFX1,              &
    &            Vt, Vq, th, sgm)

!-------------------------------------------------------------------

    INTEGER, INTENT(IN)   :: kts,kte, bl_mynn_cloudpdf
    REAL, INTENT(IN)      :: dx,PBLH1,HFX1
    REAL, DIMENSION(kts:kte), INTENT(IN) :: dz
    REAL, DIMENSION(kts:kte), INTENT(IN) :: p,exner, thl, qw, &
         &tsq, qsq, cov, th

    REAL, DIMENSION(kts:kte), INTENT(INOUT) :: vt,vq,sgm

    REAL, DIMENSION(kts:kte) :: qmq,alp,a,bet,b,ql,q1,cld,RH
    REAL, DIMENSION(kts:kte), INTENT(OUT) :: qc_bl1D,cldfra_bl1D
    DOUBLE PRECISION :: t3sq, r3sq, c3sq

    REAL :: qsl,esat,qsat,tlk,qsat_tl,dqsl,cld0,q1k,eq1,qll,&
         &q2p,pt,rac,qt,t,xl,rsl,cpm,cdhdz,Fng,qww,alpha,beta,bb,ls_min,ls,wt
    INTEGER :: i,j,k

    REAL :: erf

    !JOE: NEW VARIABLES FOR ALTERNATE SIGMA
    REAL::dth,dtl,dqw,dzk
    REAL, DIMENSION(kts:kte), INTENT(IN) :: Sh,el

    !JOE: variables for BL clouds
    REAL::zagl,cld9,damp,edown,RHcrit,RHmean,RHsum,RHnum,Hshcu,PBLH2,ql_limit
    REAL, PARAMETER :: Hfac = 3.0     !cloud depth factor for HFX (m^3/W)
    REAL, PARAMETER :: HFXmin = 50.0  !min W/m^2 for BL clouds
    REAL            :: RH_00L, RH_00O, phi_dz, lfac
    REAL, PARAMETER :: cdz = 2.0
    REAL, PARAMETER :: mdz = 1.5

    !JAYMES:  variables for tropopause-height estimation
    REAL            :: theta1, theta2, ht1, ht2
    INTEGER         :: k_tropo

    REAL, PARAMETER :: rr2=0.7071068, rrp=0.3989423

    k_tropo=5

    zagl = 0.

    SELECT CASE(bl_mynn_cloudpdf)

      CASE (0) ! ORIGINAL MYNN PARTIAL-CONDENSATION SCHEME

        DO k = kts,kte-1
           t  = th(k)*exner(k)

           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           qsl=ep_2*esat/(p(k)-ep_3*esat)
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*ev/( rd*t**2 )
           !RH (0 to 1.0)
           RH(k)=MAX(MIN(1.0,qw(k)/MAX(1.E-8,qsl)),0.001)

           alp(k) = 1.0/( 1.0+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           !NOTE: negative bl_mynn_cloudpdf will zero-out the stratus subgrid clouds
           !      at the end of this subroutine.
           !Sommeria and Deardorff (1977) scheme, as implemented
           !in Nakanishi and Niino (2009), Appendix B
           t3sq = MAX( tsq(k), 0.0 )
           r3sq = MAX( qsq(k), 0.0 )
           c3sq =      cov(k)
           c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )
           r3sq = r3sq +bet(k)**2*t3sq -2.0*bet(k)*c3sq
           !DEFICIT/EXCESS WATER CONTENT
           qmq(k) = qw(k) -qsl
           !ORIGINAL STANDARD DEVIATION: limit e-6 produces ~10% more BL clouds
           !than e-10
           sgm(k) = SQRT( MAX( r3sq, 1.0d-10 ))
           !NORMALIZED DEPARTURE FROM SATURATION
           q1(k)   = qmq(k) / sgm(k)
           !CLOUD FRACTION. rr2 = 1/SQRT(2) = 0.707
           cld(k) = 0.5*( 1.0+erf( q1(k)*rr2 ) )

        END DO

      CASE (1, -1) !ALTERNATIVE FORM (Nakanishi & Niino 2004 BLM, eq. B6, and
                       !Kuwano-Yoshida et al. 2010 QJRMS, eq. 7):
        DO k = kts,kte-1
           t  = th(k)*exner(k)
           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           qsl=ep_2*esat/(p(k)-ep_3*esat)
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*ev/( rd*t**2 )
           !RH (0 to 1.0)
           RH(k)=MAX(MIN(1.0,qw(k)/MAX(1.E-8,qsl)),0.001)

           alp(k) = 1.0/( 1.0+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           if (k .eq. kts) then
             dzk = 0.5*dz(k)
           else
             dzk = 0.5*( dz(k) + dz(k-1) )
           end if
           dth = 0.5*(thl(k+1)+thl(k)) - 0.5*(thl(k)+thl(MAX(k-1,kts)))
           dqw = 0.5*(qw(k+1) + qw(k)) - 0.5*(qw(k) + qw(MAX(k-1,kts)))
           sgm(k) = SQRT( MAX( (alp(k)**2 * MAX(el(k)**2,0.1) * &
                             b2 * MAX(Sh(k),0.03))/4. * &
                      (dqw/dzk - bet(k)*(dth/dzk ))**2 , 1.0e-10) )
           qmq(k) = qw(k) -qsl
           q1(k)   = qmq(k) / sgm(k)
           cld(k) = 0.5*( 1.0+erf( q1(k)*rr2 ) )
        END DO

      CASE (2, -2)
          !Diagnostic statistical scheme of Chaboureau and Bechtold (2002), JAS
          !JAYMES- this added 27 Apr 2015
        DO k = kts,kte-1
           t  = th(k)*exner(k)
           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           qsl=ep_2*esat/(p(k)-ep_3*esat)
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*ev/( rd*t**2 )
           !RH (0 to 1.0)
           RH(k)=MAX(MIN(1.0,qw(k)/MAX(1.E-8,qsl)),0.001)

           alp(k) = 1.0/( 1.0+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           xl = xl_blend(t)                    ! obtain latent heat

           tlk = thl(k)*(p(k)/p1000mb)**rcp    ! recover liquid temp (tl) from thl
           qsat_tl = qsat_blend(tlk,p(k))      ! get saturation water vapor mixing ratio
                                               !   at tl and p

           rsl = xl*qsat_tl / (r_v*tlk**2)     ! slope of C-C curve at t = tl
                                               ! CB02, Eqn. 4

           cpm = cp + qw(k)*cpv                ! CB02, sec. 2, para. 1

           a(k) = 1./(1. + xl*rsl/cpm)         ! CB02 variable "a"

           qmq(k) = a(k) * (qw(k) - qsat_tl) ! saturation deficit/excess;
                                               !   the numerator of Q1

           b(k) = a(k)*rsl                     ! CB02 variable "b"

           dtl =    0.5*(thl(k+1)*(p(k+1)/p1000mb)**rcp + tlk) &
               & - 0.5*(tlk + thl(MAX(k-1,kts))*(p(MAX(k-1,kts))/p1000mb)**rcp)

           dqw = 0.5*(qw(k+1) + qw(k)) - 0.5*(qw(k) + qw(MAX(k-1,kts)))

           if (k .eq. kts) then
             dzk = 0.5*dz(k)
           else
             dzk = 0.5*( dz(k) + dz(k-1) )
           end if

           cdhdz = dtl/dzk + (g/cpm)*(1.+qw(k))  ! expression below Eq. 9
                                                 ! in CB02
           zagl = zagl + dz(k)
           ls_min = MIN(MAX(zagl,25.),300.) ! Let this be the minimum possible length scale:
                                        !   25 m < ls_min(=zagl) < 300 m
           lfac=MIN(4.25+dx/4000.,6.)   ! A dx-dependent multiplier for the master length scale:
                                        !   lfac(750 m) = 4.4
                                        !   lfac(3 km)  = 5.0
                                        !   lfac(13 km) = 6.0

           ls = MAX(MIN(lfac*el(k),900.),ls_min)  ! Bounded:  ls_min < ls < 900 m
                   ! Note: CB02 use 900 m as a constant free-atmosphere length scale.
                   ! Above 300 m AGL, ls_min remains 300 m.  For dx = 3 km, the
                   ! MYNN master length scale (el) must exceed 60 m before ls
                   ! becomes responsive to el, otherwise ls = ls_min = 300 m.

           sgm(k) = MAX(1.e-10, 0.225*ls*SQRT(MAX(0., & ! Eq. 9 in CB02:
                   & (a(k)*dqw/dzk)**2              & ! < 1st term in brackets,
                   & -2*a(k)*b(k)*cdhdz*dqw/dzk     & ! < 2nd term,
                   & +b(k)**2 * cdhdz**2)))           ! < 3rd term
                   ! CB02 use a multiplier of 0.2, but 0.225 is chosen
                   ! based on tests

           q1(k) = qmq(k) / sgm(k)  ! Q1, the normalized saturation

           cld(k) = MAX(0., MIN(1., 0.5+0.36*ATAN(1.55*q1(k)))) ! Eq. 7 in CB02

         END DO
    END SELECT

    zagl = 0.
    RHsum=0.
    RHnum=0.
    RHmean=0.1 !initialize with small value for small PBLH cases
    damp =0
    PBLH2=MAX(10.,PBLH1)

    SELECT CASE(bl_mynn_cloudpdf)

      CASE (-1 : 1) ! ORIGINAL MYNN PARTIAL-CONDENSATION SCHEME
                    ! OR KUWANO ET AL.
        DO k = kts,kte-1
           t    = th(k)*exner(k)
           q1k  = q1(k)
           zagl = zagl + dz(k)
           !q1=0.
           !cld(k)=0.

           !COMPUTE MEAN RH IN PBL (NOT PRESSURE WEIGHTED). 
           IF (zagl < PBLH2 .AND. PBLH2 > 400.) THEN
              RHsum=RHsum+RH(k)
              RHnum=RHnum+1.0
              RHmean=RHsum/RHnum
           ENDIF
           RHcrit = 1. - 0.35*(1.0 - (MAX(250.- MAX(HFX1,HFXmin),0.0)/200.)**2)
           if (HFX1 > HFXmin) then
              cld9=MIN(MAX(0., (rh(k)-RHcrit)/(1.1-RHcrit)), 1.)**2
           else
              cld9=0.0
           endif

           edown=PBLH2*.1
           !Vary BL cloud depth (Hshcu) by mean RH in PBL and HFX
           !(somewhat following results from Zhang and Klein (2013, JAS))
           Hshcu=200. + (RHmean+0.5)**1.5*MAX(HFX1,0.)*Hfac
           if (zagl < PBLH2-edown) then
              damp=MIN(1.0,exp(-ABS(((PBLH2-edown)-zagl)/edown)))
           elseif(zagl >= PBLH2-edown .AND. zagl < PBLH2+Hshcu)then
              damp=1.
           elseif (zagl >= PBLH2+Hshcu)then
              damp=MIN(1.0,exp(-ABS((zagl-(PBLH2+Hshcu))/500.)))
           endif
           cldfra_bl1D(k)=cld9*damp
           !cldfra_bl1D(k)=cld(k) ! JAYMES: use this form to retain the Sommeria-Deardorff value

           !use alternate cloud fraction to estimate qc for use in BL clouds-radiation
           eq1  = rrp*EXP( -0.5*q1k*q1k )
           qll  = MAX( cldfra_bl1D(k)*q1k + eq1, 0.0 )
           !ESTIMATED LIQUID WATER CONTENT (UNNORMALIZED)
           ql (k) = alp(k)*sgm(k)*qll
           if(cldfra_bl1D(k)>0.01 .and. ql(k)<1.E-6)ql(k)=1.E-6
           qc_bl1D(k)=ql(k)*damp
           !now recompute estimated lwc for PBL scheme's use
           !qll IS THE NORMALIZED LIQUID WATER CONTENT (Sommeria and
           !Deardorff (1977, eq 29a). rrp = 1/(sqrt(2*pi)) = 0.3989
           eq1  = rrp*EXP( -0.5*q1k*q1k )
           qll  = MAX( cld(k)*q1k + eq1, 0.0 )
           !ESTIMATED LIQUID WATER CONTENT (UNNORMALIZED)
           ql (k) = alp(k)*sgm(k)*qll

           q2p = xlvcp/exner(k)
           pt = thl(k) +q2p*ql(k) ! potential temp

           !qt is a THETA-V CONVERSION FOR TOTAL WATER (i.e., THETA-V = qt*THETA)
           qt   = 1.0 +p608*qw(k) -(1.+p608)*ql(k)
           rac  = alp(k)*( cld(k)-qll*eq1 )*( q2p*qt-(1.+p608)*pt )

           !BUOYANCY FACTORS: wherever vt and vq are used, there is a
           !"+1" and "+tv0", respectively, so these are subtracted out here.
           !vt is unitless and vq has units of K.
           vt(k) =      qt-1.0 -rac*bet(k)
           vq(k) = p608*pt-tv0 +rac

           !To avoid FPE in radiation driver, when these two quantities are multiplied by eachother,
           ! add limit to qc_bl and cldfra_bl:
           IF (QC_BL1D(k) < 1E-6 .AND. ABS(CLDFRA_BL1D(k)) > 0.01) QC_BL1D(k)= 1E-6
           IF (CLDFRA_BL1D(k) < 1E-2)THEN
              CLDFRA_BL1D(k)=0.
              QC_BL1D(k)=0.
           ENDIF

        END DO
      CASE ( 2, -2)
        ! JAYMES- this option added 8 May 2015
        ! The cloud water formulations are taken from CB02, Eq. 8.
        ! "fng" represents the non-Gaussian contribution to the liquid
        ! water flux; these formulations are from Cuijpers and Bechtold
        ! (1995), Eq. 7.  CB95 also draws from Bechtold et al. 1995,
        ! hereafter BCMT95
        DO k = kts,kte-1
           t    = th(k)*exner(k)
           q1k  = q1(k)
           zagl = zagl + dz(k)
           IF (q1k < 0.) THEN
              ql (k) = sgm(k)*EXP(1.2*q1k-1)
           ELSE IF (q1k > 2.) THEN
              ql (k) = sgm(k)*q1k
           ELSE
              ql (k) = sgm(k)*(EXP(-1.) + 0.66*q1k + 0.086*q1k**2)
           ENDIF

           !Next, adjust our initial estimates of cldfra and ql based
           !on tropopause-height and PBLH considerations
           !JAYMES:  added 4 Nov 2016
           if ((cld(k) .gt. 0.) .or. (ql(k) .gt. 0.))  then
              if (k .le. k_tropo) then
                 !At and below tropopause: impose an upper limit on ql; assume that
                 !a maximum of 0.5 percent supersaturation in water vapor can be
                 !available for cloud production
                 ql_limit = 0.005 * qsat_blend( th(k)*exner(k), p(k) )
                 ql(k) = MIN( ql(k), ql_limit )
              else
                 !Above tropopause:  eliminate subgrid clouds from CB scheme
                 cld(k) = 0.
                 ql(k) = 0.
              endif
           endif

           !Buoyancy-flux-related calculations follow...
           ! "Fng" represents the non-Gaussian transport factor
           ! (non-dimensional) from from Bechtold et al. 1995
           ! (hereafter BCMT95), section 3(c).  Their suggested
           ! forms for Fng (from their Eq. 20) are:
           ! For purposes of the buoyancy flux in stratus, we will use Fng = 1
           Fng = 1.

           xl    = xl_blend(t)
           bb = b(k)*t/th(k) ! bb is "b" in BCMT95.  Their "b" differs from
                             ! "b" in CB02 (i.e., b(k) above) by a factor
                             ! of T/theta.  Strictly, b(k) above is formulated in
                             ! terms of sat. mixing ratio, but bb in BCMT95 is
                             ! cast in terms of sat. specific humidity.  The
                             ! conversion is neglected here.
           qww   = 1.+0.61*qw(k)
           alpha = 0.61*th(k)
           beta  = (th(k)/t)*(xl/cp) - 1.61*th(k)

           vt(k) = qww   - cld(k)*beta*bb*Fng   - 1.
           vq(k) = alpha + cld(k)*beta*a(k)*Fng - tv0
           ! vt and vq correspond to beta-theta and beta-q, respectively,
           ! in NN09, Eq. B8.  They also correspond to the bracketed
           ! expressions in BCMT95, Eq. 15, since (s*ql/sigma^2) = cldfra*Fng
           ! The "-1" and "-tv0" terms are included for consistency with
           ! the legacy vt and vq formulations (above).

           ! increase the cloud fraction estimate below PBLH+1km
           if (zagl .lt. PBLH2+1000.) cld(k) = MIN( 1., 1.8*cld(k) )
           ! return a cloud condensate and cloud fraction for icloud_bl option:
           cldfra_bl1D(k) = cld(k)
           qc_bl1D(k) = ql(k)

           !To avoid FPE in radiation driver, when these two quantities are multiplied by eachother,
           ! add limit to qc_bl and cldfra_bl:
           IF (QC_BL1D(k) < 1E-6 .AND. ABS(CLDFRA_BL1D(k)) > 0.01) QC_BL1D(k)= 1E-6
           IF (CLDFRA_BL1D(k) < 1E-2)THEN
              CLDFRA_BL1D(k)=0.
              QC_BL1D(k)=0.
           ENDIF

         END DO

      END SELECT !end cloudPDF option

      !FOR TESTING PURPOSES ONLY, ISOLATE ON THE MASS-CLOUDS.
      IF (bl_mynn_cloudpdf .LT. 0) THEN
         DO k = kts,kte-1
              cldfra_bl1D(k) = 0.0
              qc_bl1D(k) = 0.0
         END DO
      ENDIF

      cld(kte) = cld(kte-1)
      ql(kte) = ql(kte-1)
      vt(kte) = vt(kte-1)
      vq(kte) = vq(kte-1)
      qc_bl1D(kte)=0.
      cldfra_bl1D(kte)=0.

    RETURN

  END SUBROUTINE mym_condensation

! ==================================================================


END MODULE module_sf_mynn

