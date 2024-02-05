module driver

  use alf_vars; use alf_utils
  USE nr, ONLY : locate
  
  implicit none
  save

  !f2py intent(hide) pset
  type(PARAMS) :: pset

  integer :: is_setup=0

contains
  
  ! Call setup and update dlstep/inlam
  subroutine setup_alf
        
    implicit none
    integer :: i
    
    call setup()
    is_setup = 1
    
    !Need to do these here, which would normally be done within alf.f90
    
    !define the log wavelength grid used in velbroad.f90
    dlstep = (LOG(sspgrid%lam(nl))-LOG(sspgrid%lam(1)))/nl
    DO i=1,nl
       lnlam(i) = i*dlstep+LOG(sspgrid%lam(1))
    ENDDO
    
    ! write(*,*) "dlstep, nl_fit:", dlstep, nl
    
  end subroutine

  ! Set fit type
  subroutine set_fit_type(ft)
    implicit none
    integer, intent(in) :: ft
    fit_type = ft
    WRITE(*,*) 'Fit_type:', fit_type
  end subroutine
  
  ! Set wave limits
  subroutine set_wave_limits(il1, il2)
    implicit none
    double precision, intent(in) :: il1, il2
    
    l1 = il1
    l2 = il2
    nlint = 1
    
    WRITE(*,*) 'Wavelength limits l1, l2, nlint:', l1(1), l2(nlint), nlint
  end subroutine
  
  ! Set emline mask
  subroutine set_maskem(em)
    implicit none
    integer, intent(in) :: em
    maskem = em
    WRITE(*,*) 'maskem:', maskem
  end subroutine
  
  ! Set IMF type
  subroutine set_imf(mw, type)
    implicit none
    integer, intent(in) :: mw, type
    mwimf = mw
    imf_type = type 
    WRITE(*,*) 'mwimf:', mwimf
    WRITE(*,*) 'imf_type:', imf_type
  end subroutine

  ! Set fit_hermite flag
  subroutine set_fit_hermite(herm)
    implicit none
    integer, intent(in) :: herm
    fit_hermite = herm
    WRITE(*,*) 'fit_hermite:', fit_hermite
  end subroutine

  ! Set fit_two_ages flag
  subroutine set_fit_two_ages(two_ages)
    implicit none
    integer, intent(in) :: two_ages
    fit_two_ages = two_ages
    WRITE(*,*) 'fit_two_ages:', fit_two_ages
  end subroutine

  ! Set velbroad_simple flag
  subroutine set_velbroad_simple(velb)
    implicit none
    integer, intent(in) :: velb
    velbroad_simple = velb
    WRITE(*,*) 'velbroad_simple:', velbroad_simple
  end subroutine
    
  ! Get dimensions of the SSP spectra
  subroutine get_nspec(ns)
    implicit none
    integer, intent(out) :: ns
    ns = nl
  end subroutine
  
  ! Get number of parameters
  subroutine get_npar(np)
    implicit none
    integer, intent(out) :: np
    np = npar
  end subroutine
  
  ! Get number of filters (R, i, K)
  subroutine get_nfil(nf)
    implicit none
    integer, intent(out) :: nf
    nf = nfil
  end subroutine
  
  ! Get the SSP wavelength grid
  subroutine get_grid_lam(nl, lam)
    implicit none
    integer, intent(in) :: nl
    double precision, dimension(nl), intent(out) :: lam
    lam = sspgrid%lam    
  end subroutine
  
  ! Get number of ssp ages
  subroutine get_nage(nages)
    implicit none
    integer, intent(out) :: nages
    nages = nage
  end subroutine
  
  ! Get the SSP age grid
  subroutine get_ssp_logagegrid(nage, logagegrid)
    implicit none
    integer, intent(in) :: nage
    double precision, dimension(nage), intent(out) :: logagegrid
    logagegrid = sspgrid%logagegrid    
  end subroutine
  
  subroutine get_default_parameters(np, posarr)
    implicit none
    integer, intent(in) :: np
    double precision, dimension(np), intent(out) :: posarr
  
    TYPE(PARAMS) :: pos
  
    posarr(1) = pos%velz
    posarr(2) = pos%sigma
    posarr(3) = pos%logage
    posarr(4) = pos%zh
    !end of the super-simple and Powell-mode parameters

    posarr(5) = pos%feh
    posarr(6) = pos%ah
    posarr(7) = pos%ch
    posarr(8) = pos%nh
    posarr(9) = pos%nah
    posarr(10) = pos%mgh
    posarr(11) = pos%sih
    posarr(12) = pos%kh
    posarr(13) = pos%cah
    posarr(14) = pos%tih
    !end of the simple model parameters

    posarr(15) = pos%vh
    posarr(16) = pos%crh
    posarr(17) = pos%mnh
    posarr(18) = pos%coh
    posarr(19) = pos%nih
    posarr(20) = pos%cuh
    posarr(21) = pos%srh
    posarr(22) = pos%bah
    posarr(23) = pos%euh

    posarr(24) = pos%teff
    posarr(25) = pos%imf1
    posarr(26) = pos%imf2
    posarr(27) = pos%logfy
    posarr(28) = pos%sigma2
    posarr(29) = pos%velz2
    posarr(30) = pos%logm7g
    posarr(31) = pos%hotteff
    posarr(32) = pos%loghot
    posarr(33) = pos%fy_logage
    posarr(34) = pos%logemline_h
    posarr(35) = pos%logemline_oii
    posarr(36) = pos%logemline_oiii
    posarr(37) = pos%logemline_sii
    posarr(38) = pos%logemline_ni
    posarr(39) = pos%logemline_nii
    posarr(40) = pos%logtrans
    posarr(41) = pos%jitter
    posarr(42) = pos%logsky
    posarr(43) = pos%imf3
    posarr(44) = pos%imf4
    posarr(45) = pos%h3
    posarr(46) = pos%h4
    
  end subroutine
  
  subroutine velbroad_spec(lambda, spec, sigma, minl, maxl)
    implicit none
    double precision, dimension(:), intent(in) :: lambda
    double precision, dimension(:), intent(inout) :: spec
    double precision, intent(in) :: sigma,minl,maxl
    
    CALL VELBROAD(lambda, spec, sigma, minl, maxl)
    
  end subroutine
  
  subroutine velbroad_spec_hermite(lambda, spec, sigma, minl, maxl, h34)
    implicit none
    double precision, dimension(:), intent(in) :: lambda
    double precision, dimension(:), intent(inout) :: spec
    double precision, intent(in) :: sigma,minl,maxl
    double precision, dimension(:), intent(in), optional :: h34
    
    CALL VELBROAD(lambda, spec, sigma, minl, maxl, h34)
    
  end subroutine
  
  ! Get model spectrum for array of parameter values `posarr` and flexible IMF params
  ! Call GETMODEL(pos, mspec) (without mw=1)
  subroutine get_spec(nl, np, posarr, mspec)
    implicit none
    integer, intent(in) :: nl
    integer, intent(in) :: np
    double precision, dimension(np), intent(in) :: posarr
    double precision, dimension(nl), intent(out) :: mspec
    
    TYPE(PARAMS) :: pos

    !WRITE(*,*) 'ParamsA: ', posarr
    
    !arr->str

    pos%velz   = posarr(1)
    pos%sigma  = posarr(2)
    pos%logage = posarr(3)
    pos%zh     = posarr(4)
    !end of the super-simple and Powell-mode parameters

    pos%feh    = posarr(5)
    pos%ah     = posarr(6)
    pos%ch     = posarr(7)
    pos%nh     = posarr(8)
    pos%nah    = posarr(9)
    pos%mgh    = posarr(10)
    pos%sih    = posarr(11)
    pos%kh     = posarr(12)
    pos%cah    = posarr(13)
    pos%tih    = posarr(14)
    !end of the simple model parameters

    pos%vh     = posarr(15)
    pos%crh    = posarr(16)
    pos%mnh    = posarr(17)
    pos%coh    = posarr(18)
    pos%nih    = posarr(19)
    pos%cuh    = posarr(20)
    pos%srh    = posarr(21)
    pos%bah    = posarr(22)
    pos%euh    = posarr(23)

    pos%teff      = posarr(24)
    pos%imf1      = posarr(25)
    pos%imf2      = posarr(26)
    pos%logfy     = posarr(27)
    pos%sigma2    = posarr(28)
    pos%velz2     = posarr(29)
    pos%logm7g    = posarr(30)
    pos%hotteff   = posarr(31)
    pos%loghot    = posarr(32)
    pos%fy_logage = posarr(33)
    pos%logemline_h    = posarr(34)
    pos%logemline_oii  = posarr(35)
    pos%logemline_oiii = posarr(36)
    pos%logemline_sii  = posarr(37)
    pos%logemline_ni   = posarr(38)
    pos%logemline_nii  = posarr(39)
    pos%logtrans  = posarr(40)
    pos%jitter = posarr(41)
    pos%logsky = posarr(42)
    pos%imf3   = posarr(43)
    pos%imf4   = posarr(44)
    pos%h3     = posarr(45)
    pos%h4     = posarr(46)
    
    !WRITE(*,*) 'sigma, age, FeH, MgH: ', pos%sigma, pos%logage, pos%feh, pos%mgh
    !WRITE(*,*) 'Fit_type, powell_fitting: ', fit_type, powell_fitting
    
    CALL GETMODEL(pos, mspec)
    
  end subroutine

  ! Get model spectrum for array of parameter values `posarr` and MW (Kroupa) IMF
  ! Call with GETMODEL(pos, mspec, mw=1)
  subroutine get_spec_mw(nl, np, posarr, mspec)
    implicit none
    integer, intent(in) :: nl
    integer, intent(in) :: np
    double precision, dimension(np), intent(in) :: posarr
    double precision, dimension(nl), intent(out) :: mspec
    
    TYPE(PARAMS) :: pos

    !WRITE(*,*) 'ParamsA: ', posarr
    
    !arr->str

    pos%velz   = posarr(1)
    pos%sigma  = posarr(2)
    pos%logage = posarr(3)
    pos%zh     = posarr(4)
    !end of the super-simple and Powell-mode parameters

    pos%feh    = posarr(5)
    pos%ah     = posarr(6)
    pos%ch     = posarr(7)
    pos%nh     = posarr(8)
    pos%nah    = posarr(9)
    pos%mgh    = posarr(10)
    pos%sih    = posarr(11)
    pos%kh     = posarr(12)
    pos%cah    = posarr(13)
    pos%tih    = posarr(14)
    !end of the simple model parameters

    pos%vh     = posarr(15)
    pos%crh    = posarr(16)
    pos%mnh    = posarr(17)
    pos%coh    = posarr(18)
    pos%nih    = posarr(19)
    pos%cuh    = posarr(20)
    pos%srh    = posarr(21)
    pos%bah    = posarr(22)
    pos%euh    = posarr(23)

    pos%teff      = posarr(24)
    pos%imf1      = posarr(25)
    pos%imf2      = posarr(26)
    pos%logfy     = posarr(27)
    pos%sigma2    = posarr(28)
    pos%velz2     = posarr(29)
    pos%logm7g    = posarr(30)
    pos%hotteff   = posarr(31)
    pos%loghot    = posarr(32)
    pos%fy_logage = posarr(33)
    pos%logemline_h    = posarr(34)
    pos%logemline_oii  = posarr(35)
    pos%logemline_oiii = posarr(36)
    pos%logemline_sii  = posarr(37)
    pos%logemline_ni   = posarr(38)
    pos%logemline_nii  = posarr(39)
    pos%logtrans  = posarr(40)
    pos%jitter = posarr(41)
    pos%logsky = posarr(42)
    pos%imf3   = posarr(43)
    pos%imf4   = posarr(44)
    pos%h3     = posarr(45)
    pos%h4     = posarr(46)
    
    !WRITE(*,*) 'sigma, age, FeH, MgH: ', pos%sigma, pos%logage, pos%feh, pos%mgh
    !WRITE(*,*) 'Fit_type, powell_fitting: ', fit_type, powell_fitting
    
    CALL GETMODEL(pos, mspec, mw=1)
    
  end subroutine
  
  subroutine get_m2l(nl, np, nf, posarr, m2l)
    implicit none
    integer, intent(in) :: nl, np, nf
    !integer, intent(in) :: np
    double precision, dimension(np), intent(in) :: posarr
    double precision, dimension(nf), intent(out) :: m2l

    double precision, dimension(nl) :: lam, mspec    
    double precision :: msto

    TYPE(PARAMS) :: pos

    CALL get_grid_lam(nl, lam)
    
    !WRITE(*,*) 'ParamsA: ', posarr
    
    !arr->str

    pos%velz   = posarr(1)
    pos%sigma  = posarr(2)
    pos%logage = posarr(3)
    pos%zh     = posarr(4)
    !end of the super-simple and Powell-mode parameters

    pos%feh    = posarr(5)
    pos%ah     = posarr(6)
    pos%ch     = posarr(7)
    pos%nh     = posarr(8)
    pos%nah    = posarr(9)
    pos%mgh    = posarr(10)
    pos%sih    = posarr(11)
    pos%kh     = posarr(12)
    pos%cah    = posarr(13)
    pos%tih    = posarr(14)
    !end of the simple model parameters

    pos%vh     = posarr(15)
    pos%crh    = posarr(16)
    pos%mnh    = posarr(17)
    pos%coh    = posarr(18)
    pos%nih    = posarr(19)
    pos%cuh    = posarr(20)
    pos%srh    = posarr(21)
    pos%bah    = posarr(22)
    pos%euh    = posarr(23)

    pos%teff      = posarr(24)
    pos%imf1      = posarr(25)
    pos%imf2      = posarr(26)
    pos%logfy     = posarr(27)
    pos%sigma2    = posarr(28)
    pos%velz2     = posarr(29)
    pos%logm7g    = posarr(30)
    pos%hotteff   = posarr(31)
    pos%loghot    = posarr(32)
    pos%fy_logage = posarr(33)
    pos%logemline_h    = posarr(34)
    pos%logemline_oii  = posarr(35)
    pos%logemline_oiii = posarr(36)
    pos%logemline_sii  = posarr(37)
    pos%logemline_ni   = posarr(38)
    pos%logemline_nii  = posarr(39)
    pos%logtrans  = posarr(40)
    pos%jitter = posarr(41)
    pos%logsky = posarr(42)
    pos%imf3   = posarr(43)
    pos%imf4   = posarr(44)
    pos%h3     = posarr(45)
    pos%h4     = posarr(46)
    
    !turn off various parameters for computing M/L
    pos%logemline_h    = -8.0
    pos%logemline_oiii = -8.0
    pos%logemline_nii  = -8.0
    pos%logemline_sii  = -8.0
    pos%logemline_ni   = -8.0
    pos%logtrans       = -8.0
    
    ! Is this safe in order to be a bit faster? 
    pos%sigma       = 0.
    
    !WRITE(*,*) 'sigma, age, FeH, MgH: ', pos%sigma, pos%logage, pos%feh, pos%mgh
    !WRITE(*,*) 'Fit_type, powell_fitting: ', fit_type, powell_fitting
    
    CALL GETMODEL(pos,mspec,mw=1)
    
    !compute the main sequence turn-off mass vs. t and Z
    msto = MAX(MIN(10**(msto_t0+msto_t1*pos%logage) * &
         (msto_z0+msto_z1*pos%zh+msto_z2*pos%zh**2),3.0),0.75)
    
    CALL GETM2L(lam,mspec,pos,m2l,mw=1) !compute M/L_MW
    
  end subroutine
  
end module
