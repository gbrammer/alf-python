module driver

  use alf_vars; use alf_utils
  USE nr, ONLY : locate
  
  implicit none
  save

  !f2py intent(hide) pset
  type(PARAMS) :: pset

  integer :: is_setup=0

contains

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
    
    write(*,*) "dlstep, nl_fit:", dlstep, nl
    
  end subroutine
  
  subroutine get_nspec(ns)
    implicit none
    integer, intent(out) :: ns
    ns = nl
  end subroutine
  
  subroutine get_npar(np)
    implicit none
    integer, intent(out) :: np
    np = npar
  end subroutine
    
  subroutine get_grid_lam(nl, lam)
    implicit none
    integer, intent(in) :: nl
    double precision, dimension(nl), intent(out) :: lam
    lam = sspgrid%lam    
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
    posarr(34) = pos%logtrans
    posarr(35) = pos%logemline_h
    posarr(36) = pos%logemline_oiii
    posarr(37) = pos%logemline_sii
    posarr(38) = pos%logemline_ni
    posarr(39) = pos%logemline_nii
    posarr(40) = pos%jitter
    posarr(41) = pos%imf3
    posarr(42) = pos%logsky
    posarr(43) = pos%imf4
    posarr(44) = pos%h3
    posarr(45) = pos%h4
    
  end subroutine
  
  subroutine velbroad_spec(lambda, spec, sigma, minl, maxl)
    implicit none
    double precision, dimension(:), intent(in) :: lambda
    double precision, dimension(:), intent(inout) :: spec
    double precision, intent(in) :: sigma,minl,maxl
    
    CALL VELBROAD(lambda, spec, sigma, minl, maxl)
    
  end subroutine
  
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
    pos%logtrans  = posarr(34)
    pos%logemline_h    = posarr(35)
    pos%logemline_oiii = posarr(36)
    pos%logemline_sii  = posarr(37)
    pos%logemline_ni   = posarr(38)
    pos%logemline_nii  = posarr(39)
    pos%jitter = posarr(40)
    pos%imf3   = posarr(41)
    pos%logsky = posarr(42)
    pos%imf4   = posarr(43)
    pos%h3     = posarr(44)
    pos%h4     = posarr(45)
    
    !WRITE(*,*) 'sigma, age, FeH, MgH: ', pos%sigma, pos%logage, pos%feh, pos%mgh
    !WRITE(*,*) 'Fit_type, powell_fitting: ', fit_type, powell_fitting
    
    CALL GETMODEL(pos,mspec,mw=1)
    
  end subroutine
  
end module
