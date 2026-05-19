module patmo_photo
contains

  !**************
  subroutine loadPhotoMetric(fname)
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    character::aout
    real*8::rout(4)
    integer::ios,i

    !open metric file (energy: mid/eV, span/eV, left/eV, right/eV)
    open(22,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: problem loading "//trim(fname)
      stop
    end if

    read(22,'(a)') aout

    !loop on bins to read (already interpolated by python)
    do i=1,photoBinsNumber
      read(22,*,iostat=ios) rout(:)
      if(ios/=0) then
        print *,"ERROR: problem while reading "//trim(fname)
        stop
      end if
      !load energy mid, span, left, right (eV)
      energyMid(i) = rout(1)
      energySpan(i) = rout(2)
      energyLeft(i) = rout(3)
      energyRight(i) = rout(4)
    end do

    close(22)

  end subroutine loadPhotoMetric

  !***************
  subroutine loadAllPhotoXsecs()
    use patmo_parameters
    implicit none

    !O2 -> O + O
    call loadPhotoXsec("xsecs/O2__O_O.dat",1)
    !O3 -> O + O2
    call loadPhotoXsec("xsecs/O3__O_O2.dat",2)
    !O3 -> O(1D) + O2
    call loadPhotoXsec("xsecs/O3__O(1D)_O2.dat",3)
    !N2 -> N + N
    call loadPhotoXsec("xsecs/N2__N_N.dat",4)
    !OH -> O + H
    call loadPhotoXsec("xsecs/OH__O_H.dat",5)
    !OH -> O(1D) + H
    call loadPhotoXsec("xsecs/OH__O(1D)_H.dat",6)
    !HO2 -> OH + O
    call loadPhotoXsec("xsecs/HO2__OH_O.dat",7)
    !H2O -> OH + H
    call loadPhotoXsec("xsecs/H2O__OH_H.dat",8)
    !H2O -> H2 + O
    call loadPhotoXsec("xsecs/H2O__H2_O.dat",9)
    !H2 -> H + H
    call loadPhotoXsec("xsecs/H2__H_H.dat",10)
    !N2O -> N2 + O(1D)
    call loadPhotoXsec("xsecs/N2O__N2_O(1D).dat",11)
    !NO -> N + O
    call loadPhotoXsec("xsecs/NO__N_O.dat",12)
    !NO2 -> NO + O
    call loadPhotoXsec("xsecs/NO2__NO_O.dat",13)
    !NO3 -> NO + O2
    call loadPhotoXsec("xsecs/NO3__NO_O2.dat",14)
    !NO3 -> O + NO2
    call loadPhotoXsec("xsecs/NO3__O_NO2.dat",15)
    !N2O5 -> NO2 + NO3
    call loadPhotoXsec("xsecs/N2O5__NO2_NO3.dat",16)
    !N2O5 -> O + NO + NO3
    call loadPhotoXsec("xsecs/N2O5__O_NO_NO3.dat",17)
    !HNO3 -> OH + NO2
    call loadPhotoXsec("xsecs/HNO3__OH_NO2.dat",18)
    !HNO3 -> H + NO3
    call loadPhotoXsec("xsecs/HNO3__H_NO3.dat",19)
    !CH4 -> CH3 + H
    call loadPhotoXsec("xsecs/CH4__CH3_H.dat",20)
    !CH3OOH -> CH3O + OH
    call loadPhotoXsec("xsecs/CH3OOH__CH3O_OH.dat",21)
    !CH2O -> H + HCO
    call loadPhotoXsec("xsecs/CH2O__H_HCO.dat",22)
    !CH2O -> H2 + CO
    call loadPhotoXsec("xsecs/CH2O__H2_CO.dat",23)
    !HCO -> H + CO
    call loadPhotoXsec("xsecs/HCO__H_CO.dat",24)
    !CO2 -> CO + O
    call loadPhotoXsec("xsecs/CO2__CO_O.dat",25)
    !H2O2 -> OH + OH
    call loadPhotoXsec("xsecs/H2O2__OH_OH.dat",26)
    !H2O2 -> H + HO2
    call loadPhotoXsec("xsecs/H2O2__H_HO2.dat",27)
    !COS -> CO + S
    call loadPhotoXsec("xsecs/COS__CO_S.dat",28)
    !SO -> S + O
    call loadPhotoXsec("xsecs/SO__S_O.dat",29)
    !CS2 -> CS + S
    call loadPhotoXsec("xsecs/CS2__CS_S.dat",30)
    !H2S -> SH + H
    call loadPhotoXsec("xsecs/H2S__SH_H.dat",31)
    !SO2 -> SO + O
    call loadPhotoXsec("xsecs/SO2__SO_O.dat",32)
    !SO3 -> SO2 + O
    call loadPhotoXsec("xsecs/SO3__SO2_O.dat",33)
    !H2SO4 -> SO2 + OH + OH
    call loadPhotoXsec("xsecs/H2SO4__SO2_OH_OH.dat",34)
    !CH3OH -> CH3 + OH
    call loadPhotoXsec("xsecs/CH3OH__CH3_OH.dat",35)
    !CH3OH -> CH3O + H
    call loadPhotoXsec("xsecs/CH3OH__CH3O_H.dat",36)
    !S2O2 -> SO + SO
    call loadPhotoXsec("xsecs/S2O2__SO_SO.dat",37)
    !S2O -> SO + S
    call loadPhotoXsec("xsecs/S2O__SO_S.dat",38)
    !N2H4 -> H + N2H3
    call loadPhotoXsec("xsecs/N2H4__H_N2H3.dat",39)
    !NH3 -> H + NH2
    call loadPhotoXsec("xsecs/NH3__H_NH2.dat",40)
    !NH3 -> H2 + NH
    call loadPhotoXsec("xsecs/NH3__H2_NH.dat",41)
    !HCN -> CN + H
    call loadPhotoXsec("xsecs/HCN__CN_H.dat",42)
    !C2H4 -> C2H2 + H2
    call loadPhotoXsec("xsecs/C2H4__C2H2_H2.dat",43)
    !CHOCHO -> CO + CO + H2
    call loadPhotoXsec("xsecs/CHOCHO__CO_CO_H2.dat",44)
    !CHOCHO -> CH2O + CO
    call loadPhotoXsec("xsecs/CHOCHO__CH2O_CO.dat",45)
    !HCOOH -> CO2 + H2
    call loadPhotoXsec("xsecs/HCOOH__CO2_H2.dat",46)
    !HCOOH -> HCO + OH
    call loadPhotoXsec("xsecs/HCOOH__HCO_OH.dat",47)
    !CH3CHO -> CH4 + CO
    call loadPhotoXsec("xsecs/CH3CHO__CH4_CO.dat",48)
    !CH3CHO -> CH3 + HCO
    call loadPhotoXsec("xsecs/CH3CHO__CH3_HCO.dat",49)
    !C2H6 -> CH3 + CH3
    call loadPhotoXsec("xsecs/C2H6__CH3_CH3.dat",50)
    !C2H6 -> C2H5 + H
    call loadPhotoXsec("xsecs/C2H6__C2H5_H.dat",51)
    !C2H6 -> C2H4 + H2
    call loadPhotoXsec("xsecs/C2H6__C2H4_H2.dat",52)

  end subroutine loadAllPhotoXsecs

  !***************
  subroutine loadPhotoXsec(fname,index)
    use patmo_commons
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    character(len=100)::aout
    integer,intent(in)::index
    integer::ios,i
    real*8::rout(2)

    !open xsec file
    open(22,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: problem loading "//trim(fname)
      stop
    end if

    !skip header (one line)
    read(22,'(a)') aout

    !loop on bins to read (already interpolated by python)
    do i=1,photoBinsNumber
      read(22,*,iostat=ios) rout(:)
      if(ios/=0) then
        print *,"ERROR: problem while reading "//trim(fname)
        stop
      end if
      !load xsecs for all cells (cm2)
      xsecAll(i,index) = rout(2)
    end do

    close(22)
  end subroutine loadPhotoXsec

  !**************************
  function fluxBB(x,Tbb)
    use patmo_constants
    implicit none
    real*8,intent(in)::x,Tbb
    real*8::fluxBB,xexp

    !exponent
    xexp = x/kboltzmann_eV/Tbb

    !default value
    fluxBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
      fluxBB = 2d0*x**3/planck_eV**2/clight**2 &
          / (exp(xexp)-1d0)
    end if

  end function fluxBB

end module patmo_photo
