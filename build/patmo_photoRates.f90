module patmo_photoRates
contains

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

    !O2 -> O + O
    krate(:,269) = integrateXsec(1, tau(:,:))

    !O3 -> O + O2
    krate(:,270) = integrateXsec(2, tau(:,:))

    !O3 -> O(1D) + O2
    krate(:,271) = integrateXsec(3, tau(:,:))

    !N2 -> N + N
    krate(:,272) = integrateXsec(4, tau(:,:))

    !OH -> O + H
    krate(:,273) = integrateXsec(5, tau(:,:))

    !OH -> O(1D) + H
    krate(:,274) = integrateXsec(6, tau(:,:))

    !HO2 -> OH + O
    krate(:,275) = integrateXsec(7, tau(:,:))

    !H2O -> OH + H
    krate(:,276) = integrateXsec(8, tau(:,:))

    !H2O -> H2 + O
    krate(:,277) = integrateXsec(9, tau(:,:))

    !H2 -> H + H
    krate(:,278) = integrateXsec(10, tau(:,:))

    !N2O -> N2 + O(1D)
    krate(:,279) = integrateXsec(11, tau(:,:))

    !NO -> N + O
    krate(:,280) = integrateXsec(12, tau(:,:))

    !NO2 -> NO + O
    krate(:,281) = integrateXsec(13, tau(:,:))

    !NO3 -> NO + O2
    krate(:,282) = integrateXsec(14, tau(:,:))

    !NO3 -> O + NO2
    krate(:,283) = integrateXsec(15, tau(:,:))

    !N2O5 -> NO2 + NO3
    krate(:,284) = integrateXsec(16, tau(:,:))

    !N2O5 -> O + NO + NO3
    krate(:,285) = integrateXsec(17, tau(:,:))

    !HNO3 -> OH + NO2
    krate(:,286) = integrateXsec(18, tau(:,:))

    !HNO3 -> H + NO3
    krate(:,287) = integrateXsec(19, tau(:,:))

    !CH4 -> CH3 + H
    krate(:,288) = integrateXsec(20, tau(:,:))

    !CH3OOH -> CH3O + OH
    krate(:,289) = integrateXsec(21, tau(:,:))

    !CH2O -> H + HCO
    krate(:,290) = integrateXsec(22, tau(:,:))

    !CH2O -> H2 + CO
    krate(:,291) = integrateXsec(23, tau(:,:))

    !HCO -> H + CO
    krate(:,292) = integrateXsec(24, tau(:,:))

    !CO2 -> CO + O
    krate(:,293) = integrateXsec(25, tau(:,:))

    !H2O2 -> OH + OH
    krate(:,294) = integrateXsec(26, tau(:,:))

    !H2O2 -> H + HO2
    krate(:,295) = integrateXsec(27, tau(:,:))

    !COS -> CO + S
    krate(:,296) = integrateXsec(28, tau(:,:))

    !SO -> S + O
    krate(:,297) = integrateXsec(29, tau(:,:))

    !CS2 -> CS + S
    krate(:,298) = integrateXsec(30, tau(:,:))

    !H2S -> SH + H
    krate(:,299) = integrateXsec(31, tau(:,:))

    !SO2 -> SO + O
    krate(:,300) = integrateXsec(32, tau(:,:))

    !SO3 -> SO2 + O
    krate(:,301) = integrateXsec(33, tau(:,:))

    !H2SO4 -> SO2 + OH + OH
    krate(:,302) = integrateXsec(34, tau(:,:))

    !CH3OH -> CH3 + OH
    krate(:,303) = integrateXsec(35, tau(:,:))

    !CH3OH -> CH3O + H
    krate(:,304) = integrateXsec(36, tau(:,:))

    !S2O2 -> SO + SO
    krate(:,305) = integrateXsec(37, tau(:,:))

    !S2O -> SO + S
    krate(:,306) = integrateXsec(38, tau(:,:))

    !N2H4 -> H + N2H3
    krate(:,307) = integrateXsec(39, tau(:,:))

    !NH3 -> H + NH2
    krate(:,308) = integrateXsec(40, tau(:,:))

    !NH3 -> H2 + NH
    krate(:,309) = integrateXsec(41, tau(:,:))

    !HCN -> CN + H
    krate(:,310) = integrateXsec(42, tau(:,:))

    !C2H4 -> C2H2 + H2
    krate(:,311) = integrateXsec(43, tau(:,:))

    !CHOCHO -> CO + CO + H2
    krate(:,312) = integrateXsec(44, tau(:,:))

    !CHOCHO -> CH2O + CO
    krate(:,313) = integrateXsec(45, tau(:,:))

    !HCOOH -> CO2 + H2
    krate(:,314) = integrateXsec(46, tau(:,:))

    !HCOOH -> HCO + OH
    krate(:,315) = integrateXsec(47, tau(:,:))

    !CH3CHO -> CH4 + CO
    krate(:,316) = integrateXsec(48, tau(:,:))

    !CH3CHO -> CH3 + HCO
    krate(:,317) = integrateXsec(49, tau(:,:))

    !C2H6 -> CH3 + CH3
    krate(:,318) = integrateXsec(50, tau(:,:))

    !C2H6 -> C2H5 + H
    krate(:,319) = integrateXsec(51, tau(:,:))

    !C2H6 -> C2H4 + H2
    krate(:,320) = integrateXsec(52, tau(:,:))

  end subroutine computePhotoRates

  !*************
  function integrateXsec(index,tau)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    integer,intent(in)::index
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)
    real*8::integrateXsec(cellsNumber), dE, mu, coef
    integer::j

    ! !loop on cells (stride photobins)
    ! do j=1,cellsNumber
    !    integrateXsec(j) = sum(xsecAll(:,index)*photoFlux(:) &
        !         /energyMid(:)*energySpan(:)*exp(-tau(:,j))) / planck_eV
    ! end do

    !dE = (wavelengMax-wavelengMin)/photoBinsNumber (nm)
    dE = 0.1
    !mu =cosine(zenith_angle)
    mu = 0.500000
    !Parameter of incident solar flux = coef (Cronin, 2014- DOI:10.1175/JAS-D-13-0392.1)
    coef = 0.500000
    !loop on cells (stride photobins)
    do j=1,cellsNumber
      integrateXsec(j) = sum(xsecAll(:,index)*coef*photoFlux(:)*exp(-tau(:,j)/mu)*dE)
    end do

  end function integrateXsec

end module patmo_photoRates
