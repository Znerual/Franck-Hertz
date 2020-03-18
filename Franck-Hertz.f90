!  FranckHertz.f90 
!
!  FUNCTIONS:
!  FranckHertz - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: FranckHertz
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
!Typischer Druck 20 mbar -> 20mbar = n * kb T -> n =20*10^-3 (mBar)* 10^5 (in Pascal) / (1,38 *10^-23 (Boltzmann) * 300 (Raumtemp.)) = 5 * 10^23
!Stoﬂquerschnitt (Hartkugelmodell) etwa 2.8*10^-19 m^2
!emfp = 7.1 * 10^-2 = 7.1cm
    
    !Erwartungswert <x^N> = Integral(x * p(x) dx) bei bekanntem p(x)
    !               <x^N> = 1/N*Sum((x - <x>)^N)
    
    !       hier: p(x) = 1/ emfp * exp(-x/emfp)
    
    
    
    program FranckHertz
    !for using the random function in Intel Compiler
    use ifport
    implicit none
    
    
    !Functions
    double precision :: p
    
    ! Variables
    integer, parameter :: N_trajektorien = 4000 !Anzahl der Trajektorien
    double precision, parameter :: distanz_kg = 2d-2 !Abstand von Kathode zu Gitter
    double precision, parameter :: emfp = 7.1d-2 !emfp (Elastic mean free pathlength) = 1 / n * sigma (n= Dichte, Sigma = Stoﬂquerschnitt)
    double precision, parameter :: E_anregung = 4.9d0 !eV
    double precision, parameter :: Ug = 1d0, Ub_max = 2d1, Ub_min = 0.5d0 !Volt
    double precision, parameter :: U_step = 1d-2
    double precision, parameter :: c = 1.686d5 / distanz_kg
    double precision :: scaleX, scaleY, x, fx !For creating the rectangle for the monte carlo simulation and the random guesses
    double precision :: E, Ub , Uh, Ia !Energie des Elektrons (eV), Beschleunigungsspannung, Gegenspannung, Katheodenspannung, Anodenstrom
    integer :: i,n=0, t
    double precision :: dummy, max = 0
    
   
    
    scaleX = distanz_kg
    scaleY = 1d0 /  emfp !Because 1/Lambda is the max of the function
    
    print*, "Scale X: " ,scaleX, " ScaleY: ", scaleY
    open(1, file="output.dat", action="write")
    
    !Initialisiere den Zufallszahlengenerator
    call System_Clock(count=t)
    dummy = rand(t)
    
    do Ub = Ub_min, Ub_max, U_step
        do i = 0, N_trajektorien
            x = rand() * scaleX
            fx = rand() * scaleY
            E = Ub/distanz_kg * x
            if (fx > p(x, emfp) .and. E > E_anregung) E = E - E_anregung !Freie Wegl‰nge ist zu klein, genug Energie zur Anregung vorhanden, Teilchen stˆﬂt
            E = E + (scaleX - x) * Ub/distanz_kg - Ug!Energie (nach Stoﬂ) weiter erhˆht und durch Gegenspannung erniedrigt
            if (E >= 0) n = n + 1!Teilchen ist an der Anode angekommen
            !print*, "Teilchen simuliert bei x: " , x ," f(x): ", fx, " mit Energie " , (Ub-Ug), " hat jetzt Energie ", E
        end do
        !Ia = dble(n) / dble(N_trajektorien) * scaleX * scaleY
        Ia = c * sqrt(Ub) * n/N_trajektorien
        n = 0
        write(1,*) Ub, Ia 
    end do
    
    close(1)
    
    
    end program FranckHertz

    double precision function p(x,lambda)
        implicit none
        double precision x, lambda
        p = 1d0 / lambda * exp(-x/lambda)
    end function
