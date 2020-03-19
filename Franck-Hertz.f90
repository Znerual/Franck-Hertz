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
    integer, parameter :: N_simulation = 20 !Anzahl der Simulierten Stoﬂmˆglichkeiten
    double precision, parameter :: U_step = 1d-2 !^-2!
    double precision, parameter :: distanz_kg = 2d-1 !Abstand von Kathode zu Gitter
    double precision, parameter :: emfp = 7.1d-2 !emfp (Elastic mean free pathlength) = 1 / n * sigma (n= Dichte, Sigma = Stoﬂquerschnitt)
    double precision, parameter :: E_anregung = 4.9d0 !eV
    double precision, parameter :: Ug = 1d0, Ub_max = 8d1, Ub_min = 0.5d0 !Volt, eV 
    double precision, parameter :: c = 3.7d-9 / distanz_kg**2 !4/9 e0 sqrt(2 q / me) *S /d^2    !S mit 4cm^2 angenommen (Anodenfl‰che)
    double precision, parameter :: simulations_breite = distanz_kg / N_simulation
    double precision, parameter :: W = 4d0, Temp = 5d4 , kB = 1.38064852d-23, beta = Temp * kb, Uh = 1d3,U_sat= 1d2, q_elem= 1.6d-19 !Zur bestimmung der Exponentialverteilung der Startenergie und Austrittselektronen, !Temp = in etwa Austrittstem
    double precision, parameter :: Strom_faktor = 10d0
    double precision :: scaleX, scaleY, x,delta_x, fx !For creating the rectangle for the monte carlo simulation and the random guesses
    double precision :: E, Ub , Ia !Energie des Elektrons (eV), Beschleunigungsspannung, Gegenspannung, Katheodenspannung, Anodenstrom
    integer :: i,n=0,j, t
    integer :: N_trajektorien = 1000 !Anzahl an austretenden Elektronen
    double precision :: dummy
    double precision :: N_Verteilung !Variablen f¸r Startenergieverteilung
   
    
    scaleX = distanz_kg
    scaleY = 1d0 /  emfp !Because 1/Lambda is the max of the function
    
    print*, "Scale X: " ,scaleX, " ScaleY: ", scaleY
    open(1, file="output.dat", action="write")
    
    !Initialisiere den Zufallszahlengenerator
    call System_Clock(count=t)
    dummy = rand(t)
    
    do Ub = Ub_min, Ub_max, U_step
       
        do i = 0, N_trajektorien
            !Generieren der Exponentialverteilten Startenergie 
            N_Verteilung = rand()
            !E= (1d0 - exp(-1d0/N_Verteilung)) * ((kb*Temp)/q_elem - W)
            E = 0
            !print* ,"E: ", E
            do j = 0,N_simulation
                delta_x = rand() * simulations_breite
                x = delta_x + j * simulations_breite 
                fx = rand() * scaleY
                E = E + Ub/distanz_kg * delta_x
                if (fx > p(x, emfp) .and. E > E_anregung) E = E - E_anregung !Freie Wegl‰nge ist zu klein, genug Energie zur Anregung vorhanden, Teilchen stˆﬂt
            end do
            E = E - Ug
            if (E >= 0) n = n + 1!Teilchen ist an der Anode angekommen
            !print*, "Teilchen simuliert bei x: " , x ," f(x): ", fx, " mit Energie " , (Ub-Ug), " hat jetzt Energie ", E
        end do
        if (Ub > U_sat) then
            Ia = c * sqrt(Ub)**(3)  * n /N_trajektorien
            print*,"‹ber der S‰ttigungsspannung!"
        else 
            Ia = c * sqrt(Ub)** 3 * n /N_trajektorien
        end if
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
