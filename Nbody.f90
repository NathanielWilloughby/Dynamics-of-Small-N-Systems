! Main program - automates multiple N-body simulations
PROGRAM n_body_master
    IMPLICIT NONE
    INTEGER :: nsims, i
    DOUBLE PRECISION :: munit, runit, vunit, tunit

    ! Some adjustable values
    nsims = 1 ! Number of simulations to run
    munit = 1.99e30 ! Equal to one N-body mass unit (which is currently set to one solar mass)
    runit = 3.086e16 ! Equal to one N-body length unit (which is currently set to one parsec)

    ! Get velocity unit in km/s and time unit in years
    CALL n_body_conversion(munit, runit, tunit, vunit)

    ! Run N-body simulations by repeatedly calling n_body_simulation(n, ndims, initialConditions, tunit, vunit)
    ! n - number of objects
    ! ndims - number of dimensions (probably should always be 3)
    ! initialConditions - an integer used to determine which setup subroutine to use
    ! tunit - N-body time unit in years
    ! vunit - N-body velocity unit in km/s
    ! Currently available initial conditions are:
    ! 1 - Plummer sphere
    ! 2 - Figure of eight test (assumes n=3)
    ! If any integer that is not listed above is used, the program will default to the Plummer sphere setup
    DO i=1,nsims
        CALL n_body_simulation(10, 3, 1, tunit, vunit)
    END DO

END PROGRAM n_body_master

! One full simulation
SUBROUTINE n_body_simulation(n, ndims, initialConditions, tunit, vunit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, ndims, initialConditions
    DOUBLE PRECISION, INTENT(IN) :: tunit, vunit
    INTEGER :: i, j, k, count, encounters(1:n, 1:n)
    DOUBLE PRECISION :: t, tcount, tLim, tWrite, dtOpt(1:n), dt, dt2, dt3, dt4, dt5, recipdt2, recipdt3, recip6, recip24, recip120, E, E0, Ek, Ep, Eerr, EerrLim
    DOUBLE PRECISION :: dist, dist2, dist3, dist5, dr(1:ndims), dv(1:ndims), vmag2, amag2, vdotr, mTot, m(1:n), comPos(1:ndims), comVel(1:ndims)
    DOUBLE PRECISION :: pos(1:ndims, 1:n), vel(1:ndims, 1:n), acc(1:ndims, 1:n), jrk(1:ndims, 1:n), snp(1:ndims, 1:n), crk(1:ndims, 1:n)
    DOUBLE PRECISION :: accPred(1:ndims, 1:n), jrkPred(1:ndims, 1:n), eta, origDists(1:n), vmags(1:n), closest(1:2, 1:n), closestDist, closeThreshold

    ! Variable definitions:
    ! n - number of objects
    ! ndims - number of dimensions (should pretty much always be 3)
    ! initialConditions - an integer specifying which setup subroutine should be used (see n_body_master for more information)
    ! tunit - N-body time unit in years
    ! vunit - N-body velocity unit in km/s
    ! i - iteration variable used when looping over all objects
    ! j - iteration variable used when looping over pairs of objects
    ! k - iteration variable used when looping over dimensions
    ! count - debug variable that tracks the current timestep number
    ! encounters - 2D array in which every entry corresponds to a pair of objects. Stores a 1 if these objects are currently in a close encounter, stores a zero otherwise
    ! Note that around half of the encounters array is redundant as it will be symmetrical, but I am leaving it as an n x n array for ease of coding
    ! t - current simulation time
    ! tcount - time since previous file write (not exactly, since file writes are followed by tcount = tcount - tWrite rather than tcount = 0.)
    ! tLim - simulation time limit
    ! tWrite - time between successive file writes
    ! dtOpt - optimal timestep lengths for each object
    ! dt, dt2, dt3, dt4, dt5, recipdt2, recipdt3 - timestep length and various powers of it
    ! recip6 - 1/6
    ! recip24 - 1/24
    ! recip120 - 1/120
    ! E - total energy
    ! E0 - initial (t=0) total energy
    ! Ek - kinetic energy
    ! Ep - potential energy
    ! Eerr - fractional difference between total energy and initial total energy
    ! EerrLim - maximum allowed value of Eerr
    ! dist, dist2, dist3, dist5 - distance between two objects and various powers of it (note that dist4 is absent because it is not needed)
    ! dist2 also stores the distance from the origin squared at times
    ! dr - displacement vector between two objects (pos_i - pos_j)
    ! dv - difference in velocity between two objects (vel_i - vel_j)
    ! vmag2 - velocity magnitude squared (vel . vel)
    ! amag2 - acceleration magnitude squared (acc . acc)
    ! vdotr - dv . dr
    ! mTot - total mass of all objects
    ! m - mass of each object
    ! comPos - centre of mass
    ! comVel - centre of velocity
    ! pos - positions
    ! vel - velocities
    ! acc - accelerations
    ! jrk - jerks
    ! snp - snaps
    ! crk - crackles
    ! accPred - accelerations at predicted end-of-timestep positions and velocities (before accounting for snap and crackle)
    ! jrkPred - jerks at predicted end-of-timestep positions and velocities (before accounting for snap and crackle)
    ! eta - accuracy parameter. Optimal timestep lengths are scaled by SQRT(eta)
    ! origDists - distances from the origin
    ! vmags - velocity magnitudes
    ! closest - array that stores the distance to the nearest other object and the ID of that object (which is just the index of that object in the pos array)
    ! closestDist - keeps track of the distance to the nearest other object out of those checked so far
    ! closeThreshold - maximum distance between two objects for something to be considered a "close encounter"

    ! Opening the files that will be written to - will sort this out properly later as for now I don't think that what I've done will allow automation
    OPEN(7, file='xpositions.txt', status='replace')
    OPEN(8, file='ypositions.txt', status='replace')
    OPEN(9, file='zpositions.txt', status='replace')
    OPEN(10, file='xvelocities.txt', status='replace')
    OPEN(11, file='yvelocities.txt', status='replace')
    OPEN(12, file='zvelocities.txt', status='replace')
    OPEN(13, file='origin_dists.txt', status='replace')
    OPEN(14, file='velocity_magnitudes.txt', status='replace')
    OPEN(15, file='closest_objects.txt', status='replace')
    OPEN(16, file='energy_error.txt', status='replace')
    OPEN(17, file='close_encounters.txt', status='replace') ! Format of this will be "object i, object j, distance, time)



    ! Constants:
    recip6 = 1./6.
    recip24 = 1./24.
    recip120 = 1./120.

    ! Adjustable values
    tLim = 1.e8 / tunit
    tWrite = 1.e5 / tunit
    eta = 1.e-12
    EerrLim = 1.e-5
    closeThreshold = 1.e-3

    ! Set up initial conditions
    IF (initialConditions == 2) THEN
        CALL setup_fig_eight(pos, vel, m) ! No n as n=3 is assumed
    ELSE
        CALL setup_plummer(n, pos, vel, m)
    END IF

    ! Conversion to COM frame
    mTot = 0. ! mTot is assumed to be zero before the following code begins to run
    comPos = 0. ! comPos is assumed to be zero before the following code begins to run
    comVel = 0. ! comVel is assumed to be zero before the following code begins to run
    DO i=1,n
        mTot = mTot + m(i)
        comPos(:) = comPos(:) + (m(i) * pos(:,i))
        comVel(:) = comVel(:) + (m(i) * vel(:,i))
    END DO
    
    comPos(:) = comPos(:) / mTot
    comVel(:) = comVel(:) / mTot

    DO i=1,n
        pos(:,i) = pos(:,i) - comPos(:)
        vel(:,i) = vel(:,i) - comVel(:)
    END DO

    ! Initial total energy
    Ek = 0. ! Ek is assumed to be zero before the following code begins to run
    Ep = 0. ! Ep is assumed to be zero before the following code begins to run
    encounters = 0 ! The following code assumes this is the case before it is run
    DO i=1,n
        vmag2 = 0. ! Need to reset to zero for each new i
        DO k=1,ndims
            vmag2 = vmag2 + (vel(k,i)**2)
        END DO
        Ek = Ek + (0.5 * m(i) * vmag2)
    END DO

    DO i=1,n-1
        DO j=i+1,n ! These boundaries avoid double-counting pairs of objects
            dist2 = 0. ! Need to reset to zero for each new j
            DO k=1,ndims
                dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
            END DO
            dist = SQRT(dist2)
            Ep = Ep - ((m(i) * m(j)) / dist)

            ! Checking for a close encounter - this fits nicely into the energy calculation as it would also ideally avoid double-counting pairs of objects
            IF (dist <= closeThreshold) THEN
                encounters(i,j) = 1
                WRITE(17,*) i, j, dist, t*tunit
            END IF ! Don't need to check whether the two objects are already in an encounter here as this runs at t=0, meaning everything in encounters will be zero anyway   
        END DO
    END DO
    
    E0 = Ek + Ep

    DO i=1,n
        ! Getting distances from origin and velocity magnitudes
        dist2 = 0. ! Need to reset to zero for each new i
        vmag2 = 0. ! Need to reset to zero for each new i
        DO k=1,ndims
            dist2 = dist2 + (pos(k,i)**2)
            vmag2 = vmag2 + (vel(k,i)**2)
        END DO
        origDists(i) = SQRT(dist2)
        vmags(i) = SQRT(vmag2)
    END DO

    DO i=1,n
        ! Finding the closest other object to each object
        closestDist = 99999. ! Make this arbitrarily large at the start of each loop
        DO j=1,n
            IF (i==j) CYCLE ! Avoids having the closest object to a given object be itself
            dist2 = 0. ! Reset for each j
            DO k=1,ndims
                dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
            END DO
            dist = SQRT(dist2)
            IF (dist<closestDist) THEN
                closest(1,i) = j ! Index of the current closest object to object i
                closest(2,i) = dist
                closestDist = dist ! The closest distance to object i should now be dist
            END IF
        END DO
    END DO


    ! Writing initial conditions to files 
    WRITE(7,*) pos(1,:)
    WRITE(8,*) pos(2,:)
    WRITE(9,*) pos(3,:)
    WRITE(10,*) vel(1,:)*vunit
    WRITE(11,*) vel(2,:)*vunit
    WRITE(12,*) vel(3,:)*vunit
    WRITE(13,*) origDists(:)
    WRITE(14,*) vmags(:)*vunit
    WRITE(15,*) closest(:,:) ! Entries 1, 3, 5, etc in each line are the IDs of the nearest objects. Entries 2, 4, 6, etc are the distances to those objects
    WRITE(16,*) Eerr, t*tunit ! Note that Eerr and t*tunit will both be zero here. They are included so that all output files have the same number of lines

    ! t and tcount need to be zero at the start of the simulation
    t = 0.
    tcount = 0.

    count = 0 ! Set to zero before the first timestep, but increase by 1 at the start of every timestep
    ! Main loop
    DO
        count = count + 1 ! Tracking the current timestep number

        acc = 0. ! Needs to be reset to zero for each new timestep
        jrk = 0. ! Needs to be reset to zero for each new timestep
        ! Other variables also need to be reset to zero, but this will occur later in the timestep when those variables are actually needed

        ! Accelerations, jerks and adaptive timestep length
        DO i=1,n
            amag2 = 0. ! Needs to be reset to zero for each new i
            DO j=1,n
                IF (i==j) CYCLE ! A given object is unaffected by its own gravity
                dist2 = 0. ! Needs to be reset to zero for each new j
                vdotr = 0. ! Needs to be reset to zero for each new j
                dr(:) = pos(:,i) - pos(:,j)
                dv(:) = vel(:,i) - vel(:,j)
                DO k=1,ndims
                    vdotr = vdotr + (dv(k) * dr(k))
                    dist2 = dist2 + (dr(k)**2)
                END DO
                dist = SQRT(dist2)
                dist3 = dist**3
                dist5 = dist**5
                acc(:,i) = acc(:,i) - ((m(j) * dr(:)) / dist3)
                jrk(:,i) = jrk(:,i) + (m(j) * ((dv(:) / dist3) - ((3. * vdotr * dr(:)) / dist5)))
                DO k=1,ndims
                    amag2 = amag2 + (acc(k,i)**2)
                END DO
            END DO
            dtOpt(i) = SQRT(eta / amag2)
        END DO
        dt = MINVAL(dtOpt)
        dt2 = dt**2
        dt3 = dt**3
        dt4 = dt**4
        dt5 = dt**5
        recipdt2 = 1./dt2
        recipdt3 = 1./dt3

        ! Predicted end-of-timestep positions and velocities (before snap and crackle are accounted for)
        ! Note that the predicted positions and velocities do NOT need to be in separate arrays to the pre-timestep positions and celocities, as the
        ! pre-timestep positions and velocities are not needed for the remainder of the timestep after the predicted positions and velocities have been calculated
        DO i=1,n
            pos(:,i) = pos(:,i) + (vel(:,i) * dt) + (0.5 * acc(:,i) * dt2) + (recip6 * jrk(:,i) * dt3)
            vel(:,i) = vel(:,i) + (acc(:,i) * dt) + (0.5 * jrk(:,i) * dt2) ! Velocity MUST be updated AFTER position, as updating the position requires the
                                                                           ! pre-timestep velocity
        END DO

        ! Predicted accelerations and jerks at the predicted positions and velocities
        ! Note that the predicted accelerations and jerks DO need separate arrays to the accelerations and jerks used in predicting the end-of-timestep positions
        ! and velocities, as both pairs of accelerations and jerks are needed to estimate the snaps and crackles
        accPred = 0. ! Needs to be reset to zero in every timestep
        jrkPred = 0. ! Needs to be reset to zero in every timestep
        DO i=1,n
            ! amag2 does not need to be reset to zero as it is not used in this loop
            DO j=1,n
                IF (i==j) CYCLE
                dist2 = 0.
                vdotr = 0.
                dr(:) = pos(:,i) - pos(:,j)
                dv(:) = vel(:,i) - vel(:,j)
                DO k=1,ndims
                    vdotr = vdotr + (dv(k) * dr(k))
                    dist2 = dist2 + (dr(k)**2)
                END DO
                dist = SQRT(dist2)
                dist3 = dist**3
                dist5 = dist**5
                accPred(:,i) = accPred(:,i) - ((m(j) * dr(:)) / dist3)
                jrkPred(:,i) = jrkPred(:,i) + (m(j) * ((dv(:) / dist3) - ((3. * vdotr * dr(:)) / dist5)))
            END DO
        END DO

        ! Estimating the snap and crackle of each object. This is in a separate DO loop to the calculation of accPred and jrkPred because the estimation of the
        ! snap and crackle for a given object requires knowledge of accPred and jrkPred for all other objects
        snp = 0. ! Needs to be reset to zero in every timestep
        crk = 0. ! Needs to be reset to zero in every timestep
        DO i=1,n
            DO j=1,n
                IF (i==j) CYCLE
                ! Note that the subtraction in the snap estimation is outside the parentheses so the jerk term is added rather than subtracted
                snp(:,i) = snp(:,i) - (recipdt2 * ((6. * (acc(:,j) - accPred(:,j))) + (dt * ((4. * jrk(:,j)) + (2. * jrkPred(:,j))))))
                crk(:,i) = crk(:,i) + (recipdt3 * ((12. * (acc(:,j) - accPred(:,j))) + (6. * dt * (jrk(:,j) + jrkPred(:,j)))))
            END DO
        END DO


        ! Position and velocity corrections. This COULD be including in the previous DO loop after cycling through all j's for a given i as the snap/crackle
        ! calculations in the actual loop do not involve the predicted positions/velocities of the relevant objects (the contributions of the predicted positions/velocities
        ! are contained within accPred and jrkPred, which were calculated beforehand). Updating the positions and velocities during the same DO loop as
        ! the snaps and crackles are esimated shouldn't affect the snap and crackle estimations, and won't have any effect on the calculated positions and velocities
        ! as long as they are updated AFTER cycling through all j's for a given i. I have decided against putting them in the same DO loop
        ! in order to make the code easier to follow
        DO i=1,n
            ! The contributions from orders up to the jerk are already included within the current values of pos and vel for a given object before they are updated here,
            ! so only the contributions from the snap and crackle need to be accounted for here
            pos(:,i) = pos(:,i) + (recip24 * snp(:,i) * dt4) + (recip120 * crk(:,i) * dt5)
            vel(:,i) = vel(:,i) + (recip6 * snp(:,i) * dt3) + (recip24 * crk(:,i) * dt4)
        END DO

        ! Time incrementation
        t = t + dt
        tcount = tcount + dt

        ! Checking energy error (this should be done every timestep because the code should stop running if the energy error gets too large)
        ! The outputs from before the energy error exceeds the limit are fine to keep, but anything after the energy error exceeds the limit should be discarded
        ! so time is saved by just ending the simulation if the energy error exceeds the limit. Data will be written to the files for the timestep in which the energy
        ! error becomes too large because it could be useful to see if the outputs give some insight as to why the energy error became so large (especially if it is sudden)
        Ek = 0. ! Needs to be reset to zero in every timestep
        Ep = 0. ! Needs to be reset to zero in every timestep

        DO i=1,n
            vmag2 = 0.
            DO k=1,ndims
                vmag2 = vmag2 + (vel(k,i)**2)
            END DO
            Ek = Ek + (0.5 * m(i) * vmag2)
        END DO

        DO i=1,n-1
            DO j=i+1,n ! These boundaries avoid double-counting pairs of objects
                dist2 = 0. ! Need to reset to zero for each new j
                DO k=1,ndims
                    dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
                END DO
                dist = SQRT(dist2)
                Ep = Ep - ((m(i) * m(j)) / dist)

                ! Checking for a close encounter - this fits nicely into the energy calculation as it would also ideally avoid double-counting pairs of objects
                IF (dist <= closeThreshold) THEN
                    IF (encounters(i,j)==0) THEN
                        encounters(i,j) = 1
                        WRITE(17,*) i, j, dist, t*tunit
                    END IF
            
                ELSE
                    IF (encounters(i,j)==1) THEN
                        encounters(i,j) = 0
                        WRITE(17,*) i, j, dist, t*tunit ! This is technically after the encounter ends but it feels easier to code it this way. I may change it in the future
                    END IF
                END IF
            
            END DO
        END DO
        
        E = Ek + Ep
        Eerr = ABS((E-E0)/E0)

        IF (Eerr >= EerrLim) THEN ! Stop the code if the energy error is too large. This takes priority over the regular file write conditions that depend on either
                                  ! the simulated time elapsed since the beginning of the simulation (file write if t>=tlim, then end the simulation) or on the 
                                  ! file write interval being reached (file write if tcount>=tWrite, then subtract tWrite from tcount)

            DO i=1,n
                ! Getting distances from origin and velocity magnitudes
                dist2 = 0. ! Need to reset to zero for each new i
                vmag2 = 0. ! Need to reset to zero for each new i
                DO k=1,ndims
                    dist2 = dist2 + (pos(k,i)**2)
                    vmag2 = vmag2 + (vel(k,i)**2)
                END DO
                origDists(i) = SQRT(dist2)
                vmags(i) = SQRT(vmag2)
            END DO

            DO i=1,n
                ! Finding the closest other object to each object
                closestDist = 99999. ! Make this arbitrarily large at the start of each loop
                DO j=1,n
                    IF (i==j) CYCLE ! Avoids having the closest object to a given object be itself
                    dist2 = 0. ! Reset for each j
                    DO k=1,ndims
                        dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
                    END DO
                    dist = SQRT(dist2)
                    IF (dist<closestDist) THEN
                        closest(1,i) = j ! Index of the current closest object to object i
                        closest(2,i) = dist
                        closestDist = dist ! The closest distance to object i should now be dist
                    END IF
                END DO
            END DO

            ! FILE WRITES
            WRITE(7,*) pos(1,:)
            WRITE(8,*) pos(2,:)
            WRITE(9,*) pos(3,:)
            WRITE(10,*) vel(1,:)*vunit
            WRITE(11,*) vel(2,:)*vunit
            WRITE(12,*) vel(3,:)*vunit
            WRITE(13,*) origDists(:)
            WRITE(14,*) vmags(:)*vunit
            WRITE(15,*) closest(:,:) ! Entries 1, 3, 5, etc in each line are the IDs of the nearest objects. Entries 2, 4, 6, etc are the distances to those objects
            WRITE(16,*) Eerr, t*tunit
            EXIT
        END IF

        ! File writes
        IF (t >= tLim) THEN ! This takes priority over the coundition (tcount >= tWrite) in the event that both are satisfied

            DO i=1,n
                ! Getting distances from origin and velocity magnitudes
                dist2 = 0. ! Need to reset to zero for each new i
                vmag2 = 0. ! Need to reset to zero for each new i
                DO k=1,ndims
                    dist2 = dist2 + (pos(k,i)**2)
                    vmag2 = vmag2 + (vel(k,i)**2)
                END DO
                origDists(i) = SQRT(dist2)
                vmags(i) = SQRT(vmag2)
            END DO

            DO i=1,n
                ! Finding the closest other object to each object
                closestDist = 99999. ! Make this arbitrarily large at the start of each loop
                DO j=1,n
                    IF (i==j) CYCLE ! Avoids having the closest object to a given object be itself
                    dist2 = 0. ! Reset for each j
                    DO k=1,ndims
                        dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
                    END DO
                    dist = SQRT(dist2)
                    IF (dist<closestDist) THEN
                        closest(1,i) = j ! Index of the current closest object to object i
                        closest(2,i) = dist
                        closestDist = dist ! The closest distance to object i should now be dist
                    END IF
                END DO
            END DO

            ! FILE WRITES
            WRITE(7,*) pos(1,:)
            WRITE(8,*) pos(2,:)
            WRITE(9,*) pos(3,:)
            WRITE(10,*) vel(1,:)*vunit
            WRITE(11,*) vel(2,:)*vunit
            WRITE(12,*) vel(3,:)*vunit
            WRITE(13,*) origDists(:)
            WRITE(14,*) vmags(:)*vunit
            WRITE(15,*) closest(:,:) ! Entries 1, 3, 5, etc in each line are the IDs of the nearest objects. Entries 2, 4, 6, etc are the distances to those objects
            WRITE(16,*) Eerr, t*tunit
            EXIT

        ELSE IF (tcount >= tWrite) THEN

            DO i=1,n
                ! Getting distances from origin and velocity magnitudes
                dist2 = 0. ! Need to reset to zero for each new i
                vmag2 = 0. ! Need to reset to zero for each new i
                DO k=1,ndims
                    dist2 = dist2 + (pos(k,i)**2)
                    vmag2 = vmag2 + (vel(k,i)**2)
                END DO
                origDists(i) = SQRT(dist2)
                vmags(i) = SQRT(vmag2)
            END DO

            DO i=1,n
                ! Finding the closest other object to each object
                closestDist = 99999. ! Make this arbitrarily large at the start of each loop
                DO j=1,n
                    IF (i==j) CYCLE ! Avoids having the closest object to a given object be itself
                    dist2 = 0. ! Reset for each j
                    DO k=1,ndims
                        dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
                    END DO
                    dist = SQRT(dist2)
                    IF (dist<closestDist) THEN
                        closest(1,i) = j ! Index of the current closest object to object i
                        closest(2,i) = dist
                        closestDist = dist ! The closest distance to object i should now be dist
                    END IF
                END DO
            END DO

            ! FILE WRITES
            WRITE(7,*) pos(1,:)
            WRITE(8,*) pos(2,:)
            WRITE(9,*) pos(3,:)
            WRITE(10,*) vel(1,:)*vunit
            WRITE(11,*) vel(2,:)*vunit
            WRITE(12,*) vel(3,:)*vunit
            WRITE(13,*) origDists(:)
            WRITE(14,*) vmags(:)*vunit
            WRITE(15,*) closest(:,:) ! Entries 1, 3, 5, etc in each line are the IDs of the nearest objects. Entries 2, 4, 6, etc are the distances to those objects
            WRITE(16,*) Eerr, t*tunit
            tcount = tcount - tWrite ! Not just resetting tcount to zero because this could cause the file outputs to drift from the intended times over the course
                                     ! of the simulation. Subtracting twrite keeps the file outputs at times as close to multiples of twrite as possible
        END IF
         
    END DO

    RETURN
    
END SUBROUTINE n_body_simulation

! Gives the N-body velocity unit in km/s and the time unit in years
SUBROUTINE n_body_conversion(munit, runit, tunit, vunit)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: munit, runit
    DOUBLE PRECISION, INTENT(OUT) :: tunit, vunit
    DOUBLE PRECISION :: G, year

    G = 6.67e-11 ! Gravitational constant in SI units
    year = 365.24 * 86400. ! 86400 seconds in a day, 365.24 days in a year
    tunit = SQRT((runit**3) / (G * munit)) / year ! N-body time unit in years
    vunit = runit / (1.e3 * tunit * year)
    
    
    
    ! N-body velocity unit in km/s (the factor of 1.e3 on the denominator converts from m/s to km/s)

    RETURN
    
END SUBROUTINE n_body_conversion

! Plummer sphere (3 dimensions only). Method from Aarseth, Hénon and Wielen (1974)
SUBROUTINE setup_plummer(n, pos, vel, m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i
    DOUBLE PRECISION, INTENT(OUT) :: pos(1:3, 1:n), vel(1:3, 1:n), m(1:n)
    DOUBLE PRECISION :: ranNum, ranNum2, twopi, twothirds, rad, sqrtr2z2, vesc, root2, q2, gq, vmag


    twopi = 8. * ATAN(1.)   ! 2 * pi
    ranNum = 0.             ! Random number in a uniform distribution between 0 and 1
    ranNum2 = 0.            ! Another random number in a uniform distribution between 0 and 1, used when a new random number needs to be generated while the value stored
                            ! in ranNum is also still needed
    twothirds = 2./3.       ! 2/3
    rad = 0.                ! Distance of a given object from the origin. Not an array as the generation of initial conditions for a given object only needs the value
                            ! of rad for the object itself
    sqrtr2z2 = 0.           ! SQRT(r^2 - z^2). Also stores SQRT(v^2 - v_z^2) as this is value is only needed for a given object once SQRT(r^2 - z^2) is no longer needed
    vesc = 0.               ! Escape velocity at a given object's position, which represents the maximum 'allowed' velocity at that position. Not an array for the same
                            ! reasons as rad
    vmag = 0.               ! Velocity magnitude of a given object. Not an array for the same reasons as rad
    q2 = 0.                 ! q^2, where q = vmag/vesc
    gq = 0.                 ! g(q) = q^2 - (1-q^2)^7/2
    root2 = SQRT(2.)        ! Square root of 2

    m = 1. ! All objects will be 1 solar mass (for now)
    DO i=1,n
        CALL RANDOM_NUMBER(ranNum)
        rad = SQRT(1./((ranNum**(-1.*twothirds)) - 1.))
        CALL RANDOM_NUMBER(ranNum)
        pos(3,i) = rad * (1. - (2.*ranNum))
        CALL RANDOM_NUMBER(ranNum) ! Need a different random number for the x and y coordinates to the one used for the z coordinates
        sqrtr2z2 = SQRT((rad**2) - (pos(3,i)**2))
        pos(1,i) = COS(twopi * ranNum) * sqrtr2z2
        pos(2,i) = SIN(twopi * ranNum) * sqrtr2z2
        vesc = root2 * ((1. + (rad**2))**(-0.25))
        DO
            ! This will run continuously until a certain condition is met (this avoids having to use GOTO)
            CALL RANDOM_NUMBER(ranNum)
            q2 = ranNum**2
            gq = q2 * ((1. - q2)**3.5)
            CALL RANDOM_NUMBER(ranNum2)
            IF ((0.1 * ranNum2) < gq) THEN
                EXIT ! This ends the current DO loop (which would otherwise run infinitely)
            END IF
        END DO
        vmag = ranNum * vesc ! Note this uses ranNum, NOT ranNum2
        CALL RANDOM_NUMBER(ranNum)
        vel(3,i) = vmag * (1. - (2.*ranNum))
        CALL RANDOM_NUMBER(ranNum) ! Need a different random number for the x and y velocity components to the one used for the z velocity component
        sqrtr2z2 = SQRT((vmag**2) - (vel(3,i)**2)) ! Reusing this variable name as the value it was set to earlier in the loop is no longer needed at this point
        vel(1,i) = COS(twopi * ranNum) * sqrtr2z2
        vel(2,i) = SIN(twopi * ranNum) * sqrtr2z2
    END DO

    RETURN

END SUBROUTINE setup_plummer

SUBROUTINE setup_fig_eight(pos, vel, m)
    DOUBLE PRECISION, INTENT(OUT) :: pos(1:3, 1:3), vel(1:3, 1:3), m(1:3)

    m = 1.

    pos(1,1) = 0.9700436
    pos(2,1) = -0.24308753
    pos(1,2) = -1. * pos(1,1)
    pos(2,2) = -1. * pos(2,1)
    pos(1,3) = 0. ! Should already be the case on calling this subroutine, but this line is included on the off chance it isn't the case
    pos(2,3) = 0. ! See above comment
    pos(3,:) = 0. ! See comment two lines above

    vel(1,1) = 0.466203685
    vel(2,1) = 0.43236573
    vel(1,2) = vel(1,1)
    vel(2,2) = vel(2,1)
    vel(1,3) = -2. * vel(1,1)
    vel(2,3) = -2. * vel(2,1)
    vel(3,:) = 0. ! Should already be the case on calling this subroutine, but this line is included on the off chance it isn't the case

    RETURN

END SUBROUTINE setup_fig_eight