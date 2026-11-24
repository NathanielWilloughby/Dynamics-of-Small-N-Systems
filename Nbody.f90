! Main program - automates multiple N-body simulations
PROGRAM n_body_master
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    INTEGER :: nsims, i
    REAL(dp) :: munit, runit, vunit, tunit

    ! Some adjustable values
    nsims = 1 ! Number of simulations to run
    munit = 1.99e30_dp ! Equal to one N-body mass unit (which is currently set to one solar mass)
    runit = 3.086e16_dp ! Equal to one N-body length unit (which is currently set to one parsec)

    ! Get velocity unit in km/s and time unit in Myr
    CALL n_body_conversion(munit, runit, tunit, vunit)

    ! Run N-body simulations by repeatedly calling n_body_simulation(n, ndims, initialConditions, tunit, vunit)
    ! n - number of objects
    ! ndims - number of dimensions (probably should always be 3)
    ! initialConditions - an integer used to determine which setup subroutine to use
    ! tunit - N-body time unit in Myr
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
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    INTEGER, INTENT(IN) :: n, ndims, initialConditions
    REAL(dp), INTENT(IN) :: tunit, vunit
    INTEGER :: i, j, k, timestepNum, encounters(1:n, 1:n), maxPower, timestepMultipliers(1:n), sync, numPossibleTimestepLengths, trialdt
    REAL(dp) :: dtBase, t, tcount, tLim, tWrite, dtOpt(1:n), dt, dt2, dt3, recipdt2, recipdt3, recip6, recip24, recip120, E, E0, Ek, Ep, Eerr, EerrLim
    REAL(dp) :: dist, dist2, dist3, dist5, dr(1:ndims), dv(1:ndims), vmag2, amag2, jmag2, smag2, cmag2, vdotr, mTot, m(1:n), comPos(1:ndims), comVel(1:ndims)
    REAL(dp) :: pos(1:ndims, 1:n), vel(1:ndims, 1:n), acc(1:ndims, 1:n), jrk(1:ndims, 1:n), snp(1:ndims, 1:n), crk(1:ndims, 1:n), endPos(1:ndims, 1:n), endVel(1:ndims, 1:n)
    REAL(dp) :: accPred(1:ndims, 1:n), jrkPred(1:ndims, 1:n), eta, origDists(1:n), vmags(1:n), closest(1:2, 1:n), closestDist, closeThreshold
    REAL(dp) :: desiredVirRat, virRat, timeSinceLast(1:n), reciplog2, dtMin, exponent
    LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: doThisObject
    INTEGER, ALLOCATABLE, DIMENSION(:) :: dtAllowedMultipliers

    ! Variable definitions:
    ! n - number of objects
    ! ndims - number of dimensions (should pretty much always be 3)
    ! initialConditions - an integer specifying which setup subroutine should be used (see n_body_master for more information)
    ! tunit - N-body time unit in years
    ! vunit - N-body velocity unit in km/s
    ! i - iteration variable used when looping over all objects
    ! j - iteration variable used when looping over pairs of objects
    ! k - iteration variable used when looping over dimensions
    ! timestepNum - tracks the current timestep number, used for handling multi-particle timesteps and for debugging
    ! encounters - 2D array in which every entry corresponds to a pair of objects. Stores a 1 if these objects are currently in a close encounter, stores a zero otherwise
    ! Note that around half of the encounters array is redundant as it will be symmetrical, but I am leaving it as an n x n array for ease of coding
    ! maxPower - Limit on how many powers of 2 the timestep length can cover
    ! timestepMultipliers - timestep length for object i is timestepMultipliers(i) * dtBase. Ranges from 2^-maxPower to 2^maxPower
    ! sync - number of timesteps between times when ALL possible timestep lengths will be synchronised (will be equal to (2^maxPower)^2)
    ! numPossibleTimestepLengths - number of allowed timestep lengths, will be equal to 1 + (2*maxPower)
    ! trialdt - temporarily stores the optimal timestep length multiplier for an object whilst checking if a move to a new timestep length is allowed
    ! dtBase - base (default) timestep length
    ! t - current simulation time
    ! tcount - time since previous file write (not exactly, since file writes are followed by tcount = tcount - tWrite rather than tcount = 0.)
    ! tLim - simulation time limit
    ! tWrite - time between successive file writes
    ! dtOpt - optimal timestep lengths for each object
    ! dt, dt2, dt3, recipdt2, recipdt3 - timestep length and various powers of it
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
    ! jmag2 - jerk magnitude squared (jrk . jrk)
    ! smag2 - snap magnitude squared (snp . snp)
    ! cmag2 - crackle magnitude squared (crk . crk)
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
    ! endPos - end-of-timestep positions (used for interpolating when handling short-timestep objects)
    ! endVel - end-of-timestep velocities (used for interpolating when handling short-timestep objects)
    ! accPred - accelerations at predicted end-of-timestep positions and velocities (before accounting for snap and crackle)
    ! jrkPred - jerks at predicted end-of-timestep positions and velocities (before accounting for snap and crackle)
    ! eta - accuracy parameter. Optimal timestep lengths are scaled by SQRT(eta)
    ! origDists - distances from the origin
    ! vmags - velocity magnitudes
    ! closest - array that stores the distance to the nearest other object and the ID of that object (which is just the index of that object in the pos array)
    ! closestDist - keeps track of the distance to the nearest other object out of those checked so far
    ! closeThreshold - maximum distance between two objects for something to be considered a "close encounter"
    ! desiredVirRat - desired virial ratio - may make this a parameter for this subroutine in the future
    ! virRat - actual virial ratio (at t=0 before the velocities are scaled to make the actual virial ratio equal to the desired virial ratio)
    ! dtAllowedMultipliers - allowed timestep length multipliers
    ! timeSinceLast - amount of simulated time since each object last had its position and velocity updated
    ! reciplog2 - 1 / log(2). Here 'log' is the base 10 logarithm. This is used in converting to a base 2 logarithm
    ! dtMin - minimum possible timestep length - shorthand for dtBase * dtAllowedMultipliers(1)
    ! exponent - a given timestep multiplier is 2^exponent
    ! doThisObject - 2D array indicating which objects have each timestep length. A chunk of this array that corresponds to a given timestep length will be passed
    ! into the timestep subroutine


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


    ! Adjustable values
    tLim = 1.e2_dp / tunit
    tWrite = 1.e-1_dp / tunit
    eta = 1.e-6_dp
    EerrLim = 1.e-5_dp
    closeThreshold = 1.e-3_dp
    desiredVirRat = 1._dp
    dtBase = 1.e-7_dp ! This value is going to be hardcoded for now, but I may make it able to change based on the initial conditions at some point
    maxPower = 6

    ! Constants
    recip6 = 1._dp/6._dp
    recip24 = 1._dp/24._dp
    recip120 = 1._dp/120._dp
    reciplog2 = 1._dp / log(2._dp)
    numPossibleTimestepLengths = 1 + (2 * maxPower)
    sync = (2 ** maxPower) ** 2

    ! Setting up doThisObject
    ALLOCATE(doThisObject(1:n, 1:numPossibleTimestepLengths))
    doThisObject = .false. ! All false by default

    ! Setting up dtAllowedMultipliers
    ALLOCATE(dtAllowedMultipliers(1:numPossibleTimestepLengths))
    DO i=(-1 * maxPower),maxPower
        dtAllowedMultipliers(i+maxPower+1) = 2 ** i ! i+maxPower+1 should convert to the correct position in the dtAllowedMultipliers array
    END DO
    dtMin = dtBase * dtAllowedMultipliers(1)


    ! Set up initial conditions
    IF (initialConditions == 2) THEN
        CALL setup_fig_eight(pos, vel, m) ! No n as n=3 is assumed
    ELSE
        CALL setup_plummer(n, pos, vel, m)
    END IF

    ! Conversion to COM frame
    mTot = 0._dp ! mTot is assumed to be zero before the following code begins to run
    comPos = 0._dp ! comPos is assumed to be zero before the following code begins to run
    comVel = 0._dp ! comVel is assumed to be zero before the following code begins to run
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
    Ek = 0._dp ! Ek is assumed to be zero before the following code begins to run
    Ep = 0._dp ! Ep is assumed to be zero before the following code begins to run
    encounters = 0 ! The following code assumes this is the case before it is run
    DO i=1,n
        vmag2 = 0._dp ! Need to reset to zero for each new i
        DO k=1,ndims
            vmag2 = vmag2 + (vel(k,i)**2)
        END DO
        Ek = Ek + (0.5_dp * m(i) * vmag2)
    END DO

    DO i=1,n-1
        DO j=i+1,n ! These boundaries avoid double-counting pairs of objects
            dist2 = 0._dp ! Need to reset to zero for each new j
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
    
    ! Scaling velocities to make the virial ratio of the system be equal to the desired value
    virRat = (2._dp * Ek) / ABS(Ep)
    vel(:,:) = vel(:,:) * SQRT(desiredVirRat / virRat)

    ! Recalculating kinetic energy
    Ek = 0._dp
    DO i=1,n
        vmag2 = 0._dp
        DO k=1,ndims
            vmag2 = vmag2 + (vel(k,i)**2)
        END DO
        Ek = Ek + (0.5_dp * m(i) * vmag2)
    END DO

    E0 = Ek + Ep

    DO i=1,n
        ! Getting distances from origin and velocity magnitudes
        dist2 = 0._dp ! Need to reset to zero for each new i
        vmag2 = 0._dp ! Need to reset to zero for each new i
        DO k=1,ndims
            dist2 = dist2 + (pos(k,i)**2)
            vmag2 = vmag2 + (vel(k,i)**2)
        END DO
        origDists(i) = SQRT(dist2)
        vmags(i) = SQRT(vmag2)
    END DO

    DO i=1,n
        ! Finding the closest other object to each object
        closestDist = 99999._dp ! Make this arbitrarily large at the start of each loop
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

    ! t, tcount and timeSinceLast need to be zero at the start of the simulation
    t = 0._dp
    tcount = 0._dp
    timeSinceLast = 0._dp

    ! All objects will initially have the same timestep length
    timestepMultipliers = 0

    timestepNum = 0 ! Set to zero before the first timestep, but increase by 1 at the end of every shortest-length timestep
                    ! IMPORTANT: This means that when using timestepNum for debugging, the first timestep has timestepNum = 0

    ! Main loop

    ! Calculating initial optimal timestep lengths based on the accelerations, jerks, snaps and crackles of each object after one base timestep has elapsed
    ! This is not a "real" timestep, and instead is used to assign initial values to each index in the timesteps array. This should be treated as part of the setup
    ! Comments will be mostly omitted from this section (except where something is significantly different to the main DO loop)
    ! See the equivalent parts in the main DO loop for comments

    acc = 0._dp
    jrk = 0._dp
    endPos(:,:) = pos(:,:) ! Using endPos and endVel in this section so that the values stored in pos and vel can be kept for the first "real" timestep
    endVel(:,:) = vel(:,:)
    DO i=1,n
        DO j=1,n
            IF (i==j) CYCLE
            dist2 = 0._dp
            vdotr = 0._dp
            dr(:) = endPos(:,i) - endPos(:,j)
            dv(:) = endVel(:,i) - endVel(:,j)
            DO k=1,ndims
                vdotr = vdotr + (dv(k) * dr(k))
                dist2 = dist2 + (dr(k)**2)
            END DO
            dist = SQRT(dist2)
            dist3 = dist**3
            dist5 = dist**5
            acc(:,i) = acc(:,i) - ((m(j) * dr(:)) / dist3)
            jrk(:,i) = jrk(:,i) + (m(j) * ((dv(:) / dist3) - ((3._dp * vdotr * dr(:)) / dist5)))
        END DO
    END DO

    ! No need to calculate optimal timestep length as this section is assuming that dt=dtBase
    ! dt4 and dt5 are not needed here as the positions and velocities are not being fully updated since this is not a "real" timestep
    dt = dtBase
    dt2 = dt**2
    dt3 = dt**3
    recipdt2 = 1._dp/dt2
    recipdt3 = 1._dp/dt3
    DO i=1,n
        endPos(:,i) = endPos(:,i) + (endVel(:,i) * dt) + (0.5 * acc(:,i) * dt2) + (recip6 * jrk(:,i) * dt3)
        endVel(:,i) = endVel(:,i) + (acc(:,i) * dt) + (0.5 * jrk(:,i) * dt2)
    END DO

    accPred = 0._dp
    jrkPred = 0._dp
    DO i=1,n
        DO j=1,n
            IF (i==j) CYCLE
            dist2 = 0._dp
            vdotr = 0._dp
            dr(:) = endPos(:,i) - endPos(:,j)
            dv(:) = endVel(:,i) - endVel(:,j)
            DO k=1,ndims
                vdotr = vdotr + (dv(k) * dr(k))
                dist2 = dist2 + (dr(k)**2)
            END DO
            dist = SQRT(dist2)
            dist3 = dist**3
            dist = dist**5
            accPred(:,i) = accPred(:,i) - ((m(j) * dr(:)) / dist3)
            jrkPred(:,i) = jrkPred(:,i)+ (m(j) * ((dv(:) / dist3) - ((3._dp * vdotr * dr(:)) / dist5)))
        END DO
    END DO

    snp = 0._dp
    crk = 0._dp
    DO i=1,n
        DO j=1,n
            IF (i==j) CYCLE
            snp(:,i) = snp(:,i) - (recipdt2 * ((6._dp * (acc(:,j) - accPred(:,j))) + (dt * ((4._dp * jrk(:,j)) + (2._dp * jrkPred(:,j))))))
            crk(:,i) = crk(:,i) + (recipdt3 * ((12._dp * (acc(:,j) - accPred(:,j))) + (6._dp * dt * (jrk(:,j) + jrkPred(:,j)))))
        END DO
    END DO

    ! Do not need to update endPos or endVel here as this is not a "real" timestep

    DO i=1,n
        amag2 = 0._dp
        jmag2 = 0._dp
        smag2 = 0._dp
        cmag2 = 0._dp
        DO k=1,ndims
            amag2 = amag2 + (acc(k,i)**2)
            jmag2 = jmag2 + (jrk(k,i)**2)
            smag2 = smag2 + (snp(k,i)**2)
            cmag2 = cmag2 + (crk(k,i)**2)
        END DO
        dtOpt(i) = eta * (SQRT(amag2 * smag2) + jmag2)
        dtOpt(i) = dtOpt(i) / (SQRT(jmag2 * cmag2) + smag2)
        dtOpt(i) = SQRT(dtOpt(i))
        CALL closest_timestep(maxPower, dtBase, dtAllowedMultipliers, dtOpt(i), timestepMultipliers(i), numPossibleTimestepLengths)
    END DO

    ! TODO: REWRITE THE MAIN LOOP ACCOUNTING FOR MULTI-PARTICLE TIMESTEPS

    ! The following code assumes this is the case at the start of the first timestep
    endPos(:,:) = pos(:,:)
    endVel(:,:) = vel(:,:)

    DO

        acc = 0._dp ! Needs to be reset to zero for each new timestep
        jrk = 0._dp ! Needs to be reset to zero for each new timestep
        ! Other variables also need to be reset to zero, but this will occur later in the timestep when those variables are actually needed

        ! Subroutine for performing a timestep of a given length for some or all of the objects is
        ! timestep(n, ndims, m, pos, vel, acc, jrk, snp, crk, endPos, endVel, dt, dtBase, doThisObject, timestepMultipliers, timeSinceLast)

        ! Filling up doThisObject
        DO i=1,n
            exponent = LOG(REAL(timestepMultipliers(i))) * reciplog2 ! Base 2 logarithm of timestepMultipliers(i)
            doThisObject(i, (INT(exponent) + maxPower + 1)) = .true. ! exponent should be x.0000000, so converting to an integer should be fine
                                                                     ! Adding maxPower + 1 means that the minimum possible timestep length has a second index of 1
                                                                     ! in the doThisObject array
        END DO

        DO i=1,numPossibleTimestepLengths
            IF (MODULO(timestepNum, (sync / (2**(i-1)))) == 0) THEN
                ! The above condition means that a longer timestep of a given length synchronises with the shortest possible timestep length on this timestep
                ! Objects with the given longer timestep length will have their positions and velocities updated
                ! Note that the case where i=numPossibleTimestepLengths corresponds to the shortest possible timestep length, and so the condition in the above
                ! IF statement will always be true
                DO j=1,n
                    IF (doThisObject(j,(numPossibleTimestepLengths + 1 - i)) .eqv. .true.) THEN
                        ! numPossibleTimestepLengths + 1 - i means that for i=1, doThisObject(j,i), which corresponds to the longest possible timestep length, will
                        ! be considered

                        ! endPos and endVel are the arrays that are actually updated by the timestep subroutine. The start-of-timestep values (which will be the values
                        ! stored in endPos and endVel before they are updated - as the unupdated values are those from the end of the previous timestep)
                        ! are needed for interpolation of positions and velocities when considering other objects
                        ! The values in endPos(:,j) and endVel(:,j) before updating are not needed when updating the positions and velocities of object j
                        pos(:,j) = endPos(:,j)
                        vel(:,j) = endVel(:,j)

                        timeSinceLast(j) = 0._dp ! Setting this to zero here as the value in timeSinceLast(j) is not needed when object j is having its
                                                 ! position and velocity updated by the timestep subroutine. Setting to zero here therefore avoids having a
                                                 ! separate DO loop after calling the timestep subroutine to reset timeSinceLast to zero for the objects that were updated
                    END IF
                END DO
                CALL timestep(n, ndims, m, pos, vel, acc, jrk, snp, crk, endPos, endVel, (dtBase * dtAllowedMultipliers(i)), dtBase, &
                doThisObject(:, (numPossibleTimestepLengths + 1 - i)), timestepMultipliers, timeSinceLast)
            END IF
        END DO

        ! Updating timeSinceLast for all objects that aren't on the shortest possible timestep length
        DO i=1,n
            IF (doThisObject(i,1) .eqv. .false.) THEN
                timeSinceLast(i) = timeSinceLast(i) + dtMin
            END IF
        END DO

        ! Updating timestep lengths
        ! NOTE: Objects can only move to a longer timestep length in cases where the timestep length they are moving to synchronises with the one they are currently on
        ! at the start of the next timestep
        DO i=1,n
            amag2 = 0._dp
            jmag2 = 0._dp
            smag2 = 0._dp
            cmag2 = 0._dp
            DO k=1,ndims
                amag2 = amag2 + (acc(k,i)**2)
                jmag2 = jmag2 + (jrk(k,i)**2)
                smag2 = smag2 + (snp(k,i)**2)
                cmag2 = cmag2 + (crk(k,i)**2)
            END DO
            dtOpt(i) = eta * (SQRT(amag2 * smag2) + jmag2)
            dtOpt(i) = dtOpt(i) / (SQRT(jmag2 * cmag2) + smag2)
            dtOpt(i) = SQRT(dtOpt(i))
            CALL closest_timestep(maxPower, dtBase, dtAllowedMultipliers, dtOpt(i), trialdt, numPossibleTimestepLengths)

            IF (trialdt < timestepMultipliers(i)) THEN
                ! Moving to a shorter timestep is always allowed as they will always synchronise with the current timestep length for a given object
                timestepMultipliers(i) = trialdt
            ELSE IF (trialdt > timestepMultipliers(i)) THEN
                ! Not using an ELSE statement because I don't the following to happen if trialdt==timestepMultipliers(i). When that is the case, no change will be made
                ! to timestepMultipliers(i) regardless

                ! Moving to a longer timestep is only allowed when the target timestep length synchronises with the current timestep length for that object at the start of
                ! the next timestep
                ! Due to the way that the timestep multipliers have been defined (they have consistent multiplication factors between them), if the shortest possible
                ! timestep length synchronises with the target timestep length, then the current timestep length will also synchronise (this is actually
                ! the case if any timestep length that is shorter than the current one synchronises with the target timestep length). This means that the code can just
                ! check whether the target timestep length synchronises with the shortest possible timestep length, preventing me from having to work out how to check
                ! whether two arbitrary timestep lengths synchronise

                IF (MODULO((timestepNum+1), (trialdt * (2**maxPower))) == 0) THEN
                    ! The above condition checks whether the target timestep length and the shortest possible timestep length synchronise at the start
                    ! of the next timestep (hence the timestepNum+1). The number of shortest possible timestep lengths per target timestep length is just
                    ! the ratio between their multipliers, which is trialdt/2^-maxPower, which is equivalent to trialdt * 2^maxPower

                    ! timestepNum+1 will mean that objects only move up to a longer timestep if the target and current timestep lengths synchronise at the start of
                    ! the next shortest possible timestep, which is what I want to be the case (this comment is here so that I don't forget about this and panic
                    ! about it in the future)

                    timestepMultipliers(i) = trialdt
                END IF
            END IF
        END DO

        ! File outputs should only be done when all timestep lengths are synchronised at the start of the next timestep of shortest possible length

        IF (MODULO((timestepNum+1), sync) == 0) THEN

            ! Calculating total energy
            Ek = 0._dp
            Ep = 0._dp

            DO i=1,n
                vmag2 = 0._dp
                DO k=1,ndims
                    vmag2 = vmag2 + (endVel(k,i)**2)
                END DO
                Ek = Ek + (0.5_dp * m(i) * vmag2)
            END DO

            DO i=1,n-1
                DO j=i+1,n
                    dist2 = 0._dp
                    DO k=1,ndims
                        dist2 = dist2 + ((endPos(k,i) - endPos(k,j))**2)
                    END DO
                    dist = SQRT(dist2)
                    Ep = Ep - ((m(i) * m(j)) / dist)

                    ! Could also check for a close encounter, but not doing so right now
                END DO
            END DO

            E = Ek + Ep
            Eerr = ABS((E-E0)/E0)

            IF (Eerr >= EerrLim) THEN ! Stop the code if the energy error is too large. This takes priority over the regular file write conditions that depend on either
                                      ! the simulated time elapsed since the beginning of the simulation (file write if t>=tlim, then end the simulation) or on the 
                                      ! file write interval being reached (file write if tcount>=tWrite, then subtract tWrite from tcount)

                DO i=1,n
                    ! Getting distances from origin and velocity magnitudes
                    dist2 = 0._dp ! Need to reset to zero for each new i
                    vmag2 = 0._dp ! Need to reset to zero for each new i
                    DO k=1,ndims
                        dist2 = dist2 + (endPos(k,i)**2)
                        vmag2 = vmag2 + (endVel(k,i)**2)
                    END DO
                    origDists(i) = SQRT(dist2)
                    vmags(i) = SQRT(vmag2)
                END DO

                DO i=1,n
                    ! Finding the closest other object to each object
                    closestDist = 99999._dp ! Make this arbitrarily large at the start of each loop
                    DO j=1,n
                        IF (i==j) CYCLE ! Avoids having the closest object to a given object be itself
                        dist2 = 0._dp ! Reset for each j
                        DO k=1,ndims
                            dist2 = dist2 + ((endPos(k,i) - endPos(k,j))**2)
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
                WRITE(7,*) endPos(1,:)
                WRITE(8,*) endPos(2,:)
                WRITE(9,*) endPos(3,:)
                WRITE(10,*) endVel(1,:)*vunit
                WRITE(11,*) endVel(2,:)*vunit
                WRITE(12,*) endVel(3,:)*vunit
                WRITE(13,*) origDists(:)
                WRITE(14,*) vmags(:)*vunit
                WRITE(15,*) closest(:,:) ! Entries 1, 3, 5, etc in each line are the IDs of the nearest objects. Entries 2, 4, 6, etc are the distances to those objects
                WRITE(16,*) Eerr, t*tunit
                EXIT
            END IF

            IF (t >= tLim) THEN ! This takes priority over the coundition (tcount >= tWrite) in the event that both are satisfied

                DO i=1,n
                    ! Getting distances from origin and velocity magnitudes
                    dist2 = 0._dp ! Need to reset to zero for each new i
                    vmag2 = 0._dp ! Need to reset to zero for each new i
                    DO k=1,ndims
                        dist2 = dist2 + (endPos(k,i)**2)
                        vmag2 = vmag2 + (endVel(k,i)**2)
                    END DO
                    origDists(i) = SQRT(dist2)
                    vmags(i) = SQRT(vmag2)
                END DO

                DO i=1,n
                    ! Finding the closest other object to each object
                    closestDist = 99999._dp ! Make this arbitrarily large at the start of each loop
                    DO j=1,n
                        IF (i==j) CYCLE ! Avoids having the closest object to a given object be itself
                        dist2 = 0._dp ! Reset for each j
                        DO k=1,ndims
                            dist2 = dist2 + ((endPos(k,i) - endPos(k,j))**2)
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
                WRITE(7,*) endPos(1,:)
                WRITE(8,*) endPos(2,:)
                WRITE(9,*) endPos(3,:)
                WRITE(10,*) endVel(1,:)*vunit
                WRITE(11,*) endVel(2,:)*vunit
                WRITE(12,*) endVel(3,:)*vunit
                WRITE(13,*) origDists(:)
                WRITE(14,*) vmags(:)*vunit
                WRITE(15,*) closest(:,:) ! Entries 1, 3, 5, etc in each line are the IDs of the nearest objects. Entries 2, 4, 6, etc are the distances to those objects
                WRITE(16,*) Eerr, t*tunit
                EXIT

            ELSE IF (tcount >= tWrite) THEN

                DO i=1,n
                    ! Getting distances from origin and velocity magnitudes
                    dist2 = 0._dp ! Need to reset to zero for each new i
                    vmag2 = 0._dp ! Need to reset to zero for each new i
                    DO k=1,ndims
                        dist2 = dist2 + (endPos(k,i)**2)
                        vmag2 = vmag2 + (endVel(k,i)**2)
                    END DO
                    origDists(i) = SQRT(dist2)
                    vmags(i) = SQRT(vmag2)
                END DO

                DO i=1,n
                    ! Finding the closest other object to each object
                    closestDist = 99999._dp ! Make this arbitrarily large at the start of each loop
                    DO j=1,n
                        IF (i==j) CYCLE ! Avoids having the closest object to a given object be itself
                        dist2 = 0._dp ! Reset for each j
                        DO k=1,ndims
                            dist2 = dist2 + ((endPos(k,i) - endPos(k,j))**2)
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
                WRITE(7,*) endPos(1,:)
                WRITE(8,*) endPos(2,:)
                WRITE(9,*) endPos(3,:)
                WRITE(10,*) endVel(1,:)*vunit
                WRITE(11,*) endVel(2,:)*vunit
                WRITE(12,*) endVel(3,:)*vunit
                WRITE(13,*) origDists(:)
                WRITE(14,*) vmags(:)*vunit
                WRITE(15,*) closest(:,:) ! Entries 1, 3, 5, etc in each line are the IDs of the nearest objects. Entries 2, 4, 6, etc are the distances to those objects
                WRITE(16,*) Eerr, t*tunit
                tcount = tcount - tWrite ! Not just resetting tcount to zero because this could cause the file outputs to drift from the intended times over the course
                                         ! of the simulation. Subtracting twrite keeps the file outputs at times as close to multiples of twrite as possible
            END IF

        END IF
    
    END DO
    
    RETURN
    
END SUBROUTINE n_body_simulation

! Subroutine that performs a single timestep of a given length for some or all of the objects
SUBROUTINE timestep(n, ndims, m, pos, vel, acc, jrk, snp, crk, endPos, endVel, dt, dtBase, doThisObject, timestepMultipliers, timeSinceLast)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    INTEGER, INTENT(IN) :: n, ndims, timestepMultipliers(1:n)
    LOGICAL, INTENT(IN) :: doThisObject(1:n)
    REAL(dp), INTENT(INOUT) :: pos(1:ndims, 1:n), vel(1:ndims, 1:n), acc(1:ndims, 1:n), jrk(1:ndims, 1:n), snp(1:ndims, 1:n), crk(1:ndims, 1:n)
    REAL(dp), INTENT(INOUT) :: endPos(1:ndims, 1:n), endVel(1:ndims, 1:n)
    REAL(dp), INTENT(IN) :: dt, dtBase, m(1:n), timeSinceLast(1:n)
    INTEGER :: i, j, k
    REAL(dp) :: interPos(1:ndims, 1:n), interVel(1:ndims, 1:n), dtFrac, dr(1:ndims), dv(1:ndims), dt2, dt3, dt4, dt5, recipdt2, recipdt3, recip6, recip24, recip120
    REAL(dp) :: vdotr, dist, dist2, dist3, dist5, accPred(1:ndims, 1:n), jrkPred(1:ndims, 1:n)

    dt2 = dt**2
    dt3 = dt**3
    dt4 = dt**4
    dt5 = dt**5
    recipdt2 = 1._dp / dt2
    recipdt3 = 1._dp / dt3
    recip6 = 1._dp/6._dp
    recip24 = 1._dp/24._dp
    recip120 = 1._dp/120._dp

    DO i=1,n ! Interpolation of positions and velocities - needs to be done BEFORE the rest of the subroutine (other than setting the values of constants)
        ! First order approximation - assuming constant velocity throughout a long timestep when interpolating positions, and assuming constant acceleration when
        ! interpolating velocities. Neither of these are true in reality
        IF (doThisObject(i) .eqv. .true.) THEN
            ! These objects dont need to have their positions or velocities interpolated
            interPos(:,i) = pos(:,i)
            interVel(:,i) = vel(:,i)
        ELSE
            IF ((dtBase * timestepMultipliers(i)) < dt) THEN
                ! Short-timestep objects - no need to interpolate anything as the current timestep does not begin during one of their timesteps
                ! Just set interPos and interVel to pos and vel respectively. Due to the timestep lengths being handled in descending order, these objects will not
                ! have had their positions updated yet and so pos and vel store the start-of-timestep positions and velocities
                interPos(:,i) = pos(:,i)
                interVel(:,i) = vel(:,i)
            ELSE
                ! Long-timestep objects, need to interpolate
                dtFrac = timeSinceLast(i) / (dtBase * timestepMultipliers(i))
                dr(:) = endPos(:,i) - pos(:,i) ! Here dr and dv are the end-of-timestep values minus the start-of-timestep values, rather than the difference between
                                               ! the values for a pair of objects
                dv(:) = endVel(:,i) - vel(:,i)
                interPos(:,i) = pos(:,i) + (dr(:) * dtFrac)
                interVel(:,i) = vel(:,i) + (dv(:) * dtFrac)
            END IF
        END IF
    END DO
            
    DO i=1,n
        IF (doThisObject(i) .eqv. .true.) THEN ! Only calculate new derivatives of the motion for certain objects
            acc(:,i) = 0._dp
            jrk(:,i) = 0._dp
            
            DO j=1,n
                IF (i==j) CYCLE
                dist2 = 0._dp
                vdotr = 0._dp
                dr(:) = interPos(:,i) - interPos(:,j)
                dv(:) = interVel(:,i) - interVel(:,j)
                DO k=1,ndims
                    vdotr = vdotr + (dv(k) * dr(k))
                    dist2 = dist2 + (dr(k)**2)
                END DO
                dist = SQRT(dist2)
                dist3 = dist**3
                dist5 = dist**5
                acc(:,i) = acc(:,i) - ((m(j) * dr(:)) / dist3)
                jrk(:,i) = jrk(:,i) + (m(j) * ((dv(:) / dist3) - ((3._dp * vdotr * dr(:)) / dist5)))
            END DO
        END IF
    END DO

    ! Predicted end-of-timestep positions and velocities (before snap and crackle are accounted for)
    ! Note that the predicted positions and velocities do NOT need to be in separate arrays to the pre-timestep positions and celocities, as the
    ! pre-timestep positions and velocities are not needed for the remainder of the timestep after the predicted positions and velocities have been calculated
    DO i=1,n
        IF (doThisObject(i) .eqv. .true.) THEN
            interPos(:,i) = interPos(:,i) + (interVel(:,i) * dt) + (0.5_dp * acc(:,i) * dt2) + (recip6 * jrk(:,i) * dt3)
            interVel(:,i) = interVel(:,i) + (acc(:,i) * dt) + (0.5_dp * jrk(:,i) * dt2) ! Velocity MUST be updated AFTER position, as updating the position requires the
                                                                                        ! pre-timestep velocity
        END IF
    END DO

    ! Need to re-interpolate to get end-of-timestep positions and velocities for the other objects
    DO i=1,n
        ! No need to reassign interPos or interVel for the objects that are actually having their positions and velocities updated this timestep
        IF (doThisObject(i) .eqv. .false.) THEN
            IF ((dtBase * timestepMultipliers(i)) < dt) THEN
                CONTINUE ! THIS IS NOT CORRECT RIGHT NOW, UNSURE HOW I WANT TO APPROACH THIS AT THE MOMENT, WILL DISCUSS WITH SIMON - USING CONTINUE AS A PLACEHOLDER
                         ! Currently this leaves interPos and interVel unchanged for the short-timestep objects, which is more than likely not correct
            ELSE
                dtFrac = (timeSinceLast(i) + dt) / (dtBase * timestepMultipliers(i))
                dr(:) = endPos(:,i) - pos(:,i)                              
                dv(:) = endVel(:,i) - vel(:,i)
                interPos(:,i) = pos(:,i) + (dr(:) * dtFrac) ! It is correct to use pos and vel on the right-hand sides of these two lines because of the way dtFrac has
                                                            ! been defined
                interVel(:,i) = vel(:,i) + (dv(:) * dtFrac)
            END IF
        END IF
    END DO

    ! Predicted accelerations and jerks at the predicted positions and velocities
    DO i=1,n
        IF (doThisObject(i) .eqv. .true.) THEN
            accPred(:,i) = 0._dp
            jrkPred(:,i) = 0._dp
            
            DO j=1,n
                IF (i==j) CYCLE
                dist2 = 0._dp
                vdotr = 0._dp
                dr(:) = interPos(:,i) - interPos(:,j)
                dv(:) = interVel(:,i) - interVel(:,j)
                DO k=1,ndims
                    vdotr = vdotr + (dv(k) * dr(k))
                    dist2 = dist2 + (dr(k)**2)
                END DO
                dist = SQRT(dist2)
                dist3 = dist**3
                dist5 = dist**5
                accPred(:,i) = accPred(:,i) - ((m(j) * dr(:)) / dist3)
                jrkPred(:,i) = jrkPred(:,i) + (m(j) * ((dv(:) / dist3) - ((3._dp * vdotr * dr(:)) / dist5)))
            END DO
        END IF
    END DO

    ! Snaps and crackles
    DO i=1,n
        IF (doThisObject(i) .eqv. .true.) THEN
            snp(:,i) = 0._dp
            crk(:,i) = 0._dp
            
            DO j=1,n
                IF (i==j) CYCLE

                ! No idea if the dt's (and powers of dt) in the following 2 lines need to be the dt for that specific object or if I can juse use the value of dt
                ! for the current object(s) being considered - I will ask Simon about this when I next meet with him
                snp(:,i) = snp(:,i) - (recipdt2 * ((6._dp * (acc(:,j) - accPred(:,j))) + (dt * ((4._dp * jrk(:,j)) + (2._dp * jrkPred(:,j))))))
                crk(:,i) = crk(:,i) + (recipdt3 * ((12._dp * (acc(:,j) - accPred(:,j))) + (6._dp * dt * (jrk(:,j) + jrkPred(:,j)))))
            END DO
        END IF
    END DO

    ! Position and velocity corrections. This COULD be included in the previous DO loop after cycling through all j's for a given i as the snap/crackle
    ! calculations in the actual loop do not involve the predicted positions/velocities of the relevant objects (the contributions of the predicted positions/velocities
    ! are contained within accPred and jrkPred, which were calculated beforehand). Updating the positions and velocities during the same DO loop as
    ! the snaps and crackles are esimated shouldn't affect the snap and crackle estimations, and won't have any effect on the calculated positions and velocities
    ! as long as they are updated AFTER cycling through all j's for a given i. I have decided against putting them in the same DO loop
    ! in order to make the code easier to follow
    DO i=1,n
        IF (doThisObject(i) .eqv. .true.) THEN
            interPos(:,i) = interPos(:,i) + (recip24 * snp(:,i) * dt4) + (recip120 * crk(:,i) * dt5)
            interVel(:,i) = interVel(:,i) + (recip6 * snp(:,i) * dt3) + (recip24 * crk(:,i) * dt4)

            ! Also setting endPos and endVel to be equal to interPos and interVel respectively for the objects that are having their positions and velocities
            ! updated on this specific timestep. This allows both the start- and end-of-timestep positions and velocities to be stored for interpolation purposes
            endPos(:,i) = interPos(:,i)
            endVel(:,i) = interVel(:,i)
        END IF
    END DO

    RETURN

END SUBROUTINE timestep

                


    
! Gives the N-body velocity unit in km/s and the time unit in Myr
SUBROUTINE n_body_conversion(munit, runit, tunit, vunit)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    REAL(dp), INTENT(IN) :: munit, runit
    REAL(dp), INTENT(OUT) :: tunit, vunit
    REAL(dp) :: G, year

    G = 6.67e-11_dp ! Gravitational constant in SI units
    year = 365.24_dp * 86400._dp ! 86400 seconds in a day, 365.24 days in a year
    tunit = SQRT((runit**3) / (G * munit)) / (1.e6_dp * year) ! N-body time unit in Myr (should be about 14.9 Myr)
    vunit = runit / (1.e9_dp * tunit * year) ! N-body velocity unit in km/s. 1 vunit should be (1 pc)/(14.9 Myr), or about (6.56e-2 km/s).
                                             ! The factor of runit/tunit converts from N-body velocity units to m/Myr. The factor of 1/year then converts to m/megasecond.
                                             ! Finally, the factor of 1.e9_dp on the denominator converts to km/s

    RETURN
    
END SUBROUTINE n_body_conversion

! Plummer sphere (3 dimensions only). Method from Aarseth, Hénon and Wielen (1974)
SUBROUTINE setup_plummer(n, pos, vel, m)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, seed
    REAL(dp), INTENT(OUT) :: pos(1:3, 1:n), vel(1:3, 1:n), m(1:n)
    REAL(dp) :: ranNum, ranNum2, twopi, twothirds, rad, sqrtr2z2, vesc, root2, q2, gq, vmag
    DOUBLE PRECISION :: ran3


    seed = 242134               ! Seed for the pseudo-random number generator
    twopi = 8._dp * ATAN(1._dp) ! 2 * pi
    ranNum = 0._dp              ! Random number in a uniform distribution between 0 and 1
    ranNum2 = 0._dp             ! Another random number in a uniform distribution between 0 and 1, used when a new random number needs to be generated while the value stored
                                ! in ranNum is also still needed
    twothirds = 2._dp/3._dp     ! 2/3
    rad = 0._dp                 ! Distance of a given object from the origin. Not an array as the generation of initial conditions for a given object only needs the value
                                ! of rad for the object itself
    sqrtr2z2 = 0._dp            ! SQRT(r^2 - z^2). Also stores SQRT(v^2 - v_z^2) as this is value is only needed for a given object once SQRT(r^2 - z^2) is no longer needed
    vesc = 0._dp                ! Escape velocity at a given object's position, which represents the maximum 'allowed' velocity at that position. Not an array for the same
                                ! reasons as rad
    vmag = 0._dp                ! Velocity magnitude of a given object. Not an array for the same reasons as rad
    q2 = 0._dp                  ! q^2, where q = vmag/vesc
    gq = 0._dp                  ! g(q) = q^2 - (1-q^2)^7/2
    root2 = SQRT(2._dp)         ! Square root of 2

    m = 1._dp ! All objects will be 1 solar mass (for now)
    DO i=1,n
        ranNum = ran3(seed)
        rad = SQRT(1._dp/((ranNum**(-1._dp*twothirds)) - 1._dp))
        ranNum = ran3(seed)
        pos(3,i) = rad * (1._dp - (2._dp*ranNum))
        ranNum = ran3(seed) ! Need a different random number for the x and y coordinates to the one used for the z coordinates
        sqrtr2z2 = SQRT((rad**2) - (pos(3,i)**2))
        pos(1,i) = COS(twopi * ranNum) * sqrtr2z2
        pos(2,i) = SIN(twopi * ranNum) * sqrtr2z2
        vesc = root2 * ((1._dp + (rad**2))**(-0.25_dp))
        DO
            ! This will run continuously until a certain condition is met (this avoids having to use GOTO)
            ranNum = ran3(seed)
            q2 = ranNum**2
            gq = q2 * ((1._dp - q2)**3.5_dp)
            ranNum2 = ran3(seed)
            IF ((0.1_dp * ranNum2) < gq) THEN
                EXIT ! This ends the current DO loop (which would otherwise run infinitely)
            END IF
        END DO
        vmag = ranNum * vesc ! Note this uses ranNum, NOT ranNum2
        ranNum = ran3(seed)
        vel(3,i) = vmag * (1._dp - (2._dp*ranNum))
        ranNum = ran3(seed) ! Need a different random number for the x and y velocity components to the one used for the z velocity component
        sqrtr2z2 = SQRT((vmag**2) - (vel(3,i)**2)) ! Reusing this variable name as the value it was set to earlier in the loop is no longer needed at this point
        vel(1,i) = COS(twopi * ranNum) * sqrtr2z2
        vel(2,i) = SIN(twopi * ranNum) * sqrtr2z2
    END DO

    RETURN

END SUBROUTINE setup_plummer

SUBROUTINE setup_fig_eight(pos, vel, m)
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    REAL(dp), INTENT(OUT) :: pos(1:3, 1:3), vel(1:3, 1:3), m(1:3)

    m = 1._dp

    pos(1,1) = 0.9700436_dp
    pos(2,1) = -0.24308753_dp
    pos(1,2) = -1._dp * pos(1,1)
    pos(2,2) = -1._dp * pos(2,1)
    pos(1,3) = 0._dp ! Should already be the case on calling this subroutine, but this line is included on the off chance it isn't the case
    pos(2,3) = 0._dp ! See above comment
    pos(3,:) = 0._dp ! See comment two lines above

    vel(1,1) = 0.466203685_dp
    vel(2,1) = 0.43236573_dp
    vel(1,2) = vel(1,1)
    vel(2,2) = vel(2,1)
    vel(1,3) = -2._dp * vel(1,1)
    vel(2,3) = -2._dp * vel(2,1)
    vel(3,:) = 0._dp ! Should already be the case on calling this subroutine, but this line is included on the off chance it isn't the case

    RETURN

END SUBROUTINE setup_fig_eight

SUBROUTINE closest_timestep(maxPower, dtBase, dtAllowedMultipliers, dtOpt, timestepMultiplier, numPossibleTimestepLengths)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    INTEGER :: i
    INTEGER, INTENT(IN) :: numPossibleTimestepLengths
    INTEGER, INTENT(IN) :: maxPower, dtAllowedMultipliers(1:numPossibleTimestepLengths)
    INTEGER, INTENT(OUT) :: timestepMultiplier
    REAL(dp), INTENT(IN) :: dtBase, dtOpt
    REAL(dp) :: dtRatio

    dtRatio = dtOpt / dtBase

    IF (dtRatio < dtAllowedMultipliers(1)) THEN
        timestepMultiplier = 2 ** (-1 * maxPower)
    ELSE IF (dtRatio > dtAllowedMultipliers(numPossibleTimestepLengths)) THEN
        timestepMultiplier = 2 ** maxPower
    ELSE
        ! The following assumes that dtAllowedMultipliers is sorted in ascending order (which it should be)
        DO i=2,numPossibleTimestepLengths ! Start from i=2 because the case for i=1 has already been checked earlier in the IF statement
            IF (dtRatio < dtAllowedMultipliers(i)) THEN
                timestepMultiplier = dtAllowedMultipliers(i-1) ! Rounding down to the nearest allowed multiplier
                EXIT
            END IF
        END DO
    END IF

    RETURN
    
END SUBROUTINE closest_timestep









! Pseudo-random number generator (NOT MY OWN CODE)
DOUBLE PRECISION FUNCTION RAN3(IDUM)                                   
  IMPLICIT NONE
  INTEGER idum,iff,k,i,ii,inext,inextp
  DOUBLE PRECISION FAC,mbig,mseed,mz,ma,mj,mk
  PARAMETER (MBIG=4000000., MSEED=1618033., MZ=0., FAC=1./MBIG)
  DIMENSION MA(55)
  SAVE
  DATA IFF/0/
  IF(IDUM<0 .OR. IFF==0)THEN
     IFF=1
     MJ=MSEED-IABS(IDUM)
     MJ=MOD(MJ,MBIG)
     MA(55)=MJ
     MK=1
     DO 100 I=1,54
        II=MOD(21*I,55)
        MA(II)=MK
        MK=MJ-MK
        IF(MK.LT.MZ)MK=MK+MBIG
        MJ=MA(II)
100     CONTINUE
        DO 200 K=1,4
           DO 300 I=1,55
              MA(I)=MA(I)-MA(1+MOD(I+30,55))
              IF(MA(I)<MZ)MA(I)=MA(I)+MBIG
300           CONTINUE
200           CONTINUE
              INEXT=0
              INEXTP=31
              IDUM=1
           ENDIF
           INEXT=INEXT+1
           IF(INEXT==56)INEXT=1
           INEXTP=INEXTP+1
           IF(INEXTP==56)INEXTP=1
           MJ=MA(INEXT)-MA(INEXTP)
           IF(MJ<MZ)MJ=MJ+MBIG
           MA(INEXT)=MJ
           RAN3=MJ*FAC
           RETURN
         END FUNCTION ran3