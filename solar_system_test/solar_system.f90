PROGRAM solar_system
    IMPLICIT NONE
    INTEGER :: n, i, j, k, ndims
    DOUBLE PRECISION :: t, tcount, twrite, tlim, dt, dt2, dt3, dt4, dt5, eta, G, pi, year, GReal, mReal, lReal, tReal, vReal, initPhase, mTot, E0, Ek, Ep, E, Eerr
    DOUBLE PRECISION :: vmag2, dist, dist2, dist3, dist5, amag2, vdotr, recip6, recip24, recip120, recipdt2, recipdt3
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m, dtOpt, acRad, acVel, comPos, comVel, dr, dv, sunDist
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pos, vel, acc, jrk, snp, crk, accPred, jrkPred

    OPEN(7, file='sun_data.txt', status='replace')
    OPEN(8, file='mercury_data.txt', status='replace')
    OPEN(9, file='venus_data.txt', status='replace')
    OPEN(10, file='earth_data.txt', status='replace')
    OPEN(11, file='mars_data.txt', status='replace')
    OPEN(12, file='jupiter_data.txt', status='replace')
    OPEN(13, file='saturn_data.txt', status='replace')
    OPEN(14, file='uranus_data.txt', status='replace')
    OPEN(15, file='neptune_data.txt', status='replace')
    OPEN(16, file='energy_error.txt', status='replace')

    n = 9               ! Number of objects
    ndims = 3           ! Number of dimensions

    ALLOCATE(m(1:n))
    ALLOCATE(pos(1:ndims, 1:n))
    ALLOCATE(vel(1:ndims, 1:n))
    ALLOCATE(acc(1:ndims, 1:n))
    ALLOCATE(jrk(1:ndims, 1:n))
    ALLOCATE(snp(1:ndims, 1:n))
    ALLOCATE(crk(1:ndims, 1:n))
    ALLOCATE(dtOpt(1:n))
    ALLOCATE(acRad(1:n))
    ALLOCATE(acVel(1:n))
    ALLOCATE(comPos(1:ndims))
    ALLOCATE(comVel(1:ndims))
    ALLOCATE(accPred(1:ndims, 1:n))
    ALLOCATE(jrkPred(1:ndims, 1:n))
    ALLOCATE(dr(1:ndims))
    ALLOCATE(dv(1:ndims))
    ALLOCATE(sunDist(1:n))

    t = 0.              ! Time
    tcount = 0.         ! Tracks the time between file writes
    twrite = 0.         ! Time between file writes, to be defined in terms of tReal later
    tlim = 0.           ! Total simulated time, to be defined in terms of tReal later
    dt = 0.             ! Timestep length, will vary each timestep according to the acceleration of each object, and will be made non-zero before 
                        ! the first time it is used in calculations
    dt2 = 0.            ! dt^2
    dt3 = 0.            ! dt^3
    dt4 = 0.            ! dt^4
    dt5 = 0.            ! dt^5
    recipdt2 = 0.       ! dt^-2
    recipdt3 = 0.       ! dt^-3
    dtOpt = 0.          ! Optimal timestep lengths for each object, which depends on the current acceleration of that object
    eta = 1e-10         ! Used in calculating the optimal timestep lengths for each object, all of the timestep lengths are scaled by SQRT(eta)
    G = 1.              ! Gravitational constant, which is set to 1 in N-body units
    pi = 4. * ATAN(1.)  ! Pi to maximum computational precision
    year = 3.15e7       ! Number of seconds in a year
    acRad = 0.          ! Actual semi-major axes for each planet. Position 1 in this array is 0 as it corresponds to the Sun, this could be completely omitted 
                        ! but has been included for consistency with the sizes of other arrays
    acVel = 0.          ! Same as acRad but for orbital velocities
    pos = 0.            ! Positions
    vel = 0.            ! Velocities
    acc = 0.            ! Accelerations
    jrk = 0.            ! Jerks
    snp = 0.            ! Snaps
    crk = 0.            ! Crackles
    m = 0.              ! Masses
    GReal = 6.67e-11    ! SI value of the gravitational constant, used for scaling to N-body units
    mReal = 1.99e30     ! Solar mass, used for scaling to N-body units
    lReal = 1.50e11     ! 1 AU, used for scaling to N-body units
    tReal = 0.          ! N-body time unit, will be calculated later
    vReal = 0.          ! N-body velocity unit, will be calculated later
    initPhase = 0.      ! Initial position of a given planet in its orbit. Does not need to be an array as the value can be re-randomised for each planet after the initial
                        ! position and velocity of a given planet are adjusted according to the value of initPhase
    mTot = 0.           ! Total mass, will be calculated when calculating the position and velocity of the origin of the COM frame
    comPos = 0.         ! COM origin position, will be calculated at the same time as mTot
    comVel = 0.         ! COM origin velocity, will be calculated at the same time as mTot
    E0 = 0.             ! Initial total energy, will be calculated after moving to COM frame
    Ek = 0.             ! Tracks the total kinetic energy
    Ep = 0.             ! Tracks the total potential energy
    E = 0.              ! Tracks the total energy
    Eerr = 0.           ! Fractional difference between E0 and E
    vmag2 = 0.          ! Velocity magnitude squared (vel.vel)
    dist = 0.           ! Temporarily stores the distance between two objects
    dist2 = 0.          ! Temporarily stores dist^2
    dist3 = 0.          ! Temporarily stores dist^3
    dist5 = 0.          ! Temporarily stores dist^5. Note that dist^4 is not used in the simulation at all so there is no dist4 variable
    accPred = 0.        ! Acceleration at initially predicted end-of-timestep positions and velocities (before accounting for snap and crackle)
    jrkPred = 0.        ! Same as accPred but for the jerk
    amag2 = 0.          ! Acceleration magnitude squared (acc.acc)
    dr = 0.             ! Disrance vector between two objects (pos_i - pos_j)
    dv = 0.             ! Difference in velocity between two objects (vel_i - vel_j)
    vdotr = 0.          ! Dot product of dv and dr
    recip6 = 1./6.      ! Multiply by this instead of dividing by 6
    recip24 = 1./24.    ! Multiply by this instead of dividing by 24
    recip120 = 1./120.  ! Multiply by this instead of dividing by 120
    sunDist = 0.        ! Distance to the Sun for each object. The value of sunDist(1) is trivially zero throughout as this is the distance from the Sun to itself
                        ! This value is included for the sake of consistency in which object corresponds to which index in various arrays

    ! N-body scaling
    tReal = ((lReal**3)/(mReal*GReal))**0.5
    vReal = lReal / tReal

    ! Time limits
    twrite = (0.01 * year) / tReal ! Write to files 100 times each SI year
    tlim = (1000. * year) / tReal  ! Run for 1000 SI years in total

    ! Initial conditions (most information is taken from the NASA planetary fact sheet)

    ! Sun
    m(1) = 1.99e30 / mReal

    ! Mercury
    m(2) = 0.33e24 / mReal
    acRad(2) = 5.79e10 / lReal
    acVel(2) = 4.74e4 / vReal

    ! Venus
    m(3) = 4.87e24 / mReal
    acRad(3) = 1.08e11 / lReal
    acVel(3) = 3.50e4 / vReal

    ! Earth
    m(4) = 5.97e24 / mReal
    acRad(4) = 1.50e11 / lReal
    acVel(4) = 2.98e4 / vReal

    ! Mars
    m(5) = 0.64e24 / mReal
    acRad(5) = 2.20e11 / lReal
    acVel(5) = 2.41e4 / vReal

    ! Jupiter
    m(6) = 1898e24 / mReal
    acRad(6) = 7.79e11 / lReal
    acVel(6) = 1.31e4 / vReal

    ! Saturn
    m(7) = 568e24 / mReal
    acRad(7) = 1.43e12 / lReal
    acVel(7) = 0.97e4 / vReal

    ! Uranus
    m(8) = 86.8e24 / mReal
    acRad(8) = 2.87e12 / lReal
    acVel(8) = 0.68e4 / vReal

    ! Neptune
    m(9) = 102e24 / mReal
    acRad(9) = 4.52e12 / lReal
    acVel(9) = 0.54e4 / vReal

    ! Randomisation
    DO i=2,n ! Starts at n=2 as the initial position and velocity of the Sun do not need to be randomised (it starts at rest at the origin before moving to the COM frame)
        CALL RANDOM_NUMBER(initPhase)
        initPhase = 2. * pi * initPhase         ! Initial angle of the planet's position in its orbit above y=0
        pos(1,i) = acRad(i) * COS(initPhase)    ! rcos(theta). This gives the correct sign on the x-position for all orbital phase angles
        pos(2,i) = acRad(i) * SIN(initPhase)    ! rsin(theta). This gives the correct sign on the y-position for all orbital phase angles
        vel(1,i) = -acVel(i) * SIN(initPhase)   ! -vsin(theta). Orbits are anti-clockwise so the x-component of velocity is in the opposite direction to the x-component
                                                ! of the position relative to the Sun, and the equation used gives the correct sign on the x-velocity for all 
                                                ! orbital phase angles
        vel(2,i) = acVel(i) * COS(initPhase)    ! vcos(theta). This gives the correct sign on the y-velocity for all orbital phase angles
    END DO

    ! COM frame calculation
    DO i=1,n
        mTot = mTot + m(i)
        DO k=1,ndims
            comPos(k) = comPos(k) + (m(i) * pos(k,i))
            comVel(k) = comVel(k) + (m(i) * vel(k,i))
        END DO
    END DO
    comPos(:) = comPos(:) / mTot
    comVel(:) = comVel(:) / mTot

    ! Moving into COM frame
    DO i=1,n
        DO k=1,ndims
            pos(k,i) = pos(k,i) - comPos(k)
            vel(k,i) = vel(k,i) - comVel(k)
        END DO
    END DO

    ! Initial total energy. Note Ek and Ep do not need to be reset to zero here as they will not have been changed from zero as of this point in the simulation
    DO i=1,n
        vmag2 = 0. ! Needs to be reset to zero for each iteration
        DO k=1,ndims
            vmag2 = vmag2 + (vel(k,i)**2)
        END DO
        Ek = Ek + (0.5 * m(i) * vmag2)
        
        DO j=i+1,n ! These boundaries avoid double-counting pairs of objects in the potential energy calculations
            dist2 = 0. ! Needs to be reset to zero for each iteration
            DO k=1,ndims
                dist2 = dist2 + ((pos(k,i)-pos(k,j))**2)
            END DO
            dist = SQRT(dist2) ! Do not need to reset dist to zero each iteration here as it is not calculated directly through incrementation
            Ep = Ep - ((G * m(i) * m(j)) / dist)
        END DO
    END DO
    
    E0 = Ek + Ep

    ! Main Loop
    DO
        ! Resetting values that are calculated through incrementation over all objects
        Ek = 0.
        Ep = 0.
        acc = 0.
        jrk = 0.
        snp = 0.
        crk = 0.
        accPred = 0.
        jrkPred = 0.

        ! Accelerations, jerks and adaptive timestep length
        DO i=1,n
            amag2 = 0. ! Needs to be reset to zero for each iteration
            DO j=1,n
                IF (i==j) CYCLE ! Finding the contribution of particle i to the acceleration and jerk of particle i is redundant
                dist2 = 0.
                vdotr = 0.
                DO k=1,ndims
                    dr(k) = pos(k,i) - pos(k,j)
                    dv(k) = vel(k,i) - vel(k,j)
                    vdotr = vdotr + (dv(k) * dr(k))
                    dist2 = dist2 + (dr(k)**2)
                END DO
                dist = SQRT(dist2)
                dist3 = dist**3
                dist5 = dist**5
                DO k=1,ndims
                    acc(k,i) = acc(k,i) - ((G * m(j) * dr(k)) / dist3)
                    amag2 = amag2 + (acc(k,i)**2)
                    jrk(k,i) = jrk(k,i) + (G * m(j) * ((dv(k) / dist3) - ((3. * vdotr * dr(k)) / dist5)))
                END DO
            END DO
            dtOpt(i) = SQRT(eta / amag2)
        END DO
        dt = MINVAL(dtOpt(:))
        dt2 = dt**2
        dt3 = dt**3
        dt4 = dt**4
        dt5 = dt**5
        recipdt2 = 1. / dt2
        recipdt3 = 1. / dt3

        ! Predicted end-of-timestep positions and velocities (before snap and crackle are accounted for)
        DO i=1,n
            ! Predicted positions/velocities do not to be in separate arrays to the pre-timestep positions/velocities.
            ! This is because the pre-timestep positions/velocities are not needed for the remainder of the timestep after the predicted
            ! positions/velocities have been calculated
            pos(:,i) = pos(:,i) + (vel(:,i) * dt) + (0.5 * acc(:,i) * dt2) + (recip6 * jrk(:,i) * dt3)
            vel(:,i) = vel(:,i) + (acc(:,i) * dt) + (0.5 * jrk(:,i) * dt2) ! Velocity MUST be updated AFTER position as updating the position requires
                                                                           ! the pre-timestep velocity
        END DO

        ! Predicted accelerations and jerks at the predicted positions and velocities
        DO i=1,n
            ! Predicted accelerations/jerks do need separate arrays to the accelerations/jerks used in predicting the end-of-timestep positions/velocities
            ! because both pairs are needed to estimate the snaps/crackles
            DO j=1,n
                IF(i==j) CYCLE
                dist2 = 0.
                vdotr = 0.
                DO k=1,ndims
                    dr(k) = pos(k,i) - pos(k,j)
                    dv(k) = vel(k,i) - vel(k,j)
                    vdotr = vdotr + (dv(k) * dr(k))
                    dist2 = dist2 + (dr(k)**2)
                END DO
                dist = SQRT(dist2)
                dist3 = dist**3
                dist5 = dist**5
                DO k=1,ndims
                    accPred(k,i) = accPred(k,i) - ((G * m(j) * dr(k)) / dist3)
                    jrkPred(k,i) = jrkPred(k,i) + (G * m(j) * ((dv(k) / dist3) - ((3. * vdotr * dr(k)) / dist5)))
                END DO
            END DO
        END DO

        ! Estimating the snap and crackle of each object. Separate DO loop to the calculation of accPred and jrkPred because the estimation of the snap and crackle
        ! for a given object requires knowledge of accPred and jrkPred for all other objects
        DO i=1,n
            DO j=1,n
                IF (i==j) CYCLE
                DO k=1,ndims
                    ! Note that the subtraction in the snap estimation is outside the parentheses so the jerk term is added rather than subtracted
                    snp(k,i) = snp(k,i) - (recipdt2 * ((6. * (acc(k,j) - accPred(k,j))) + (dt * ((4. * jrk(k,j)) + (2. * jrkPred(k,j))))))
                    crk(k,i) = crk(k,i) + (recipdt3 * ((12. * (acc(k,j) - accPred(k,j))) + (6. * dt * (jrk(k,j) + jrkPred(k,j)))))
                END DO
            END DO
        END DO

        ! Position and velocity corrections. This COULD be including in the previous DO loop after cycling through all j's for a given i as the snap/crackle
        ! calculations in the actual loop do not involve the predicted positions/velocities of the relevant objects (the contributions of the predicted positions/velocities
        ! are contained within accPred and jrkPred, which were calculated beforehand). Updating the positions and velocities during the same DO loop as
        ! the snaps and crackles are esimated shouldn't affect the snap and crackle estimations, and won't have any effect on the calculated positions and velocities
        ! as long as they are updated AFTER cycling through all j's for a given i. I have decided against putting them in the same DO loop
        ! in order to make the code easier to follow
        DO i=1,n
            ! The contributions from orders up to the jerk are already included within the current values of pos and vel for a given object before they are updated here
            ! so only the contributions from the snap and crackle need to be accounted for here
            pos(:,i) = pos(:,i) + (recip24 * snp(:,i) * dt4) + (recip120 * crk(:,i) * dt5)
            vel(:,i) = vel(:,i) + (recip6 * snp(:,i) * dt3) + (recip24 * crk(:,i) * dt4)
        END DO

        ! Time incrementation

        t = t + dt
        tcount = tcount + dt

        ! File writes
        IF (t >= tlim) THEN ! This takes priority over the condition (tcount >= twrite) in the event that both are satisfied
            
            ! Energy error calculation. This will be re-used in the ELSE IF segment after this, but I don't want it to be calculated on every single timestep
            ! so increasing the number of lines of code by re-writing this block of code again later will decrease the runtime
            DO i=1,n
                vmag2 = 0.
                DO k=1,ndims
                    vmag2 = vmag2 + (vel(k,i)**2)
                END DO
                Ek = Ek + (0.5 * m(i) * vmag2)

                DO j=i+1,n
                    dist2 = 0.
                    DO k=1,ndims
                        dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
                    END DO
                    dist = SQRT(dist2)
                    Ep = Ep - ((G * m(i) * m(j)) / dist)
                END DO
            END DO
            E = Ek + Ep
            Eerr = ABS((E-E0)/E0)

            ! Distance to the Sun
            sunDist = 0.
            DO i=2,n ! i=1 is redundant as this would calculate the distance from the Sun to itself
                DO k=1,ndims
                    sunDist(i) = sunDist(i) + ((pos(k,i) - pos(k,1))**2) ! It doesn't matter which way around the subtraction is as the square is taken anyway
                END DO
                sunDist(i) = SQRT(sunDist(i))
            END DO

            ! File writes
            WRITE(7,*) pos(:,1)
            WRITE(8,*) pos(:,2), sunDist(2)
            WRITE(9,*) pos(:,3), sunDist(3)
            WRITE(10,*) pos(:,4), sunDist(4)
            WRITE(11,*) pos(:,5), sunDist(5)
            WRITE(12,*) pos(:,6), sunDist(6)
            WRITE(13,*) pos(:,7), sunDist(7)
            WRITE(14,*) pos(:,8), sunDist(8)
            WRITE(15,*) pos(:,9), sunDist(9)
            WRITE(16,*) Eerr, t
            EXIT

        ELSE IF (tcount >= twrite) THEN

            ! Energy error
            DO i=1,n
                vmag2 = 0.
                DO k=1,ndims
                    vmag2 = vmag2 + (vel(k,i)**2)
                END DO
                Ek = Ek + (0.5 * m(i) * vmag2)

                DO j=i+1,n
                    dist2 = 0.
                    DO k=1,ndims
                        dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
                    END DO
                    dist = SQRT(dist2)
                    Ep = Ep - ((G * m(i) * m(j)) / dist)
                END DO
            END DO
            E = Ek + Ep
            Eerr = ABS((E-E0)/E0)

            ! Distance to the Sun
            sunDist = 0.
            DO i=2,n
                DO k=1,ndims
                    sunDist(i) = sunDist(i) + ((pos(k,i) - pos(k,1))**2)
                END DO
                sunDist(i) = SQRT(sunDist(i))
            END DO

            ! File writes
            WRITE(7,*) pos(:,1)
            WRITE(8,*) pos(:,2), sunDist(2)
            WRITE(9,*) pos(:,3), sunDist(3)
            WRITE(10,*) pos(:,4), sunDist(4)
            WRITE(11,*) pos(:,5), sunDist(5)
            WRITE(12,*) pos(:,6), sunDist(6)
            WRITE(13,*) pos(:,7), sunDist(7)
            WRITE(14,*) pos(:,8), sunDist(8)
            WRITE(15,*) pos(:,9), sunDist(9)
            WRITE(16,*) Eerr, t

            tcount = tcount - twrite ! Not just resetting tcount to zero because this could cause the file outputs to drift from the intended times over the course
                                     ! of the simulation. Subtracting twrite keeps the file outputs at times as close to multiples of twrite as possible
        END IF

    END DO

END PROGRAM solar_system