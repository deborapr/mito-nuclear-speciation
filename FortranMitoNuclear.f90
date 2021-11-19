! program FortranMitoNuclear.f90

! f2py -c --fcompiler=gnu95 -m FortranMitoNuclear FortranMitoNuclear.f90

! Debora Princepe & Marcus A.M. de Aguiar  - 21/Dec/2019

! Uses hard walls boundary conditions
! Number of individuals grows until N(t) >= N
! sex separation in vector s()
! g(i,j) = nuclear gene j of individual i
! gm(i,j) = mitochondrial gene j o individual i

! let S increases up to S + nrmax to find surrugate mother
! if no new mother is available within S+nrmax, individual dies.

! module defining global variables
MODULE globals
INTEGER(1), SAVE :: g(5000,15000),gm(5000,5000),s(5000)
INTEGER, SAVE :: iseed(12),input_seed
INTEGER, SAVE :: ndim,nf1,nf2,nb,nbm,mnpm,ineighborg,ineighbor,ineighbort,nrmax,iitime,ntime
INTEGER, SAVE :: igt,igt_old,igtfull,igtfull_old,igtm,iibound
INTEGER, SAVE :: inisbit,nc,nct,nco,ineighborden
INTEGER, SAVE :: iini,ifin,jini,jfin
INTEGER, SAVE :: x(5000),y(5000),neig(5000),neigsp(5000),neigspt(5000),ispecies(5000,500),ispeciesm(5000,500)
INTEGER, SAVE :: worg(500,500,20),nworg(500,500)
INTEGER, SAVE :: ndq(5),idqx(5,10000),idqy(5,10000),ispv(500)
INTEGER, SAVE :: t(5000,5000),tp(5000,5000),p1(5000),p2(5000),ispidx(5000),ispidxold(5000)
INTEGER, SAVE :: ispvm(500),ispidxm(5000)
INTEGER, SAVE :: max_species
INTEGER, SAVE :: ext_time_old(500),ext_time(500),t_sp(500,500),t_sp_old(500,500)
INTEGER, SAVE :: size_ext(500),size_ext_old(500),ispv_old(500)
INTEGER, SAVE :: sister_sp(500),sister_sp_old(500)
REAL, SAVE :: rg,qmat,aux,rho0,radius,mut,diff,mutmit,fmax,fdelta,width
REAL, SAVE :: fitness(5000),distmitonuc(5000)
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module that reads and writes the populations
MODULE readwrite
INTEGER, SAVE :: iread
CHARACTER*25, SAVE :: name
CHARACTER*30, SAVE :: name1
CONTAINS
	    
SUBROUTINE READPOP
USE globals
! initialize or read population
IF(iread==1) THEN
    name = 'new.dat'
    name1 = 'pop-'//name
    OPEN(UNIT=10,FILE=name1,STATUS='old')
    READ(10,*) iitime,nct,nf1,nf2
    READ(10,*) amutd,amutmitd,diffd !dumb only
    READ(10,*) radiusd,rgd,nbd,nbmd,mnpmd !dumb only
    DO i=1,nct
        READ (10,901) x(i),y(i),(g(i,j),j=1,nb)
    END DO
    DO i=1,nct
        READ (10,902) (gm(i,j),j=1,nbm)
    END DO
    READ (10,903) (s(i),i=1,nct)
    CLOSE(10)
    name1 = 'time-'//name  !!write file with time matrix
    OPEN(UNIT=10,FILE=name1,STATUS='old')
    READ(10,*) nct
    DO i=1,nct-1
        READ(10,904) (t(i,jj),jj=i+1,nct)
    END DO
    t = t + transpose(t)
    CLOSE(10)
ELSE
    nct = nco
    DO i=1,nct   ! random distribution
        ictr = 0
        DO WHILE (ictr == 0)
            CALL RANDOM_NUMBER(aux)
            ii = int(aux*nf1)+1
            CALL RANDOM_NUMBER(aux)
            jj = int(aux*nf2)+1
            x(i) = ii
            y(i) = jj
            ictr = 1
        END DO
        CALL RANDOM_NUMBER(aux)
        s(i) = 1
        IF(aux < 0.5) s(i) = 0
    END DO
    iitime = 0
END IF
901 FORMAT(i4,1x,i4,1x,200000i1)
902 FORMAT(200000i1)
903 FORMAT(200000(i1,1x))
904 FORMAT(5000(i7,1x)) !time matrix
END SUBROUTINE READPOP

SUBROUTINE WRITEPOP
USE GLOBALS
name='new.dat'
name1 = 'pop-'//name

OPEN(UNIT=10,FILE=name1,STATUS='UNKNOWN',POSITION='REWIND')
WRITE(10,*) iitime,nct,nf1,nf2
WRITE(10,*) mut,mutmit,diff
WRITE(10,*) radius,rg,nb,nbm,mnpm
DO i=1,nct
    WRITE (10,901) x(i),y(i),(g(i,j),j=1,nb)
END DO
DO i=1,nct
    WRITE (10,902) (gm(i,j),j=1,nbm)
END DO
WRITE (10,903) (s(i),i=1,nct)
CLOSE(10)

name1 = 'time-'//name  !!write file with time matrix
OPEN(UNIT=9,FILE=name1,STATUS='unknown')
WRITE(9,*) nct
DO i=1,nct-1
    WRITE(9,904) (t(i,jj),jj=i+1,nct)
END DO
CLOSE(9)

901 FORMAT(i4,1x,i4,1x,200000i1)
902 FORMAT(200000i1)
903 FORMAT(200000(i1,1x))
904 FORMAT(5000(i7,1x))
END SUBROUTINE WRITEPOP


SUBROUTINE WRITEPHYLOGENY
USE GLOBALS
OPEN(unit=15,file='matrix-phy.dat',status='unknown')
WRITE(15,*) igt,iitime,nb
DO ii=1,igt
    write(15,902) (t(ispecies(ii,1),ispecies(jj,1)),jj=1,igt)
END DO
CLOSE(15)

OPEN(unit=15,file='matrix-phy-color.dat',status='unknown')
do ii=1,igt
    write(15,*) ispidx(ispecies(ii,1))
end do
CLOSE(15)

OPEN(unit=15,file='matrix-phy-genetics.dat',status='unknown')
do ii=1,igt
    WRITE (15,903) (g(ispecies(ii,1),jj),jj=1,nb)
end do
CLOSE(15)

OPEN(unit=15,file='matrix-phy-pop.dat',status='unknown')
do ii=1,igt
    write(15,*) x(ispecies(ii,1)),y(ispecies(ii,1)),ispidx(ispecies(ii,1))
end do
CLOSE(15)

902 FORMAT(1000(1x,i7))
903 FORMAT(200000i1)
END SUBROUTINE WRITEPHYLOGENY

SUBROUTINE WRITEFULLPHYLOGENY
USE GLOBALS
OPEN(unit=15,file='matrix-fullphy.dat',status='unknown')
WRITE(15,*) igtfull,iitime,nb
do ii=1,igtfull
    write(15,902) (t_sp(ii,jj),jj=1,igtfull)
end do
CLOSE(15)

OPEN(unit=15,file='ext-times.dat',status='unknown')
WRITE(15,*) igtfull,iitime,nb
do ii=1,igtfull
    write(15,902) ext_time(ii)
end do
CLOSE(15)

OPEN(unit=15,file='ext-sizes.dat',status='unknown')
WRITE(15,*) igtfull,iitime,nb
do ii=1,igtfull
    write(15,902) size_ext(ii),ext_time(ii),sister_sp(ii)
end do
CLOSE(15)

902 FORMAT(1000(1x,i7))
END SUBROUTINE WRITEFULLPHYLOGENY

SUBROUTINE WRITEFITNESS
USE GLOBALS
OPEN(unit=15,file='fitness.dat',status='unknown')
do ii=1,nct
    write(15,902) fitness(ii),distmitonuc(ii),ispidx(ii)
end do
CLOSE(15)
902 FORMAT(2(2x,f7.5),2x,i5)
END SUBROUTINE WRITEFITNESS

END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module random seed
MODULE randomseed
CONTAINS

SUBROUTINE initialize
USE globals, ONLY: aux,iseed,input_seed
iseed = 0
iseed(12) = input_seed
CALL RANDOM_SEED(put=iseed)
CALL RANDOM_NUMBER(aux)
END SUBROUTINE initialize
    
END MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! begin main program
SUBROUTINE submitonuclear(vec_int,vec_real)
USE globals
USE readwrite
USE randomseed

INTEGER, INTENT(IN), DIMENSION(12) :: vec_int
REAL, INTENT(IN), DIMENSION(6) :: vec_real

INTEGER(1) :: one,gl(5000,15000),gml(5000,5000),sl(5000)
INTEGER :: indxorg(5000),xl(5000),yl(5000)
INTEGER xj(32),yj(32),ichoose(1)
INTEGER iix,iiy,idensity
REAL :: neighfit(5000)


! input data
nf1 = vec_int(1)   ! Note: nf1 >= nf2
nf2 = vec_int(2)
nco = vec_int(3)
nrmax = vec_int(4)
mnpm = vec_int(5)
njump = vec_int(6)
nb = vec_int(7)
nbm = vec_int(8)
rg = vec_int(9)
ntime = vec_int(10)
iread = vec_int(11)
input_seed = vec_int(12)

mut = vec_real(1)
mutmit = vec_real(2)
diff = vec_real(3)
qmat = vec_real(4)
radius = vec_real(5)
width = vec_real(6)

nrmax = nrmax + 1
one = 1


! diffusion table
CALL JUMPTABLE(xj,yj)

! initialize random number generator
CALL initialize

! open species plot abundance output files
OPEN(UNIT=8,FILE='speciesplot.dat',STATUS='unknown')
OPEN(UNIT=15,FILE='number0.dat',STATUS='unknown')
OPEN(UNIT=19,FILE='hist0.dat',STATUS='unknown')


! model parameters and initializations
idensity = 15   ! maximum number of individuals per site before  mating is idensity/3
rho0 = float(nco)/float(nf1*nf2)
max_species = 500

g = 0  ! identical genomes
gl = 0
gm = 0 ! identical mitochondria
gml = 0
x = 0
y = 0
xl = 0
yl = 0
t = 0  ! initial times to common ancestor
tp = 0 ! aux variables
p1 = 0 ! first parent
p2 = 0 ! second parent
ispidx = 1
igt_old = 1
size_ext = 0  ! population sizes at time of extinction
sister_sp = 0 ! sister species for extinct species
sister_sp_old = 0

! initialize populations
CALL READPOP

WRITE(6,*) 'total area =',nf1*nf2
write(6,*) 'initial pop =',nct
WRITE(6,*) 'average density =', rho0
WRITE(6,*) 'average number of individuals in S =',3.1416*rho0*radius**2

! Place the organisms by location in the poster space
nworg =0
worg = 0
indxorg = 0
DO i=1,nct     
    nworg(x(i),y(i)) = nworg(x(i),y(i)) + 1
    worg(x(i),y(i),nworg(x(i),y(i))) = i
    indxorg(i) = nworg(x(i),y(i))
END DO

! Set up template for the locations that are part of the neighborhood for a radius
ndim = int(4.0*(radius+nrmax)**2)
idqx = 0
idqy = 0
ndq = 0
DO irad = 1,nrmax 
    idq = 0
    rs = radius + irad - 1
    rs2 = rs**2
    nj = int(rs)+1
    DO itx=-nj,nj
        DO ity=-nj,nj
            IF(itx*itx+ity*ity <= rs2) THEN
                idq = idq + 1
                idqx(irad,idq) = itx
                idqy(irad,idq) = ity
            END IF
        END DO
    END DO
ndq(irad) = idq
END DO

write(6,*)
write(6,*) 'partial time   total time   total pop.  species   ind./species' 
write(6,*)

!
! Time evolution: mating, mutation and diffusion
!
ntest = nco
DO j=1,ntime
	!Mating
    knext = 0   ! count individuals of the next generation
    CALL POPFITNESS
    looppop: DO k=1,nct
        CALL FINDNEIG(k,1)
        nlocal = int(rho0*iibound*0.7)
        ! check if extra offspring is produced
        iextra = 0
        IF(nct < ntest .AND. ineighbort < nlocal) THEN
            nthere = nworg(x(k),y(k))
            IF(ineighborg >= mnpm .AND. nthere < idensity/3) THEN
                iextra = 1
            END IF
        END IF

        IF(iextra == 1) THEN
            knext = knext + 1
            !
            ! choose mate according to fitness
            neighfit = 0.0
            neighfit(1) = fitness(neig(1))
            DO i=2,ineighborg
                neighfit(i) = neighfit(i-1)+fitness(neig(i))
            END DO
            neighfit = neighfit/neighfit(ineighborg)
            CALL RANDOM_NUMBER(aux)
            ichoose = minloc(neighfit-aux,mask=(neighfit - aux > 0))
            kmate = neig(ichoose(1))
            ! end choosing mate
            !
            DO kc=1,nb
                CALL RANDOM_NUMBER(aux)
                IF(aux < 0.5) THEN
                    gl(knext,kc) = g(k,kc)
                ELSE
                    gl(knext,kc) = g(kmate,kc)
                END IF
                CALL RANDOM_NUMBER(aux)
                IF(aux < mut) gl(knext,kc) = one-gl(knext,kc)   ! Mutation
            END DO
            ifemale = k
            imale = kmate
            IF(s(k) == 1) THEN
                ifemale = kmate
                imale = k
            END IF
            gml(knext,1:nbm) = gm(ifemale,1:nbm)
            DO kc=1,nbm
                CALL RANDOM_NUMBER(aux)
                IF(aux < mutmit) gml(knext,kc) = one-gml(knext,kc)   ! Mutation
            END DO
            xl(knext) = x(k)
            yl(knext) = y(k)
            sl(knext) = 1
            CALL RANDOM_NUMBER(aux)
            IF(aux < 0.5) sl(knext) = 0
            ! save parents
            p1(knext) = ifemale
            p2(knext) = imale

            ! regular offspring - same pair k-kmate
            knext = knext + 1
            DO kc=1,nb
                CALL RANDOM_NUMBER(aux)
                IF(aux < 0.5) THEN
                    gl(knext,kc) = g(k,kc)
                ELSE
                    gl(knext,kc) = g(kmate,kc)
                END IF
                CALL RANDOM_NUMBER(aux)
                IF(aux < mut) gl(knext,kc) = one-gl(knext,kc)   ! Mutation
            END DO
            gml(knext,1:nbm) = gm(ifemale,1:nbm)
            DO kc=1,nbm
                CALL RANDOM_NUMBER(aux)
                IF(aux < mutmit) gml(knext,kc) = one-gml(knext,kc)   ! Mutation
            END DO
            xl(knext) = x(k)
            yl(knext) = y(k)
            sl(knext) = 1
            CALL RANDOM_NUMBER(aux)
            IF(aux < 0.5) sl(knext) = 0
            ! save parents
            p1(knext) = ifemale
            p2(knext) = imale
        ELSE
			! only regular offspring
            irad = 1
            DO WHILE(ineighbor < 2)
                irad = irad + 1
                IF (irad == nrmax + 1) THEN
                    CYCLE looppop
                END IF
                CALL FINDNEIG(k,irad)
            END DO
            kmother = k
            irad = 1
            CALL FINDMATE(k,kmother,kmate,irad)
            IF(ineighborg < 2) THEN
                CYCLE looppop
            ELSE
                knext = knext + 1
                DO kc=1,nb
                    CALL RANDOM_NUMBER(aux)
                    IF(aux < 0.5) THEN
                        gl(knext,kc) = g(kmother,kc)
                    ELSE
                        gl(knext,kc) = g(kmate,kc)
                    END IF
                    CALL RANDOM_NUMBER(aux)
                    IF(aux < mut) gl(knext,kc) = one-gl(knext,kc)   ! Mutation
                END DO
                ifemale = kmother
                imale = kmate
                IF(s(kmother) == 1) THEN
                    ifemale = kmate
                    imale = kmother
                END IF

                gml(knext,1:nbm) = gm(ifemale,1:nbm)
                DO kc=1,nbm
                    CALL RANDOM_NUMBER(aux)
                    IF(aux < mutmit) gml(knext,kc) = one-gml(knext,kc)   ! Mutation
                END DO
                xl(knext) = x(k)
                yl(knext) = y(k)
                sl(knext) = 1
                CALL RANDOM_NUMBER(aux)
                IF(aux < 0.5) sl(knext) = 0
                ! save parents
                p1(knext) = ifemale
                p2(knext) = imale
            END IF
        END IF
    END DO looppop

	!update number of individuals, genomes and sex
    nct = knext
    g = gl
    gl = 0
    gm = gml
    gml = 0
    s = sl
    sl = 0

    ! update time matrix
    do ii=1,nct-1
        do jj=ii+1,nct
            it11 = t(p1(ii),p1(jj))
            it12 = t(p1(ii),p2(jj))
            it21 = t(p2(ii),p1(jj))
            it22 = t(p2(ii),p2(jj))
            itlm = min(it11,it12,it21,it22)
            tp(ii,jj) = itlm + 1
            tp(jj,ii) = itlm + 1
        end do
    end do
    t = tp
    tp = 0

	!Diffusion
    IF(diff /= 0.0) THEN
        DO i=1,nct
            CALL RANDOM_NUMBER(aux)
            IF(aux < diff) THEN
                ibound = 0
                DO WHILE(ibound == 0)
                    CALL RANDOM_NUMBER(aux)
                    jjump = INT(njump*aux)+1
                    iix = xl(i)+xj(jjump)
                    iiy = yl(i)+yj(jjump)
                    IF (iix > nf1) CYCLE
                    IF (iix < 1) CYCLE
                    IF (iiy > nf2) CYCLE
                    IF (iiy < 1) CYCLE
                    xl(i) = iix
                    yl(i) = iiy
                    ibound = 1
                END DO
            END IF
        END DO
    END IF

    x = xl
    y = yl


    ! Place the organisms by location in the poster space
    nworg =0
    worg = 0
    indxorg = 0
    DO i=1,nct
        nworg(x(i),y(i)) = nworg(x(i),y(i)) + 1
        worg(x(i),y(i),nworg(x(i),y(i))) = i
        indxorg(i) = nworg(x(i),y(i))
    END DO

! calculate species every generation starting at 2
    IF (j >= 2) THEN
        ispidxold = ispidx
        igtfull_old = igtfull
        ext_time_old = ext_time
        t_sp_old = t_sp
        ispv_old = ispv
        size_ext_old = size_ext
        sister_sp_old = sister_sp
        CALL FINDSPECIES
        write(15,*) j+iitime, igt,  sum(distmitonuc)/float(nct)
        DO k=1,igt
            write(19,*) ispv(k)
        END DO
        write(6,*) j,j+iitime,nct,igt,(ispv(k),k=1,igt)
        CALL MRCAT_SPECIES(j+iitime)
        write(6,*) igt,igtfull
        !write(6,*) (ext_time(l),l=1,igtfull)
        write(6,*)
    END IF

END DO  ! end loop in time
CLOSE(15)
CLOSE(19)

CALL FINDSPECIESMIT
write(6,*)
write(6,*) 'mitochondrial species'
write(6,100) igtm
write(6,100) (ispvm(k),k=1,igtm)

CALL FINDSPECIES
write(6,*)
write(6,*) 'nuclear species'
write(6,100) igt
write(6,100) (ispv(k),k=1,igt)

100 FORMAT(50(i4,2x))


DO i=1,igt
    ij = ispv(i)
    DO jjj=1,ij
        jj = ispecies(i,jjj)
        WRITE(8,*) x(jj),y(jj),i
    END DO
END DO

OPEN(unit=19,file='abund0.dat',status='unknown')
DO k=1,igt
    write(19,*) ispv(k)
END DO
CLOSE(19)

iptime = ntime
iitime = iitime + iptime
CALL WRITEPOP
CALL WRITEPHYLOGENY
CALL WRITEFULLPHYLOGENY
CALL WRITEFITNESS

END SUBROUTINE submitonuclear


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MRCAT_SPECIES(current_time)
USE GLOBALS
INTEGER :: parent_sp(500),sphist(500),sister1(500),sister2(500),sister3(500),true_ext(500)
INTEGER current_time,kmax(1),kmaxs(1),kmaxt(1),kmaxq(1)

parent_sp = 0
ext_time = 0
size_ext = 0
sister1 = 0
sister2 = 0
sister3 = 0
sister_sp = 0

! save parent species of extant species
! parent_sp(i) = parent species of leaving species i
!
DO i=1,igt
    sphist = 0
    DO j=1,ispv(i)
        k = ispecies(i,j)              ! jth individual from species i
        kp = p1(k)                     ! parent
        kps = ispidxold(kp)            ! species of the parent
        sphist(kps) = sphist(kps) + 1  ! increment histogram
    END DO
    kmax = maxloc(sphist)
    parent_sp(i) =  kmax(1)     ! take most frequent ancestor species
    sphist(kmax(1)) = 0
    kmaxs = maxloc(sphist)       ! take the second most frequent ancestor
    if(sphist(kmaxs(1)) /= 0) then
        sister1(i) = kmaxs(1)
        sphist(kmaxs(1)) = 0
        kmaxt = maxloc(sphist)       ! take the third most frequent ancestor
        if(sphist(kmaxt(1)) /= 0) then
            sister2(i) = kmaxt(1)
            sphist(kmaxt(1)) = 0
            kmaxq = maxloc(sphist)       ! take the fourth most frequent ancestor
            if(sphist(kmaxq(1)) /= 0) then
                sister3(i) = kmaxq(1)
            end if
        end if
    end if
END DO

! check for true extinctions
true_ext = 0
do i=1,igt_old
    iext = 0
    do j=1,igt
        if(parent_sp(j) == i) iext = 1
        if(sister1(j) == i) iext = 1
        if(sister2(j) == i) iext = 1
        if(sister3(j) == i) iext = 1
    end do
    if(iext == 0) then
        true_ext(i) = 1
    end if
end do

! check for extinctions / reversals / fusions
!
igtfull = igt
DO i=1,igtfull_old
    DO j=1,igt
        IF(i == parent_sp(j)) THEN
            EXIT
        ELSE IF(j == igt) THEN
            igtfull = igtfull + 1     ! extinct species: copy but set extinction time to /= 0
            parent_sp(igtfull) = i
            ext_time(igtfull) = ext_time_old(i)
            size_ext(igtfull) = size_ext_old(i)
            if(sister_sp_old(i) > 0) then
                sister_sp(igtfull) = igtfull
            else
                sister_sp(igtfull) = sister_sp_old(i)
            end if
            if(ext_time_old(i) == 0) ext_time(igtfull) = iitime + ntime - current_time + 1
            if(size_ext_old(i) == 0) size_ext(igtfull) = ispv_old(i)
            if(i <= igt_old) then
                if(true_ext(i) == 1) then
                    sister_sp(igtfull) = -1
                else
                    sister_sp(igtfull) = igtfull
                end if
            end if
            EXIT
        END IF
    END DO
END DO

igt_old = igt

! new list of species is
! 1,2,...,igt,iext1,iext2,..., = extant + extinct species

! ext_times are zero for 1,2,...,igt and ext_time_old for extinct


! craete new time matrix
!
t_sp = 0
DO i1=1,igtfull
    k1 = parent_sp(i1)
    DO i2=i1+1,igtfull
        k2 = parent_sp(i2)
        t_sp(i1,i2) = t_sp_old(k1,k2) + 1
        t_sp(i2,i1) = t_sp(i1,i2)
    END DO
END DO

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE JUMPTABLE(xj,yj)
integer xj(32),yj(32)
xj(1)=0;    yj(1)=1;
xj(2)=1;    yj(2)=0;
xj(3)=0;    yj(3)=-1;
xj(4)=-1;   yj(4)=0;  !first 4 neighbors, dist-max = 1

xj(5)=1;    yj(5)=1;
xj(6)=1;    yj(6)=-1;
xj(7)=-1;   yj(7)=-1;
xj(8)=-1;   yj(8)=1;
xj(9)=0;    yj(9)=2;
xj(10)=2;   yj(10)=0;
xj(11)=0;   yj(11)=-2;
xj(12)=-2;  yj(12)=0;
xj(13)=1;   yj(13)=2;
xj(14)=2;   yj(14)=1;
xj(15)=2;   yj(15)=-1;
xj(16)=1;   yj(16)=-2;
xj(17)=-1;  yj(17)=-2;
xj(18)=-2;  yj(18)=-1;
xj(19)=-2;  yj(19)=1;
xj(20)=-1;  yj(20)=2; !first 20 neighbors, dist-max = 2

xj(21)=0;   yj(21)=3;
xj(22)=3;   yj(22)=0;
xj(23)=0;   yj(23)=-3;
xj(24)=-3;  yj(24)=0;
xj(25)=1;   yj(25)=3;
xj(26)=3;   yj(26)=1;
xj(27)=3;   yj(27)=-1;
xj(28)=1;   yj(28)=-3;
xj(29)=-1;  yj(29)=-3;
xj(30)=-3;  yj(30)=-1;
xj(31)=-3;  yj(31)=1;
xj(32)=-1;  yj(32)=3; !first 32 neighbors, dist-max = 3

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For each individual k, find a kmother and its mate      !
! Mother is usually k, but not always. 
! At the beginning NEIG was called with irad = 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDMATE(k,kmother,kmate,irad)
USE globals
INTEGER, INTENT(IN) :: k
INTEGER, INTENT(OUT) :: kmate
INTEGER, INTENT(INOUT) :: kmother
INTEGER ichoose(1)
REAL :: neighfit(5000)

neighfit = 0.0

! choose a new mother if aux < fitness or if potential mates < P
! the choice may involve increasing irad until new mother has at least 2 potential mates
!
CALL RANDOM_NUMBER(aux)
qnorm = 1.0
IF(fdelta /= 0.0) qnorm = 2.0*qmat*(fmax-fitness(k))/fdelta
IF(aux < qnorm .OR. ineighborg < mnpm) THEN
    ! CHOOSE MOTHER ACCORDING TO FITNESS
    ! compute cumulative normalized fitness vector
    neighfit(1) = fitness(neigspt(1))
    DO i=2,ineighbort
        neighfit(i) = neighfit(i-1)+fitness(neigspt(i))
    END DO
    neighfit = neighfit/neighfit(ineighbort)
    CALL RANDOM_NUMBER(aux)
    ichoose = minloc(neighfit-aux,mask=(neighfit - aux > 0))
    kmother = neigspt(ichoose(1))  ! choose a new mother
    CALL FINDNEIG(kmother,irad) ! find mother's neighbors
    icc = 0
    DO WHILE (ineighborg < 2 .AND. irad < nrmax)  ! choose a neighbor mother with ineighborg > 1
        IF(icc > 1) THEN
            irad = irad + 1  ! increase search radius
            icc = 0
            CALL FINDNEIG(k,irad)
            ! compute cumulative normalized fitness vector
            IF(ineighbort /= 0) THEN
                neighfit(1) = fitness(neigspt(1))
                DO i=2,ineighbort
                    neighfit(i) = neighfit(i-1)+fitness(neigspt(i))
                END DO
                neighfit = neighfit/neighfit(ineighbort)
            END IF
        END IF
        CALL RANDOM_NUMBER(aux)  !try a different neighbor within the same irad if icc=0
        ichoose = minloc(neighfit-aux,mask=(neighfit - aux > 0))
        kmother = neigspt(ichoose(1))
        ! if there are no spatial neighbors kmother=0 and need to increase neighborhood size
        IF(kmother /= 0) THEN
            CALL FINDNEIG(kmother,irad) ! find mother's neighbors
        END IF
        icc = icc + 1
    END DO
    IF(irad > nrmax) WRITE(6,*) 'picking last individual in search'
END IF


! choose a mate for the mother
!
IF (ineighborg >= 2) THEN
    ! choose mate according to fitness
    neighfit = 0.0
    neighfit(1) = fitness(neig(1))
    DO i=2,ineighborg
        neighfit(i) = neighfit(i-1)+fitness(neig(i))
    END DO
    neighfit = neighfit/neighfit(ineighborg)
    CALL RANDOM_NUMBER(aux)
    ichoose = minloc(neighfit-aux,mask=(neighfit - aux > 0))
    kmate = neig(ichoose(1))
    ! end choosing mate
END IF

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Find spatial and genetic neighbors of kmother !
!   only neighbors of oposite sex are selected   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDNEIG(kmother,irad)
USE globals
INTEGER, INTENT(IN) :: kmother

ineighbor = 0
ineighborg = 0
ineighbort = 0
iibound = 0
neig = 0    ! genetic neighbors
neigsp = 0  ! spatial neighbors of oposite sex
neigspt = 0 ! total spatial neighbors
! the focal individual kmother is not included as a neighbor of itself
ix = x(kmother)
iy = y(kmother)
loop1: DO isite = 1,ndq(irad) 
    ix1 = ix+idqx(irad,isite)
    iy1 = iy+idqy(irad,isite)
    IF (ix1 > nf1) CYCLE
    IF (ix1 < 1) CYCLE
    IF (iy1 > nf2) CYCLE
    IF (iy1 < 1) CYCLE
    iibound = iibound + 1 ! number of sites in the neighborhood
    loop2: DO iworg = 1,nworg(ix1,iy1)
        korg = worg(ix1,iy1,iworg)
        if(korg /= kmother) then
            ineighbort = ineighbort + 1
            neigspt(ineighbort) = korg
            if(s(korg) /= s(kmother)) then
                ineighbor = ineighbor + 1
                neigsp(ineighbor) = korg
                dista = 0
                DO l=1,nb
                    dista = dista + ABS(g(kmother,l)-g(neigsp(ineighbor),l))
                    IF(dista > rg) CYCLE loop2
                END DO
                ineighborg = ineighborg + 1
                neig(ineighborg) = neigsp(ineighbor)
            end if
        end if
    END DO loop2
END DO loop1

END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find species using nuclear chromosome         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDSPECIES
USE globals, ONLY: nct,nb,rg,g,igt,ispecies,ispv,ispidx
INTEGER :: species(5000),auxy1(5000),auxy2(5000)

itot = 0  ! count total population in groups
igt = 0   ! count number of groups
i2 = nct
DO i=1,i2  ! initialize aux2 -- contains all individuals not yet classified
    auxy2(i) = i
END DO 
	
DO WHILE (itot < nct)
    icr = auxy2(1)   !take first individual and find its species
    isp = 0
    ispold = 1
    i1 = 0
    auxy1 = 0
    loop1: DO i=1,i2
        ii = auxy2(i)
        dista = 0
        DO l=1,nb
            IF(g(icr,l) /= g(ii,l)) dista = dista + 1
            IF(dista > rg) THEN
                i1 = i1 + 1
                auxy1(i1) = ii      !put creatures with dist > rg into aux1
                CYCLE loop1
            END IF
        END DO
        isp = isp + 1
        species(isp) = ii   !collect individuals with dist <= rg from icr
    END DO loop1

    !check if individuals in aux1 have to be included; put the rest in aux2
    itest = 1
    DO WHILE(itest /= 0)
        i2 = 0
        auxy2 = 0
        itest = 0
        isp0 = isp
        IF(i1 /= 0) THEN
            loop2: DO i=1,i1
                DO ji=ispold+1,isp0
                    dista = 0
                    DO l=1,nb
                        IF(g(auxy1(i),l) /= g(species(ji),l)) dista = dista + 1  !HERE
                        IF(dista > rg) EXIT
                    END DO
                    IF(dista <= rg) THEN
                        isp = isp + 1
                        species(isp) = auxy1(i)   ! colect the aux1 individual
                        itest = 1                 ! indicates that the process has to be repeated
                        CYCLE loop2
                    END IF
                END DO
                i2 = i2 + 1
                auxy2(i2) = auxy1(i)  ! put individual in aux2
            END DO loop2
        END IF
        auxy1 = auxy2   ! aux1 contains the creatures not in the species
        i1 = i2
        ispold = isp0
    END DO

    itot = itot + isp    !total number of individuals classified into species
    igt = igt + 1        !number of species

	! save species info
    DO i=1,isp
        ispecies(igt,i) = species(i)
        ispidx(species(i)) = igt
    END DO
    ispv(igt) = isp          ! number of individuals in species

END DO

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find species using mitochondrial chromosome         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDSPECIESMIT
USE globals, ONLY: nct,nbm,rg,gm,igtm,ispeciesm,ispvm,ispidxm
INTEGER :: species(5000),auxy1(5000),auxy2(5000)

itot = 0  ! count total population in groups
igtm = 0   ! count number of groups
i2 = nct
DO i=1,i2  ! initialize aux2 -- contains all individuals not yet classified
    auxy2(i) = i
END DO

DO WHILE (itot < nct)

    icr = auxy2(1)   !take first individual and find its species
    isp = 0
    ispold = 1
    i1 = 0
    auxy1 = 0
    loop1: DO i=1,i2
        ii = auxy2(i)
        dista = 0
        DO l=1,nbm
            IF(gm(icr,l) /= gm(ii,l)) dista = dista + 1
            IF(dista > rg) THEN
                i1 = i1 + 1
                auxy1(i1) = ii      !put creatures with dist > rg into aux1
                CYCLE loop1
            END IF
        END DO
        isp = isp + 1
        species(isp) = ii   !collect individuals with dist <= rg from icr
    END DO loop1

    !check if individuals in aux1 have to be included; put the rest in aux2
    itest = 1
    DO WHILE(itest /= 0)
        i2 = 0
        auxy2 = 0
        itest = 0
        isp0 = isp
        IF(i1 /= 0) THEN
            loop2: DO i=1,i1
                DO ji=ispold+1,isp0
                    dista = 0
                    DO l=1,nbm
                        IF(gm(auxy1(i),l) /= gm(species(ji),l)) dista = dista + 1  !HERE
                        IF(dista > rg) EXIT
                    END DO
                    IF(dista <= rg) THEN
                        isp = isp + 1
                        species(isp) = auxy1(i)   ! colect the aux1 individual
                        itest = 1                 ! indicates that the process has to be repeated
                        CYCLE loop2
                    END IF
                END DO
                i2 = i2 + 1
                auxy2(i2) = auxy1(i)  ! put individual in aux2
            END DO loop2
        END IF
        auxy1 = auxy2   ! aux1 contains the creatures not in the species
        i1 = i2
        ispold = isp0
    END DO

    itot = itot + isp    !total number of individuals classified into species
    igtm = igtm + 1        !number of species

	! save species info
    DO i=1,isp
        ispeciesm(igtm,i) = species(i)
        ispidxm(species(i)) = igtm
    END DO
    ispvm(igtm) = isp          ! number of individuals in species

END DO

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the fitness of kmother's neighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE POPFITNESS
USE globals

fitness = 0.0
anbm = float(nbm)

w2 = width**2
DO k=1,nct
    dista = 0
    DO l=1,nbm
        dista = dista + ABS(g(k,l)-gm(k,l))
    END DO
    dista = dista/anbm
    distmitonuc(k) = dista
    fitness(k) = exp(-0.5*dista**2/w2)
END DO
fmax = maxval(fitness)
fmin = minval(fitness)
fdelta = fmax-fmin

END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE init_random_seed()
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)
END SUBROUTINE
