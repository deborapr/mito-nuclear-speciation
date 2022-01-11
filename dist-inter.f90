! program dist-inter.f90

! calculates nuclear and mitochondrial genetic distances between pair of individuals of different species

! module defining global variables
MODULE globals
INTEGER(1), ALLOCATABLE, SAVE :: gm(:,:), g(:,:), ix(:), iy(:)
INTEGER, SAVE :: nb,nbm,nct,igt
INTEGER, ALLOCATABLE, SAVE :: ispecies(:,:),ispv(:),ispidx(:)
REAL, SAVE :: rg
CHARACTER*100 :: stname,folder,path
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
    name1 = trim(path)//'pop-new.dat'
    OPEN(UNIT=10,FILE=name1,STATUS='old')
		READ(10,*) iitime,nct,nf1,nf2
		READ(10,*) amutd,amutmitd,diffd !dumb only
		READ(10,*) radiusd,rg,nb,nbm,mnpmd !dumb only
        ALLOCATE (g(nct,nb),gm(nct,nbm),ix(nct), iy(nct))
		DO i=1,nct
			READ (10,901) ix(i),iy(i),(g(i,j),j=1,nb)
		END DO
		DO i=1,nct
			READ (10,902) (gm(i,j),j=1,nbm)
		END DO
    CLOSE(10)
901 FORMAT(i4,1x,i4,1x,200000i1)
902 FORMAT(200000i1)
	END SUBROUTINE READPOP

END MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! begin main program
PROGRAM simmatrix
USE globals
USE readwrite

read(*,*) i
write(folder, '( "run_",I2.2)') i
path = trim('./')//trim(folder)//trim('/')

CALL READPOP
ALLOCATE (ispidx(nct),ispv(nct),ispecies(nct,nct))
CALL FINDSPECIES
write(6,*) 'species sizes'
DO i=1,igt
    write(6,*) i,ispv(i)
END DO


stname = trim(path)//'nuclear_dist_sp.dat'
OPEN(unit=15,file=stname)

!dist = 0
ii = 0
DO i1=1,igt
    DO j1=i1+1,igt
        DO i2=1,ispv(i1)
            ind1 = ispecies(i1,i2)
            DO j2=1,ispv(j1)
                ind2 = ispecies(j1,j2)
                idist = 0
                DO l=1,nb
                    IF(g(ind1,l) /= g(ind2,l)) idist = idist + 1
                END DO
                write(15,*) i1, j1, idist
            END DO
        END DO
    END DO
END DO
CLOSE(15)

stname = trim(path)//'mito_dist_sp.dat'
OPEN(unit=15,file=stname)

!dist = 0
ii = 0
DO i1=1,igt
    DO j1=i1+1,igt
        DO i2=1,ispv(i1)
            ind1 = ispecies(i1,i2)
            DO j2=1,ispv(j1)
                ind2 = ispecies(j1,j2)
                idist = 0
                DO l=1,nbm
                    IF(gm(ind1,l) /= gm(ind2,l)) idist = idist + 1
                END DO
                write(15,*) i1, j1, idist
            END DO
        END DO
    END DO
END DO
CLOSE(15)

END PROGRAM simmatrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find species using nuclear chromosome         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDSPECIES
USE globals
INTEGER, ALLOCATABLE :: species(:),auxy1(:),auxy2(:)
INTEGER iix,iiy

ALLOCATE (species(nct),auxy1(nct),auxy2(nct))

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
			loop2:	DO i=1,i1
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

