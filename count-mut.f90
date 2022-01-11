! program count-mut.f90

! counts the number of mutations in the nuclear and mitochondrial dna's

! module defining global variables
MODULE globals
INTEGER(1), ALLOCATABLE, SAVE :: g(:,:),gm(:,:)
INTEGER, SAVE :: nb,nbm,nct
CHARACTER*100 :: stname,folder,path
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module that reads and writes the populations
MODULE readwrite
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
        ALLOCATE (g(nct,nb),gm(nct,nbm))
		DO i=1,nct
			READ (10,901) ix,iy,(g(i,j),j=1,nb)
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
PROGRAM countmut
USE globals
USE readwrite
INTEGER, ALLOCATABLE :: mutmit(:),mutnuc(:),mutnucneutral(:),mutnuctot(:)


read(*,*) i
write(folder, '( "run_",I2.2)') i
path = trim('./')//trim(folder)//trim('/')

CALL READPOP

ALLOCATE (mutmit(nct),mutnuc(nct),mutnucneutral(nct),mutnuctot(nct))

mutnuc = 0
mutmit = 0
mutnuctot = 0
mutnucneutral = 0
DO j=1,nct
    DO k=1,nbm
        mutnuc(j) = mutnuc(j) + g(j,k)
    END DO
    DO k=nbm+1,nb
        mutnucneutral(j) = mutnucneutral(j) + g(j,k)
    END DO
    DO k=1,nbm
        mutmit(j) = mutmit(j) + gm(j,k)
    END DO
    mutnuctot(j) = mutnuc(j) + mutnucneutral(j)
!    write(6,*) j,mutmit(j),mutnuc(j),mutnucneutral(j),mutnuctot(j)
END DO


OPEN(unit=15,file=trim(path)//'mit-mutations.dat')
OPEN(unit=16,file=trim(path)//'nuc-mutations.dat')
OPEN(unit=17,file=trim(path)//'nuc-neutral-mutations.dat')
OPEN(unit=18,file=trim(path)//'nuc-tot-mutations.dat')
DO i=1,nct
    write(15,*) mutmit(i)
    write(16,*) mutnuc(i)
    write(17,*) mutnucneutral(i)
    write(18,*) mutnuctot(i)
END DO
CLOSE(15)
CLOSE(16)
CLOSE(17)
CLOSE(18)

t = 2000.0
anorm = t*nbm*nct
anorm1 = t*nct*(nb-nbm)
anorm2 = t*nb*nct
write(6,*)
write(6,*) 'average effective mitochondrial mutation rate =',sum(mutmit)/anorm
write(6,*)
write(6,*) 'average effective nuclear coupled mutation rate =',sum(mutnuc)/anorm
write(6,*)
write(6,*) 'average effective nuclear uncoupled mutation rate =',sum(mutnucneutral)/anorm1
write(6,*)
write(6,*) 'average effective total nuclear mutation rate =',sum(mutnuctot)/anorm2

END PROGRAM countmut


