! tree-ngas.f90

! f2py -c --fcompiler=gnu95 -m FortranTree FortranTree.f90

! generates genealogical tree from the time matrix
! uses average distances to groups
! computes gamma, alpha and sackin indexes

! input files:
! matrix.dat - genetic distance matrix
! matrix-color.dat - species of each individual, to be used as colors

! output files:
! tree.dat - genealogic tree
! tree-species.dat - colors of individuals at the base of tree according to their species
! ngas.dat - n, gamma, alpha, sackin

! includes genetic distances and colors

SUBROUTINE subtree
REAL :: t(100,100),tori(100,100)      ! common ancestor matrix
REAL:: gtime(0:100)  ! groups variables
INTEGER :: gsize(0:100),g(100,100),idxg(100)  ! groups variables
INTEGER :: k1(100),k2(100),color(100)  ! branch type variables
INTEGER :: nc,itmin(2),p
REAL :: x1(100),x2(100),xav(100)
REAL :: times(100)

OPEN(unit=12,file='matrix-phy.dat',status='old')
color = 0
read(12,*) nc,nt,nb
do i=1,nc
    read(12,*) (t(i,j),j=1,nc)
    color(i) = i
end do
close(12)
tori = t

g = 0
gsize = 0
gtime = 0.0
gsize(0) = 1

! Construct groups and find branch type of each new
! individual added to a group
! group type is
! (i)   if k1 = k2 = 0
! (iia) if k1 = 0, k2 /= 0
! (iib) if k1 /= 0, k2 = 0
! (iii) if k1 /= 0, k2 /= 0

! construct first group
! branch type is always (i)
p = 1
itmin = minloc(t,mask = t > 0)
ind1 = itmin(1)
ind2 = itmin(2)
g(p,1) = ind1
g(p,2) = ind2
gsize(p) = 2
gtime(p) = t(ind1,ind2)
k1(p) = 0
k2(p) = 0

! set distances to ind2 to average values between ind1 and ind2
! for tori also set distances to ind1
DO k=1,nc
    IF(t(k,ind2) /= 0) THEN
        t(k,ind2) = 0.5*(t(k,ind2)+t(k,ind1))
        t(ind2,k) = t(k,ind2)
        tori(k,ind2) = t(k,ind2)
        tori(ind2,k) = t(k,ind2)
        tori(k,ind1) = t(k,ind2)
        tori(ind1,k) = t(k,ind2)
    END IF
END DO
tori(ind1,ind2) = gtime(p)
tori(ind2,ind1) = gtime(p)
tori(ind1,ind1) = 0.0
tori(ind2,ind2) = 0.0

! set line and column ind1 to zero
t(1:nc,ind1) = 0
t(ind1,1:nc) = 0


DO WHILE (gsize(p) < nc)
    p = p + 1
    itmin = minloc(t,mask = t > 0)
    ind1 = itmin(1)
    ind2 = itmin(2)
    gtime(p) = t(ind1,ind2)
    ki = 0
    kj = 0

    ! find branch type

    !check if ind1 is already in a group
    loop1: do k=p-1,1,-1
        do l=1,gsize(k)
            if(ind1 == g(k,l)) then
                ki = k
                exit loop1
            end if
        end do
    end do loop1

    !check if ind2 is already in a group
    loop2: do k=p-1,1,-1
        do l=1,gsize(k)
            if(ind2 == g(k,l)) then
                kj = k
                exit loop2
            end if
        end do
    end do loop2

    k1(p) = ki
    k2(p) = kj

    if( ki /= 0 .and. kj /= 0) then
        ! both ind1 and ind2 belong to groups: merge lists
        do k=1,gsize(ki)
            g(p,k) = g(ki,k)
        end do
        do k=1,gsize(kj)
            g(p,gsize(ki)+k) = g(kj,k)
        end do
        gsize(p) = gsize(ki) + gsize(kj)
    else if(ki /= 0) then
        ! only ind1 belongs to a group: new group = old group + ind2
        do k=1,gsize(ki)
            g(p,k) = g(ki,k)
        end do
        g(p,gsize(ki)+1) = ind2
        gsize(p) = gsize(ki) + 1
    else if(kj /= 0) then
        ! only ind2 belongs to a group: new group = old group + ind1
        do k=1,gsize(kj)
            g(p,k) = g(kj,k)
        end do
        g(p,gsize(kj)+1) = ind1
        gsize(p) = gsize(kj) + 1
    else
        ! create new group (cherry)
        g(p,1) = ind1
        g(p,2) = ind2
        gsize(p) = 2
    end if
    ! set distances to ind2 to average values between ind1 and ind2
    DO k=1,nc
        IF(t(k,ind2) /= 0 .and. t(k,ind1) /= 0) THEN
            gtot = float(gsize(ki)+gsize(kj))
            wi = float(gsize(ki))/gtot
            wj = float(gsize(kj))/gtot
            t(k,ind2) = wj*t(k,ind2)+ wi*t(k,ind1)
            t(ind2,k) = t(k,ind2)
            tori(k,ind2) = t(k,ind2)
            tori(ind2,k) = t(k,ind2)
            tori(k,ind1) = t(k,ind2)
            tori(ind1,k) = t(k,ind2)
        END IF
    END DO
    tori(ind1,ind2) = gtime(p)
    tori(ind2,ind1) = gtime(p)
    tori(ind1,ind1) = 0.0
    tori(ind2,ind2) = 0.0

    ! set line and column to zero
    t(1:nc,ind1) = 0
    t(ind1,1:nc) = 0
END DO

write(6,*)
write(6,*) 'order of species at base of the tree'
write(6,100) (g(p,i),i=1,gsize(p))
write(6,*)

! re-order locations on x-axis according to final order
do i=1,gsize(p)
    idxg(g(p,i)) = i  
end do

! position of groups as averages over individuals
xav = 0.0
do i=1,p
    do j=1,gsize(i)
        xav(i) = xav(i) + idxg(g(i,j))
    end do
    xav(i) = xav(i)/float(gsize(i))
end do

OPEN(unit=15,file='tree.dat',status='unknown')

!write tree
!write(6,*) 'group        x1        x2         t01       t02         time '
!write(6,*)
do i=1,p
    if( k1(i) /= 0 .and. k2(i) /= 0) then
        x1(i) = xav(k1(i))
        x2(i) = xav(k2(i))
    else if(k1(i) /= 0) then
        x1(i) = xav(k1(i))
        x2(i) = idxg(g(i,gsize(i)))
    else if(k2(i) /= 0) then
        x1(i) = idxg(g(i,gsize(i)))
        x2(i) = xav(k2(i))
    else
        x1(i) = idxg(g(i,1))
        x2(i) = idxg(g(i,2))
    end if
!    write(6,102) i,x1(i),x2(i),gtime(k1(i)),gtime(k2(i)),gtime(i)
    write(15,103) x1(i),gtime(k1(i))
    write(15,103) x1(i),gtime(i)
    write(15,103) x2(i),gtime(i)
    write(15,103) x2(i),gtime(k2(i))
    write(15,100)  -1,-1
end do
    write(15,103) xav(p),gtime(p)
    write(15,103) xav(p),gtime(p) + 1
    write(15,100) -1,-1
close(15)

OPEN(unit=16,file='order.dat',status='unknown')
do i=1,nc
    write(16,*) idxg(g(p,i)),color(g(p,i))
end do
close(16)

DO i=1,p
    times(i) = tori(g(p,i),g(p,i+1))
END DO

! get matrix in canonical form
t = 0
do i=1,nc-1
    t(i,i+1) = times(i)
end do
do j=2,nc-1
    do k=1,nc-j
        t(k,k+j) = max(t(k,k+j-1),t(k+1,k+j))
    end do
end do
t = t + transpose(t)
tori = t

CALL gamma(t,nc,gm)

CALL sackin(tori,nc,sa,sn,sy)

CALL alpha_calc(nc,gm,alpha)

OPEN(unit=15,file='ntgasns-upgma.dat')
    write(15,107) nt,nc,gm,alpha,sn,sa
CLOSE(15)

100 format(100i5)
103 format(1x,f10.3,1x,f10.3)
107 format(i5,2x,i4,3(2x,f7.3),2x,f6.0)

END SUBROUTINE subtree



!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gamma(t,nc,gm)
    REAL, SAVE :: xx1, xxd, TT, SE, SI, DD
    INTEGER :: ind1, ind2
    REAL :: gtime(1:100), gt(1:100), itmin(2), t(100,100)

    gtime = 0
    do i=1,nc-1
        itmin = minloc(t,mask = t > 0)
        ind1 = itmin(1)
        ind2 = itmin(2)
        gtime(i) = t(ind1,ind2)

        ! set line and column to zero
        t(1:nc,ind1) = 0
        t(ind1,1:nc) = 0
    end do

    xxd=gtime(1)
    gt(1)=0
    gt(nc)=xxd

    do i=2,nc-1
        xx1=gtime(i)
        gt(nc-i+1)=xx1-xxd
        xxd=xx1
    end do

    TT=0
    do i=2,nc
        TT=TT+i*gt(i)
    end do

    SE=0
    do i=2,nc-1
        SI=0
        do j=2,i
            SI=SI+j*gt(j)
        end do
        SE=SE+SI
    end do

    DD=TT*sqrt(1.0/(12*(nc-2)))
    gm=((1.0/(nc-2))*SE - TT/2.0)/DD

RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sackin(t,n,sa,sn,sy)
REAL :: row(100)
REAL t(100,100)

index = 0
do i=1,n
    row = 0
    irow = 0
    do j=1,n
        if(j == 1) then
            irow = 1
            row(1) = t(i,j)
        else
            ieq = 0
            do k=1,irow
                if(t(i,j) == row(k)) ieq = 1
            end do
            if(ieq == 0) then
                irow = irow + 1
                row(irow) = t(i,j)
            end if
        end if
    end do
    irow = irow - 1
    index = index + irow
end do
write(6,*)

ys = 0.0
do i=2,n
    ys = ys + 1.0/float(i)
end do
sa = index
sy = 2.0*float(n)*ys

! Compute mean square deviation for Sackin

h1 = 0.0
h2 = 0.0
do j=1,n
    h1 = h1 + 1.0/float(j)
    h2 = h2 + 1.0/float(j*j)
end do
an = float(n)
an2 = an**2
sigma = sqrt((7.0 - 4.0*h2)*an2 -(2.0*h1 + 1.0)*an)

sn = (index-sy)/sigma

RETURN
END

SUBROUTINE alpha_calc(nc,gm,alpha)
error = 0.001
alpha0 = 0.0
deltaalpha = 1.0

DO WHILE(abs(deltaalpha) > error)
    CALL gammafunc(nc,alpha0,gamma0,gammap0)
    deltaalpha = (gm - gamma0)/gammap0
    alpha0 = alpha0 + deltaalpha
END DO
alpha = alpha0

RETURN
END

SUBROUTINE gammafunc(nc,alpha0,gamma0,gammap0)
REAL :: aux,aux1,aux2,ai,deno
REAL :: tk(2:nc), tkp(2:nc)

    tk = 0.0
    tkp = 0.0
    tk(2) = 2.0**(1.0-alpha0)
    tkp(2) = -log(2.0)*tk(2)
    do i=3,nc
        ai = float(i)
        aux = ai**(1.0-alpha0)
        tk(i) = tk(i-1) + aux
        tkp(i) = tkp(i-1) - log(ai)*aux
    end do

    deno = tk(nc)/sqrt(float(nc-2)*12.0)

    aux1 = 0.0
    aux2 = 0.0
    do i=2,nc-1
        aux1 = aux1 + tk(i)
        aux2 = aux2 + tkp(i)
    end do

    gamma0 = (aux1/float(nc-2) - 0.5*tk(nc))/deno
    gammap0 = (aux2 - aux1*tkp(nc)/tk(nc))/(deno*float(nc-2))

RETURN
END
