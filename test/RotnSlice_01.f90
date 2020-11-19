!     Last change:  AL   17 Jul 2008    1:50 pm
module global_var
   REAL(KIND=8) :: v_luce=299792458d0,ev_joule=1.60217733d-19,m0=510998.9d0
   INTEGER :: O_stat,R_stat,Npartic
   character(len=19) :: Astra_fmt='(8(x,es11.4e2),2i4)'
   character(len=20) :: Astra_fmt_d='(8(x,es19.12e2),2i4)'
   type part
           real(kind=8) :: x,y,z,px,py,pz,t,q
           integer :: f1,f2
   end type
   type (part),dimension(:), allocatable :: S_bunch !They are the Full bunch and the Sliced bunch
end module

module math_tools
        contains
        function var (array)
           implicit none
           real(kind=8) :: var,x
           real(kind=8),intent(in),dimension(:) :: array
           x=sum(array)/size(array)
           var=sum((array-x)**2)/size(array)
           end function var
        function mean (array)
           implicit none
           real(kind=8) :: mean,x
           real(kind=8),intent(in),dimension(:) :: array
           mean=sum(array)/size(array)
           end function mean
end module

module subr
  USE global_var
  contains
  subroutine bunch_slicer(zmin_slice,zmax_slice,F_bunch,S_bunch,last)
  IMPLICIT NONE
  type (part),dimension(:), allocatable,intent(in) :: F_bunch
  type (part),dimension(:), allocatable,intent(out) :: S_bunch
  REAL(KIND=8) :: zmin_slice,zmax_slice
  INTEGER      :: i,np
  logical      :: last
  np=0
  open (2000,status='scratch')
  !do i=1,Npartic !ecco l'errore, scorreva anche valori di bunch che non esistemvano più ma erano ancora valorizzati
  do i=1,size(F_bunch)
      if ((.not.last) .and. ((F_bunch(i)%z >= zmin_slice).AND.(F_bunch(i)%z < zmax_slice))) then
         np=np+1
         write (2000,*) F_bunch(i)

         else if ((last) .and. ((F_bunch(i)%z >= zmin_slice).AND.(F_bunch(i)%z <= zmax_slice))) then
         np=np+1
         write (2000,*) F_bunch(i)
         end if
         end do
  ALLOCATE (S_bunch(np))
  rewind(2000)
  do i=1,np
     read(2000,*) S_bunch(i)
  end do
  close (2000)
  end subroutine

  subroutine  xyz_rot(Ibunch,ang_x,ang_y)
      use global_var
      implicit none
      type (part),dimension(:),allocatable :: Ibunch,Rbunch
      real(kind=8),intent(in) :: ang_x,ang_y
      !real(kind=8),dimension(:),allocatable :: modulo_xyz,modulo_pxpypz
      !print *,'size=',size(Ibunch)
      Rbunch=Ibunch
      !***stesso formalismo delle matrici di rotazioni presenti su wikipedia*** 

      !***rotazione solo in theta_x o attorno all'asse y***
      !Rbunch%x=Ibunch%x*dcos(ang_x)+Ibunch%z*dsin(ang_x)
      !Rbunch%y=Ibunch%y
      !Rbunch%z=-Ibunch%x*dsin(ang_x)+Ibunch%z*dcos(ang_x)
      !!queste sono le matirici corrette
      !Rbunch%x= dcos(ang_x)  0    dsin(ang_x)
      !Rbunch%y=      0       1       0
      !Rbunch%z=-dsin(ang_x)  0    dcos(ang_x)
      
      !***rotazione solo in theta_y o attorno all'asse x***
      !********************************************************************************************
      !***se ruoto intorno all'asse x, con la regola della mano destra, un vettore              ***
      !***che ruota partedo sovrapposto all'asse z e che va verso y apre un angolo negativo,    ***
      !***visto che nel mio formalismo questo ancgolo ï¿½ positivo, moltiplico x -1 i seni        ***
      !********************************************************************************************
      !Rbunch%x=Ibunch%x
      !Rbunch%y=Ibunch%y*dcos(ang_y)-Ibunch%z*dsin(ang_y)
      !Rbunch%z=Ibunch%y*dsin(ang_y)+Ibunch%z*dcos(ang_y)
      !Rbunch%x=     1       0       0
      !Rbunch%y=     0   dcos(ang_y)  -dsin(ang_y)
      !Rbunch%z=     0   dsin(ang_y)   dcos(ang_y)

      !prodotto delle due matrici di rotazione, applicato al vettore colonna (x,y,z)
      !in precedenza ruotavo per angoli negativi, applicando una correzione ad un tilt,
      !adesso ruoto per con angoli positivi, quando correggo moltiplico x -1 l'angolo.
      Rbunch%x=Ibunch%x*dcos(ang_x)-Ibunch%y*dsin(ang_x)*dsin(ang_y)+Ibunch%z*dsin(ang_x)*dcos(ang_y)
      Rbunch%y=Ibunch%y*dcos(ang_y)+Ibunch%z*dsin(ang_y)
      Rbunch%z=-1.0d0*Ibunch%x*dsin(ang_x)-Ibunch%y*dcos(ang_x)*dsin(ang_y)+Ibunch%z*dcos(ang_x)*dcos(ang_y)
      !prodotto delle due matrici di rotazione, applicato al vettore colonna (px,py,pz)
      Rbunch%px=Ibunch%px*dcos(ang_x)-Ibunch%py*dsin(ang_x)*dsin(ang_y)+Ibunch%pz*dsin(ang_x)*dcos(ang_y)
      Rbunch%py=Ibunch%py*dcos(ang_y)+Ibunch%pz*dsin(ang_y)
      Rbunch%pz=-1.0d0*Ibunch%px*dsin(ang_x)-Ibunch%py*dcos(ang_x)*dsin(ang_y)+Ibunch%pz*dcos(ang_x)*dcos(ang_y)
      
      

      Ibunch=Rbunch
  end subroutine
end module      

module inizi_and_messages
        USE global_var
        CHARACTER(LEN=50) :: astra_file,astra_out_file,arg1,arg2,arg3
        real :: real_val
        integer :: N_s !numero di slice
        complex :: theta
        character(len=8) :: mode !"one_file" or "listfile"
        contains
               subroutine iniz_and_msg(arg1,arg2,arg3,N_s,mode)
                 implicit none
                 CHARACTER(LEN=50) :: arg1,arg2,arg3
                 integer :: N_s
                 character(len=8) :: mode !"one_file" or "listfile"
                 if (trim(arg1)=="") then
                       print *,"ERROR: <Bunch file> or <file list> is missing"
                       stop
                       end if
                 if (trim(arg1)=="help") then
                         print *
                         print *,'**************************  INSTRUCTIONS  ***************************************'
                         print *,"- Performs bunch slice analysis or to compute Trace3D code inputs"
                         print *,'- Computes eta_x, etaP_x, eta_y, etaP_y and other beam paramiters'
                         print *,'- Straightens tilted bunches (e.g. passed throught dipole fields)'
                         print *
                         print *,'Angles are computed by <Px>/<Pz> and <Py>/<Pz> OR imposed by the User' 
                         print *
                         print *,'This Code useful to compute fitness functions, based on eta, '
                         print *,'in the GIOTTO software (e.g. dispersion closing after dispersive paths)' 
                         print *
                         print *,'Different WORKING MODES:'
                         print *,'** 1 ** Bunch slice analysis:'
                         print *,'>Bunch_xx.exe <astra_file> [N_slice]'
                         print *
                         print *,'   -Generates <an_out_xx.txt> beam parameters file with eta, eta_prime'
                         print *,'    xx is the astra run index'
                         print *,'    -If N_slice = 1 it returns values for the whole bunch and'
                         print *,'     also values in Trace3D  format'
                         print *,'    -If N_slice is not specified it means N_slice=1'
                         print *,'    -If N_slice > 1 the slice analysis is performed'
                         print *
                         print *,"** 2 ** Elaborate a list of Astra's bunch distributions:"
                         print *,'>Bunch_xx.exe <List_of_file.txt>'
                         print *
                         print *,'   -Generates <an_out_xx.txt> for every file into the list'
                         print *
                         print *,"** 3 ** Bunch rotation by an User defined angle:"
                         print *,'>Bunch_xx.exe <astra_file> <output_file> <(theta_x,theta_y)>'
                         print *
                         print *,'   -It generates <output_file>, an Astra distribution that propagate in the'
                         print *,'    direction (theta_x,theta_y) [rad]. It is a fortran complex format '
                         print *,'    theta_x>0 means z_versor rotating versus x-axis'
                         print *,'    theta_y>0 means z_versor rotating versus y-axis'
                         print *
                         print *,'** FURTHER CAPABILITY **:'
                         print *,'Optional input file <analyzer_in.txt> to be added in working folders' 
                         print *,'   -it contains the NAMELIST "&optional_par" to change following variables' 
                         print *
                         print *,'VARIABLE       |type    |default  |meaning'
                         print *,'freq [MHz]      real     2856.0d6  RF for long-emit [rad-keV] (Trace3D format)'
                         print *,'rotate          logical  TRUE      rotates the bunches'
                         print *,'print_rot_files logical  FALSE     generates or not rotated bunches (Astra format)'
                         print *,'del_lost_part   logical  TRUE      if TRUE delete lost particles (id<-10)'
                         print *
                         print *,' A.Bacci: alberto.bacci@mi.infn.it & M.Rossetti Conti: marcello.rossetti@mi.infn.it'
                         print *,'***********************Milano -october- 2020:*************************************'
                         stop
                         end if

                         if (trim(arg3)=="") then ! analysis of Astra files modes
                                open(300,file=trim(arg1),status='old',iostat=O_stat)
                                if (O_stat /=0) then
                                        print *,"ERROR: The file in input does not exist" 
                                        stop

                                        else
                                        read(300,*,iostat=R_stat) real_val
                                        close(300)
                                        end if

                                if (R_stat==0) then !the file contains real number, it is an Astra file
                                       mode="one_file"
                                       astra_file=arg1  
                                       if (trim(arg2)=="") arg2='1'
                                       read(arg2,*) N_s

                                       else          !the file contains text, it is a list of files
                                       mode="listfile"
                                       N_s=1
                                       end if
                                end if
                          if (trim(arg3)/="") then ! bunch rotation mode with external angles
                                  print *,"mode rotation from user"
                                  mode="one_file"
                                  N_s=1
                                  astra_file=arg1
                                  astra_out_file=trim(arg2)
                                  read(arg3,*) theta
                                  end if
                 end subroutine
               end module inizi_and_messages

!---------------------------------------------------
program main
   USE subr
   USE global_var
   USE math_tools
   USE File_analysis
   USE inizi_and_messages
   IMPLICIT NONE
   type (part),dimension(:),allocatable :: Bunch,bunch_buffer
   REAL(KIND=8),DIMENSION(:),ALLOCATABLE ::pr,p,xp,yp
   
   REAL(KIND=8) :: p_med,x_max,y_max,x_med,y_med,p_j,sig_z,xp_med,yp_med
   REAL(KIND=8) :: x_med_rot,y_med_rot,z_med_rot,px_med_rot,py_med_rot,pz_med_rot !valori rotoluti del pacchetto che non ha ruotato
   REAL(KIND=8) :: zmin_slice,zmax_slice,pz_med,g,px_med,py_med,eta_x,eta_y,etaP_x,etaP_y
   REAL(KIND=8) :: E_var
   real(kind=8) :: x_var,y_var,z_var,xp_var,yp_var,px_var,py_var,pz_var,x_cor,y_cor,z_cor,emit_x,emit_y,emit_z   ,s_pos
   REAL(KIND=8) :: curr,beta,dz,z_med
   real(kind=8) :: ti,tf
   REAL(kind=8) :: m2deg,freq=2856.0d6
   real(kind=8) :: theta_x,theta_y,beam_energy
   logical      :: print_rot_files=.false.,rotate=.true.,del_lost_part=.true.
   INTEGER :: i,k,j,np,Q_perc,Nfile,N_file,n
   integer :: myid
   CHARACTER(LEN=8) :: bunch_name='bunchxxx'
   character(len=17) :: out_file='an_out_xx.txt'
   character(len=14) :: fmt_out_file='(txx,tl2,i3.3)'
   logical :: last=.false.
   namelist /optional_par/ freq,rotate,print_rot_files,del_lost_part
   
   ! Parametri per la namelist di comunicazione RotnSlice-GIOTTO, il prefisso "r_" sta per rot, indica il fatto che sono variaili temporanee
   character(len=17) :: nml_file='an_out_xx.nml'
   character(len=14) :: fmt_nml_file='(txx,tl2,i3.3)'
   integer :: r_np
   REAL(kind=8) :: r_Den,r_Xmax,r_Ymax,r_SigX,r_SigY,r_SigZ,r_SigPx,r_SigPy,r_SigPz,r_EmitX,r_EmitY,r_Curr,r_EtaX,&
   r_EtaPx,r_EtaY,r_EtaPy,r_AlphaX_Twiss,r_AlphaY_Twiss,r_BetaX_Twiss,r_BetaY_Twiss,r_GammaX_Twiss,r_GammaY_Twiss
   namelist /OutParam/ r_np,r_Den,r_Xmax,r_Ymax,r_SigX,r_SigY,r_SigZ,r_SigPx,r_SigPy,r_SigPz,r_EmitX,r_EmitY,r_Curr,r_EtaX,&
   r_EtaPx,r_EtaY,r_EtaPy,r_AlphaX_Twiss,r_AlphaY_Twiss,r_BetaX_Twiss,r_BetaY_Twiss,r_GammaX_Twiss,r_GammaY_Twiss

   call cpu_time(ti)

   call get_command_argument(1,arg1)
   call get_command_argument(2,arg2)
   call get_command_argument(3,arg3)

   call iniz_and_msg(arg1,arg2,arg3,N_s,mode)

   Print *
   Print *
   Print *,"************* Bunch Rotation or Bunch Analysis is running ***********"
   !if ((mode.eq.0).and.(arg3/="")) then
   if (mode.eq."one_file") then
        N_file=1
        !l=len(trim(astra_file))
        write(fmt_out_file(3:4),'(i2)') len(trim(astra_file))
        write(fmt_nml_file(3:4),'(i2)') len(trim(astra_file))
        read (astra_file,fmt=fmt_out_file,iostat=R_stat) myid
        if (R_stat/=0) myid=99 !non ha trovato nessun suffisso numerico nel nome originale del file di Astra in input
        write(out_file(8:9),'(i2.2)') myid
        write(nml_file(8:9),'(i2.2)') myid
       else 
        N_file=N_record(trim(arg1))
        end if

   OPEN (100,FILE='analyzer_in.txt',status='old',iostat=O_stat)
   if (O_stat==0) read (100,NML=optional_par,iostat=R_stat)
   if (R_stat/=0) then
           print *,'ATTENTION:'
           print *,'The Namelist <optional_par> could be spelled incorrectly'
           print *,'Variables: freq, rotate, print_rot_files'
           print *,'keep default values: 2856.0D6, true, false'
           end if
   close (100)

   !-----alcuni preliminari----
   m2deg=(360*freq/v_luce)
   !---------------------------
   
   
   do Nfile=1,N_file
      if ((N_s==1).and.(Nfile==1))   write(*,fmt='(2x,A)') 'MODE: whole bunch analysis'
      if ((N_s/=1).and.(Nfile==1))   write(*,fmt='(2x,A)') 'MODE: Slice analysis'
      if ((arg3/="").and.(Nfile==1)) write(*,fmt='(2x,A)') 'MODE: bunch rotation by User defined angle @ prompt'
      if (mode=="listfile") then
              open (300, file=arg1,status='old')
              read (300,*) astra_file
              end if

      Npartic=N_record(trim(astra_file))
      write (*,fmt='(I3,3A,I9)')Nfile,'  Working on file: ',trim(astra_file),'   N. particles=',Npartic 
      !print *,Nfile,'Working on file = ',trim(astra_file),'   N. particles=',Npartic 
      allocate (bunch(Npartic))

      !-----costanti di normalizzazione e preliminari----
      !-----qui devo caricare il bunch da Astra dentro le matrici appana allocate-----
      open (200,file=astra_file,status='old')
      read (200,*) bunch(1)
      do i=2,Npartic
            read (200,*) bunch(i) 
            Bunch(i)%z=Bunch(i)%z+Bunch(1)%z
            Bunch(i)%pz=Bunch(i)%pz+Bunch(1)%pz
            !Bunch(i)%t=Bunch(i)%t+Bunch(1)%t !lascio i valori relativi, che tanto non tocco
            end do
      close (200)

      if (del_lost_part) then
           write (*,fmt='(2x,2(A,2x,i9))') 'Removed lost particles:'&
           &,size(bunch)-size(pack(bunch(:)%x,bunch%f2>-10)),' on',size(bunch)
           allocate (bunch_buffer(size(pack(bunch(:)%x,bunch%f2>-10))))
           bunch_buffer=pack(bunch,bunch(:)%f2>-10)
           !print *,"dim buffer and original= ",size(bunch_buffer),size(bunch)
           !bunch_buffer=bunch
           deallocate (bunch)
           allocate (bunch(size(bunch_buffer(:)%x)))
           bunch=bunch_buffer
           print *,size(bunch)
           deallocate (bunch_buffer)
           end if

      x_med= mean(bunch%x)
      y_med= mean(bunch%y)
      z_med= mean(bunch%z)
      px_med=mean(bunch%px)
      py_med=mean(bunch%py)
      pz_med=mean(bunch%pz)
      !s_pos= bunch(1)%t*1.0d-9*v_luce !questo va concolato con il beta non con v_luce, in particolare se non sono ultra relativistico
      
      !--- mi sposto in 0,0,0 per POSIZIONI E MOMENTI per RUOTARE
      bunch%x=bunch%x-x_med
      bunch%y=bunch%y-y_med
      bunch%z=bunch%z-z_med
      bunch%px=bunch%px-px_med
      bunch%py=bunch%py-py_med
      bunch%pz=bunch%pz-pz_med

      !se ruoto imposto gli angoli
      if (rotate) then 
         theta_x=-1.0d0*datan(px_med/pz_med) !angoli negativi per radrizzare
         theta_y=-1.0d0*datan(py_med/pz_med)
         if (arg3/="") then
                 theta_x=real(theta) !angoli letti nella linea di lancio
                 theta_y=aimag(theta)
                 end if
         else
         theta_x=0.0d0
         theta_y=0.0d0
         end if
        
      call xyz_rot(bunch,theta_x,theta_y) !rotazione del pacchetto

      !***rotazione dei valori medi, stessa matrice di rot. della "CALL"
      px_med_rot=px_med*dcos(theta_x)-py_med*dsin(theta_x)*dsin(theta_y)+pz_med*dsin(theta_x)*dcos(theta_y)
      py_med_rot=py_med*dcos(theta_y)+pz_med*dsin(theta_y)
      pz_med_rot=-1.0d0*px_med*dsin(theta_x)-py_med*dcos(theta_x)*dsin(theta_y)+pz_med*dcos(theta_x)*dcos(theta_y)
      !ATTENZINE: NO RUOTO I VALORI MEDI DELLE POSIZIONI. i.e. se ruoto di 90ï¿½ la z di riferimento non si trasformerï¿½
      !nella posizione di x medio. i.e. se sono sull'orbita di un LINAC, a 4m dal cathodo e ruoto, fascio rimane a 4 metri
      !dal catodo in Z, non a x_medio=4 m. cioï¿½ non ruoto i sistemi di riferimento, ma solo le posizioni. Il discorso ï¿½ diverso
      !per i momenti.
      !************************************************************************************************************************

      !***ANALISI SLICE (N_s>1) analisi FULL BUNCH (N_s==1)
      zmin_slice=minval(Bunch%z)
      dz=(maxval(Bunch%z)-zmin_slice)/REAL(N_s,8)
            
      do k=1,N_s
            if (k==N_s) last=.true. 
            if (N_s/=1) write(*,fmt='(2x,A,i3.3)') "Slice N.",k
            zmax_slice=zmin_slice+dz+epsilon(dz)
            call bunch_slicer(zmin_slice,zmax_slice,bunch,S_bunch,last) !give in output S_bunch
            zmin_slice=zmax_slice

            np=size(S_bunch)
            if (print_rot_files) then
               if (mode=="listfile") write (bunch_name(6:8),fmt='(I3.3)') Nfile
               if (mode=="one_file") write (bunch_name(6:8),fmt='(I3.3)') k
               print *,'Print rotated bunch/slice centered in pos.(0,0,0) & momenta(0,0,0) -FILENAME--> ',bunch_name,'  Part#=',np
               open (unit=500,file=bunch_name)
               write (500,fmt=Astra_fmt_d) S_bunch(1)
               write (500,fmt=Astra_fmt_d) (S_bunch(n), n=2,np)
               close (500)
               end if
            !*** dopo la rotazione riaggiungo i valori medi
            !**** POSIZIONI: tengo quelli del fascio originario 
            !***** (tipicamente il fascio dritto che voglio ruotare ma tenere nella posizione di partenza)
            !**** MOMENTI: aggiungo i valori medi dei momenti ruotati, infatti se fuoto di 90ï¿½ tutto il 
            !***** pz si deve trasformare in px
            S_bunch%x=S_bunch%x+x_med
            S_bunch%y=S_bunch%y+y_med
            S_bunch%z=S_bunch%z+z_med
            S_bunch%px=S_bunch%px+px_med_rot
            S_bunch%py=S_bunch%py+py_med_rot
            S_bunch%pz=S_bunch%pz+pz_med_rot

            x_med_rot= mean(S_bunch%x)
            y_med_rot= mean(S_bunch%y)
            z_med_rot= mean(S_bunch%z)

            if (arg3/="") then !questo serve per stampare il file ruotato con angolo richiesto dall'esterno
                  open (unit=600,file=astra_out_file)
                  write (600,fmt=Astra_fmt_d) S_bunch(1)
                  do i=2,np
                      S_Bunch(i)%z=S_Bunch(i)%z-S_Bunch(1)%z
                      S_Bunch(i)%pz=S_Bunch(i)%pz-S_Bunch(1)%pz
                      !S_Bunch(i)%t=S_Bunch(i)%t-s_pos/(1.e-9*v_luce)
                      write (600,fmt=Astra_fmt_d) S_bunch(i)
                      end do
                  close (600)
                  write(*,fmt='(2(2x,A))') "The output file is:",astra_out_file

            else 
                  !Modalitï¿½ in cui si raddrizza il fascio per analizzarlo
                  allocate (p(np),pr(np),xp(np),yp(np))
                  p=dsqrt(S_bunch%px**2+S_bunch%py**2+S_bunch%pz**2)
                  p_med=mean(p)
                  xp=(S_bunch%px/S_bunch%pz)
                  yp=(S_bunch%py/S_bunch%pz)
                  xp_med=mean(xp) 
                  yp_med=mean(yp)
             
                  g=dsqrt((p_med/m0)**2+1.0D0)
                  beta=DSQRT((g**2)-1.0d0)/g
                  s_pos= bunch(1)%t*1.0d-9*v_luce*beta !questo va concolato con il beta non con v_luce, in particolare se non sono ultra relativistico
                  x_max=MAXVAL(S_bunch%x)
                  y_max=MAXVAL(S_bunch%y)
                  E_var=0.0D0
                  
                  x_var=0.0D0  
                  y_var=0.0D0  
                  z_var=0.0D0
                  
                  xp_var=0.0D0 
                  yp_var=0.0D0

                  px_var=0.0D0
                  py_var=0.0D0                  
                  pz_var=0.0D0
             
                  z_cor=0.0D0
                  do j=1,np
                     x_var=x_var+(S_bunch(j)%x-x_med_rot)**2   !leggere app1 sopra, metodo per velocizzare
                     y_var=y_var+(S_bunch(j)%y-y_med_rot)**2
                     z_var=z_var+(S_bunch(j)%z-z_med_rot)**2  
                     xp_var=xp_var+(xp(j)-xp_med)**2
                     yp_var=yp_var+(yp(j)-yp_med)**2
                     px_var=px_var+((S_bunch(j)%px-px_med_rot)**2)
                     py_var=py_var+((S_bunch(j)%py-py_med_rot)**2)
                     pz_var=pz_var+((S_bunch(j)%pz-pz_med_rot)**2)
                     z_cor=z_cor+(S_bunch(j)%z*(S_bunch(j)%pz-pz_med_rot))
                     !**********************************************************
                     !quando il delta_p ï¿½ diviso per p_med, l'emittanza long.
                     !ï¿½ in [mm-mrad] che ï¿½ un vecchio formalismo (quello che si 
                     !trova in Trace3D premendo il tasto W
                     !z_cor=z_cor+(S_bunch(j)%z*((S_bunch(j)%pz-pz_med))/pz_med)
                     !pz_var=pz_var+(((S_bunch(j)%pz-pz_med)/pz_med)**2)  
                     !**********************************************************
                     E_var=E_var+(dsqrt(p(j)**2+m0**2)-dsqrt(p_med**2+m0**2))**2
                     end do
             
                  !dispersion evaluation
                  pr=(p-p_med)/p_med
                  !***queste varianze e valori medi, se inseriti nel ciclo sopra
                  !***molto probabilmente pormetterebbero un esecuzione parecchio piï¿½ rapida
                  eta_x=mean(S_bunch%x*pr)/var(pr)
                  eta_y=mean(S_bunch%y*pr)/var(pr)
                  etaP_x=mean(xp*pr)/var(pr)
                  etaP_y=mean(yp*pr)/var(pr)
             
                  beam_energy=m0*1.0D-6*(g-1.0d0)
                  E_var=E_var/np
                  
                  x_var=x_var/np
                  y_var=y_var/np
                  z_var=z_var/np
                  sig_z=DSQRT(z_var)
             
                  xp_var=xp_var/np
                  yp_var=yp_var/np
                  
                  px_var=px_var/np
                  py_var=py_var/np
                  pz_var=pz_var/np
                  !la covarianza come nota klauss Floettmann o Tesi PhD Marcello Ressitti Conti
                  !x_cor=mean(S_bunch%x* (S_bunch%px/(S_bunch%pz))) - (mean(S_bunch%x)*mean((S_bunch%px/(S_bunch%pz))))
                  !y_cor=mean(S_bunch%y* (S_bunch%py/(S_bunch%pz))) - (mean(S_bunch%y)*mean((S_bunch%py/(S_bunch%pz))))
                  x_cor=mean(S_bunch%x* (S_bunch%px)) - (mean(S_bunch%x)*mean((S_bunch%px)))
                  y_cor=mean(S_bunch%y* (S_bunch%py)) - (mean(S_bunch%y)*mean((S_bunch%py)))
                  z_cor=z_cor/np
             
                  !emit_x=DSQRT(x_var*xp_var-x_cor**2)
                  !emit_y=DSQRT(y_var*yp_var-y_cor**2)
                  !emit_z=DSQRT(z_var*pz_var-z_cor**2)

                  emit_x=DSQRT(x_var*px_var-x_cor**2)
                  emit_y=DSQRT(y_var*py_var-y_cor**2)
                  emit_z=DSQRT(z_var*pz_var-z_cor**2)
             
                  if (N_s.eq.1) curr=((-1.0d-9)*real(np,8)*S_bunch(1)%q*v_luce*beta)/(sig_z*DSQRT(12.0d0))
                  if (N_s.gt.1) curr=((-1.0d-9)*real(np,8)*S_bunch(1)%q*v_luce*beta)/(dz)
                  Q_perc=INT(100*np/Npartic)
                  
                  !Calcolo una volta sola le variabili per l'output usando i valori temporanei "r_XXXX"
                  r_np=np
                  r_Den=1.0d-3*dsqrt(E_var)
                  r_Xmax=x_max*1.0d3
                  r_Ymax=y_max*1.0d3
                  r_SigX=1.0d3*DSQRT(x_var)
                  r_SigY=1.0d3*DSQRT(y_var)
                  r_SigZ=sig_z*1.d3
                  r_SigPx=1.0d3*DSQRT(px_var)
                  r_SigPy=1.0d3*DSQRT(py_var)
                  r_SigPz=1.0d3*DSQRT(pz_var)
                  !r_EmitX=1.d6*beta*g*emit_x
                  !r_EmitY=1.d6*beta*g*emit_y
                  r_EmitX=(1.d6/m0)*emit_x
                  r_EmitY=(1.d6/m0)*emit_y
                  r_Curr=curr
                  r_EtaX=eta_x
                  r_EtaPx=etaP_x
                  r_EtaY=eta_y
                  r_EtaPy=etaP_y
                  r_AlphaX_Twiss=(-1.0d0)*x_cor/emit_x
                  r_AlphaY_Twiss=(-1.0d0)*y_cor/emit_y
                  r_BetaX_Twiss=x_var/emit_x
                  r_BetaY_Twiss=y_var/emit_y
                  r_GammaX_Twiss=xp_var/emit_x
                  r_GammaY_Twiss=yp_var/emit_y
                  select case (mode)
                     case ("one_file")
                        if (k.eq.1) then 
                            OPEN (1000,FILE=out_file,STATUS='replace')
                            WRITE(1000,FMT='(A7,A13,14A13)') "slice","n_partic","dz[mm]","x_max[mm]"&
                            &,"y_max[mm]","dE [keV]","sig_x[mm]","sig_y[mm]","sig_z[mm]","I_med[A]","e_x[mm*mrad]","e_y[mm*mrad]"&
                            &,"etaX","etaX_p","etaY","etaY_p"
                            open(314,file=nml_file,status='replace')
                            write(nml=OutParam, unit=314)
                        end if

                        write (1000,1112) k,r_np,'(',Q_perc,'%)',dz*1.0d3,r_Xmax,r_Ymax,r_Den,r_SigX,&
                        r_SigY,r_SigZ,r_Curr,r_EmitX,r_EmitY,r_EtaX,r_EtaPx,r_EtaY,r_EtaPy
                        if (N_s.eq.1) then
                           write (1000,*)
                           write (1000,*) 'bunch average ref. values: x_med[m]=',x_med,' y_med[m]=',y_med
                           write (1000,*)
                           write (1000,*) 'beam energy=',beam_energy,'  <P_tot> [MeV/c]= ',p_med/1.0D6
                           write (1000,*) 'unorm. emit_x,y=',r_emitx/(beta*g),', ',r_emity/(beta*g)
                           write (1000,*) 'z_emit[keV-mm]=',emit_z,'  z_emit[deg-keV]=', 1.0d-3*m2deg*emit_z
                           write (1000,*) 'Courant-Snyder parameters:'
                           write (1000,*) 'beta_x,y [mm/mrad]= ',r_BetaX_Twiss,',  ',r_BetaY_Twiss
                           write (1000,*) 'beta_z   = ',(z_var*m2deg*1.0d3)/emit_z,' [deg/KeV]'
                           write (1000,*) 'alfa_x,y = ',r_AlphaX_Twiss,',  ',r_AlphaY_Twiss
                           write (1000,*) 'alfa_z = ',(-1.0d0)*z_cor/emit_z
                           write (1000,*) 'gamma_x=',r_GammaX_Twiss,' gamma_y=',r_GammaX_Twiss
                           end if
                 
                     case ("listfile")
                        if (Nfile.eq.1) then 
                                OPEN (1001,FILE='file_analysis.txt',STATUS='replace')
                                write (1001,fmt='(13A13)') "s_pos[m]","z_lab[m]","sig_x[mm]","sig_y[mm]","sig_z[mm]","e_x[um]",&
                                "e_y[um]","Beta_x[m]","Beta_y[m]","etaX[m]","etaX_p[rad]","etaY[m]","etaY_p[rad]"
                        end if
                        write (1001,fmt='(13E13.5)') s_pos,z_med,r_SigX,r_SigY,r_SigZ,&
                        r_EmitX,r_EmitY,r_BetaX_Twiss,r_BetaY_Twiss,r_EtaX,r_EtaPx,r_EtaY,r_EtaPy
                  end select
                  deallocate (S_bunch,p,pr,xp,yp)
            end if
      end do
      deallocate (bunch)
   end do
   close (300)
   call cpu_time(tf)
   write (*,fmt='(2x,A,es11.4e2)') 'computation time [s] = ',tf-ti
   Print *,"******************************* END *********************************"
   close (1000)
   close (1001)
   1112 FORMAT(2i7,a1,i3,a2,14E13.4)
end
