module File_analysis
      implicit none
      contains
        function N_record(nomefile) result(Nrec)
           implicit none
           integer :: Nrec
           character(len=*),intent(in) :: nomefile
           integer ::stat_100
           
           Nrec=0
           open (100,file=nomefile,status='old')
           Nrec=0
           do
              read (100,*,iostat=stat_100)
              Nrec=Nrec+1
              if (stat_100==-1) exit
           end do
           Nrec=Nrec-1
           close(100)
        end function N_record

        function N_columns(nomefile) result(Ncol)
           implicit none
           integer:: Ncol
           character(len=*),intent(in) :: nomefile
           integer :: stat_100
           CHARACTER :: OldCh, Ch
           CHARACTER :: BL = " "

           stat_100=0
           Ncol=0
           open (100,file=nomefile,status='old')
           OldCh = BL
           DO WHILE (stat_100/=-2) ! check for EOR
              READ (100,'(A1)',iostat=stat_100,ADVANCE='NO') Ch
              IF (stat_100==0) THEN ! protect against EOR
                 IF ((Ch==BL.or.Ch==achar(9)).AND.(OldCh/=BL.or.OldCh/=achar(9))) THEN !end of column
                    Ncol=Ncol+1 ! ... signals end of word
                 END IF
                 OldCh=Ch
              END IF
           END DO
           IF (OldCh /= BL.or.OldCh/=achar(9)) THEN ! if last char actually read is non-blank ..
              Ncol=Ncol+1 ! ... count another word
           END IF
           close (100)
           end function N_columns
        end module File_analysis
