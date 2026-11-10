subroutine probe(l,mbra,lp,mket,lpot,mpot,radial_geom,angular_geom)
implicit none
integer,intent(in)::l,mbra,lp,mket,lpot,mpot
complex*16,intent(out)::radial_geom(3,3),angular_geom(3,3)
integer,parameter::lmax2n=8
complex*16,parameter::zeroc=(0d0,0d0)
real*8,parameter::hf=0.5d0
integer::ih,sp,tp,t,u,v1,v,alpha,beta,lint,nu
integer::lv(1),lpv(1),llv(1),mmv(1)
complex*16::c_am(3,-1:1),ca_st_lm(3,-1:1,-1:1),radial_factor,angular_factor
complex*16::vpuu(1),vpdd(1),vpud(1),vpdu(1),wpuu(1),wpdd(1),wpud(1),wpdu(1)
logical::duplicate_radial_channel
real*8::gaunt1,gaunt_pot,gaunt_couple,channel_diff,channel_scale
real*8,external::gaunt
integer,external::notri
ih=1;lv=lp;lpv=l;llv=lpot;mmv=mpot
vpuu=1;vpdd=1;vpud=1;vpdu=1;wpuu=1;wpdd=1;wpud=1;wpdu=1
call c_alpha_m(c_am)
call c_factor(l,mbra,ca_st_lm)
radial_geom=zeroc;angular_geom=zeroc

       ! ATPAR/FORFHS writes all (m,m') replicas of one radial channel
       ! consecutively, so comparison with the preceding record is sufficient
       ! and avoids an O(ihmx**2) search at every k point.
       ! Fortran .AND. need not short circuit: guard ih-1 separately.
       IF (ih.GT.1) THEN
       IF (lpv(ih-1).EQ.l .AND. lv(ih-1).EQ.lp .AND. &
           llv(ih-1).EQ.llv(ih) .AND. mmv(ih-1).EQ.mmv(ih)) THEN
          channel_diff = ABS(vpuu(ih-1)-vpuu(ih))+ &
                         ABS(vpdd(ih-1)-vpdd(ih))+ &
                         ABS(vpud(ih-1)-vpud(ih))+ &
                         ABS(vpdu(ih-1)-vpdu(ih))+ &
                         ABS(wpuu(ih-1)-wpuu(ih))+ &
                         ABS(wpdd(ih-1)-wpdd(ih))+ &
                         ABS(wpud(ih-1)-wpud(ih))+ &
                         ABS(wpdu(ih-1)-wpdu(ih))
          channel_scale = ABS(vpuu(ih))+ABS(vpdd(ih))+ABS(vpud(ih))+ &
                          ABS(vpdu(ih))+ABS(wpuu(ih))+ABS(wpdd(ih))+ &
                          ABS(wpud(ih))+ABS(wpdu(ih))
          duplicate_radial_channel = &
               channel_diff.LE.1.d-12*(1.d0+channel_scale)
       ENDIF
       ENDIF

             DO sp = 0,2
                IF (NOTRI(sp,1,1).LT.0) CYCLE
                DO tp = -sp,sp
                   DO t = -1,1
                      DO u = -1,1
                         ! GAUNT does not reliably return zero for a
                         ! magnetically forbidden triplet.  Every historical
                         ! stress loop tests -m1+m2+m3=0 before calling it;
                         ! do the same here.  Omitting this guard fed
                         ! undefined values into the experimental G2 term.
                         IF (-tp+t+u.NE.0) CYCLE
                         gaunt1 = GAUNT(sp,1,1,tp,t,u)
                         IF (ABS(gaunt1).LT.1.d-30) CYCLE
                         ! The recoupled angular momentum is not a basis
                         ! expansion index.  It may exceed lmax2n; the exact
                         ! finite upper bound follows from L'' x l'.
                         DO v1 = 0,llv(ih)+lp
                            IF (NOTRI(v1,llv(ih),lp).LT.0) CYCLE
                            IF (NOTRI(l,sp,v1).LT.0) CYCLE
                            DO v = -v1,v1
                               IF (-v+mmv(ih)+mket.NE.0) CYCLE
                               gaunt_pot = GAUNT(v1,llv(ih),lp,v, &
                                                mmv(ih),mket)
                               IF (ABS(gaunt_pot).LT.1.d-30) CYCLE
                               IF (-mbra+tp+v.NE.0) CYCLE
                               gaunt_couple = GAUNT(l,sp,v1,mbra,tp,v)
                               IF (ABS(gaunt_couple).LT.1.d-30) CYCLE
                               DO alpha = 1,3
                                  DO beta = 1,3
                                     ! This Gaunt ordering expands n_alpha*n_beta
                                     ! on the KET side of <Y_bra|...|Y_ket>.
                                     ! Conjugating c_am here reverses xy/yz signs.
                                     ! The spherical block uses another ordering
                                     ! and must NOT be changed in the same way.
                                     radial_factor = &
                                          c_am(alpha,t)* &
                                          c_am(beta,u)*gaunt1* &
                                          gaunt_pot*gaunt_couple
                                     radial_geom(alpha,beta) = &
                                          radial_geom(alpha,beta)+radial_factor
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO

             DO v = -1,1,2
                ! grad_Omega Y_l contains l+1 even for l=lmax2n.
                ! This is an intermediate angular harmonic, not an extra
                ! radial basis channel.  Cutting it spoils the exact zero
                ! trace n dot grad_Omega=0, including in TRACE mode.
                IF (l+v.LT.0) CYCLE
                DO t = -1,1
                   IF (mbra+t.LT.-(l+v) .OR. &
                       mbra+t.GT. (l+v)) CYCLE
                   DO tp = -1,1
                      ! Likewise this intermediate momentum belongs to the
                      ! four-harmonic product, not to the LAPW expansion.
                      ! Truncating it at lmax2n violates the hydrostatic trace
                      ! identity of the angular derivative.
                      DO lint = 0,llv(ih)+lp
                         IF (NOTRI(lint,llv(ih),lp).LT.0) CYCLE
                         IF (NOTRI(lint,1,l+v).LT.0) CYCLE
                         DO nu = -lint,lint
                            IF (-nu+mmv(ih)+mket.NE.0) CYCLE
                            gaunt_pot = GAUNT(lint,llv(ih),lp,nu, &
                                             mmv(ih),mket)
                            IF (ABS(gaunt_pot).LT.1.d-30) CYCLE
                            IF (-nu+tp+mbra+t.NE.0) CYCLE
                            gaunt_couple = GAUNT(lint,1,l+v,nu,tp, &
                                                mbra+t)
                            IF (ABS(gaunt_couple).LT.1.d-30) CYCLE
                            DO alpha = 1,3
                               DO beta = 1,3
                                  angular_factor = hf*( &
                                       conjg(c_am(alpha,tp))* &
                                       conjg(ca_st_lm(beta,v,t)) + &
                                       conjg(c_am(beta,tp))* &
                                       conjg(ca_st_lm(alpha,v,t)) )
                                  angular_geom(alpha,beta) = &
                                       angular_geom(alpha,beta) + &
                                       angular_factor*gaunt_pot*gaunt_couple
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
end subroutine
