	module reallocate
	  !     nur 1 (generischer) Name wird von aussen angesprochen
	  interface doreallocate
	    module procedure doreallocate_r8_d1
	    module procedure doreallocate_r8_d2
	    module procedure doreallocate_c16_d1
	    module procedure doreallocate_c16_d2
	    module procedure doreallocate_i4_d1
	    module procedure doreallocate_i4_d2
	    module procedure hugo     !   ;)
	  end interface
	contains

	  !     leider sind mehrere subroutines notwendig fuer verschiedene Typen
	  subroutine doreallocate_r8_d1(tf, newdimension)
	    real*8, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
	    !     nur 1 mal kopieren reicht
	    !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    !     der Zeiger wird nur auf das neue Feld umgebogen, nicht neu alloziert
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_r8_d2(tf, newdimension1, newdimension2)
	    real*8, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_c16_d1(tf, newdimension)
	    complex*16, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_c16_d2(tf, newdimension1, newdimension2)
	    complex*16, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d1(tf, newdimension)
	    integer*4, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d2(tf, newdimension1, newdimension2)
	    integer*4, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	!     Es gibt auch Methoden, um das Programm unleserlich zu machen :-)
	!     das sollten wir besser vermeiden ;-)
	  subroutine hugo
	    write(*,*) " Hier ist Hugo"
	  end subroutine
	end module 
