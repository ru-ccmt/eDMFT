! Minimal data modules for the unmodified native reference routines.
module TestCases
logical :: SimpInt=.false., ConOpt=.true.
end module
module struct
real(8) :: br1(3,3),vi
end module
module varyatoms
real(8) :: ztot,nec1,nec2
logical :: useatoms=.false.
end module
