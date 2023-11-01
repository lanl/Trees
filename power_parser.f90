!----------------------------------------------------------------
! Queries gridlist for input reals
!----------------------------------------------------------------
subroutine QueryFuellist_real(variableName,variableValue, &
    fileunit,variableDefault)
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(in)  :: fileunit
  real,intent(in)  :: variableDefault
  real,intent(out) :: variableValue

  integer :: ierror
  character(len=1)    :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  variableValue=variableDefault

  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    read(text,*) word
    if(word.eq.variableName)then
      read(text,*) word,equal,variableValue
      exit
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryFuellist_real

!----------------------------------------------------------------
! Queries gridlist for input real arrays
!----------------------------------------------------------------
subroutine QueryFuellist_real_array(variableName,array, &
    arraysize,fileunit,variableDefault)
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(in)  :: arraysize
  integer,intent(in)  :: fileunit
  real,intent(in)  :: variableDefault
  real,intent(out) :: array(arraysize)

  integer :: ierror
  character(len=1)    :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  array(:)=variableDefault

  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    read(text,*) word
    if(word.eq.variableName)then
      read(text,*) word,equal,array
      exit
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryFuellist_real_array

!----------------------------------------------------------------
! Queries gridlist for input integers
!----------------------------------------------------------------
subroutine QueryFuellist_integer(variableName,variableValue, &
    fileunit,variableDefault)
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(in)  :: fileunit
  integer,intent(in)  :: variableDefault
  integer,intent(out) :: variableValue

  integer :: ierror
  character(len=1)    :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  variableValue=variableDefault

  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    read(text,*) word
    if(word.eq.variableName)then
      read(text,*) word,equal,variableValue
      exit
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryFuellist_integer

!----------------------------------------------------------------
! Queries gridlist for input integer arrays
!----------------------------------------------------------------
subroutine QueryFuellist_integer_array(variableName,array, &
    arraysize,fileunit,variableDefault)
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(in)  :: arraysize
  integer,intent(in)  :: fileunit
  integer,intent(in)  :: variableDefault
  integer,intent(out) :: array(arraysize)

  integer :: ierror
  character(len=1)    :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  array(:)=variableDefault

  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    read(text,*) word
    if(word.eq.variableName)then
      read(text,*) word,equal,array
      exit
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryFuellist_integer_array

!----------------------------------------------------------------
! Queries gridlist for input string arrays
!----------------------------------------------------------------
subroutine QueryFuellist_string(variableName,variableValue, &
    fileunit,variableDefault)
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(in)  :: fileunit
  character(len=*),intent(in)  :: variableDefault
  character(len=*),intent(out) :: variableValue

  integer :: i,ierror
  character(len=20)   :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  variableValue=variableDefault

  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    read(text,*) word
    if(word.eq.variableName)then
      read(text,*) word,equal,variableValue
      exit
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryFuellist_string

!----------------------------------------------------------------
! Queries gridlist for input string arrays
!----------------------------------------------------------------
subroutine QueryFuellist_string_array(variableName,array, &
    arraysize,fileunit,variableDefault)
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(in)  :: arraysize
  integer,intent(in)  :: fileunit
  character(len=*),intent(in)  :: variableDefault
  character(len=*),intent(out) :: array(arraysize)

  integer :: i,ierror
  character(len=20)   :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  array(:)=variableDefault

  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    read(text,*) word
    if(word.eq.variableName)then
      read(text,*) word,equal,array
      exit
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryFuellist_string_array
