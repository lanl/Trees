module QueryGridlist_mod
  use iso_fortran_env, only : int64, real64
  Implicit None

  contains
    
    !-------------------------------------------------------------------
    ! Queries gridlist for input variable values regardless of type
    !-------------------------------------------------------------------
    subroutine QueryFuellist_variable(variableName,variableValue, &
        fileunit,variableDefault)
      Implicit None
    
      ! Local Variables
      character(len=*),intent(in) :: variableName
      class(*),intent(inout) :: variableValue
      class(*),intent(in) :: variableDefault
      integer,intent(in)  :: fileunit
    
      integer :: ierror,eqpos,nfound
      character(len=20)   :: word
      character(len=1000) :: queryName
      character(len=1000) :: text,valueText
    
      ! Executable Code
      ! Set default value using select type
      select type (x=>variableValue)
        type is (integer)
          select type (y=>variableDefault)
            type is (integer)
              x = y
          end select
        type is (integer(int64))
          select type (y=>variableDefault)
            type is (integer(int64))
              x = y
          end select
        type is (character(len=*))
          select type (y=>variableDefault)
            type is (character(len=*))
              x = y
          end select
      end select
      
      queryName = variableName
      call lowercase_inplace(queryName)
      nfound = 0
      do ! Iterate through all lines of file
        read(fileunit,"(a)",iostat=ierror) text ! Read line into text
        if (ierror/=0) exit
        if(text/='')then
          eqpos=index(text,'=')
          word =trim(adjustl(text(:eqpos-1)))
          call lowercase_inplace(word)
          if(trim(word).eq.trim(queryName))then
            nfound = nfound + 1
            if(nfound.gt.1)then
              print*,'ERROR: duplicate gridlist key: ', &
                trim(variableName)
              STOP
            endif
            valueText=trim(text(eqpos+1:))
            select type (x=>variableValue)
              type is (integer)
                read(valueText,*) x
              type is (integer(int64))
                read(valueText,*) x
              type is (character(len=*))
                read(valueText,*) x
              class default
                print*,"Unexpected type in reading fuellist"
            end select
          endif
        endif
      enddo
      rewind(fileunit)

    end subroutine QueryFuellist_variable
    
    !-------------------------------------------------------------------
    ! Queries gridlist for input array values regardless of type
    !-------------------------------------------------------------------
    subroutine QueryFuellist_array(variableName,arrayValue,size, &
        fileunit,variableDefault)
      Implicit None

      ! Local Variables
      character(len=*),intent(in) :: variableName
      integer,intent(in)  :: size
      integer,intent(in)  :: fileunit
      class(*),intent(inout) :: arrayValue(size)
      class(*),intent(in) :: variableDefault

      integer :: ierror,eqpos,nfound
      character(len=20)   :: word
      character(len=1000) :: queryName
      character(len=1000) :: text,valueText

      ! Executable Code
      ! Set default value using select type
      select type (x=>arrayValue)
        type is (integer)
          select type (y=>variableDefault)
            type is (integer)
              x(:) = y
          end select
        type is (integer(int64))
          select type (y=>variableDefault)
            type is (integer(int64))
              x(:) = y
          end select
        type is (character(len=*))
          select type (y=>variableDefault)
            type is (character(len=*))
              x(:) = y
          end select
      end select
      
      queryName = variableName
      call lowercase_inplace(queryName)
      nfound = 0
      do ! Iterate through all lines of file
        read(fileunit,"(a)",iostat=ierror) text ! Read line into text
        if (ierror/=0) exit
        if(text/='')then
          eqpos=index(text,'=')
          word =trim(adjustl(text(:eqpos-1)))
          call lowercase_inplace(word)
          if(trim(word).eq.trim(queryName))then
            nfound = nfound + 1
            if(nfound.gt.1)then
              print*,'ERROR: duplicate gridlist key: ', &
                trim(variableName)
              STOP
            endif
            valueText=trim(text(eqpos+1:))
            select type (x=>arrayValue)
              type is (integer)
                read(valueText,*) x
              type is (integer(int64))
                read(valueText,*) x
              type is (character(len=*))
                read(valueText,*) x
              class default
                print*,"Unexpected type in reading fuellist"
            end select
          endif
        endif
      enddo
      rewind(fileunit)
    end subroutine QueryFuellist_array

    !-------------------------------------------------------------------
    ! Modifies text to include both uppercase and lowercase characters
    ! in string
    !-------------------------------------------------------------------
    subroutine lowercase_inplace(text)
      Implicit None

      ! Local Variables
      character(len=*),intent(inout) :: text
      integer :: i,c

      ! Executable Code
      do i=1,len_trim(text)
        c = iachar(text(i:i))
        if(c.ge.iachar('A').and.c.le.iachar('Z'))then
          text(i:i) = achar(c + iachar('a') - iachar('A'))
        endif
      enddo
    end subroutine lowercase_inplace

end module QueryGridlist_mod
