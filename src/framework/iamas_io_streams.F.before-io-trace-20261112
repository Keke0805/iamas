! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module iamas_io_streams

#define COMMA ,
#define STREAM_DEBUG_WRITE(M) call iamas_log_write(M)

#ifdef SINGLE_PRECISION
#define REAL_IO_TYPE iAMAS_IO_REAL
#else
#define REAL_IO_TYPE iAMAS_IO_DOUBLE
#endif

   use iamas_attlist
   use iamas_derived_types
   use iamas_pool_routines
   use iamas_kind_types, only : StrKIND
   use iamas_timekeeping
   use iamas_io
   use iamas_log
   use iamas_dmpar

   interface iAMAS_streamAddField
      module procedure iAMAS_streamAddField_0dInteger
      module procedure iAMAS_streamAddField_1dInteger
      module procedure iAMAS_streamAddField_2dInteger
      module procedure iAMAS_streamAddField_3dInteger
      module procedure iAMAS_streamAddField_0dReal
      module procedure iAMAS_streamAddField_1dReal
      module procedure iAMAS_streamAddField_2dReal
      module procedure iAMAS_streamAddField_3dReal
      module procedure iAMAS_streamAddField_4dReal
      module procedure iAMAS_streamAddField_5dReal
      module procedure iAMAS_streamAddField_0dChar
      module procedure iAMAS_streamAddField_1dChar
   end interface iAMAS_streamAddField

   ! For multiple files @hxy
   interface iAMAS_streamAddField_mul
      module procedure iAMAS_streamAddField_0dInteger_mul
      module procedure iAMAS_streamAddField_1dInteger_mul
      module procedure iAMAS_streamAddField_2dInteger_mul
      module procedure iAMAS_streamAddField_3dInteger_mul
      module procedure iAMAS_streamAddField_0dReal_mul
      module procedure iAMAS_streamAddField_1dReal_mul
      module procedure iAMAS_streamAddField_2dReal_mul
      module procedure iAMAS_streamAddField_3dReal_mul
      module procedure iAMAS_streamAddField_4dReal_mul
      module procedure iAMAS_streamAddField_5dReal_mul
      module procedure iAMAS_streamAddField_0dChar_mul
      module procedure iAMAS_streamAddField_1dChar_mul
   end interface iAMAS_streamAddField_mul
   
   interface iAMAS_streamUpdateField
      module procedure iAMAS_streamUpdateField_0dInteger
      module procedure iAMAS_streamUpdateField_1dInteger
      module procedure iAMAS_streamUpdateField_2dInteger
      module procedure iAMAS_streamUpdateField_3dInteger
      module procedure iAMAS_streamUpdateField_0dReal
      module procedure iAMAS_streamUpdateField_1dReal
      module procedure iAMAS_streamUpdateField_2dReal
      module procedure iAMAS_streamUpdateField_3dReal
      module procedure iAMAS_streamUpdateField_4dReal
      module procedure iAMAS_streamUpdateField_5dReal
      module procedure iAMAS_streamUpdateField_0dChar
      module procedure iAMAS_streamUpdateField_1dChar
   end interface iAMAS_streamUpdateField

   interface iAMAS_readStreamAtt
      module procedure iAMAS_readStreamAtt_0dInteger
      module procedure iAMAS_readStreamAtt_1dInteger
      module procedure iAMAS_readStreamAtt_0dReal
      module procedure iAMAS_readStreamAtt_1dReal
      module procedure iAMAS_readStreamAtt_text
   end interface iAMAS_readStreamAtt

   interface iAMAS_writeStreamAtt
      module procedure iAMAS_writeStreamAtt_0dInteger
      module procedure iAMAS_writeStreamAtt_1dInteger
      module procedure iAMAS_writeStreamAtt_0dReal
      module procedure iAMAS_writeStreamAtt_1dReal
      module procedure iAMAS_writeStreamAtt_text
   end interface iAMAS_writeStreamAtt


   private mergeArrays


   contains


   subroutine iAMAS_createStream(stream, ioContext, fileName, ioFormat, ioDirection, precision, &
                                clobberRecords, clobberFiles, truncateFiles, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(out) :: stream
      type (iamas_io_context_type), pointer :: ioContext
      character (len=*), intent(in) :: fileName
      integer, intent(in) :: ioFormat
      integer, intent(in) :: ioDirection
      integer, intent(in), optional :: precision
      logical, intent(in), optional :: clobberRecords
      logical, intent(in), optional :: clobberFiles
      logical, intent(in), optional :: truncateFiles
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR


      stream % fileHandle = iAMAS_io_open(fileName, ioDirection, ioFormat, ioContext, clobber_file=clobberFiles, truncate_file=truncateFiles, &
                                         ierr=io_err) 
      !
      ! Catch a few special errors
      !
      if (io_err == iAMAS_IO_ERR_NOEXIST_READ) then
         if (present(ierr)) ierr = iAMAS_IO_ERR_NOEXIST_READ
         return
      end if
      if (io_err == iAMAS_IO_ERR_WOULD_CLOBBER) then
         if (present(ierr)) ierr = iAMAS_STREAM_CLOBBER_FILE
         return
      end if

      ! General error
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

      stream % ioDirection = ioDirection
      stream % ioFormat = ioFormat
      stream % filename = fileName
      if (present(clobberRecords)) then
         stream % clobberRecords = clobberRecords
      else
         stream % clobberRecords = .false.
      end if
      if (present(precision)) stream % defaultPrecision = precision

      stream % isInitialized = .true.

   end subroutine iAMAS_createStream

   ! For multiple file @ftao
   subroutine iAMAS_createStream_mul(stream, ioContext, fileName, ioFormat, ioDirection, precision, &
                                clobberRecords, clobberFiles, truncateFiles, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(out) :: stream
      type (iamas_io_context_type), pointer :: ioContext
      character (len=*), intent(in) :: fileName
      integer, intent(in) :: ioFormat
      integer, intent(in) :: ioDirection
      integer, intent(in), optional :: precision
      logical, intent(in), optional :: clobberRecords
      logical, intent(in), optional :: clobberFiles
      logical, intent(in), optional :: truncateFiles
      integer, intent(out), optional :: ierr
      type(iAMAS_IO_Handle_type) :: temp_Handle

      integer :: io_err

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR


      ! if (ioDirection == iAMAS_IO_WRITE) then
      ! temp_Handle = iAMAS_io_open_mul(fileName, ioDirection, ioFormat, ioContext, clobber_file=clobberFiles, truncate_file=truncateFiles, &
      !                                    ierr=io_err)
      ! endif
      stream % fileHandle = iAMAS_io_open_mul(fileName, ioDirection, ioFormat, ioContext, clobber_file=clobberFiles, truncate_file=truncateFiles, &
                                         ierr=io_err) 
      !
      ! Catch a few special errors
      !
      if (io_err == iAMAS_IO_ERR_NOEXIST_READ) then
         if (present(ierr)) ierr = iAMAS_IO_ERR_NOEXIST_READ
         return
      end if
      if (io_err == iAMAS_IO_ERR_WOULD_CLOBBER) then
         if (present(ierr)) ierr = iAMAS_STREAM_CLOBBER_FILE
         return
      end if

      ! General error
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

      stream % ioDirection = ioDirection
      stream % ioFormat = ioFormat
      stream % filename = fileName
      ! call iamas_log_write("createStream_mul === fileName = "//trim(fileName))
      if (present(clobberRecords)) then
         stream % clobberRecords = clobberRecords
      else
         stream % clobberRecords = .false.
      end if
      if (present(precision)) stream % defaultPrecision = precision

      stream % isInitialized = .true.

   end subroutine iAMAS_createStream_mul

   integer function iAMAS_seekStream(stream, seekTime, seekPosition, actualTime, maxRecords, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: seekTime
      integer, intent(in) :: seekPosition
      character (len=*), intent(out) :: actualTime
      integer, intent(out), optional :: maxRecords
      integer, intent(out), optional :: ierr
      
      integer :: io_err
      integer :: i
      integer :: timeDim
      character (len=ShortStrKIND), dimension(:), pointer :: xtimes
      character (len=ShortStrKIND) :: strTemp
      type (iAMAS_Time_type) :: sliceTime, startTime
      type (iAMAS_TimeInterval_type) :: timeDiff, minTimeDiff

      character (len=ShortStrKIND) :: xtime0, xtime1, xtime2, xtimeGuess
      type (iAMAS_Time_type) :: time0, time1, time2, timeGuess, timeGuessData
      type (iAMAS_TimeInterval_type) :: timeInterval

!      call iamas_log_write('Called iAMAS_seekStream')

      !
      ! Initialize output arguments
      !
      if (present(ierr)) ierr = iAMAS_STREAM_NOERR
      iAMAS_seekStream = 0

      if (present(maxRecords)) maxRecords = 0

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

      call iAMAS_io_inq_dim(stream % fileHandle, 'Time', timeDim, io_err)
      if (timeDim <= 0 .or. io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

      if (present(maxRecords)) maxRecords = timeDim

!call iamas_log_write('Found $i times in file', intArgs=(/timeDim/))

      call iAMAS_io_inq_var(stream % fileHandle, 'xtime', ierr=io_err)
      if (io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         if (present(maxRecords)) maxRecords = 0
         return
      end if

!call iamas_log_write('Found xtime variable')

      ! the following section of code improves performance for forcing files
      ! that have large numbers of evenly spaced times by reducing the number
      ! of PIO character string reads needed. The following section of code is
      ! only exercised for an extact time whence in the stream manager routine.
      ! In that case the first two times in the stream file are read, their
      ! interval calculated and the desired time index inferred. This time is
      ! read and checked to be identical with the desired time. This way for
      ! evenly spaced input files with n times only 3 reads instead of n is needed.
      ! Should this method fail to produce the right time index or the seek is
      ! desired for a non-exact whence the original method is used. The limit of
      ! timeDIm > 20 ensures optimal performance for the sea ice core by avoiding
      ! using the new method for monthly forcing files which dont have even time
      ! spacing.
      if (seekPosition == iAMAS_STREAM_EXACT_TIME .and. timeDim > 20) then

         ! first time
         call iAMAS_io_set_frame(stream % fileHandle, 1, io_err)
         call iAMAS_io_get_var(stream % fileHandle, 'xtime', xtime1, io_err)
         call iamas_set_time(curr_time=time1, dateTimeString=trim(xtime1(1:32)))

         ! second time
         call iAMAS_io_set_frame(stream % fileHandle, 2, io_err)
         call iAMAS_io_get_var(stream % fileHandle, 'xtime', xtime2, io_err)
         call iamas_set_time(curr_time=time2, dateTimeString=trim(xtime2(1:32)))

         ! get interval of first two times
         timeInterval = time2 - time1

         ! get time we desire to find
         call iamas_set_time(curr_time=time0, dateTimeString=trim(seekTime))

         ! guess assuming constant intervals what index we want
         timeGuess = time1
         do i = 1, timeDim

            if (timeGuess == time0) then

               call iAMAS_io_set_frame(stream % fileHandle, i, io_err)
               call iAMAS_io_get_var(stream % fileHandle, 'xtime', xtimeGuess, io_err)

               call iamas_set_time(curr_time=timeGuessData, dateTimeString=trim(xtimeGuess(1:32)))
               if (timeGuessData == timeGuess) then
                  ! found the correct time
                  iAMAS_seekStream = i
                  write(actualTime, '(a)') trim(xtimeGuess(1:32))
                  return
               endif

            endif

            timeGuess = timeGuess + timeInterval
         enddo ! i

      endif

      allocate(xtimes(timeDim))

      do i=1,timeDim
         call iAMAS_io_set_frame(stream % fileHandle, i, io_err)
         call iAMAS_io_get_var(stream % fileHandle, 'xtime', xtimes(i), io_err)
!call iamas_log_write('... just read in xtime='//xtimes(i))
      end do

      if (len(seekTime) > 32) then
         write(strTemp, '(a)') seekTime(1:32)
      else
         write(strTemp, '(a)') trim(seekTime)
      end if

      call iamas_set_timeInterval(interval=minTimeDiff, DD=10000)
      call iamas_set_time(curr_time=startTime, dateTimeString=strTemp)

      do i=1,timeDim
         write(strTemp, '(a)') trim(xtimes(i)(1:32))
!call iamas_log_write('... converted strTemp='//strTemp)
         call iamas_set_time(curr_time=sliceTime, dateTimeString=strTemp)
         if (seekPosition == iAMAS_STREAM_EXACT_TIME) then
            if (sliceTime == startTime) then
               minTimeDiff = timeDiff
               iAMAS_seekStream = i
            end if
         else if (seekPosition == iAMAS_STREAM_NEAREST) then
            timeDiff = abs(sliceTime - startTime)
            if (timeDiff < minTimeDiff) then
               minTimeDiff = timeDiff
               iAMAS_seekStream = i
            end if
         else if (seekPosition == iAMAS_STREAM_LATEST_BEFORE) then
            if (sliceTime <= startTime) then
               timeDiff = abs(sliceTime - startTime)
               if (timeDiff < minTimeDiff) then
                  minTimeDiff = timeDiff
                  iAMAS_seekStream = i
               end if
            end if
         else if (seekPosition == iAMAS_STREAM_EARLIEST_AFTER) then
            if (sliceTime >= startTime) then
               timeDiff = abs(sliceTime - startTime)
               if (timeDiff < minTimeDiff) then
                  minTimeDiff = timeDiff
                  iAMAS_seekStream = i
               end if
            end if
         else if (seekPosition == iAMAS_STREAM_LATEST_STRICTLY_BEFORE) then
            if (sliceTime < startTime) then
               timeDiff = abs(sliceTime - startTime)
               if (timeDiff < minTimeDiff) then
                  minTimeDiff = timeDiff
                  iAMAS_seekStream = i
               end if
            end if
         else if (seekPosition == iAMAS_STREAM_EARLIEST_STRICTLY_AFTER) then
            if (sliceTime > startTime) then
               timeDiff = abs(sliceTime - startTime)
               if (timeDiff < minTimeDiff) then
                  minTimeDiff = timeDiff
                  iAMAS_seekStream = i
               end if
            end if
         else
            call iamas_log_write('Error in iAMAS_seekStream: unrecognized seekPosition', iAMAS_LOG_ERR)
            deallocate(xtimes)
            if (present(ierr)) ierr = iAMAS_IO_ERR
            return
         end if
      end do

      if (iAMAS_seekStream == 0) then
         deallocate(xtimes)
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

      write(actualTime, '(a)') trim(xtimes(iAMAS_seekStream)(1:32))

      deallocate(xtimes)

   end function iAMAS_seekStream


   subroutine iAMAS_streamTime(stream, frame, frameValidTime, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      integer, intent(in) :: frame
      character (len=*), intent(out) :: frameValidTime
      integer, intent(out), optional :: ierr
      
      integer :: io_err
      integer :: timeDim


!      call iamas_log_write('Called iAMAS_streamTime')

      !
      ! Initialize output arguments
      !
      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

      call iAMAS_io_inq_var(stream % fileHandle, 'xtime', ierr=io_err)
      if (io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

!call iamas_log_write('Found xtime variable')

      call iAMAS_io_set_frame(stream % fileHandle, frame, io_err)
      call iAMAS_io_get_var(stream % fileHandle, 'xtime', frameValidTime, io_err)
      if (io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

   end subroutine iAMAS_streamTime

   subroutine iAMAS_streamAddField_0dInteger_mul(stream, field, ierr)
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field0DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field0dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(:), pointer :: dimNames
      integer, dimension(:), pointer :: indices
      integer, dimension(:), pointer :: dimSizes
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = 0

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(0))
      allocate(dimSizes(0))
      allocate(dimNames(0))
      localDimSize = 0
      totalLocalDimSize = 0

      call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), iAMAS_IO_INT, dimNames, dimSizes, &
                                       field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, ierr=io_err)

      deallocate(indices)
      deallocate(dimSizes)
      deallocate(dimNames)
      if (io_err /= iAMAS_STREAM_NOERR) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_0D_INT
      new_field_list_node % int0dField => field

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_0dInteger_mul

   subroutine iAMAS_streamAddField_0dInteger(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field0DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field0dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(:), pointer :: dimNames
      integer, dimension(:), pointer :: indices
      integer, dimension(:), pointer :: dimSizes
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = 0

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(0))
      allocate(dimSizes(0))
      allocate(dimNames(0))
      globalDimSize = 0
      totalDimSize = 0

      call iAMAS_streamAddField_generic(stream, trim(field % fieldName), iAMAS_IO_INT, dimNames, dimSizes, &
                                       field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)

      deallocate(indices)
      deallocate(dimSizes)
      deallocate(dimNames)
      if (io_err /= iAMAS_STREAM_NOERR) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_0D_INT
      new_field_list_node % int0dField => field

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_0dInteger

   subroutine iAMAS_streamAddField_1dInteger_mul(stream, field, ierr)

      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field1DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr
      
      integer :: mpi_ierr
      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field1dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      character (len=StrKIND), dimension(0) :: dimNames0
      integer, dimension(0) :: dimSizes0
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if(.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"
      
      ! if (trim(field%fieldName) == 'indexToCellID') then
      !    call iamas_log_write('indexToCellID totalDimSize == $i ', intArgs=(/totalDimSize/))
      ! end if
      ! if (trim(field%fieldName) == 'indexToEdgeID') then
      !    call iamas_log_write('indexToEdgeID totalDimSize== $i ', intArgs=(/totalDimSize/))
      ! end if
      ! if (trim(field%fieldName) == 'indexToVertexID') then
      !    call iamas_log_write('indexToVertexID totalDimSize== $i ', intArgs=(/totalDimSize/))
      ! end if
      
      
      ! if(trim(field % fieldName)=='indexToCellID') then
      !     do i=1, totalDimSize
      !        call iamas_log_write('indexToCellID == $i ', intArgs=(/indices(i)/) ) 
      !     end do
      ! end if
      ! if(trim(field % fieldName)=='indexToEdgeID') then
      !     do i=1, totalDimSize
      !        call iamas_log_write('indexToEdgeID == $i ', intArgs=(/indices(i)/) ) 
      !     end do
      ! end if
      ! if(trim(field % fieldName)=='indexToVertexID') then
      !     do i=1, totalDimSize
      !        call iamas_log_write('indexToVertexID == $i ', intArgs=(/indices(i)/) ) 
      !     end do
      ! end if
      

      localDimSize = totalDimSize
      if (field % isDecomposed) then
         call MPI_Allreduce(localDimSize, totalLocalDimSize, 1, MPI_INTEGER, MPI_SUM, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
      else
         totalLocalDimSize = localDimSize
      end if

      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), iAMAS_IO_INT, dimNames0, &
                                             dimSizes0, field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                             ierr=io_err)
            if (io_err == iAMAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), iAMAS_IO_INT, field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_1D_INT
      new_field_list_node % int1dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_1dInteger_mul


   subroutine iAMAS_streamAddField_1dInteger(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field1DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field1dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      character (len=StrKIND), dimension(0) :: dimNames0
      integer, dimension(0) :: dimSizes0
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if(.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), iAMAS_IO_INT, dimNames0, &
                                             dimSizes0, field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                             indices, ierr=io_err)
            if (io_err == iAMAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), iAMAS_IO_INT, field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_1D_INT
      new_field_list_node % int1dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_1dInteger

   subroutine iAMAS_streamAddField_2dInteger_mul(stream, field, ierr)

      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field2DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: mpi_ierr
      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field2dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if(.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"
      localDimSize = totalDimSize
      if (field % isDecomposed) then
         call MPI_Allreduce(localDimSize, totalLocalDimSize, 1, MPI_INTEGER, MPI_SUM, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
      else
         totalLocalDimSize = localDimSize
      end if

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), iAMAS_IO_INT, field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, localDimSize,totalLocalDimSize, &
                                             ierr=io_err)
            if (io_err == iAMAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), iAMAS_IO_INT, field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_2D_INT
      new_field_list_node % int2dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_2dInteger_mul


   subroutine iAMAS_streamAddField_2dInteger(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field2DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field2dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if(.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), iAMAS_IO_INT, field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                             indices, ierr=io_err)
            if (io_err == iAMAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), iAMAS_IO_INT, field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_2D_INT
      new_field_list_node % int2dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_2dInteger

   subroutine iAMAS_streamAddField_3dInteger_mul(stream, field, ierr)

      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field3DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: mpi_ierr
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field3dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if(.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      localDimSize = totalDimSize
      if (field % isDecomposed) then
         call MPI_Allreduce(localDimSize, totalLocalDimSize, 1, MPI_INTEGER, MPI_SUM, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
      else
         totalLocalDimSize = localDimSize
      end if
      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), iAMAS_IO_INT, field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, localDimSize,totalLocalDimSize, &
                                             ierr=io_err)
            if (io_err == iAMAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), iAMAS_IO_INT, field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_3D_INT
      new_field_list_node % int3dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_3dInteger_mul

   subroutine iAMAS_streamAddField_3dInteger(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field3DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field3dInteger), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if(.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"
      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), iAMAS_IO_INT, field % dimNames(2:ndims), &
                                             field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                             indices, ierr=io_err)
            if (io_err == iAMAS_STREAM_NOERR) then
               isAvailable(i) = .true.
               any_success = .true.
            end if
         end do
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), iAMAS_IO_INT, field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_3D_INT
      new_field_list_node % int3dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_3dInteger

   subroutine iAMAS_streamAddField_0dReal_mul(stream, field, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field0DReal), intent(in), target :: field
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field0dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(:), pointer :: dimNames
      integer, dimension(:), pointer :: indices
      integer, dimension(:), pointer :: dimSizes
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      integer :: local_precision

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = 0

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(0))
      allocate(dimSizes(0))
      allocate(dimNames(0))
      globalDimSize = 0
      totalDimSize = 0

      localDimSize = 0
      totalLocalDimSize = 0

      call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), REAL_IO_TYPE , dimNames, dimSizes, &
                                       field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, precision=local_precision, ierr=io_err)

      deallocate(indices)
      deallocate(dimSizes)
      deallocate(dimNames)
      if (io_err /= iAMAS_STREAM_NOERR) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_0D_REAL
      new_field_list_node % real0dField => field

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_0dReal_mul

   subroutine iAMAS_streamAddField_0dReal(stream, field, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field0DReal), intent(in), target :: field
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field0dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(:), pointer :: dimNames
      integer, dimension(:), pointer :: indices
      integer, dimension(:), pointer :: dimSizes
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      integer :: local_precision

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = 0

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(0))
      allocate(dimSizes(0))
      allocate(dimNames(0))
      globalDimSize = 0
      totalDimSize = 0

      call iAMAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , dimNames, dimSizes, &
                                       field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)

      deallocate(indices)
      deallocate(dimSizes)
      deallocate(dimNames)
      if (io_err /= iAMAS_STREAM_NOERR) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_0D_REAL
      new_field_list_node % real0dField => field

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_0dReal

   subroutine iAMAS_streamAddField_1dReal_mul(stream, field, arrayFieldList, precision, ierr)

      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field1DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: mpi_ierr
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field1dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      character (len=StrKIND), dimension(0) :: dimNames0
      integer, dimension(0) :: dimSizes0
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"
      localDimSize = totalDimSize
      if (field % isDecomposed) then
         call MPI_Allreduce(localDimSize, totalLocalDimSize, 1, MPI_INTEGER, MPI_SUM, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
      else
         totalLocalDimSize = localDimSize
      end if
!call iamas_log_write('... local dim size is $i ', intArgs=(/localDimSize/))
!call iamas_log_write('... total local dim size is $i ', intArgs=(/totalLocalDimSize/))

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , dimNames0, &
                                                      dimSizes0, field % hasTimeDimension, field % isDecomposed, localDimSize,totalLocalDimSize, &
                                                      precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , dimNames0, &
                                                dimSizes0, field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                                precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_1D_REAL
      new_field_list_node % real1dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_1dReal_mul

   subroutine iAMAS_streamAddField_1dReal(stream, field, arrayFieldList, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field1DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field1dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      character (len=StrKIND), dimension(0) :: dimNames0
      integer, dimension(0) :: dimSizes0
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , dimNames0, &
                                                      dimSizes0, field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                      indices, precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , dimNames0, &
                                                dimSizes0, field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                indices, precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_1D_REAL
      new_field_list_node % real1dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_1dReal

   subroutine iAMAS_streamAddField_2dReal_mul(stream, field, arrayFieldList, precision, ierr)

      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field2DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: mpi_ierr
      integer :: io_err
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field2dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"
      localDimSize = totalDimSize
      if (field % isDecomposed) then
         call MPI_Allreduce(localDimSize, totalLocalDimSize, 1, MPI_INTEGER, MPI_SUM, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
      else
         totalLocalDimSize = localDimSize
      end if

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                      field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                                      precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed,localDimSize, totalLocalDimSize, &
                                                precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize,  precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_2D_REAL
      new_field_list_node % real2dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_2dReal_mul

   subroutine iAMAS_streamAddField_2dReal(stream, field, arrayFieldList, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field2DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field2dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                      field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                      indices, precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                indices, precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_2D_REAL
      new_field_list_node % real2dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_2dReal

   subroutine iAMAS_streamAddField_3dReal_mul(stream, field, arrayFieldList, precision, ierr)

      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field3DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: mpi_ierr
      integer :: io_err
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field3dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"
      localDimSize = totalDimSize
      ! call iamas_log_write('localDimSize is $i', intArgs=(/localDimSize/))
      if (field % isDecomposed) then
         call MPI_Allreduce(localDimSize, totalLocalDimSize, 1, MPI_INTEGER, MPI_SUM, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
      else
         totalLocalDimSize = localDimSize
      end if
     ! call iamas_log_write('totalLocalDimSize is $i', intArgs=(/totalLocalDimSize/))

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                      field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                                      precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, &
                                            localDimSize,totalLocalDimSize, precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_3D_REAL
      new_field_list_node % real3dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')
!call iamas_log_write('DEBUGGING : Finished adding 3d real field '//trim(field % fieldName))

   end subroutine iAMAS_streamAddField_3dReal_mul


   subroutine iAMAS_streamAddField_3dReal(stream, field, arrayFieldList, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field3DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field3dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                      field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                      indices, precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                indices, precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_3D_REAL
      new_field_list_node % real3dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')
!call iamas_log_write('DEBUGGING : Finished adding 3d real field '//trim(field % fieldName))

   end subroutine iAMAS_streamAddField_3dReal

   subroutine iAMAS_streamAddField_4dReal_mul(stream, field, arrayFieldList, precision, ierr)

      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field4DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: mpi_ierr
      integer :: io_err
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field4dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"
      localDimSize = totalDimSize
      if (field % isDecomposed) then
         call MPI_Allreduce(localDimSize, totalLocalDimSize, 1, MPI_INTEGER, MPI_SUM, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
      else
         totalLocalDimSize = localDimSize
      end if

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                      field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                                      precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                                precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize,  precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_4D_REAL
      new_field_list_node % real4dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')
!call iamas_log_write('DEBUGGING : Finished adding 4d real field '//trim(field % fieldName))

   end subroutine iAMAS_streamAddField_4dReal_mul

   subroutine iAMAS_streamAddField_4dReal(stream, field, arrayFieldList, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field4DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field4dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                      field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                      indices, precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                indices, precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_4D_REAL
      new_field_list_node % real4dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')
!call iamas_log_write('DEBUGGING : Finished adding 4d real field '//trim(field % fieldName))

   end subroutine iAMAS_streamAddField_4dReal

   subroutine iAMAS_streamAddField_5dReal_mul(stream, field, arrayFieldList, precision, ierr)

      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field5DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: mpi_ierr
      integer :: io_err
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field5dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"
      localDimSize = totalDimSize
      if (field % isDecomposed) then
         call MPI_Allreduce(localDimSize, totalLocalDimSize, 1, MPI_INTEGER, MPI_SUM, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
      else
         totalLocalDimSize = localDimSize
      end if

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                      field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                                      precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                                precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_5D_REAL
      new_field_list_node % real5dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')
!call iamas_log_write('DEBUGGING : Finished adding 5d real field '//trim(field % fieldName))

   end subroutine iAMAS_streamAddField_5dReal_mul

   subroutine iAMAS_streamAddField_5dReal(stream, field, arrayFieldList, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field5DReal), intent(in), target :: field
      ! jwfeng
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i,j,iVarArray
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field5dReal), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      logical :: any_success
      logical, dimension(:), pointer :: isAvailable
      integer :: local_precision

      type (iamas_pool_type), pointer :: meshPool
      integer, dimension(:), pointer :: indexArray
      integer, pointer :: indexDimension
      character (len=StrKIND) :: elementName, elementNamePlural
      logical :: meshFieldDim
      ! jwfeng: For debug info
      character (len=StrKIND) :: debugMessage


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

     
      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      if (.not. field % isPersistent .or. .not. field % isActive) then
        return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = size(field % dimSizes)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
#include "add_field_indices.inc"

      
      any_success = .false.
      if (field % isVarArray) then
         allocate(isAvailable(size(field % constituentNames)))
         isAvailable(:) = .false.

         ! jwfeng: Check that if this varArray has specific vars in stream_list
         iVarArray = 0
         do i = 1, maxVarArrays
            if ( trim(field % fieldname) .eq. trim(arrayFieldList%arrayFields(i)) ) then
                iVarArray = i
            end if
         end do

         if (iVarArray .ne. 0) then
            do i=1,size(field % constituentNames)
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                    !  call iamas_log_write("Adding "//trim(field%fieldname)//"%"//trim(field%constituentNames(i))//" to array field list in io_streams")
                     call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                      field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                      indices, precision=local_precision, ierr=io_err)
                     if (io_err == iAMAS_STREAM_NOERR) then
                        isAvailable(i) = .true.
                        any_success = .true.
                     end if
                  end if
               end do
            end do
         else
            do i=1,size(field % constituentNames)
               call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), REAL_IO_TYPE , field % dimNames(2:ndims), &
                                                field % dimSizes(2:ndims), field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                                indices, precision=local_precision, ierr=io_err)
               if (io_err == iAMAS_STREAM_NOERR) then
                  isAvailable(i) = .true.
                  any_success = .true.
               end if
            end do
         end if
      else
         nullify(isAvailable)
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), REAL_IO_TYPE , field % dimNames, field % dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, precision=local_precision, ierr=io_err)
         if (io_err == iAMAS_STREAM_NOERR) then
            any_success = .true.
         end if
      end if

      deallocate(indices)
      if (.not. any_success) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            if ( iVarArray .ne. 0) then
               do j = 1, maxVarsInArray
                  if ( trim(field % constituentNames(i)) .eq. trim(arrayFieldList%varsInArray(j,iVarArray)) ) then
                     call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
                  end if
               end do
            else
               call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
            end if
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_5D_REAL
      new_field_list_node % real5dField => field
      new_field_list_node % isAvailable => isAvailable

!call iamas_log_write('... done adding field')
!call iamas_log_write('DEBUGGING : Finished adding 5d real field '//trim(field % fieldName))

   end subroutine iAMAS_streamAddField_5dReal

   subroutine iAMAS_streamAddField_0dChar_mul(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field0DChar), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field0dChar), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: dimSizes
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = 1

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(0))
      allocate(dimSizes(1))
      dimSizes(1) = ShortStrKIND
      dimNames(1) = 'StrLen'
      globalDimSize = ShortStrKIND
      totalDimSize = ShortStrKIND

      localDimSize = ShortStrKIND
      totalLocalDimSize = ShortStrKIND
      
      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), iAMAS_IO_CHAR, dimNames(1:1), &
                                             dimSizes, field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                             ierr=io_err)
         end do
      else
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), iAMAS_IO_CHAR, dimNames(1:1), dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, ierr=io_err)
      end if

      deallocate(indices)
      deallocate(dimSizes)
      if (io_err /= iAMAS_STREAM_NOERR) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_0D_CHAR
      new_field_list_node % char0dField => field

!call iamas_log_write('... done adding field'

   end subroutine iAMAS_streamAddField_0dChar_mul

   subroutine iAMAS_streamAddField_0dChar(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field0DChar), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field0dChar), pointer :: field_ptr
      character (len=StrKIND), dimension(5) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: dimSizes
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = 1

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(0))
      allocate(dimSizes(1))
      dimSizes(1) = ShortStrKIND
      dimNames(1) = 'StrLen'
      globalDimSize = ShortStrKIND
      totalDimSize = ShortStrKIND

      
      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), iAMAS_IO_CHAR, dimNames(1:1), &
                                             dimSizes, field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                             indices, ierr=io_err)
         end do
      else
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), iAMAS_IO_CHAR, dimNames(1:1), dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)
      end if

      deallocate(indices)
      deallocate(dimSizes)
      if (io_err /= iAMAS_STREAM_NOERR) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_0D_CHAR
      new_field_list_node % char0dField => field

!call iamas_log_write('... done adding field'

   end subroutine iAMAS_streamAddField_0dChar


   subroutine iAMAS_streamAddField_1dChar_mul(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field1DChar), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: localDimSize, totalLocalDimSize
      integer :: ndims
      type (field1dChar), pointer :: field_ptr
      character (len=StrKIND), dimension(2) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: dimSizes
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = 2

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(1))
      allocate(dimSizes(2))
      dimSizes(1) = ShortStrKIND
      dimNames(1) = 'StrLen'
      dimSizes(2) = field % dimSizes(1)
      dimNames(2) = field % dimNames(1)
      globalDimSize = dimSizes(2)
      totalDimSize = dimSizes(2)
      
      localDimSize = dimSizes(2)
      totalLocalDimSize = dimSizes(2)

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic_mul(stream, trim(field % constituentNames(i)), iAMAS_IO_CHAR, dimNames, &
                                             dimSizes, field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, &
                                            ierr=io_err)
         end do
      else
         call iAMAS_streamAddField_generic_mul(stream, trim(field % fieldName), iAMAS_IO_CHAR, dimNames, dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, localDimSize, totalLocalDimSize, ierr=io_err)
      end if

      deallocate(indices)
      deallocate(dimSizes)
      if (io_err /= iAMAS_STREAM_NOERR) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_1D_CHAR
      new_field_list_node % char1dField => field

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_1dChar_mul

   subroutine iAMAS_streamAddField_1dChar(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field1DChar), intent(in), target :: field
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: totalDimSize, globalDimSize
      integer :: ndims
      type (field1dChar), pointer :: field_ptr
      character (len=StrKIND), dimension(2) :: dimNames
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: dimSizes
      integer, dimension(:), pointer :: indices
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

!call iamas_log_write('... Adding field '//trim(field % fieldName)//' to stream')

      ndims = 2

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      ! 
      ! Determine whether the field is decomposed, the indices that are owned by this task's blocks,
      !    and the total number of outer-indices owned by this task
      ! 
      idim = ndims
      allocate(indices(1))
      allocate(dimSizes(2))
      dimSizes(1) = ShortStrKIND
      dimNames(1) = 'StrLen'
      dimSizes(2) = field % dimSizes(1)
      dimNames(2) = field % dimNames(1)
      globalDimSize = dimSizes(2)
      totalDimSize = dimSizes(2)
      
      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call iAMAS_streamAddField_generic(stream, trim(field % constituentNames(i)), iAMAS_IO_CHAR, dimNames, &
                                             dimSizes, field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, &
                                             indices, ierr=io_err)
         end do
      else
         call iAMAS_streamAddField_generic(stream, trim(field % fieldName), iAMAS_IO_CHAR, dimNames, dimSizes, &
                                          field % hasTimeDimension, field % isDecomposed, totalDimSize, globalDimSize, indices, ierr=io_err)
      end if

      deallocate(indices)
      deallocate(dimSizes)
      if (io_err /= iAMAS_STREAM_NOERR) then
          if (present(ierr)) ierr = iAMAS_IO_ERR
          return
      end if

      if (field % isVarArray) then
         do i=1,size(field % constituentNames)
            call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % constituentNames(i)), field % attLists(i) % attList)
         end do
      else
         call put_get_field_atts(stream % fileHandle, stream % ioDirection, trim(field % fieldname), field % attLists(1) % attList)
      end if


      !
      ! Set field pointer and type in fieldList
      !
      new_field_list_node => stream % fieldList
      do while (associated(new_field_list_node % next))
         new_field_list_node => new_field_list_node % next
      end do
      new_field_list_node % field_type = FIELD_1D_CHAR
      new_field_list_node % char1dField => field

!call iamas_log_write('... done adding field')

   end subroutine iAMAS_streamAddField_1dChar

   subroutine iAMAS_streamAddField_generic_mul(stream, fieldName, fieldType, dimNames, dimSizes, hasTimeDimension, isDecomposed, &
                                          localDimSize, totalLocalDimSize, precision, ierr)

      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: fieldName
      integer, intent(in) :: fieldType
      character (len=StrKIND), dimension(:), intent(in) :: dimNames
      integer, dimension(:), intent(in) :: dimSizes
      logical, intent(in) :: hasTimeDimension
      logical, intent(in) :: isDecomposed
      integer, intent(in) :: localDimSize
      integer, intent(in) :: totalLocalDimSize
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: mpi_ierr
      integer :: io_err
      integer :: i
      integer :: idim
      integer :: ndims
      integer :: dimTemp
      character (len=StrKIND), dimension(5) :: dimNamesLocal
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: dimSizesInq
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node
      
      integer, allocatable :: indices(:)
      integer :: elements_per_pe
      integer :: local_ntasks
      integer :: local_rank
      integer :: start1
      integer, allocatable :: dimcount(:)

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      ndims = size(dimNames)

!      call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      allocate(new_field_list_node)
      nullify(new_field_list_node % next)

      if (stream % ioDirection == iAMAS_IO_WRITE) then
!         call iamas_log_write('Call iAMAS_streamAddField_generic_mul write')
!call iamas_log_write('... defining field')

         !
         ! Define inner dimensions
         !
         elements_per_pe = 1
         do idim = 1, ndims-1
!call iamas_log_write('... defining dimension '// trim(dimNames(idim))//" $i", intArgs=(/ dimSizes(idim)/))
            write(dimNamesLocal(idim),'(a)') dimNames(idim)
            call iAMAS_io_def_dim(stream % fileHandle, trim(dimNames(idim)), dimSizes(idim), io_err)
            elements_per_pe = elements_per_pe * dimSizes(idim)
      !      call iamas_log_write('dimSizes $i ', intArgs=(/ dimSizes(idim)/))

            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
               ! call iamas_log_write('wrong at 3479')
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               return
            end if
         end do

         !
         ! Define outer-most dimension
         !
         idim = ndims
         if (idim > 0) write(dimNamesLocal(idim),'(a)') dimNames(idim)

         if (isDecomposed) then
             ! Used when allocating temp array
            new_field_list_node % totalDimSize = localDimSize
         else
            new_field_list_node % totalDimSize = totalLocalDimSize
         end if

         new_field_list_node % isDecomposed = isDecomposed

         if (ndims > 0) then
            call iAMAS_io_def_dim(stream % fileHandle, trim(dimNames(idim)), totalLocalDimSize, io_err)
            
            elements_per_pe = elements_per_pe * localDimSize

            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
               ! call iamas_log_write('wrong at 3508')
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               return
            end if
         end if

         !
         ! Define time dimension if necessary
         !
         if (hasTimeDimension) then
!call iamas_log_write('... defining Time dimension ')
            call iAMAS_io_def_dim(stream % fileHandle, 'Time', iAMAS_IO_UNLIMITED_DIM, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
               ! call iamas_log_write('wrong at 3525')
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               return
            end if
            ndims = ndims + 1
            write(dimNamesLocal(ndims),'(a)') 'Time'
         end if

         !
         ! Define variable to low-level interface
         !
!call iamas_log_write('... defining var to low-level interface with ndims $i', intArgs=(/ndims/))

         call iAMAS_io_def_var(stream % fileHandle, trim(fieldName), fieldType, dimNamesLocal(1:ndims), precision=precision, ierr=io_err)
         call iAMAS_io_err_mesg(io_err, .false.)
         if (io_err /= iAMAS_IO_NOERR) then
            ! call iamas_log_write('wrong at 3542')
            if (present(ierr)) ierr = iAMAS_IO_ERR
            deallocate(new_field_list_node)
            return
         end if
!         call iamas_log_write('Call iAMAS_streamAddField_generic_mul finish def')
        
      else if (stream % ioDirection == iAMAS_IO_READ) then
         ! call iamas_log_write('... inquiring about')
!         call iamas_log_write('... inquiring about '//trim(fieldName))

         call iAMAS_io_inq_var(stream % fileHandle, trim(fieldName), dimnames=dimNamesInq, dimsizes=dimSizesInq, ierr=io_err)
         ! If the field does not exist in the input file, we should handle this situation gracefully at higher levels
         !   without printing disconcerting error messages
         !call iAMAS_io_err_mesg(io_err, .false.)
         if (io_err /= iAMAS_IO_NOERR) then
!            call iamas_log_write('... inquiring '//trim(fieldName)//' error')
            if (present(ierr)) ierr = iAMAS_IO_ERR
            deallocate(new_field_list_node)
            return
         end if

! Here, we should probably do a check to make sure the file agrees with what iAMAS expects for the field
         do i=1,ndims
!            call iamas_log_write('Comparing '//trim(dimNames(i))//' '//trim(dimNamesInq(i)))
            if (trim(dimNames(i)) /= trim(dimNamesInq(i))) then
               ! call iamas_log_write('Mismatched dimension name in field')
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               deallocate(dimNamesInq)
               deallocate(dimSizesInq)
               return
            end if
            if (i < ndims) then
               dimTemp = dimSizes(i)
            else
               if ( isDecomposed ) then
                  dimTemp = totalLocalDimSize
               else
                  dimTemp = dimSizes(i)
               end if
            end if
!            call iamas_log_write('Comparing $i $i', intArgs=(/dimTemp, dimSizesInq(i)/))
            if (dimTemp /= dimSizesInq(i)) then
!               call iamas_log_write('Mismatched dimension size in field')
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               deallocate(dimNamesInq)
               deallocate(dimSizesInq)
               return
            end if
         end do 

         ! Set outer dimension sizes depending on whether the field is decomposed
         if (isDecomposed) then
!            call iamas_log_write(' is Decomposed' //trim(fieldName))
            new_field_list_node % totalDimSize = localDimSize
         else
!            call iamas_log_write(' is not Decomposed' //trim(fieldName))
            new_field_list_node % totalDimSize = totalLocalDimSize
         end if

         new_field_list_node % isDecomposed = isDecomposed

         deallocate(dimNamesInq)
         deallocate(dimSizesInq)

      end if
!      
!      if (isDecomposed) then
!         call iamas_log_write(' is Decomposed' //trim(fieldName))
!      else
!         call iamas_log_write(' is not Decomposed' //trim(fieldName))
!      end if

      ! Generate indices
      call MPI_Comm_size(stream % fileHandle % ioContext % dminfo % local_comm, local_ntasks, mpi_ierr)
      call MPI_Comm_rank(stream % fileHandle % ioContext % dminfo % local_comm, local_rank, mpi_ierr)
      
      allocate(indices(1:localDimSize))
      allocate(dimcount(1:local_ntasks))
    
      call MPI_Allgather(localDimSize, 1, MPI_INT, dimcount, 1, MPI_INT, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
!      call iamas_log_write('iAMAS_streamAddField_generic_mul finish MPI_Allgather')
      
!      do i = 1, local_ntasks
!         call iamas_log_write('dimcount $i', intArgs=(/dimcount(i)/))
!      end do
      
      start1 = 1
      do i = 1, local_rank
         start1 = start1 + dimcount(i)
      end do
      
      indices(1:localDimSize) = (/(i, i=start1, start1 + localDimSize - 1, 1)/)
      
      ! do i = 1, size(indices)
      !    call iamas_log_write('indices $i', intArgs=(/indices(i)/) )
      ! end do
      !
      ! Set variable indices
      !
      if (ndims > 0 .and. isDecomposed) then
      !   call iamas_log_write('call iAMAS_io_set_var_indices_mul')
         call iAMAS_io_set_var_indices_mul(stream % fileHandle, trim(fieldName), indices, io_err)
      !   call iamas_log_write('iAMAS_streamAddField_generic_mul finish set var indices mul')
         call iAMAS_io_err_mesg(io_err, .false.)
         if (io_err /= iAMAS_IO_NOERR) then
            if (present(ierr)) ierr = iAMAS_IO_ERR
            deallocate(new_field_list_node)
            return
         end if
      end if


      !
      ! Add field pointer to the list
      !
      if (.not. associated(stream % fieldList)) then
!call iamas_log_write('Adding field to the head of the list')
         stream % fieldList => new_field_list_node
      else
!call iamas_log_write('Adding field to the tail of the list')
         field_list_cursor => stream % fieldList
         do while (associated(field_list_cursor % next))
            field_list_cursor => field_list_cursor % next
         end do
         field_list_cursor % next => new_field_list_node
      end if
      
      deallocate(indices)
      deallocate(dimcount)
!      call iamas_log_write('iAMAS_streamAddField_generic_mul finish')

   end subroutine iAMAS_streamAddField_generic_mul

   subroutine iAMAS_streamAddField_generic(stream, fieldName, fieldType, dimNames, dimSizes, hasTimeDimension, isDecomposed, &
                                          totalDimSize, globalDimSize, indices, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: fieldName
      integer, intent(in) :: fieldType
      character (len=StrKIND), dimension(:), intent(in) :: dimNames
      integer, dimension(:), intent(in) :: dimSizes
      logical, intent(in) :: hasTimeDimension
      logical, intent(in) :: isDecomposed
      integer, intent(in) :: totalDimSize
      integer, intent(in) :: globalDimSize
      integer, dimension(:), intent(in) :: indices
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i
      integer :: idim
      integer :: ndims
      integer :: dimTemp
      character (len=StrKIND), dimension(5) :: dimNamesLocal
      character (len=StrKIND), dimension(:), pointer :: dimNamesInq
      integer, dimension(:), pointer :: dimSizesInq
      type (field_list_type), pointer :: field_list_cursor
      type (field_list_type), pointer :: new_field_list_node

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      ndims = size(dimNames)

!call iamas_log_write('... field has $i dimensions', intArgs=(/ndims/))

      allocate(new_field_list_node)
      nullify(new_field_list_node % next)

      if (stream % ioDirection == iAMAS_IO_WRITE) then
!call iamas_log_write('... defining field')

         !
         ! Define inner dimensions
         !
         do idim = 1, ndims-1
!call iamas_log_write('... defining dimension '// trim(dimNames(idim))//" $i", intArgs=(/ dimSizes(idim)/))
            write(dimNamesLocal(idim),'(a)') dimNames(idim)
            call iAMAS_io_def_dim(stream % fileHandle, trim(dimNames(idim)), dimSizes(idim), io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               return
            end if
         end do

         !
         ! Define outer-most dimension
         !
         idim = ndims
         if (idim > 0) write(dimNamesLocal(idim),'(a)') dimNames(idim)

         if (isDecomposed) then
            new_field_list_node % totalDimSize = totalDimSize
         else
            new_field_list_node % totalDimSize = globalDimSize
         end if

         new_field_list_node % isDecomposed = isDecomposed
  
         if (ndims > 0) then
!call iamas_log_write('... defining dimension '// trim(dimNames(idim))//" $i", intArgs=(/ globalDimSize/))
            call iAMAS_io_def_dim(stream % fileHandle, trim(dimNames(idim)), globalDimSize, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               return
            end if
         end if

         !
         ! Define time dimension if necessary
         !
         if (hasTimeDimension) then
!call iamas_log_write('... defining Time dimension ')
            call iAMAS_io_def_dim(stream % fileHandle, 'Time', iAMAS_IO_UNLIMITED_DIM, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               return
            end if
            ndims = ndims + 1
            write(dimNamesLocal(ndims),'(a)') 'Time'
         end if

         !
         ! Define variable to low-level interface
         !
!call iamas_log_write('... defining var to low-level interface with ndims $i', intArgs=(/ndims/))

         call iAMAS_io_def_var(stream % fileHandle, trim(fieldName), fieldType, dimNamesLocal(1:ndims), precision=precision, ierr=io_err)
         call iAMAS_io_err_mesg(io_err, .false.)
         if (io_err /= iAMAS_IO_NOERR) then
            if (present(ierr)) ierr = iAMAS_IO_ERR
            deallocate(new_field_list_node)
            return
         end if
        
      else if (stream % ioDirection == iAMAS_IO_READ) then
!call iamas_log_write('... inquiring about')

         call iAMAS_io_inq_var(stream % fileHandle, trim(fieldName), dimnames=dimNamesInq, dimsizes=dimSizesInq, ierr=io_err)
         ! If the field does not exist in the input file, we should handle this situation gracefully at higher levels
         !   without printing disconcerting error messages
         !call iAMAS_io_err_mesg(io_err, .false.)
         if (io_err /= iAMAS_IO_NOERR) then
            if (present(ierr)) ierr = iAMAS_IO_ERR
            deallocate(new_field_list_node)
            return
         end if

! Here, we should probably do a check to make sure the file agrees with what iAMAS expects for the field
         do i=1,ndims
!call iamas_log_write('Comparing '//trim(dimNames(i))//' '//trim(dimNamesInq(i)))
            if (trim(dimNames(i)) /= trim(dimNamesInq(i))) then
!call iamas_log_write('Mismatched dimension name in field')
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               deallocate(dimNamesInq)
               deallocate(dimSizesInq)
               return
            end if
            if (i < ndims) then
               dimTemp = dimSizes(i)
            else
               if ( isDecomposed ) then
                  dimTemp = globalDimSize
               else
                  dimTemp = dimSizes(i)
               end if
            end if
!call iamas_log_write('Comparing $i $i', intArgs=(/dimTemp, dimSizesInq(i)/))
            if (dimTemp /= dimSizesInq(i)) then
!call iamas_log_write('Mismatched dimension size in field')
               if (present(ierr)) ierr = iAMAS_IO_ERR
               deallocate(new_field_list_node)
               deallocate(dimNamesInq)
               deallocate(dimSizesInq)
               return
            end if
         end do 

         ! Set outer dimension sizes depending on whether the field is decomposed
         if (isDecomposed) then
            new_field_list_node % totalDimSize = totalDimSize
         else
            new_field_list_node % totalDimSize = globalDimSize
         end if

         new_field_list_node % isDecomposed = isDecomposed

         deallocate(dimNamesInq)
         deallocate(dimSizesInq)

      end if

      
      !
      ! Set variable indices
      !
      if (ndims > 0 .and. isDecomposed) then
         call iAMAS_io_set_var_indices(stream % fileHandle, trim(fieldName), indices, io_err)
         call iAMAS_io_err_mesg(io_err, .false.)
         if (io_err /= iAMAS_IO_NOERR) then
            if (present(ierr)) ierr = iAMAS_IO_ERR
            deallocate(new_field_list_node)
            return
         end if
      end if


      !
      ! Add field pointer to the list
      !
      if (.not. associated(stream % fieldList)) then
!call iamas_log_write('Adding field to the head of the list')
         stream % fieldList => new_field_list_node
      else
!call iamas_log_write('Adding field to the tail of the list')
         field_list_cursor => stream % fieldList
         do while (associated(field_list_cursor % next))
            field_list_cursor => field_list_cursor % next
         end do
         field_list_cursor % next => new_field_list_node
      end if

   end subroutine iAMAS_streamAddField_generic


   subroutine iAMAS_streamUpdateField_0dInteger(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field0DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 0d integer field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_0D_INT) then
            if (field_cursor % int0dField % fieldname == field % fieldname) then
               STREAM_DEBUG_WRITE( '... found 0d integer named '//trim(field_cursor % int0dField % fieldname) )
               field_cursor % int0dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 0d integer field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_0dInteger


   subroutine iAMAS_streamUpdateField_1dInteger(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field1DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 1d integer field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_1D_INT) then
            if (field_cursor % int1dField % fieldname == field % fieldname .and. &
                field_cursor % int1dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % int1dField % dimNames(1) == field % dimNames(1)) then
               STREAM_DEBUG_WRITE( '... found 1d integer named '//trim(field_cursor % int1dField % fieldname) )
               field_cursor % int1dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 1d integer field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_1dInteger


   subroutine iAMAS_streamUpdateField_2dInteger(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field2DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 2d integer field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_2D_INT) then
            if (field_cursor % int2dField % fieldname == field % fieldname .and. &
                field_cursor % int2dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % int2dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % int2dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % int2dField % dimNames(2) == field % dimNames(2)) then
               STREAM_DEBUG_WRITE( '... found 2d integer named '//trim(field_cursor % int2dField % fieldname) )
               field_cursor % int2dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 2d integer field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_2dInteger


   subroutine iAMAS_streamUpdateField_3dInteger(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field3DInteger), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 3d integer field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_3D_INT) then
            if (field_cursor % int3dField % fieldname == field % fieldname .and. &
                field_cursor % int3dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % int3dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % int3dField % dimSizes(3) == field % dimSizes(3) .and. &
                field_cursor % int3dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % int3dField % dimNames(2) == field % dimNames(2) .and. &
                field_cursor % int3dField % dimNames(3) == field % dimNames(3)) then
               STREAM_DEBUG_WRITE( '... found 3d integer named '//trim(field_cursor % int3dField % fieldname) )
               field_cursor % int3dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 3d integer field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_3dInteger


   subroutine iAMAS_streamUpdateField_0dReal(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field0DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 0d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_0D_REAL) then
            if (field_cursor % real0dField % fieldname == field % fieldname) then
               STREAM_DEBUG_WRITE( '... found 0d real named '//trim(field_cursor % real0dField % fieldname) )
               field_cursor % real0dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 0d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_0dReal


   subroutine iAMAS_streamUpdateField_1dReal(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field1DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 1d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_1D_REAL) then
            if (field_cursor % real1dField % fieldname == field % fieldname .and. &
                field_cursor % real1dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real1dField % dimNames(1) == field % dimNames(1)) then
               STREAM_DEBUG_WRITE( '... found 1d real named '//trim(field_cursor % real1dField % fieldname) )
               field_cursor % real1dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 1d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_1dReal


   subroutine iAMAS_streamUpdateField_2dReal(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field2DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 2d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_2D_REAL) then
            if (field_cursor % real2dField % fieldname == field % fieldname .and. &
                field_cursor % real2dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real2dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % real2dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % real2dField % dimNames(2) == field % dimNames(2)) then
               STREAM_DEBUG_WRITE( '... found 2d real named '//trim(field_cursor % real2dField % fieldname) )
               field_cursor % real2dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 2d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_2dReal


   subroutine iAMAS_streamUpdateField_3dReal(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field3DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 3d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_3D_REAL) then
            if (field_cursor % real3dField % fieldname == field % fieldname .and. &
                field_cursor % real3dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real3dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % real3dField % dimSizes(3) == field % dimSizes(3) .and. &
                field_cursor % real3dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % real3dField % dimNames(2) == field % dimNames(2) .and. &
                field_cursor % real3dField % dimNames(3) == field % dimNames(3)) then
               STREAM_DEBUG_WRITE( '... found 3d real named '//trim(field_cursor % real3dField % fieldname) )
               field_cursor % real3dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 3d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_3dReal


   subroutine iAMAS_streamUpdateField_4dReal(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field4DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 4d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_4D_REAL) then
            if (field_cursor % real4dField % fieldname == field % fieldname .and. &
                field_cursor % real4dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real4dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % real4dField % dimSizes(3) == field % dimSizes(3) .and. &
                field_cursor % real4dField % dimSizes(4) == field % dimSizes(4) .and. &
                field_cursor % real4dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % real4dField % dimNames(2) == field % dimNames(2) .and. &
                field_cursor % real4dField % dimNames(3) == field % dimNames(3) .and. &
                field_cursor % real4dField % dimNames(4) == field % dimNames(4)) then
               STREAM_DEBUG_WRITE( '... found 4d real named '//trim(field_cursor % real4dField % fieldname) )
               field_cursor % real4dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 4d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_4dReal


   subroutine iAMAS_streamUpdateField_5dReal(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field5DReal), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 5d real field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_5D_REAL) then
            if (field_cursor % real5dField % fieldname == field % fieldname .and. &
                field_cursor % real5dField % dimSizes(1) == field % dimSizes(1) .and. &
                field_cursor % real5dField % dimSizes(2) == field % dimSizes(2) .and. &
                field_cursor % real5dField % dimSizes(3) == field % dimSizes(3) .and. &
                field_cursor % real5dField % dimSizes(4) == field % dimSizes(4) .and. &
                field_cursor % real5dField % dimSizes(5) == field % dimSizes(5) .and. &
                field_cursor % real5dField % dimNames(1) == field % dimNames(1) .and. &
                field_cursor % real5dField % dimNames(2) == field % dimNames(2) .and. &
                field_cursor % real5dField % dimNames(3) == field % dimNames(3) .and. &
                field_cursor % real5dField % dimNames(4) == field % dimNames(4) .and. &
                field_cursor % real5dField % dimNames(5) == field % dimNames(5)) then
               STREAM_DEBUG_WRITE( '... found 5d real named '//trim(field_cursor % real5dField % fieldname) )
               field_cursor % real5dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 5d real field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_5dReal


   subroutine iAMAS_streamUpdateField_0dChar(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field0DChar), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 0d char field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_0D_CHAR) then
            if (field_cursor % char0dField % fieldname == field % fieldname) then
               STREAM_DEBUG_WRITE( '... found 0d char named '//trim(field_cursor % char0dField % fieldname) )
               field_cursor % char0dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 0d char field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_0dChar


   subroutine iAMAS_streamUpdateField_1dChar(stream, field, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      type (field1DChar), intent(in), target :: field
      integer, intent(out), optional :: ierr

      type (field_list_type), pointer :: field_cursor


      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      STREAM_DEBUG_WRITE( '... Updating 1d char field '//trim(field % fieldName)//' in stream' )

      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_1D_CHAR) then
            if (field_cursor % char1dField % fieldname == field % fieldname) then
               STREAM_DEBUG_WRITE( '... found 1d char named '//trim(field_cursor % char0dField % fieldname) )
               field_cursor % char1dField => field
               exit
            end if
         end if
         field_cursor => field_cursor % next
      end do
      if (.not. associated(field_cursor)) then
         STREAM_DEBUG_WRITE( '... 1d char field '//trim(field % fieldname)//' not found in stream' )
         if (present(ierr)) ierr = iAMAS_STREAM_FIELD_NOT_FOUND
      end if

      STREAM_DEBUG_WRITE( '... done updating field' )

   end subroutine iAMAS_streamUpdateField_1dChar


   subroutine iAMAS_readStream_mul(stream, frame, ierr)

      use pio
      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      integer, intent(in) :: frame
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: mpi_ierr
      integer :: i, j
      integer :: ncons
      integer, pointer :: ownedSize
      type (field0dInteger), pointer :: field_0dint_ptr
      type (field1dInteger), pointer :: field_1dint_ptr
      type (field2dInteger), pointer :: field_2dint_ptr
      type (field3dInteger), pointer :: field_3dint_ptr
      type (field0dReal), pointer :: field_0dreal_ptr
      type (field1dReal), pointer :: field_1dreal_ptr
      type (field2dReal), pointer :: field_2dreal_ptr
      type (field3dReal), pointer :: field_3dreal_ptr
      type (field4dReal), pointer :: field_4dreal_ptr
      type (field5dReal), pointer :: field_5dreal_ptr
      type (field0dChar), pointer :: field_0dchar_ptr
      type (field1dChar), pointer :: field_1dchar_ptr
      type (field_list_type), pointer :: field_cursor
      character (len=StrKIND), dimension(:), pointer :: char1d_temp
      integer                            :: int0d_temp
      integer, dimension(:),     pointer :: int1d_temp
      integer, dimension(:,:),   pointer :: int2d_temp
      integer, dimension(:,:,:), pointer :: int3d_temp
      real (kind=RKIND)                            :: real0d_temp
      real (kind=RKIND), dimension(:),     pointer :: real1d_temp
      real (kind=RKIND), dimension(:,:),   pointer :: real2d_temp
      real (kind=RKIND), dimension(:,:,:), pointer :: real3d_temp
      real (kind=RKIND), dimension(:,:,:,:), pointer :: real4d_temp
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: real5d_temp


      real (kind=R8KIND) :: timer_start_time, timer_stop_time!jwfeng
      real (kind=R8KIND) :: sub_timer_start_time, sub_timer_stop_time!jwfeng
      real (kind=R8KIND) :: sub_timer_sum, reading_time

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR
      reading_time = 0.0

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if


      !
      ! Set time frame to real
      !
      call iAMAS_io_set_frame(stream % fileHandle, frame, io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if


      !
      ! Loop over fields in the stream
      !
      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_0D_INT) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % int0dField % fieldName))
call iamas_log_write('Reading in field '//trim(field_cursor % int0dField % fieldName))
! call iamas_log_write('FIELD_0D_INT Reading in field '//trim(field_cursor % int0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('MGD calling iAMAS_io_get_var now...')
            call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % int0dField % fieldName, int0d_temp, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
                if (present(ierr)) ierr = iAMAS_IO_ERR
                return
            end if

!call iamas_log_write('Distributing and Copying field to other blocks')

            ! call iamas_dmpar_bcast_int( iamas_io_handle_dminfo(stream % fileHandle), int0d_temp)
            call MPI_Bcast(int0d_temp, 1, MPI_INTEGER, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
            field_0dint_ptr => field_cursor % int0dField
            do while (associated(field_0dint_ptr))
               field_0dint_ptr % scalar = int0d_temp
               field_0dint_ptr => field_0dint_ptr % next
            end do

         else if (field_cursor % field_type == FIELD_1D_INT) then
!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % int1dField % fieldName))
 call iamas_log_write('Reading in field '//trim(field_cursor % int1dField % fieldName))

            if (field_cursor % int1dField % isVarArray) then
               ncons = size(field_cursor % int1dField % constituentNames)
            else
               ncons = 1
               allocate(int1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % int1dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % int1dField % constituentNames(j), int0d_temp, io_err)
               else
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % int1dField % fieldName, int1d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (.not. field_cursor % int1dField % isVarArray) then
                      deallocate(int1d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_1dint_ptr => field_cursor % int1dField
                  i = 1
                  do while (associated(field_1dint_ptr))
                     if (trim(field_1dint_ptr % dimNames(1)) == 'nCells') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, field_1dint_ptr % dimNames(1), ownedSize)
                     end if

                     if (field_cursor % int1dField % isVarArray) then
                        field_1dint_ptr % array(j) = int0d_temp
                     else
                        field_1dint_ptr % array(1:ownedSize) = int1d_temp(i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_1dint_ptr => field_1dint_ptr % next
                  end do

               else
   
                  if (field_cursor % int1dField % isVarArray) then
                     ! call iamas_dmpar_bcast_int( iamas_io_handle_dminfo(stream % fileHandle), int0d_temp)
                     call MPI_Bcast(int0d_temp, 1, MPI_INTEGER, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_1dint_ptr => field_cursor % int1dField
                     do while (associated(field_1dint_ptr))
                        field_1dint_ptr % array(j) = int0d_temp
                        field_1dint_ptr => field_1dint_ptr % next
                     end do
                  else
                     ! call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int1d_temp), int1d_temp(:))
                     call MPI_Bcast(int1d_temp(:), size(int1d_temp), MPI_INTEGER, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_1dint_ptr => field_cursor % int1dField
                     do while (associated(field_1dint_ptr))
                        field_1dint_ptr % array(:) = int1d_temp(:)
                        field_1dint_ptr => field_1dint_ptr % next
                     end do
                  end if
               end if
            end do

            if (.not. field_cursor % int1dField % isVarArray) then
               deallocate(int1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_INT) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % int2dField % fieldName))
call iamas_log_write('FIELD_2D_INT Reading in field '//trim(field_cursor % int2dField % fieldName))
            if (field_cursor % int2dField % isVarArray) then
               ncons = size(field_cursor % int2dField % constituentNames)
               allocate(int1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int2d_temp(field_cursor % int2dField % dimSizes(1), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % int2dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % int2dField % constituentNames(j), int1d_temp, io_err)
               else
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % int2dField % fieldName, int2d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % int2dField % isVarArray) then
                      deallocate(int1d_temp)
                   else
                      deallocate(int2d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_2dint_ptr => field_cursor % int2dField
                  i = 1
                  do while (associated(field_2dint_ptr))
                     if (trim(field_2dint_ptr % dimNames(2)) == 'nCells') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, field_2dint_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % int2dField % isVarArray) then
                        field_2dint_ptr % array(j,1:ownedSize) = int1d_temp(i:i+ownedSize-1)
                     else
                        field_2dint_ptr % array(:,1:ownedSize) = int2d_temp(:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_2dint_ptr => field_2dint_ptr % next
                  end do

               else
   
                  if (field_cursor % int2dField % isVarArray) then
                     ! call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int1d_temp), int1d_temp(:))
                     call MPI_Bcast(int1d_temp(:), size(int1d_temp), MPI_INTEGER, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_2dint_ptr => field_cursor % int2dField
                     do while (associated(field_2dint_ptr))
                        field_2dint_ptr % array(j,:) = int1d_temp(:)
                        field_2dint_ptr => field_2dint_ptr % next
                     end do
                  else
                     ! call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int2d_temp), int2d_temp(:,1))
                     call MPI_Bcast(int2d_temp(:,1), size(int2d_temp), MPI_INTEGER, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_2dint_ptr => field_cursor % int2dField
                     do while (associated(field_2dint_ptr))
                        field_2dint_ptr % array(:,:) = int2d_temp(:,:)
                        field_2dint_ptr => field_2dint_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % int2dField % isVarArray) then
               deallocate(int1d_temp)
            else
               deallocate(int2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_INT) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % int3dField % fieldName))
 call iamas_log_write('FIELD_3D_INT Reading in field '//trim(field_cursor % int3dField % fieldName))
            if (field_cursor % int3dField % isVarArray) then
               ncons = size(field_cursor % int3dField % constituentNames)
               allocate(int2d_temp(field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int3d_temp(field_cursor % int3dField % dimSizes(1), &
                                   field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % int3dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % int3dField % constituentNames(j), int2d_temp, io_err)
               else
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % int3dField % fieldName, int3d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % int3dField % isVarArray) then
                      deallocate(int2d_temp)
                   else
                      deallocate(int3d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_3dint_ptr => field_cursor % int3dField
                  i = 1
                  do while (associated(field_3dint_ptr))
                     if (trim(field_3dint_ptr % dimNames(3)) == 'nCells') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, field_3dint_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % int3dField % isVarArray) then
                        field_3dint_ptr % array(j,:,1:ownedSize) = int2d_temp(:,i:i+ownedSize-1)
                     else
                        field_3dint_ptr % array(:,:,1:ownedSize) = int3d_temp(:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_3dint_ptr => field_3dint_ptr % next
                  end do

               else
   
                  if (field_cursor % int3dField % isVarArray) then
                     ! call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int2d_temp), int2d_temp(:,1))
                     call MPI_Bcast(int2d_temp(:,1), size(int2d_temp), MPI_INTEGER, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_3dint_ptr => field_cursor % int3dField
                     do while (associated(field_3dint_ptr))
                        field_3dint_ptr % array(j,:,:) = int2d_temp(:,:)
                        field_3dint_ptr => field_3dint_ptr % next
                     end do
                  else
                     ! call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int3d_temp), int3d_temp(:,1,1))
                     call MPI_Bcast(int3d_temp(:,1,1), size(int3d_temp), MPI_INTEGER, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_3dint_ptr => field_cursor % int3dField
                     do while (associated(field_3dint_ptr))
                        field_3dint_ptr % array(:,:,:) = int3d_temp(:,:,:)
                        field_3dint_ptr => field_3dint_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % int3dField % isVarArray) then
               deallocate(int2d_temp)
            else
               deallocate(int3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_0D_REAL) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real0dField % fieldName))
call iamas_log_write('Reading in field '//trim(field_cursor % real0dField % fieldName))
! call iamas_log_write('FIELD_0D_REAL Reading in field '//trim(field_cursor % real0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/) )

!call iamas_log_write('MGD calling iAMAS_io_get_var now...')
            call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real0dField % fieldName, real0d_temp, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
                if (present(ierr)) ierr = iAMAS_IO_ERR
                return
            end if

!call iamas_log_write('Distributing and Copying field to other blocks')

            ! call iamas_dmpar_bcast_real( iamas_io_handle_dminfo(stream % fileHandle), real0d_temp)
            call MPI_Bcast(real0d_temp, 1, MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
            field_0dreal_ptr => field_cursor % real0dField
            do while (associated(field_0dreal_ptr))
               field_0dreal_ptr % scalar = real0d_temp
               field_0dreal_ptr => field_0dreal_ptr % next
            end do

         else if (field_cursor % field_type == FIELD_1D_REAL) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real1dField % fieldName))
 call iamas_log_write('FIELD_1D_REAL Reading in field '//trim(field_cursor % real1dField % fieldName))
            if (field_cursor % real1dField % isVarArray) then
               ! call iamas_log_write('isVararray')
               ncons = size(field_cursor % real1dField % constituentNames)
            else
               ncons = 1
               ! call iamas_log_write('is not Vararray')
               allocate(real1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real1dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real1dField % constituentNames(j), real0d_temp, io_err)
               else
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real1dField % fieldName, real1d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (.not. field_cursor % real1dField % isVarArray) then
                      deallocate(real1d_temp)
                   end if
                   return
               end if
               ! call iamas_log_write('After get var mul ')
   
               if (field_cursor % isDecomposed) then
                  ! call iamas_log_write(' is decomposed ')
                  ! Distribute field to multiple blocks
                  field_1dreal_ptr => field_cursor % real1dField
                  i = 1

                  do while (associated(field_1dreal_ptr))
                     if (trim(field_1dreal_ptr % dimNames(1)) == 'nCells') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, field_1dreal_ptr % dimNames(1), ownedSize)
                     end if
                     ! call iamas_log_write('ownedSize   $i', intArgs=(/ownedSize/))

                     if (field_cursor % real1dField % isVarArray) then
                        field_1dreal_ptr % array(j) = real0d_temp
                     else
                        field_1dreal_ptr % array(1:ownedSize) = real1d_temp(i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_1dreal_ptr => field_1dreal_ptr % next
                  end do

               else
                  ! call iamas_log_write('start bcast')
   
                  if (field_cursor % real1dField % isVarArray) then
                     ! call iamas_dmpar_bcast_real( iamas_io_handle_dminfo(stream % fileHandle), real0d_temp)
                     call MPI_Bcast(real0d_temp, 1, MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_1dreal_ptr => field_cursor % real1dField
                     do while (associated(field_1dreal_ptr))
                        field_1dreal_ptr % array(j) = real0d_temp
                        field_1dreal_ptr => field_1dreal_ptr % next
                     end do
                  else
                     ! call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real1d_temp), real1d_temp(:))
                     call MPI_Bcast(real1d_temp(:), size(real1d_temp), MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_1dreal_ptr => field_cursor % real1dField
                     do while (associated(field_1dreal_ptr))
                        field_1dreal_ptr % array(:) = real1d_temp(:)
                        field_1dreal_ptr => field_1dreal_ptr % next
                     end do
                  end if
                  ! call iamas_log_write('finish bcast')
               end if
            end do

            if (.not. field_cursor % real1dField % isVarArray) then
               deallocate(real1d_temp)
            end if
            ! call iamas_log_write('FINISH  FIELD_1D_REAL Reading in field '//trim(field_cursor % real1dField % fieldName))

         else if (field_cursor % field_type == FIELD_2D_REAL) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real2dField % fieldName))
call iamas_dmpar_get_time(timer_start_time)!jwfeng
 call iamas_log_write('FIELD_2D_REAL Reading in field '//trim(field_cursor % real2dField % fieldName))
            if (field_cursor % real2dField % isVarArray) then
               ncons = size(field_cursor % real2dField % constituentNames)
               allocate(real1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real2d_temp(field_cursor % real2dField % dimSizes(1), &
                                    field_cursor % totalDimSize))
            end if

            sub_timer_sum = 0
            do j=1,ncons
               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % real2dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iamas_log_write('--- reading in constitutent '//trim(field_cursor % real2dField % constituentNames(j)))!jwfeng
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real2dField % constituentNames(j), real1d_temp, io_err)
               else
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real2dField % fieldName, real2d_temp, io_err)
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               sub_timer_sum = sub_timer_sum + (sub_timer_stop_time - sub_timer_start_time)
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % real2dField % isVarArray) then
                      deallocate(real1d_temp)
                   else
                      deallocate(real2d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_2dreal_ptr => field_cursor % real2dField
                  i = 1
                  do while (associated(field_2dreal_ptr))
                     if (trim(field_2dreal_ptr % dimNames(2)) == 'nCells') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, field_2dreal_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % real2dField % isVarArray) then
                        field_2dreal_ptr % array(j,1:ownedSize) = real1d_temp(i:i+ownedSize-1)
                     else
                        field_2dreal_ptr % array(:,1:ownedSize) = real2d_temp(:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_2dreal_ptr => field_2dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real2dField % isVarArray) then
                     ! call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real1d_temp), real1d_temp(:))
                     call MPI_Bcast(real1d_temp(:), size(real1d_temp), MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_2dreal_ptr => field_cursor % real2dField
                     do while (associated(field_2dreal_ptr))
                        field_2dreal_ptr % array(j,:) = real1d_temp(:)
                        field_2dreal_ptr => field_2dreal_ptr % next
                     end do
                  else
                     ! call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real2d_temp), real2d_temp(:,1))
                     call MPI_Bcast(real2d_temp(:,1), size(real2d_temp), MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_2dreal_ptr => field_cursor % real2dField
                     do while (associated(field_2dreal_ptr))
                        field_2dreal_ptr % array(:,:) = real2d_temp(:,:)
                        field_2dreal_ptr => field_2dreal_ptr % next
                     end do
                  end if
               end if
            end do
            reading_time = reading_time + sub_timer_sum

            if (field_cursor % real2dField % isVarArray) then
               deallocate(real1d_temp)
            else
               deallocate(real2d_temp)
            end if
call iamas_dmpar_get_time(timer_stop_time)!jwfeng
call iamas_log_write('Timing for reading field '//trim(field_cursor % real2dField % fieldName)// ': $r s', &
     realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng
call iamas_log_write('Timing for reading field kernel '//trim(field_cursor % real2dField % fieldName)// ': $r s', &
     realArgs=(/real(sub_timer_sum, kind = RKIND)/))!jwfeng

         else if (field_cursor % field_type == FIELD_3D_REAL) then

call iamas_dmpar_get_time(timer_start_time)!jwfeng
!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real3dField % fieldName))
!call iamas_log_write('DEBUGGING : reading a 3d real array')
! call iamas_log_write('FIELD_3D_REAL Reading in field '//trim(field_cursor % real3dField % fieldName))
            if (field_cursor % real3dField % isVarArray) then
!call iamas_log_write('DEBUGGING : reading a 3d real var-array')
               ncons = size(field_cursor % real3dField % constituentNames)
               allocate(real2d_temp(field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real3d_temp(field_cursor % real3dField % dimSizes(1), &
                                    field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            end if

            sub_timer_sum = 0
            do j=1,ncons
               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % real3dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real3dField % constituentNames(j), real2d_temp, io_err)
               else
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real3dField % fieldName, real3d_temp, io_err)
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               sub_timer_sum = sub_timer_sum + (sub_timer_stop_time - sub_timer_start_time)
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % real3dField % isVarArray) then
                      deallocate(real2d_temp)
                   else
                      deallocate(real3d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_3dreal_ptr => field_cursor % real3dField
                  i = 1
                  do while (associated(field_3dreal_ptr))
                     if (trim(field_3dreal_ptr % dimNames(3)) == 'nCells') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, field_3dreal_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % real3dField % isVarArray) then
!call iamas_log_write('DEBUGGING : copying the temporary array')
                        field_3dreal_ptr % array(j,:,1:ownedSize) = real2d_temp(:,i:i+ownedSize-1)
                     else
                        field_3dreal_ptr % array(:,:,1:ownedSize) = real3d_temp(:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_3dreal_ptr => field_3dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real3dField % isVarArray) then
                     ! call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real2d_temp), real2d_temp(:,1))
                     call MPI_Bcast(real2d_temp(:,1), size(real2d_temp), MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_3dreal_ptr => field_cursor % real3dField
                     do while (associated(field_3dreal_ptr))
                        field_3dreal_ptr % array(j,:,:) = real2d_temp(:,:)
                        field_3dreal_ptr => field_3dreal_ptr % next
                     end do
                  else
                     ! call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real3d_temp), real3d_temp(:,1,1))
                     call MPI_Bcast(real3d_temp(:,1,1), size(real3d_temp), MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_3dreal_ptr => field_cursor % real3dField
                     do while (associated(field_3dreal_ptr))
                        field_3dreal_ptr % array(:,:,:) = real3d_temp(:,:,:)
                        field_3dreal_ptr => field_3dreal_ptr % next
                     end do
                  end if
               end if
            end do
            reading_time = reading_time + sub_timer_sum

            if (field_cursor % real3dField % isVarArray) then
               deallocate(real2d_temp)
            else
               deallocate(real3d_temp)
            end if
call iamas_dmpar_get_time(timer_stop_time)!jwfeng
call iamas_log_write('Timing for reading field '//trim(field_cursor % real3dField % fieldName)// ': $r s', &
     realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng
call iamas_log_write('Timing for reading field kernel '//trim(field_cursor % real3dField % fieldName)// ': $r s', &
     realArgs=(/real(sub_timer_sum, kind = RKIND)/))!jwfeng
         else if (field_cursor % field_type == FIELD_4D_REAL) then

call iamas_dmpar_get_time(timer_start_time)!jwfeng
!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real3dField % fieldName))
!call iamas_log_write('DEBUGGING : reading a 4d real array')
 call iamas_log_write('FIELD_4D_REAL Reading in field '//trim(field_cursor % real4dField % fieldName))
            if (field_cursor % real4dField % isVarArray) then
!call iamas_log_write('DEBUGGING : reading a 4d real var-array')
               ncons = size(field_cursor % real4dField % constituentNames)
               allocate(real3d_temp(field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real4d_temp(field_cursor % real4dField % dimSizes(1), &
                                    field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real4dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
!call iamas_log_write('DEBUGGING : calling get_var for a constitutent')
                  call iamas_log_write('--- reading in constitutent '//trim(field_cursor % real4dField % constituentNames(j)))!jwfeng
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real4dField % constituentNames(j), real3d_temp, io_err)
               else
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real4dField % fieldName, real4d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % real4dField % isVarArray) then
                      deallocate(real3d_temp)
                   else
                      deallocate(real4d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_4dreal_ptr => field_cursor % real4dField
                  i = 1
                  do while (associated(field_4dreal_ptr))
                     if (trim(field_4dreal_ptr % dimNames(4)) == 'nCells') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, field_4dreal_ptr % dimNames(4), ownedSize)
                     end if

                     if (field_cursor % real4dField % isVarArray) then
!call iamas_log_write('DEBUGGING : copying the temporary array')
                        field_4dreal_ptr % array(j, :,:,1:ownedSize) = real3d_temp(:,:,i:i+ownedSize-1)
                     else
                        field_4dreal_ptr % array(:,:,:,1:ownedSize) = real4d_temp(:,:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_4dreal_ptr => field_4dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real4dField % isVarArray) then
                     ! call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real3d_temp), real3d_temp(:,1,1))
                     call MPI_Bcast(real3d_temp(:,1,1), size(real3d_temp), MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_4dreal_ptr => field_cursor % real4dField
                     do while (associated(field_4dreal_ptr))
                        field_4dreal_ptr % array(j,:,:,:) = real3d_temp(:,:,:)
                        field_4dreal_ptr => field_4dreal_ptr % next
                     end do
                  else
                     ! call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real4d_temp), real4d_temp(:,1,1,1))
                     call MPI_Bcast(real4d_temp(:,1,1,1), size(real4d_temp), MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_4dreal_ptr => field_cursor % real4dField
                     do while (associated(field_4dreal_ptr))
                        field_4dreal_ptr % array(:,:,:,:) = real4d_temp(:,:,:,:)
                        field_4dreal_ptr => field_4dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real4dField % isVarArray) then
               deallocate(real3d_temp)
            else
               deallocate(real4d_temp)
            end if

call iamas_dmpar_get_time(timer_stop_time)!jwfeng
call iamas_log_write('Timing for reading field '//trim(field_cursor % real4dField % fieldName)// ': $r s', &
     realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

         else if (field_cursor % field_type == FIELD_5D_REAL) then

call iamas_dmpar_get_time(timer_start_time)!jwfeng
!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real3dField % fieldName))
!call iamas_log_write('DEBUGGING : reading a 4d real array')
call iamas_log_write('FIELD_5D_REAL Reading in field '//trim(field_cursor % real5dField % fieldName))
            if (field_cursor % real5dField % isVarArray) then
!call iamas_log_write('DEBUGGING : reading a 4d real var-array')
               ncons = size(field_cursor % real5dField % constituentNames)
               allocate(real4d_temp(field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real5d_temp(field_cursor % real5dField % dimSizes(1), &
                                    field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real5dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
!call iamas_log_write('DEBUGGING : calling get_var for a constitutent')
                  call iamas_log_write('--- reading in constitutent '//trim(field_cursor % real5dField % constituentNames(j)))!jwfeng
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
               else
                  call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % real5dField % fieldName, real5d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % real5dField % isVarArray) then
                      deallocate(real4d_temp)
                   else
                      deallocate(real5d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_5dreal_ptr => field_cursor % real5dField
                  i = 1
                  do while (associated(field_5dreal_ptr))
                     if (trim(field_5dreal_ptr % dimNames(5)) == 'nCells') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, field_5dreal_ptr % dimNames(5), ownedSize)
                     end if

                     if (field_cursor % real5dField % isVarArray) then
!call iamas_log_write('DEBUGGING : copying the temporary array')
                        field_5dreal_ptr % array(j,:,:,:,1:ownedSize) = real4d_temp(:,:,:,i:i+ownedSize-1)
                     else
                        field_5dreal_ptr % array(:,:,:,:,1:ownedSize) = real5d_temp(:,:,:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_5dreal_ptr => field_5dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real5dField % isVarArray) then
                     ! call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real4d_temp), real4d_temp(:,1,1,1))
                     call MPI_Bcast(real4d_temp(:,1,1,1), size(real4d_temp), MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_5dreal_ptr => field_cursor % real5dField
                     do while (associated(field_5dreal_ptr))
                        field_5dreal_ptr % array(j,:,:,:,:) = real4d_temp(:,:,:,:)
                        field_5dreal_ptr => field_5dreal_ptr % next
                     end do
                  else
                     ! call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real5d_temp), real5d_temp(:,1,1,1,1))
                     call MPI_Bcast(real5d_temp(:,1,1,1,1), size(real5d_temp), MPI_REAL, 0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
                     field_5dreal_ptr => field_cursor % real5dField
                     do while (associated(field_5dreal_ptr))
                        field_5dreal_ptr % array(:,:,:,:,:) = real5d_temp(:,:,:,:,:)
                        field_5dreal_ptr => field_5dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real5dField % isVarArray) then
               deallocate(real4d_temp)
            else
               deallocate(real5d_temp)
            end if

call iamas_dmpar_get_time(timer_stop_time)!jwfeng
call iamas_log_write('Timing for reading field '//trim(field_cursor % real5dField % fieldName)// ': $r s', &
     realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng


         else if (field_cursor % field_type == FIELD_0D_CHAR) then

!call iamas_log_write('Reading in field '//trim(field_cursor % char0dField % fieldName))
 call iamas_log_write('FIELD_0D_CHAR Reading in field '//trim(field_cursor % char0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('MGD calling iAMAS_io_get_var now...')
            call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % char0dField % fieldName, field_cursor % char0dField % scalar, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
                if (present(ierr)) ierr = iAMAS_IO_ERR
                return
            end if

!call iamas_log_write('Distributing and Copying field to other blocks')

            ! call iamas_dmpar_bcast_char( iamas_io_handle_dminfo(stream % fileHandle), field_cursor % char0dField % scalar)
            call MPI_Bcast(field_cursor % char0dField % scalar, len(field_cursor % char0dField % scalar), MPI_CHARACTER, &
                        0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
            field_0dchar_ptr => field_cursor % char0dField
            do while (associated(field_0dchar_ptr))
               field_0dchar_ptr % scalar = field_cursor % char0dField % scalar
               field_0dchar_ptr => field_0dchar_ptr % next
            end do

         else if (field_cursor % field_type == FIELD_1D_CHAR) then
call iamas_log_write('Reading in field '//trim(field_cursor % char1dField % fieldName))
! call iamas_log_write('FIELD_1D_CHAR Reading in field '//trim(field_cursor % char1dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('MGD calling iAMAS_io_get_var now...')
            allocate(char1d_temp(field_cursor % char1dField % dimSizes(1)))
            call iAMAS_io_get_var_mul(stream % fileHandle, field_cursor % char1dField % fieldName, char1d_temp, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
                if (present(ierr)) ierr = iAMAS_IO_ERR
                deallocate(char1d_temp)
                return
            end if

!call iamas_log_write('Distributing and Copying field to other blocks')

            ! call iamas_dmpar_bcast_chars( iamas_io_handle_dminfo(stream % fileHandle), &
            !                             field_cursor % char1dField % dimSizes(1), char1d_temp)
            call MPI_Bcast(char1d_temp, field_cursor % char1dField % dimSizes(1)*len(char1d_temp(1)), MPI_CHARACTER, &
                                     0, stream % fileHandle % ioContext % dminfo % local_comm, mpi_ierr)
            field_1dchar_ptr => field_cursor % char1dField
            do while (associated(field_1dchar_ptr))
               field_1dchar_ptr % array(:)(:) = char1d_temp(:)(:)
               field_1dchar_ptr => field_1dchar_ptr % next
            end do

            deallocate(char1d_temp)
         end if
         field_cursor => field_cursor % next
      end do
!call iamas_log_write('Finished fieldlist loop...')
      call iamas_log_write('!!! Timing for reading  kernel: $r s', realArgs=(/real(reading_time, kind = RKIND)/))!jwfeng

   end subroutine iAMAS_readStream_mul


   subroutine iAMAS_readStream(stream, frame, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      integer, intent(in) :: frame
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i, j
      integer :: ncons
      integer, pointer :: ownedSize
      type (field0dInteger), pointer :: field_0dint_ptr
      type (field1dInteger), pointer :: field_1dint_ptr
      type (field2dInteger), pointer :: field_2dint_ptr
      type (field3dInteger), pointer :: field_3dint_ptr
      type (field0dReal), pointer :: field_0dreal_ptr
      type (field1dReal), pointer :: field_1dreal_ptr
      type (field2dReal), pointer :: field_2dreal_ptr
      type (field3dReal), pointer :: field_3dreal_ptr
      type (field4dReal), pointer :: field_4dreal_ptr
      type (field5dReal), pointer :: field_5dreal_ptr
      type (field0dChar), pointer :: field_0dchar_ptr
      type (field1dChar), pointer :: field_1dchar_ptr
      type (field_list_type), pointer :: field_cursor
      character (len=StrKIND), dimension(:), pointer :: char1d_temp
      integer                            :: int0d_temp
      integer, dimension(:),     pointer :: int1d_temp
      integer, dimension(:,:),   pointer :: int2d_temp
      integer, dimension(:,:,:), pointer :: int3d_temp
      real (kind=RKIND)                            :: real0d_temp
      real (kind=RKIND), dimension(:),     pointer :: real1d_temp
      real (kind=RKIND), dimension(:,:),   pointer :: real2d_temp
      real (kind=RKIND), dimension(:,:,:), pointer :: real3d_temp
      real (kind=RKIND), dimension(:,:,:,:), pointer :: real4d_temp
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: real5d_temp


      real (kind=R8KIND) :: timer_start_time, timer_stop_time!jwfeng
      real (kind=R8KIND) :: sub_timer_start_time, sub_timer_stop_time!jwfeng

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if


      !
      ! Set time frame to real
      !
      call iAMAS_io_set_frame(stream % fileHandle, frame, io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if


      !
      ! Loop over fields in the stream
      !
      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (field_cursor % field_type == FIELD_0D_INT) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % int0dField % fieldName))
!call iamas_log_write('Reading in field '//trim(field_cursor % int0dField % fieldName))
call iamas_log_write('Reading in field '//trim(field_cursor % int0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('MGD calling iAMAS_io_get_var now...')
            call iAMAS_io_get_var(stream % fileHandle, field_cursor % int0dField % fieldName, int0d_temp, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
                if (present(ierr)) ierr = iAMAS_IO_ERR
                return
            end if

!call iamas_log_write('Distributing and Copying field to other blocks')

            call iamas_dmpar_bcast_int( iamas_io_handle_dminfo(stream % fileHandle), int0d_temp)
            field_0dint_ptr => field_cursor % int0dField
            do while (associated(field_0dint_ptr))
               field_0dint_ptr % scalar = int0d_temp
               field_0dint_ptr => field_0dint_ptr % next
            end do

         else if (field_cursor % field_type == FIELD_1D_INT) then
!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % int1dField % fieldName))

            if (field_cursor % int1dField % isVarArray) then
               ncons = size(field_cursor % int1dField % constituentNames)
            else
               ncons = 1
               allocate(int1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % int1dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % int1dField % constituentNames(j), int0d_temp, io_err)
               else
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % int1dField % fieldName, int1d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (.not. field_cursor % int1dField % isVarArray) then
                      deallocate(int1d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_1dint_ptr => field_cursor % int1dField
                  i = 1
                  do while (associated(field_1dint_ptr))
                     if (trim(field_1dint_ptr % dimNames(1)) == 'nCells') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, field_1dint_ptr % dimNames(1), ownedSize)
                     end if

                     if (field_cursor % int1dField % isVarArray) then
                        field_1dint_ptr % array(j) = int0d_temp
                     else
                        field_1dint_ptr % array(1:ownedSize) = int1d_temp(i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_1dint_ptr => field_1dint_ptr % next
                  end do

               else
   
                  if (field_cursor % int1dField % isVarArray) then
                     call iamas_dmpar_bcast_int( iamas_io_handle_dminfo(stream % fileHandle), int0d_temp)
                     field_1dint_ptr => field_cursor % int1dField
                     do while (associated(field_1dint_ptr))
                        field_1dint_ptr % array(j) = int0d_temp
                        field_1dint_ptr => field_1dint_ptr % next
                     end do
                  else
                     call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int1d_temp), int1d_temp(:))
                     field_1dint_ptr => field_cursor % int1dField
                     do while (associated(field_1dint_ptr))
                        field_1dint_ptr % array(:) = int1d_temp(:)
                        field_1dint_ptr => field_1dint_ptr % next
                     end do
                  end if
               end if
            end do

            if (.not. field_cursor % int1dField % isVarArray) then
               deallocate(int1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_INT) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % int2dField % fieldName))
call iamas_log_write('Reading in field '//trim(field_cursor % int2dField % fieldName))
            if (field_cursor % int2dField % isVarArray) then
               ncons = size(field_cursor % int2dField % constituentNames)
               allocate(int1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int2d_temp(field_cursor % int2dField % dimSizes(1), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % int2dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % int2dField % constituentNames(j), int1d_temp, io_err)
               else
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % int2dField % fieldName, int2d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % int2dField % isVarArray) then
                      deallocate(int1d_temp)
                   else
                      deallocate(int2d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_2dint_ptr => field_cursor % int2dField
                  i = 1
                  do while (associated(field_2dint_ptr))
                     if (trim(field_2dint_ptr % dimNames(2)) == 'nCells') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, field_2dint_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % int2dField % isVarArray) then
                        field_2dint_ptr % array(j,1:ownedSize) = int1d_temp(i:i+ownedSize-1)
                     else
                        field_2dint_ptr % array(:,1:ownedSize) = int2d_temp(:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_2dint_ptr => field_2dint_ptr % next
                  end do

               else
   
                  if (field_cursor % int2dField % isVarArray) then
                     call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int1d_temp), int1d_temp(:))
                     field_2dint_ptr => field_cursor % int2dField
                     do while (associated(field_2dint_ptr))
                        field_2dint_ptr % array(j,:) = int1d_temp(:)
                        field_2dint_ptr => field_2dint_ptr % next
                     end do
                  else
                     call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int2d_temp), int2d_temp(:,1))
                     field_2dint_ptr => field_cursor % int2dField
                     do while (associated(field_2dint_ptr))
                        field_2dint_ptr % array(:,:) = int2d_temp(:,:)
                        field_2dint_ptr => field_2dint_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % int2dField % isVarArray) then
               deallocate(int1d_temp)
            else
               deallocate(int2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_INT) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % int3dField % fieldName))
call iamas_log_write('Reading in field '//trim(field_cursor % int3dField % fieldName))
            if (field_cursor % int3dField % isVarArray) then
               ncons = size(field_cursor % int3dField % constituentNames)
               allocate(int2d_temp(field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int3d_temp(field_cursor % int3dField % dimSizes(1), &
                                   field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % int3dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % int3dField % constituentNames(j), int2d_temp, io_err)
               else
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % int3dField % fieldName, int3d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % int3dField % isVarArray) then
                      deallocate(int2d_temp)
                   else
                      deallocate(int3d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_3dint_ptr => field_cursor % int3dField
                  i = 1
                  do while (associated(field_3dint_ptr))
                     if (trim(field_3dint_ptr % dimNames(3)) == 'nCells') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, field_3dint_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % int3dField % isVarArray) then
                        field_3dint_ptr % array(j,:,1:ownedSize) = int2d_temp(:,i:i+ownedSize-1)
                     else
                        field_3dint_ptr % array(:,:,1:ownedSize) = int3d_temp(:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_3dint_ptr => field_3dint_ptr % next
                  end do

               else
   
                  if (field_cursor % int3dField % isVarArray) then
                     call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int2d_temp), int2d_temp(:,1))
                     field_3dint_ptr => field_cursor % int3dField
                     do while (associated(field_3dint_ptr))
                        field_3dint_ptr % array(j,:,:) = int2d_temp(:,:)
                        field_3dint_ptr => field_3dint_ptr % next
                     end do
                  else
                     call iamas_dmpar_bcast_ints( iamas_io_handle_dminfo(stream % fileHandle), size(int3d_temp), int3d_temp(:,1,1))
                     field_3dint_ptr => field_cursor % int3dField
                     do while (associated(field_3dint_ptr))
                        field_3dint_ptr % array(:,:,:) = int3d_temp(:,:,:)
                        field_3dint_ptr => field_3dint_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % int3dField % isVarArray) then
               deallocate(int2d_temp)
            else
               deallocate(int3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_0D_REAL) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real0dField % fieldName))
!call iamas_log_write('Reading in field '//trim(field_cursor % real0dField % fieldName))
call iamas_log_write('Reading in field '//trim(field_cursor % real0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/) )

!call iamas_log_write('MGD calling iAMAS_io_get_var now...')
            call iAMAS_io_get_var(stream % fileHandle, field_cursor % real0dField % fieldName, real0d_temp, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
                if (present(ierr)) ierr = iAMAS_IO_ERR
                return
            end if

!call iamas_log_write('Distributing and Copying field to other blocks')

            call iamas_dmpar_bcast_real( iamas_io_handle_dminfo(stream % fileHandle), real0d_temp)
            field_0dreal_ptr => field_cursor % real0dField
            do while (associated(field_0dreal_ptr))
               field_0dreal_ptr % scalar = real0d_temp
               field_0dreal_ptr => field_0dreal_ptr % next
            end do

         else if (field_cursor % field_type == FIELD_1D_REAL) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real1dField % fieldName))
call iamas_log_write('Reading in field '//trim(field_cursor % real1dField % fieldName))
            if (field_cursor % real1dField % isVarArray) then
               ncons = size(field_cursor % real1dField % constituentNames)
            else
               ncons = 1
               allocate(real1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real1dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real1dField % constituentNames(j), real0d_temp, io_err)
               else
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real1dField % fieldName, real1d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (.not. field_cursor % real1dField % isVarArray) then
                      deallocate(real1d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_1dreal_ptr => field_cursor % real1dField
                  i = 1

                  do while (associated(field_1dreal_ptr))
                     if (trim(field_1dreal_ptr % dimNames(1)) == 'nCells') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, field_1dreal_ptr % dimNames(1), ownedSize)
                     end if

                     if (field_cursor % real1dField % isVarArray) then
                        field_1dreal_ptr % array(j) = real0d_temp
                     else
                        field_1dreal_ptr % array(1:ownedSize) = real1d_temp(i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_1dreal_ptr => field_1dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real1dField % isVarArray) then
                     call iamas_dmpar_bcast_real( iamas_io_handle_dminfo(stream % fileHandle), real0d_temp)
                     field_1dreal_ptr => field_cursor % real1dField
                     do while (associated(field_1dreal_ptr))
                        field_1dreal_ptr % array(j) = real0d_temp
                        field_1dreal_ptr => field_1dreal_ptr % next
                     end do
                  else
                     call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real1d_temp), real1d_temp(:))
                     field_1dreal_ptr => field_cursor % real1dField
                     do while (associated(field_1dreal_ptr))
                        field_1dreal_ptr % array(:) = real1d_temp(:)
                        field_1dreal_ptr => field_1dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (.not. field_cursor % real1dField % isVarArray) then
               deallocate(real1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_REAL) then

!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real2dField % fieldName))
call iamas_dmpar_get_time(timer_start_time)!jwfeng
call iamas_log_write('Reading in field '//trim(field_cursor % real2dField % fieldName))
            if (field_cursor % real2dField % isVarArray) then
               ncons = size(field_cursor % real2dField % constituentNames)
               allocate(real1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real2d_temp(field_cursor % real2dField % dimSizes(1), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real2dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
                  call iamas_log_write('--- reading in constitutent '//trim(field_cursor % real2dField % constituentNames(j)))!jwfeng
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real2dField % constituentNames(j), real1d_temp, io_err)
               else
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real2dField % fieldName, real2d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % real2dField % isVarArray) then
                      deallocate(real1d_temp)
                   else
                      deallocate(real2d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_2dreal_ptr => field_cursor % real2dField
                  i = 1
                  do while (associated(field_2dreal_ptr))
                     if (trim(field_2dreal_ptr % dimNames(2)) == 'nCells') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, field_2dreal_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % real2dField % isVarArray) then
                        field_2dreal_ptr % array(j,1:ownedSize) = real1d_temp(i:i+ownedSize-1)
                     else
                        field_2dreal_ptr % array(:,1:ownedSize) = real2d_temp(:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_2dreal_ptr => field_2dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real2dField % isVarArray) then
                     call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real1d_temp), real1d_temp(:))
                     field_2dreal_ptr => field_cursor % real2dField
                     do while (associated(field_2dreal_ptr))
                        field_2dreal_ptr % array(j,:) = real1d_temp(:)
                        field_2dreal_ptr => field_2dreal_ptr % next
                     end do
                  else
                     call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real2d_temp), real2d_temp(:,1))
                     field_2dreal_ptr => field_cursor % real2dField
                     do while (associated(field_2dreal_ptr))
                        field_2dreal_ptr % array(:,:) = real2d_temp(:,:)
                        field_2dreal_ptr => field_2dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real2dField % isVarArray) then
               deallocate(real1d_temp)
            else
               deallocate(real2d_temp)
            end if
call iamas_dmpar_get_time(timer_stop_time)!jwfeng
call iamas_log_write('Timing for reading field '//trim(field_cursor % real2dField % fieldName)// ': $r s', &
     realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

         else if (field_cursor % field_type == FIELD_3D_REAL) then

call iamas_dmpar_get_time(timer_start_time)!jwfeng
!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real3dField % fieldName))
!call iamas_log_write('DEBUGGING : reading a 3d real array')
call iamas_log_write('Reading in field '//trim(field_cursor % real3dField % fieldName))
            if (field_cursor % real3dField % isVarArray) then
!call iamas_log_write('DEBUGGING : reading a 3d real var-array')
               ncons = size(field_cursor % real3dField % constituentNames)
               allocate(real2d_temp(field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real3d_temp(field_cursor % real3dField % dimSizes(1), &
                                    field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real3dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
!call iamas_log_write('DEBUGGING : calling get_var for a constitutent')
                  call iamas_log_write('--- reading in constitutent '//trim(field_cursor % real3dField % constituentNames(j)))!jwfeng
                  call iamas_dmpar_get_time(sub_timer_start_time)!jwfeng
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real3dField % constituentNames(j), real2d_temp, io_err)
                  call iamas_dmpar_get_time(sub_timer_stop_time)!jwfeng
                  call iamas_log_write('--- Timing for reading constitutent '//trim(field_cursor % real3dField % constituentNames(j))// ': $r s', &
                       realArgs=(/real(sub_timer_stop_time - sub_timer_start_time, kind = RKIND)/))!jwfeng
               else
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real3dField % fieldName, real3d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % real3dField % isVarArray) then
                      deallocate(real2d_temp)
                   else
                      deallocate(real3d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_3dreal_ptr => field_cursor % real3dField
                  i = 1
                  do while (associated(field_3dreal_ptr))
                     if (trim(field_3dreal_ptr % dimNames(3)) == 'nCells') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, field_3dreal_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % real3dField % isVarArray) then
!call iamas_log_write('DEBUGGING : copying the temporary array')
                        field_3dreal_ptr % array(j,:,1:ownedSize) = real2d_temp(:,i:i+ownedSize-1)
                     else
                        field_3dreal_ptr % array(:,:,1:ownedSize) = real3d_temp(:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_3dreal_ptr => field_3dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real3dField % isVarArray) then
                     call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real2d_temp), real2d_temp(:,1))
                     field_3dreal_ptr => field_cursor % real3dField
                     do while (associated(field_3dreal_ptr))
                        field_3dreal_ptr % array(j,:,:) = real2d_temp(:,:)
                        field_3dreal_ptr => field_3dreal_ptr % next
                     end do
                  else
                     call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real3d_temp), real3d_temp(:,1,1))
                     field_3dreal_ptr => field_cursor % real3dField
                     do while (associated(field_3dreal_ptr))
                        field_3dreal_ptr % array(:,:,:) = real3d_temp(:,:,:)
                        field_3dreal_ptr => field_3dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real3dField % isVarArray) then
               deallocate(real2d_temp)
            else
               deallocate(real3d_temp)
            end if
call iamas_dmpar_get_time(timer_stop_time)!jwfeng
call iamas_log_write('Timing for reading field '//trim(field_cursor % real3dField % fieldName)// ': $r s', &
     realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng
         else if (field_cursor % field_type == FIELD_4D_REAL) then

call iamas_dmpar_get_time(timer_start_time)!jwfeng
!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real3dField % fieldName))
!call iamas_log_write('DEBUGGING : reading a 4d real array')
call iamas_log_write('Reading in field '//trim(field_cursor % real4dField % fieldName))
            if (field_cursor % real4dField % isVarArray) then
!call iamas_log_write('DEBUGGING : reading a 4d real var-array')
               ncons = size(field_cursor % real4dField % constituentNames)
               allocate(real3d_temp(field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real4d_temp(field_cursor % real4dField % dimSizes(1), &
                                    field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real4dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
!call iamas_log_write('DEBUGGING : calling get_var for a constitutent')
                  call iamas_log_write('--- reading in constitutent '//trim(field_cursor % real4dField % constituentNames(j)))!jwfeng
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real4dField % constituentNames(j), real3d_temp, io_err)
               else
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real4dField % fieldName, real4d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % real4dField % isVarArray) then
                      deallocate(real3d_temp)
                   else
                      deallocate(real4d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_4dreal_ptr => field_cursor % real4dField
                  i = 1
                  do while (associated(field_4dreal_ptr))
                     if (trim(field_4dreal_ptr % dimNames(4)) == 'nCells') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, field_4dreal_ptr % dimNames(4), ownedSize)
                     end if

                     if (field_cursor % real4dField % isVarArray) then
!call iamas_log_write('DEBUGGING : copying the temporary array')
                        field_4dreal_ptr % array(j, :,:,1:ownedSize) = real3d_temp(:,:,i:i+ownedSize-1)
                     else
                        field_4dreal_ptr % array(:,:,:,1:ownedSize) = real4d_temp(:,:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_4dreal_ptr => field_4dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real4dField % isVarArray) then
                     call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real3d_temp), real3d_temp(:,1,1))
                     field_4dreal_ptr => field_cursor % real4dField
                     do while (associated(field_4dreal_ptr))
                        field_4dreal_ptr % array(j,:,:,:) = real3d_temp(:,:,:)
                        field_4dreal_ptr => field_4dreal_ptr % next
                     end do
                  else
                     call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real4d_temp), real4d_temp(:,1,1,1))
                     field_4dreal_ptr => field_cursor % real4dField
                     do while (associated(field_4dreal_ptr))
                        field_4dreal_ptr % array(:,:,:,:) = real4d_temp(:,:,:,:)
                        field_4dreal_ptr => field_4dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real4dField % isVarArray) then
               deallocate(real3d_temp)
            else
               deallocate(real4d_temp)
            end if

call iamas_dmpar_get_time(timer_stop_time)!jwfeng
call iamas_log_write('Timing for reading field '//trim(field_cursor % real4dField % fieldName)// ': $r s', &
     realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

         else if (field_cursor % field_type == FIELD_5D_REAL) then

call iamas_dmpar_get_time(timer_start_time)!jwfeng
!call iamas_log_write('DEBUGGING : *************** '//trim(field_cursor % real3dField % fieldName))
!call iamas_log_write('DEBUGGING : reading a 4d real array')
call iamas_log_write('Reading in field '//trim(field_cursor % real5dField % fieldName))
            if (field_cursor % real5dField % isVarArray) then
!call iamas_log_write('DEBUGGING : reading a 4d real var-array')
               ncons = size(field_cursor % real5dField % constituentNames)
               allocate(real4d_temp(field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real5d_temp(field_cursor % real5dField % dimSizes(1), &
                                    field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % real5dField % isVarArray) then
                  if (.not. field_cursor % isAvailable(j)) cycle
!call iamas_log_write('DEBUGGING : calling get_var for a constitutent')
                  call iamas_log_write('--- reading in constitutent '//trim(field_cursor % real5dField % constituentNames(j)))!jwfeng
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
               else
                  call iAMAS_io_get_var(stream % fileHandle, field_cursor % real5dField % fieldName, real5d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR) then
                   if (present(ierr)) ierr = iAMAS_IO_ERR
                   if (field_cursor % real5dField % isVarArray) then
                      deallocate(real4d_temp)
                   else
                      deallocate(real5d_temp)
                   end if
                   return
               end if
   
               if (field_cursor % isDecomposed) then
                  ! Distribute field to multiple blocks
                  field_5dreal_ptr => field_cursor % real5dField
                  i = 1
                  do while (associated(field_5dreal_ptr))
                     if (trim(field_5dreal_ptr % dimNames(5)) == 'nCells') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, field_5dreal_ptr % dimNames(5), ownedSize)
                     end if

                     if (field_cursor % real5dField % isVarArray) then
!call iamas_log_write('DEBUGGING : copying the temporary array')
                        field_5dreal_ptr % array(j,:,:,:,1:ownedSize) = real4d_temp(:,:,:,i:i+ownedSize-1)
                     else
                        field_5dreal_ptr % array(:,:,:,:,1:ownedSize) = real5d_temp(:,:,:,:,i:i+ownedSize-1)
                     end if
                     i = i + ownedSize
                     field_5dreal_ptr => field_5dreal_ptr % next
                  end do

               else
   
                  if (field_cursor % real5dField % isVarArray) then
                     call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real4d_temp), real4d_temp(:,1,1,1))
                     field_5dreal_ptr => field_cursor % real5dField
                     do while (associated(field_5dreal_ptr))
                        field_5dreal_ptr % array(j,:,:,:,:) = real4d_temp(:,:,:,:)
                        field_5dreal_ptr => field_5dreal_ptr % next
                     end do
                  else
                     call iamas_dmpar_bcast_reals( iamas_io_handle_dminfo(stream % fileHandle), size(real5d_temp), real5d_temp(:,1,1,1,1))
                     field_5dreal_ptr => field_cursor % real5dField
                     do while (associated(field_5dreal_ptr))
                        field_5dreal_ptr % array(:,:,:,:,:) = real5d_temp(:,:,:,:,:)
                        field_5dreal_ptr => field_5dreal_ptr % next
                     end do
                  end if
               end if
            end do

            if (field_cursor % real5dField % isVarArray) then
               deallocate(real4d_temp)
            else
               deallocate(real5d_temp)
            end if

call iamas_dmpar_get_time(timer_stop_time)!jwfeng
call iamas_log_write('Timing for reading field '//trim(field_cursor % real5dField % fieldName)// ': $r s', &
     realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng


         else if (field_cursor % field_type == FIELD_0D_CHAR) then

!call iamas_log_write('Reading in field '//trim(field_cursor % char0dField % fieldName))
call iamas_log_write('Reading in field '//trim(field_cursor % char0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('MGD calling iAMAS_io_get_var now...')
            call iAMAS_io_get_var(stream % fileHandle, field_cursor % char0dField % fieldName, field_cursor % char0dField % scalar, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
                if (present(ierr)) ierr = iAMAS_IO_ERR
                return
            end if

!call iamas_log_write('Distributing and Copying field to other blocks')

            call iamas_dmpar_bcast_char( iamas_io_handle_dminfo(stream % fileHandle), field_cursor % char0dField % scalar)
            field_0dchar_ptr => field_cursor % char0dField
            do while (associated(field_0dchar_ptr))
               field_0dchar_ptr % scalar = field_cursor % char0dField % scalar
               field_0dchar_ptr => field_0dchar_ptr % next
            end do

         else if (field_cursor % field_type == FIELD_1D_CHAR) then
!call iamas_log_write('Reading in field '//trim(field_cursor % char1dField % fieldName))
call iamas_log_write('Reading in field '//trim(field_cursor % char1dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('MGD calling iAMAS_io_get_var now...')
            allocate(char1d_temp(field_cursor % char1dField % dimSizes(1)))
            call iAMAS_io_get_var(stream % fileHandle, field_cursor % char1dField % fieldName, char1d_temp, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR) then
                if (present(ierr)) ierr = iAMAS_IO_ERR
                deallocate(char1d_temp)
                return
            end if

!call iamas_log_write('Distributing and Copying field to other blocks')

            call iamas_dmpar_bcast_chars( iamas_io_handle_dminfo(stream % fileHandle), &
                                        field_cursor % char1dField % dimSizes(1), char1d_temp)
            field_1dchar_ptr => field_cursor % char1dField
            do while (associated(field_1dchar_ptr))
               field_1dchar_ptr % array(:)(:) = char1d_temp(:)(:)
               field_1dchar_ptr => field_1dchar_ptr % next
            end do

            deallocate(char1d_temp)
         end if
         field_cursor => field_cursor % next
      end do
!call iamas_log_write('Finished fieldlist loop...')

   end subroutine iAMAS_readStream


   subroutine iAMAS_writeStream_mul(stream, frame, arrayFieldList, ierr)

      use pio
      use mpi
      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      integer, intent(in) :: frame
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(out), optional :: ierr

      character (len=StrKIND) :: dimensionSuffix
      integer :: io_err
      integer :: i, j, k, iVarArray
      integer :: ncons
      integer, pointer :: ownedSize
      type (field0dInteger), pointer :: field_0dint_ptr
      type (field1dInteger), pointer :: field_1dint_ptr
      type (field2dInteger), pointer :: field_2dint_ptr
      type (field3dInteger), pointer :: field_3dint_ptr
      type (field0dReal), pointer :: field_0dreal_ptr
      type (field1dReal), pointer :: field_1dreal_ptr
      type (field2dReal), pointer :: field_2dreal_ptr
      type (field3dReal), pointer :: field_3dreal_ptr
      type (field4dReal), pointer :: field_4dreal_ptr
      type (field5dReal), pointer :: field_5dreal_ptr
      type (field0dChar), pointer :: field_0dchar_ptr
      type (field1dChar), pointer :: field_1dchar_ptr
      type (field_list_type), pointer :: field_cursor
      integer                            :: int0d_temp
      integer, dimension(:),     pointer :: int1d_temp
      integer, dimension(:,:),   pointer :: int2d_temp
      integer, dimension(:,:,:), pointer :: int3d_temp
      real (kind=RKIND)                                :: real0d_temp
      real (kind=RKIND), dimension(:),         pointer :: real1d_temp
      real (kind=RKIND), dimension(:,:),       pointer :: real2d_temp
      real (kind=RKIND), dimension(:,:,:),     pointer :: real3d_temp
      real (kind=RKIND), dimension(:,:,:,:),   pointer :: real4d_temp
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: real5d_temp
      
      integer :: local_rank,mpi_ierr
      real (kind=R8KIND) :: timer_start_time, timer_stop_time
      real (kind=R8KIND) :: sub_timer_start_time, sub_timer_stop_time
      real (kind=R8KIND) :: sub_timer_sum, left_timer_sum
      real (kind=R8KIND) :: sync_timer_start_time, sync_timer_stop_time

      call MPI_Comm_rank(stream % fileHandle % ioContext % dminfo % local_comm, local_rank, mpi_ierr)
      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      !
      ! Set time frame to write
      !
      call iAMAS_io_set_frame(stream % fileHandle, frame, io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

      !
      ! Check whether we will clobber any records
      !
      if (iAMAS_io_would_clobber_records(stream % fileHandle, io_err)) then
         if (.not. stream % clobberRecords) then
            if (present(ierr)) ierr = iAMAS_STREAM_CLOBBER_RECORD
            return
         else
            call iamas_log_write('iAMAS I/O: Overwriting existing record $i in output file '// &
                  trim(stream % filename), iAMAS_LOG_WARN, intArgs=(/frame/))
         end if
      end if

      if ( stream % blockWrite ) then
         dimensionSuffix = ''
      else
         dimensionSuffix = 'Solve'
      end if

      !
      ! Loop over fields in the stream
      !
      field_cursor => stream % fieldList
      sub_timer_sum = 0
      left_timer_sum = 0
      do while (associated(field_cursor))

         if (field_cursor % field_type == FIELD_0D_INT) then

!call iamas_log_write('Writing out field '//trim(field_cursor % int0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('Copying field from first block')
            int0d_temp = field_cursor % int0dField % scalar

!call iamas_log_write('MGD calling iAMAS_io_put_var_mul now...')
            call iamas_dmpar_get_time(sub_timer_start_time)
            call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % int0dField % fieldName, int0d_temp, io_err)
            call iamas_dmpar_get_time(sub_timer_stop_time)
            sub_timer_sum = sub_timer_sum + (sub_timer_stop_time - sub_timer_start_time)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_INT) then

            if (field_cursor % int1dField % isVarArray) then
               ncons = size(field_cursor % int1dField % constituentNames)
            else
               ncons = 1
               allocate(int1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_1dint_ptr => field_cursor % int1dField
                  i = 1
                  do while (associated(field_1dint_ptr))
                     if (trim(field_1dint_ptr % dimNames(1)) == 'nCells') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, trim(field_1dint_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, trim(field_1dint_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, trim(field_1dint_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, field_1dint_ptr % dimNames(1), ownedSize)
                     end if

                     if (field_cursor % int1dField % isVarArray) then
! I suspect we will never hit this code, as it doesn't make sense, really
                        int0d_temp = field_1dint_ptr % array(j)
                     else
                        int1d_temp(i:i+ownedSize-1) = field_1dint_ptr % array(1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_1dint_ptr => field_1dint_ptr % next
                     else
                        nullify(field_1dint_ptr)
                     end if
                  end do
               else
                  if (field_cursor % int1dField % isVarArray) then
                     int0d_temp = field_cursor % int1dField % array(j)
                  else
                     int1d_temp(:) = field_cursor % int1dField % array(:)
                  end if
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)

               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % int1dField % isVarArray) then
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % int1dField % constituentNames(j), int0d_temp, io_err)
               else
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % int1dField % fieldName, int1d_temp, io_err)
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               sub_timer_sum = sub_timer_sum + (sub_timer_stop_time - sub_timer_start_time)
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (.not. field_cursor % int1dField % isVarArray) then
               deallocate(int1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_INT) then

            if (field_cursor % int2dField % isVarArray) then
               ncons = size(field_cursor % int2dField % constituentNames)
               allocate(int1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int2d_temp(field_cursor % int2dField % dimSizes(1), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_2dint_ptr => field_cursor % int2dField
                  i = 1
                  do while (associated(field_2dint_ptr))
                     if (trim(field_2dint_ptr % dimNames(2)) == 'nCells') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, trim(field_2dint_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, trim(field_2dint_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, trim(field_2dint_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, field_2dint_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % int2dField % isVarArray) then
                        int1d_temp(i:i+ownedSize-1) = field_2dint_ptr % array(j,1:ownedSize)
                     else
                        int2d_temp(:,i:i+ownedSize-1) = field_2dint_ptr % array(:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_2dint_ptr => field_2dint_ptr % next
                     else
                        nullify(field_2dint_ptr)
                     end if
                  end do
               else
                  if (field_cursor % int2dField % isVarArray) then
                     int1d_temp(:) = field_cursor % int2dField % array(j,:)
                  else
                     int2d_temp(:,:) = field_cursor % int2dField % array(:,:)
                  end if
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)

               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % int2dField % isVarArray) then
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % int2dField % constituentNames(j), int1d_temp, io_err)
               else
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % int2dField % fieldName, int2d_temp, io_err)
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               sub_timer_sum = sub_timer_sum + (sub_timer_stop_time - sub_timer_start_time)
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % int2dField % isVarArray) then
               deallocate(int1d_temp)
            else
               deallocate(int2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_INT) then

            call iamas_dmpar_get_time(sub_timer_start_time)
            if (field_cursor % int3dField % isVarArray) then
               ncons = size(field_cursor % int3dField % constituentNames)
               allocate(int2d_temp(field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int3d_temp(field_cursor % int3dField % dimSizes(1), &
                                   field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_3dint_ptr => field_cursor % int3dField
                  i = 1
                  do while (associated(field_3dint_ptr))
                     if (trim(field_3dint_ptr % dimNames(3)) == 'nCells') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, trim(field_3dint_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, trim(field_3dint_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, trim(field_3dint_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, field_3dint_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % int3dField % isVarArray) then
                        int2d_temp(:,i:i+ownedSize-1) = field_3dint_ptr % array(j,:,1:ownedSize)
                     else
                        int3d_temp(:,:,i:i+ownedSize-1) = field_3dint_ptr % array(:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_3dint_ptr => field_3dint_ptr % next
                     else
                        nullify(field_3dint_ptr)
                     end if
                  end do
               else
                  if (field_cursor % int3dField % isVarArray) then
                     int2d_temp(:,:) = field_cursor % int3dField % array(j,:,:)
                  else
                     int3d_temp(:,:,:) = field_cursor % int3dField % array(:,:,:)
                  end if
               end if

               if (field_cursor % int3dField % isVarArray) then
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % int3dField % constituentNames(j), int2d_temp, io_err)
               else
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % int3dField % fieldName, int3d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % int3dField % isVarArray) then
               deallocate(int2d_temp)
            else
               deallocate(int3d_temp)
            end if
            call iamas_dmpar_get_time(sub_timer_stop_time)
            left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)

         else if (field_cursor % field_type == FIELD_0D_REAL) then

!call iamas_log_write('Writing out field '//trim(field_cursor % real0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('Copying field from first block')
            call iamas_dmpar_get_time(sub_timer_start_time)
            real0d_temp = field_cursor % real0dField % scalar

!call iamas_log_write('MGD calling iAMAS_io_put_var_mul now...')
            call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real0dField % fieldName, real0d_temp, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            call iamas_dmpar_get_time(sub_timer_stop_time)
            left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)

         else if (field_cursor % field_type == FIELD_1D_REAL) then

            if (field_cursor % real1dField % isVarArray) then
               ncons = size(field_cursor % real1dField % constituentNames)
            else
               ncons = 1
               allocate(real1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_1dreal_ptr => field_cursor % real1dField
                  i = 1
                  do while (associated(field_1dreal_ptr))
                     if (trim(field_1dreal_ptr % dimNames(1)) == 'nCells') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, trim(field_1dreal_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, trim(field_1dreal_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, trim(field_1dreal_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, field_1dreal_ptr % dimNames(1), ownedSize)
                     end if

                     if (field_cursor % real1dField % isVarArray) then
! I suspect we will never hit this code, as it doesn't make sense, really
                        real0d_temp = field_1dreal_ptr % array(j)
                     else
                        real1d_temp(i:i+ownedSize-1) = field_1dreal_ptr % array(1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_1dreal_ptr => field_1dreal_ptr % next
                     else
                        nullify(field_1dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real1dField % isVarArray) then
                     real0d_temp = field_cursor % real1dField % array(j)
                  else
                     real1d_temp(:) = field_cursor % real1dField % array(:)
                  end if
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)

               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % real1dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real1dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real1dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real1dField % constituentNames(j), real0d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real1dField % constituentNames(j), real0d_temp, io_err)
                  end if

               else
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real1dField % fieldName, real1d_temp, io_err)
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               sub_timer_sum = sub_timer_sum + (sub_timer_stop_time - sub_timer_start_time)
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (.not. field_cursor % real1dField % isVarArray) then
               deallocate(real1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_REAL) then

            if (field_cursor % real2dField % isVarArray) then
               ncons = size(field_cursor % real2dField % constituentNames)
               allocate(real1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real2d_temp(field_cursor % real2dField % dimSizes(1), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_2dreal_ptr => field_cursor % real2dField
                  i = 1
                  do while (associated(field_2dreal_ptr))
                     if (trim(field_2dreal_ptr % dimNames(2)) == 'nCells') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, trim(field_2dreal_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, trim(field_2dreal_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, trim(field_2dreal_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, field_2dreal_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % real2dField % isVarArray) then
                        real1d_temp(i:i+ownedSize-1) = field_2dreal_ptr % array(j,1:ownedSize)
                     else
                        real2d_temp(:,i:i+ownedSize-1) = field_2dreal_ptr % array(:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_2dreal_ptr => field_2dreal_ptr % next
                     else
                        nullify(field_2dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real2dField % isVarArray) then
                     real1d_temp(:) = field_cursor % real2dField % array(j,:)
                  else
                     real2d_temp(:,:) = field_cursor % real2dField % array(:,:)
                  end if
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)

               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % real2dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real2dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real2dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real2dField % constituentNames(j), real1d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real2dField % constituentNames(j), real1d_temp, io_err)
                  end if

               else
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real2dField % fieldName, real2d_temp, io_err)
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               sub_timer_sum = sub_timer_sum + (sub_timer_stop_time - sub_timer_start_time)
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % real2dField % isVarArray) then
               deallocate(real1d_temp)
            else
               deallocate(real2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_REAL) then
            ! call iamas_log_write('====== field 3d real, rank  $i', intArgs=(/local_rank/))

            if (field_cursor % real3dField % isVarArray) then
               ncons = size(field_cursor % real3dField % constituentNames)
               allocate(real2d_temp(field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
               ! call iamas_log_write('====== 3d real dimSizes(1) is   $i', intArgs=(/field_cursor % real3dField % dimSizes(1)/))
               ! call iamas_log_write('====== 3d real dimSizes(2) is   $i', intArgs=(/field_cursor % real3dField % dimSizes(2)/))
               ! call iamas_log_write('====== 3d real totalDimSize is   $i', intArgs=(/field_cursor % totalDimSize/))
            else
               ncons = 1
               allocate(real3d_temp(field_cursor % real3dField % dimSizes(1), &
                                    field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_3dreal_ptr => field_cursor % real3dField
                  i = 1
                  do while (associated(field_3dreal_ptr))
                     if (trim(field_3dreal_ptr % dimNames(3)) == 'nCells') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, trim(field_3dreal_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, trim(field_3dreal_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, trim(field_3dreal_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, field_3dreal_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % real3dField % isVarArray) then
                        real2d_temp(:,i:i+ownedSize-1) = field_3dreal_ptr % array(j,:,1:ownedSize)
                        ! do k = i, i + ownedSize - 1
                        !    real2d_temp(:, k) = local_rank
                        ! end do
                     else
                        real3d_temp(:,:,i:i+ownedSize-1) = field_3dreal_ptr % array(:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_3dreal_ptr => field_3dreal_ptr % next
                     else
                        nullify(field_3dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real3dField % isVarArray) then
                     real2d_temp(:,:) = field_cursor % real3dField % array(j,:,:)
                  else
                     real3d_temp(:,:,:) = field_cursor % real3dField % array(:,:,:)
                  end if
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)

               call iamas_dmpar_get_time(sub_timer_start_time)
               if (field_cursor % real3dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real3dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real3dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real3dField % constituentNames(j), real2d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real3dField % constituentNames(j), real2d_temp, io_err)
                  end if

               else
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real3dField % fieldName, real3d_temp, io_err)
               end if
               call iamas_dmpar_get_time(sub_timer_stop_time)
               sub_timer_sum = sub_timer_sum + (sub_timer_stop_time - sub_timer_start_time)
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % real3dField % isVarArray) then
               deallocate(real2d_temp)
            else
               deallocate(real3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_4D_REAL) then

            call iamas_dmpar_get_time(sub_timer_start_time)
            if (field_cursor % real4dField % isVarArray) then
               ncons = size(field_cursor % real4dField % constituentNames)
               allocate(real3d_temp(field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real4d_temp(field_cursor % real4dField % dimSizes(1), &
                                    field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_4dreal_ptr => field_cursor % real4dField
                  i = 1
                  do while (associated(field_4dreal_ptr))
                     if (trim(field_4dreal_ptr % dimNames(4)) == 'nCells') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, trim(field_4dreal_ptr % dimNames(4)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, trim(field_4dreal_ptr % dimNames(4)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, trim(field_4dreal_ptr % dimNames(4)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, field_4dreal_ptr % dimNames(4), ownedSize)
                     end if

                     if (field_cursor % real4dField % isVarArray) then
                        real3d_temp(:,:,i:i+ownedSize-1) = field_4dreal_ptr % array(j,:,:,1:ownedSize)
                     else
                        real4d_temp(:,:,:,i:i+ownedSize-1) = field_4dreal_ptr % array(:,:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_4dreal_ptr => field_4dreal_ptr % next
                     else
                        nullify(field_4dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real4dField % isVarArray) then
                     real3d_temp(:,:,:) = field_cursor % real4dField % array(j,:,:,:)
                  else
                     real4d_temp(:,:,:,:) = field_cursor % real4dField % array(:,:,:,:)
                  end if
               end if

               if (field_cursor % real4dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real4dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real4dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real4dField % constituentNames(j), real3d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real4dField % constituentNames(j), real3d_temp, io_err)
                  end if

               else
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real4dField % fieldName, real4d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % real4dField % isVarArray) then
               deallocate(real3d_temp)
            else
               deallocate(real4d_temp)
            end if
            call iamas_dmpar_get_time(sub_timer_stop_time)
            left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)

         else if (field_cursor % field_type == FIELD_5D_REAL) then

            call iamas_dmpar_get_time(sub_timer_start_time)
            if (field_cursor % real5dField % isVarArray) then
               ncons = size(field_cursor % real5dField % constituentNames)
               allocate(real4d_temp(field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real5d_temp(field_cursor % real5dField % dimSizes(1), &
                                    field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_5dreal_ptr => field_cursor % real5dField
                  i = 1
                  do while (associated(field_5dreal_ptr))
                     if (trim(field_5dreal_ptr % dimNames(5)) == 'nCells') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, trim(field_5dreal_ptr % dimNames(5)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, trim(field_5dreal_ptr % dimNames(5)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, trim(field_5dreal_ptr % dimNames(5)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, field_5dreal_ptr % dimNames(5), ownedSize)
                     end if

                     if (field_cursor % real5dField % isVarArray) then
                        real4d_temp(:,:,:,i:i+ownedSize-1) = field_5dreal_ptr % array(j,:,:,:,1:ownedSize)
                     else
                        real5d_temp(:,:,:,:,i:i+ownedSize-1) = field_5dreal_ptr % array(:,:,:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_5dreal_ptr => field_5dreal_ptr % next
                     else
                        nullify(field_5dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real5dField % isVarArray) then
                     real4d_temp(:,:,:,:) = field_cursor % real5dField % array(j,:,:,:,:)
                  else
                     real5d_temp(:,:,:,:,:) = field_cursor % real5dField % array(:,:,:,:,:)
                  end if
               end if

               if (field_cursor % real5dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real5dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real5dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
                  end if

                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
               else
                  call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % real5dField % fieldName, real5d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % real5dField % isVarArray) then
               deallocate(real4d_temp)
            else
               deallocate(real5d_temp)
            end if
            call iamas_dmpar_get_time(sub_timer_stop_time)
            left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)

         else if (field_cursor % field_type == FIELD_0D_CHAR) then

!call iamas_log_write('Writing out field '//trim(field_cursor % char0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('Copying field from first block')
!call iamas_log_write('MGD calling iAMAS_io_put_var_mul now...')
            call iamas_dmpar_get_time(sub_timer_start_time)
            call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % char0dField % fieldName, field_cursor % char0dField % scalar, io_err)
            call iamas_dmpar_get_time(sub_timer_stop_time)
            sub_timer_sum = sub_timer_sum + (sub_timer_stop_time - sub_timer_start_time)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_CHAR) then

!call iamas_log_write('Writing out field '//trim(field_cursor % char1dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('Copying field from first block')
!call iamas_log_write('MGD calling iAMAS_io_put_var_mul now...')
            call iamas_dmpar_get_time(sub_timer_start_time)
            call iAMAS_io_put_var_mul(stream % fileHandle, field_cursor % char1dField % fieldName, field_cursor % char1dField % array, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            call iamas_dmpar_get_time(sub_timer_stop_time)
            left_timer_sum = left_timer_sum + (sub_timer_stop_time - sub_timer_start_time)


         end if
         field_cursor => field_cursor % next
      end do


      !
      ! Sync all fields with disk
      !
      call iamas_dmpar_get_time(sync_timer_start_time)
      call iAMAS_io_sync(stream % fileHandle)
      call iamas_dmpar_get_time(sync_timer_stop_time)

   end subroutine iAMAS_writeStream_mul

   subroutine iAMAS_writeStream(stream, frame, arrayFieldList, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      integer, intent(in) :: frame
      type (iAMAS_arrayField_list_type), intent(in):: arrayFieldList
      integer, intent(out), optional :: ierr

      character (len=StrKIND) :: dimensionSuffix
      integer :: io_err
      integer :: i, j, k, iVarArray
      integer :: ncons
      integer, pointer :: ownedSize
      type (field0dInteger), pointer :: field_0dint_ptr
      type (field1dInteger), pointer :: field_1dint_ptr
      type (field2dInteger), pointer :: field_2dint_ptr
      type (field3dInteger), pointer :: field_3dint_ptr
      type (field0dReal), pointer :: field_0dreal_ptr
      type (field1dReal), pointer :: field_1dreal_ptr
      type (field2dReal), pointer :: field_2dreal_ptr
      type (field3dReal), pointer :: field_3dreal_ptr
      type (field4dReal), pointer :: field_4dreal_ptr
      type (field5dReal), pointer :: field_5dreal_ptr
      type (field0dChar), pointer :: field_0dchar_ptr
      type (field1dChar), pointer :: field_1dchar_ptr
      type (field_list_type), pointer :: field_cursor
      integer                            :: int0d_temp
      integer, dimension(:),     pointer :: int1d_temp
      integer, dimension(:,:),   pointer :: int2d_temp
      integer, dimension(:,:,:), pointer :: int3d_temp
      real (kind=RKIND)                                :: real0d_temp
      real (kind=RKIND), dimension(:),         pointer :: real1d_temp
      real (kind=RKIND), dimension(:,:),       pointer :: real2d_temp
      real (kind=RKIND), dimension(:,:,:),     pointer :: real3d_temp
      real (kind=RKIND), dimension(:,:,:,:),   pointer :: real4d_temp
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: real5d_temp

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      !
      ! Set time frame to write
      !
      call iAMAS_io_set_frame(stream % fileHandle, frame, io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR) then
         if (present(ierr)) ierr = iAMAS_IO_ERR
         return
      end if

      !
      ! Check whether we will clobber any records
      !
      if (iAMAS_io_would_clobber_records(stream % fileHandle, io_err)) then
         if (.not. stream % clobberRecords) then
            if (present(ierr)) ierr = iAMAS_STREAM_CLOBBER_RECORD
            return
         else
            call iamas_log_write('iAMAS I/O: Overwriting existing record $i in output file '// &
                  trim(stream % filename), iAMAS_LOG_WARN, intArgs=(/frame/))
         end if
      end if

      if ( stream % blockWrite ) then
         dimensionSuffix = ''
      else
         dimensionSuffix = 'Solve'
      end if

      !
      ! Loop over fields in the stream
      !
      field_cursor => stream % fieldList
      do while (associated(field_cursor))

         if (field_cursor % field_type == FIELD_0D_INT) then

!call iamas_log_write('Writing out field '//trim(field_cursor % int0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('Copying field from first block')
            int0d_temp = field_cursor % int0dField % scalar

!call iamas_log_write('MGD calling iAMAS_io_put_var now...')
            call iAMAS_io_put_var(stream % fileHandle, field_cursor % int0dField % fieldName, int0d_temp, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_INT) then

            if (field_cursor % int1dField % isVarArray) then
               ncons = size(field_cursor % int1dField % constituentNames)
            else
               ncons = 1
               allocate(int1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_1dint_ptr => field_cursor % int1dField
                  i = 1
                  do while (associated(field_1dint_ptr))
                     if (trim(field_1dint_ptr % dimNames(1)) == 'nCells') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, trim(field_1dint_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, trim(field_1dint_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, trim(field_1dint_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_1dint_ptr % block % dimensions, field_1dint_ptr % dimNames(1), ownedSize)
                     end if

                     if (field_cursor % int1dField % isVarArray) then
! I suspect we will never hit this code, as it doesn't make sense, really
                        int0d_temp = field_1dint_ptr % array(j)
                     else
                        int1d_temp(i:i+ownedSize-1) = field_1dint_ptr % array(1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_1dint_ptr => field_1dint_ptr % next
                     else
                        nullify(field_1dint_ptr)
                     end if
                  end do
               else
                  if (field_cursor % int1dField % isVarArray) then
                     int0d_temp = field_cursor % int1dField % array(j)
                  else
                     int1d_temp(:) = field_cursor % int1dField % array(:)
                  end if
               end if

               if (field_cursor % int1dField % isVarArray) then
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % int1dField % constituentNames(j), int0d_temp, io_err)
               else
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % int1dField % fieldName, int1d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (.not. field_cursor % int1dField % isVarArray) then
               deallocate(int1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_INT) then

            if (field_cursor % int2dField % isVarArray) then
               ncons = size(field_cursor % int2dField % constituentNames)
               allocate(int1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int2d_temp(field_cursor % int2dField % dimSizes(1), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_2dint_ptr => field_cursor % int2dField
                  i = 1
                  do while (associated(field_2dint_ptr))
                     if (trim(field_2dint_ptr % dimNames(2)) == 'nCells') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, trim(field_2dint_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, trim(field_2dint_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, trim(field_2dint_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_2dint_ptr % block % dimensions, field_2dint_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % int2dField % isVarArray) then
                        int1d_temp(i:i+ownedSize-1) = field_2dint_ptr % array(j,1:ownedSize)
                     else
                        int2d_temp(:,i:i+ownedSize-1) = field_2dint_ptr % array(:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_2dint_ptr => field_2dint_ptr % next
                     else
                        nullify(field_2dint_ptr)
                     end if
                  end do
               else
                  if (field_cursor % int2dField % isVarArray) then
                     int1d_temp(:) = field_cursor % int2dField % array(j,:)
                  else
                     int2d_temp(:,:) = field_cursor % int2dField % array(:,:)
                  end if
               end if

               if (field_cursor % int2dField % isVarArray) then
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % int2dField % constituentNames(j), int1d_temp, io_err)
               else
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % int2dField % fieldName, int2d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % int2dField % isVarArray) then
               deallocate(int1d_temp)
            else
               deallocate(int2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_INT) then

            if (field_cursor % int3dField % isVarArray) then
               ncons = size(field_cursor % int3dField % constituentNames)
               allocate(int2d_temp(field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int3d_temp(field_cursor % int3dField % dimSizes(1), &
                                   field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_3dint_ptr => field_cursor % int3dField
                  i = 1
                  do while (associated(field_3dint_ptr))
                     if (trim(field_3dint_ptr % dimNames(3)) == 'nCells') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, trim(field_3dint_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, trim(field_3dint_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, trim(field_3dint_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_3dint_ptr % block % dimensions, field_3dint_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % int3dField % isVarArray) then
                        int2d_temp(:,i:i+ownedSize-1) = field_3dint_ptr % array(j,:,1:ownedSize)
                     else
                        int3d_temp(:,:,i:i+ownedSize-1) = field_3dint_ptr % array(:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_3dint_ptr => field_3dint_ptr % next
                     else
                        nullify(field_3dint_ptr)
                     end if
                  end do
               else
                  if (field_cursor % int3dField % isVarArray) then
                     int2d_temp(:,:) = field_cursor % int3dField % array(j,:,:)
                  else
                     int3d_temp(:,:,:) = field_cursor % int3dField % array(:,:,:)
                  end if
               end if

               if (field_cursor % int3dField % isVarArray) then
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % int3dField % constituentNames(j), int2d_temp, io_err)
               else
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % int3dField % fieldName, int3d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % int3dField % isVarArray) then
               deallocate(int2d_temp)
            else
               deallocate(int3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_0D_REAL) then

!call iamas_log_write('Writing out field '//trim(field_cursor % real0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('Copying field from first block')
            real0d_temp = field_cursor % real0dField % scalar

!call iamas_log_write('MGD calling iAMAS_io_put_var now...')
            call iAMAS_io_put_var(stream % fileHandle, field_cursor % real0dField % fieldName, real0d_temp, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_REAL) then

            if (field_cursor % real1dField % isVarArray) then
               ncons = size(field_cursor % real1dField % constituentNames)
            else
               ncons = 1
               allocate(real1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_1dreal_ptr => field_cursor % real1dField
                  i = 1
                  do while (associated(field_1dreal_ptr))
                     if (trim(field_1dreal_ptr % dimNames(1)) == 'nCells') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, trim(field_1dreal_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, trim(field_1dreal_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, trim(field_1dreal_ptr % dimNames(1)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_1dreal_ptr % block % dimensions, field_1dreal_ptr % dimNames(1), ownedSize)
                     end if

                     if (field_cursor % real1dField % isVarArray) then
! I suspect we will never hit this code, as it doesn't make sense, really
                        real0d_temp = field_1dreal_ptr % array(j)
                     else
                        real1d_temp(i:i+ownedSize-1) = field_1dreal_ptr % array(1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_1dreal_ptr => field_1dreal_ptr % next
                     else
                        nullify(field_1dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real1dField % isVarArray) then
                     real0d_temp = field_cursor % real1dField % array(j)
                  else
                     real1d_temp(:) = field_cursor % real1dField % array(:)
                  end if
               end if

               if (field_cursor % real1dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real1dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real1dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var(stream % fileHandle, field_cursor % real1dField % constituentNames(j), real0d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var(stream % fileHandle, field_cursor % real1dField % constituentNames(j), real0d_temp, io_err)
                  end if

               else
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % real1dField % fieldName, real1d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (.not. field_cursor % real1dField % isVarArray) then
               deallocate(real1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_REAL) then

            if (field_cursor % real2dField % isVarArray) then
               ncons = size(field_cursor % real2dField % constituentNames)
               allocate(real1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real2d_temp(field_cursor % real2dField % dimSizes(1), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_2dreal_ptr => field_cursor % real2dField
                  i = 1
                  do while (associated(field_2dreal_ptr))
                     if (trim(field_2dreal_ptr % dimNames(2)) == 'nCells') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, trim(field_2dreal_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, trim(field_2dreal_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, trim(field_2dreal_ptr % dimNames(2)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_2dreal_ptr % block % dimensions, field_2dreal_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % real2dField % isVarArray) then
                        real1d_temp(i:i+ownedSize-1) = field_2dreal_ptr % array(j,1:ownedSize)
                     else
                        real2d_temp(:,i:i+ownedSize-1) = field_2dreal_ptr % array(:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_2dreal_ptr => field_2dreal_ptr % next
                     else
                        nullify(field_2dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real2dField % isVarArray) then
                     real1d_temp(:) = field_cursor % real2dField % array(j,:)
                  else
                     real2d_temp(:,:) = field_cursor % real2dField % array(:,:)
                  end if
               end if

               if (field_cursor % real2dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real2dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real2dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var(stream % fileHandle, field_cursor % real2dField % constituentNames(j), real1d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var(stream % fileHandle, field_cursor % real2dField % constituentNames(j), real1d_temp, io_err)
                  end if

               else
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % real2dField % fieldName, real2d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % real2dField % isVarArray) then
               deallocate(real1d_temp)
            else
               deallocate(real2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_REAL) then

            if (field_cursor % real3dField % isVarArray) then
               ncons = size(field_cursor % real3dField % constituentNames)
               allocate(real2d_temp(field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real3d_temp(field_cursor % real3dField % dimSizes(1), &
                                    field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_3dreal_ptr => field_cursor % real3dField
                  i = 1
                  do while (associated(field_3dreal_ptr))
                     if (trim(field_3dreal_ptr % dimNames(3)) == 'nCells') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, trim(field_3dreal_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, trim(field_3dreal_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, trim(field_3dreal_ptr % dimNames(3)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_3dreal_ptr % block % dimensions, field_3dreal_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % real3dField % isVarArray) then
                        real2d_temp(:,i:i+ownedSize-1) = field_3dreal_ptr % array(j,:,1:ownedSize)
                     else
                        real3d_temp(:,:,i:i+ownedSize-1) = field_3dreal_ptr % array(:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_3dreal_ptr => field_3dreal_ptr % next
                     else
                        nullify(field_3dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real3dField % isVarArray) then
                     real2d_temp(:,:) = field_cursor % real3dField % array(j,:,:)
                  else
                     real3d_temp(:,:,:) = field_cursor % real3dField % array(:,:,:)
                  end if
               end if

               if (field_cursor % real3dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real3dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real3dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var(stream % fileHandle, field_cursor % real3dField % constituentNames(j), real2d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var(stream % fileHandle, field_cursor % real3dField % constituentNames(j), real2d_temp, io_err)
                  end if

               else
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % real3dField % fieldName, real3d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % real3dField % isVarArray) then
               deallocate(real2d_temp)
            else
               deallocate(real3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_4D_REAL) then

            if (field_cursor % real4dField % isVarArray) then
               ncons = size(field_cursor % real4dField % constituentNames)
               allocate(real3d_temp(field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real4d_temp(field_cursor % real4dField % dimSizes(1), &
                                    field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_4dreal_ptr => field_cursor % real4dField
                  i = 1
                  do while (associated(field_4dreal_ptr))
                     if (trim(field_4dreal_ptr % dimNames(4)) == 'nCells') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, trim(field_4dreal_ptr % dimNames(4)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, trim(field_4dreal_ptr % dimNames(4)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, trim(field_4dreal_ptr % dimNames(4)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_4dreal_ptr % block % dimensions, field_4dreal_ptr % dimNames(4), ownedSize)
                     end if

                     if (field_cursor % real4dField % isVarArray) then
                        real3d_temp(:,:,i:i+ownedSize-1) = field_4dreal_ptr % array(j,:,:,1:ownedSize)
                     else
                        real4d_temp(:,:,:,i:i+ownedSize-1) = field_4dreal_ptr % array(:,:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_4dreal_ptr => field_4dreal_ptr % next
                     else
                        nullify(field_4dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real4dField % isVarArray) then
                     real3d_temp(:,:,:) = field_cursor % real4dField % array(j,:,:,:)
                  else
                     real4d_temp(:,:,:,:) = field_cursor % real4dField % array(:,:,:,:)
                  end if
               end if

               if (field_cursor % real4dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real4dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real4dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var(stream % fileHandle, field_cursor % real4dField % constituentNames(j), real3d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var(stream % fileHandle, field_cursor % real4dField % constituentNames(j), real3d_temp, io_err)
                  end if

               else
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % real4dField % fieldName, real4d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % real4dField % isVarArray) then
               deallocate(real3d_temp)
            else
               deallocate(real4d_temp)
            end if

         else if (field_cursor % field_type == FIELD_5D_REAL) then

            if (field_cursor % real5dField % isVarArray) then
               ncons = size(field_cursor % real5dField % constituentNames)
               allocate(real4d_temp(field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real5d_temp(field_cursor % real5dField % dimSizes(1), &
                                    field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_5dreal_ptr => field_cursor % real5dField
                  i = 1
                  do while (associated(field_5dreal_ptr))
                     if (trim(field_5dreal_ptr % dimNames(5)) == 'nCells') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, trim(field_5dreal_ptr % dimNames(5)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nEdges') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, trim(field_5dreal_ptr % dimNames(5)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nVertices') then
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, trim(field_5dreal_ptr % dimNames(5)) &
                                                     // trim(dimensionSuffix), ownedSize)
                     else
                        call iamas_pool_get_dimension(field_5dreal_ptr % block % dimensions, field_5dreal_ptr % dimNames(5), ownedSize)
                     end if

                     if (field_cursor % real5dField % isVarArray) then
                        real4d_temp(:,:,:,i:i+ownedSize-1) = field_5dreal_ptr % array(j,:,:,:,1:ownedSize)
                     else
                        real5d_temp(:,:,:,:,i:i+ownedSize-1) = field_5dreal_ptr % array(:,:,:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     if ( .not. stream % blockWrite ) then
                        field_5dreal_ptr => field_5dreal_ptr % next
                     else
                        nullify(field_5dreal_ptr)
                     end if
                  end do
               else
                  if (field_cursor % real5dField % isVarArray) then
                     real4d_temp(:,:,:,:) = field_cursor % real5dField % array(j,:,:,:,:)
                  else
                     real5d_temp(:,:,:,:,:) = field_cursor % real5dField % array(:,:,:,:,:)
                  end if
               end if

               if (field_cursor % real5dField % isVarArray) then
                  ! jwfeng: Check that if this varArray has specific vars in stream_list
                  iVarArray = 0
                  do k = 1, maxVarArrays
                     if ( trim(field_cursor % real5dField % fieldname) .eq. trim(arrayFieldList%arrayFields(k)) ) then
                         iVarArray = k
                     end if
                  end do
                  if (iVarArray .ne. 0) then
                     do k = 1, maxVarsInArray
                        if ( trim(field_cursor % real5dField % constituentNames(j)) .eq. trim(arrayFieldList%varsInArray(k,iVarArray)) ) then
                           call iAMAS_io_put_var(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
                        end if
                     end do
                  else
                     call iAMAS_io_put_var(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
                  end if

                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
               else
                  call iAMAS_io_put_var(stream % fileHandle, field_cursor % real5dField % fieldName, real5d_temp, io_err)
               end if
               call iAMAS_io_err_mesg(io_err, .false.)
               if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR
            end do

            if (field_cursor % real5dField % isVarArray) then
               deallocate(real4d_temp)
            else
               deallocate(real5d_temp)
            end if

         else if (field_cursor % field_type == FIELD_0D_CHAR) then

!call iamas_log_write('Writing out field '//trim(field_cursor % char0dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('Copying field from first block')
!call iamas_log_write('MGD calling iAMAS_io_put_var now...')
            call iAMAS_io_put_var(stream % fileHandle, field_cursor % char0dField % fieldName, field_cursor % char0dField % scalar, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_CHAR) then

!call iamas_log_write('Writing out field '//trim(field_cursor % char1dField % fieldName))
!call iamas_log_write('   > is the field decomposed? $l', logicArgs=(/field_cursor % isDecomposed/))
!call iamas_log_write('   > outer dimension size $i', intArgs=(/field_cursor % totalDimSize/))

!call iamas_log_write('Copying field from first block')
!call iamas_log_write('MGD calling iAMAS_io_put_var now...')
            call iAMAS_io_put_var(stream % fileHandle, field_cursor % char1dField % fieldName, field_cursor % char1dField % array, io_err)
            call iAMAS_io_err_mesg(io_err, .false.)
            if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR


         end if
         field_cursor => field_cursor % next
      end do


      !
      ! Sync all fields with disk
      !
      call iAMAS_io_sync(stream % fileHandle)

   end subroutine iAMAS_writeStream


   subroutine iAMAS_readStreamAtt_0dInteger(stream, attName, attValue, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      integer, intent(out) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      call iAMAS_io_get_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_readStreamAtt_0dInteger


   subroutine iAMAS_readStreamAtt_1dInteger(stream, attName, attValue, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      integer, dimension(:), pointer :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      call iAMAS_io_get_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_readStreamAtt_1dInteger


   subroutine iAMAS_readStreamAtt_0dReal(stream, attName, attValue, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(out) :: attValue
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: local_precision

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      call iAMAS_io_get_att(stream % fileHandle, attName, attValue, precision=local_precision, ierr=io_err)
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_readStreamAtt_0dReal


   subroutine iAMAS_readStreamAtt_1dReal(stream, attName, attValue, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), pointer :: attValue
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: local_precision

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      call iAMAS_io_get_att(stream % fileHandle, attName, attValue, precision=local_precision, ierr=io_err)
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_readStreamAtt_1dReal


   subroutine iAMAS_readStreamAtt_text(stream, attName, attValue, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      character (len=*), intent(out) :: attValue
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      call iAMAS_io_get_att(stream % fileHandle, attName, attValue, ierr=io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_readStreamAtt_text


   subroutine iAMAS_writeStreamAtt_0dInteger(stream, attName, attValue, syncVal, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      integer, intent(in) :: attValue
      logical, intent(in), optional :: syncVal
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      call iAMAS_io_put_att(stream % fileHandle, attName, attValue, syncVal=syncVal, ierr=io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_writeStreamAtt_0dInteger


   subroutine iAMAS_writeStreamAtt_1dInteger(stream, attName, attValue, syncVal, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      integer, dimension(:), intent(in) :: attValue
      logical, intent(in), optional :: syncVal
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      call iAMAS_io_put_att(stream % fileHandle, attName, attValue, syncVal=syncVal, ierr=io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_writeStreamAtt_1dInteger


   subroutine iAMAS_writeStreamAtt_0dReal(stream, attName, attValue, syncVal, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      real (kind=RKIND), intent(in) :: attValue
      logical, intent(in), optional :: syncVal
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: local_precision

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      call iAMAS_io_put_att(stream % fileHandle, attName, attValue, syncVal=syncVal, precision=local_precision, ierr=io_err)
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_writeStreamAtt_0dReal


   subroutine iAMAS_writeStreamAtt_1dReal(stream, attName, attValue, syncVal, precision, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      real (kind=RKIND), dimension(:), intent(in) :: attValue
      logical, intent(in), optional :: syncVal
      integer, intent(in), optional :: precision
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: local_precision

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      if (present(precision)) then
         local_precision = precision
      else
         local_precision = stream % defaultPrecision
      end if

      call iAMAS_io_put_att(stream % fileHandle, attName, attValue, syncVal=syncVal, precision=local_precision, ierr=io_err)
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_writeStreamAtt_1dReal


   subroutine iAMAS_writeStreamAtt_text(stream, attName, attValue, syncVal, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      character (len=*), intent(in) :: attName
      character (len=*), intent(in) :: attValue
      logical, intent(in), optional :: syncVal
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         return
      end if

      call iAMAS_io_put_att(stream % fileHandle, attName, attValue, syncVal=syncVal, ierr=io_err) 
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) ierr = iAMAS_IO_ERR

   end subroutine iAMAS_writeStreamAtt_text


   subroutine iAMAS_closeStream(stream, ierr)

      implicit none

      type (iAMAS_Stream_type), intent(inout) :: stream
      integer, intent(out), optional :: ierr

      integer :: io_err
      type (field_list_type), pointer :: field_cursor

      if (present(ierr)) ierr = iAMAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = iAMAS_STREAM_NOT_INITIALIZED
         if (.not. stream % valid)then
             if (present(ierr))ierr = iAMAS_IO_NOERR
         end if
         return
      end if

      call iAMAS_io_close(stream % fileHandle, io_err)
      call iAMAS_io_err_mesg(io_err, .false.)
      if (io_err /= iAMAS_IO_NOERR .and. present(ierr)) then
          ierr = iAMAS_IO_ERR
      end if
      
      if (.not. stream % valid)then
          if (present(ierr))ierr = iAMAS_IO_NOERR
          return 
      end if

!call iamas_log_write('Deallocating global attribute list')
      call iamas_deallocate_attlist(stream % attList)

!call iamas_log_write('Deallocating field list')
      field_cursor => stream % fieldList
      do while (associated(field_cursor))
         if (associated(field_cursor % isAvailable)) then
            deallocate(field_cursor % isAvailable)
!call iamas_log_write('Deallocating isAvailable array')
         end if
         stream % fieldList => stream % fieldList % next
         deallocate(field_cursor)
         field_cursor => stream % fieldList
      end do

      stream % isInitialized = .false.

   end subroutine iAMAS_closeStream


   subroutine mergeArrays(array1, array2)

      implicit none

      integer, dimension(:), pointer :: array1
      integer, dimension(:), intent(in) :: array2

      integer :: n1, n2
      integer, dimension(:), pointer :: newArray

      n1 = size(array1)
      n2 = size(array2)

      allocate(newArray(n1+n2))

      newArray(1:n1) = array1(:)
      newArray(n1+1:n1+n2) = array2(:)

      deallocate(array1)
      array1 => newArray

   end subroutine mergeArrays


   subroutine put_get_field_atts(handle, ioDirection, fieldname, attList)

      implicit none

      type (iAMAS_IO_Handle_type), intent(inout) :: handle
      integer, intent(in) :: ioDirection
      character (len=*), intent(in) :: fieldname
      type (att_list_type), pointer :: attList

      type (att_list_type), pointer :: att_cursor

      if (.not. associated(attList)) return

      att_cursor => attList
      if (ioDirection == iAMAS_IO_WRITE) then
         do while (associated(att_cursor))
            select case (att_cursor % attType)
               case (iAMAS_ATT_INT)
                  call iAMAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueInt, fieldname)
               case (iAMAS_ATT_INTA)
                  call iAMAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueIntA, fieldname)
               case (iAMAS_ATT_REAL)
                  call iAMAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueReal, fieldname)
               case (iAMAS_ATT_REALA)
                  call iAMAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueRealA, fieldname)
               case (iAMAS_ATT_TEXT)
                  call iAMAS_io_put_att(handle, trim(att_cursor % attName), att_cursor % attValueText, fieldname)
            end select 
            att_cursor => att_cursor % next
         end do
      else
         do while (associated(att_cursor))
            select case (att_cursor % attType)
               case (iAMAS_ATT_INT)
                  call iAMAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueInt, fieldname)
               case (iAMAS_ATT_INTA)
                  call iAMAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueIntA, fieldname)
               case (iAMAS_ATT_REAL)
                  call iAMAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueReal, fieldname)
               case (iAMAS_ATT_REALA)
                  call iAMAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueRealA, fieldname)
               case (iAMAS_ATT_TEXT)
                  call iAMAS_io_get_att(handle, trim(att_cursor % attName), att_cursor % attValueText, fieldname)
            end select 
            att_cursor => att_cursor % next
         end do
      end if

   end subroutine put_get_field_atts

end module iamas_io_streams
