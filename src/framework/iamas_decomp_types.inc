   integer, parameter :: iAMAS_DECOMP_NOERR = 1000, &
                         iAMAS_DECOMP_ERROR = 1001

   abstract interface
      function iamas_decomp_function(block, manager, globalDimSize, numBlocks, ownedIndices) result(iErr)
         import block_type
         import iamas_streamManager_type

         type (block_type), intent(in) :: block
         type (iamas_streamManager_type), intent(inout) :: manager
         integer, intent(in) :: globalDimSize
         integer, intent(in) :: numBlocks
         integer, dimension(:), pointer :: ownedIndices
         integer :: iErr
      end function
   end interface

   type iamas_decomp_list
      integer :: nameLen
      character (len=StrKIND) :: decompName
      procedure (iamas_decomp_function), pointer, nopass :: decompFunc => null()
      type (iamas_decomp_list), pointer :: next => null()
   end type iamas_decomp_list
