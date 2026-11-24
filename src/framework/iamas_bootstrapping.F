! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://git.ustc.edu.cn/ustc_lacar/development/iamas/-/blob/iAMAS/LICENSE
!
module iamas_bootstrapping

   use iamas_derived_types
   use iamas_field_routines
   use iamas_pool_routines
   use iamas_dmpar
   use iamas_block_decomp
   use iamas_block_creator
   use iamas_sort
   use iamas_timekeeping
   use iamas_io_streams
   use iamas_log


   integer :: readCellStart, readCellEnd, nReadCells
   integer :: readEdgeStart, readEdgeEnd, nReadEdges
   integer :: readVertexStart, readVertexEnd, nReadVertices


   contains


   !***********************************************************************
   !
   !  routine iamas_bootstrap_framework_phase1
   !
   !> \brief   Obtains mesh partition, builds halos, and allocates blocks.
   !> \author  Michael Duda, Doug Jacobsen
   !> \date    03/04/2015
   !> \details
   !>  This routine is responsible for reading basic decomposed field dimensions,
   !>  determining a block decomposition of the mesh, building halos for fields,
   !>  and, finally, allocating blocklist for a domain.
   !>
   !>  Phase1 does not currently allocate fields, or define any other
   !>  dimensions aside from the ones listed below.
   !>
   !>  Dimensions required to be present in the grid file:
   !>      nCells
   !>      nEdges
   !>      nVertices
   !>      vertexDegree
   !>      maxEdges
   !>
   !>  Fields required to be present in the grid file:
   !>      indexToCellID
   !>      indexToEdgeID
   !>      indexToVertexID
   !>      {x,y,z}Cell
   !>      {x,y,z}Edge
   !>      {x,y,z}Vertex
   !>      nEdgesOnCell
   !>      cellsOnCell
   !>      edgesOnCell
   !>      verticesOnCell
   !>      cellsOnEdge
   !>      cellsOnVertex
   !>
   !>  Attributes required to be present in the grid file:
   !>      on_a_sphere
   !>      sphere_radius
   !>      parent_id
   !>      mesh_spec
   !>      is_periodic
   !>      x_period
   !>      y_period
   !>      ***** these are needed by iamas_block_creator_finalize_block_init()
   !>            so they can be set in the mesh pool and queried by, e.g.,
   !>            iamas_initialize_vectors()
   !
   !-----------------------------------------------------------------------
   subroutine iamas_bootstrap_framework_phase1(domain, mesh_filename, mesh_iotype) !{{{

      use mpi
      implicit none

      type (domain_type), pointer :: domain
      character(len=*), intent(in) :: mesh_filename
      integer, intent(in) :: mesh_iotype

      type (block_type), pointer :: readingBlock

      integer :: nCells, nEdges, maxEdges, nVertices, vertexDegree
      integer :: nCells_mul, nEdges_mul, nVertices_mul, maxEdges_mul, vertexDegree_mul

      integer :: ierr
      type (iAMAS_IO_Handle_type) :: inputHandle, inputHandle_mul

      character (len=StrKIND) :: c_on_a_sphere, c_is_periodic, c_parent_id, c_file_id, c_mesh_spec
      real (kind=RKIND) :: r_sphere_radius, r_x_period, r_y_period

      type (field1dInteger), pointer :: indexToCellIDField
      type (field1dInteger), pointer :: indexToEdgeIDField
      type (field1dInteger), pointer :: indexToVertexIDField
      type (field1dInteger), pointer :: nEdgesOnCellField
      type (field2dInteger), pointer :: cellsOnCellField
      type (field2dInteger), pointer :: edgesOnCellField
      type (field2dInteger), pointer :: verticesOnCellField
      type (field2dInteger), pointer :: cellsOnEdgeField
      type (field2dInteger), pointer :: cellsOnVertexField

      type (field1dReal), pointer :: xCellField,   yCellField,   zCellField
      type (field1dReal), pointer :: xEdgeField,   yEdgeField,   zEdgeField
      type (field1dReal), pointer :: xVertexField, yVertexField, zVertexField

      type (field1dInteger), pointer :: nCellsSolveField
      type (field1dInteger), pointer :: nVerticesSolveField
      type (field1dInteger), pointer :: nEdgesSolveField

      type (field1DInteger), pointer :: indexToCellID_Block
      type (field1DInteger), pointer :: nEdgesOnCell_Block
      type (field2DInteger), pointer :: cellsOnCell_Block
      type (field2DInteger), pointer :: verticesOnCell_Block
      type (field2DInteger), pointer :: edgesOnCell_Block

      type (field1DInteger), pointer :: indexToVertexID_Block
      type (field2DInteger), pointer :: cellsOnVertex_Block

      type (field1DInteger), pointer :: indexToEdgeID_Block
      type (field2DInteger), pointer :: cellsOnEdge_Block

      integer, dimension(:), pointer :: local_cell_list
      integer, dimension(:), pointer :: block_id, block_start, block_count
      type (graph) :: partial_global_graph_info

      type (iAMAS_Time_type) :: startTime
      character(len=StrKIND) :: timeStamp, restartTimeStamp
      character(len=StrKIND) :: filename

      integer, pointer :: config_num_halos, config_number_of_blocks
      logical, pointer :: config_explicit_proc_decomp
      character (len=StrKIND), pointer :: config_block_decomp_file_prefix, config_proc_decomp_file_prefix
      integer :: nHalos
      real (kind=R8KIND) :: timer_start_time, timer_stop_time!jwfeng

      integer :: mpi_ierr
      logical :: io_use_mulfiles, is_init_atm


      call iamas_pool_get_config(domain % configs, 'config_num_halos', config_num_halos)
      call iamas_pool_get_config(domain % configs, 'config_number_of_blocks', config_number_of_blocks)
      call iamas_pool_get_config(domain % configs, 'config_explicit_proc_decomp', config_explicit_proc_decomp)
      call iamas_pool_get_config(domain % configs, 'config_block_decomp_file_prefix', config_block_decomp_file_prefix)
      call iamas_pool_get_config(domain % configs, 'config_proc_decomp_file_prefix', config_proc_decomp_file_prefix)

      nHalos = config_num_halos
      io_use_mulfiles = domain % ioContext % io_use_mulfiles
      is_init_atm = domain % ioContext % is_init_atm
      if( io_use_mulfiles .and. (.not. is_init_atm)) then
         inputHandle_mul = iAMAS_io_open_mul(trim(mesh_filename), iAMAS_IO_READ, mesh_iotype, domain % ioContext, ierr=ierr)
      else
         inputHandle = iAMAS_io_open(trim(mesh_filename), iAMAS_IO_READ, mesh_iotype, domain % ioContext, ierr=ierr)
      end if

      if (ierr /= iAMAS_IO_NOERR) then
         call iamas_log_write('Could not open input file '''//trim(mesh_filename)//''' to read mesh fields', iAMAS_LOG_CRIT)
      else
         call iamas_log_write('Bootstrapping framework with mesh fields from input file '''//trim(mesh_filename)//'''')
      end if


      !
      ! Read global number of cells/edges/vertices
      !
      call iamas_log_write('--- Reading global number of cells/edges/vertices')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      if(io_use_mulfiles .and. (.not. is_init_atm)) then
         call iamas_io_inq_dim(inputHandle_mul, 'nCells', nCells_mul, ierr)
         call iamas_io_inq_dim(inputHandle_mul, 'nEdges', nEdges_mul, ierr)
         call iamas_io_inq_dim(inputHandle_mul, 'nVertices', nVertices_mul, ierr)
         call iamas_io_inq_dim(inputHandle_mul, 'maxEdges', maxEdges_mul, ierr)
         call iamas_io_inq_dim(inputHandle_mul, 'vertexDegree', vertexDegree_mul, ierr)
         call MPI_Allreduce(nCells_mul, nCells, 1, MPI_INTEGER, MPI_SUM, domain%dminfo%io_comm, mpi_ierr)
         call MPI_Allreduce(nEdges_mul, nEdges, 1, MPI_INTEGER, MPI_SUM, domain%dminfo%io_comm, mpi_ierr)
         call MPI_Allreduce(nVertices_mul, nVertices, 1, MPI_INTEGER, MPI_SUM, domain%dminfo%io_comm, mpi_ierr)

         maxEdges = maxEdges_mul
         vertexDegree = vertexDegree_mul
      else
         call iamas_io_inq_dim(inputHandle, 'nCells', nCells, ierr)
         call iamas_io_inq_dim(inputHandle, 'nEdges', nEdges, ierr)
         call iamas_io_inq_dim(inputHandle, 'nVertices', nVertices, ierr)
         call iamas_io_inq_dim(inputHandle, 'maxEdges', maxEdges, ierr)
         call iamas_io_inq_dim(inputHandle, 'vertexDegree', vertexDegree, ierr)
      end if
      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

      !
      ! Determine the range of cells/edges/vertices that a processor will initially read
      !   from the input file
      !
      call iamas_log_write('---Getting index of cells/edges/vertices')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      if(io_use_mulfiles .and. (.not. is_init_atm)) then
         call iamas_dmpar_get_index_range_mul(domain % dminfo, 1, nCells_mul, readCellStart, readCellEnd)
      else
          call iamas_dmpar_get_index_range(domain % dminfo, 1, nCells, readCellStart, readCellEnd)
      end if
      nReadCells = readCellEnd - readCellStart + 1

      if(io_use_mulfiles .and. (.not. is_init_atm)) then
         call iamas_dmpar_get_index_range_mul(domain % dminfo, 1, nEdges_mul, readEdgeStart, readEdgeEnd)
      else
         call iamas_dmpar_get_index_range(domain % dminfo, 1, nEdges, readEdgeStart, readEdgeEnd)
      end if
      nReadEdges = readEdgeEnd - readEdgeStart + 1

      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iamas_dmpar_get_index_range_mul(domain % dminfo, 1, nVertices_mul, readVertexStart, readVertexEnd)
      else
          call iamas_dmpar_get_index_range(domain % dminfo, 1, nVertices, readVertexStart, readVertexEnd)
      end if
      nReadVertices = readVertexEnd - readVertexStart + 1

      allocate(readingBlock)
      readingBlock % domain => domain
      readingBlock % blockID = domain % dminfo % my_proc_id
      readingBlock % localBlockID = 0

      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

      !
      ! Allocate and read fields that we will need in order to ultimately work out
      !   which cells/edges/vertices are owned by each block, and which are ghost
      !

      call iamas_log_write('--- Setting up cell block fields')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iamas_io_setup_cell_block_fields_mul(inputHandle_mul, nreadCells, readCellStart, readingBlock, maxEdges, &
                                           indexTocellIDField, xCellField, yCellField, zCellField, nEdgesOnCellField, &
                                           cellsOnCellField, edgesOnCellField, verticesOnCellField, nHalos)
      else
          call iamas_io_setup_cell_block_fields(inputHandle, nreadCells, readCellStart, readingBlock, maxEdges, &
                                           indexTocellIDField, xCellField, yCellField, zCellField, nEdgesOnCellField, &
                                           cellsOnCellField, edgesOnCellField, verticesOnCellField, nHalos)
      end if
      call iamas_dmpar_get_time(timer_stop_time)!jwfeng

      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng
      call iamas_log_write('--- Setting up edge block fields')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iamas_io_setup_edge_block_fields_mul(inputHandle_mul, nReadEdges, readEdgeStart, readingBlock, indexToEdgeIDField, &
                                            xEdgeField, yEdgeField, zEdgeField, cellsOnEdgeField, nHalos)
      else
          call iamas_io_setup_edge_block_fields(inputHandle, nReadEdges, readEdgeStart, readingBlock, indexToEdgeIDField, & 
              xEdgeField, yEdgeField, zEdgeField, cellsOnEdgeField, nHalos)
      end if
      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

      call iamas_log_write('--- Setting up vertex block fields')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iamas_io_setup_vertex_block_fields_mul(inputHandle_mul, nReadVertices, readVertexStart, readingBlock, vertexDegree, &
                                                 indexToVertexIDField, xVertexField, yVertexField, zVertexField, cellsOnVertexField, &
                                                 nHalos)
      else
          call iamas_io_setup_vertex_block_fields(inputHandle, nReadVertices, readVertexStart, readingBlock, vertexDegree, &
                                                 indexToVertexIDField, xVertexField, yVertexField, zVertexField, cellsOnVertexField, &
                                                 nHalos)
      end if
      
      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng
      !
      ! Set up a graph derived data type describing the connectivity for the cells
      !   that were read by this process
      ! A partial description is passed to the block decomp module by each process,
      !   and the block decomp module returns with a list of global cell indices
      !   that belong to the block on this process
      !
      partial_global_graph_info % nVertices = nReadCells
      partial_global_graph_info % nVerticesTotal = nCells
      partial_global_graph_info % maxDegree = maxEdges
      partial_global_graph_info % ghostStart = nVertices+1
      allocate(partial_global_graph_info % vertexID(nReadCells))
      allocate(partial_global_graph_info % nAdjacent(nReadCells))
      allocate(partial_global_graph_info % adjacencyList(maxEdges, nReadCells))


      partial_global_graph_info % vertexID(:) = indexToCellIDField % array(:)
      partial_global_graph_info % nAdjacent(:) = nEdgesOnCellField % array(:)
      partial_global_graph_info % adjacencyList(:,:) = cellsOnCellField % array(:,:)

      ! TODO: Ensure (by renaming or exchanging) that initial cell range on each proc is contiguous
      !       This situation may occur when reading a restart file with cells/edges/vertices written
      !       in a scrambled order

      !
      ! Determine which cells are owned by this process. Currently, this routine just reads a graph.info.part.*
      !    file, but in principle it could call some online, distributed mesh partitioning library.
      !
      call iamas_log_write('--- Reading graph.info.part.* file')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      call iamas_block_decomp_cells_for_proc(domain % dminfo, partial_global_graph_info, local_cell_list, block_id, block_start, &
                                            block_count, config_number_of_blocks, config_explicit_proc_decomp, &
                                            config_block_decomp_file_prefix, config_proc_decomp_file_prefix)

      deallocate(partial_global_graph_info % vertexID)
      deallocate(partial_global_graph_info % nAdjacent)
      deallocate(partial_global_graph_info % adjacencyList)

      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

      !
      ! Transfer owned cells from local_cell_list into the indexToCellID_Block field (which may contain multiple blocks)
      !
      call iamas_log_write('--- Building 0 and 1 halo of cell/edge/vertex block fields')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      call iamas_block_creator_setup_blocks_and_0halo_cells(nHalos, domain, indexToCellID_Block, local_cell_list, block_id, &
                                                           block_start, block_count)

      !
      ! The *Field arguments contain fields read by iamas_io_setup_*_block_fields on the simple, contiguous decompositions
      !    that were set up by iamas_dmpar_get_index_range. These next calls communicate the *Field mesh fields from these
      !    simple decompositions to the native iAMAS decomposition in the *_Block fields.
      !
      call iamas_block_creator_build_0halo_cell_fields(nHalos, indexToCellIDField, nEdgesOnCellField, cellsOnCellField, &
                                                      verticesOnCellField, edgesOnCellField, indexToCellID_Block, &
                                                      nEdgesOnCell_Block, cellsOnCell_Block, verticesOnCell_Block, &
                                                      edgesOnCell_Block)

      call iamas_block_creator_build_0_and_1halo_edge_fields(nHalos, indexToEdgeIDField, cellsOnEdgeField, indexToCellID_Block, &
                                                            nEdgesOnCell_Block, edgesOnCell_Block, indexToEdgeID_Block, &
                                                            cellsOnEdge_Block, nEdgesSolveField)
      call iamas_block_creator_build_0_and_1halo_edge_fields(nHalos, indexToVertexIDField, cellsOnVertexField, indexToCellID_Block, &
                                                            nEdgesOnCell_Block, verticesOnCell_Block, indexToVertexID_Block, &
                                                            cellsOnVertex_Block, nVerticesSolveField)
      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

      !
      ! The following three calls are responsible for taking the *_Block fields that upon entry
      !    to the routines are defined only for owned elements and adding nHalos layers on to them so that
      !    upon exit they contain all owned and halo elements. A side effect of this is that the exchange
      !    lists for these fields are valid (and can be used, e.g., by iamas_block_creator_finalize_block_phase1).
      !
      call iamas_log_write('--- Building halos of cell/edge/vertex block fields')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      call iamas_block_creator_build_cell_halos(nHalos, indexToCellID_Block, nEdgesOnCell_Block, cellsOnCell_Block, &
                                               verticesOnCell_Block, edgesOnCell_Block, nCellsSolveField)

      call iamas_block_creator_build_edge_halos(nHalos, indexToCellID_Block, nEdgesOnCell_Block, nCellsSolveField, &
                                               edgesOnCell_Block, indexToEdgeID_Block, cellsOnEdge_Block, nEdgesSolveField)
      call iamas_block_creator_build_edge_halos(nHalos, indexToCellID_Block, nEdgesOnCell_Block, nCellsSolveField, &
                                               verticesOnCell_Block, indexToVertexID_Block, cellsOnVertex_Block, &
                                               nVerticesSolveField)
      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

      !
      ! Before we can allocate blocks, we need the attributes on_a_sphere and sphere_radius so
      !    they can be propagated as configs to subpools
      !
      call iamas_log_write('--- Getting attributes ... ')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iAMAS_io_get_att(inputHandle_mul, 'on_a_sphere', c_on_a_sphere, ierr=ierr)
      else
          call iAMAS_io_get_att(inputHandle, 'on_a_sphere', c_on_a_sphere, ierr=ierr)
      end if
     
      if (ierr /= iAMAS_IO_NOERR) then
         call iamas_log_write('Attribute on_a_sphere not found in '//trim(mesh_filename), iAMAS_LOG_WARN)
         call iamas_log_write('   Setting on_a_sphere to ''YES''', iAMAS_LOG_WARN)
         domain % on_a_sphere = .true.
      else
         if (index(c_on_a_sphere, 'YES') /= 0) then
            domain % on_a_sphere = .true.
         else
            domain % on_a_sphere = .false.
         end if
      end if

      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iAMAS_io_get_att(inputHandle_mul, 'sphere_radius', r_sphere_radius, ierr=ierr)
      else
          call iAMAS_io_get_att(inputHandle, 'sphere_radius', r_sphere_radius, ierr=ierr)
      end if
      if (ierr /= iAMAS_IO_NOERR) then
         call iamas_log_write('Attribute sphere_radius not found in '//trim(mesh_filename), iAMAS_LOG_WARN)
         call iamas_log_write('   Setting sphere_radius to 1.0', iAMAS_LOG_WARN)
         domain % sphere_radius = 1.0
      else
         domain % sphere_radius = r_sphere_radius
      end if

      !
      ! We need to read in the file_id, parent_id, and mesh spec attributes as well, to
      ! set them appropriately in the output files.
      !
      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iAMAS_io_get_att(inputHandle_mul, 'file_id', c_file_id, ierr=ierr)
      else 
          call iAMAS_io_get_att(inputHandle, 'file_id', c_file_id, ierr=ierr)
      end if
      if ( ierr /= iAMAS_IO_NOERR) then
         call iamas_log_write('Attribute file_id not found in '//trim(mesh_filename), iAMAS_LOG_WARN)
         call iamas_log_write('  Using '''' for the previous file_id.', iAMAS_LOG_WARN)
         c_file_id = ''
      end if

      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iAMAS_io_get_att(inputHandle_mul, 'parent_id', c_parent_id, ierr=ierr)
      else
          call iAMAS_io_get_att(inputHandle, 'parent_id', c_parent_id, ierr=ierr)
      end if
      if ( ierr /= iAMAS_IO_NOERR) then
         call iamas_log_write('Attribute parent_id not found in '//trim(mesh_filename), iAMAS_LOG_WARN)
         call iamas_log_write('  Setting parent_id to ''''', iAMAS_LOG_WARN)
         c_parent_id = ''
      end if

      if ( trim(c_file_id) == '' ) then
         domain % parent_id = c_parent_id
      else
         domain % parent_id = trim(c_file_id) // NEW_LINE('A') // trim(c_parent_id)
      end if

      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iAMAS_io_get_att(inputHandle_mul, 'mesh_spec', c_mesh_spec, ierr=ierr)
      else
          call iAMAS_io_get_att(inputHandle, 'mesh_spec', c_mesh_spec, ierr=ierr)
      end if
      if ( ierr /= iAMAS_IO_NOERR) then
         call iamas_log_write('Attribute mesh_spec not found in '//trim(mesh_filename), iAMAS_LOG_WARN)
         call iamas_log_write('  Setting mesh_spec to ''0.0''', iAMAS_LOG_WARN)
         domain % mesh_spec = '0.0'
      else
        domain % mesh_spec = c_mesh_spec
      end if

      !
      ! Get attributes based on planar mesh
      !
      if ( .not. domain % on_a_sphere ) then
         if(io_use_mulfiles .and. (.not. is_init_atm)) then
              call iAMAS_io_get_att(inputHandle_mul, 'is_periodic', c_is_periodic, ierr=ierr)
         else
              call iAMAS_io_get_att(inputHandle, 'is_periodic', c_is_periodic, ierr=ierr)
         end if
         if ( ierr /= iAMAS_IO_NOERR ) then
            call iamas_log_write('is_periodic attribute not found in '//trim(mesh_filename), iAMAS_LOG_WARN)
            call iamas_log_write('   Setting is_periodic to .false.', iAMAS_LOG_WARN)
            domain % is_periodic = .false.
         else
            if (index(c_is_periodic, 'YES') /= 0) then
               domain % is_periodic = .true.
            else
               domain % is_periodic = .false.
            end if
         end if
      else
        domain % is_periodic = .false.
      end if

      if ( domain % is_periodic ) then
         if(io_use_mulfiles .and. (.not. is_init_atm)) then
             call iAMAS_io_get_att(inputHandle_mul, 'x_period', r_x_period, ierr=ierr)
         else
             call iAMAS_io_get_att(inputHandle, 'x_period', r_x_period, ierr=ierr)
         end if
         if ( ierr /= iAMAS_IO_NOERR ) then
            call iamas_log_write('Required attribute ''x_period'' missing from periodic input mesh. Exiting...')
         else
            domain % x_period = r_x_period
         end if

         if(io_use_mulfiles .and. (.not. is_init_atm)) then
             call iAMAS_io_get_att(inputHandle_mul, 'y_period', r_y_period, ierr=ierr)
         else
             call iAMAS_io_get_att(inputHandle, 'y_period', r_y_period, ierr=ierr)
         end if
         if ( ierr /= iAMAS_IO_NOERR ) then
            call iamas_log_write('Required attribute ''y_period'' missing from periodic input mesh. Exiting...')
         else
            domain % y_period = r_y_period
         end if
      else
         domain % x_period = 0
         domain % y_period = 0
      end if
      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng


     ! Allocate blocks, and copy indexTo arrays into blocks
     call iamas_log_write('--- Creating final blocks-phase1 ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     call iamas_block_creator_finalize_block_phase1(nHalos, domain % blocklist, &
             nCellsSolveField, nEdgesSolveField, nVerticesSolveField, vertexDegree, maxEdges, &
             indexToCellID_Block, indexToEdgeID_Block, indexToVertexID_Block)
      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

      call iamas_log_write('--- Closing/Destroying/Deallocating files halos and fields ...')!jwfeng
      call iamas_dmpar_get_time(timer_start_time) !jwfeng
      if(io_use_mulfiles .and. (.not. is_init_atm)) then
          call iAMAS_io_close(inputHandle_mul, ierr)
      else
          call iAMAS_io_close(inputHandle, ierr)
      end if 

      call iamas_dmpar_destroy_mulithalo_exchange_list(indexToCellIDField % sendList)
      call iamas_dmpar_destroy_mulithalo_exchange_list(indexToCellIDField % recvList)
      call iamas_dmpar_destroy_mulithalo_exchange_list(indexToCellIDField % copyList)

      call iamas_dmpar_destroy_mulithalo_exchange_list(indexToEdgeIDField % sendList)
      call iamas_dmpar_destroy_mulithalo_exchange_list(indexToEdgeIDField % recvList)
      call iamas_dmpar_destroy_mulithalo_exchange_list(indexToEdgeIDField % copyList)

      call iamas_dmpar_destroy_mulithalo_exchange_list(indexToVertexIDField % sendList)
      call iamas_dmpar_destroy_mulithalo_exchange_list(indexToVertexIDField % recvList)
      call iamas_dmpar_destroy_mulithalo_exchange_list(indexToVertexIDField % copyList)

      call iamas_deallocate_field(indexToCellIDField)
      call iamas_deallocate_field(indexToEdgeIDField)
      call iamas_deallocate_field(indexToVertexIDField)

      call iamas_deallocate_field(nEdgesOnCellField)
      call iamas_deallocate_field(cellsOnCellField)
      call iamas_deallocate_field(edgesOnCellField)
      call iamas_deallocate_field(verticesOnCellField)
      call iamas_deallocate_field(cellsOnEdgeField)
      call iamas_deallocate_field(cellsOnVertexField)

      call iamas_deallocate_field(indexToCellID_Block)
      call iamas_deallocate_field(nEdgesOnCell_Block)
      call iamas_deallocate_field(cellsOnCell_Block)
      call iamas_deallocate_field(verticesOnCell_Block)
      call iamas_deallocate_field(edgesOnCell_Block)
      call iamas_deallocate_field(indexToVertexID_Block)
      call iamas_deallocate_field(cellsOnVertex_Block)
      call iamas_deallocate_field(indexToEdgeID_Block)
      call iamas_deallocate_field(cellsOnEdge_Block)

      call iamas_deallocate_field(nCellsSolveField)
      call iamas_deallocate_field(nVerticesSolveField)
      call iamas_deallocate_field(nEdgesSolveField)

      deallocate(local_cell_list)
      deallocate(block_id)
      deallocate(block_start)
      deallocate(block_count)
      deallocate(readingBlock)
      call iamas_dmpar_get_time(timer_stop_time)!jwfeng
      call iamas_log_write('--- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

   end subroutine iamas_bootstrap_framework_phase1 !}}}


   !***********************************************************************
   !
   !  routine iamas_bootstrap_framework_phase2
   !
   !> \brief   Defines dimensions and allocates fields
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    03/04/2015
   !> \details
   !>  This routine assumes blocks have already been setup correctly (some
   !>  fields have the correct indices, and exchange lists are created).
   !>  It will finish the setup of blocks by defining all remaining dimensions,
   !>  and allocating all fields and structs.
   !
   !-----------------------------------------------------------------------
   subroutine iamas_bootstrap_framework_phase2(domain) !{{{

      use iamas_stream_manager
      use iamas_stream_list
      use mpi

      implicit none

      type (domain_type), pointer :: domain

      type (iamas_pool_type), pointer :: readableDimensions
      type (iamas_pool_type), pointer :: streamDimensions

      character (len=StrKIND) :: streamName, streamFilenameTemplate, streamFilenameInterval, streamReferenceTime, streamFilename
      character (len=StrKIND) :: fieldName
      integer :: streamDirection
      logical :: streamActive
      logical :: fieldActive
      integer :: ioType
      integer :: mpi_ierr

      type (iamas_pool_iterator_type) :: poolItr

      type (iAMAS_IO_Handle_type) :: inputHandle, inputHandle_mul

      character (len=StrKIND), dimension(:), pointer :: dimNames
      integer, pointer :: dimValue
      integer :: tempDim, tempDim_mul
      integer :: i, err_level, err_local

      logical :: io_use_mulfiles, is_init_atm

      io_use_mulfiles = domain % ioContext % io_use_mulfiles 
      is_init_atm = domain % ioContext % is_init_atm
      call iamas_pool_create_pool(readableDimensions)

      err_level = iamas_pool_get_error_level()
      call iamas_pool_set_error_level(iAMAS_POOL_SILENT)

      call iamas_log_write(' ')
      call iamas_log_write(' ')

      ! Reading dimensions from streams
      call iamas_log_write('Reading dimensions from input streams ...')
      call iamas_stream_mgr_begin_iteration(domain % streamManager)
      do while ( iamas_stream_mgr_get_next_stream(domain % streamManager, streamID = streamName, directionProperty = streamDirection, &
                                                 activeProperty = streamActive) )

         if ( streamActive .and. ( streamDirection == iAMAS_STREAM_INPUT .or. streamDirection == iAMAS_STREAM_INPUT_OUTPUT ) ) then

            call iamas_stream_mgr_begin_iteration(domain % streamManager, streamID=streamName)

            ! Build stream dimension pool from the list of fields
            call iamas_pool_create_pool(streamDimensions)

            do while ( iamas_stream_mgr_get_next_field(domain % streamManager, streamName, fieldName, isActive=fieldActive) )

               if (fieldActive) then
                  call get_dimlist_for_field(domain % blocklist % allFields, fieldName, dimNames)

                  do i=1,size(dimNames)
                     call iamas_pool_get_dimension(streamDimensions, dimNames(i), dimValue)
                     if ( .not. associated(dimValue) ) then
                        call iamas_pool_add_dimension(streamDimensions, dimNames(i), iAMAS_MISSING_DIM)
                     end if
                  end do
                  deallocate(dimNames)
               end if

            end do

            ! Determine stream filename
            call iamas_get_stream_filename(domain % streamManager, streamID = streamName, filename = streamFilename, ierr = err_local)

            ! Determine stream io_type
            call iAMAS_stream_mgr_get_property(domain % streamManager, streamName, &
                                              iAMAS_STREAM_PROPERTY_IOTYPE, ioType, ierr = err_local)

            ! Try to open file
            if(trim(streamName) == 'input' .or. trim(streamName) == 'restart' .or. trim(streamName) == 'chem_input')then
               if (io_use_mulfiles .and. (.not. is_init_atm)) then
                  inputHandle_mul = iAMAS_io_open_mul(trim(streamFilename), iAMAS_IO_READ, ioType, domain % ioContext, ierr = err_local)
               else
                  inputHandle = iAMAS_io_open(trim(streamFilename), iAMAS_IO_READ, ioType, domain % ioContext, ierr = err_local)
               end if
            else
               inputHandle = iAMAS_io_open(trim(streamFilename), iAMAS_IO_READ, ioType, domain % ioContext, ierr = err_local)
            end if

            ! If to determine if file was opened or not.
            if ( err_local == iAMAS_IO_NOERR ) then

               call iamas_log_write(' ')
               call iamas_log_write('----- reading dimensions from stream '''//trim(streamName)//''' using file ' &
                                    //trim(streamFilename) )

               ! Iterate over list of dimensions we determined we need from the above loop
               call iamas_pool_begin_iteration(streamDimensions)
               do while ( iamas_pool_get_next_member(streamDimensions, poolItr) )
                  if ( poolItr % memberType == iAMAS_POOL_DIMENSION ) then
                     ! Try to read the dimension
                     if(trim(streamName) == 'input' .or. trim(streamName) == 'restart' .or. trim(streamName) == 'chem_input') then
                        if (io_use_mulfiles .and. (.not. is_init_atm)) then
                           call iamas_io_inq_dim(inputHandle_mul, trim(poolItr % memberName), tempDim_mul, ierr = err_local)
                           if (trim(poolItr%memberName) == 'nCells' .or. &
                               trim(poolItr%memberName) == 'nVertices' .or. &
                               trim(poolItr%memberName) == 'nEdges') then
                               call MPI_Allreduce(tempDim_mul, tempDim, 1, MPI_INTEGER, MPI_SUM, domain%dminfo%io_comm, mpi_ierr)
                           else
                              tempDim = tempDim_mul
                           end if
                        else
                           call iamas_io_inq_dim(inputHandle, trim(poolItr % memberName), tempDim, ierr = err_local)
                        end if
                     else
                        call iamas_io_inq_dim(inputHandle, trim(poolItr % memberName), tempDim, ierr = err_local)
                     end if

                     ! Check to see if the dimension has already been defined
                     call iamas_pool_get_dimension(readableDimensions, poolItr % memberName, dimValue)

                     ! If to see if dimension was read or not
                     if ( err_local == iAMAS_IO_NOERR ) then
                        call iamas_log_write('       ' // trim(poolItr % memberName) // ' = $i', intArgs=(/tempDim/) )

                        if ( .not. associated(dimValue) ) then
                           call iamas_pool_add_dimension(readableDimensions, poolItr % memberName, tempDim)
                        else if ( dimValue /= tempDim .and. dimValue == iAMAS_MISSING_DIM ) then
                           dimValue = tempDim
                        else if ( dimValue /= tempDim ) then
                           call iamas_log_write('Dimension ' // trim(poolItr % membername) &
                                                // ' was read with an inconsistent value.', iAMAS_LOG_CRIT)
                        end if
                     else
                        call iamas_log_write('       ' // trim(poolItr % memberName) // ' *** not found in stream ***')
                     end if

                  end if
               end do

               ! Close file
               if(trim(streamName) == 'input' .or. trim(streamName) == 'restart' .or. trim(streamName) == 'chem_input') then
                  if (io_use_mulfiles .and. (.not. is_init_atm)) then
                     call iamas_io_close(inputHandle_mul)
                  else
                     call iamas_io_close(inputHandle)
                  end if
               else
                  call iamas_io_close(inputHandle)
               end if
            else
               call iamas_log_write(' ')
               call iamas_log_write('  *** unable to open input file '//trim(streamFilename)//' for stream ''' &
                                    //trim(streamName)//'''')
            end if

            ! Destroy pool that contains list of streams dimensions
            call iamas_pool_destroy_pool(streamDimensions)

         else if ( .not. streamActive .and. ( streamDirection == iAMAS_STREAM_INPUT .or. streamDirection == iAMAS_STREAM_INPUT_OUTPUT ) ) then

            call iamas_log_write(' ')
            call iamas_log_write('----- skipping inactive stream '''//trim(streamName)//'''')
            
         end if

      end do

      call iamas_log_write(' ')
      call iamas_log_write('----- done reading dimensions from input streams -----')
      call iamas_log_write(' ')
      call iamas_log_write(' ')

      call iamas_pool_set_error_level(err_level)

      ! Allocate blocks, and copy indexTo arrays into blocks
      call iamas_block_creator_finalize_block_phase2(domain % streamManager, domain % blocklist, readableDimensions)
      call iamas_link_fields(domain)

      call iamas_pool_destroy_pool(readableDimensions)

   end subroutine iamas_bootstrap_framework_phase2 !}}}

   subroutine iamas_io_setup_cell_block_fields_mul(inputHandle, nReadCells, readCellStart, readingBlock, maxEdges, indexToCellID, &
                                              xCell, yCell, zCell, nEdgesOnCell, cellsOnCell, edgesOnCell, verticesOnCell, nHalos) !{{{

     implicit none

     type (iAMAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadCells
     integer, intent(in) :: readCellStart
     integer, intent(in) :: maxEdges
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToCellID
     type (field1dReal), pointer :: xCell
     type (field1dReal), pointer :: yCell
     type (field1dReal), pointer :: zCell
     type (field1dInteger), pointer :: nEdgesOnCell
     type (field2dInteger), pointer :: cellsOnCell
     type (field2dInteger), pointer :: edgesOnCell
     type (field2dInteger), pointer :: verticesOnCell
     integer, intent(in) :: nHalos

     integer :: i, ierr
     integer, dimension(:), pointer :: readIndices
     real (kind=R8KIND) :: timer_start_time, timer_stop_time!jwfeng
     real (kind=R8KIND) :: sub_timer_start_time, sub_timer_stop_time!jwfeng

     !
     ! Allocate and read fields that we will need in order to ultimately work out
     !   which cells/edges/vertices are owned by each block, and which are ghost
     !

     ! Global cell indices
     call iamas_log_write('      --- Reading and processing indexToCellID ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(indexToCellID)
     allocate(indexToCellID % array(nReadCells))
     allocate(readIndices(nReadCells))
     do i=1,nReadCells
        readIndices(i) = i + readCellStart - 1
     end do
     call iamas_log_write('          --- Calling iAMAS_io_inq_var ...')!jwfeng
     call iAMAS_io_inq_var(inputHandle, 'indexToCellID', ierr=ierr)
     call iamas_log_write('          --- Calling iAMAS_io_set_var_indices_mul ...')!jwfeng
     call iAMAS_io_set_var_indices_mul(inputHandle, 'indexToCellID', readIndices, ierr=ierr)
     call iamas_log_write('          --- Calling iamas_io_get_var ...')!jwfeng
     call iamas_io_get_var_mul(inputHandle, 'indexToCellID', indexToCellID % array, ierr)
     indexToCellID % dimSizes(1) = nReadCells
     indexToCellID % block => readingBlock
     call iamas_log_write('          --- Initializing multihalo exchange list ...')!jwfeng
     call iamas_dmpar_get_time(sub_timer_start_time) !jwfeng
     call iamas_dmpar_init_multihalo_exchange_list(indexToCellID % sendList, nHalos)
     call iamas_dmpar_init_multihalo_exchange_list(indexToCellID % recvList, nHalos)
     call iamas_dmpar_init_multihalo_exchange_list(indexToCellID % copyList, nHalos)
     call iamas_dmpar_get_time(sub_timer_stop_time)!jwfeng
     call iamas_log_write('         --- Cost: $r s', &
          realArgs=(/real(sub_timer_stop_time - sub_timer_start_time, kind = RKIND)/))!jwfeng
     nullify(indexToCellID % next)
     call iamas_dmpar_get_time(timer_stop_time)  !jwfeng
     call iamas_log_write('     --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng


     ! Number of cell/edges/vertices adjacent to each cell
     call iamas_log_write('      --- Reading and processing nEdgesOnCell ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(nEdgesOnCell)
     allocate(nEdgesOnCell % array(nReadCells))
     call iAMAS_io_inq_var(inputHandle, 'nEdgesOnCell', ierr=ierr)
     call iAMAS_io_set_var_indices_mul(inputHandle, 'nEdgesOnCell', readIndices, ierr=ierr)
     call iamas_io_get_var_mul(inputHandle, 'nEdgesOnCell', nEdgesOnCell % array, ierr)
     nEdgesOnCell % dimSizes(1) = nReadCells
     nEdgesOnCell % block => readingBlock
     nEdgesOnCell % sendList => indexToCellID % sendList
     nEdgesOnCell % recvList => indexToCellID % recvList
     nEdgesOnCell % copyList => indexToCellID % copyList
     nullify(nEdgesOnCell % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     ! Global indices of cells adjacent to each cell
     call iamas_log_write('      --- Reading and processing cellsOnCell ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(cellsOnCell)
     allocate(cellsOnCell % array(maxEdges,nReadCells))
     call iAMAS_io_inq_var(inputHandle, 'cellsOnCell', ierr=ierr)
     call iAMAS_io_set_var_indices_mul(inputHandle, 'cellsOnCell', readIndices, ierr=ierr)
     call iamas_io_get_var_mul(inputHandle, 'cellsOnCell', cellsOnCell % array, ierr)
     cellsOnCell % dimSizes(1) = maxEdges
     cellsOnCell % dimSizes(2) = nReadCells
     cellsOnCell % block => readingBlock
     cellsOnCell % sendList => indexToCellID % sendList
     cellsOnCell % recvList => indexToCellID % recvList
     cellsOnCell % copyList => indexToCellID % copyList
     nullify(cellsOnCell % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     ! Global indices of edges adjacent to each cell
     call iamas_log_write('      --- Reading and processing edgesOnCell ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(edgesOnCell)
     allocate(edgesOnCell % array(maxEdges,nReadCells))
     call iAMAS_io_inq_var(inputHandle, 'edgesOnCell', ierr=ierr)
     call iAMAS_io_set_var_indices_mul(inputHandle, 'edgesOnCell', readIndices, ierr=ierr)
     call iamas_io_get_var_mul(inputHandle, 'edgesOnCell', edgesOnCell % array, ierr)
     edgesOnCell % dimSizes(1) = maxEdges
     edgesOnCell % dimSizes(2) = nReadCells
     edgesOnCell % block => readingBlock
     edgesOnCell % sendList => indexToCellID % sendList
     edgesOnCell % recvList => indexToCellID % recvList
     edgesOnCell % copyList => indexToCellID % copyList
     nullify(edgesOnCell % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     ! Global indices of vertices adjacent to each cell
     call iamas_log_write('      --- Reading and processing verticesOnCell ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(verticesOnCell)
     allocate(verticesOnCell % array(maxEdges,nReadCells))
     call iAMAS_io_inq_var(inputHandle, 'verticesOnCell', ierr=ierr)
     call iAMAS_io_set_var_indices_mul(inputHandle, 'verticesOnCell', readIndices, ierr=ierr)
     call iamas_io_get_var_mul(inputHandle, 'verticesOnCell', verticesOnCell % array, ierr)
     verticesOnCell % dimSizes(1) = maxEdges
     verticesOnCell % dimSizes(2) = nReadCells
     verticesOnCell % block => readingBlock
     verticesOnCell % sendList => indexToCellID % sendList
     verticesOnCell % recvList => indexToCellID % recvList
     verticesOnCell % copyList => indexToCellID % copyList
     nullify(verticesOnCell % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     deallocate(readIndices)

   end subroutine iamas_io_setup_cell_block_fields_mul !}}}

   subroutine iamas_io_setup_cell_block_fields(inputHandle, nReadCells, readCellStart, readingBlock, maxEdges, indexToCellID, &
                                              xCell, yCell, zCell, nEdgesOnCell, cellsOnCell, edgesOnCell, verticesOnCell, nHalos) !{{{

     implicit none

     type (iAMAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadCells
     integer, intent(in) :: readCellStart
     integer, intent(in) :: maxEdges
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToCellID
     type (field1dReal), pointer :: xCell
     type (field1dReal), pointer :: yCell
     type (field1dReal), pointer :: zCell
     type (field1dInteger), pointer :: nEdgesOnCell
     type (field2dInteger), pointer :: cellsOnCell
     type (field2dInteger), pointer :: edgesOnCell
     type (field2dInteger), pointer :: verticesOnCell
     integer, intent(in) :: nHalos

     integer :: i, ierr
     integer, dimension(:), pointer :: readIndices
     real (kind=R8KIND) :: timer_start_time, timer_stop_time!jwfeng
     real (kind=R8KIND) :: sub_timer_start_time, sub_timer_stop_time!jwfeng

     !
     ! Allocate and read fields that we will need in order to ultimately work out
     !   which cells/edges/vertices are owned by each block, and which are ghost
     !

     ! Global cell indices
     call iamas_log_write('      --- Reading and processing indexToCellID ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(indexToCellID)
     allocate(indexToCellID % array(nReadCells))
     allocate(readIndices(nReadCells))
     do i=1,nReadCells
        readIndices(i) = i + readCellStart - 1
     end do
     call iamas_log_write('          --- Calling iAMAS_io_inq_var ...')!jwfeng
     call iAMAS_io_inq_var(inputHandle, 'indexToCellID', ierr=ierr)
     call iamas_log_write('          --- Calling iAMAS_io_set_var_indices ...')!jwfeng
     call iAMAS_io_set_var_indices(inputHandle, 'indexToCellID', readIndices, ierr=ierr)
     call iamas_log_write('          --- Calling iamas_io_get_var ...')!jwfeng
     call iamas_io_get_var(inputHandle, 'indexToCellID', indexToCellID % array, ierr)
     indexToCellID % dimSizes(1) = nReadCells
     indexToCellID % block => readingBlock
     call iamas_log_write('          --- Initializing multihalo exchange list ...')!jwfeng
     call iamas_dmpar_get_time(sub_timer_start_time) !jwfeng
     call iamas_dmpar_init_multihalo_exchange_list(indexToCellID % sendList, nHalos)
     call iamas_dmpar_init_multihalo_exchange_list(indexToCellID % recvList, nHalos)
     call iamas_dmpar_init_multihalo_exchange_list(indexToCellID % copyList, nHalos)
     call iamas_dmpar_get_time(sub_timer_stop_time)!jwfeng
     call iamas_log_write('          --- Cost: $r s', &
          realArgs=(/real(sub_timer_stop_time - sub_timer_start_time, kind = RKIND)/))!jwfeng
     nullify(indexToCellID % next)
     call iamas_dmpar_get_time(timer_stop_time) !jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng


     ! Number of cell/edges/vertices adjacent to each cell
     call iamas_log_write('      --- Reading and processing nEdgesOnCell ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(nEdgesOnCell)
     allocate(nEdgesOnCell % array(nReadCells))
     call iAMAS_io_inq_var(inputHandle, 'nEdgesOnCell', ierr=ierr)
     call iAMAS_io_set_var_indices(inputHandle, 'nEdgesOnCell', readIndices, ierr=ierr)
     call iamas_io_get_var(inputHandle, 'nEdgesOnCell', nEdgesOnCell % array, ierr)
     nEdgesOnCell % dimSizes(1) = nReadCells
     nEdgesOnCell % block => readingBlock
     nEdgesOnCell % sendList => indexToCellID % sendList
     nEdgesOnCell % recvList => indexToCellID % recvList
     nEdgesOnCell % copyList => indexToCellID % copyList
     nullify(nEdgesOnCell % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     ! Global indices of cells adjacent to each cell
     call iamas_log_write('      --- Reading and processing cellsOnCell ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(cellsOnCell)
     allocate(cellsOnCell % array(maxEdges,nReadCells))
     call iAMAS_io_inq_var(inputHandle, 'cellsOnCell', ierr=ierr)
     call iAMAS_io_set_var_indices(inputHandle, 'cellsOnCell', readIndices, ierr=ierr)
     call iamas_io_get_var(inputHandle, 'cellsOnCell', cellsOnCell % array, ierr)
     cellsOnCell % dimSizes(1) = maxEdges
     cellsOnCell % dimSizes(2) = nReadCells
     cellsOnCell % block => readingBlock
     cellsOnCell % sendList => indexToCellID % sendList
     cellsOnCell % recvList => indexToCellID % recvList
     cellsOnCell % copyList => indexToCellID % copyList
     nullify(cellsOnCell % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     ! Global indices of edges adjacent to each cell
     call iamas_log_write('      --- Reading and processing edgesOnCell ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(edgesOnCell)
     allocate(edgesOnCell % array(maxEdges,nReadCells))
     call iAMAS_io_inq_var(inputHandle, 'edgesOnCell', ierr=ierr)
     call iAMAS_io_set_var_indices(inputHandle, 'edgesOnCell', readIndices, ierr=ierr)
     call iamas_io_get_var(inputHandle, 'edgesOnCell', edgesOnCell % array, ierr)
     edgesOnCell % dimSizes(1) = maxEdges
     edgesOnCell % dimSizes(2) = nReadCells
     edgesOnCell % block => readingBlock
     edgesOnCell % sendList => indexToCellID % sendList
     edgesOnCell % recvList => indexToCellID % recvList
     edgesOnCell % copyList => indexToCellID % copyList
     nullify(edgesOnCell % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     ! Global indices of vertices adjacent to each cell
     call iamas_log_write('      --- Reading and processing verticesOnCell ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(verticesOnCell)
     allocate(verticesOnCell % array(maxEdges,nReadCells))
     call iAMAS_io_inq_var(inputHandle, 'verticesOnCell', ierr=ierr)
     call iAMAS_io_set_var_indices(inputHandle, 'verticesOnCell', readIndices, ierr=ierr)
     call iamas_io_get_var(inputHandle, 'verticesOnCell', verticesOnCell % array, ierr)
     verticesOnCell % dimSizes(1) = maxEdges
     verticesOnCell % dimSizes(2) = nReadCells
     verticesOnCell % block => readingBlock
     verticesOnCell % sendList => indexToCellID % sendList
     verticesOnCell % recvList => indexToCellID % recvList
     verticesOnCell % copyList => indexToCellID % copyList
     nullify(verticesOnCell % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     deallocate(readIndices)

   end subroutine iamas_io_setup_cell_block_fields !}}}

   subroutine iamas_io_setup_edge_block_fields_mul(inputHandle, nReadEdges, readEdgeStart, readingBlock, indexToEdgeID, &
                                              xEdge, yEdge, zEdge, cellsOnEdge, nHalos) !{{{

     implicit none

     type (iAMAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadEdges
     integer, intent(in) :: readEdgeStart
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToEdgeID
     type (field1dReal), pointer :: xEdge
     type (field1dReal), pointer :: yEdge
     type (field1dReal), pointer :: zEdge
     type (field2dInteger), pointer :: cellsOnEdge
     integer, intent(in) :: nHalos

     integer :: i, ierr
     integer, dimension(:), pointer :: readIndices
     real (kind=R8KIND) :: timer_start_time,timer_stop_time!jwfeng
     real (kind=R8KIND) :: sub_timer_start_time,sub_timer_stop_time!jwfeng

     !
     ! Allocate and read fields that we will need in order to ultimately work out
     !   which cells/edges/vertices are owned by each block, and which are ghost
     !

     allocate(readIndices(nReadEdges))

     ! Global edge indices
     call iamas_log_write('      --- Reading and processing indexToEdgeID ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(indexToEdgeID)
     allocate(indexToEdgeID % array(nReadEdges))
     !allocate(indexToEdgeID % array(nReadEdges))
     do i=1,nReadEdges
        readIndices(i) = i + readEdgeStart - 1
     end do
     call iAMAS_io_inq_var(inputHandle, 'indexToEdgeID', ierr=ierr)
     call iAMAS_io_set_var_indices_mul(inputHandle, 'indexToEdgeID', readIndices, ierr=ierr)
     call iamas_io_get_var_mul(inputHandle, 'indexToEdgeID', indexToEdgeID % array, ierr)
     indexToEdgeID % dimSizes(1) = nREadEdges
     indexToEdgeID % block => readingBlock
     call iamas_log_write('          --- Initializing multihalo exchange list ...')!jwfeng
     call iamas_dmpar_get_time(sub_timer_start_time) !jwfeng
     call iamas_dmpar_init_multihalo_exchange_list(indexToEdgeID % sendList, nHalos+1)
     call iamas_dmpar_init_multihalo_exchange_list(indexToEdgeID % recvList, nHalos+1)
     call iamas_dmpar_init_multihalo_exchange_list(indexToEdgeID % copyList, nHalos+1)
     call iamas_dmpar_get_time(sub_timer_stop_time) !jwfeng
     call iamas_log_write('          --- Cost: $r s', &
          realArgs=(/real(sub_timer_stop_time - sub_timer_start_time, kind = RKIND)/))!jwfeng
     nullify(indexToEdgeID % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng


     ! Global indices of cells adjacent to each edge
     !    used for determining which edges are owned by a block, where
     !    iEdge is owned iff cellsOnEdge(1,iEdge) is an owned cell
     call iamas_log_write('      --- Reading and processing cellsOnEdge ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(cellsOnEdge)
     allocate(cellsOnEdge % array(2,nReadEdges))
     call iAMAS_io_inq_var(inputHandle, 'cellsOnEdge', ierr=ierr)
     call iAMAS_io_set_var_indices_mul(inputHandle, 'cellsOnEdge', readIndices, ierr=ierr)
     call iamas_io_get_var_mul(inputHandle, 'cellsOnEdge', cellsOnEdge % array, ierr)
     cellsOnEdge % dimSizes(1) = 2
     cellsOnEdge % dimSizes(2) = nReadEdges
     cellsOnEdge % block => readingBlock
     cellsOnEdge % sendList => indexToEdgeID % sendList
     cellsOnEdge % recvList => indexToEdgeID % recvList
     cellsOnEdge % copyList => indexToEdgeID % copyList
     nullify(cellsOnEdge % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     deallocate(readIndices)

   end subroutine iamas_io_setup_edge_block_fields_mul !}}}


   subroutine iamas_io_setup_edge_block_fields(inputHandle, nReadEdges, readEdgeStart, readingBlock, indexToEdgeID, &
                                              xEdge, yEdge, zEdge, cellsOnEdge, nHalos) !{{{

     implicit none

     type (iAMAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadEdges
     integer, intent(in) :: readEdgeStart
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToEdgeID
     type (field1dReal), pointer :: xEdge
     type (field1dReal), pointer :: yEdge
     type (field1dReal), pointer :: zEdge
     type (field2dInteger), pointer :: cellsOnEdge
     integer, intent(in) :: nHalos

     integer :: i, ierr
     integer, dimension(:), pointer :: readIndices
     real (kind=R8KIND) :: timer_start_time,timer_stop_time!jwfeng
     real (kind=R8KIND) :: sub_timer_start_time,sub_timer_stop_time!jwfeng

     !
     ! Allocate and read fields that we will need in order to ultimately work out
     !   which cells/edges/vertices are owned by each block, and which are ghost
     !

     allocate(readIndices(nReadEdges))

     ! Global edge indices
     call iamas_log_write('      --- Reading and processing indexToEdgeID ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(indexToEdgeID)
     allocate(indexToEdgeID % array(nReadEdges))
     !allocate(indexToEdgeID % array(nReadEdges))
     do i=1,nReadEdges
        readIndices(i) = i + readEdgeStart - 1
     end do
     call iAMAS_io_inq_var(inputHandle, 'indexToEdgeID', ierr=ierr)
     call iAMAS_io_set_var_indices(inputHandle, 'indexToEdgeID', readIndices, ierr=ierr)
     call iamas_io_get_var(inputHandle, 'indexToEdgeID', indexToEdgeID % array, ierr)
     indexToEdgeID % dimSizes(1) = nREadEdges
     indexToEdgeID % block => readingBlock
     call iamas_log_write('          --- Initializing multihalo exchange list ...')!jwfeng
     call iamas_dmpar_get_time(sub_timer_start_time) !jwfeng
     call iamas_dmpar_init_multihalo_exchange_list(indexToEdgeID % sendList, nHalos+1)
     call iamas_dmpar_init_multihalo_exchange_list(indexToEdgeID % recvList, nHalos+1)
     call iamas_dmpar_init_multihalo_exchange_list(indexToEdgeID % copyList, nHalos+1)
     call iamas_dmpar_get_time(sub_timer_stop_time) !jwfeng
     call iamas_log_write('          --- Cost: $r s', &
          realArgs=(/real(sub_timer_stop_time - sub_timer_start_time, kind = RKIND)/))!jwfeng
     nullify(indexToEdgeID % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng


     ! Global indices of cells adjacent to each edge
     !    used for determining which edges are owned by a block, where
     !    iEdge is owned iff cellsOnEdge(1,iEdge) is an owned cell
     call iamas_log_write('      --- Reading and processing cellsOnEdge ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(cellsOnEdge)
     allocate(cellsOnEdge % array(2,nReadEdges))
     call iAMAS_io_inq_var(inputHandle, 'cellsOnEdge', ierr=ierr)
     call iAMAS_io_set_var_indices(inputHandle, 'cellsOnEdge', readIndices, ierr=ierr)
     call iamas_io_get_var(inputHandle, 'cellsOnEdge', cellsOnEdge % array, ierr)
     cellsOnEdge % dimSizes(1) = 2
     cellsOnEdge % dimSizes(2) = nReadEdges
     cellsOnEdge % block => readingBlock
     cellsOnEdge % sendList => indexToEdgeID % sendList
     cellsOnEdge % recvList => indexToEdgeID % recvList
     cellsOnEdge % copyList => indexToEdgeID % copyList
     nullify(cellsOnEdge % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     deallocate(readIndices)

   end subroutine iamas_io_setup_edge_block_fields !}}}

   subroutine iamas_io_setup_vertex_block_fields_mul(inputHandle, nReadVertices, readVertexStart, readingBlock, vertexDegree, &
                                                indexToVertexID, xVertex, yVertex, zVertex, cellsOnVertex, nHalos) !{{{

     implicit none

     type (iAMAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadVertices
     integer, intent(in) :: readVertexStart
     integer, intent(in) :: vertexDegree
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToVertexID
     type (field1dReal), pointer :: xVertex
     type (field1dReal), pointer :: yVertex
     type (field1dReal), pointer :: zVertex
     type (field2dInteger), pointer :: cellsOnVertex
     integer, intent(in) :: nHalos

     integer :: i, ierr
     integer, dimension(:), pointer :: readIndices
     real (kind=R8KIND) :: timer_start_time, timer_stop_time!jwfeng
     real (kind=R8KIND) :: sub_timer_start_time, sub_timer_stop_time!jwfeng

     ! Global vertex indices
     call iamas_log_write('      --- Reading and processing indexToVertexID ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(indexToVertexID)
     allocate(indexToVertexID % array(nReadVertices))
     allocate(readIndices(nReadVertices))
     do i=1,nReadVertices
        readIndices(i) = i + readVertexStart - 1
     end do
     call iAMAS_io_inq_var(inputHandle, 'indexToVertexID', ierr=ierr)
     call iAMAS_io_set_var_indices_mul(inputHandle, 'indexToVertexID', readIndices, ierr=ierr)
     call iamas_io_get_var_mul(inputHandle, 'indexToVertexID', indexToVertexID % array, ierr)
     indexToVertexID % dimSizes(1) = nReadVertices
     indexToVertexID % block => readingBlock
     call iamas_log_write('          --- Initializing multihalo exchange list ...')!jwfeng
     call iamas_dmpar_get_time(sub_timer_start_time) !jwfeng
     call iamas_dmpar_init_multihalo_exchange_list(indexToVertexID % sendList, nHalos+1)
     call iamas_dmpar_init_multihalo_exchange_list(indexToVertexID % recvList, nHalos+1)
     call iamas_dmpar_init_multihalo_exchange_list(indexToVertexID % copyList, nHalos+1)
     call iamas_dmpar_get_time(sub_timer_stop_time)!jwfeng
     call iamas_log_write('          --- Cost: $r s', &
          realArgs=(/real(sub_timer_stop_time - sub_timer_start_time, kind = RKIND)/))!jwfeng
     nullify(indexToVertexID % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng


     ! Global indices of cells adjacent to each vertex
     !    used for determining which vertices are owned by a block, where
     !    iVtx is owned iff cellsOnVertex(1,iVtx) is an owned cell
     call iamas_log_write('      --- Reading and processing cellsOnVertex ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(cellsOnVertex)
     allocate(cellsOnVertex % array(vertexDegree,nReadVertices))
     call iAMAS_io_inq_var(inputHandle, 'cellsOnVertex', ierr=ierr)
     call iAMAS_io_set_var_indices_mul(inputHandle, 'cellsOnVertex', readIndices, ierr=ierr)
     call iamas_io_get_var_mul(inputHandle, 'cellsOnVertex', cellsOnVertex % array, ierr)
     cellsOnVertex % dimSizes(1) = vertexDegree
     cellsOnVertex % dimSizes(2) = nReadVertices
     cellsOnVertex % block => readingBlock
     cellsOnVertex % sendList => indexToVertexID % sendList
     cellsOnVertex % recvList => indexToVertexID % recvList
     cellsOnVertex % copyList => indexToVertexID % copyList
     nullify(cellsOnVertex % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     deallocate(readIndices)

   end subroutine iamas_io_setup_vertex_block_fields_mul !}}}

   subroutine iamas_io_setup_vertex_block_fields(inputHandle, nReadVertices, readVertexStart, readingBlock, vertexDegree, &
                                                indexToVertexID, xVertex, yVertex, zVertex, cellsOnVertex, nHalos) !{{{

     implicit none

     type (iAMAS_IO_Handle_type) :: inputHandle
     integer, intent(in) :: nReadVertices
     integer, intent(in) :: readVertexStart
     integer, intent(in) :: vertexDegree
     type (block_type), pointer :: readingBlock
     type (field1dInteger), pointer :: indexToVertexID
     type (field1dReal), pointer :: xVertex
     type (field1dReal), pointer :: yVertex
     type (field1dReal), pointer :: zVertex
     type (field2dInteger), pointer :: cellsOnVertex
     integer, intent(in) :: nHalos

     integer :: i, ierr
     integer, dimension(:), pointer :: readIndices
     real (kind=R8KIND) :: timer_start_time, timer_stop_time!jwfeng
     real (kind=R8KIND) :: sub_timer_start_time, sub_timer_stop_time!jwfeng

     ! Global vertex indices
     call iamas_log_write('      --- Reading and processing indexToVertexID ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(indexToVertexID)
     allocate(indexToVertexID % array(nReadVertices))
     allocate(readIndices(nReadVertices))
     do i=1,nReadVertices
        readIndices(i) = i + readVertexStart - 1
     end do
     call iAMAS_io_inq_var(inputHandle, 'indexToVertexID', ierr=ierr)
     call iAMAS_io_set_var_indices(inputHandle, 'indexToVertexID', readIndices, ierr=ierr)
     call iamas_io_get_var(inputHandle, 'indexToVertexID', indexToVertexID % array, ierr)
     indexToVertexID % dimSizes(1) = nReadVertices
     indexToVertexID % block => readingBlock
     call iamas_log_write('          --- Initializing multihalo exchange list ...')!jwfeng
     call iamas_dmpar_get_time(sub_timer_start_time) !jwfeng
     call iamas_dmpar_init_multihalo_exchange_list(indexToVertexID % sendList, nHalos+1)
     call iamas_dmpar_init_multihalo_exchange_list(indexToVertexID % recvList, nHalos+1)
     call iamas_dmpar_init_multihalo_exchange_list(indexToVertexID % copyList, nHalos+1)
     call iamas_dmpar_get_time(sub_timer_stop_time)!jwfeng
     call iamas_log_write('          --- Cost: $r s', &
          realArgs=(/real(sub_timer_stop_time - sub_timer_start_time, kind = RKIND)/))!jwfeng
     nullify(indexToVertexID % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng


     ! Global indices of cells adjacent to each vertex
     !    used for determining which vertices are owned by a block, where
     !    iVtx is owned iff cellsOnVertex(1,iVtx) is an owned cell
     call iamas_log_write('      --- Reading and processing cellsOnVertex ...')!jwfeng
     call iamas_dmpar_get_time(timer_start_time) !jwfeng
     allocate(cellsOnVertex)
     allocate(cellsOnVertex % array(vertexDegree,nReadVertices))
     call iAMAS_io_inq_var(inputHandle, 'cellsOnVertex', ierr=ierr)
     call iAMAS_io_set_var_indices(inputHandle, 'cellsOnVertex', readIndices, ierr=ierr)
     call iamas_io_get_var(inputHandle, 'cellsOnVertex', cellsOnVertex % array, ierr)
     cellsOnVertex % dimSizes(1) = vertexDegree
     cellsOnVertex % dimSizes(2) = nReadVertices
     cellsOnVertex % block => readingBlock
     cellsOnVertex % sendList => indexToVertexID % sendList
     cellsOnVertex % recvList => indexToVertexID % recvList
     cellsOnVertex % copyList => indexToVertexID % copyList
     nullify(cellsOnVertex % next)
     call iamas_dmpar_get_time(timer_stop_time)!jwfeng
     call iamas_log_write('      --- Cost: $r s', &
          realArgs=(/real(timer_stop_time - timer_start_time, kind = RKIND)/))!jwfeng

     deallocate(readIndices)

   end subroutine iamas_io_setup_vertex_block_fields !}}}


   !***********************************************************************
   !
   !  routine get_dimlist_for_field
   !
   !> \brief   Returns an array of dimension names for the specified field
   !> \author  Michael Duda
   !> \date    06 March 2015
   !> \details
   !>  Given a pool containing the field indicated by the 'fieldName' argument,
   !>  this routine allocates the output array, 'dimNames', to the
   !>  dimensionality of the field and fills the array with the dimension names.
   !>  
   !>  If the specified field is a var_array, the dimNames array will have one
   !>  fewer dimension than that of the var_array, and it will be filled with 
   !>  the outer-most dimension names.
   !>
   !>  If the specified field is is inactive (in a package sense) at the time 
   !>  of the call to this routine, the dimNames array will be returned with 
   !>  zero size.
   !
   !-----------------------------------------------------------------------
   subroutine get_dimlist_for_field(allFields, fieldName, dimNames)!{{{

      implicit none

      type (iamas_pool_type), intent(in) :: allFields
      character (len=*), intent(in) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: dimNames

      type (iamas_pool_field_info_type) :: fieldInfo
      type (field1DReal), pointer :: real1DField
      type (field2DReal), pointer :: real2DField
      type (field3DReal), pointer :: real3DField
      type (field4DReal), pointer :: real4DField
      type (field5DReal), pointer :: real5DField
      type (field1DInteger), pointer :: int1DField
      type (field2DInteger), pointer :: int2DField
      type (field3DInteger), pointer :: int3DField
      type (field1DChar), pointer :: char1DField


      call iamas_pool_get_field_info(allFields, fieldName, fieldInfo)

      if ( .not. fieldInfo % isActive ) then
         allocate(dimNames(0))
         return
      end if

      if ( fieldInfo % fieldType == iAMAS_POOL_REAL ) then

         if ( fieldInfo % nDims == 1 ) then
            call iamas_pool_get_field(allFields, fieldName, real1DField)
            if ( .not. real1DField % isVarArray ) then
               allocate(dimNames(1))
               dimNames(:) = real1DField % dimNames(:)
            else
               allocate(dimNames(0))
            end if

         else if ( fieldInfo % nDims == 2 ) then
            call iamas_pool_get_field(allFields, fieldName, real2DField)
            if ( .not. real2DField % isVarArray ) then
               allocate(dimNames(2))
               dimNames(:) = real2DField % dimNames(:)
            else
               allocate(dimNames(1))
               dimNames(1:1) = real2DField % dimNames(2:2)
            end if

         else if ( fieldInfo % nDims == 3 ) then
            call iamas_pool_get_field(allFields, fieldName, real3DField)
            if ( .not. real3DField % isVarArray ) then
               allocate(dimNames(3))
               dimNames(:) = real3DField % dimNames(:)
            else
               allocate(dimNames(2))
               dimNames(1:2) = real3DField % dimNames(2:3)
            end if

         else if ( fieldInfo % nDims == 4 ) then
            call iamas_pool_get_field(allFields, fieldName, real4DField)
            if ( .not. real4DField % isVarArray ) then
               allocate(dimNames(4))
               dimNames(:) = real4DField % dimNames(:)
            else
               allocate(dimNames(3))
               dimNames(1:3) = real4DField % dimNames(2:4)
            end if

         else if ( fieldInfo % nDims == 5 ) then
            call iamas_pool_get_field(allFields, fieldName, real5DField)
            if ( .not. real5DField % isVarArray ) then
               allocate(dimNames(5))
               dimNames(:) = real5DField % dimNames(:)
            else
               allocate(dimNames(4))
               dimNames(1:4) = real5DField % dimNames(2:5)
            end if

         else
            allocate(dimNames(0))
         end if

      else if ( fieldInfo % fieldType == iAMAS_POOL_INTEGER ) then

         if ( fieldInfo % nDims == 1 ) then
            call iamas_pool_get_field(allFields, fieldName, int1DField)
            if ( .not. int1DField % isVarArray ) then
               allocate(dimNames(1))
               dimNames(:) = int1DField % dimNames(:)
            else
               allocate(dimNames(0))
            end if

         else if ( fieldInfo % nDims == 2 ) then
            call iamas_pool_get_field(allFields, fieldName, int2DField)
            if ( .not. int2DField % isVarArray ) then
               allocate(dimNames(2))
               dimNames(:) = int2DField % dimNames(:)
            else
               allocate(dimNames(1))
               dimNames(1:1) = int2DField % dimNames(2:2)
            end if

         else if ( fieldInfo % nDims == 3 ) then
            call iamas_pool_get_field(allFields, fieldName, int3DField)
            if ( .not. int3DField % isVarArray ) then
               allocate(dimNames(3))
               dimNames(:) = int3DField % dimNames(:)
            else
               allocate(dimNames(2))
               dimNames(1:2) = int3DField % dimNames(2:3)
            end if

         else
            allocate(dimNames(0))
         end if

      else if ( fieldInfo % fieldType == iAMAS_POOL_CHARACTER ) then

         if ( fieldInfo % nDims == 1 ) then
            call iamas_pool_get_field(allFields, fieldName, char1DField)
            if ( .not. char1DField % isVarArray ) then
               allocate(dimNames(1))
               dimNames(:) = char1DField % dimNames(:)
            else
               allocate(dimNames(0))
            end if

         else
            allocate(dimNames(0))
         end if

      else

         allocate(dimNames(0))

      end if

   end subroutine get_dimlist_for_field!}}}

end module iamas_bootstrapping
