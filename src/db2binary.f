c----------------------------------------------------------------------
c
c  Binary DB2 I/O module for high-performance ligand reading
c  Provides gzip-compatible binary format (db2bin.gz) 
c  ~5-20x faster than text parsing
c  
c  Phong Lam, Uppsala University, January 2026
c
c----------------------------------------------------------------------
      module db2binary

      use db2type
      use optionstype
      use filenums
      use spheres

      implicit none

c Binary format magic and version
      character(6), parameter :: DB2BIN_MAGIC = 'DB2BIN'
      integer, parameter :: DB2BIN_VERSION = 1

      contains

c----------------------------------------------------------------------
c Helper: read integer from gzipped binary file
c----------------------------------------------------------------------
      subroutine read_int_gz(fdlig, val, istat)
      integer(kind=8), intent(in) :: fdlig
      integer, intent(out) :: val
      integer, intent(out) :: istat
      character(4) :: buf
      
      call gzbread(fdlig, buf, 4, istat)
      val = transfer(buf, val)
      end subroutine read_int_gz

c----------------------------------------------------------------------
c Helper: read real from gzipped binary file
c----------------------------------------------------------------------
      subroutine read_real_gz(fdlig, val, istat)
      integer(kind=8), intent(in) :: fdlig
      real, intent(out) :: val
      integer, intent(out) :: istat
      character(4) :: buf
      
      call gzbread(fdlig, buf, 4, istat)
      val = transfer(buf, val)
      end subroutine read_real_gz

c----------------------------------------------------------------------
c Helper: read integer array from gzipped binary file  
c----------------------------------------------------------------------
      subroutine read_int_array_gz(fdlig, arr, n, istat)
      integer(kind=8), intent(in) :: fdlig
      integer, intent(in) :: n
      integer, intent(out) :: arr(n)
      integer, intent(out) :: istat
      character, allocatable :: buf(:)
      integer :: nbytes
      
      nbytes = n * 4
      allocate(buf(nbytes))
      call gzbread(fdlig, buf, nbytes, istat)
      arr = transfer(buf, arr)
      deallocate(buf)
      end subroutine read_int_array_gz

c----------------------------------------------------------------------
c Helper: read real array from gzipped binary file
c----------------------------------------------------------------------
      subroutine read_real_array_gz(fdlig, arr, n, istat)
      integer(kind=8), intent(in) :: fdlig
      integer, intent(in) :: n
      real, intent(out) :: arr(n)
      integer, intent(out) :: istat
      character, allocatable :: buf(:)
      integer :: nbytes
      
      nbytes = n * 4
      allocate(buf(nbytes))
      call gzbread(fdlig, buf, nbytes, istat)
      arr = transfer(buf, arr)
      deallocate(buf)
      end subroutine read_real_array_gz

c----------------------------------------------------------------------
c read_db2_binary: Fast binary reader for db2bin.gz files
c Replaces ligread2 when use_binary_db2 = .true.
c----------------------------------------------------------------------
      subroutine read_db2_binary(db2lig, options0, 
     &    inputend, rigid_heavy_coords, 
     &    rigid_colors, rigid_heavy_count, fdlig, ligc, db2c)

      type(db2), intent(inout) :: db2lig
      type(options), intent(inout) :: options0
      logical, intent(inout) :: inputend
      real, intent(inout), dimension(3, MAXPTS) :: rigid_heavy_coords
      integer, intent(inout), dimension(MAXPTS) :: rigid_colors
      integer, intent(inout) :: rigid_heavy_count
      integer(kind=8), intent(inout) :: fdlig
      integer, intent(inout) :: ligc
      integer, intent(inout) :: db2c

c temporary variables
      character(6) :: magic
      integer :: version, inputstatus
      integer :: total_atoms, total_bonds, total_coords
      integer :: total_confs, total_sets, total_clusters
      integer :: mlines, i, j, curset, curconf
      integer :: arb_count
      real :: total_charge, total_polar_solv, total_apolar_solv
      real :: total_solv, total_surfarea
      real :: set_energy_val, set_max_strain_val
      integer :: conf_len, broken_val, conf_val
      integer :: set_start, set_end, match_start, match_end
      integer :: match_color
      real :: match_x, match_y, match_z
      integer, allocatable :: temp_int_arr(:)
      real, allocatable :: temp_real_arr(:)

c read magic number (6 bytes)
      call gzbread(fdlig, magic, 6, inputstatus)
      if (inputstatus .le. 0) then
        inputend = .true.
        return
      endif
      
      if (magic .ne. DB2BIN_MAGIC) then
        write(OUTDOCK, *) 'ERROR: Invalid binary DB2 magic: ', magic
        write(OUTDOCK, *) 'Expected: ', DB2BIN_MAGIC
        inputend = .true.
        return
      endif

c read version
      call read_int_gz(fdlig, version, inputstatus)
      if (version .ne. DB2BIN_VERSION) then
        write(OUTDOCK, *) 'WARNING: Binary DB2 version mismatch'
        write(OUTDOCK, *) 'File version: ', version
        write(OUTDOCK, *) 'Expected: ', DB2BIN_VERSION
      endif

c read header integers
      call read_int_gz(fdlig, total_atoms, inputstatus)
      call read_int_gz(fdlig, total_bonds, inputstatus)
      call read_int_gz(fdlig, total_coords, inputstatus)
      call read_int_gz(fdlig, total_confs, inputstatus)
      call read_int_gz(fdlig, total_sets, inputstatus)
      call read_int_gz(fdlig, rigid_heavy_count, inputstatus)
      call read_int_gz(fdlig, total_clusters, inputstatus)
      call read_int_gz(fdlig, mlines, inputstatus)

c read header reals
      call read_real_gz(fdlig, total_charge, inputstatus)
      call read_real_gz(fdlig, total_polar_solv, inputstatus)
      call read_real_gz(fdlig, total_apolar_solv, inputstatus)
      call read_real_gz(fdlig, total_solv, inputstatus)
      call read_real_gz(fdlig, total_surfarea, inputstatus)

c read refcod (16 bytes, padded) and trim NUL characters
      call gzbread(fdlig, db2lig%refcod, 16, inputstatus)
c replace NUL chars in refcod with spaces
      do i = 1, 16
        if (db2lig%refcod(i:i) .eq. char(0)) then
          db2lig%refcod(i:i) = ' '
        endif
      enddo

c read protcod (9 bytes)
      call gzbread(fdlig, db2lig%protcod, 9, inputstatus)

c read smiles (78 bytes)
      call gzbread(fdlig, db2lig%smiles, 78, inputstatus)

c read molname (78 bytes)
      call gzbread(fdlig, db2lig%molname, 78, inputstatus)

c read arbitrary lines count and data
      call read_int_gz(fdlig, arb_count, inputstatus)
      
c store header values (update mlines to include arbitrary)
      db2lig%total_atoms = total_atoms
      db2lig%total_bonds = total_bonds
      db2lig%total_coords = total_coords
      db2lig%total_confs = total_confs
      db2lig%total_sets = total_sets
      db2lig%total_clusters = total_clusters
      db2lig%mlines = mlines
      db2lig%total_charge = total_charge
      db2lig%total_polar_solv = total_polar_solv
      db2lig%total_apolar_solv = total_apolar_solv
      db2lig%total_solv = total_solv
      db2lig%total_surfarea = total_surfarea

c allocate arrays
      call allocate_db2(db2lig, OUTDOCK, total_sets,
     &    total_coords, total_atoms, total_bonds,
     &    total_confs, mlines, total_clusters,
     &    options0%dockovalent, options0%nsav)

c initialize set_conf_len to zero
      db2lig%set_conf_len(:) = 0
      
c initialize arbitrary array to spaces
      if (allocated(db2lig%arbitrary)) then
          db2lig%arbitrary = ' '
      endif

c read arbitrary lines (after allocate since it allocates the array)
      do i = 1, arb_count
        call gzbread(fdlig, db2lig%arbitrary(i), 78, inputstatus)
      enddo

c----------------------------------------------------------------------
c read atom arrays
c----------------------------------------------------------------------
      allocate(temp_int_arr(total_atoms))
      allocate(temp_real_arr(total_atoms))
      
c read atom_num (int32)
      call read_int_array_gz(fdlig, temp_int_arr, total_atoms, 
     &                       inputstatus)
      db2lig%atom_num(1:total_atoms) = temp_int_arr
      
c read atom_name (4 bytes each)
      do i = 1, total_atoms
        call gzbread(fdlig, db2lig%atom_name(i), 4, inputstatus)
      enddo
      
c read atom_type (5 bytes each)
      do i = 1, total_atoms
        call gzbread(fdlig, db2lig%atom_type(i), 5, inputstatus)
      enddo
      
c read atom_vdwtype
      call read_int_array_gz(fdlig, temp_int_arr, total_atoms, 
     &                       inputstatus)
      db2lig%atom_vdwtype(1:total_atoms) = temp_int_arr
      
c read atom_color
      call read_int_array_gz(fdlig, temp_int_arr, total_atoms,
     &                       inputstatus)
      db2lig%atom_color(1:total_atoms) = temp_int_arr
      
      call read_real_array_gz(fdlig, temp_real_arr, total_atoms,
     &                        inputstatus)
      db2lig%atom_charge(1:total_atoms) = temp_real_arr
      
      call read_real_array_gz(fdlig, temp_real_arr, total_atoms,
     &                        inputstatus)
      db2lig%atom_polsolv(1:total_atoms) = temp_real_arr
      
      call read_real_array_gz(fdlig, temp_real_arr, total_atoms,
     &                        inputstatus)
      db2lig%atom_apolsolv(1:total_atoms) = temp_real_arr
      
      call read_real_array_gz(fdlig, temp_real_arr, total_atoms,
     &                        inputstatus)
      db2lig%atom_totalsolv(1:total_atoms) = temp_real_arr
      
      call read_real_array_gz(fdlig, temp_real_arr, total_atoms,
     &                        inputstatus)
      db2lig%atom_surfarea(1:total_atoms) = temp_real_arr
      
      deallocate(temp_int_arr)
      deallocate(temp_real_arr)

c----------------------------------------------------------------------
c read bond data
c----------------------------------------------------------------------
      allocate(temp_int_arr(total_bonds))
      
c read bond_num (int32)
      call read_int_array_gz(fdlig, temp_int_arr, total_bonds,
     &                       inputstatus)
      db2lig%bond_num(1:total_bonds) = temp_int_arr
      
c read bond_start (int32)
      call read_int_array_gz(fdlig, temp_int_arr, total_bonds,
     &                       inputstatus)
      db2lig%bond_start(1:total_bonds) = temp_int_arr
      
c read bond_end (int32)
      call read_int_array_gz(fdlig, temp_int_arr, total_bonds,
     &                       inputstatus)
      db2lig%bond_end(1:total_bonds) = temp_int_arr
      
c read bond_type (2 bytes each)
      do i = 1, total_bonds
        call gzbread(fdlig, db2lig%bond_type(i), 2, inputstatus)
      enddo
      
      deallocate(temp_int_arr)

c----------------------------------------------------------------------
c read coordinate data
c----------------------------------------------------------------------
      do i = 1, total_coords
        call read_int_gz(fdlig, db2lig%coord_index(1,i), inputstatus)
        call read_int_gz(fdlig, db2lig%coord_index(2,i), inputstatus)
        call read_int_gz(fdlig, db2lig%coord_index(3,i), inputstatus)
        call read_real_gz(fdlig, db2lig%coords(1,i), inputstatus)
        call read_real_gz(fdlig, db2lig%coords(2,i), inputstatus)
        call read_real_gz(fdlig, db2lig%coords(3,i), inputstatus)
      enddo

c----------------------------------------------------------------------
c read rigid matching spheres
c----------------------------------------------------------------------
      do i = 1, rigid_heavy_count
        call read_int_gz(fdlig, rigid_colors(i), inputstatus)
        call read_real_gz(fdlig, rigid_heavy_coords(1,i), inputstatus)
        call read_real_gz(fdlig, rigid_heavy_coords(2,i), inputstatus)
        call read_real_gz(fdlig, rigid_heavy_coords(3,i), inputstatus)
      enddo

c----------------------------------------------------------------------
c read conformation coord indices
c----------------------------------------------------------------------
      do i = 1, total_confs
        call read_int_gz(fdlig, db2lig%conf_coord(1,i), inputstatus)
        call read_int_gz(fdlig, db2lig%conf_coord(2,i), inputstatus)
      enddo

c----------------------------------------------------------------------
c read set data
c----------------------------------------------------------------------
      do curset = 1, total_sets
        call read_int_gz(fdlig, conf_len, inputstatus)
        call read_int_gz(fdlig, broken_val, inputstatus)
        call read_real_gz(fdlig, set_energy_val, inputstatus)
        call read_real_gz(fdlig, set_max_strain_val, inputstatus)
        
        db2lig%set_conf_len(curset) = conf_len
        db2lig%set_broken(curset) = broken_val
        db2lig%set_energy(curset) = set_energy_val * options0%iscale
        db2lig%set_total_strain(curset) = set_energy_val
        db2lig%set_max_strain(curset) = set_max_strain_val
        
        if ((set_energy_val .ge. options0%total_strain) .or.
     &      (set_max_strain_val .ge. options0%max_strain)) then
          db2lig%set_above_strain(curset) = 1
        else
          db2lig%set_above_strain(curset) = 0
        endif
        
        do curconf = 1, conf_len
          call read_int_gz(fdlig, conf_val, inputstatus)
          db2lig%set_conf(curset, curconf) = conf_val
        enddo
      enddo

c----------------------------------------------------------------------
c read cluster data
c----------------------------------------------------------------------
      do i = 1, total_clusters
        call read_int_gz(fdlig, db2lig%cluster_set_start(i),inputstatus)
        call read_int_gz(fdlig, db2lig%cluster_set_end(i), inputstatus)
        call read_int_gz(fdlig, db2lig%cluster_match_start(i),
     &                   inputstatus)
        call read_int_gz(fdlig, db2lig%cluster_match_end(i),inputstatus)
      enddo

c read additional matching spheres for each cluster
      do i = 1, total_clusters
        do j = db2lig%cluster_match_start(i), 
     &         db2lig%cluster_match_end(i)
          call read_int_gz(fdlig, db2lig%addmatch_color(j), inputstatus)
          call read_real_gz(fdlig, db2lig%addmatch_coord(1,j),
     &                      inputstatus)
          call read_real_gz(fdlig, db2lig%addmatch_coord(2,j),
     &                      inputstatus)
          call read_real_gz(fdlig, db2lig%addmatch_coord(3,j),
     &                      inputstatus)
        enddo
      enddo

c----------------------------------------------------------------------
c post-processing: compute derived values like ligread2 does
c----------------------------------------------------------------------
c compute total heavy atom count (atoms where vdwtype != 6 and != 7)
      db2lig%total_heavy_atoms = 0
      db2lig%ligcharge = 0.0
      do i = 1, total_atoms
        if ((db2lig%atom_vdwtype(i) .ne. 6) .and. 
     &      (db2lig%atom_vdwtype(i) .ne. 7)) then
          db2lig%total_heavy_atoms = db2lig%total_heavy_atoms + 1
        endif
        db2lig%ligcharge = db2lig%ligcharge + db2lig%atom_charge(i)
      enddo

c initialize atom_heavy_valence to 0
      do i = 1, total_atoms
        db2lig%atom_heavy_valence(i) = 0
      enddo

      ligc = ligc + 1
      inputend = .false.

      return
      end subroutine read_db2_binary

      end module db2binary
