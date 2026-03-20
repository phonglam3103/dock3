c----------------------------------------------------------------------
c
c  --calculates force field scores.                      ECMeng  11/91
c  --removed eescore, just calculates van der waal scores
c  --removed all hierarchy stuff from scoring
c  --removed global variables. rgc 12/2011.
c------------------------------------------------------------------------
      module chemscore

      ! XXX: GFortran chokes on this when handling abval_precomp types
      !implicit none
      contains

      subroutine calc_chemscore(atomstart, atomend, 
     &    coord_index, transfm_coords, offset, invgrid, grdpts,
     &    sra, srb, atom_vdwtype,
     &    abval_precomp, ngrd, 
     &    maxatm, maxatmpos, maxtyv, 
     &    bmpvdw, vdwasum, vdwbsum, numscored, bumped)

      use vdwconsts

c XXX: GFortran needs these first (TS)
c max constants from max.h
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in) :: maxtyv !how many vdw atom types there are
c input variables
      integer, intent(in) :: atomstart, atomend !atoms to score
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: transfm_coords(3, maxatm) !ligand atoms coordinates moved into bs
      real, intent(in) :: offset(3) !grid xmin, ymin, zmin in angstroms
      real, intent(in) :: invgrid !angstroms per grid
      integer, intent(in) :: grdpts(3) !number of grid points in each dimenion (x,y,z)
      real, intent(in) :: sra(maxtyv), srb(maxtyv) !square roots of vdw parameters per vdwtype
      integer, intent(in) :: atom_vdwtype(maxatm) !the vdw type of each ligand atom
      real, intent(in) :: abval_precomp(2, 8, ngrd) !actual grid values of vdw scores,
          ! precomputed a la the Carchia optimizations
      integer, intent(in) :: ngrd !how many grid points there are, dynamically changes each run
      real, intent(in) :: bmpvdw !bump threshold for early exit
c return values
      real, intent(inout) :: vdwasum, vdwbsum !where to put scores
      integer, intent(inout) :: numscored !keeps track of how many atom positions are evaluated
      logical, intent(out) :: bumped !true if early exit due to bump

      real avdw, bvdw
c       avdw: vdw repulsion energy of a ligand atom
c       bvdw: vdw dispersion energy of a ligand atom

      integer :: i, n_in_block, block_idx, atomindex, count, array_index
      real :: x_tmp(8), y_tmp(8), z_tmp(8)
      real :: xgr_tmp(8), ygr_tmp(8), zgr_tmp(8)
      integer :: nx_tmp(8), ny_tmp(8), nz_tmp(8)
      integer :: array_idx_tmp(8)
      integer :: atom_idx_tmp(8)
      real :: a_val(8, 8), b_val(8, 8) ! [coeff_idx, atom_in_block]
      real :: apt_tmp(8), bpt_tmp(8)

      vdwasum = 0.0
      vdwbsum = 0.0
      bumped = .false.

      do block_idx = atomstart, atomend, 8
        n_in_block = min(8, atomend - block_idx + 1)

        ! 1. Gather coordinates and compute grid indices
        do i = 1, n_in_block
          count = block_idx + i - 1
          atomindex = coord_index(2, count)
          atom_idx_tmp(i) = atomindex
          
          x_tmp(i) = (transfm_coords(1, atomindex) - offset(1)) 
     &               * invgrid
          y_tmp(i) = (transfm_coords(2, atomindex) - offset(2)) 
     &               * invgrid
          z_tmp(i) = (transfm_coords(3, atomindex) - offset(3)) 
     &               * invgrid

          nx_tmp(i) = int(x_tmp(i))
          ny_tmp(i) = int(y_tmp(i))
          nz_tmp(i) = int(z_tmp(i))

          xgr_tmp(i) = x_tmp(i) - float(nx_tmp(i))
          ygr_tmp(i) = y_tmp(i) - float(ny_tmp(i))
          zgr_tmp(i) = z_tmp(i) - float(nz_tmp(i))

          array_idx_tmp(i) = grdpts(1)*grdpts(2)*nz_tmp(i) + 
     &                       grdpts(1)*ny_tmp(i) + nx_tmp(i) + 1
        enddo

        ! 2. Gather interpolation coefficients
        do i = 1, n_in_block
          array_index = array_idx_tmp(i)
          ! AVAL coefficients
          a_val(1, i) = abval_precomp(AVAL_INDEX, 1, array_index)
          a_val(2, i) = abval_precomp(AVAL_INDEX, 2, array_index)
          a_val(3, i) = abval_precomp(AVAL_INDEX, 3, array_index)
          a_val(4, i) = abval_precomp(AVAL_INDEX, 4, array_index)
          a_val(5, i) = abval_precomp(AVAL_INDEX, 5, array_index)
          a_val(6, i) = abval_precomp(AVAL_INDEX, 6, array_index)
          a_val(7, i) = abval_precomp(AVAL_INDEX, 7, array_index)
          a_val(8, i) = abval_precomp(AVAL_INDEX, 8, array_index)
          
          ! BVAL coefficients
          b_val(1, i) = abval_precomp(BVAL_INDEX, 1, array_index)
          b_val(2, i) = abval_precomp(BVAL_INDEX, 2, array_index)
          b_val(3, i) = abval_precomp(BVAL_INDEX, 3, array_index)
          b_val(4, i) = abval_precomp(BVAL_INDEX, 4, array_index)
          b_val(5, i) = abval_precomp(BVAL_INDEX, 5, array_index)
          b_val(6, i) = abval_precomp(BVAL_INDEX, 6, array_index)
          b_val(7, i) = abval_precomp(BVAL_INDEX, 7, array_index)
          b_val(8, i) = abval_precomp(BVAL_INDEX, 8, array_index)
        enddo

        ! 3. Compute interpolated values using FMA Factorization
        do i = 1, n_in_block
          apt_tmp(i) = zgr_tmp(i) * (ygr_tmp(i) * 
     &                 (xgr_tmp(i) * a_val(1,i) + a_val(4,i)) + 
     &                 (xgr_tmp(i) * a_val(3,i) + a_val(7,i))) + 
     &                 (ygr_tmp(i) * (xgr_tmp(i) * a_val(2,i) + 
     &                 a_val(6,i)) + (xgr_tmp(i) * a_val(5,i) + 
     &                 a_val(8,i)))

          bpt_tmp(i) = zgr_tmp(i) * (ygr_tmp(i) * 
     &                 (xgr_tmp(i) * b_val(1,i) + b_val(4,i)) + 
     &                 (xgr_tmp(i) * b_val(3,i) + b_val(7,i))) + 
     &                 (ygr_tmp(i) * (xgr_tmp(i) * b_val(2,i) + 
     &                 b_val(6,i)) + (xgr_tmp(i) * b_val(5,i) + 
     &                 b_val(8,i)))
        enddo

        ! 4. Sum up the results
        do i = 1, n_in_block
          atomindex = atom_idx_tmp(i)
          avdw = sra(atom_vdwtype(atomindex)) * apt_tmp(i) 
          bvdw = srb(atom_vdwtype(atomindex)) * bpt_tmp(i)
          vdwasum = vdwasum + avdw
          vdwbsum = vdwbsum + bvdw
        enddo

        ! EARLY EXIT: Check if bumped after each block
        if ((vdwasum - vdwbsum) .gt. bmpvdw) then
          bumped = .true.
          numscored = numscored + 1
          return
        endif
      enddo

      numscored = numscored + 1
      return
      end subroutine calc_chemscore

      end module chemscore

