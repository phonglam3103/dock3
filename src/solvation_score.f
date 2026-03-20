      module solvation_score

      implicit none
      contains

c----------------------------------------------------------------------
c  Prorate desolvation using vol. occluded by receptor
c-----------------------------------------------------------------------
      subroutine calc_solvation_score(atomstart, atomend,
     &    atom_vdwtype, atom_polsolv, atom_apolsolv, coord_index, 
     &    transfm_coords, sperang, smax, scadif, 
     &    solgrd_precomp, hsolgrd_precomp, hsolflag, ldesolv,
     &    maxatm, maxatmpos,
     &    solvp, solva)
      
c XXX: GFortran needs this to come first (TS)
c maximum values, from max.h
      integer, intent(in) :: maxatm !how many atoms are in the ligand
      integer, intent(in) :: maxatmpos !max number of positions all atoms can take
c input values
      integer, intent(in) :: atomstart, atomend !which atoms to calculate
      integer, intent(in) :: atom_vdwtype(maxatm) !atom dock vdw type
      real, intent(in) :: atom_polsolv(maxatm) !atom polar solvation
      real, intent(in) :: atom_apolsolv(maxatm) !atom apolar solvation
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: transfm_coords(3, maxatm) !ligand atoms coordinates moved into bs
      integer, intent(in) :: sperang ! grid points/angstrom for desolvation grid.
      real, intent(in) :: smax(3) !maximum real space xyz values
      integer, intent(in) :: scadif(3) !dimensions of the solvation grid
      real, intent(in) :: solgrd_precomp(8, 
     &    0:scadif(1), 0:scadif(2), 0:scadif(3))
      real, intent(in) :: hsolgrd_precomp(8, 
     &    0:scadif(1), 0:scadif(2), 0:scadif(3))
      logical, intent(in) :: hsolflag !whether or not there is a separate hydrogen solvmap
      integer, intent(in) :: ldesolv !0 = none, 1 = full, 2 = partial
c return values
      real, intent(inout) :: solvp, solva !polar and apolar desolvation

      real    volocl ! % desolvation based on sol. volume occluded by receptor.

c     use precomputed values, removed hierarchy knowledge
c     REMOVED OUTSIDE GRID CHECKS

      integer :: i, n_in_block, block_idx, atomindex, count, nx, ny, nz
      real :: slx_tmp(8), sly_tmp(8), slz_tmp(8)
      real :: xgr_tmp(8), ygr_tmp(8), zgr_tmp(8)
      integer :: nx_tmp(8), ny_tmp(8), nz_tmp(8)
      integer :: atom_idx_tmp(8)
      real :: a_val(8, 8) ! [coeff_idx, atom_in_block]
      real :: volocl_tmp(8)

c     use precomputed values, removed hierarchy knowledge
c     REMOVED OUTSIDE GRID CHECKS

      volocl = 0.
      if (ldesolv .eq. 0) then
        volocl = 0. !no ligand desolvation at all
      else if (ldesolv .eq. 1) then
        volocl = 1. !full ligand desolvation
      endif !otherwise ldesolv = 2, so use partial
      solvp = 0.0
      solva = 0.0

      do block_idx = atomstart, atomend, 8
        n_in_block = min(8, atomend - block_idx + 1)

        ! 1. Gather coordinates and compute grid indices (only for ldesolv=2)
        if (ldesolv .eq. 2) then
          do i = 1, n_in_block
            count = block_idx + i - 1
            atomindex = coord_index(2, count)
            atom_idx_tmp(i) = atomindex
            
            slx_tmp(i) = smax(1) - 
     &                   transfm_coords(1, atomindex) * float(sperang)
            sly_tmp(i) = smax(2) - 
     &                   transfm_coords(2, atomindex) * float(sperang)
            slz_tmp(i) = smax(3) - 
     &                   transfm_coords(3, atomindex) * float(sperang)
            
            nx_tmp(i) = int(slx_tmp(i))
            ny_tmp(i) = int(sly_tmp(i))
            nz_tmp(i) = int(slz_tmp(i))

            xgr_tmp(i) = slx_tmp(i) - float(nx_tmp(i))
            ygr_tmp(i) = sly_tmp(i) - float(ny_tmp(i))
            zgr_tmp(i) = slz_tmp(i) - float(nz_tmp(i))
          enddo

          ! 2. Gather interpolation coefficients
          do i = 1, n_in_block
            atomindex = atom_idx_tmp(i)
            nx = nx_tmp(i)
            ny = ny_tmp(i)
            nz = nz_tmp(i)
            
            if (hsolflag .and. ((atom_vdwtype(atomindex) .eq. 6) .or.
     &          (atom_vdwtype(atomindex) .eq. 7))) then !hydrogen atom 
              a_val(1, i) = hsolgrd_precomp(1, nx, ny, nz)
              a_val(2, i) = hsolgrd_precomp(2, nx, ny, nz)
              a_val(3, i) = hsolgrd_precomp(3, nx, ny, nz)
              a_val(4, i) = hsolgrd_precomp(4, nx, ny, nz)
              a_val(5, i) = hsolgrd_precomp(5, nx, ny, nz)
              a_val(6, i) = hsolgrd_precomp(6, nx, ny, nz)
              a_val(7, i) = hsolgrd_precomp(7, nx, ny, nz)
              a_val(8, i) = hsolgrd_precomp(8, nx, ny, nz)
            else !heavy atom
              a_val(1, i) = solgrd_precomp(1, nx, ny, nz)
              a_val(2, i) = solgrd_precomp(2, nx, ny, nz)
              a_val(3, i) = solgrd_precomp(3, nx, ny, nz)
              a_val(4, i) = solgrd_precomp(4, nx, ny, nz)
              a_val(5, i) = solgrd_precomp(5, nx, ny, nz)
              a_val(6, i) = solgrd_precomp(6, nx, ny, nz)
              a_val(7, i) = solgrd_precomp(7, nx, ny, nz)
              a_val(8, i) = solgrd_precomp(8, nx, ny, nz)
            endif
          enddo

          ! 3. Compute volocl using FMA Factorization
          do i = 1, n_in_block
            volocl_tmp(i) = zgr_tmp(i) * (ygr_tmp(i) * 
     &                      (xgr_tmp(i) * a_val(1,i) + a_val(4,i)) + 
     &                      (xgr_tmp(i) * a_val(3,i) + a_val(7,i))) + 
     &                      (ygr_tmp(i) * (xgr_tmp(i) * a_val(2,i) + 
     &                      a_val(6,i)) + (xgr_tmp(i) * a_val(5,i) + 
     &                      a_val(8,i)))
            
            if (volocl_tmp(i) .ge. 1.0 ) then
              volocl_tmp(i) = 1.0
            elseif ( volocl_tmp(i) .lt. 0.0) then
              volocl_tmp(i) = 0.0
            endif
          enddo
        else ! ldesolv is 0 or 1
          do i = 1, n_in_block
            count = block_idx + i - 1
            atom_idx_tmp(i) = coord_index(2, count)
            volocl_tmp(i) = volocl
          enddo
        endif

        ! 4. Sum up results
        do i = 1, n_in_block
          atomindex = atom_idx_tmp(i)
          solvp = solvp-(atom_polsolv(atomindex)*volocl_tmp(i))
          solva = solva-(atom_apolsolv(atomindex)*volocl_tmp(i))
        enddo
      enddo

      return
      end subroutine calc_solvation_score

      end module solvation_score
