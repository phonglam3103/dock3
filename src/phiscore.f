!module for phimap scoring, contains phipdb & charge (also desolvation routines)
      module phiscore

      implicit none
      contains     

************************************************************************
c   This is Honig et al.'s phitopdb program, turned into a subroutine
c   for dock2 scoring.  Modifications by BKS.
c   Brian K. Shoichet, 08/90.
c   removed group/hierarchy knowledge rgc 12/2009
c   removed global variables rgc 12/2011
************************************************************************
      subroutine phipdb(atomstart, atomend, 
     &    coord_index, transfm_coords, oldmid, phiscale, goff,
     &    phimap_precomp, maxatm, maxatmpos, nsize, potent)

c XXX: These need to come first or Gfortran complains (TS)
c max constants, see max.h
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in) :: nsize !cubic grid size of phimap
c input variables
      integer, intent(in) :: atomstart, atomend !where the atoms and other things are
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: transfm_coords(3, maxatm) !ligand atoms coordinates moved into bs
      real, intent(in) :: oldmid(3) !the center of the phimap
      real, intent(in) :: phiscale !angstroms per grid point on phimap
      real, intent(in) :: goff !distance from center to lower grid edges, cubic so all same
      real, intent(in) :: phimap_precomp(8, nsize, nsize, nsize) !actual phimap values
c return variable
      real, intent(inout) :: potent(maxatm) !electrostatic potential at each atom


      integer :: i, n_in_block, block_idx, atomindex, count, nx, ny, nz
      real :: x_tmp(8), y_tmp(8), z_tmp(8)
      real :: xgr_tmp(8), ygr_tmp(8), zgr_tmp(8)
      integer :: nx_tmp(8), ny_tmp(8), nz_tmp(8)
      integer :: atom_idx_tmp(8)
      real :: a_val(8, 8) ! [coeff_idx, atom_in_block]

c     loop over all atoms in blocks of 8
      do block_idx = atomstart, atomend, 8
        n_in_block = min(8, atomend - block_idx + 1)

        ! 1. Gather coordinates and compute grid indices
        do i = 1, n_in_block
          count = block_idx + i - 1
          atomindex = coord_index(2, count)
          atom_idx_tmp(i) = atomindex
          
          x_tmp(i) = (transfm_coords(1, atomindex) - oldmid(1)) 
     &               * phiscale + goff
          y_tmp(i) = (transfm_coords(2, atomindex) - oldmid(2)) 
     &               * phiscale + goff
          z_tmp(i) = (transfm_coords(3, atomindex) - oldmid(3)) 
     &               * phiscale + goff

          nx_tmp(i) = int(x_tmp(i))
          ny_tmp(i) = int(y_tmp(i))
          nz_tmp(i) = int(z_tmp(i))

          xgr_tmp(i) = x_tmp(i) - float(nx_tmp(i))
          ygr_tmp(i) = y_tmp(i) - float(ny_tmp(i))
          zgr_tmp(i) = z_tmp(i) - float(nz_tmp(i))
        enddo

        ! 2. Gather interpolation coefficients
        do i = 1, n_in_block
          nx = nx_tmp(i)
          ny = ny_tmp(i)
          nz = nz_tmp(i)
          ! These 8 values are contiguous in phimap_precomp
          a_val(1, i) = phimap_precomp(1, nx, ny, nz)
          a_val(2, i) = phimap_precomp(2, nx, ny, nz)
          a_val(3, i) = phimap_precomp(3, nx, ny, nz)
          a_val(4, i) = phimap_precomp(4, nx, ny, nz)
          a_val(5, i) = phimap_precomp(5, nx, ny, nz)
          a_val(6, i) = phimap_precomp(6, nx, ny, nz)
          a_val(7, i) = phimap_precomp(7, nx, ny, nz)
          a_val(8, i) = phimap_precomp(8, nx, ny, nz)
        enddo

        ! 3. Compute potential using FMA Factorization
        do i = 1, n_in_block
          atomindex = atom_idx_tmp(i)
          potent(atomindex) = zgr_tmp(i) * (ygr_tmp(i) * 
     &                        (xgr_tmp(i) * a_val(1,i) + a_val(4,i)) + 
     &                        (xgr_tmp(i) * a_val(3,i) + a_val(7,i))) + 
     &                        (ygr_tmp(i) * (xgr_tmp(i) * a_val(2,i) + 
     &                        a_val(6,i)) + (xgr_tmp(i) * a_val(5,i) + 
     &                        a_val(8,i)))
        enddo
      enddo

      return
      end subroutine phipdb

c  --calculates electrostatic interaction energy using DelPhi
c    potential and ligand atomic charges.  Units are kT.  BKS
c removed global variables, etc, rgc 12/2011
      function charge(atomstart, atomend, coord_index, atom_charge, 
     &    potent, maxatm, maxatmpos) result(chargesum)

c XXX: These need to come first or Gfortran complains (TS)
c max constants
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take

c input variables
      integer, intent(in) :: atomstart, atomend
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: atom_charge(maxatm) !partial charge at each atom
      real, intent(in) :: potent(maxatm) !electrostatic potential at each atom

c output values
      !real, intent(out) :: chargesum ! GFortran Chokes on this
      real :: chargesum

      integer atomindex, count !used to find correct atoms

      chargesum = 0.0 !the returned value
      do count = atomstart, atomend !for every atom
        atomindex = coord_index(2, count) !gets actual atom non-consecutive #
        chargesum = chargesum + (potent(atomindex) * 
     &      atom_charge(atomindex))
      enddo
      chargesum = chargesum * 0.5924 !const to convert to kcal/mole
      return
      end function charge

c  --calculates electrostatic interaction energy using DelPhi
c    potential and ligand atomic charges.  Units are kT.  
c this uses charge squared, to scale ligand desolvation terms correctly
c removed global variables, etc, rgc 12/2011
      function charge_squared(atomstart, atomend, coord_index, 
     &    atom_charge, potent, maxatm, maxatmpos) result(chargesum)

c XXX: These need to come first or Gfortran complains (TS)
c max constants
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take

c input variables
      integer, intent(in) :: atomstart, atomend
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: atom_charge(maxatm) !partial charge at each atom
      real, intent(in) :: potent(maxatm) !electrostatic potential at each atom

c output values
      !real, intent(out) :: chargesum ! GFortran Chokes on this
      real :: chargesum

      integer atomindex, count !used to find correct atoms

      chargesum = 0.0 !the returned value
      do count = atomstart, atomend !for every atom
        atomindex = coord_index(2, count) !gets actual atom non-consecutive #
        chargesum = chargesum + (potent(atomindex) * 
     &      atom_charge(atomindex) * atom_charge(atomindex)) !this is where
                           !the charged squared term is applied
      enddo
      chargesum = chargesum * 0.5924 !const to convert to kcal/mole
      return
      end function charge_squared

!function basically converts kT units to kcal/mol, used for receptor desolvation
!uses atom_heavy_valence to reduce the effect by 1/3 for each heavy atom
!connected. atom_heavy_valence should be 0..3 inclusive, above 3 makes no sense
!and neither does below 0. no need to reference options0%recdes_valence_penalty
!since if it isn't set then everything will be 0.
      function convert_sum(atomstart, atomend, coord_index,  
     &    kT_energy, atom_vdwtype, atom_heavy_valence, 
     &    maxatm, maxatmpos) 
     &    result(kcal_energy_sum)

c XXX: These need to come first or Gfortran complains (TS)
c max constants
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take

c input variables
      integer, intent(in) :: atomstart, atomend
      integer, intent(in) :: coord_index(3, maxatmpos) !coord#, atom#, conf# for each coord
      real, intent(in) :: kT_energy(maxatm) !energy in kT at each atom
      integer, intent(in) :: atom_vdwtype(maxatm) !atom dock vdw type
      integer, intent(in) :: atom_heavy_valence(maxatm) !atom heavy valence
c output values
      !real, intent(out) :: kcal_energy_sum ! GFortran chokes on this (TS)
      real :: kcal_energy_sum

      integer atomindex, count !used to find correct atoms
      real this_energy !used to temporarily store energy

      kcal_energy_sum = 0.0 !the returned value
      do count = atomstart, atomend !for every atom
        atomindex = coord_index(2, count) !gets actual atom non-consecutive #
        !ignore hydrogen atoms for receptor desolvation, they don't have enough
        !bulk to desolvate, only the heavy atoms do
        if ((atom_vdwtype(atomindex) .ne. 6) .and. 
     &      (atom_vdwtype(atomindex) .ne. 7)) then
          this_energy = kT_energy(atomindex) 
!          write (6,*) 'xxxpre', this_energy !debugging
          this_energy = this_energy * 
     &        (5 - atom_heavy_valence(atomindex))/5.0
!          write (6,*) 'xxxpost', this_energy, 
!     &        atom_heavy_valence(atomindex)
          kcal_energy_sum = kcal_energy_sum + this_energy
        endif
      enddo
      kcal_energy_sum = kcal_energy_sum * 0.5924 !const to convert to kcal/mole
      return
      end function convert_sum

      end module phiscore
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 David M. Lorber and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c       Copyright (C) 1992 Regents of the University of California
c                      Brian K. Shoichet and I.D. Kuntz
c                         All Rights Reserved.
c
c----------------------------------------------------------------------
