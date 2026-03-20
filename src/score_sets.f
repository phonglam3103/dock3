c just score the conformations directly, no hierarchy, no tree necessary
c this is the most complex subroutine, since it deals with clouds/clusters
c and flexible receptor scoring
      module score_sets
 
      implicit none
      contains

      subroutine calc_score_sets(current_match, dock_status, 
     &    setstart, setend, db2lig, options0, ligscore, ligscoreeach, 
     &    grids0, vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0,
     &    rec_des_r0, match,
     &    allminlimits, allmaxlimits,
     &    sra, srb, 
     &    maxor, maxconfs, maxsets, maxatm, maxatmpos,
     &    maxtyv)
 
      use status
      use db2type
      use optionstype
      use ligscoretype
      use gridstype
      use phimaptype
      use gisttype
      use solvmaptype
      use transfm_conf
      use vdwmaptype
      use matchtype
      use doall_conf
      use filenums

c constants, set in max.h and passed in
      integer, intent(in) :: maxor !the size of many arrays that store information about
                    !orientations, is checked several times
      integer, intent(in) :: maxconfs !the size of any array having to do with conformations
      integer, intent(in) :: maxsets !how many max sets can be in a ligand file
      integer, intent(in) :: maxatm !how many atoms max can be in a ligand
      integer, intent(in) :: maxatmpos !how many positions all ligand atoms can take
      integer, intent(in) :: maxtyv !how many vdw atom types there are
c variables
      integer, intent(in) :: current_match !which match number this orientation is from
      integer, intent(inout) :: dock_status !dock_status is a enumerated type, see status.h
      integer, intent(in) :: setstart, setend !controls only doing a certain number of sets
      type(db2), intent(inout) :: db2lig !ligand information
      type(options), intent(in) :: options0 !useful options
      type(ligscoret), intent(inout) :: ligscore
      type(ligscoret), intent(inout) :: ligscoreeach
      type(flexgrids), intent(inout) :: grids0 !precomp values
      type(vdwmap), intent(in) :: vdwmap0
      type(phimap), intent(in) :: phimap0, recdes0, ligdes0
      type(gistparm), intent(in) :: gist0, rec_des_r0
      type(solvmap), intent(in) :: solvmap0
      type(matcht), intent(in) :: match
      real, intent(in), dimension(3) :: allminlimits, allmaxlimits !grid limits in xyz
      real, intent(in), dimension(maxtyv) :: sra, srb !square roots of vdw parameters per vdwtype

      integer :: maingrid, flexgrid !which grid to use
      integer :: tempset !keeps track of which set we are scoring
      integer :: loopconf, tempconf, tempconfnum, oriconf !keeps track of which conf we are scoring, oriconf is the origin conformation, which conf preseids the current conf
      integer :: bestset !keeps track of the best set in terms of energy
      integer :: set_status !status.h status of the current set
      integer :: place_to_save !insertion point into list
      !not accessed outside of this function so made here
      logical, dimension(maxsets, options0%total_receptors) :: goodset !for each set overall, used for flexible receptor
      logical :: goodwhole
      logical :: onegood = .false. ! look to see if none are good. 
      real, dimension(maxsets, grids0%total_combinations) :: 
     &     wholevs, wholees, wholegi, wholeas, wholeps, 
     &     wholeds,wholerd, wholehs, wholers, wholescore !score for all the grids (or just 1 if only 1) 
      real, dimension(maxsets, options0%total_receptors) ::
     &    setscore !scores per set (whole position)
c     Packed array for better locality (Array of Structures)
c     1:vdw, 2:es, 3:gi, 4:as, 5:ps, 6:ds, 7:rd, 8:hs
      real, dimension(8, maxconfs, options0%total_receptors) :: 
     &    conf_scores 
      integer, dimension(maxconfs, options0%total_receptors) :: 
     &    conf_status !track status of each conf
      real, dimension(maxconfs, options0%total_receptors) ::
     &    conf_weighted_total !pre-calculated weighted total score
      real, dimension(maxsets, options0%total_receptors) ::
     &    set_weighted_total !accumulated weighted total per set
      character (len=255) :: reccode
      integer, dimension(grids0%group_len) :: temp_combo !list of grids
      integer :: temp_combo_num !temporary counter for combinations
      integer :: temp_combo_num_grid !temporary variable for grid index in combo
      integer :: tempcoord,tempatom,count
c     ! this will store the list grid point already asigned to the molecule 
c     integer,dimension(options0%total_receptors,10000) :: gistgridindex
      integer,dimension(5000,db2lig%total_confs,
     &     options0%total_receptors) :: 
     &     gistgridindex
c     integer,dimension(100000) :: gistgridindex_cur_rec
      integer,dimension(db2lig%total_confs,
     &    options0%total_receptors) :: lengistgi ! this is the curent length of list grid point already asigned to the molecule
c     just for debugging, jklyu, 2019.03.28
      character (len=*), parameter :: PDBFORM =
     &  '(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,11X,A2,A2)'
      integer :: atomcount,conf,confcount
      integer atomcoord!used to find the right place to put things

      logical logged_tconf, logged_oconf
c     2D bump cache for O(1) lookup: is_bumped(conf, grid)
      logical, dimension(maxconfs, options0%total_receptors) ::
     &    is_bumped
      logical :: skip_set

c     Initialize bump cache for all grids
      is_bumped(:,:) = .false.
      conf_weighted_total(:,:) = 0.0

      maingrid = grids0%group_part(1, 1) !do the invariant or only grid first

c     intialize to zero the lenth of the gist index list.  
      do tempconf = 1, db2lig%total_confs ! loop over all 
        lengistgi(tempconf,maingrid) = 0
        if (options0%flexible_receptor) then !
          do flexgrid = 2, options0%total_receptors !
             lengistgi(tempconf,flexgrid) = 0
          enddo
        endif
      enddo
      !first do the first conf, which is rigid component
      !tempconf = 1 ! db2lig%conf_sorted(1) will equal 1
      tempconf = 1
c     write(*,*) "starting tempconf = ", tempconf
      !gistgridindex_cur_rec = gistgridindex(:,maingrid)
c      write(6,*) "I AM in score_sets"
c      write(6,*) "conf_status before is",conf_status(tempconf, maingrid)
      call doflush(6)
      conf_status(tempconf, maingrid) = 
     &    run_conf(tempconf, current_match,
     &    db2lig, options0, grids0, maingrid,
     &    vdwmap0, phimap0, recdes0, ligdes0, gist0, solvmap0,
     &    rec_des_r0, match,
     &    conf_scores(1, tempconf, maingrid), 
     &    conf_scores(2, tempconf, maingrid), 
     &    conf_scores(3, tempconf, maingrid), 
     &    conf_scores(4, tempconf, maingrid), 
     &    conf_scores(5, tempconf, maingrid), 
     &    conf_scores(6, tempconf, maingrid),
     &    conf_scores(7, tempconf, maingrid),
     &    conf_scores(8, tempconf, maingrid),
     &    allminlimits, allmaxlimits,
     &    sra, srb, 
     &    maxor, maxconfs, maxatm, maxatmpos, maxtyv, 
     &    ligscore%numscored,
     &    gistgridindex(:,tempconf,
     &        maingrid),
     &    lengistgi(tempconf,maingrid),
     &    gistgridindex(:,
c     &                      db2lig%conf_conection(tempconf),
     &                      tempconf,
     &        maingrid),
     &    lengistgi(
c     &                      db2lig%conf_conection(tempconf),
     &                      tempconf,
     &        maingrid))


c      write(6,*) "conf_status after is", conf_status(tempconf, maingrid)
c      just for debugging turn off the rigidbump check
      if ((conf_status(tempconf, maingrid) .ne. ALLOKAY) .or. 
     &    (conf_scores(1, tempconf, maingrid) .gt. 
     &    options0%bmpvdwrigid)) then
        dock_status = BUMPED !the new version of bumping, rely on the vdw score
c       Mark rigid conf as bumped in cache
        is_bumped(tempconf, maingrid) = .true.
      else !the rigid component isn't terrible, go ahead and score
c       LAZY EVALUATION: Only score confs on-demand during set iteration
c       Mark all non-rigid confs as NOTCHECKED - they will be scored
c       when first encountered in the set loop below
        do loopconf = 2, db2lig%total_confs
          do flexgrid = 1, options0%total_receptors
            conf_status(loopconf, flexgrid) = NOTCHECKED
          enddo
        enddo
        do flexgrid = 1, options0%total_receptors !go through every grid & 
                            !set, assemble scores for each combination
c         Note: is_bumped is now initialized at start, not reset per grid
c         Skip entire set loop if rigid conf (conf 1) bumped on this grid
          if (is_bumped(1, flexgrid)) then
            dock_status = BUMPED
c           Mark all sets as not good for this grid
            do tempset = setstart, setend
              goodset(tempset, flexgrid) = .false.
            enddo
            cycle
          endif
            set_status = NOTCHECKED !temporary for this set
          do tempset = setstart, setend !going through every set as well
            skip_set = .false.

c           Pre-scan: Check if any conformation in this set is already known to be bumped
            do tempconfnum = 2, db2lig%set_conf_len(tempset)
              tempconf = db2lig%set_conf(tempset, tempconfnum)
              if (is_bumped(tempconf, flexgrid)) then
                 skip_set = .true.
                 exit
              endif
            enddo

            if (skip_set) then
              set_status = BUMPED
              dock_status = BUMPED
              goodset(tempset, flexgrid) = .false.
              cycle ! Skip to next set
            endif

            !add up scores from each conf to find the best one.
c           Init weighted total
            set_weighted_total(tempset, flexgrid) = 0.
            !broken/clashed check on next line, if necessary
            if ((db2lig%set_broken(tempset) .eq. 0) .or. 
     &            (options0%check_clash .eqv. .false.)) then 
              if ((db2lig%set_above_strain(tempset) .eq. 0) .or.
     &            (options0%check_strain .eqv. .false.)) then 
                set_status = ALLOKAY
                do tempconfnum = 1, db2lig%set_conf_len(tempset)
                  tempconf = db2lig%set_conf(tempset, tempconfnum)
c                 Check bump cache first (O(1) lookup)
                  if (is_bumped(tempconf, flexgrid)) then
                    set_status = BUMPED
                    exit
                  endif
                  !build on the fly if necessary
                  if (conf_status(tempconf, flexgrid) .eq. 
     &                NOTCHECKED) then !if it hasn't been checked already
                    conf_status(tempconf, flexgrid) = run_conf(
     &                  tempconf, current_match,
     &                  db2lig, options0, grids0, flexgrid,
     &                  vdwmap0, phimap0, recdes0, ligdes0, gist0,
     &                  solvmap0, rec_des_r0, match,
     &                  conf_scores(1, tempconf, flexgrid), 
     &                  conf_scores(2, tempconf, flexgrid), 
     &                  conf_scores(3, tempconf, flexgrid), 
     &                  conf_scores(4, tempconf, flexgrid), 
     &                  conf_scores(5, tempconf, flexgrid), 
     &                  conf_scores(6, tempconf, flexgrid),
     &                  conf_scores(7, tempconf, flexgrid),
     &                  conf_scores(8, tempconf, flexgrid),
     &                  allminlimits, allmaxlimits,
     &                  sra, srb, 
     &                  maxor, maxconfs, maxatm, maxatmpos, maxtyv, 
c    &                  ligscore%numscored,gistgridindex_cur_rec,
     &                  ligscore%numscored,gistgridindex(:,tempconf,
     &                      flexgrid),
     &                  lengistgi(tempconf,flexgrid),
     &                  gistgridindex(:,
c     &                      db2lig%conf_conection(tempconf),
     &                      tempconf,
     &                      flexgrid),
     &                  lengistgi(
c     &                      db2lig%conf_conection(tempconf),
     &                      tempconf,
     &                      flexgrid))
                  endif             
                  
                  if (conf_status(tempconf, flexgrid) .eq. ALLOKAY) then
c                   Algorithm 6: Pre-calculate weighted total
c                   Split lines to respect 72-col limit
                    conf_weighted_total(tempconf, flexgrid) = 
     &               conf_scores(1, tempconf, flexgrid)*options0%vscale+
     &               conf_scores(2, tempconf, flexgrid)*options0%escale+
     &               conf_scores(3, tempconf, flexgrid)*
     &                                          options0%gistscale+
     &               conf_scores(5, tempconf, flexgrid)*
     &                                          options0%solvscale+
     &               conf_scores(4, tempconf, flexgrid)*
     &                                          options0%solvscale
                    conf_weighted_total(tempconf, flexgrid) = 
     &                 conf_weighted_total(tempconf, flexgrid) +
     &                 conf_scores(6, tempconf, flexgrid)*
     &                                          options0%rdscale +
     &                 conf_scores(7, tempconf, flexgrid)* 
     &                                          options0%rec_des_scale +
     &                 conf_scores(8, tempconf, flexgrid)*
     &                                          options0%rhscale
                  
                  else if (conf_status(tempconf, flexgrid) .ne. 
     &                 ALLOKAY) then
                    set_status = conf_status(tempconf, flexgrid)
c                   Mark conf as bumped for O(1) cache lookup
                    if (conf_status(tempconf,flexgrid).eq.BUMPED) then
                      is_bumped(tempconf, flexgrid) = .true.
                    endif
                    exit !break out of loop for this set
                  endif
c                 Optimized inner loop: accumulate only the weighted total
                  set_weighted_total(tempset, flexgrid) = 
     &                 set_weighted_total(tempset, flexgrid) + 
     &                 conf_weighted_total(tempconf, flexgrid)
                enddo
              else !above strain
                set_status = STRAIN
                dock_status = STRAIN
              endif      
            else !broken conf
              set_status = CLASHES !nothing else needs checked, no need to sum
              dock_status = CLASHES ! this is so that if everything clashes then it reports as a clash.
            endif
            goodset(tempset, flexgrid) = .false.
            if (set_status .eq. ALLOKAY) then !each conf scored
              goodset(tempset, flexgrid) = .true.
              setscore(tempset, flexgrid) = 
     &            set_weighted_total(tempset, flexgrid) + 
     &            db2lig%set_energy(tempset) + 
     &            options0%rec_energy(flexgrid)
            endif
          enddo      
        enddo
        do temp_combo_num = 1, grids0%total_combinations
          temp_combo = grids0%all_combos(temp_combo_num, :)
          reccode = getcode(grids0, temp_combo_num)        
          do tempset = setstart, setend !going through every set as well
            goodwhole = .true.
            wholevs(tempset, temp_combo_num) = 0.
            wholees(tempset, temp_combo_num) = 0.
            wholegi(tempset, temp_combo_num) = 0.
            wholeas(tempset, temp_combo_num) = 0.
            wholeps(tempset, temp_combo_num) = 0.
            wholeds(tempset, temp_combo_num) = 0.
            wholerd(tempset, temp_combo_num) = 0.
            wholehs(tempset, temp_combo_num) = 0.
            wholers(tempset, temp_combo_num) = 0.
            wholescore(tempset, temp_combo_num) = 0.
            do flexgrid = 1, grids0%group_len
              if (goodset(tempset, temp_combo(flexgrid))) then
c               Use weighted total for quick check
c               Note: We removed setvs, setes etc intermediate arrays.
c               Only accumulate wholescore here for checking.
c               Components will be reconstructed later if good.
                wholescore(tempset, temp_combo_num) = 
     &              wholescore(tempset, temp_combo_num) + 
     &              set_weighted_total(tempset, temp_combo(flexgrid)) +
     &              options0%rec_energy(temp_combo(flexgrid))
              else
                goodwhole = .false.
                exit !break out, no need to do more
              endif
            enddo
c           Add set_energy to complete the wholescore
            if (goodwhole) then
                 wholescore(tempset, temp_combo_num) = 
     &              wholescore(tempset, temp_combo_num) +
     &              db2lig%set_energy(tempset)

c             Lazy check against cutoff
c             Note: ligscore%pose_score is sorted best (lowest) to worst.
c             ligscore%savedcount is the number of saved poses.
c             We want to save if buffer is not full OR if score is better than worst.
              if ((ligscore%savedcount .lt. options0%nsav) .or. 
     &            (wholescore(tempset, temp_combo_num) .lt. 
     &             ligscore%pose_score(ligscore%savedcount))) then

c               It's a top hit! Reconstruct component scores now.
                do flexgrid = 1, grids0%group_len
                   temp_combo_num_grid = temp_combo(flexgrid)

                   do tempconfnum=1, db2lig%set_conf_len(tempset)
                      tempconf = db2lig%set_conf(tempset, tempconfnum)
                      
                      wholevs(tempset, temp_combo_num) = 
     &                   wholevs(tempset, temp_combo_num) + 
     &                   conf_scores(1, tempconf, temp_combo_num_grid) * 
     &                   options0%vscale

                      wholees(tempset, temp_combo_num) = 
     &                   wholees(tempset, temp_combo_num) + 
     &                   conf_scores(2, tempconf, temp_combo_num_grid) * 
     &                   options0%escale

                      wholegi(tempset, temp_combo_num) =
     &                   wholegi(tempset, temp_combo_num) +
     &                   conf_scores(3, tempconf, temp_combo_num_grid) * 
     &                   options0%gistscale

                      wholeas(tempset, temp_combo_num) = 
     &                   wholeas(tempset, temp_combo_num) +
     &                   conf_scores(4, tempconf, temp_combo_num_grid) * 
     &                   options0%solvscale

                      wholeps(tempset, temp_combo_num) = 
     &                   wholeps(tempset, temp_combo_num) + 
     &                   conf_scores(5, tempconf, temp_combo_num_grid) * 
     &                   options0%solvscale

                      wholeds(tempset, temp_combo_num) = 
     &                   wholeds(tempset, temp_combo_num) + 
     &                   conf_scores(6, tempconf, temp_combo_num_grid) * 
     &                   options0%rdscale

                      wholerd(tempset, temp_combo_num) = 
     &                   wholerd(tempset, temp_combo_num) + 
     &                   conf_scores(7, tempconf, temp_combo_num_grid) * 
     &                   options0%rec_des_scale

                      wholehs(tempset, temp_combo_num) = 
     &                   wholehs(tempset, temp_combo_num) + 
     &                   conf_scores(8, tempconf, temp_combo_num_grid) * 
     &                   options0%rhscale
                   enddo
                   
                   ! Add receptor energy for this grid
                   wholers(tempset, temp_combo_num) = 
     &                wholers(tempset, temp_combo_num) + 
     &                options0%rec_energy(temp_combo_num_grid)

                enddo 
                
              else
                 goodwhole = .false. ! Mark as bad to skip saving block below
              endif
            endif ! end lazy check block

            if (goodwhole) then !don't save if not good
              onegood = .true. ! there is at least one that is good. 

              !all logic about whether or not to save (if score good enough,
              ! etc.) is in the functions called here

              !if dockovalent supply the ligands coords to save
               if (options0%dockovalent .eqv. .true.) then
                  do tempconfnum = 1, db2lig%set_conf_len(tempset)
                     tempconf = db2lig%set_conf(tempset, tempconfnum)
                     do tempcoord = db2lig%conf_coord(1, tempconf), 
     &                    db2lig%conf_coord(2, tempconf)
                        tempatom = db2lig%coord_index(2, tempcoord)
                        do count = 1, 3
                           db2lig%transfm_coords(count, tempatom) = 
     &                          db2lig%coords(count,tempcoord)
                        enddo
                     enddo
                  enddo
               endif
              call save_set_score(ligscore, current_match, tempset, 
     &            wholevs(tempset, temp_combo_num), 
     &            wholees(tempset, temp_combo_num), 
     &            wholegi(tempset, temp_combo_num), 
     &            wholeas(tempset, temp_combo_num), 
     &            wholeps(tempset, temp_combo_num), 
c     &            db2lig%set_energy(tempset), 
     &            db2lig%set_total_strain(tempset), 
c     &            wholers(tempset, temp_combo_num), 
     &            db2lig%set_max_strain(tempset), 
     &            wholeds(tempset, temp_combo_num), 
     &            wholerd(tempset, temp_combo_num), 
     &            wholehs(tempset, temp_combo_num), 
     &            wholescore(tempset, temp_combo_num), 
     &            reccode, temp_combo_num, options0%nsav,
     &            options0%dockovalent,db2lig%save_cov_coords,
     &            db2lig%total_atoms,db2lig%transfm_coords)
              if (options0%flexible_receptor .and.
     &            options0%score_each_flex) then
                call save_set_score_each_flex(ligscoreeach, 
     &              current_match, tempset, 
     &              wholevs(tempset, temp_combo_num), 
     &              wholees(tempset, temp_combo_num), 
     &              wholegi(tempset, temp_combo_num), 
     &              wholeas(tempset, temp_combo_num), 
     &              wholeps(tempset, temp_combo_num), 
c     &              db2lig%set_energy(tempset), 
     &              db2lig%set_total_strain(tempset), 
c     &              wholers(tempset, temp_combo_num), 
     &              db2lig%set_max_strain(tempset), 
     &              wholeds(tempset, temp_combo_num), 
     &              wholerd(tempset, temp_combo_num), 
     &              wholehs(tempset, temp_combo_num), 
     &              wholescore(tempset, temp_combo_num), 
     &              reccode, temp_combo_num)
              endif
            endif
          enddo
        enddo
        if (onegood) then
           dock_status = ALLOKAY !at least one conformation so this run is good
        endif
      endif
      return
      end subroutine calc_score_sets

      end module score_sets

