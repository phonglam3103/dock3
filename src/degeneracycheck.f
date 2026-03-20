c makes a hash of the lig spheres and rec spheres. checks this against
c previous hashes, if new (not seen before) add to list and return true. 
c if seen before, return false.
c OPTIMIZED: Uses hash buckets for O(1) average lookup instead of O(n) linear scan

      logical function degeneracycheck(nodetotal, sphpairs, 
     &          hash, hashcount,
     &          maxpts, maxctr, maxor)

      implicit none

c XXX: GFortran needs these to come first (TS)
c bunch of maximum parameters
      integer maxpts !maximum number of sphere points
      integer maxctr !maximum number of sphere centers
      integer maxor !maximum number of orientations

c first are all the input parameters
      integer nodetotal !match sphere parameters, how many sphere pairs to use
      integer sphpairs(maxctr, 2) !match sphere parameters
      integer hash(maxor, maxctr) !hash of sph numbers
      integer hashcount !how many hashes have been stored

c hash bucket parameters - using prime number for bucket count
      integer, parameter :: NBUCKETS = 65537
      integer, save :: bucket_heads(NBUCKETS)
      integer, save :: bucket_next(100000)
      logical, save :: buckets_initialized = .false.
      integer, save :: last_hashcount = 0

c temporary variables
      integer nodeindex, hashindex, centerindex, tmpindex !temporary counters
      integer tmphash(maxctr) !the temporarily computed hash
      logical match_found, all_matched
      integer bucket_idx, hash_val, idx

c initialize or reset buckets if hashcount was reset (new ligand)
      if (.not. buckets_initialized .or. hashcount .lt. last_hashcount) 
     &    then
        do idx = 1, NBUCKETS
          bucket_heads(idx) = 0
        enddo
        buckets_initialized = .true.
        last_hashcount = 0
      endif

      degeneracycheck = .true. !default is true
      
c compute hash for current match
      hash_val = 0
      do nodeindex = 1, nodetotal !for every node
        tmphash(nodeindex) = sphpairs(nodeindex, 1) + 
     &                       sphpairs(nodeindex, 2) * maxpts !makes uniq number
        hash_val = hash_val + tmphash(nodeindex) * (7919 ** nodeindex)
      enddo
      
c compute bucket index (1 to NBUCKETS)
      bucket_idx = mod(abs(hash_val), NBUCKETS) + 1
      
c search only within this bucket instead of all entries
      idx = bucket_heads(bucket_idx)
      do while (idx .gt. 0)
        !check if this entry matches - assume order doesn't matter
        all_matched = .true.
        do tmpindex = 1, maxctr
          match_found = .false.
          do centerindex = 1, maxctr
            if (tmphash(tmpindex) .eq. hash(idx, centerindex)) then
              match_found = .true.
              exit
            endif
          enddo
          if (.not. match_found) then
            all_matched = .false.
            exit !break out of the do loop
          endif       
        enddo       
        if (all_matched) then !whole tmphash matched
          degeneracycheck = .false.
          return !quit and return, match found, nothing to do
        endif
        idx = bucket_next(idx)
      enddo
      
c not found - add new entry
      hashcount = hashcount + 1
      last_hashcount = hashcount
      do centerindex = 1, maxctr !copy them all, even 0s
        hash(hashcount, centerindex) = tmphash(centerindex)
      enddo
c insert at head of bucket's linked list
      bucket_next(hashcount) = bucket_heads(bucket_idx)
      bucket_heads(bucket_idx) = hashcount

      return !degeneracycheck as the return value, is still true if here
      end
