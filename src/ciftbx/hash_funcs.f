C
C       hash_funcs.f -- a library of hash table management routines
C
C                                      by
C
C                              Herbert J. Bernstein
C                                Bernstein + Sons
C                    5 Brewster Lane, Bellport, NY 11713-0177, USA
C                       email: yaya@bernstein-plus-sons.com
C
C       work on these routines done in part at Brookhaven National
C       Laboratory, under contract to the U.S. Department of Energy
C
C       work on these routines as been done in part at Dowling College 
C       under contract to the International Union of Crystallography 
C       and under grants from the National Science Foundation and 
C       the U.S. Department of Energy.
C
C       Copyright (C) Herbert J. Bernstein 1997 -- 2006
C
C       hash_funcs.f is free software; you can redistribute this software 
C       and/or modify this software under the terms of the GNU General 
C       Public License as published by the Free Software Foundation; 
C       either version 2 of the License, or (at your option) any later version.
C
C       Alternatively you may reditribute and/or modify hash_funcs.f
C       under the terms of the GNU Lesser General Public 
C       License as published by the Free Software Foundation; either 
C       version 2.1 of the License, or (at your option) any later version.
C
C       This software is distributed in the hope that it will be useful,
C       but WITHOUT ANY WARRANTY; without even the implied warranty of
C       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C       GNU General Public License for more details.
C
C       You should have received a copy of the GNU General Public License
C       along with this software; if not, write to the Free Software
C       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C       You should have received a copy of the GNU Lesser General Public License
C       along with this software; if not, write to the Free Software
C       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C-------------------------------------------------------------------------------
C
C       Routines
C
C       hash_init          Initializes a hash table controlled list
C                          call hash_init(data_structure_args)
C
C       hash_find          Searches for a string in a list
C                          call hash_find(name,data_structure_args,ifind)
C
C       hash_fnext         Searches for the next matching string in a list
C                          call hash_fnext(name,data_structure_args,ifind,icurr)
C
C       hash_store         Inserts a new string in a list
C                          call hash_store(name,data_structure_args,ifind)
C
C       hash_snext         Inserts a string in a list, allowing duplicates
C                          call hash_store(name,data_structure_args,ifind,icurr)
C
C       hash_value         Integer function returns index into hash_list
C                          ih = hash_value(name,hash_length)
C
C       The necessary data_structure_args for these routines are
C          name_list   -- an array of character strings
C                         character*(*) name_list(list_length)
C          chain_list  -- chain pointers for searches
C                         integer chain_list(list_length)
C          list_length -- the size of the list arrays
C                         integer list_length
C          num_list    -- number of entries in the list
C                         integer num_list
C          hash_table  -- the initial hashed pointers
C                         integer hash_table
C          hash_length -- the size of the hash table
C                         integer hash_length
C
C
C       The three remaining arguments are
C          name        -- string to search for
C                         character*(*) name
C          ifind       -- return value, 0 for not found (hash_find)
C                         or list full (hash_store), otherwise
C                         the index in name_list of the entry
C          icurr       -- the prior matching index
C
C       The relationship among the arrays used is:
C
C       hash_table is an array (preferably of a modest prime
C       dimension) which starts containing all zeros, which are
C       replaced by pointers to entries in name_list, based
C       values returned by hash_value ranging from 1 to hash_length.
C       Each name is placed in name_list.  A initial zero is placed
C       in the matching entry in chain_list, when the first entry
C       is made.  When a new entry with the same hash_value must be
C       placed a pointer is inserted into chain_list to hook the
C       values together.
C
        subroutine hash_init(name_list,chain_list,list_length,num_list,
     *                       hash_table,hash_length)
C
C       initialization routine for a hash table controlled list
C          name_list   -- a list of character strings
C          chain_list  -- chain pointers for searches
C          list_length -- the size of the list arrays
C          num_list    -- number of entries in the list
C          hash_table  -- the initial hashed pointers
C          hash_length -- the size of the hash table
C
           integer list_length
           character*(*) name_list(list_length)
           integer hash_length,num_list,i
           integer chain_list(list_length)
           integer hash_table(hash_length)
           num_list=0
           do i = 1,hash_length
           hash_table(i)=0
           enddo
           return
           end          

        subroutine
     *  hash_find(name,name_list,chain_list,list_length,num_list,
     *                       hash_table,hash_length,ifind)
C
C       search routine for a hash table controlled list
C          name        -- string to find
C          name_list   -- a list of character strings
C          chain_list  -- chain pointers for searches
C          list_length -- the size of the list arrays
C          num_list    -- number of entries in the list
C          hash_table  -- the initial hashed pointers
C          hash_length -- the size of the hash table
C          ifind       -- returned index or 0
C
           character*(*) name
           integer list_length, hash_length
           integer lenn
           character*(*) name_list(list_length)
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,num_list,ih,ip
           integer lastnb
           ifind=0
           ih=hash_value(name,hash_length)
           ip=hash_table(ih)
           lenn = lastnb(name)
 100       if (ip.eq.0) return
           if (name_list(ip).eq.name(1:lenn)) then
             ifind=ip
             return
           else
             ip=chain_list(ip)
             go to 100
           endif
           end
           
        subroutine
     *  hash_fnext(name,name_list,chain_list,list_length,num_list,
     *                       hash_table,hash_length,ifind,icurr)
C
C       search routine for a hash table controlled list
C          name        -- string to find
C          name_list   -- a list of character strings
C          chain_list  -- chain pointers for searches
C          list_length -- the size of the list arrays
C          num_list    -- number of entries in the list
C          hash_table  -- the initial hashed pointers
C          hash_length -- the size of the hash table
C          ifind       -- returned index or 0
C          icurr       -- current match or 0
C
           character*(*) name
           integer hash_length
           integer list_length
           integer lenn
           character*(*) name_list(list_length)
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,num_list,ih,ip
           integer lastnb
           integer icurr
           
           ifind=0
           if (icurr.eq.0) then
             ih = hash_value(name,hash_length)
             ip = hash_table(ih)
           else
             ip = chain_list(icurr)
           endif
           lenn = lastnb(name)
 100       if (ip.eq.0) return
           if (name_list(ip).eq.name(1:lenn)) then
             ifind=ip
             return
           else
             ip=chain_list(ip)
             go to 100
           endif
           end

        subroutine
     *  hash_store(name,name_list,chain_list,list_length,num_list,
     *                       hash_table,hash_length,ifind)
C
C       store routine for a hash table controlled list
C          name        -- string to find
C          name_list   -- a list of character strings
C          chain_list  -- chain pointers for searches
C          list_length -- the size of the list arrays
C          num_list    -- number of entries in list
C          hash_table  -- the initial hashed pointers
C          hash_length -- the size of the hash table
C          ifind       -- index of entry or 0 (table full)
C
           integer list_length
           character*(*) name
           character*(*) name_list(list_length)
           integer hash_length
           integer lenn
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,num_list,ih,ip,iq
           integer lastnb

           ifind=0
           ih = hash_value(name,hash_length)
           ip=hash_table(ih)
           iq=0
           lenn = lastnb(name)
 100       if (ip.eq.0) go to 200
           if (name_list(ip).eq.name(1:lenn)) then
             ifind=ip
             return
           else
             iq=ip
             ip=chain_list(ip)
             go to 100
           endif
 200       if (num_list.lt.list_length) then
             num_list=num_list+1
             name_list(num_list)=name(1:lenn)
             chain_list(num_list)=0
             if (iq.eq.0) then
               hash_table(ih)=num_list
             else
               chain_list(iq)=num_list
             endif
             ifind=num_list
             return
           else
             ifind = 0
             return
           endif
           end


        subroutine
     *  hash_snext(name,name_list,chain_list,list_length,num_list,
     *                       hash_table,hash_length,ifind,icurr)
C
C       store routine for a hash table controlled list
C          name        -- string to find
C          name_list   -- a list of character strings
C          chain_list  -- chain pointers for searches
C          list_length -- the size of the list arrays
C          num_list    -- number of entries in list
C          hash_table  -- the initial hashed pointers
C          hash_length -- the size of the hash table
C          ifind       -- index of entry or 0 (table full)
C          icurr       -- current match or 0
C
           integer list_length
           character*(*) name
           character*(*) name_list(list_length)
           integer hash_length
           integer lenn
           integer chain_list(list_length)
           integer hash_table(hash_length)
           integer hash_value
           integer ifind,num_list,ih,ip,iq
           integer lastnb
           integer icurr

           ifind=0
           ih = 0
           if (icurr.eq.0) then
             ih = hash_value(name,hash_length)
             ip=hash_table(ih)
             iq = 0
           else
             ip = chain_list(icurr)
             iq = icurr
           endif
           lenn = lastnb(name)
 100       if (ip.eq.0) go to 200
           if (name_list(ip).eq.name(1:lenn)) then
             ifind=ip
             return
           else
             iq=ip
             ip=chain_list(ip)
             go to 100
           endif
 200       if (num_list.lt.list_length) then
             num_list=num_list+1
             name_list(num_list)=name(1:lenn)
             chain_list(num_list)=0
             if (iq.eq.0) then
               hash_table(ih)=num_list
             else
               chain_list(iq)=num_list
             endif
             ifind=num_list
             return
           else
             ifind = 0
             return
           endif
           end

      integer function hash_value(name,hash_length)
C
C     function to return a hash value of string name to fit
C     a hash table of length hash_length
      character*(*) name
      integer lastnb
      
      integer hash_length,id,ii,i,ic,lenn
      lenn = lastnb(name)
      hash_value=1
      id = 0
      do ii = 1,lenn
        i = 1+lenn-ii
        ic = ichar(name(i:i))
        if (ic.ge.65) then
          hash_value=mod(hash_value*(ic-64),hash_length)+1
          id = id+1
          if (id.gt.3) return
        endif
      enddo
      return
      end
        
