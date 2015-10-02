#! /usr/bin/awk -f
#
# basinmerge - Script to automatize the plot of combined basins in a
#              crystal.
#
# Author: Victor Lua~na (1998.07.17)
#
# Uses: tessel, off2off, basin2off
#
# The program interprets a collection of orders/instructions written
# according to the next language:
#
# Ins. 01:   FILETYPE [ OFF | COFF | BASIN | DBASIN ]
#
# Ins. 02:   NEQION  name.s  file.s  x.r y.r z.r
#    Enter basin description for each non-equivalent ion
#
# Ins. 03:   CRYS2CART
#              c11.r  c12.r  c13.r  c14.r
#              c21.r  c22.r  c23.r  c24.r
#              c31.r  c32.r  c33.r  c34.r
#              c41.r  c42.r  c43.r  c44.r
#    Matrix that transforms from crystal to cartesian coordinates.
#    It is written in the critic output.
#
# Ins. 04:   CART2CRYST
#              d11.r  d12.r  d13.r  d14.r
#              d21.r  d22.r  d23.r  d24.r
#              d31.r  d32.r  d33.r  d34.r
#              d41.r  d42.r  d43.r  d44.r
#    Matrix that transforms from cartesian to crystal coordinates.
#    It is written in the critic output.
#
# Ins. 05:   COPY  fromname.s  toname.s  fileto.s
#              t11.r  t12.r  t13.r  t14.r
#              t21.r  t22.r  t23.r  t24.r
#              t31.r  t32.r  t33.r  t34.r
#              t41.r  t42.r  t43.r  t44.r
#    Create a copy of the basin file of a NEQ ion into an equivalent
#    position. The transformation matrix that is entered should act
#    on the crystallographic coordinates, and it is produced by
#    critic.
#
# Ins. 06:   MERGETO  destfile.s  file1.s  [ file2.s ... ]
#    Merge a collection of files into a single one.
#
# Ins. 07:   BASIN2OFF  order.s
#    Run directly the basin2off program.
#
# Ins. 08:   OFF2OFF  order.s
#    Run directly the off2off program.
#

BEGIN {
   filetypes["coff"]   = 1
   filetypes["off"]    = 1
   filetypes["basin"]  = 1
   filetypes["dbasin"] = 1
   pid = "123"
   }

{ print "DBG(basinmerge) Order -> ", $0 }

tolower($1) == "filetype" {
   if (tolower($2) in filetypes) { type = tolower($2) }
   else { print "Wrong file type"; exit(1) }
   }

tolower($1) == "neqion" {
   neq++
   neq_name[neq] = $2
   neq_file[neq] = $3
   neq_x[neq] = $4
   neq_y[neq] = $5
   neq_z[neq] = $6
   name_neq[$2] = neq
   }

tolower($1) == "crys2cart" {
   getline; cy2ca[1,1] = $1; cy2ca[1,2] = $2; cy2ca[1,3] = $3; cy2ca[1,4] = $4
   getline; cy2ca[2,1] = $1; cy2ca[2,2] = $2; cy2ca[2,3] = $3; cy2ca[2,4] = $4
   getline; cy2ca[3,1] = $1; cy2ca[3,2] = $2; cy2ca[3,3] = $3; cy2ca[3,4] = $4
   getline; cy2ca[4,1] = $1; cy2ca[4,2] = $2; cy2ca[4,3] = $3; cy2ca[4,4] = $4
   }

tolower($1) == "cart2crys" {
   getline; ca2cy[1,1] = $1; ca2cy[1,2] = $2; ca2cy[1,3] = $3; ca2cy[1,4] = $4
   getline; ca2cy[2,1] = $1; ca2cy[2,2] = $2; ca2cy[2,3] = $3; ca2cy[2,4] = $4
   getline; ca2cy[3,1] = $1; ca2cy[3,2] = $2; ca2cy[3,3] = $3; ca2cy[3,4] = $4
   getline; ca2cy[4,1] = $1; ca2cy[4,2] = $2; ca2cy[4,3] = $3; ca2cy[4,4] = $4
   }

tolower($1) == "copy" {
   from = $2
   to = $3
   tofile = $4
   if (from in name_neq) {
      ifrom = name_neq[from]
      fromfile = neq_file[ifrom]
      }
   else {
      print "Wrong ion to copy!"
      print "Source ions known:"
      printf "->"
      for (i=1; i<=neq; i++) {
         printf " %s", neq_name[i]
         }
      printf "\n"
      exit (1)
      }
   getline; tmat[1,1] = $1; tmat[1,2] = $2; tmat[1,3] = $3; tmat[1,4] = $4
   getline; tmat[2,1] = $1; tmat[2,2] = $2; tmat[2,3] = $3; tmat[2,4] = $4
   getline; tmat[3,1] = $1; tmat[3,2] = $2; tmat[3,3] = $3; tmat[3,4] = $4
   getline; tmat[4,1] = $1; tmat[4,2] = $2; tmat[4,3] = $3; tmat[4,4] = $4
   close("tmp"pid)
   print ca2cy[1,1], ca2cy[1,2], ca2cy[1,3], ca2cy[1,4]  >  "tmp"pid
   print ca2cy[2,1], ca2cy[2,2], ca2cy[2,3], ca2cy[2,4]  >  "tmp"pid
   print ca2cy[3,1], ca2cy[3,2], ca2cy[3,3], ca2cy[3,4]  >  "tmp"pid
   print ca2cy[4,1], ca2cy[4,2], ca2cy[4,3], ca2cy[4,4]  >  "tmp"pid
   print tmat[1,1], tmat[1,2], tmat[1,3], tmat[1,4]  >  "tmp"pid
   print tmat[2,1], tmat[2,2], tmat[2,3], tmat[2,4]  >  "tmp"pid
   print tmat[3,1], tmat[3,2], tmat[3,3], tmat[3,4]  >  "tmp"pid
   print tmat[4,1], tmat[4,2], tmat[4,3], tmat[4,4]  >  "tmp"pid
   print cy2ca[1,1], cy2ca[1,2], cy2ca[1,3], cy2ca[1,4]  >  "tmp"pid
   print cy2ca[2,1], cy2ca[2,2], cy2ca[2,3], cy2ca[2,4]  >  "tmp"pid
   print cy2ca[3,1], cy2ca[3,2], cy2ca[3,3], cy2ca[3,4]  >  "tmp"pid
   print cy2ca[4,1], cy2ca[4,2], cy2ca[4,3], cy2ca[4,4]  >  "tmp"pid
   if ((type == "off")  ||  (type == "coff")) {
      order = "off2off -c -t tmp" pid " -o " tofile " " fromfile
      }
   else if (type == "basin") {
      order = "off2off -b -t tmp" pid " -o " tofile " " fromfile
      }
   else if (type == "dbasin") {
      order = "off2off -d -t tmp" pid " -o " tofile " " fromfile
      }
   else { print "Wrong file type to copy!"; exit(1) }
   print order
   system (order)
   if (to in name_neq) { }
   else {
      neq++
      neq_name[neq] = to
      neq_file[neq] = tofile
      neq_x[neq] = -1
      neq_y[neq] = -1
      neq_z[neq] = -1
      name_neq[to] = neq
      }
   }

tolower($1) == "mergeto" {
   tofile = $2
   $1 = ""
   $2 = ""
   fromfile = $0
   if ((type == "off")  ||  (type == "coff")) {
      order = "off2off -o " tofile " " fromfile
      }
   else if (type == "basin") {
      order = "off2off -b -o " tofile " " fromfile
      }
   else if (type == "dbasin") {
      order = "off2off -d -o " tofile " " fromfile
      }
   else { print "Wrong file type to copy!"; exit(1) }
   print order
   system (order)
   }

tolower($1) == "basin2off" {
   print
   system ($0)
   }

tolower($1) == "off2off" {
   print
   system ($0)
   }
