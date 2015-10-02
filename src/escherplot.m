#! /usr/bin/octave -q

args = argv();
source(args{1});
mol = cr_crystalbox(cr);
rep = mol_ball(mol);
mol = mol_getfragment(mol,find(mol.atnumber < 104));
rep = mol_stick(mol,rep);
rep = cr_unitcell(cr,rep);
rep_write_obj(rep,strcat(substr(args{1},1,index(args{1},".m")-1),".obj"));
