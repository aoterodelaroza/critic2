#! /usr/bin/octave -q

## Requires modifying ld1.x to:
## - Increase ndmx in Modules/radial_grids.f90
## - Increase the maximum Z in ld1_readin.f90 (search for "wrong nuclear charge")
## - Increase the precision of the wfc (E28.14), in atomic/src/write_files.f90
## - Increase the number of known atomic species to 118 (Modules/atomic_number.f90)

atoms={"pu","am","cm","bk","cf","es","fm","md","no","lr","rf","db","sg",...
       "bh","hs","mt","ds","rg","cn","nh","fl","mc","lv","ts","og"}; 
atnums=[94  95  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118];

confs={...
        "[Rn] 5f6 7s2",...
        "[Rn] 5f7 7s2",...
        "[Rn] 5f7 6d1 7s2",...
        "[Rn] 5f9 7s2",...
        "[Rn] 5f10 7s2",...
        "[Rn] 5f11 7s2",...
        "[Rn] 5f12 7s2",...
        "[Rn] 5f13 7s2",...
        "[Rn] 5f14 7s2",...
        "[Rn] 5f14 6d1 7s2",...
        "[Rn] 5f14 6d2 7s2",...
        "[Rn] 5f14 6d3 7s2",...
        "[Rn] 5f14 6d4 7s2",...
        "[Rn] 5f14 6d5 7s2",...
        "[Rn] 5f14 6d6 7s2",...
        "[Rn] 5f14 6d7 7s2",...
        "[Rn] 5f14 6d9 7s1",...
        "[Rn] 5f14 6d10 7s1",...
        "[Rn] 5f14 6d10 7s2",...
        "[Rn] 5f14 6d10 7s2 7p1",...
        "[Rn] 5f14 6d10 7s2 7p2",...
        "[Rn] 5f14 6d10 7s2 7p3",...
        "[Rn] 5f14 6d10 7s2 7p4",...
        "[Rn] 5f14 6d10 7s2 7p5",...
        "[Rn] 5f14 6d10 7s2 7p6",...
      };

aocc = struct();
ll = {"S","P","D","F"};
loc = [2 6 10 14];
for i = 1:8
  for j = 1:length(ll)
    aocc = setfield(aocc,sprintf("%d%s",i,ll{j}),loc(j));
  endfor
endfor

for i = 1:length(atoms)
  tocc = aocc;
  str = confs{1};
  aa = strfields(str);
  for j = 1:length(aa)
    if (aa{j}(1:1) == "[")
      continue
    endif
    n = str2num(aa{j}(1:1));
    chan = upper(aa{j}(2:2));
    occ = str2num(aa{j}(3:end));
    tocc = setfield(tocc,sprintf("%d%s",n,chan),occ);
  endfor

  fid = fopen(sprintf("%s.in",atoms{i}),"w");
  fprintf(fid,"&input\n");
  fprintf(fid,"  title='%s',\n",atoms{i});
  fprintf(fid,"  zed=%d.,\n",atnums(i));
  fprintf(fid,"  rel=1,\n");
  fprintf(fid,"  config='%s',\n",confs{i});
  fprintf(fid,"  iswitch=1,\n");
  fprintf(fid,"  dft='PBE'\n");
  fprintf(fid,"  max_out_wfc=99,\n");
  fprintf(fid,"  xmin=-6,\n");
  fprintf(fid,"  dx=0.002,\n");
  fprintf(fid,"  rmax=100.,\n");
  fprintf(fid,"/\n");
  fclose(fid);
  system(sprintf("~/src/espresso-6.1/bin/ld1.x < %s.in",atoms{i}));
  
  fid = fopen("ld1.wfc","r");
  line = fgetl(fid);
  aa = strfields(line);
  wfn = zeros(9000,length(aa)-1);

  n = 0;
  while (!feof(fid))
    n++;
    line = fgetl(fid);
    wfn(n,:) = str2num(line);
  endwhile
  fclose(fid);

  nwfn = length(aa) - 2;
  fid = fopen(sprintf("%s_pbe.wfc",atoms{i}),"w");
  fprintf(fid,"%d\n",nwfn);
  for j = nwfn+2:-1:3
    fprintf(fid,"%s ",aa{j});
  endfor
  fprintf(fid,"\n");
  for j = nwfn+2:-1:3
    occ = getfield(tocc,upper(aa{j}));
    fprintf(fid,"%d ",occ);
  endfor
  fprintf(fid,"\n");
  fprintf(fid,"%.10f %.10f %.10f %d\n",-6,atnums(i),0.002,n);
  for j = 1:n
    fprintf(fid,"%.14e ",wfn(j,1));
    for k = nwfn+1:-1:2
      fprintf(fid,"%.14e ",wfn(j,k));
    endfor
    fprintf(fid,"\n");
  endfor
  fclose(fid);
endfor
