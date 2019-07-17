#! /usr/bin/octave -q

## 1. Requires modifying ld1.x to:
## - Increase ndmx in Modules/radial_grids.f90
## - Increase the maximum Z in ld1_readin.f90 (search for "wrong nuclear charge")
## - Increase the precision of the wfc (E28.14), in atomic/src/write_files.f90
## - Increase the number of known atomic species to 118 (Modules/atomic_number.f90)

atoms={...
        "h" ,"he","li","be","b" ,"c" ,"n" ,"o" ,"f" ,"ne",... ## 1-10
        "na","mg","al","si","p" ,"s" ,"cl","ar","k" ,"ca",... ## 11-20
        "sc","ti","v" ,"cr","mn","fe","co","ni","cu","zn",... ## 21-30
        "ga","ge","as","se","br","kr","rb","sr","y" ,"zr",... ## 31-40
        "nb","mo","tc","ru","rh","pd","ag","cd","in","sn",... ## 41-50
        "sb","te","i" ,"xe","cs","ba","la","ce","pr","nd",... ## 51-60
        "pm","sm","eu","gd","tb","dy","ho","er","tm","yb",... ## 61-70
        "lu","hf","ta","w" ,"re","os","ir","pt","au","hg",... ## 71-80
        "tl","pb","bi","po","at","rn","fr","ra","ac","th",... ## 81-90
        "pa","u" ,"np","pu","am","cm","bk","cf","es","fm",... ## 91-100
        "md","no","lr","rf","db","sg","bh","hs","mt","ds",... ## 101-110
        "rg","cn","nh","fl","mc","lv","ts","og",...           ## 111-108
      };
atnums=1:118;
confs={... ## from wikipedia - see notes on nickel and heavy elements
        "1s1",...                        ## h  (1) 
        "1s2",...                        ## he (2) 
        "[He] 2s1",...                   ## li (3) 
        "[He] 2s2",...                   ## be (4) 
        "[He] 2s2 2p1",...               ## b  (5) 
        "[He] 2s2 2p2",...               ## c  (6) 
        "[He] 2s2 2p3",...               ## n  (7) 
        "[He] 2s2 2p4",...               ## o  (8) 
        "[He] 2s2 2p5",...               ## f  (9) 
        "[He] 2s2 2p6",...               ## ne (10) 
        "[Ne] 3s1",...                   ## na (11) 
        "[Ne] 3s2",...                   ## mg (12) 
        "[Ne] 3s2 3p1",...               ## al (13) 
        "[Ne] 3s2 3p2",...               ## si (14) 
        "[Ne] 3s2 3p3",...               ## p  (15) 
        "[Ne] 3s2 3p4",...               ## s  (16) 
        "[Ne] 3s2 3p5",...               ## cl (17) 
        "[Ne] 3s2 3p6",...               ## ar (18) 
        "[Ar] 4s1",...                   ## k  (19) 
        "[Ar] 4s2",...                   ## ca (20) 
        "[Ar] 3d1 4s2",...               ## sc (21) 
        "[Ar] 3d2 4s2",...               ## ti (22) 
        "[Ar] 3d3 4s2",...               ## v  (23) 
        "[Ar] 3d5 4s1",...               ## cr (24) 
        "[Ar] 3d5 4s2",...               ## mn (25) 
        "[Ar] 3d6 4s2",...               ## fe (26) 
        "[Ar] 3d7 4s2",...               ## co (27) 
        "[Ar] 3d8 4s2",...               ## ni (28) 
        "[Ar] 3d10 4s1",...              ## cu (29) 
        "[Ar] 3d10 4s2",...              ## zn (30) 
        "[Ar] 3d10 4s2 4p1",...          ## ga (31) 
        "[Ar] 3d10 4s2 4p2",...          ## ge (32) 
        "[Ar] 3d10 4s2 4p3",...          ## as (33) 
        "[Ar] 3d10 4s2 4p4",...          ## se (34) 
        "[Ar] 3d10 4s2 4p5",...          ## br (35) 
        "[Ar] 3d10 4s2 4p6",...          ## kr (36) 
        "[Kr] 5s1",...                   ## rb (37) 
        "[Kr] 5s2",...                   ## sr (38) 
        "[Kr] 4d1 5s2",...               ## y  (39) 
        "[Kr] 4d2 5s2",...               ## zr (40) 
        "[Kr] 4d4 5s1",...               ## nb (41) 
        "[Kr] 4d5 5s1",...               ## mo (42) 
        "[Kr] 4d5 5s2",...               ## tc (43) 
        "[Kr] 4d7 5s1",...               ## ru (44) 
        "[Kr] 4d8 5s1",...               ## rh (45) 
        "[Kr] 4d10",...                  ## pd (46) 
        "[Kr] 4d10 5s1",...              ## ag (47) 
        "[Kr] 4d10 5s2",...              ## cd (48) 
        "[Kr] 4d10 5s2 5p1",...          ## in (49) 
        "[Kr] 4d10 5s2 5p2",...          ## sn (50) 
        "[Kr] 4d10 5s2 5p3",...          ## sb (51) 
        "[Kr] 4d10 5s2 5p4",...          ## te (52) 
        "[Kr] 4d10 5s2 5p5",...          ## i  (53) 
        "[Kr] 4d10 5s2 5p6",...          ## xe (54) 
        "[Xe] 6s1",...                   ## cs (55) 
        "[Xe] 6s2",...                   ## ba (56) 
        "[Xe] 5d1 6s2",...               ## la (57) 
        "[Xe] 4f1 5d1 6s2",...           ## ce (58) 
        "[Xe] 4f3 6s2",...               ## pr (59) 
        "[Xe] 4f4 6s2",...               ## nd (60) 
        "[Xe] 4f5 6s2",...               ## pm (61) 
        "[Xe] 4f6 6s2",...               ## sm (62) 
        "[Xe] 4f7 6s2",...               ## eu (63) 
        "[Xe] 4f7 5d1 6s2",...           ## gd (64) 
        "[Xe] 4f9 6s2",...               ## tb (65) 
        "[Xe] 4f10 6s2",...              ## dy (66) 
        "[Xe] 4f11 6s2",...              ## ho (67) 
        "[Xe] 4f12 6s2",...              ## er (68) 
        "[Xe] 4f13 6s2",...              ## tm (69) 
        "[Xe] 4f14 6s2",...              ## yb (70) 
        "[Xe] 4f14 5d1 6s2",...          ## lu (71) 
        "[Xe] 4f14 5d2 6s2",...          ## hf (72) 
        "[Xe] 4f14 5d3 6s2",...          ## ta (73) 
        "[Xe] 4f14 5d4 6s2",...          ## w  (74) 
        "[Xe] 4f14 5d5 6s2",...          ## re (75) 
        "[Xe] 4f14 5d6 6s2",...          ## os (76) 
        "[Xe] 4f14 5d7 6s2",...          ## ir (77) 
        "[Xe] 4f14 5d9 6s1",...          ## pt (78) 
        "[Xe] 4f14 5d10 6s1",...         ## au (79) 
        "[Xe] 4f14 5d10 6s2",...         ## hg (80) 
        "[Xe] 4f14 5d10 6s2 6p1",...     ## tl (81) 
        "[Xe] 4f14 5d10 6s2 6p2",...     ## pb (82) 
        "[Xe] 4f14 5d10 6s2 6p3",...     ## bi (83) 
        "[Xe] 4f14 5d10 6s2 6p4",...     ## po (84) 
        "[Xe] 4f14 5d10 6s2 6p5",...     ## at (85) 
        "[Xe] 4f14 5d10 6s2 6p6",...     ## rn (86) 
        "[Rn] 7s1",...                   ## fr (87) 
        "[Rn] 7s2",...                   ## ra (88) 
        "[Rn] 6d1 7s2",...               ## ac (89) 
        "[Rn] 6d2 7s2",...               ## th (90) 
        "[Rn] 5f2 6d1 7s2",...           ## pa (91) 
        "[Rn] 5f3 6d1 7s2",...           ## u  (92) 
        "[Rn] 5f4 6d1 7s2",...           ## np (93) 
        "[Rn] 5f6 7s2",...               ## pu (94) 
        "[Rn] 5f7 7s2",...               ## am (95) 
        "[Rn] 5f7 6d1 7s2",...           ## cm (96) 
        "[Rn] 5f9 7s2",...               ## bk (97) 
        "[Rn] 5f10 7s2",...              ## cf (98) 
        "[Rn] 5f11 7s2",...              ## es (99) 
        "[Rn] 5f12 7s2",...              ## fm (100) 
        "[Rn] 5f13 7s2",...              ## md (101) 
        "[Rn] 5f14 7s2",...              ## no (102) 
        "[Rn] 5f14 7s2 7p1",...          ## lr (103) 
        "[Rn] 5f14 6d2 7s2",...          ## rf (104) 
        "[Rn] 5f14 6d3 7s2",...          ## db (105) 
        "[Rn] 5f14 6d4 7s2",...          ## sg (106) 
        "[Rn] 5f14 6d5 7s2",...          ## bh (107) 
        "[Rn] 5f14 6d6 7s2",...          ## hs (108) 
        "[Rn] 5f14 6d7 7s2",...          ## mt (109) 
        "[Rn] 5f14 6d8 7s2",...          ## ds (110) 
        "[Rn] 5f14 6d9 7s2",...          ## rg (111) 
        "[Rn] 5f14 6d10 7s2",...         ## cn (112) 
        "[Rn] 5f14 6d10 7s2 7p1",...     ## nh (113) 
        "[Rn] 5f14 6d10 7s2 7p2",...     ## fl (114) 
        "[Rn] 5f14 6d10 7s2 7p3",...     ## mc (115) 
        "[Rn] 5f14 6d10 7s2 7p4",...     ## lv (116) 
        "[Rn] 5f14 6d10 7s2 7p5",...     ## ts (117) 
        "[Rn] 5f14 6d10 7s2 7p6",...     ## og (118) 
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
  system(sprintf("~/git/espresso/bin/ld1.x < %s.in",atoms{i}));
  
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
  xocc = zeros(1,nwfn+2);
  for j = nwfn+2:-1:3
    occ = getfield(tocc,upper(aa{j}));
    fprintf(fid,"%d ",occ);
    xocc(j) = occ;
  endfor
  fprintf(fid,"\n");
  fprintf(fid,"%.10f %.10f %.10f %d\n",-6,atnums(i),0.002,n);
  rr = rho = zeros(1,n);
  for j = 1:n
    fprintf(fid,"%.14e ",wfn(j,1));
    for k = nwfn+1:-1:2
      fprintf(fid,"%.14e ",wfn(j,k));
    endfor
    rr(j) = wfn(j,1);
    rho(j) += sum(xocc(3:nwfn+2) .* wfn(j,2:nwfn+1).^2) / (4*pi*rr(j)^2);
    fprintf(fid,"\n");
  endfor
  fclose(fid);

  fid = fopen(sprintf("%s.rho",atoms{i}),"w");
  for j = 1:n
    fprintf(fid,"%.14e %.14e\n",rr(j),rho(j));
  endfor
  fclose(fid);
endfor
