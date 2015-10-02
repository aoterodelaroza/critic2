#! /bin/sh

# Use:
#   makemake.sh
#
# Assumes files and modules share name and general common sense.

gawk '
{ 
    if (FILENAME != file[fs]){
        fs++
        file[fs] = FILENAME
        f[fs] = FILENAME
        gsub(/\.[fF]90$/,"",f[fs])
    }
}
/^( |\t)*module( |\t)*[^ \t\n]*( |\t)*$/{
    ismodule[fs] = 1
}
/^( |\t)*use( |\t)*[^ \t\n]*/{
    nm = tolower($2)
    idx = index(nm,",")
    if (idx != 0)
       nm = substr(nm,0,idx-1)
    for (i=1;i<=uses[nm];i++){
       if (use[nm,i] == f[fs])
          next
    }
    uses[nm]++
    use[nm,uses[nm]] = f[fs]
}

END{
    for (i=1;i<=fs;i++){
        if (ismodule[i] && uses[tolower(f[i])]){
            str = sprintf(": %s.$(MODEXT)",f[i])
            for (j=1;j<=uses[tolower(f[i])];j++)
                str = sprintf("%s.$(OBJEXT) %s",use[f[i],j],str)
            print str 
        }
    }
}
' *.f90 *.F90
