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
    issubmodule[f[fs]] = 0
}
/^( |\t)*submodule( |\t)*/{
    ismodule[fs] = 1
    issubmodule[f[fs]] = 1
    gsub(/^.*\(/,"")
    gsub(/\).*$/,"")
    split($0,a,",")
    for (i in a)
        hassubmodule[a[i]] = 1
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
            str = sprintf(": %s.mod",f[i])
            for (j=1;j<=uses[tolower(f[i])];j++){
	    	nm = use[f[i],j]
		if (issubmodule[nm])
                   str = sprintf("%s.o %s.smod %s",nm,nm,str)
		else if (hassubmodule[nm])
                   str = sprintf("%s.o %s.mod %s.smod %s",nm,nm,nm,str)
                else
                   str = sprintf("%s.o %s.mod %s",nm,nm,str)
	    }
            print str 
        }
    }
}
' *.f90 *.F90
