#! /usr/bin/gawk -f

/\\begin{quote}\\begin{verbatim}/ {
   getline
   if ($0 ~ /runwien/) {
       print "\\runwienlist"
       print "\\begin{lstlisting}"
       inq = 1
       }
   else if ($0 ~ /ascii/) {
       print "\\asciilist"
       print "\\begin{lstlisting}"
       inq = 1
       }
   else if ($0 ~ /tessel/) {
       print "\\tessellist"
       print "\\begin{lstlisting}"
       inq = 1
       }
   else if ($0 ~ /critic/) {
       print "\\criticlist"
       print "\\begin{lstlisting}"
       inq = 1
       }
   else if ($0 ~ /octave/) {
       print "\\octavelist"
       print "\\begin{lstlisting}"
       inq = 1
       }
   else {
       print "\\criticlist"
       print "\\begin{lstlisting}"
       print
       inq = 1
       }
   next
   }

/\\end{verbatim}/ && inq {
   print "\\end{lstlisting}"
   inq = ""
   for(;$0 !~ /quote/;getline)
       ;
   next
   }

{ print }
