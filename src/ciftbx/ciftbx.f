C
C
C    \ | /            /##|    @@@@  @   @@@@@   |      |             @    @
C     \|/ STAR       /###|   @      @   @     __|__    |             @    @
C  ----*----        /####|  @       @   @@@@    |      |___  __  __  @@@@@@
C     /|\          /#####|   @      @   @       |      |   \   \/         @
C    / | \         |#####|    @@@@  @   @       \___/  \___/ __/\__       @
C                  |#####|________________________________________________
C                 ||#####|                 ___________________            |
C        __/|_____||#####|________________|&&&&&&&&&&&&&&&&&&&||          |
C<\\\\\\\\_ |_____________________________|&&& 29 Nov 2009  &&||          |
C          \|     ||#####|________________|&&&&&&&&&&&&&&&&&&&||__________|
C                  |#####|
C                  |#####|                Version 4.1.0 Release
C                  |#####|
C                 /#######\
C                |#########|
C                    ====
C                     ||
C           An extended tool box of fortran routines for manipulating CIF data.
C                     ||
C                     ||  CIFtbx Version 4
C                     ||        by
C                     ||
C                     ||  Sydney R. Hall (syd@crystal.uwa.edu.au)
C                     ||  Crystallography Centre
C                     ||  University of Western Australia
C                     ||  Nedlands 6009, AUSTRALIA
C                     ||
C                     ||       and
C                     ||
C                     ||  Herbert J. Bernstein (yaya@bernstein-plus-sons.com)
C                     ||  Bernstein + Sons
C                     ||  5 Brewster Lane
C                     ||  Bellport, NY 11713, U.S.A.
C                     ||
C The latest program source and information is available from:
C                     ||
C Em: syd@crystal.uwa.edu.au       ,-_|\      Sydney R. Hall
C sendcif@crystal.uwa.edu.au      /     \     Crystallography Centre
C Fx: +61 9 380 1118  ||      --> *_,-._/     University of Western Australia
C Ph: +61 9 380 2725  ||               v      Nedlands 6009, AUSTRALIA
C                     ||
C                     ||
C_____________________||_____________________________________________________
C
C This is a version of CIFtbx which has been extended to work with CIF2, DDLm,
C DDL 2 and mmCIF as well as with DDL 1.4 and core CIF dictionaries.  CIFtbx
C version 1 was written by Sydney R. Hall (see Hall, S. R., "CIF Applications
C IV.  CIFtbx: a Tool Box for Manipulating CIFs,"  J. Appl. Cryst (1993). 26,
C 482-494.  The revisions for version 2 were done by Herbert J. Bernstein
C and Sydney R. Hall (see Hall, S. R. and Bernstein, H. J., "CIFtbx 2:
C Extended Tool Box for Manipulating CIFs," J. Appl. Cryst.(1996). 29,
C 598-603)
C
C The revisions for releases 3 and 4 were done by Herbert J. Bernstein, work
C funded in part by the International Union of Crystallography
C
C___________________________________________________________________________
C
C
C    GENERAL TOOLS
C
C
C    init_      Sets the device numbers of files.   (optional)
C               [logical function always returned .true.]
C
C               <input CIF dev number> Set input CIF device     (def=1)
C
C               <output CIF dev number>Set output CIF device    (def=2)
C
C               <diracc dev number>    Set direct access formatted
C                                      scratch device number    (def=3)
C
C               <error  dev number>    Set error message device (def=6)
C
C
C
C    dict_      Requests a CIF dictionary be used for various data checks.
C               [logical function returned as .true. if the name dictionary
C               was opened and if the check codes are recognisable.  The
C               data item names used in the first dictionary loaded are
C               considered to be preferred by the user to aliases found
C               in dictionaries loaded in later calls.  On exit from dict_
C               the variable dicname_ is either equal to the filename, or,
C               if the dictionary had a value for the tag dictionary_name
C               of dictionary.title, dicname_ is set to that value.
C               The variable dicver_ is blank or set to the value of
C               _dictionary_version or of _dictionary.version  The check codes
C               'catck' and 'catno' turn on and off checking of dictionary
C               catgeory conventions.  The default is 'catck'.  The check
C               codes 'parck' and 'parno' turn on and off checking of
C               parent-child relationships.  The default of 'parck'.  Three check
C               codes control the handling of tags from the current dictionary
C               which duplicate tags from a dictionary loaded earlier.  These
C               codes ('first', 'final' and 'nodup') have effect only for the
C               current call to dict_  The default is 'first'.]
C
C               <dictionary filename>  A CIF dictionary in DDL format
C                                      or blank if just setting flags
C                                      or resetting the dictionary
C
C               <check code string>    The codes specifying the types of
C                                      checks to be applied to the CIF.
C
C                                      'valid'  data name validation check.
C                                      'dtype'  data item data type check.
C                                      'catck'  check datanames against
C                                               categories
C                                      'catno'  don't check datanames against
C                                               categories
C                                      'parck'  check datanames against
C                                               parent-child relationships
C                                      'parno'  don't check datanames against
C                                               parent-child relationships
C                                      'first'  accept first dictionary's
C                                               definitions of duplicate tags
C                                      'final'  accept final dictionary's
C                                               definitions of duplicate tags
C                                      'nodup'  do not accept duplicate tag
C                                               definitions
C                                      'parck'  check datanames against parent-
C                                               child relationahips
C                                      'parno'  don't check datanames against
C                                               parent-child relationships
C                                      'reset'  switch off checking flags
C                                      'close'  close existing dictionaries
C
C___________________________________________________________________________
C
C
C   CIF ACCESS TOOLS  ("the get_ing commands")
C
C
C
C    ocif_      Opens the CIF containing the required data.
C               [logical function returned .true. if CIF opened]
C
C               <CIF filename>        A blank name signals that the
C                                     currently open input CIF file
C                                     will be read.
C
C
C
C    data_      Identifies the data block containing the data to be requested.
C               [logical function returned .true. if block found]
C
C               <data block name>     A blank name signals that the next
C                                     encountered block is used (the block
C                                     name is stored in the variable bloc_).
C
C
C    bkmrk_     Saves or restores the current position so that data from
C               elsewhere in the cif can be examined.
C               [logical function returned as .true. on save if there was
C               room in internal storage to hold the current position, .true.
C               on restore if the bookmark number used was valid.  If the
C               argument is zero, the call is to save the position and return
C               the bookmark number in the argument.  If the argument is
C               non-zero, the call is to restore the position saved for the
C               bookmark number given.  The bookmark and the argument are
C               cleared.  The position set on return allow reprocessing of
C               the data item or loop row last processed when the bookmark
C               was placed.
C
C               NOTE:  All bookmarks are cleared by a call to data_]
C
C               <integer variable>    Bookmark number
C
C
C    find_      Find the location of the requested item in the CIF.
C               [The argument "name" may be a data item name, blank
C               for the next such item.  The argument "type" may be
C               blank for unrestricted acceptance of any non-comment
C               string (use cmnt_ to see comments), including loop headers,
C               "name" to accept only the name itself and "valu"
C               to accept only the value, or "head" to position to the
C               head of the CIF.  Except when the "head" is requested,
C               the position is left after the data item provided.  If the
C               item found is of type "name", posnam_ is set, otherwise,
C               posval_]
C
C               <data item name>      A blank name signals that the next
C                                     item of the type specified is needed
C
C               <data item type>      blank, 'head', 'name' or 'valu'
C
C               <character variable>  Returned string is of length long_.
C
C
C
C    test_      Identify the data attributes of the named data item.
C               [logical function returned as .true. if the item is present or
C               .false. if it is not. The data attributes are stored in the
C               common variables list_, type_, dictype_, diccat_ and dicname_.
C               The list, array, tuple or table attribites are stored in
C               ttype_, depth_ index_.
C
C               The values in dictype_, diccat_ and dicname_ are valid
C               whether or not the data item is found in the input CIF, as
C               long as the named data item is found in the dictionaries
C               declared by calls to dict_.  The data item name found
C               in the input CIF is stored in tagname_.  The appropriate
C               column numbers are stored in posnam_, posval_, posend_ and (for
C               numbers) in posdec_.  The quotation mark, if any, used is
C               stored in quote_.
C
C               list_ is an integer variable containing the sequential number
C               of the loop block in the data block. If the item is not within
C               a loop structure this value will be zero.
C
C               type_ is a character*4 variable with the possible values:
C                      'numb'  for number data
C                      'char'  for character data
C                      'text'  for text data
C                      'null'  if data missing or '?' or '.'
C                              also used for blank quoted fields if
C                              nblank_ is true
C
C               ttype_ is a character*4 variable with the container type:
C                      'list'  for list or array data   [item,...]
C                      'tupl'  for tuple data           (item,...)
C                      'tabl'  for table data           {item,...}
C               The meanings change if rdbkt_, rdbrc_ or rdprn_ are
C               false.  If rdbkt_ is false, the meanings are
C                      'tupl'  for tuple data           (item,...)
C                      'list'  for list or table data   {item,...}
C               If rdprn_ is false, then 'list' is used for all
C               container types.  If depth_ is 0, then ttype_ is not
C               valid and will contain '    '
C
C               depth_ is an integer variable with the depth into a
C               list, array, tuple or table.  A depth of zero means that
C               no list, array, tuple or table is being processed.
C
C               index_ is an integer variable with the index (from 1)
C               across a list, array, tuple or table.  An index of zero
C               means that no list, array, tuple or table is being processed.
C
C               dictype_ is a character*(NUMCHAR) variable with the type code
C               given in the dictionary entry for the named data item.  If
C               no dictionary was used, or no type code was specified, this
C               field will simply agree with type_.  If a dictionary was used,
C               this type may be more specific than the one given by type_.
C
C               diccat_ is a character*(NUMCHAR) variable with the category
C               of the named data item, or '(none)'
C
C               dicname_ is a character*(NUMCHAR) variable with the name of
C               the data item which is found in the dictionary for the
C               named data item.  If alias_ is .true., this name may
C               differ from the name given in the call to test_.  If alias_
C               is .false. or no preferred alias is found, dicname_ agrees with
C               the data item name.
C
C               tagname_ is a character*(NUMCHAR) variable with the name
C               of the data item as found in the input CIF.  It will be
C               blank if the data item name requested is not found in the
C               input CIF and may differ from the data item name provided
C               by the user if the name used in the input CIF is an
C               alias of the data item name and alias_ is .true.
C
C               posnam_, posval_, posend_  and posdec_ are integer variables
C               which may be examined if information about the horizontal
C               position of the name and data read are needed.  posnam_ is the
C               starting column of the data name found (most often 1).
C               posval_ is the starting column of the data value.  If the
C               field is numeric, then posdec_ will contain the effective
C               column number of the decimal point.  For whole numbers, the
C               effective position of the decimal point is one column to the
C               right of the field.  posend_ contains the ending column of the
C               data value.
C
C               quote_ is a character*3 variable which may be examined to
C               determine if a quotation character was used on character data.]
C
C               <data name>           Name of the data item to be tested.
C
C
C    dtype_     Return the dictionary type of a data name, if any.
C               [logical function returned as .true. if the item has a type
C               in the dctionary, .false. if not.  The type returned is
C               one of the base type used by type_ (see above), if possible]
C
C               <data name>          Name of the item for which a type is needed
C               <data type>          Returned type from the dictionary
C
C
C    name_      Get the NEXT data name in the current data block.
C               [logical function returned as .true. if a new data name exists
C               in the current data block, and .false. when the end of the data
C               block is reached.]
C
C               <data name>           Returned name of next data item in block.
C
C
C
C    numb_      Extracts the number and its standard deviation (if appended).
C               [logical function returned as .true. if number present. If
C               .false. arguments 2 and 3 are unaltered. If the esd is not
C               attached to the number argument 3 is unaltered.]
C
C               <data name>           Name of the number sought.
C
C               <real variable>       Returned number.
C
C               <real variable>       Returned standard deviation.
C
C
C
C    numd_      Extracts the number and its standard deviation (if appended)
C               as double precision variables.
C               [logical function returned as .true. if number present. If
C               .false. arguments 2 and 3 are unaltered. If the esd is not
C               attached to the number argument 3 is unaltered.]
C
C               <data name>           Name of the number sought.
C
C               <double precision variable>
C                                     Returned number.
C
C               <double precision variable>
C                                     Returned standard deviation.
C
C
C
C    char_      Extracts character and text strings.
C               [logical function returned as .true. if the string is present.
C               Note that if the character string is text this function is
C               called repeatedly until the logical variable text_ is .false.
C               Non-text blank (quoted blanks) or empty ('' or "") fields
C               are converted by char to a null field, if nblank_ is true.]
C
C               <data name>           Name of the string sought.
C
C               <character variable>  Returned string is of length long_.
C
C    charnp_    Extracts character and text strings.
C               [logical function returned as .true. if the string is present.
C               Note that if the character string is text this function is
C               called repeatedly until the logical variable text_ is .false.
C               If the value is found in a container, then charnp_ should
C               be called repeatedly until both text_ is false and depth_
C               is zero. 
C
C               Non-text blank (quoted blanks) or empty ('' or "") fields
C               are converted by char to a null field, if nblank_ is true.
C               Only the number of characters returned in the third argument
C               are set.  This value is never less than 1, but may be less
C               than the allocated length of the returned string.]
C
C               <data name>           Name of the string sought.
C
C               <character variable>  Returned string is of length long_.
C
C               <integer variable>    Returned length of valid characters.
C
C
C    cmnt_      Extracts the next comment from the data block.
C               [logical function returned as .true. if a comment is present.
C               The initial comment character "#" is _not_ included in the
C               returned string.  A completely blank line is treated as
C               a comment.  A comment may be extracted while reading a list,
C               array, tuple or table]
C
C               <character variable>  Returned string is of length long_.
C
C
C    delim_     Reports the most recently seen delimiter prior to the
C               most recently extracted tag or value at the specified
C               depth.  Outside of bracketed constructs, only delimiters
C               at depth 0 (top level) can be seen.  This is not the
C               quoting character for a quoted string or text field.
C               See the variable quote_.
C               [logical function returned as .true. if the depth is
C               not negative and greater than or equal to the current
C               depth.  At depth 0, in a correctly formatted CIF, the
C               delimiter returned is always a blank,]
C
C               <integer variable>    Depth
C                    
C               <character variable>  Returned string is of length 1
C
C               <integer variable>    column position of delimiter
C
C               <integer variable>    record position of delimiter
C
C
C
C    purge_     Closes existing data files and clears tables and pointers.
C               [subroutine call]
C
C____________________________________________________________________________
C
C
C
C   CIF CREATION TOOLS ("the put_ing commands")
C
C
C
C    pfile_     Create a file with the specified file name.
C               [logical function returned as .true. if the file is opened.
C               The value will be .false. if the file already exists.]
C
C               <file name>           Blank for use of currently open file
C
C
C
C    pdata_     Put a data block command into the created CIF.
C               [logical function returned as .true. if the block is created.
C               The value will be .false. if the block name already exists.
C               Produces a save frame instead of a data block if the
C               variable saveo_ is true during the call.  No block duplicate
C               check is made for a save frame.]
C
C               <block name>
C
C    pdelim_    Emit a specific delimiter
C               [logical function returned as .true. if the delimiter is
C               appropriate to the context.  Emitting a '(', '{' or '['
C               increases the output depth by one.  Emitting a ')', '}'
C               or ']' decreases the output depth by one.  Emitting a ' ',
C               ',' or ':' does not change the depth.  Emitting a ','
C               or ':' at depth_ 0 is an error that can be overridden
C               by the second argument being .true..  Emitting a ' ' at
C               a depth_ greater than 0 is an error that can be overridden
C               by the second argument being .true.. ]
C
C               <character variable>   The one-character delimiter string
C
C               <logical variable>     .true. if an invalid delimiter is
C                                      to be forced out
C
C               <integer variable>     Column position at which to write
C                                      the delimiter or 0 if not specified
C
C
C
C    ploop_     Put a loop_ data name into the created CIF.
C               [logical function returned as .true. if the invocation
C               conforms with the CIF logical structure.  If pposval_
C               is non-zero, the "loop_" header is positioned to
C               that column.  If pposnam_ is non-zero, the data name is
C               positioned to that column.]
C
C               <data name>         If the name is blank on the first call
C                                   of a loop, only the "loop_" is placed.
C
C
C
C    pchar_     Put a character string into the created CIF.
C               [logical function returned as .true. if the name is unique,
C               AND, if dict_ is invoked, is a name defined in the dictionary,
C               AND, if the invocation conforms to the CIF logical structure.
C               The action of pchar_ is modified by the variables pquote_ and
C               nblanko_.  If pquote_ is non-blank, it is used as a quotation
C               character for the string written by pchar_.  The valid values
C               are '''', '"', ';', '(', '{', '[', '''''''', and '"""'.
C               In the last six cases a text field, bracketed construct or
C               multi-line triple-quoted string is written.  If the string
C               contains a matching character to the value of quote_, or if
C               quote_ is not one of the valid quotation characters, a valid,
C               non-conflicting quotation character is used or the line-folding
C               conventions are used to prevent the close-quote from being
C               followed by white space.  Except when writing a text field, if
C               nblanko_ is true, pchar_ converts a blank string to
C               an unquoted period.]
C
C               <data name>         If the name is blank, do not output name.
C
C               <character string>  A character string of MAXBUF chars or less.
C
C
C
C    pcmnt_     Puts a comment into the created CIF.
C               [logical function returned as .true.  The comment character
C               "#" should not be included in the string.  A blank comment
C               is presented as a blank line without the leading "#"].
C
C               <character string>  A character string of MAXBUF chars or less.
C
C
C    pnumb_     Put a single precision number and its esd into the created CIF.
C               [logical function returned as .true. if the name is unique,
C               AND, if dict_ is invoked, is a name defined in the dictionary,
C               AND, if the invocation conforms to the CIF logical structure.
C               The number of esd digits is controlled by the variable
C               esdlim_]
C
C               <data name>         If the name is blank, do not output name.
C
C               <real variable>     Number to be inserted.
C
C               <real variable>     Esd number to be appended in parentheses.
C
C
C    pnumd_     Put a double precision number and its esd into the created CIF.
C               [logical function returned as .true. if the name is unique,
C               AND, if dict_ is invoked, is a name defined in the dictionary,
C               AND, if the invocation conforms to the CIF logical structure.
C               The number of esd digits is controlled by the variable
C               esdlim_]
C
C               <data name>         If the name is blank, do not output name.
C
C               <double precision variable>
C                                   Number to be inserted.
C
C               <double precision variable>
C                                   Esd number to be appended in parentheses.
C
C
C
C    ptext_     Put a character string into the created CIF.
C               [logical function returned as .true. if the name is unique,
C               AND, if dict_ is invoked, is a name defined in the dictionary,
C               AND, if the invocation conforms to the CIF logical structure.
C               ptext_ is invoked repeatedly until the text is finished. Only
C               the first invocation will insert a data name.
C
C               If used when pclipt_ is .true. if the first character of the
C               text field is blank, it is removed.
C
C               If used when pfold_ is non-zero, the text field will be marked
C               as folded even if the first line is small enough to fit.
C               In order to produce a non-folded text field in the midst
C               of generally folded items, pfold_ should be set to 0 before
C               calling ptext_ and then restored to the previous value.]
C
C               <data name>         If the name is blank, do not output name.
C
C               <character string>  A character string of MAXBUF chars or less.
C
C
C    prefx_     Puts a prefix onto subsequent lines of the created CIF.
C               [logical function returned as .true.  The second argument
C               may be zero to suppress a previously used prefix, or
C               greater than the non-blank length of the string to force
C               a left margin.  Any change in the length of the prefix
C               string flushes pending partial output lines, but does _not_
C               force completion of pending text blocks or loops.
C               This function allows the CIF output functions to be used
C               within what appear to be text fields to support annotation
C               of a CIF. ]
C
C               <character string>  A character string of MAXBUF chars or less.
C
C               <integer variable>  The length of the prefix string to use.
C
C
C
C
C    close_     Close the creation CIF. MUST be used if pfile_ is used.
C               [subroutine call]
C
C
C____________________________________________________________________________
C
C
C
C....The CIF tool box also provides variables for data access control:
C
C
C    alias_      Logical variable: if left .true. then all calls to
C                CIFtbx functions may use aliases of data item names.
C                The preferred synonym from the dictionary will be
C                subsituted internally, provided aliased data names were
C                supplied by an input dictionary (via dict_).  The
C                default is .true., but alias_ may be set to .false.
C                in an application.
C
C    aliaso_     Logical variable: if set .true. then cif output
C                routines will convert aliases to the names to preferred
C                synonyms from the dictionary.  The default is .false., but
C                aliaso_ may be set to .true. in an application.  The
C                setting of aliaso_ is independent of the setting of
C                alias_.
C
C    align_      Logical variable signals alignment of loop_ lists during
C                the creation of a CIF. The default is .true.
C
C    append_     Logical variable:  if set .true. each call to ocif_ will
C                append the information found to the current cif.  The default
C                is .false.
C
C    bloc_       Character*(NUMCHAR) variable: the current block name.
C
C    clipt_      Logical variable: if set .true., when reading text fields,
C                an extra blank is inserted before the character string
C                returned for the first line of a text field, emulating
C                the behavior of CIFtbx versions prior to version 4.
C
C    decp_       Logical variable: set when processing numeric input, .true.
C                if there is a decimal point in the numeric value, .false.
C                otherwise
C
C    depth_      Integer variable: set to the depth within a list, array, tuple
C                or table
C
C    dictype_    Character*(NUMCHAR) variable: the precise data type code
C                (see test_)
C
C    diccat_     Character*(NUMCHAR) variable: the category (see test_)
C
C    dicname_    Character*(NUMCHAR) variable: the root alias (see test_) or
C                the name of the dictionary just loaded (see dict_)
C
C    dicpname_   Character*(NUMCHAR) variable: the parent (see test_)
C
C    dicver_     Character*(NUMCHAR) variable: the version of the dictionary
C                just loaded (see dict_)
C
C    esdlim_     Integer variable:  Specifies the upper limit of esd's
C                produced by pnumb_, and, implicitly, the lower limit.
C                The default value is 19, which limits esd's to the range
C                2-19.  Typical values of esdlim_ might be 9 (limiting
C                esd's to the range 1-9), 19, or 29 (limiting esd's
C                to the range 3-29).  If esdlim_ is given as a negative
C                value, the upper limit of esd's is the absolute value
C                of esdlim_ and the lower limit is 1.
C
C    esddig_     Integer variable:  The number of esd digits in the last
C                number read from a CIF.  Will be zero if no esd
C                was given.
C
C    file_       Character*(MAXBUF) variable: the filename of the current file.
C                Warning:  only file_(1:longf_) is valid
C
C    fold_       Logical variable signals that the current text block
C                began with the ';\' fold indicator. Only meaningful
C                when text_ is .true. and type_ is 'text'.
C                (fold_ is .true. if the indicator is present)
C
C    glob_       Logical variable signals that the current data block
C                is actually a global block (.true. for a global block).
C
C    globo_      Logical variable signals that the output data block from
C                pdata_ is actually a global block (.true. for a global block).
C
C    index_      Integer variable: Specifies the one-based index of the current
C                item in a list, array, tuple or table
C
C    line_       Integer variable: Specifies the input/output line limit
C                for processing a CIF. The default value is 80 characters.
C                This may be set by the program. The max value is MAXBUF
C                which has a default value of 2048.  In order to use
C                the CIF 1.1 line folding protocol for lines that
C                cannot be fit into line_ characters, the variable
C                pfold_ must be set to a non-zero value less than
C                or equal to line_
C
C    list_       Integer variable: the loop block number (see test_).
C
C    long_       Integer variable: the length of the data string in strg_.
C
C    longf_      Integer variable: the length of the filename in file_.
C
C    loop_       Logical variable signals if another loop packet is present.
C
C    lzero_      Logical variable: set when processing numeric input, .true.
C                if the numeric value is of the form [sign]0.nnnn rather than
C                [sign].nnnn, .false. otherwise
C
C    nblank_     Logical variable: if set .true. then all calls to
C                to char_ or test_ which encounter a non-text quoted blank
C                will return the type as 'null' rather than 'char'.
C
C    nblanko_    Logical variable: if set .true. then cif output
C                routines will convert quoted blank strings to an
C                unquoted period (i.e. to a data item of type null).
C
C    pclipt_     Logical variable: if set .true., when writing text fields,
C                if there is a blank as the first character to be written,
C                it is removed, emulating the behavior of CIFtbx versions
C                prior to version 4.
C
C    pdecp_      Logical variable: if set .true. then cif numeric output
C                routines will insert a decimal point in all numbers written by
C                pnumb_ or pnumbd_.  If set .false. then a decimal point will be
C                written only when needed.  The default is .false.
C
C    pesddig_    Integer variable: if set non-zero, and esdlim_ is negative,
C                controls the number of digits for esd's produced by
C                pnumb_ and pnumd_
C
C    pfold_      Integer variable:  If set non-zero, specifies a column
C                on which output lines are to be folded.  The default is 0.
C                If pfold_ is set to a value greater than line_ the
C                value of line_ will be used instead.  Non-zero values of
C                pfold_ less than 4 are not valid and will be reset to 4.
C                Non-zero values of pfold_ less than 80 can cause conflict
C                with the syntactic requirements of creating a valid CIF.
C
C    plzero_     Logical variable: if set .true. then cif numeric output
C                routines will insert a zero before a leading decimal point,
C                The default is .false.
C
C    pposdec_    Integer variable giving the position of the decimal point
C                for the next number to be written.  This acts very much like
C                a decimal centered tab in a word processor, to help align
C                columns of number on a decimal point, if a decimal point
C                is present.
C
C    pposend_    Integer variable giving the ending column of the next
C                number or quoted character value to be written.  Used to
C                pad with zeros or blanks.
C
C    pposnam_    Integer variable giving the starting column of the next
C                name or comment or data block to be written.
C
C    pposval_    Integer variable giving the starting column of the next
C                data value to be written by pchar_, pnumb_ or pnumd_.
C                Also used to set the position of the initial "loop_"
C                in a ploop_ call or to set the position of a terminal "save_"
C                for a save frame in a pdata_ call for which saveo_ is .true.
C
C    posdec_     Integer variable giving the position of the decimal point
C                for the last number read, if a decimal point was present.
C
C    posend_     Integer variable giving the ending column of the last
C                data value read, not including a terminal quote.
C
C    posnam_     Integer variable giving the starting column of the last
C                name or comment or data block read.
C
C    posval_     Integer variable giving the starting column of the last
C                data value read.  Also reports the column of the
C                terminal "save_" of a save frame.
C
C    pquote_     Character variable giving the quotation symbol to be
C                used for the next string written, or the comment
C                flag for the next comment written.
C
C    precn_      Integer variable:  Reports the record number of the last
C                line written to the output cif.  Set to zero by init_.  Also
C                set to zero by pfile_ and close_ if the output cif file name
C                was not blank.
C
C    ptabx_      Logical variable signals tab character expansion to blanks
C                during the creation of a CIF. The default is .true.
C
C    quote_      Character variable giving the quotation symbol found
C                delimiting the last string read or the comment flag
C                for the last comment read.  The possible valid values
C                are '''', '"', ';', '''''''', and '"""'.
C                The treble quotes are recognized only if rdtq_ is .true.
C
C    rdbrc_      Logical variable:  control recognition of { ... } constructs
C                on read.  The default is .false.
C
C    rdbkt_      Logical variable:  controls recognition of [ ... ] constructs
C                on read.  The default is .false.
C
C    rdprn_      Logical variable:  controls recognition of ( ... ) constructs
C                on read.  The default is .false.
C
C    rdtq_       Logical variable:  controls recognition of """ ... """ and
C                ''' ... ''' constructs on read.  The default is .false.
C
C    rdrcqt_     Logical variable:  controls recognition of trailing punctuation
C                after a closing quote.  If .true. a closing quotation mark is
C                recognized immediately, no matter what follows the closing
C                quoation mark (the CIF 2 convention).  If .false., a closing
C                quotation mark is only effective if followed by a blank, or,
C                in bracketed constructs by a blank, a colon, a comma or 
C                the closing bracket.
C
C    recbeg_     Integer variable:  Gives the record number of the first
C                record to be used.  May be changed by the user to restrict
C                access to a CIF.
C
C    recend_     Integer variable:  Gives the record number of the last
C                record to be used.  May be changed by the user to restrict
C                access to a CIF.
C
C    recn_       Integer variable:  Reports the record number of the last
C                line read from the direct access copy of the input cif.
C
C    save_       Logical variable signals that the current data block
C                is actually a save-frame (.true. for a save-frame).
C
C    saveo_      Logical variable signals that the output data block from
C                pdata_ is actually a save-frame (.true. for a save-frame).
C
C    strg_       Character*(MAXBUF) variable: the current data item.
C
C    tabl_       Logical variable signals tab-stop alignment of output
C                during the creation of a CIF. The default is .true.
C
C    tabx_       Logical variable signals tab character expansion to blanks
C                during the reading of a CIF. The default is .true.
C
C    tbxver_     Character*32 variable: the CIFtbx version and date
C                in the form 'CIFtbx version N.N.N, DD MMM YY '
C
C    text_       Logical variable signals if another text line or is present.
C
C    type_       Character*4 variable: the data type code (see test_).
C
C    ttype_      Character*4 variable: the list, array, tuple or table type code (see test_).
C
C    unfold_     Logical variable signals that input lines are to be
C                unfolded before presentation of data.  The default
C                is .false.
C
C    xmlout_     Logical variable:  Set by the user to change the output
C                style to XML conventions.  Note that this is not a
C                cml output, but a literal translation from the input CIF.
C                The default is .false.
C
C    xmlong_     Logical variable:  Set by the user to change the style of
C                xml output if xmlout_ is .true.  When .true. (the default)
C                xml tag names are the full CIF tag names with the leading
C                '_' removed.  When .false. an attempt is made to strip
C                the leading category name as well.
C
C
C_____________________________________________________________________________
C
C
C >>>>>> Set the device numbers.
C
         function init_(devcif,devout,devdir,deverr)
C
         logical   init_
         include   'ciftbx.sys'
         integer   devcif,devout,devdir,deverr
         integer   ii,kdig
         real      ytest
         double precision ztest
         double precision tbxxdble
         real      tbxxsngl
C
         init_=.true.
         cifdev=devcif
         outdev=devout
         dirdev=devdir
         errdev=deverr

         recn_=0
         precn_=0
         plcat = ' '
         plxcat = ' '
         plhead(1) = ' '
         plxhead(1) = ' '
         pdblok = ' '
         ploopn = 0
         nstable = 0
         nivt = 0
C
C        recompute decimal single precision precision
C        This is found by computing the smallest power of
C        10 which, when added to 1, produces a change
C        and then backing off by 1
C
         decprc = .1
         do ii = 1,8
         ytest = tbxxsngl(1.+decprc/10.)
         if (ytest.eq.1.) go to 100
         decprc = decprc/10.
         enddo
100      continue
         decprc=decprc*10.
C
C        recompute decimal double precision precision
C
         kdig = 1
         dpprc = .1D0
         do ii = 1,17
         ztest = tbxxdble(1.D0+dpprc/10.)
         if (ztest.eq.1.D0) go to 200
         dpprc = dpprc/10.D0
         kdig = kdig+1
         enddo
200      continue
         dpprc=dpprc*10.D0
         write(ndpfmt,'(5h(d30.,i2,1h))') kdig-1
C
C        recompute decimal single precision minimum power of ten
C
         decmin = .1
         do ii = 1,39
         ytest = decmin/10.
         if (ytest.eq.0.) go to 300
         decmin = decmin/10.
         enddo
300      continue
C
C        recompute decimal double precision minimum power of 10
C        and its log base 10 (minexp)
C
         dpmin = .1D0
         minexp = -1
         do ii = 1,309
         ztest = dpmin/10.
         if (ztest.eq.0.D0) go to 400
         dpmin = dpmin/10.D0
         minexp = minexp-1
         enddo
400      continue
         call clearfp
         return
         end
C
C
C >>>>>> Function to defeat the optimizer
C
C     
         function tbxxdble(x)
         double precision x
         double precision tbxxdble
         tbxxdble = x
         return
         end
C
C
C >>>>>> Function to defeat the optimizer
C
C     
         function tbxxsngl(x)
         real x
         real tbxxsngl
         tbxxsngl = x
         return
         end
C
C
C
C
C
C >>>>>> Read a CIF dictionary and prepare for checks
C
         function dict_(fname,checks)
C
         logical   dict_
         logical   ocif_
         logical   data_
         logical   charnp_
         logical   test_
         integer   lastnb
         include  'ciftbx.sys'
         logical   tbxxnewd, tbxxoldd
         logical   nresult
         character fname*(*),checks*(*)
         character temp*(MAXBUF)
         character codes(11)*5,name*(MAXBUF),bxname*(NUMCHAR)
         character bpname*(NUMCHAR)
         character bcname*(NUMCHAR),bname*(NUMCHAR)
         character baname*(NUMCHAR),ganame*(NUMCHAR),btname*(NUMCHAR)
         character batag*(NUMCHAR)
         character mcstrg*(NUMCHAR)
         character riname*(NUMCHAR),rfname*(NUMCHAR)
         character xdicnam*(NUMCHAR)
         character xdicver*(NUMCHAR)
         character xmtoken*(NUMCHAR),xmtarg*(XMLCHAR),xmtyp*(NUMCHAR)
         character xxxtemp*(NUMCHAR)
         character*3 ovchk, otchk
         integer   nrecds,recends,recbegs
         integer   lchecks,lbpname,lbcname,lbaname,lbname
         integer   kdict,ifind,jfind,iafind,ick
         integer   i,j,nmatch,mycat,ksmatch,ii,jj,idstrt,kdup
         integer   nmycat,ixmtyp,nxmc,kxmc
         integer   lstrg,lxmtoken,lxmtarg,lxmtyp,kvrtp,kstrg,sindex
         integer   lbloc,kivt

C
C        Control flags for matching categories, names and types
C
C        icloop is the loop number of the block for the
C        current category
C        ictype is the type of the current category
C          0 - none found yet
C          1 - _item.category_id             (DDL2)
C          2 - _category                     (DDL1)
C          3 - _category.id                  (DDL2)
C          4 - _name.category_id             (DDLm)
C        the last ictype entry is not a type, but a tag
C        whose value may specify that this is a category
C        with the category name given under intype
C          5 - _definition.scope             (DDLm)
C        inloop is the loop number of the block for the
C        current name
C        intype is the type of the current name
C          0 - none found yet
C          1 - _item.name                    (DDL2)
C          2 - _name                         (DDL1)
C          3 - _definition.id                (DDLm)
C        ialoop is the loop number of the block for the
C        current alias
C        iatype is the type for the current alias
C          0 - none found yet
C          1 - _item_aliases.alias_name      (DDL2)
C          2 - _aliases.definition_id        (DDL2)
C        imloop is the loop number of the block for the
C        current parent
C        imtype is the type for a mandatory item
C          0 - none found yet
C          1 - _item.mandatory_code          (DDL2)
C          2 - _category_mandatory.item_id   (DDLm)
C        iptype is the type for the current parent
C          0 - none found yet
C          1 - _item_linked.parent_name      (DDL2)
C          2 - _item_link_parent             (DDL1)
C          3 - _category.parent_id           (DDLm)
C          4 - _name.linked_item_id          (DDLm)
C        itloop is the loop number of the block for the
C        current type
C        ittype is the type of the current type
C          0 - none found yet
C          1 - _item_type.code               (DDL2)
C          2 - _type                         (DDL1)
C          3 - _type.contents                (DDLm)
C        iritype is the type of the current related item
C          0 - none found yet
C          1 - _item_related.related_name    (DDL2)
C          2 - _related_item                 (DDL1)
C          3 - _type.purpose                 (DDLm)
C        irftype is the type of the current related item function
C          0 - none found yet
C          1 - _item_related.function_code   (DDL2)
C          2 - _related_function             (DDL1)
C          3 - _type.purpose                 (DDLm)
C

         integer icloop,ictype,inloop,intype,ialoop,iatype,
     * imloop,imtype,iptype,itloop,ittype,
     * iritype,irftype,icktype
C
         character*4 map_type(19),map_to(19),mapped
         character*(NUMCHAR) dt(2),dv(2),ct(5),nt(3),at(2),tt(3)
         character*(NUMCHAR) ri(3),rf(3),ck(4),pt(4),pc(2),mc(3)
         character*(NUMCHAR) ve(3),vr(4)
         data map_type
     *   /'floa','int ','yyyy','symo','ucha','ucod','name','idna',
     *    'any ','code','line','ulin','atco','fax ','phon','emai',
     *    'real','inte','coun'/
         data map_to
     *   /'numb','numb','char','char','char','char','char','char',
     *    'char','char','char','char','char','char','char','char',
     *    'numb','numb','numb'/
         data ri
     *      /'_item_related.related_name      ',
     *       '_related_item                   ',
     *       '_type.purpose                   '/
         data rf
     *      /'_item_related.function_code     ',
     *       '_related_function               ',
     *       '_type.purpose                   '/
         data dt
     *      /'_dictionary.title               ',
     *       '_dictionary_name                '/
         data dv
     *      /'_dictionary.version             ',
     *       '_dictionary_version             '/
         data ct
     *      /'_item.category_id               ',
     *       '_category                       ',
     *       '_category.id                    ',
     *       '_name.category_id               ',
     *       '_definition.scope               '/
         data nt
     *      /'_item.name                      ',
     *       '_name                           ',
     *       '_definition.id                  '/
         data at
     *      /'_item_aliases.alias_name        ',
     *       '_aliases.definition_id          '/
         data tt
     *      /'_item_type.code                 ',
     *       '_type                           ',
     *       '_type.contents                  '/
         data ck
     *      /'_category_key.name              ',
     *       '_list_reference                 ',
     *       '_category_key.generic           ',
     *       '_category_key.primitive         '/
         data pt
     *      /'_item_linked.parent_name        ',
     *       '_item_link_parent               ',
     *       '_category.parent_id             ',
     *       '_name.linked_item_id            '/
         data pc
     *      /'_item_linked.child_name         ',
     *       '_item_link_child                '/
         data mc
     *      /'_item.mandatory_code            ',
     *       '_mandatory                      ',
     *       '_category_mandatory.item_id     '/
         data ve
     *      /'_item_enumeration.value         ',
     *       '_enumeration                    ',
     *       '_enumeration_set.state          '/
         data vr
     *      /'_item_range.minimum             ',
     *       '_enumeration_range              ',
     *       '_item_range.maximum             ',
     *       '_enumeration.range              '/

C
         data codes /'valid','dtype','reset','close',
     *       'catck','catno','nodup','final','first',
     *       'parck','parno'/
C
         nrecds=nrecd
         recbegs=recbeg_
         recends=recend_
         if(append_) then
           recbeg_=nrecd
         endif
C
C        Initialize kdup to 0 ('final')
C
         kdup = 0
C
C        initialize both xdicnam and xdicver to blank
C
         xdicnam = ' '
         xdicver = ' '
C
C        preserve entry values of tcheck and vcheck in case dict fails
C
         otchk = tcheck
         ovchk = vcheck
C
C....... Are the codes OK
C
         lchecks=min(len(temp),len(checks))
         call tbxxnlc(temp(1:lchecks),checks)
         i=0
120      i=i+1
         if(i.ge.lchecks)            goto 190
         if(temp(i:i).eq.' ')        goto 120
         do 150 j=1,11
         if(temp(i:i+4).eq.codes(j)) goto 170
150      continue
         dict_=.false.
         goto 500
170      i=i+4
         if(j.eq.1) then
           vcheck='yes'
           goto 120
         endif
         if(j.eq.2) then
           tcheck='yes'
           goto 120
         endif
         if(j.eq.3) then
           vcheck = 'no '
           tcheck = 'no '
           goto 120
         endif
         if(j.eq.4) then
           vcheck = 'no '
           tcheck = 'no '
           catchk = 'yes'
           ndcname = 0
           ndict = 0
           if(nname.gt.0) then
           do 180 i = 1,nname
             dtype(i)=' '
             dxtyp(i)=' '
             cindex(i)=0
             ddict(i)=0
180        continue
           endif
           dict_=.true.
           goto 500
         endif
         if (j.eq.5) then
           catchk = 'yes'
           goto 120
         endif
         if (j.eq.6) then
           catchk = 'no '
           goto 120
         endif
         if (j.eq.10) then
           parchk = 'yes'
           goto 120
         endif
         if (j.eq.11) then
           parchk = 'no '
           goto 120
         endif
         kdup=j-8
         goto 120
C
C        if no category names have been loaded, clean up
C        the hash table for dictionary category names
C
190      if(ndcname.eq.0) then
           call hash_init(dcname,dcchain,NUMDICT,ndcname,dchash,
     *     NUMHASH)
         endif
C
C        if no dictionary names have been loaded, clean up
C        the hash table for dictionary names
C
         if(ndict.eq.0) then
           call hash_init(dicnam,dicchain,NUMDICT,ndict,dichash,
     *     NUMHASH)
         endif
         idstrt=ndict
C
C....... Open and store the dictionary
C
         dict_=.true.
         if(fname.eq.' ')            goto 500
         if(nname.gt.0) call tbxxerr(' Dict_ must precede ocif_')
         dict_=ocif_(fname)
         if(.not.dict_)              goto 500
         dictfl='yes'
C
C        At this point is is proper to update xdicnam to fname
C
         xdicnam = fname
C
C....... Loop over data blocks; extract _name's, _type etc.
C
200      if(.not.data_(' '))         goto 400
         lbloc = lastnb(bloc_)
         if(bloc_(1:1).eq.'_'.or.glob_.or.bloc_.eq.' ') then
           call tbxxclc(bname,lbname,bloc_(1:lbloc),lbloc)
         else
           call tbxxclc(bname,lbname,'_'//bloc_(1:lbloc),lbloc+1)
         endif
C
C        see if this is a dictionary defining block
C
         do i = 1,2
           if(charnp_(dt(i),name,lstrg)) then
             xdicnam = name(1:lstrg)
             do j = 1,2
               if(test_(dv(j))) then
                 xdicver = strg_(1:max(1,long_))
                 goto 200
               endif
             enddo
             goto 200
           endif
         enddo
C
Cdbg     WRITE(6,*) ndict,bloc_
C
C        Analyze loop structure for categories, names, types and parents
C
C
C        initalize loop info
C
         icloop = -1
         inloop = -1
         ialoop = -1
         imloop = -1
         itloop = -1
         ictype = 0
         intype = 0
         iatype = 0
         imtype = 0
         iptype = 0
         ittype = 0
         iritype = 0
         irftype = 0
         icktype = 0
         ixmtyp = 0
         bcname = ' '
         bpname = ' '
         lbcname = 1
         lbpname = 1
         baname = ' '
         batag = ' '
         lbaname = 1
         btname = ' '
         mycat=0
         loop_=.false.
         loopnl=0
         nmatch=0
         ksmatch=0
         riname = ' '
         rfname = ' '
C
C        Pick up category_keys and list_references
C
         do i = 1,4
210        if(charnp_(ck(i),name,lstrg)) then
             if (icktype.ne.0 .and. icktype.ne.i)
     *         call tbxxwarn
     *         (' Multiple DDL 1, 2 or m related key definitions ')
             icktype = i
             if (tbxxnewd(name(1:lstrg),ick)) then
               catkey(ick) = .true.
             else
               if(.not.catkey(ick)) then
                 ifind = aroot(ick)
215              catkey(ifind) = .true.
                 ifind = alias(ifind)
                 if (ifind.ne.0) go to 215
               endif
             endif
             if (loop_) go to 210
           endif
         enddo

C
C        Process related items
C
         do i = 1,2
           if(charnp_(ri(i),name,lstrg)) then
             if (iritype.ne.0)
     *         call tbxxwarn
     *         (' Multiple DDL 1 and 2 related item definitions ')
             iritype = i
             call tbxxnlc(riname,name(1:lstrg))
C
C            Seek the matching function, may be in the same loop or not
C
             if(charnp_(rf(i),name,lstrg)) then
               if (irftype.ne.0)
     *           call tbxxwarn
     *           (' Multiple DDL 1 and 2 related item functions ')
               irftype = i
               call tbxxnlc(rfname,name(1:lstrg))
             endif
           endif
         enddo
         loop_ = .false.
         loopnl = 0
C
C        Process categories
C
         do i = 1,5
           if(charnp_(ct(i),name,lstrg)) then
             if(i.eq.5) then
C
C              if this is a DDLm _defintion.scope with a value of
C              category, we need to get the name from _defintion.id
C
               call tbxxnlc(bcname,name(1:lstrg))
               if(bcname.eq.'category') then
                 if(.not.charnp_(nt(3),name,lstrg)) then
                   call tbxxwarn(
     *             ' DDLm category defintion without _definition.id ')
                 else
                   go to 216
                 endif
               endif
             endif

             if(ictype.ne.0)
     *         call tbxxwarn(
     *           ' Multiple DDL 1, 2 or m category definitions ')
             ictype = i
             if(loop_) icloop = loopnl
             call tbxxnlc(bcname,name(1:lstrg))
             lbcname=long_
             nmycat = ndcname+1
             call hash_store(bcname(1:long_),
     *         dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
             if(mycat.eq.0) then
               call tbxxerr(' Dictionary category names > NUMDICT ')
             endif
             if (mycat.eq.nmycat) then
               ccatkey(mycat) = 0
               xmcind(mycat)=0
             endif
C
C            if this is not a loop of categories, we expect a match
C            against the block name, unless we are doing replacements
C
             if(.not.loop_) then
               if(ictype.eq.1) then
                 if(bname(1:min(lbname,lbcname+2)).ne.
     *            '_'//bcname(1:lbcname)//'.'
     *            .and. catchk.eq.'yes'
     *            .and. (rfname(1:7).ne.'replace')) then
                 call tbxxwarn(' Category id does not match block name')
                 endif
               else
                 if(ictype.eq.2) then
                   if(bcname.ne.'dictionary_definition' .and.
     *                bcname.ne.'category_overview') then
                   if(bname(1:min(lbname,lbcname+2)).ne.
     *               '_'//bcname(1:lbcname)//'_') then
                   if(bname(1:min(lbname,lbcname+1)).ne.
     *               '_'//bcname(1:lbcname)
     *            .and. catchk.eq.'yes'
     *            .and. (rfname(1:7).ne.'replace')) then
                   call tbxxwarn(
     *               ' Category id does not match block name')
                   endif
                   endif
                   endif
                 endif
               endif
             endif
           endif
           loop_ = .false.
           loopnl = 0
         enddo
C
C        Process XML translations
C
216      loop_ = .false.
         loopnl = 0
         if(charnp_('_xml_mapping.token',xmtoken,lxmtoken)) then
230        if(charnp_('_xml_mapping.token_type',xmtyp,lxmtyp)) then
             if(charnp_('_xml_mapping.target',xmtarg,lxmtarg)) then
               if (xmnxlat.ge.XMLDEFS) then
                 call tbxxerr(' XML translations > XMLDEFS')
               else
                 xmnxlat=xmnxlat+1
                 xmlate(xmnxlat)=xmtarg(1:lxmtarg)
               endif
               if (xmtyp.eq.'data') then
                 ixmtyp = 1
                 if (xmdata.eq.0) then
                   xmdata = xmnxlat
                 else
                   call tbxxwarn(' XML duplicate DATA_ translation')
                 endif
               endif
               if (xmtyp(1:lxmtyp).eq.'category') then
                 ixmtyp = 2
                 nxmc = ndcname+1
                 call tbxxnlc(xxxtemp,xmtoken(1:lxmtoken))
                 call hash_store(xxxtemp,
     *           dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,kxmc)
                 if( kxmc.eq.nxmc) then
                   ccatkey(kxmc) = 0
                   xmcind(kxmc) = xmnxlat
                 else
                   if (xmcind(kxmc).ne.0) then
                     call tbxxwarn(
     *                 ' XML duplicate category translation')
                   else
                     xmcind(kxmc) = xmnxlat
                   endif
                 endif
               endif
               if (xmtyp.eq.'item') then
                 ixmtyp = 3
                 if (tbxxnewd(xmtoken(1:lxmtoken),ifind)) then
                   xmindex(ifind) = xmnxlat
                 else
                   if (xmindex(ifind).ne.0) then
                     call tbxxwarn(' XML duplicate item translation')
                   else
                     ifind = aroot(ifind)
 235                 xmindex(ifind) = xmnxlat
                     ifind = alias(ifind)
                     if (ifind.ne.0) go to 235
                   endif
                 endif
               endif
               if(loop_) then
                 if(charnp_('_xml_mapping.token',xmtoken,lxmtoken)) then
                   go to 230
                 else
                   call tbxxerr(' XML dictionary logic error')
                 endif
               endif
             else
               call tbxxerr(' XML target missing')
             endif
           else
             call tbxxerr(' XML token_type missing')
           endif
         else
           xmtoken = bname(1:lbname)
           lxmtoken=lbname
           if(charnp_('_xml_mapping.token_type',xmtyp,lxmtyp)) then
             if(charnp_('_xml_mapping.target',xmtarg,lxmtarg)) then
               if (xmnxlat.ge.XMLDEFS) then
                 call tbxxerr(' XML translations > XMLDEFS')
               else
                 xmnxlat=xmnxlat+1
                 xmlate(xmnxlat)=xmtarg(1:lxmtarg)
               endif
               if (xmtyp(1:lxmtyp).eq.'data') then
                 ixmtyp = 1
                 if (xmdata.eq.0) then
                   xmdata = xmnxlat
                 else
                   call tbxxwarn(' XML duplicate DATA_ translation')
                 endif
               endif
               if (xmtyp.eq.'category') then
                 ixmtyp = 2
                 nxmc = ndcname+1
                 call tbxxnlc(xxxtemp,xmtoken(1:lxmtoken))
                 call hash_store(xxxtemp,
     *           dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,kxmc)
                 if( kxmc.eq.nxmc) then
                   ccatkey(kxmc) = 0
                   xmcind(kxmc) = xmnxlat
                 else
                   if (xmcind(kxmc).ne.0) then
                     call tbxxwarn(
     *                 ' XML duplicate category translation')
                   else
                     xmcind(kxmc) = xmnxlat
                   endif
                 endif
               endif
               if (xmtyp.eq.'item') then
                 ixmtyp = 3
                 if (tbxxnewd(xmtoken(1:lxmtoken),ifind)) then
                   xmindex(ifind) = xmnxlat
                 else
                   if (xmindex(ifind).ne.0) then
                     call tbxxwarn(' XML duplicate item translation')
                   else
                     ifind = aroot(ifind)
 240                 xmindex(ifind) = xmnxlat
                     ifind = alias(ifind)
                     if (ifind.ne.0) go to 240
                     xmindex(ifind) = xmnxlat
                   endif
                 endif
               endif
               if(loop_) then
                 call tbxxerr(' XML dictionary logic error')
               endif
             else
               call tbxxerr(' XML target missing')
             endif
           endif
         endif
C
C        Process names
C
         bxname = ' '
         do i = 1,3
         if(charnp_(nt(i),name,lstrg)) then
           if(intype.ne.0)
     *       call tbxxwarn(
     *         ' Multiple DDL 1 and 2 or m name definitions ')
           intype = i
           call tbxxnlc(bxname,name(1:lstrg))
           if(loop_) inloop = loopnl
         endif
         loop_ = .false.
         loopnl=0
         enddo
         if(intype.eq.0.and.ictype.lt.3.and.(.not.glob_)
     *     .and.bname(1:lbname).ne.' '.and.ixmtyp.eq.0)
     *     call tbxxwarn (' No name defined in block')
         loop_ = .false.
         if(charnp_(at(1),name,lstrg)) then
           iatype=1
           call tbxxnlc(baname,name(1:lstrg))
           batag = name(1:lstrg)
           lbaname = lstrg
           if(loop_) ialoop = loopnl
         endif
         loop_ = .false.
         loopnl=0
         mcstrg = "no"
         if(ictype.ne.3) then
           do i=1,3
             if(charnp_(tt(i),name,lstrg)) then
               if(ittype.ne.0)
     *           call tbxxwarn(
     *             ' Multiple DDL 1 and 2 type definitions ')
               ittype = i
               call tbxxnlc(btname,name(1:lstrg))
               if(loop_) itloop = loopnl
             endif
             loop_ = .false.
             loopnl=0
           enddo
           do i = 1,2
             if(charnp_(mc(i),name,lstrg)) then
               if (imtype.ne.0)
     *           call tbxxwarn(' Multiple DDL 1 and 2 mandatory codes ')
               imtype = i
               call tbxxnlc(mcstrg,name(1:lstrg))
               if (loop_) imloop = loopnl
             endif
             loop_ = .false.
             loopnl=0
           enddo
         endif
C
C        Now test for consistent combinations
C
         if(inloop.ne.-1) then
           if(icloop.ne.-1.and.icloop.ne.inloop
     *            .and. catchk.eq.'yes')
     *       call tbxxwarn(
     *       ' Categories and names in different loops')
           if(iatype.ne.0.and.ialoop.ne.inloop) then
             if(ialoop.eq.-1) then
               if(bxname.ne.bname(1:lbname))
     *          call tbxxwarn(
     *         ' One alias, looped names, linking to first')
             else
               call tbxxwarn(
     *         ' Aliases and names in different loops '
     *         //' only using first alias ')
             endif
           endif
           if(itloop.ne.-1.and.itloop.ne.inloop)
     *       call tbxxwarn(
     *       ' Types and names in different loops')
           if(imloop.ne.-1.and.imloop.ne.inloop)
     *       call tbxxwarn(
     *       ' Mandatory codes and names in different loops')
         else
           if(icloop.ne.-1)
     *       call tbxxwarn(
     *         ' Multiple categories for one name')
           if(itloop.ne.-1)
     *       call tbxxwarn(
     *         ' Multiple types for one name')
           if(imloop.ne.-1)
     *       call tbxxwarn(
     *         ' Multiple madatory codes for one name')
         endif
C
C        Pick up parents
C
         do i = 1,2
220        if(charnp_(pt(i),name,lstrg)) then
             if (iptype.ne.0 .and. iptype.ne.i)
     *         call tbxxwarn
     *         (' Multiple DDL 1 and 2 parent definitions ')
             iptype = i
             call tbxxnlc(bpname,name(1:lstrg))
             lbpname=long_
C
C            Seek the matching child, may be in the same loop or not
C
             if (charnp_(pc(i),name,lstrg)) then
               nresult = tbxxnewd(name(1:lstrg),ifind)
               nresult = tbxxnewd(bpname(1:lbpname),dpindex(ifind))
               bpname = ' '
               lbpname = 1
             endif
             if (loop_) go to 220
           endif
         enddo

C
C        Now we need to process value enumerations and ranges
C        and load them into item value table
C
         if (tcheck .eq. 'yes' .and. bxname.ne.' ') then
         loop_ = .false.
         nresult = tbxxnewd(bxname,ifind)

         do i = 1,2
5400       if(charnp_(ve(i),name,lstrg) .and. nivt.lt.NUMIVALS) then
             call tbxxsstb(name(1:lstrg),sindex)
             if (sindex.gt.0) then
               if (deindex(ifind).eq.0) then
                 deindex(ifind)=nivt+1
               else
                 kivt = deindex(ifind)
5410             if (ivtnxt(kivt).ne.0) then
                   kivt = ivtnxt(kivt)
                   go to 5410
                 endif
                 ivtnxt(kivt)=nivt+1
               endif
               nivt = nivt+1
               ivtnxt(nivt)=0
               ivtvet(nivt)=0
               ivtsbp(nivt)=sindex
             endif
           endif
           if (loop_) go to 5400
         enddo


         do i = 1,2
         loop_ = .false.
5420     strg_=' '
         long_=1
         nresult = test_(vr(i))
         if (strg_(1:long_).ne.' '.and.type_.eq.'null')
     *     nresult = .true.
         if (nresult .and. nivt.lt.NUMIVALS) then
           nresult = charnp_(vr(i),name,lstrg)
           if (type_.ne.'char'.and.type_.ne.'numb') then
             name = '.'
             lstrg = 1
           endif
           kvrtp = -1
           if(i.eq.1 .and. lstrg<len(name)-2) then
             strg_=' '
             long_=1
             nresult = test_(vr(3))
             if (strg_(1:long_).ne.' '.and.
     *         type_.eq.'null') nresult = .true.
             if (nresult) then
               nresult = charnp_(vr(3),
     *           name(lstrg+2:len(name)),kstrg)
               if (type_.ne.'char'.and.type_.ne.'numb') then
                 name(lstrg+2:len(name)) = '.'
                 kstrg = 1
               endif
               if (name(1:lstrg).ne.name(lstrg+2:lstrg+1+kstrg))
     *           kvrtp = 1
               name(lstrg+1:lstrg+1)=':'
               lstrg = lstrg+kstrg+1
               endif
             endif
             if (name(1:lstrg).eq.'.:.') then
               loop_=.false.
             else
               call tbxxsstb(name(1:lstrg),sindex)
               if (sindex.gt.0) then
                 if (deindex(ifind).eq.0) then
                   deindex(ifind)=nivt+1
                 else
                   kivt = deindex(ifind)
5430               if (ivtnxt(kivt).ne.0) then
                     kivt = ivtnxt(kivt)
                     go to 5430
                   endif
                   ivtnxt(kivt)=nivt+1
                 endif
                 nivt = nivt+1
                 ivtnxt(nivt)=0
                 ivtvet(nivt)=kvrtp
                 ivtsbp(nivt)=sindex
               endif
               if(loop_) go to 5420
             endif
           endif
         enddo
         endif

C
C        This is the main loop
C
250      if(ictype.eq.5.or.intype.eq.0) goto 200
         if(.not.charnp_(nt(intype),name,lstrg)) goto 200
         kdict=ndict+1
251      nresult = tbxxnewd(name(1:lstrg),ifind)
         if (bpname .ne. ' ') then
           nresult=tbxxnewd(bpname(1:lbpname),dpindex(ifind))
           bpname = ' '
           lbpname = 1
         endif
         nresult = tbxxoldd(bxname,jfind)
         if (nresult.and.jfind.ne.ifind.and.deindex(ifind).eq.0)
     *     deindex(ifind) = deindex(jfind)
         if(ifind.le.idstrt) then
           if (kdup .lt. 0) then
             call tbxxerr(' Duplicate name in dictionary '//
     *       dictag(ifind)(1:lastnb(dictag(ifind))))
           endif
           if (kdup .gt.0) go to 254
           dicnam(ifind)=char(0)
           goto 251
254        continue
         endif
         if(dicnam(ifind).eq.bname(1:lbname)) nmatch=ifind
         if(dicnam(ifind)(1:lbname).eq.bname(1:lbname)) ksmatch=ifind
Cdbg     if(dicnam(ifind).ne.bname(1:lbname))
Cdbg *   call tbxxwarn (' Name mismatch: '//dicnam(ifind)//bname(1:lbname))
         if(inloop.ge.0)then
C
C          We are in a loop of names.  If it is the same loop as
C          for categories, we need to extract the matching category
C
           if(inloop.eq.icloop) then
             mycat=0
             if(charnp_(ct(ictype),name,lstrg)) then
               call tbxxnlc(bcname,name(1:lstrg))
               lbcname=lstrg
               nmycat=ndcname+1
               call hash_store(bcname,
     *         dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
               if(mycat.eq.0) then
                 call tbxxerr(' Dictionary category names > NUMDICT ')
               endif
               if(mycat.eq.nmycat) ccatkey(mycat)=0
             endif
           endif
C
C          If it is the same loop as for types, we need to extract
C          the matching type
C
           if(inloop.eq.itloop) then
             btname=' '
             if(charnp_(ct(ittype),name,lstrg)) then
               call tbxxnlc(btname,name(1:lstrg))
             endif
           endif
C
C          If it is the same loop as for mandatory codes, we need to extract
C          the matching mandatory
C
           if(inloop.eq.imloop) then
             mcstrg='no'
             if(charnp_(mc(imtype),name,lstrg)) then
               call tbxxnlc(mcstrg,name(1:lstrg))
             endif
           endif
C
C          If it is the same loop as for aliases, we need to extract
C          the matching alias
C
           if(inloop.eq.ialoop) then
             baname=' '
             batag=' '
             if(charnp_(at(1),name,lstrg)) then
               call tbxxnlc(baname,name(1:lstrg))
               batag = name(1:lstrg)
               lbaname = lstrg
             endif
           endif
         endif
C
C        now we have a name stored in dicnam at location ifind
C        the index of the category in mycat, the type in btname,
C        the alias in baname, and the mandatory code in mcstrg
C
C        First verify match between the name and category, if
C        we have one, or extract from the block name
C
         if (mycat.eq.0) then
         if (dcindex(ifind).eq.0) then
           if (dicnam(ifind).eq.bloc_) then
             call tbxxcat(dicnam(ifind),bcname,lbcname)
Cdbg         call tbxxwarn(' Extracting category name from block name '
Cdbg *       //bloc_(1:max(1,lastnb(bloc_))))
             if(bcname(1:1).ne.' ') then
               ictype = 1
               nmycat = ndcname+1
               call hash_store(bcname,
     *         dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
               if(mycat.eq.0) then
                 call tbxxerr(' Dictionary category names > NUMDICT ')
               endif
               if (mycat.eq.nmycat) then
                 ccatkey(mycat) = 0
                 xmcind(mycat) = 0
               endif
             else
               if(catchk.eq.'yes')
     *         call tbxxwarn(' No category defined in block '
     *       //bloc_(1:max(1,lastnb(bloc_)))//' and name '
     *       //dicnam(ifind)(1:max(1,lastnb(dicnam(ifind))))
     *       //' does not match')
             endif
           endif
         endif
         else
         if (bcname(1:lbcname).ne.'dictionary_definition' .and.
     *     bcname(1:lbcname).ne.'category_overview') then
           if (dicnam(ifind)(1:lbcname+1).ne.'_'//bcname(1:lbcname)
     *        .or.( dicnam(ifind)(lbcname+2:lbcname+2).ne.'_' .and.
     *          dicnam(ifind)(lbcname+2:lbcname+2).ne.'.' .and.
     *          dicnam(ifind)(lbcname+2:lbcname+2).ne.' ' )) then
                if (catchk.eq.'yes'.and.rfname(1:7).ne.'replace')
     *          call tbxxwarn(' Item name '//
     *          dicnam(ifind)(1:max(1,lastnb(dicnam(ifind))))//' '//
     *       ' does not match category name '//bcname(1:lbcname))
           endif
         endif
         endif
C
C        We will need the type in what follows.  cif_mm.dic defines
C        some higher level types.  We map them to primitive types
C
         mapped = btname(1:4)
         do i = 1,19
           if (btname(1:4).eq.map_type(i)) mapped = map_to(i)
         enddo
         if (mapped.ne.'char' .and.
     *       mapped.ne.'text' .and.
     *       mapped.ne.'null' .and.
     *       mapped.ne.'numb' .and.
     *       mapped.ne.'    ' ) then
             if (tcheck .eq. 'yes') then
               call tbxxwarn (' Item type '//
     *           btname(1:max(1,lastnb(btname)))//' not recognized')
             endif
             mapped = 'char'
         endif

C
C        There are two cases to consider, one if the name is new to
C        the dictionary, the other, if it is not
C
         if(ifind.eq.kdict) then
           aroot(ifind)=ifind
           alias(ifind)=0
           dcindex(ifind)=mycat
           dictyp(ifind)=mapped
           dicxtyp(ifind)=btname
           dmcode(ifind) = 0
           if (mcstrg .eq. 'yes') dmcode(ifind) = 1
           if (mcstrg .eq. 'implicit') dmcode(ifind) = -1
         else
           if(dcindex(ifind).ne.mycat) then
             if(dcindex(ifind).eq.0) then
               jfind=ifind
               if (aroot(ifind).ne.0) jfind=aroot(ifind)
255            continue
               dcindex(jfind)=mycat
               jfind=alias(jfind)
               if(jfind.ne.0) goto 255
             else
               if(mycat.ne.0.and.
     *           (vcheck.eq.'yes'.or.tcheck.eq.'yes')
     *           .and.catchk.eq.'yes')  then
                 if(rfname(1:7).ne.'replace')
     *           call tbxxwarn(' Attempt to redefine category for item')
                 endif
             endif
           endif
           if(dictyp(ifind).ne.mapped .or.
     *       dicxtyp(ifind).ne.btname) then
             if(dictyp(ifind).eq.' ') then
               jfind=ifind
               if (aroot(ifind).ne.0) jfind=aroot(ifind)
256            continue
               dictyp(jfind)=mapped
               dicxtyp(jfind)=btname
               jfind=alias(jfind)
               if(jfind.ne.0) go to 256
             else
               if(mapped.ne.' '.and.tcheck.eq.'yes')
     *           call tbxxwarn(' Attempt to redefine type for item')
             endif
           endif
           if(dmcode(ifind).eq.0) then
             jfind = ifind
             if (aroot(ifind).ne.0) jfind = aroot(ifind)
257          continue
             dmcode(jfind) = 0
             if (mcstrg.eq.'yes') dmcode(jfind) = 1
             if (mcstrg.eq.'implicit') dmcode(jfind) = -1
             jfind=alias(jfind)
             if(jfind.ne.0) go to 257
           else
             if((mcstrg.eq.'yes' .and. dmcode(ifind).lt.0) .or.
     *         (mcstrg.eq.'implicit' .and. dmcode(ifind).gt.0))
     *         call tbxxwarn(
     *           ' Attempt to redefine mandatory code for item')
           endif
         endif
C
C        now deal with alias, if any.
C
         if(baname.ne.' ') then
           if (tbxxnewd(baname(1:lbaname),iafind)) then
             dictag(iafind)    =batag
             aroot(iafind)     =aroot(ifind)
             if(aroot(iafind).eq.0) aroot(iafind)=ifind
             catkey(iafind)    =catkey(ifind)
             alias(ifind)      =iafind
             dcindex(iafind)   =dcindex(ifind)
             dictyp(iafind)    =dictyp(ifind)
             dicxtyp(iafind)   =dicxtyp(ifind)
             xmindex(iafind)   =xmindex(ifind)
             dmcode(iafind)    =dmcode(ifind)
             dpindex(iafind)   =dpindex(ifind)
             deindex(iafind)   =deindex(ifind)
           else
             if(aroot(iafind).ne.0 .and.
     *         aroot(iafind).ne.iafind) then
               if(aroot(iafind).eq.ifind .or.
     *           aroot(iafind).eq.aroot(ifind)) then
                 call tbxxwarn(' Duplicate definition of same alias')
               else
                 call tbxxwarn(' Conflicting definition of alias')
               endif
             else
               if((dcindex(iafind).eq.0.or.
     *           dcindex(iafind).eq.dcindex(ifind)).and.
     *           (dictyp(iafind).eq.' '.or.
     *           (dictyp(iafind).eq.dictyp(ifind) .and.
     *            dicxtyp(iafind).eq.dicxtyp(ifind)))) then
                 dcindex(iafind)   =dcindex(ifind)
                 dictyp(iafind)    =dictyp(ifind)
                 dicxtyp(iafind)   =dicxtyp(ifind)
               endif
               if(xmindex(iafind).eq.0)
     *           xmindex(iafind)=xmindex(ifind)
               if(xmindex(ifind).eq.0)
     *           xmindex(ifind)=xmindex(iafind)
               if (dmcode(iafind).eq.0)
     *           dmcode(iafind)=dmcode(ifind)
               if (dmcode(ifind).eq.0)
     *           dmcode(ifind)=dmcode(iafind)
               if (dpindex(iafind).eq.iafind
     *           .and. dpindex(ifind).ne.ifind)
     *           dpindex(iafind) = dpindex(ifind)
               if (dpindex(ifind).eq.ifind
     *           .and. dpindex(iafind).ne.iafind)
     *           dpindex(ifind) = dpindex(iafind)
               if (deindex(ifind).eq.0)
     *           deindex(ifind)=deindex(iafind)
               if (deindex(iafind).eq.0)
     *           deindex(iafind)=deindex(ifind)
               aroot(iafind)     =aroot(ifind)
               if(aroot(iafind).eq.0) aroot(iafind)=ifind
               alias(ifind)      =iafind
               if (catkey(iafind)) catkey(ifind) = .true.
               if (catkey(ifind)) catkey(iafind) = .true.
             endif
           endif
         endif
         if(inloop.ge.0) then
           baname = ' '
           batag = ' '
         endif
C
         if(inloop.ge.0.and.loop_) go to 250
         if(nmatch.eq.0) then
         if ((ksmatch.eq.0.or.inloop.lt.0)
     *     .and.(rfname(1:7).ne.'replace')) then
         call tbxxwarn(' No name in the block matches the block name')
         endif
         endif
C
C        check for aliases
C        we execute this loop only in the case of unlooped name
C        with looped alias
C
         if(inloop.lt.0.and.ialoop.ge.0) then
           loop_=.false.
           loopnl=0
           ganame=baname
260        if(.not.charnp_(at(iatype),name,lstrg)) goto 200
           call tbxxnlc(baname,name(1:lstrg))
           batag=name(1:lstrg)
           lbaname=lstrg
           if(baname.eq.ganame) then
             if(loop_) go to 260
             go to 200
           endif
           if(baname.ne.' ') then
             if (tbxxnewd(baname(1:lbaname),iafind)) then
             if(iafind.eq.0) call tbxxerr(' CIFdic names > NUMDICT')
               dictag(iafind)    =batag
               aroot(iafind)     =aroot(ifind)
               if(aroot(iafind).eq.0) aroot(iafind)=ifind
               catkey(iafind)    =catkey(ifind)
               alias(ifind)      =iafind
               dcindex(iafind)   =dcindex(ifind)
               dictyp(iafind)    =dictyp(ifind)
               dicxtyp(iafind)   =dicxtyp(ifind)
               xmindex(iafind)   =xmindex(ifind)
               dmcode(iafind)    =dmcode(ifind)
               dpindex(iafind)   =dpindex(ifind)
               deindex(iafind)   =deindex(ifind)
               ifind=iafind
             else
               if(aroot(iafind).ne.0 .and.
     *           aroot(iafind).ne.iafind) then
                 if(aroot(iafind).eq.ifind .or.
     *             aroot(iafind).eq.aroot(ifind)) then
                   call tbxxwarn(' Duplicate definition of same alias')
                 else
                   call tbxxwarn(' Conflicting definition of alias')
                 endif
               else
                 if((dcindex(iafind).eq.0.or.
     *           dcindex(iafind).eq.dcindex(ifind)).and.
     *           (dictyp(iafind).eq.' '.or.
     *           (dictyp(iafind).eq.dictyp(ifind) .and.
     *            dicxtyp(iafind).eq.dicxtyp(ifind)))) then
                 dcindex(iafind)   =dcindex(ifind)
                 dictyp(iafind)    =dictyp(ifind)
                 dicxtyp(iafind)   =dicxtyp(ifind)
                 ifind=iafind
                 endif
                 if(xmindex(iafind).eq.0)
     *             xmindex(iafind)=xmindex(ifind)
                 if(xmindex(ifind).eq.0)
     *             xmindex(ifind)=xmindex(iafind)
                 if (dmcode(iafind).eq.0)
     *             dmcode(iafind)=dmcode(ifind)
                 if (dmcode(ifind).eq.0)
     *             dmcode(ifind)=dmcode(iafind)
                 if (dpindex(iafind).eq.iafind
     *             .and. dpindex(ifind).ne.ifind)
     *             dpindex(iafind) = dpindex(ifind)
                 if (dpindex(ifind).eq.ifind
     *             .and. dpindex(iafind).ne.iafind)
     *             dpindex(ifind) = dpindex(iafind)
                 if (deindex(ifind).eq.0)
     *             deindex(ifind) = deindex(iafind)
                 if (deindex(iafind).eq.0)
     *             deindex(iafind) = deindex(ifind)
                 aroot(iafind)     =aroot(ifind)
                 if(aroot(iafind).eq.0) aroot(iafind)=ifind
                 alias(ifind)      =iafind
                 if (catkey(iafind)) catkey(ifind) = .true.
                 if (catkey(ifind)) catkey(iafind) = .true.
               endif
             endif
           endif
           if(loop_) go to 260
         endif
         go to 200
C
400      bloc_=' '
         if (ndcname.ne.0) then
         do ii = idstrt+1,ndict
         keychain(ii) = 0
         if (aroot(ii).eq.0.and.dcindex(ii).eq.0
     *     .and.catchk.eq.'yes')
     *     call tbxxwarn(' No category specified for name '//
     *       dicnam(ii)(1:max(1,lastnb(dicnam(ii)))))
         enddo
         endif
         do ii = idstrt+1,ndict
         if (dicxtyp(ii).eq.' ') then
           if (dpindex(ii).ne.ii
     *       .and. dicxtyp(dpindex(ii)).ne.' ') then
             dicxtyp(ii) = dicxtyp(dpindex(ii))
             dictyp(ii) = dicxtyp(dpindex(ii))(1:4)
           else
             dicxtyp(ii) = 'null'
             dictyp(ii) = 'null'
             if (tcheck.eq.'yes')  then
               jj = lastnb(dicnam(ii))
               if (jj.gt.0) then
               if (dicnam(ii)(jj:jj).ne.'_')
     *         call tbxxwarn(' No type specified for name '//
     *           dicnam(ii)(1:max(1,lastnb(dicnam(ii)))))
               endif
             endif
           endif
         endif
         if (catkey(ii) .or. dmcode(ii).gt.0) then
           ifind = aroot(ii)
           mycat = dcindex(ifind)
           if (mycat.ne.0) then
             jj = ccatkey(mycat)
             if (jj.eq.0) then
               ccatkey(mycat) = ifind
             else
410            if (keychain(jj).eq.0) then
                 keychain(jj) = ifind
                 keychain(ifind) = 0
               else
                 if(keychain(jj).ne.ifind) then
                   jj = keychain(jj)
                   goto 410
                 endif
               endif
             endif
           endif
         endif
         enddo
         if (.not.append_) then
           close(dirdev)
           nrecd=0
         endif
         dictfl='no '
500      continue
         if (append_) then
           nrecd=nrecds
           recend_=recends
           recbeg_=recbegs
         endif
         if(dict_) then
           dicname_=xdicnam
           dicver_ =xdicver
         else
           tcheck = otchk
           vcheck = ovchk
         endif
         if(tcheck.eq.'yes') vcheck='yes'
Cdbg     WRITE(6,'(i5,3x,a,2x,a)') (i,dicnam(i),dictyp(i),i=1,ndict)
         return
         end
C
C
C
C
C
C >>>>>> Create a new dictionary entry, or find a matching existing one
C
         function tbxxnewd(xname,ick)
         logical   tbxxnewd
         include  'ciftbx.sys'
         character xname*(*)
         character xxxtemp*(NUMCHAR)
         integer   jck, ick, ilen
         integer   lastnb
         tbxxnewd = .true.
         ilen = lastnb(xname)
         jck = ndict
         call tbxxnlc(xxxtemp,xname(1:ilen))
         call hash_store(xxxtemp,
     *     dicnam,dicchain,
     *     NUMDICT,ndict,dichash,NUMHASH,ick)
         if(ick.eq.0) call tbxxerr(' CIFdic names > NUMDICT')
         if(ick .eq. jck+1) then
           dictag(ick) = xname(1:ilen)
           dictyp(ick) = ' '
           dicxtyp(ick) = ' '
           catkey(ick) = .false.
           dpindex(ick) = ick
           deindex(ick) = 0
           alias(ick) = 0
           aroot(ick) = ick
           keychain(ick) = 0
           dcindex(ick) = 0
           xmindex(ick) = 0
           dmcode(ick) = 0
         else
           tbxxnewd = .false.
         endif
         return
         end
C
C
C
C
C
C >>>>>> Find matching existing dictionary entry if any
C
         function tbxxoldd(xname,ick)
         logical   tbxxoldd
         include  'ciftbx.sys'
         character xname*(*)
         character xxxtemp*(NUMCHAR)
         integer   ick, ilen
         integer   lastnb
         tbxxoldd = .true.
         ilen = lastnb(xname)
         call tbxxnlc(xxxtemp,xname(1:ilen))
         call hash_find(xxxtemp,
     *     dicnam,dicchain,
     *     NUMDICT,ndict,dichash,NUMHASH,ick)
         if(ick.eq.0) tbxxoldd = .false.
         return
         end
C
C
C
C
C
C >>>>>> Find position of last non_blank in a string
C        but never less than 1
C
         function lastnb(str)
C
         integer    lastnb
         include   'ciftbx.sys'
         character*(*) str
         integer lenn,ihi,itestl
         lenn = len(str)
c
         ihi = lenn
         if(ihi.eq.0) then
           ihi = 1
           go to 200
         endif
         itestl = ihi/4
         if (itestl.lt.4) go to 200
c
100      if (ihi.gt.itestl) then
         if (str(ihi-itestl+1:ihi-itestl+1).eq.' ') then
           if (str(ihi-itestl+1:ihi).eq.' ') then
             ihi = ihi-itestl
             go to 100
           endif
         endif
         endif
         itestl = itestl/2
         if (itestl.gt.3) go to 100
c
200      if (ihi.gt.1 .and. str(ihi:ihi).eq.' ') then
           ihi = ihi-1
           go to 200
         endif
         if (ihi.eq.0) ihi = 1
         lastnb = ihi
         return
         end
C
C
C
C
C
C >>>>>> Convert a character to a radix XXRADIX digit
C
C        given a character c, return a decimal value
C
         function tbxxc2dig(c)
         integer   tbxxc2dig
         character*(*) c
         include  'ciftbx.sys'
C
         tbxxc2dig = ichar(c)-ichar(' ')
C
C        The code above may not be portable, especially to non-ascii
C        computer systems.  In that case, comment out the line above
C        and uncomment the following lines.  Be sure to make the
C        matching change in tbxxd2chr.  Be certain to have at least
C        XXRADIX characters in the search string.
C
C         tbxxc2dig = index(
C     *   '+-01234567890'//
C     *   'abcdefghijlmnopqrstuvwxyz'//
C     *   'ABCDEFGHIJKLMNOPQRSTUVWXYZ',c)-1
         return
         end
C
C
C
C
C
C >>>>>> Convert a radix XXRADIX digit to a character
C
C        given an integer value, return a character
C
         function tbxxd2chr(d)
         character*1 tbxxd2chr
         integer   d
         include  'ciftbx.sys'
C
         tbxxd2chr = char(d+ichar(' '))
C
C        The code above may not be portable, especially to non-ascii
C        computer systems.  In that case, comment out the line above
C        and uncomment the following lines.  Be sure to make the
C        matching change in tbxxc2dig.  Be certain to have at least
C        XXRADIX characters in the search string.
C
C         character*(XXRADIX) digits
C         digits =
C     *   '+-01234567890'//
C     *   'abcdefghijlmnopqrstuvwxyz'//
C     *   'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
C         tbxxd2chr = digits(d+1:d+1)
         return
         end
C
C
C
C
C
C >>>>>> Convert a string to Run Length Encoded version
C
         subroutine tbxxrle(astr,bstr,mlen)
C
C        astr is the raw input string
C        bstr is the run-length-encoded string
C          beginning with the compressed length in
C            in base-XXRADIX in the first four characters
C          followed by either individual characters or run
C          flagged by XXFLAG
C        XXFLAG//tbxxd2chr(n)//c represents n copies of c
C
         character*(*) astr, bstr
         include  'ciftbx.sys'
         character*1 c
         character*1 tbxxd2chr
         integer tbxxc2dig
         integer klen, krep, ialen, iblen, mode, ii
         integer mlen
C
         ialen = len(astr)
         iblen = len(bstr)
         mode = 0
         klen = 4
         bstr(1:4) = tbxxd2chr(0)//tbxxd2chr(0)
     *     //tbxxd2chr(0)//tbxxd2chr(0)
         do ii = 1,ialen
           c = astr(ii:ii)
           if (mode .eq. -2) then
             krep = tbxxc2dig(bstr(klen-1:klen-1))
             if (c.eq.bstr(klen:klen).and.krep.lt.XXRADIX-1) then
               bstr(klen-1:klen-1) = tbxxd2chr(krep+1)
             else
               mode = 0
               if (c.eq.bstr(klen:klen)) mode=-1
             endif
           endif
           if (klen.ge.iblen) go to 100
           if (mode .ge.-1 .and. mode .le.2) then
             klen = klen+1
             bstr(klen:klen) = c
             if (klen .gt. 5) then
               if (c.eq.bstr(klen-1:klen-1)) mode=mode+1
               if (c.ne.bstr(klen-1:klen-1)) mode=0
             endif
             if (c.eq.XXFLAG .and. klen.lt.iblen-1) then
               bstr(klen+1:klen+2) = tbxxd2chr(1)//c
               mode = -2
               klen = klen+2
             endif
           endif
           if (mode.eq.2) then
             bstr(klen-2:klen-1) = XXFLAG//tbxxd2chr(3)
             mode = -2
           endif
         enddo
 100     mlen = klen
         do ii = 4,1,-1
           bstr(ii:ii) = tbxxd2chr(mod(klen,XXRADIX))
           klen = klen/XXRADIX
         enddo
         return
         end
C
C
C
C
C
C >>>>>> Decode a string from  Run Length Encoded version
C
         function tbxxrld(astr,bstr,fill)
C
C        astr is the raw output string
C        bstr is the run-length-encoded string
C          beginning with the compressed length in
C            in base-XXRADIX in the first four characters
C          followed by either individual characters or run
C          flagged by char(0)
C        char(0)//char(n)//c represents n copies of c
C        fill is a logical variable, .true. to fill astr with blanks
C        the return value is the number of valid characters in astr
C        never less than 1
C
C
         character*(*) astr, bstr
         logical fill
         integer tbxxrld
         include  'ciftbx.sys'
         character*1 c
         integer tbxxc2dig
         integer klen, krep, ialen, mode, ipos
         integer ii, jj
C
         tbxxrld = 1
         krep = 0
         ialen = len(astr)
         if (fill) then
           astr = ' '
         else
           astr(1:1) = ' '
         endif
         mode = 0
         klen = 0
         do ii = 1,4
           klen = klen*XXRADIX+tbxxc2dig(bstr(ii:ii))
         enddo
         mode = 0
         ipos = 0
         do ii = 5,klen
           c = bstr(ii:ii)
           if(mode.eq.0) then
             if(c.ne.XXFLAG) then
               if (ipos.ge.ialen) then
                 tbxxrld = ialen
                 return
               endif
               ipos = ipos+1
               astr(ipos:ipos) = c
             else
               mode = 1
             endif
           else
             if (mode.eq.1) then
               krep = tbxxc2dig(c)
               mode = -1
             else
               do jj = 1,krep
                 if (ipos.ge.ialen) return
                 ipos=ipos+1
                 astr(ipos:ipos) = c
               enddo
               mode = 0
             endif
           endif
         enddo
         if(ipos .lt. ialen) astr(ipos+1:ipos+1) = ' '
         tbxxrld = max(ipos,1)
         return
         end
C
C
C
C
C
C >>>>>> Extract the item.category_id from a save frame name
C
         subroutine tbxxcat(sfname,bcname,lbcname)
C
         character*(*) sfname,bcname
         integer lbcname,ii,ic,lastnb,lenn
C
C        Note that this logic works only for item.category_id
C        not for category.id
C
         lenn = lastnb(sfname)
         bcname = ' '
         lbcname = 1
         if (lenn.eq.0.or.sfname(1:1).ne.'_') return
         do ii = 1,lenn-2
         ic = 1+lenn-ii
         if (sfname(ic:ic).eq.'.') then
           bcname = sfname(2:ic-1)
           lbcname = ic-2
           return
         endif
         enddo
         return
         end
C
C
C
C
C
C >>>>>> Fetch line from direct access file
C
         subroutine tbxxflin(linno,lip,lipag,lipof,kip,ip,mip,mis)
C
         include   'ciftbx.sys'
         integer    linno,lip,kip,ip,mip,mis,i,mipno,miprno, kzero
         integer    lipag,lipof,kmode
C
C        linno -- the line number to locate
C        lip   -- the location of the line
C                   (page*(NUMCPP/NUMCIP)+offset)
C        lipag -- the page number (1...)
C        lipof -- the offset (1...)
C        kip   -- subindex number
C        ip    -- subindex offset
C        mip   -- master index number
C        mis   -- master index offset

         kip = (linno-1)/NUMSIP + 1
         ip = mod(linno-1,NUMSIP) + 1
         mip = (kip-1)/NUMMIP + 1
         mis = mod(kip-1,NUMMIP) + 1
C
C        test subindex page number against number in memory
C
         if (kip.ne.iabs(ipim)) then
C
C          save the current subindex page if it has been written
C
           if (ipim.lt.0) then
             do i = 1,NUMSIP
               write(scrbuf(NUMCIP*(i-1)+1:NUMCIP*i),'(i8)')
     *           ippoint(i)
             enddo
             write(dirdev,'(a)',rec=iabs(iprim)) scrbuf
             ipim = -ipim
           endif
C
C          find the appropriate master index page and slot
C
           if (mip.ne.iabs(mipim)) then
C
C            save the current master index page if it has been written
C
             if (mipim.lt.0) then
               write(scrbuf(1:NUMCIP),'(i8)')mipcp
               do i = 1,NUMMIP
                 write(scrbuf(NUMCIP*i+1:NUMCIP*(i+1)),'(i8)')
     *             mippoint(i)
               enddo
               write(dirdev,'(a)',rec=iabs(miprim))scrbuf
               mipim = -mipim
             endif
C
C            search the master index pages for a match
C
             mipno = 0
             miprno = 1
             kzero = 0
             kmode = 1
 10          read(dirdev,'(a)',rec=miprno) scrbuf
             mipno = mipno+1
             read(scrbuf(1:NUMCIP),'(i8)') mipcp
             if (mipno.ne.mip) then
               if (mipcp.eq.0) then
                 if (nfword.gt.1) then
                   nfblock = nfblock+1
                   nfword = 1
                 endif
                 mipcp = nfblock
                 nfblock = nfblock+1
                 write(scrbuf(1:NUMCIP),'(i8)') mipcp
                 write(dirdev,'(a)',rec=miprno) scrbuf
                 scrbuf = ' '
                 write(scrbuf(1:NUMCIP),'(i8)') kzero
                 write(dirdev,'(a)',rec=mipcp) scrbuf
                 kmode = -1
               endif
               miprno = mipcp
               go to 10
             endif
C
C            Have the master index in scrbuf, copy to mippoint
C
             do i = 1,NUMMIP
               if (scrbuf(NUMCIP*i+1:NUMCIP*(i+1)).eq.' ') then
                 mippoint(i) = 0
               else
                 read(scrbuf(NUMCIP*i+1:NUMCIP*(i+1)),'(i8)')
     *             mippoint(i)
               endif
             enddo
             mipim =kmode* mip
             miprim = miprno
           endif
C
C          See if the subindex page exists
C
           if (mippoint(mis).eq.0) then
             do i = 1,NUMSIP
               ippoint(i) = 0
             enddo
             if (nfword.gt.1) then
               nfblock=nfblock+1
               nfword = 1
             endif
             mippoint(mis) = nfblock
             mipim = -iabs(mipim)
             ipim = -kip
             iprim = -nfblock
             scrbuf = ' '
             write(dirdev,'(a)', rec=nfblock) scrbuf
             nfblock = nfblock+1
           else
             read(dirdev,'(a)', rec=mippoint(mis)) scrbuf
             do i = 1,NUMSIP
               if (scrbuf(NUMCIP*(i-1)+1:NUMCIP*i).eq.' ') then
                 ippoint(i) = 0
               else
               read(scrbuf(NUMCIP*(i-1)+1:NUMCIP*i),'(i8)')
     *           ippoint(i)
               endif
             enddo
             ipim = kip
             iprim = mippoint(mis)
           endif
         endif
         lip = ippoint(ip)
         lipag = (lip-1)/(NUMCPP/NUMCIP) + 1
         lipof = mod(lip-1,NUMCPP/NUMCIP) + 1
         lipof = (lipof-1)*NUMCIP + 1
         return
         end

C
C
C
C
C
C >>>>>> Store a string in the string table
C
         subroutine tbxxsstb(astrg,sindex)
C
C        store string astrg in the string table, returning the
C        index in sindex
C
         character *(*) astrg
         integer sindex
         include  'ciftbx.sys'
         character *(MAXBUF) temp
         integer mlen, ii, ibstb, icstb, ikstb, rlen
         integer iestb

         call tbxxrle(astrg,temp,mlen)
         icstb = mod(nstable,NUMCSTB)+1
         ibstb = (nstable+NUMCSTB)/NUMCSTB
         iestb = min(NUMCSTB,icstb+mlen-1)
         ikstb = iestb-icstb+1
         if (mlen+nstable .le. NUMCSTB*NUMSTB) then
           stable(ibstb)(icstb:iestb)=temp(1:ikstb)
           sindex = nstable+1
           nstable = nstable+mlen
           rlen = mlen - ikstb
           if (rlen .gt. 0) then
             do ii = ikstb+1,mlen,NUMCSTB
               ibstb = ibstb+1
               iestb = min(NUMCSTB,rlen)
               stable(ibstb)(1:iestb) = temp(ii:ii+iestb-1)
               rlen = rlen - iestb
             enddo
           endif
         else
           sindex = 0
           call tbxxwarn(
     *      ' More than NUMCSTB*NUMSTB stable characters needed')
         endif
         return
         end

C
C
C
C
C
C >>>>>> Fetch a string from the string table
C
         function tbxxfstb(astrg,sindex,fill)
C
C        fetch string astrg from the string table, starting at the
C        index in sindex, and returning the valid length.
C
C        fill is a logical variable, .true. to fill astr with blanks
C        the return value is the number of valid characters in astr
C        never less than 1, unless there is no valid string

         integer tbxxfstb
         character *(*)astrg
         integer sindex
         logical fill
         integer tbxxc2dig, tbxxrld
         integer rlen
         integer icstb, ibstb, iestb, ikstb, klen, ii

         include  'ciftbx.sys'
         character *(MAXBUF) temp

         tbxxfstb = 0

         if (sindex.le.0.or.nstable+3.gt.NUMCSTB*NUMSTB) return

         icstb = mod(sindex-1,NUMCSTB)+1
         ibstb = (sindex-1+NUMCSTB)/NUMCSTB
         iestb = min(NUMCSTB,icstb+3)
         ikstb = iestb-icstb+1
         temp(1:ikstb)=stable(ibstb)(icstb:iestb)

         rlen = 4-ikstb
         if (rlen .gt. 0) then
         temp(ikstb+1:4)=stable(ibstb+1)(1:rlen)
         endif

         klen = 0
         do ii = 1,4
           klen = klen*XXRADIX+tbxxc2dig(temp(ii:ii))
         enddo

         if (klen.gt.MAXBUF.or.klen.le.0) return
         if (sindex+klen-1.gt.NUMCSTB*NUMSTB) return

         if (klen.gt.4) then
           icstb = mod(sindex+3,NUMCSTB)+1
           ibstb = (sindex+3+NUMCSTB)/NUMCSTB
           iestb = min(NUMCSTB,icstb+klen-5)
           ikstb = iestb-icstb+1
           temp(5:ikstb+4) = stable(ibstb)(icstb:iestb)
           rlen = klen - ikstb - 4
           if (rlen .gt. 0) then
             do ii = ikstb+1,ikstb+rlen,NUMCSTB
               ibstb = ibstb+1
               iestb = min(NUMCSTB,rlen)
               temp(ii:ii+iestb-1) = stable(ibstb)(1:iestb)
               rlen = rlen - iestb
             enddo
           endif
         endif

         tbxxfstb = tbxxrld(astrg,temp(1:klen),fill)
         return
         end


C
C
C
C
C
C >>>>>> Open a CIF and copy its contents into a direct access file.
C
         function ocif_(fname)
C
         logical   ocif_
         integer   lastnb
         include  'ciftbx.sys'
         logical   test
         character fname*(*)
         integer   lfname
         integer   case,i,kp,lp,mp
         integer   klen, mlen, lip, lppag, lipof, kip, ip, mip, mis
C
         save_=.false.
         glob_=.false.
         depth_=0
         jchar=MAXBUF
         lastch=0
         if(line_.gt.MAXBUF) call tbxxerr(' Input line_ value > MAXBUF')
         if(nrecd.ne.0 .and. (.not.append_)) then
           close(dirdev)
           nrecd=0
           lrecd=0
         endif
C
C        clear the memory resident page buffer
C
         do i = 1,NUMPAGE
         mppoint(i)=0
         enddo
C
         case=ichar('a')-ichar('A')
         tab=char(05)
         if(case.lt.0) goto 100
         tab=char(09)
         bloc_=' '
C
C....... Make sure the CIF is available to open
C
100      file_(1:longf_)=' '
         lfname = len(fname)
         file_(1:lfname) = fname
         do 120 i=1,lfname
         if(file_(i:i).eq.' ' .or. file_(i:i).eq.char(0) ) goto 140
120      continue
140      longf_=i-1
         if (longf_.gt.0) then
           inquire(file=file_(1:longf_),exist=test)
           ocif_=test
           if(.not.ocif_)      goto 200
         else
           file_(1:1) = ' '
           longf_ = 1
           ocif_ = .true.
         endif
C
C....... Open up the CIF and a direct access formatted file as scratch
C
         if (file_(1:1).ne.' ')
     *   open(unit=cifdev,file=file_(1:longf_),status='OLD',
     *                    access='SEQUENTIAL',
     *                    form='FORMATTED')
         if(nrecd.eq.0)  then
           open(unit=dirdev,status='SCRATCH',access='DIRECT',
     *                    form='FORMATTED',recl=NUMCPP)
           mipim = -1
           miprim = 1
           mipcp = 0
           ipim = -1
           iprim = 2
           do i = 1,NUMPAGE
             mppoint(i) = 0
           enddo
           do i = 1,NUMMIP
             mippoint(i) = 0
           enddo
           mippoint(1)=2
           do i = 1,NUMSIP
             ippoint(i) = 0
           enddo
           nfblock = 3
           nfword = 1
         endif
         if (mppoint(1).lt.0) then
            write(dirdev,'(a)',rec=-mppoint(1)) pagebuf(1)
            mppoint(1) = 0
         endif
         if(append_ .and. nrecd.ne.0) then
           kp = 1
           lp = nfblock
           nfblock = nfblock+1
           mppoint(kp) = lp
           mp = 1
         else
           do kp = 1,NUMPAGE
             mppoint(kp)=0
           enddo
           kp = 1
           lp = 3
           nfblock = 4
           mp = 1
         endif
C
C....... Copy the CIF to the direct access file
C
160      read(cifdev,'(a)',end=180) buffer
         nrecd=nrecd+1
         irecd=nrecd
         klen = lastnb(buffer(1:MAXBUF))
         if (klen.gt.line_)
     *     call tbxxwarn(' Input line length exceeds line_')
         call tbxxrle(buffer(1:klen),scrbuf,mlen)
         if (mp+mlen-1 .gt. NUMCPP) then
           if (mp.lt.NUMCPP) pagebuf(kp)(mp:NUMCPP) = ' '
C          write(dirdev,'(a)',rec=lp) pagebuf(kp)
           mppoint(kp)=-lp
           if (nfword.gt.1) then
             nfblock = nfblock+1
             nfword = 1
           endif
           lp = nfblock
           nfblock=nfblock+1
           kp = kp+1
           if(kp.gt.NUMPAGE) kp=1
           if (mppoint(kp).lt.0) then
             write(dirdev,'(a)',rec=-mppoint(kp)) pagebuf(kp)
           endif
           mppoint(kp)=0
           mp=1
         endif
         pagebuf(kp)(mp:mp+mlen-1)=scrbuf(1:mlen)
         mppoint(kp) = -lp
         mlen = ((mlen+NUMCIP-1)/NUMCIP)
         mlen = mlen*NUMCIP
         call tbxxflin(nrecd,lip,lppag,lipof,kip,ip,mip,mis)
         ippoint(ip) = (mp-1)/NUMCIP+(lp-1)*(NUMCPP/NUMCIP)+1
         ipim = -iabs(ipim)
         mp = mp+mlen
         goto 160
C
180      if (mp.lt.NUMCPP) pagebuf(kp)(mp:NUMCPP) = ' '
         if (mp.gt.1) then
C          write(dirdev,'(a)',rec=lp) pagebuf(kp)
           mppoint(kp)=-lp
         endif
         lrecd=max(0,recbeg_-1)
         jrecd=max(0,recbeg_-1)
         jrect=-1
         irecd=max(0,recbeg_-1)
         recn_=irecd
         recend_=nrecd
         if (file_(1:1).ne.' ') close(cifdev)
200      return
         end
C
C
C
C
C
C >>>>>> Close off direct access file of the current CIF
C         and reset all data name tables and pointers
C
         subroutine purge_
C
         include   'ciftbx.sys'
C
         integer i
         if(nrecd.ne.0) close(dirdev)
         do i = 1,NUMPAGE
           mppoint(i)=0
         enddo
         do i = 1,MAXBOOK
           ibkmrk(1,i)=-1
           ibkmrk(2,i)=-1
           ibkmrk(3,i)=-1
           ibkmrk(4,i)=-1
           ibkmrk(5,i)=-1
           ibkmrk(6,i)=-1
         enddo
         recn_=0
         save_=.false.
         glob_=.false.
         jchar=MAXBUF
         depth_=0
         lastch=0
         nrecd=0
         lrecd=0
         irecd=0
         nname=0
         nhash=0
         iname=0
         loopct=0
         loopnl=0
         loop_=.false.
         text_=.false.
         textfl='no '
         append_=.false.
         recbeg_=0
         recend_=0
         nivt = 0
         nstable = 0
         return
         end
C
C
C
C
C
C >>>>>> Store the data names and pointers for the requested data block
C
         function data_(name)
C
         logical   data_
         logical   wasave
         logical   tbxxoldd
         integer   lastnb
         include  'ciftbx.sys'
         character name*(*),temp*(NUMCHAR),ltype*4
         character ctemp*(NUMCHAR)
         character xdname*(NUMCHAR)
         character ydname*(NUMCHAR)
         character isbuf*(MAXBUF),lsbuf*(MAXBUF)
         logical   ixcat(NUMDICT)
         integer   ndata,idata,nitem,npakt,i,ii,j,k,kchar,krecd
         integer   jj,icc,idd
         integer   fcatnum,lctemp,isrecd,isjchr,islast
         integer   lsrecd,lsjchr,lslast
         integer   pnname,itpos,ipp,ipj
         integer   ltemp
CDBG     if(dictfl.eq.'no ')
CDBG *     print *,' ***>>>> Entering data_ ',name
C
         jchar=MAXBUF
         depth_=0
         nname=0
         ndata=0
         nhash=0
         nitem=0
         idata=0
         iname=0
         loopct=0
         loopnl=0
         ltype=' '
         posnam_=0
         posval_=0
         posdec_=0
         posend_=0
         kchar = 0
         krecd = 0
         fcatnum = 0
         data_=.false.
         wasave=.false.
         loop_=.false.
         text_=.false.
         textfl='no '
         glob_=.false.
         do ii = 1,MAXBOOK
         ibkmrk(1,ii)=-1
         enddo
         irecd=lrecd
         lrecd=min(nrecd,recend_)
         if(name(1:1).ne.' ') irecd=max(0,recbeg_-1)
         call hash_init(dname,dchain,NUMBLOCK,nname,dhash,
     *     NUMHASH)
         call hash_init(cname,cchain,NUMBLOCK,ncname,chash,
     *     NUMHASH)
         isrecd=irecd
         isjchr=jchar
         islast=lastch
         lsrecd=isrecd
         lsjchr=isjchr
         lslast=islast
         isbuf=' '
         if(lastch.gt.0)isbuf(1:lastch)=buffer(1:lastch)
         lsbuf=' '
         if(lastch.gt.0)lsbuf(1:lastch)=isbuf(1:lastch)
         call tbxxnlc(xdname,name)
C
C....... Find the requested data block in the file
C
100      lsjchr=isjchr
         call getstr
         isjchr=jchar
         if(irecd.ne.isrecd) then
           lsrecd=isrecd
           lslast=islast
           lsbuf=' '
           if(islast.gt.0)lsbuf(1:islast)=isbuf(1:islast)
           isrecd=irecd
           islast=lastch
           isbuf=' '
           if(lastch.gt.0)isbuf(1:lastch)=buffer(1:lastch)
         endif
         if(type_.eq.'fini')           goto 500
         if(text_.or.depth_.gt.0)      goto 110
         goto 120
110      call getstr
         if (type_.eq.'fini')
     *      call tbxxerr(' Unexpected termination of file')
         if (text_.or.depth_.gt.0)     goto 100
         goto 100
120      continue
         if(type_.eq.'save') then
           if(long_.lt.6) then
             if(.not.save_)
     *         call tbxxerr(
     *           ' Save frame terminator found out of context ')
             wasave=.true.
             save_=.false.
             goto 100
           else
             if(save_)
     *         call tbxxerr(' Prior save frame not terminated ')
             save_=.true.
             if(name.eq.' ')          goto 150
             call tbxxnlc(ydname,strg_(6:long_))
             if(ydname.ne.xdname) goto 100
             goto 150
           endif
         endif
         if(type_.eq.'glob') then
           if(name.ne.' ')            goto 100
           glob_=.true.
           goto 150
         endif
         if(type_.eq.'name'.or.type_.eq.'loop') then
           if(name.ne.' ')            goto 100
           if(.not.wasave)
     *       call tbxxwarn(' Data block header missing ')
           isrecd=lsrecd
           islast=lslast
           isjchr=lsjchr
           isbuf=' '
           if(islast.gt.0)isbuf(1:islast)=lsbuf(1:islast)
           data_=.true.
           bloc_=' '
           itpos=jchar-long_
           if(tabx_) then
           itpos=0
           do ipp=1,jchar-long_
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           endif
           posnam_=itpos
           goto 204
         endif
         if(type_.ne.'data')          goto 100
         if(name.eq.' ')              goto 150
         call tbxxnlc(ydname,strg_(6:long_))
         if(ydname.ne.xdname)   goto 100
150      data_=.true.
         bloc_=strg_(6:long_)
C
CDBG     if(dictfl.eq.'no ')
CDBG *     print *, 'bloc_: '//bloc_
         itpos=jchar-long_
         if(tabx_) then
         itpos=0
         do ipp=1,jchar-long_
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         endif
         posnam_=itpos
C
C....... Get the next token and identify
C        ltype is the previous type
C
200      call getstr
CDBG     WRITE(6,*) ltype,type_,loop_,nitem,ndata,idata,iname,nname
C
         if(ltype.ne.'name')                goto 201
         if(type_.eq.'numb')                goto 203
         if(type_.eq.'char')                goto 203
         if(type_.eq.'text')                goto 203
         if(type_.eq.'null')                goto 203
         if(type_.eq.'tupl'
     *     .or.type_.eq.'tabl'
     *     .or.type_.eq.'arra')             goto 203
         if(type_.eq.'name'.and.loop_)      goto 204
CDBG     WRITE(6,*) ltype,type_,loop_,nitem,ndata,idata,iname,nname
         call tbxxerr(
     *     ' Illegal tag/value construction: tag followed by '
     *     //type_)    
C
C        The prior type was not a name (not a tag)
201      if(ltype.ne.'valu')                goto 204
C
C        The prior type was a data value
C
         if(type_.eq.'numb')                goto 202
         if(type_.eq.'char')                goto 202
         if(type_.eq.'text')                goto 202
         if(type_.eq.'null')                goto 202
         if(type_.eq.'tupl'
     *     .or.type_.eq.'tabl'
     *     .or.type_.eq.'arra')             goto 202
         goto 204
C
C        If we have a vaue followed by a value, we need to be
C        in a loop (item > 0)
C
202      if(nitem.gt.0)                     goto 205
         call tbxxerr(
     *     ' Illegal tag/value construction: value followed by '
     *     //type_)
C
C        The prior item was a tag and this is a value
C
203      ltype='valu'
CDBG     if(dictfl.eq.'no ')
CDBG *     print *, ' ***>>>>> data_ value ',strg_(1:long_)
         goto 205
C
C        Cases that get us here
C          The prior item was a tag and this is a tag in a loop
C          The prior item was neither a tag nor a value 
204      ltype=type_
C
C        We are in a loop and have a value after a value
C        or a name after a value or come from above cases
C
205      if(type_.eq.'name')           goto 206
         if(type_.eq.'loop')           goto 210
         if(type_.eq.'data')           goto 210
         if(type_.eq.'save')           goto 210
         if(type_.eq.'glob')           goto 210
         if(type_.ne.'fini')           goto 220
206      if(loop_)                     goto 270
210      if(nitem.eq.0)                goto 215
C
C....... End of loop detected; save pointers
C        loopni(loopct) -- number of items in a row
C        loopnp(loopct) -- number of rows
C
         npakt=idata/nitem
         if(npakt*nitem.ne.idata) call tbxxerr(' Item miscount in loop')
         loopni(loopct)=nitem
         loopnp(loopct)=npakt
         nitem=0
         idata=0
215      if(type_.eq.'name')           goto 270
         if(type_.eq.'data')           goto 300
         if(type_.eq.'save')           goto 300
         if(type_.eq.'glob')           goto 300
         if(type_.eq.'fini')           goto 300
C
C....... Loop_ line detected; incr loop block counter
C        record the character position in loopos(loopct)
C        record the line number in        loorec(loopct)
C        record the detabbed char pos in  loopox(loopct)
C
         loop_=.true.
CDBG     print *,' in data_ loop_ set, type_', type_
         loopct=loopct+1
         if(loopct.gt.NUMLOOP) call tbxxerr(
     *     ' Number of loop_s > NUMLOOP')
         loorec(loopct)=irecd
         loopos(loopct)=jchar-long_
         if(quote_.ne.' ') then
           if (quote_.eq.';') then
             loopos(loopct) = 1
           else
             if (quote_.eq.''''''''.or.quote_.eq.'"""') then
               loopos(loopct)=jchar-long_-3
             else
               loopos(loopct)=jchar-long_-1
             end if
           end if
         end if
         itpos=0
         do ipp=1,loopos(loopct)
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         loopox(loopct)=itpos
         goto 200
C
C....... This is the data item; store char position and length
C
220      if(loop_ .and. nitem.eq.0)
     *   call tbxxerr(' Illegal tag/value construction')
         loop_=.false.
C
         i=nname
         if(nitem.gt.0) i=i-nitem+mod(idata,nitem)+1
         if(i.lt.1) call tbxxerr(' Illegal tag/value construction')
         if(dtype(i).ne.'test')       goto 223
         if(dictfl.eq.'yes')          goto 223
         if(tcheck.eq.'no ')          goto 223
C>>>>    if(long_.eq.1.and.strg_(1:1).eq.'?') goto 223
C>>>>    if(long_.eq.1.and.strg_(1:1).eq.'.') goto 223
         if(type_.eq.'null')          goto 223
         if(type_.eq.'numb')          goto 223
         call tbxxwarn( ' Numb type violated  '//dname(i))
223      if(nitem.le.0)               goto 224
         idata=idata+1
         if(dtype(i).eq.'null') dtype(i)=type_
         if(dtype(i).eq.'numb' .and.
     *     (type_.eq.'char'.or.type_.eq.'text')) dtype(i)='char'
224      if(nname.eq.ndata)           goto 230
         ndata=ndata+1
         if(iloop(ndata).gt.1)        goto 225
         krecd=irecd
         kchar=jchar-long_-1
         if(quote_.ne.' ') then
           kchar=kchar-1
           if (quote_(2:3).ne.'  ') kchar=kchar-2
         end if
225      continue
         if(dtype(ndata).eq.'    ') dtype(ndata)=type_
         drecd(ndata)=krecd
         dchar(ndata)=kchar
         if (depth_.gt.0) then
CDBG     print *,' Setting bracket start at ',
CDBG *     'char: ', posbrkstk(1)-1, 'rec: ',srecd
         dchar(ndata) = posbrkstk(1)-1
         drecd(ndata) = srecd

         end if
         if(nloop(ndata).gt.0)        goto 230
         nloop(ndata)=0
         iloop(ndata)=long_
         if (depth_.gt.0) iloop(ndata) = 1
C
C....... Skip text lines if present
C
230      if(type_.ne.'text')           goto 250
CDBG     print *,' text field detected at 230 '
         if(nloop(ndata).eq.0.and.depth_.eq.0) dchar(ndata)=0
         if(nloop(ndata).eq.0.and.depth_.eq.0) iloop(ndata)=long_
240      call getstr
         if(type_.eq.'fini') call tbxxerr(' Unexpected end of data')
         if (type_.ne.'text'.or..not.text_) then
           if (depth_.eq.0)            goto 200
           goto 260
         endif
         goto 240
C
C....... Skip bracketed construct if present
C
250      if(depth_.eq.0)           goto 200
         
260      call getstr
         if(depth_.eq.0) goto 200
         if(type_.eq.'fini') call tbxxerr(' Unexpected end of data')
         if(type_.eq.'text') goto 240
         goto 260
C
C....... This is a data name; store name and loop parameters
C
270      call tbxxclc(temp,ltemp,strg_(1:long_),long_)
         k=0
         if(dictfl.ne.'yes' .and. ndict.gt.0) then
           tbxxrslt = tbxxoldd(temp(1:ltemp),k)
           if(k.ne.0) then
             if(alias_ .and. aroot(k).ne.0) then
               temp=dicnam(aroot(k))
               ltemp = lastnb(temp)
             endif
           endif
         endif
         pnname=nname
         call hash_store(temp(1:ltemp),
     *   dname,dchain,NUMBLOCK,nname,dhash,
     *     NUMHASH,j)
CDBG     if(dictfl.eq.'no ')
CDBG *     print *,' ***>>>>> data_ name: ',temp(1:ltemp)
         if(j.eq.pnname+1) then
           dtag(j)=strg_(1:long_)
           if(k.ne.0) dtag(j)=dictag(k)
           trecd(j)=irecd
           tchar(j)=jchar-long_
           if(quote_.ne.' '.and.quote_.ne.';') 
     *       tchar(j)=jchar-long_-1
           itpos=0
           do ipp=1,tchar(j)
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           xchar(j)=itpos
         endif
         if(j.eq.0)
     *     call tbxxerr(' Number of data names > NUMBLOCK')
         if(k.ne.0) then
           ltemp = lastnb(dicnam(k))
           temp(1:ltemp) = dicnam(k)(1:ltemp)
         endif
         if(j.ne.pnname+1) then
           call tbxxwarn(' Duplicate data item '//
     *     temp(1:ltemp))
           goto 200
         endif
         dtype(nname)=' '
         dxtyp(nname)=' '
         cindex(nname)=0
         ddict(nname)=0
         ctemp(1:6)='(none)'
         lctemp=6
C
         if(dictfl.eq.'yes' .or. vcheck.eq.'no ') goto 290
         j=k
         if(j.ne.0) then
           ddict(nname)=j
           cindex(nname)=dcindex(j)
           dxtyp(nname)=dicxtyp(j)
           dtype(nname)=dictyp(j)
           if(vcheck.eq.'no ')          goto 280
           if(dictyp(j).eq.'numb') then
             dtype(nname)='test'
           endif
           if(cindex(nname).ne.0) then
             lctemp=lastnb(dcname(cindex(nname)))
             ctemp(1:lctemp)=dcname(cindex(nname))(1:lctemp)
             goto 290
           endif
           goto  280
         endif
         call tbxxwarn(' Data name '//
     *               temp(1:ltemp)
     *               //' not in dictionary!')
280      call tbxxcat(temp(1:ltemp),ctemp,lctemp)
         if (ctemp(1:lctemp).eq.' '.or.
     *     ('_'//ctemp(1:lctemp).eq.temp(1:ltemp))) then
           ctemp = '(none)'
           lctemp= 6
           if (ndcname.ne.0.and.vcheck.eq.'yes')
     *       call tbxxwarn(' No category defined for '
     *       //temp(1:ltemp))
         else
           call hash_find(ctemp(1:lctemp),
     *       dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,j)
           if(j.ne.0) then
             cindex(nname) = j
           else
             ipj=ncname
             call hash_store(ctemp(1:lctemp),
     *         cname,cchain,NUMBLOCK,ncname,chash,NUMHASH,j)
             if (j.eq.0)
     *         call tbxxerr(' Number of categories > NUMBLOCK ')
             cindex(nname) = -j
             if (ndcname.gt.0.and.j.eq.ipj+1.and.vcheck.eq.'yes'
     *         .and.catchk.eq.'yes')
     *         call tbxxwarn(' Category '//
     *         ctemp(1:lctemp)//' first implicitly defined in cif ')
           endif
         endif
C
290      lloop(nname)=0
         nloop(nname)=0
         iloop(nname)=0
         if (nitem.eq.0) fcatnum=cindex(nname)
         if(.not.loop_)               goto 200
         nitem=nitem+1
         if(nitem.gt.NUMITEM)
     *     call tbxxerr(' Items per loop packet > NUMITEM')
         nloop(nname)=loopct
         iloop(nname)=nitem
         if (fcatnum.ne.cindex(nname)) then
           temp = '(none)'
           if (fcatnum.gt.0) temp=dcname(fcatnum)
           if (fcatnum.lt.0) temp=cname(-fcatnum)
           ltemp = lastnb(temp)
           if (ctemp(1:lctemp).ne.temp(1:ltemp)
     *     .and.catchk.eq.'yes')
     *     call tbxxwarn (' Heterogeneous categories in loop '//
     *     ctemp(1:lctemp)//' vs '//
     *     temp(1:ltemp))
           fcatnum=cindex(nname)
         endif
         goto 200
300      continue
C
C....... Are names checked against dictionary?
C
         if(dictfl.eq.'yes')          goto 500
         if(vcheck.eq.'no '.or.ndict.eq.0) goto 500
         do i=1,nname
           if(dtype(i).eq.'test') dtype(i)='numb'
         enddo

C
C        prepare for category and parent checks
C
         if ((catchk.eq.'yes'.or.parchk.eq.'yes')
     *      .and. ndict.gt.0) then
         do i = 1,ndict
           ixcat(i) = .false.
         enddo
C
C        make a pass marking all used tags and their aliases
C
         do i = 1,nname
           icc=cindex(i)
           idd=ddict(i)
           if(icc.ne.0.and.idd.ne.0) then
             icc = aroot(idd)
310          ixcat(icc) = .true.
             icc = alias(icc)
             if (icc.ne.0) goto 310
           endif
         enddo
         endif
C
C        check for category keys
C
C
C
C        now make a pass making certain the keys are
C        used
C
         if(catchk.eq.'yes' .and. ndict.gt.0) then
         do i = 1,nname
           idd=cindex(i)
           if (idd.gt.0) then
             icc=ccatkey(idd)
             if(icc.ne.0) then
             if(aroot(icc).ne.0) icc=aroot(icc)
320          if(icc.ne.0) then
               if(.not.ixcat(icc)) then
                 jj = irecd
                 irecd = drecd(i)
                 if (catkey(icc)) then
                   call tbxxwarn(' Category key '//
     *               dictag(icc)(1:lastnb(dictag(icc)))//
     *               ' not given for '//
     *               dcname(idd)(1:lastnb(dcname(idd))))
                 else
                   call tbxxwarn(' Mandatory item '//
     *               dictag(icc)(1:lastnb(dictag(icc)))//
     *               ' not given for '//
     *               dcname(idd)(1:lastnb(dcname(idd))))
                 endif
                 ixcat(icc) = .true.
                 irecd = jj
               endif
               icc = keychain(icc)
               if(icc.ne.0) go to 320
             endif
             endif
           endif
         enddo
         endif
C
C        check for parents of tags that are used
C
         if(parchk.eq.'yes' .and. ndict.gt.0) then
         do i = 1,nname
           if (ddict(i).ne.0) then
             if (dpindex(ddict(i)).ne.ddict(i)) then
               if (.not.ixcat(dpindex(ddict(i)))) then
                 call tbxxwarn(' Parent '//
     *             dicnam(dpindex(ddict(i)))
     *             (1:lastnb(dicnam(dpindex(ddict(i)))))//
     *             ' of '//
     *             dname(i)(1:lastnb(dname(i))) //
     *             ' not given')
               endif
             endif
           endif
         enddo
         endif

C
C....... End of data block; tidy up loop storage
C
500      lrecd=irecd-1
         if(type_.eq.'save'.and.long_.lt.6) then
           itpos=jchar-long_
           if(tabx_) then
           itpos=0
           do ipp=1,jchar-long_
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           endif
           posval_=itpos
         endif
         irecd=isrecd
         jchar=isjchr
         lastch=islast
         recn_=irecd
         buffer(1:1)=' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=isbuf(1:lastch)
         jrecd=irecd
         loop_=.false.
         loopct=0
         if(ndata.ne.nname) call tbxxerr(' Syntax construction error')
C
Cdbg     WRITE(6,'(a)')
Cdbg *   ' data name                       type recd char loop leng'
Cdbg     WRITE(6,'(a,1x,a,4i5)') (dname(i),dtype(i),drecd(i),dchar(i),
Cdbg *              nloop(i),iloop(i),i=1,nname)
Cdbg     WRITE(6,'(3i5)') (i,loopni(i),loopnp(i),i=1,loopct)
C
         return
         end
C
C
C
C
C
C >>>>>> Check dictionary for data name validation
C
         function dtype_(name,type)
C
         logical    dtype_, tbxxoldd
         include   'ciftbx.sys'
         integer    nln, ii
         character  name*(*),temp*(NUMCHAR),
     *              type*4
C
         character*4 map_type(19),map_to(19),mapped
         data map_type
     *   /'floa','int ','yyyy','symo','ucha','ucod','name','idna',
     *    'any ','code','line','ulin','atco','fax ','phon','emai',
     *    'real','inte','coun'/
         data map_to
     *   /'numb','numb','char','char','char','char','char','char',
     *    'char','char','char','char','char','char','char','char',
     *    'numb','numb','numb'/

C
         type = ' '
         dtype_ = .false.
         nln = min(len(name),len(temp))
         call tbxxnlc(temp(1:nln),name)
         if (ndict.eq.0) go to 200
         tbxxrslt = tbxxoldd(temp(1:nln),xdchk)
         if(xdchk.eq.0) go to 200
         mapped = dictyp(xdchk)(1:4)
         do ii = 1,19
           if (dictyp(xdchk)(1:4).eq.map_type(ii)) mapped = map_to(ii)
         enddo
         if (mapped.ne.'char'.and.mapped.ne.'numb'
     *     .and.mapped.ne.'null'.and.mapped.ne.'text') then
           call tbxxwarn(' Item type '
     *        //dictyp(xdchk)(1:max(1,lastnb(dictyp(xdchk))))//
     *       ' for item '//
     *       name(1:max(1,lastnb(name)))//' not recognized ')
           mapped = 'char'
         endif
         type = mapped
         dtype_ = .true.
200      continue
         return
         end
C
C
C
C
C
C
C >>>>>> Get the attributes of data item associated with data name
C
         function test_(temp)
C
         logical   test_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character*4 otype
         integer lname
         character  otestf*3
C
         otestf=testfl
         otype = type_
         testfl='yes'
         test_ = .false.
         call tbxxclc(name,lname,temp,len(temp))
CDBG     print *,' Entering test_ ',name(1:lname)
         if (depth_.eq.0) go to 100
         if (name(1:1).ne.' '.and.name(1:1).ne.char(0).and.
     *     name(1:lname).ne.nametb(1:lnametb))     goto 120
         call getstr
         test_=.true.
         if (type_.eq.'null') test_=.false.
         if (otype.eq.'text' .and. (.not. text_) .and.long_.eq.0) then
           quote_=' '
           textfl = 'no'
           type_ = 'null'
           test_ = .false.
           goto 200
         end if
         posval_ = jchar-long_
         posend_ = jchar-1
         if (long_.gt.0) then
            if (type_.eq.'numb') then
               call ctonum
               if(posdec_.gt.0) posdec_=posval_+posdec_-1
               jchar = posend_
            else 
              if (quote_.eq.' ') then
                 if (long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
                 if (long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
              end if
            end if
         end if
         goto 200

         
100      test_=.true.
         if(otestf.eq.'no ' .or. type_.eq.' ')  goto 120
         if(name(1:lname).eq.nametb(1:lnametb))   goto 200
120      call tbxxgitm(name(1:lname))
200      list_ =loopnl
         if(type_.eq.'null') test_=.false.
         if(type_.ne.'null'.and.type_.ne.'char'.and.
     *     type_.ne.'text'.and.type_.ne.'numb') type_='char'
CDBG     print *,' leaving test_ ', type_, depth_, strg_(1:long_)
         return
         end

C
C
C
C
C
C >>>>>> Set or Reference a bookmark
C
         function bkmrk_(mark)
C
         logical   bkmrk_
         include   'ciftbx.sys'
C
         integer   mark,ii,nitem
         character*4 flag
         bkmrk_=.true.
         if(mark.eq.0) then
           do ii=1,MAXBOOK
             if(ibkmrk(1,ii).lt.0)      goto 100
           enddo
           bkmrk_=.false.
           call tbxxwarn(' More than MAXBOOK bookmarks requested')
           return
100        mark=ii
           ibkmrk(1,ii)=iname
           ibkmrk(2,ii)=irecd
           ibkmrk(3,ii)=jchar
           if(iname.gt.0) then
             ibkmrk(2,ii) = trecd(iname)
             ibkmrk(3,ii) = tchar(iname)
           endif
           ibkmrk(4,ii)=0
           if(iname.gt.0) then
             if(nloop(iname).ne.0.and.
     *         loopnl.eq.nloop(iname).and.loopct.ne.0) then
               nitem=loopni(nloop(iname))
               ibkmrk(2,ii)=looprd(1)
               ibkmrk(3,ii)=max(0,loopch(1)-1)
               ibkmrk(4,ii)=loopct
             endif
           endif
           ibkmrk(5,ii) = depth_
           ibkmrk(6,ii) = index_
         else
           if(ibkmrk(1,mark).lt.0) then
             bkmrk_=.false.
             return
           endif
           iname=ibkmrk(1,mark)
           irecd=ibkmrk(2,mark)
           loopct=ibkmrk(4,mark)
           loop_=.false.
           text_=.false.
           textfl = 'no '
           loopnl=-1
           testfl='no '
           if(iname.gt.0) then
            if(nloop(iname).ne.0.and.loopct.ne.0) then
               nitem=loopni(nloop(iname))
               looprd(nitem+1)=ibkmrk(2,mark)
               loopch(nitem+1)=ibkmrk(3,mark)
               do ii = 1,nitem
                 lloop(ii+iname-iloop(iname))=loopct-1
               enddo
               loopct=loopct-1
               if(lloop(iname).gt.0) then
                 loop_=.true.
                 loopnl=nloop(iname)
               endif
             endif
           endif
           jchar=MAXBUF
           if(irecd.gt.0) then
             irecd=irecd-1
             call getlin(flag)
             jchar=ibkmrk(3,mark)
           endif
           depth_=0
           index_=0
           if (ibkmrk(5,mark).gt.0) then
200           call getstr
              if (depth_ .lt. 1) then
                call tbxxwarn(
     *          ' Bookmark for list, array, tuple or table corrupted')
                go to 210
              end if
              if(ibkmrk(5,mark).ne.depth_
     *          .or. ibkmrk(6,mark).ne.index_ ) go to 200
           endif
210        ibkmrk(1,mark)=-1
           mark=0
         endif
         return
         end
C
C
C
C
C
C
C >>>>>> Find the location of the requested item in the CIF
C        The argument "name" may be a data item name, blank
C        for the next such item.  The argument "type" may be
C        blank for unrestricted acceptance of any non-comment
C        string (use cmnt_ to see comments), including loop headers,
C        "name" to accept only the name itself and "valu"
C        to accept only the value, or "head" to position to the
C        head of the CIF.  Except when the "head" is requested,
C        the position is left after the data item provided.
C
         function find_(name,type,strg)
C
         logical   find_
         include   'ciftbx.sys'
         character  name*(*),type*(*),strg*(*),flag*4
         character  jjbuf*(MAXBUF)
         integer    jjchar,jjrecd,jjlast,jjlrec,jjjrec,jjdepth,jindex
C
CDBG     print *,' Entering find ', name, type
         find_  = .false.
         strg   = ' '
         long_  = 0
         jjchar = jchar
         jjrecd = lrecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjdepth = depth_
         jindex = index_
         jjbuf  = ' '
         if(lastch.gt.0) jjbuf(1:lastch)=buffer(1:lastch)
         if(type.eq.'head') then
           lrecd = min(nrecd,recend_)
           irecd = max(0,recbeg_-1)
           jchar=MAXBUF+1
           depth_=0
           call getlin(flag)
           if(flag.eq.'fini')       goto 300
           find_=.true.
           lrecd=max(0,recbeg_-1)
           return
         endif
         if(name.ne.' ') then
           testfl='no '
           call tbxxgitm(name)
           if(iname.eq.0) goto 300
           if(type.eq.'valu') then
             list_=loopnl
             strg=strg_(1:long_)
             find_=.true.
             return
           endif
           if(type.eq.'name'.or.loopnl.eq.0) then
             irecd=trecd(iname)-1
             call getlin(flag)
             jchar=tchar(iname)
             depth_=0
             posnam_=jchar+1
             call getstr
             strg=strg_(1:long_)
             recn_=irecd
             find_=.true.
             return
           endif
           if(type.eq.' ') then
             irecd=loorec(loopnl)-1
             call getlin(flag)
             jchar=loopos(loopnl)
             depth_=0
             call getstr
             posval_=loopos(loopnl)
             if(tabx_) posval_=loopox(loopnl)
             strg=strg_(1:long_)
             recn_=irecd
             find_=.true.
             return
           endif
           call tbxxerr(' Call to find_ with invalid arguments')
         endif
         if(name.eq.' ') then
           go to 200
190        if (text_.or.depth_.gt.0) then
              call getstr
              if (type_.eq.'fini')  goto 300
              if (type_.ne.'null')  goto 190  
           end if     
200        call getstr
           if(type_.eq.'fini')      goto 300
           if(type.ne.' '.and.
     *      (type_.eq.'data'.or.type_.eq.'save'.or.
     *      type_.eq.'glob'))   goto 300
           if(type.eq.'name'.and.type_.ne.'name')  goto 190
           if(type.eq.'valu'.and.
     *       type_.ne.'numb'.and.type_.ne.'text'
     *      .and.type_.ne.'char'.and.type_.ne.'null') goto 190
           find_=.true.
           strg=strg_(1:long_)
           if(type_.eq.'name') then
             posnam_=jchar-long_
           else
             posval_=jchar-long_
             if(quote_.ne.' '.and.quote_.ne.';') 
     *         posval_=posval_-1
             if(quote_.eq.'''''''' .or.quote_.eq.'"""')
     *         posval_=posval_-2
           endif
           recn_=irecd
           return
         endif

C
C        Search failed, restore pointers
C
300      irecd  = jjrecd
         lastch = jjlast
         lrecd  = jjlrec
         jchar  = jjchar
         depth_ = jjdepth
         index_ = jindex

         buffer(1:1) = ' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=jjbuf(1:lastch)
         jrecd  = jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_  = irecd
C
         return
         end
C
C
C
C
C
C
C >>>>>> Get the next data name in the data block
C
         function name_(temp)
C
         logical    name_
         include   'ciftbx.sys'
         character  temp*(*)
C
         name_=.false.
         temp=' '
         iname=iname+1
         if(iname.gt.nname)  goto 100
         name_=.true.
         temp=dtag(iname)
         if(ddict(iname).ne.0) temp=dictag(ddict(iname))
100      return
         end
C
C
C
C
C
C
C >>>>>> Extract a number data item and its standard deviation
C        This version return single precision numbers
C
         function numb_(temp,numb,sdev)
C
         logical    numb_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         integer    lname
         real       numb,sdev

CDBG     print *,'***>>> Entering numb_ for ', temp

C
         call tbxxclc(name,lname,temp,len(temp))
         if(testfl.eq.'yes')    goto 100
         if(depth_.eq.0)        goto 120
         if(name(1:1).ne.' '.and.name(1:1).ne.char(0).and.
     *     name(1:lname).ne.nametb(1:lnametb))     goto 120
     
         numb_ = .false.
         call getstr
         if (type_.ne.'numb') go to 200 
         if (type_.eq.'numb') then
            call ctonum
            if(posdec_.gt.0) posdec_=posval_+posdec_-1
            numb_ = .true.
            if (depth_.gt.0) jchar=jchar-1
         end if
CDBG     print *,'***>>> In numb_ strg_ ', strg_(1:long_)
         go to 200

100      if(name(1:lname).eq.nametb(1:lnametb))   goto 150
C
120      call tbxxgitm(name(1:lname))
C
150      continue
CDBG     print *,'***>>> In numb_ strg_ ', strg_(1:long_)
         numb_=.false.
         if(type_.ne.'numb') goto 200
         numb_=.true.
         numb =sngl(numbtb)
         if(sdevtb.ge.0.0) sdev=sngl(sdevtb)
C
200      testfl='no '
         return
         end
C
C
C
C
C
C
C >>>>>> Extract a number data item and its standard deviation
C        This version returns double precision numbers
C
         function numd_(temp,numb,sdev)
C
         logical    numd_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         integer    lname
         double precision numb,sdev
CDBG     print *,'***>>> Entering numb_ for ', temp
C
         call tbxxclc(name,lname,temp,len(temp))
         if(testfl.eq.'yes')    goto 100
         if(depth_.eq.0)        goto 120
         if(name(1:1).ne.' '.and.name(1:1).ne.char(0).and.
     *     name(1:lname).ne.nametb(1:lnametb))     goto 120
     
         numd_ = .false.
         call getstr
         if (type_.ne.'numb') go to 200 
         if (type_.eq.'numb') then
            call ctonum
            if(posdec_.gt.0) posdec_=posval_+posdec_-1
            numd_ = .true.
            if (depth_.gt.0) jchar=jchar-1
         end if
CDBG     print *,'***>>> In numd_ strg_ ', strg_(1:long_)
         go to 200

100      if(name(1:lname).eq.nametb(1:lnametb))   goto 150
C
120      call tbxxgitm(name(1:lname))
C
150      numd_=.false.
CDBG     print *,'***>>> In numd_ strg_ ', strg_(1:long_)
         if(type_.ne.'numb') goto 200
         numd_=.true.
         numb =numbtb
         if(sdevtb.ge.0.0) sdev=sdevtb
C
200      testfl='no '
         return
         end
C
C
C
C
C
C
C >>>>>> Extract a character data item.
C
         function char_(temp,strg)
C
         logical    char_, charnp_
         include   'ciftbx.sys'
         character  temp*(*), strg*(*)
         integer    lstrg,nstrg

         nstrg = len(strg)
         char_ = charnp_(temp,strg,lstrg)

         if (lstrg.lt.len(strg)) strg(lstrg+1:nstrg) = ' '
         return
         end

C
C
C
C
C
C
C >>>>>> Extract a character data item, no padding.
C
         function charnp_(temp,strg,lstrg)
C
         logical    charnp_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  strg*(*),flag*4
         integer    lstrg
         character*1 slash
         character*4 otype
         integer    ltemp, lname, klow
         integer    lastnb
C
         slash = rsolidus(1:1)
         ltemp = lastnb(temp)
         otype = type_
         call tbxxclc(name,lname,temp,ltemp)
         if(testfl.eq.'yes')    goto 100
         if(.not.text_.and.depth_.eq.0)         goto 120
         if(name(1:1).ne.' '.and.name(1:1).ne.char(0).and.
     *     name(1:lname).ne.nametb(1:lnametb))     goto 120
         charnp_=.false.
         lstrg = 1
         strg=' '
         call getstr
         if (type_.eq.'fini') goto 200
         if (otype.eq.'text' .and. (.not. text_) .and.long_.eq.0) then
           quote_=' '
           textfl = 'no'
           charnp_=.false.
           type_ = 'null'
           goto 200
         end if
         posval_ = jchar-long_
         posend_ = jchar-1
         if (long_.gt.0) then
            strg=strg_(1:long_)
            lstrg=long_
            if (type_.eq.'numb') then
               call ctonum
               if(posdec_.gt.0) posdec_=posval_+posdec_-1
            else 
              if (quote_.eq.' ') then
                 if (long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
                 if (long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
              end if
            end if
         end if
         charnp_=.true.
         goto 200
C
100      if(name(1:lname).eq.nametb(1:lnametb))     goto 150
C
120      quote_=' '
         call tbxxgitm(name(1:lname))
         text_=.false.
         if(type_.eq.'null') then
           charnp_=.false.
           text_=.false.
           textfl = 'no '
           strg_=' '
           long_=0
           goto 200
         endif
C
C        strg_(1:long_) loaded with item
C
150      charnp_=.true.
         strg(1:1)=' '
         lstrg = 1
         if(long_.gt.0) then
           strg=strg_(1:long_)
           lstrg = long_
         endif
         if(type_.eq.'char' )   goto 200
         charnp_=.false.
         if(type_.ne.'text')   goto 200
         charnp_=.true.
         call getlin(flag)
         jchar=MAXBUF+1
         if(flag.eq.'fini')    goto 200
         if(buffer(1:1).eq.';')then
           jchar=2
           textfl = 'no '
           quote_=';'
           goto 200
         endif
         irecd=irecd-1
         text_=.true.
         if (long_.gt.0) then
         if (unfold_ .and. strg(long_:long_).eq.slash) then
170        klow = long_
           long_ = long_-1
           call getlin(flag)
           if(flag.eq.'fini')    goto 210
           if(buffer(1:1).eq.';') then
             jchar=2
             textfl = 'no '
           goto 210
           endif
           quote_=' '
           jchar=lastch+1
           long_=min(len(strg_),klow+max(1,lastch)-1)
           strg_(klow:long_)=buffer(1:max(1,lastch))
           strg(long_:long_)=' '
           if(lastch.gt.0) then
             long_=min(len(strg),klow+lastch-1)
             if(long_.ge.klow) strg(klow:long_)=buffer(1:lastch)
           endif
           if( strg(long_:long_).eq.slash ) go to 170
         endif
         endif
C
200      testfl='no '
         if(long_.eq.0) strg(1:1)=' '
         lstrg = max(1,long_)
CDBG     print *,' Leaving charnp_ text_, type_, quote_: ',
CDBG *     charnp_,text_,type_,quote_
CDBG     print *, ':>>>:'//strg_(1:lstrg)
         if (type_.eq.'char' .and. quote_.eq.' ') then 
           if (strg(1:lstrg).eq.'?'.or. strg(1:lstrg).eq.'.') 
     *       type_='null'
         end if
         return
C
210      text_ = .false.
         go to 200
C
         end
C
C
C
C
C
C >>>>>> Extract a comment or terminal delimiter field
C        backing up to a prior delimiter, depth_ will not
C        be changed even when crossing a terminal delimiter
C
         function cotdb_(strg,lstrg,istd,posstart,recstart)
C
         logical   cotdb_
         logical   istd
         integer   posstart,recstart
         integer   lastnb
         include  'ciftbx.sys'
         character strg*(*),flag*4,c*1,
     *     jjbuf*(MAXBUF)
         integer   jjchar,jjrecd,jjlast,jjlrec,jjjrec
         integer   lstrg
         integer   ipp,itpos
         integer   klow
         character*1 slash
C
         jjchar = jchar
         jjrecd = irecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf=' '
         istd = .false.
         slash = rsolidus(1:1)
         if(lastch.gt.0)jjbuf(1:lastch)=buffer(1:lastch)
         lrecd = nrecd
         if (irecd.ne.recstart) then
           irecd = recstart-1
           call getlin(flag)
           if(flag.eq.'fini') then
             strg='fini'
             jchar=MAXBUF+1
             lstrg=4
             cotdb_=.false.
             posstart = jchar
             recstart = irecd
             go to 300
           endif
         end if
         jchar = posstart
         strg=' '
         lstrg=0
         cotdb_=.false.
100      jchar=jchar+1
         if (jchar.gt.jjchar.and.irecd.ge.jjrecd) go to 300
         if(jchar.le.lastch)     goto 140
C
C....... Read a new line
C
         call getlin(flag)
         if(flag.eq.'fini') then
           cotdb_=.false.
           posstart = jchar
           recstart = irecd
           strg='fini'
           lstrg=4
           go to 300
         endif
         jchar=0
         strg=char(0)
         lstrg=1
         posnam_=0
         quote_=' '
         goto 220
140      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 200
C
C....... Process this character in the line
C
         c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab.and.(.not.tabx_)) goto 190
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 200
         if(depth_.gt.0 .and.
     *    ((c.eq.']'.and.rdbkt_)
     *      .or.(c.eq.')'.and.rdprn_)
     *      .or.(c.eq.'}'.and.rdbrc_))) go to 250
         goto 300
C
C        For a tab, when not expanding to blanks, accept
C        that single character as a comment
C
190      lstrg=1
         strg=tab
         posnam_=jchar
         jchar=jchar+1
         goto 220
C
C....... Accept the remainder of the line as a comment
C
200      lstrg=lastch-jchar
         quote_=buffer(jchar:jchar)
         itpos=jchar
         if(tabx_) then
           itpos=0
           do ipp=1,jchar
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
         endif
         posnam_=itpos
         if(lstrg.gt.0) then
           strg = buffer(jchar+1:lastch)
         endif
         if(lstrg.le.0) then
           strg=' '
           lstrg=1
         endif
         if (strg.eq.slash .and. unfold_) go to 390
         jchar=MAXBUF+1
220      cotdb_=.true.
         posstart = jchar
         recstart = irecd
         go to 300
         
C
C....... Accept the next character as a terminal delimiter
C        in a bracketed construct
C
250      lstrg=1
         quote_=' '
         posval_=jchar
         strg = buffer(jchar:jchar)
         istd =.true.
         cotdb_=.true.
         posstart = jchar
         recstart = irecd

C
C....... restore pointers and exit
C
300      irecd = jjrecd
         lastch = jjlast
         lrecd = jjlrec
         jchar = jjchar
         buffer(1:1)=' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=jjbuf(1:lastch)
         jrecd=jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_=irecd
         return
C
C....... Got a comment with a folding flag
C
390      klow = 1
         lrecd=jjlrec
         cotdb_=.true.
         strg(1:1)=' '
400      jjchar = MAXBUF+1
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         lstrg=0
         go to 420
410      jchar=jchar+1
         if(jchar.le.lastch) go to 450
420      call getlin(flag)
         jchar = 1
         jjchar = 1
         if(flag.eq.'fini') then
           cotdb_=.false.
           strg='fini'
           jchar=MAXBUF+1
           lstrg=lastnb(strg)
           posstart = jchar
           recstart = irecd
           go to 300
         endif
         jchar=1
450      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 400
C
C....... Process this character in the line
C
         c=buffer(jchar:jchar)
         if(c.eq.' '.or.c.eq.tab)      goto 410
         if(c.eq.'#')       goto 470
         posstart = jchar
         recstart = irecd
         goto 300
C
C....... Accept the remainder of the line as part of the comment
C
470      lstrg=lastch-jchar
         itpos=jchar
         if(lastch.gt.jchar)
     *      strg(klow:min(len(strg),klow+lastch-2)) =
     *      buffer(jchar+1:lastch)
         klow=lastnb(strg)
         if (strg(klow:klow).eq.slash) then
           strg(klow:klow)=' '
           go to 400
         endif
         jchar=MAXBUF+1
         lstrg = klow
         lrecd=jjlrec
         posstart = jchar
         recstart = irecd
         goto 300

         end
C
C
C
C >>>>>> Extract a comment or terminal delimiter field.
C
         function cotd_(strg,istd)
C
         logical   cotd_
         logical   istd
         integer   lastnb
         include  'ciftbx.sys'
         character strg*(*),flag*4,c*1,
     *     jjbuf*(MAXBUF)
         integer   jjchar,jjrecd,jjlast,jjlrec,jjjrec
         integer   ipp,itpos
         integer   klow
         character*1 slash
C
         jjchar = jchar
         jjrecd = irecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf=' '
         istd = .false.
         slash = rsolidus(1:1)
         if(lastch.gt.0)jjbuf(1:lastch)=buffer(1:lastch)
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         strg=' '
         long_=0
         cotd_=.false.
         if (depth_.eq.0.and.jchar.gt.0) go to 105
100      jchar=jchar+1
105      if(jchar.le.lastch)     goto 140
C
C....... Read a new line
C
         call getlin(flag)
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=4
           cotd_=.false.
           return
         endif
         jchar=0
         strg=char(0)
         long_=1
         posnam_=0
         quote_=' '
         goto 220
140      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 200
C
C....... Process this character in the line
C
         c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab.and.(.not.tabx_)) goto 190
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 200
         if(depth_.gt.0 .and.
     *    ((c.eq.']'.and.rdbkt_)
     *      .or.(c.eq.')'.and.rdprn_)
     *      .or.(c.eq.'}'.and.rdbrc_))) go to 250
         goto 300
C
C        For a tab, when not expanding to blanks, accept
C        that single character as a comment
C
190      long_=1
         strg=tab
         posnam_=jchar
         jchar=jchar+1
         goto 220
C
C....... Accept the remainder of the line as a comment
C
200      long_=lastch-jchar
         quote_=buffer(jchar:jchar)
         itpos=jchar
         if(tabx_) then
           itpos=0
           do ipp=1,jchar
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
         endif
         posnam_=itpos
         if(long_.gt.0) then
           strg = buffer(jchar+1:lastch)
         endif
         if(long_.le.0) then
           strg=' '
           long_=1
         endif
         if (strg.eq.slash .and. unfold_) go to 390
         jchar=MAXBUF+1
220      lrecd=jjlrec
         cotd_=.true.
         return
         
C
C....... Accept the next character as a terminal delimiter
C        in a bracketed construct
C
250      long_=1
         quote_=' '
         depth_ = depth_-1
         posval_=jchar
         strg = buffer(jchar:jchar)
         jchar=jchar+1
         lrecd=jjlrec
         istd =.true.
         cotd_=.true.
         return

C
C....... Found a non-comment field, restore pointers
C
300      irecd = jjrecd
         lastch = jjlast
         lrecd = jjlrec
         jchar = jjchar
         buffer(1:1)=' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=jjbuf(1:lastch)
         jrecd=jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_=irecd
         return
C
C....... Got a comment with a folding flag
C
390      klow = 1
         lrecd=jjlrec
         cotd_=.true.
         strg(1:1)=' '
400      jjchar = MAXBUF+1
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         long_=0
         go to 420
410      jchar=jchar+1
         if(jchar.le.lastch) go to 450
420      call getlin(flag)
         jchar = 1
         jjchar = 1
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=lastnb(strg)
           return
         endif
         jchar=1
450      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 400
C
C....... Process this character in the line
C
         c=buffer(jchar:jchar)
         if(c.eq.' '.or.c.eq.tab)      goto 410
         if(c.eq.'#')       goto 470
         goto 500
C
C....... Accept the remainder of the line as part of the comment
C
470      long_=lastch-jchar
         itpos=jchar
         if(lastch.gt.jchar)
     *      strg(klow:min(len(strg),klow+lastch-2)) =
     *      buffer(jchar+1:lastch)
         klow=lastnb(strg)
         if (strg(klow:klow).eq.slash) then
           strg(klow:klow)=' '
           go to 400
         endif
         jchar=MAXBUF+1
         long_ = klow
         lrecd=jjlrec
         return
C
C....... Found a non-comment field, restore pointers, but return the
C        comment found so far
C
500      jchar = jjchar
         return

         end
C         
         subroutine tbxxbtab
         include  'ciftbx.sys'
         if (jchar.gt.0 .and. jchar.le.lastch) then
           if (buffer(jchar:jchar).eq.tab
     *      .and..not.tabx_) 
     *      jchar=jchar-1
         end if
         return
         end
C         
         subroutine tbxxetab
         include  'ciftbx.sys'
         if (jchar.gt.1 .and. jchar.le.lastch) then
           if (buffer(jchar-1:jchar-1).eq.tab
     *      .and..not.tabx_) then
            jchar = jchar-1 
            buffer(jchar:jchar) = ' '
            end if
         end if
         return
         end
C
C
C
C
C
C
C >>>>>> Extract a comment field.
C
         function cmnt_(strg)
C
         logical   cmnt_
         integer   lastnb
         include  'ciftbx.sys'
         character strg*(*),flag*4,c*1,
     *     jjbuf*(MAXBUF)
         integer   jjchar,jjrecd,jjlast,jjlrec,jjjrec
         integer   ipp,itpos
         integer   klow
         character*1 slash
C
         jjchar = jchar
         jjrecd = irecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf=' '
         slash = rsolidus(1:1)
         if(lastch.gt.0)jjbuf(1:lastch)=buffer(1:lastch)
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         strg=' '
         long_=0
         cmnt_=.false.
         if (depth_.eq.0 .and. jchar.gt.0) go to 105
100      jchar=jchar+1
105      if(jchar.le.lastch)     goto 140
C
C....... Read a new line
C
         call getlin(flag)
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=4
           cmnt_=.false.
           return
         endif
         jchar=0
         strg=char(0)
         long_=1
         posnam_=0
         quote_=' '
         goto 220
140      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 200
C
C....... Process this character in the line
C
         c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab.and.(.not.tabx_)) goto 190
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 200
         goto 300
C
C        For a tab, when not expanding to blanks, accept
C        that single character as a comment
C
190      long_=1
         strg=tab
         posnam_=jchar
         jchar=jchar+1
         goto 220
C
C....... Accept the remainder of the line as a comment
C
200      long_=lastch-jchar
         quote_=buffer(jchar:jchar)
         itpos=jchar
         if(tabx_) then
           itpos=0
           do ipp=1,jchar
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
         endif
         posnam_=itpos
         if(long_.gt.0) then
           strg = buffer(jchar+1:lastch)
         endif
         if(long_.le.0) then
           strg=' '
           long_=1
         endif
         if (strg.eq.slash .and. unfold_) go to 390
         jchar=MAXBUF+1
220      lrecd=jjlrec
         cmnt_=.true.
         return
C
C....... Found a non-comment field, restore pointers
C
300      irecd = jjrecd
         lastch = jjlast
         lrecd = jjlrec
         jchar = jjchar
         buffer(1:1)=' '
         if(lastch.gt.0)buffer(1:min(MAXBUF,lastch))=jjbuf(1:lastch)
         jrecd=jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_=irecd
         return
C
C....... Got a comment with a folding flag
C
390      klow = 1
         lrecd=jjlrec
         cmnt_=.true.
         strg(1:1)=' '
400      jjchar = MAXBUF+1
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         long_=0
         go to 420
410      jchar=jchar+1
         if(jchar.le.lastch) go to 450
420      call getlin(flag)
         jchar = 1
         jjchar = 1
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=lastnb(strg)
           return
         endif
         jchar=1
450      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 400
C
C....... Process this character in the line
C
         c=buffer(jchar:jchar)
         if(c.eq.' '.or.c.eq.tab)      goto 410
         if(c.eq.'#')       goto 470
         goto 500
C
C....... Accept the remainder of the line as part of the comment
C
470      long_=lastch-jchar
         itpos=jchar
         if(lastch.gt.jchar)
     *      strg(klow:min(len(strg),klow+lastch-2)) =
     *      buffer(jchar+1:lastch)
         klow=lastnb(strg)
         if (strg(klow:klow).eq.slash) then
           strg(klow:klow)=' '
           go to 400
         endif
         jchar=MAXBUF+1
         long_ = klow
         lrecd=jjlrec
         return
C
C....... Found a non-comment field, restore pointers, but return the
C        comment found so far
C
500      jchar = jjchar
         return

         end
C
C
C
C
C
C
C >>>>>> Return the delimiter prior to the most recently
C        examined value
C
         function delim_(depth,delim,posdlm,recdlm)
C
         logical   delim_
         integer   depth
         integer   posdlm
         integer   recdlm
         character*(*) delim

         include  'ciftbx.sys'
         delim = ' '
         delim_ = .false.
         posdlm = 0
         if (depth .ge.0 .and. depth .le.depth_) then
           delim = delimstack(depth+1)
           delim_ = .true.
           posdlm = posdlmstk(depth+1)
           recdlm = recdlmstk(depth+1)
         end if

         return

         end
C
C
C
C
C
C >>>>> Convert name string to lower case
C
         function tbxxlocs(name)
C
         include     'ciftbx.sys'
         character    tbxxlocs*(MAXBUF)
         character    temp*(MAXBUF),name*(*)
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
C
         temp=name
         kln = lastnb(name)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(cap,c)
         if(j.ne.0) temp(i:i)=low(j:j)
100      continue
         tbxxlocs=temp
         return
         end
C
C
C
C
C
C >>>>> Convert name string to lower case as subroutine
C
         subroutine tbxxnlc(loname, name)
C
         include     'ciftbx.sys'
         character    temp*(MAXBUF),loname*(*),name*(*)
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      lolen,olen
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
C
         lolen = len(loname)
         olen = len(name)
         kln = min(MAXBUF,lolen,olen)
         kln = lastnb(name(1:kln))
         temp(1:kln)=name(1:kln)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(cap,c)
         if(j.ne.0) temp(i:i)=low(j:j)
100      continue
         loname=temp(1:kln)
         return
         end
C
C
C
C
C
C >>>>> Convert counted name string to lower case as subroutine
C       with counts
C
         subroutine tbxxclc(loname, lloname, name, lname)
C
         include     'ciftbx.sys'
         character    temp*(MAXBUF),loname*(*),name*(*)
         integer      lloname, lname
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      lolen,olen
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
C
         lolen = len(loname)
         olen = min(len(name),lname)
         kln = min(MAXBUF,lolen,olen)
         kln = lastnb(name(1:kln))
         temp(1:kln)=name(1:kln)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(cap,c)
         if(j.ne.0) temp(i:i)=low(j:j)
100      continue
         loname(1:kln)=temp(1:kln)
         lloname = kln
         return
         end

C
C
C
C
C
C >>>>> Convert name string to upper case
C
         function tbxxupcs(name)
C
         include     'ciftbx.sys'
         character    tbxxupcs*(MAXBUF)
         character    temp*(MAXBUF),name*(*)
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
C
         temp=name
         kln = lastnb(name)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(low,c)
         if(j.ne.0) temp(i:i)=cap(j:j)
100      continue
         tbxxupcs=temp
         return
         end
C
C
C
C
C
C >>>>> Convert name string to upper case as subroutine
C
         subroutine tbxxnupc(upname, name)
C
         include     'ciftbx.sys'
         character    temp*(MAXBUF),upname*(*),name*(*)
         character    low*26,cap*26,c*1
         integer      i,j,kln
         integer      olen,uplen
         integer      lastnb
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
C
         uplen = len(upname)
         olen = len(name)
         kln = min(MAXBUF,uplen,olen)
         kln = lastnb(name(1:kln))
         temp(1:kln)=name(1:kln)
         do 100 i=1,kln
         c=temp(i:i)
         j=index(low,c)
         if(j.ne.0) temp(i:i)=cap(j:j)
100      continue
         upname=temp(1:kln)
         return
         end

C
C
C
C
C
C >>>>>> Get the data item associated with the tag.
C
         subroutine tbxxgitm(name)
C
         include   'ciftbx.sys'
         SAVE
         character name*(*)
         character flag*4
         character*1 slash
         integer   iitem,nitem,npakt
         integer   kchar,loopi,i,jdict,itpos,ipp
         integer   lastnb
C
         slash = rsolidus(1:1)
C
C....... Find requested dataname in hash list
C
         lnametb=lastnb(name)
         nametb(1:lnametb)=name(1:lnametb)
CDBG     print *,' Entering tbxxgitm: ', name(1:lnametb),' ',
CDBG *     tcheck, vcheck
         posnam_=0
         posval_=0
         posdec_=0
         posend_=0
         valid_ = .false.
         quote_=' '
         jdict = 0
         strg_= ' '
         long_=1
         if(name(1:1).eq.'_')       goto 100
         type_='null'
         ttype_='    '
         depth_=0
         index_=0
         dictype_='null'
         diccat_='(none)'
         dicname_=name
         tagname_=' '
         goto 1000
100      call hash_find(nametb(1:lnametb),
     *     dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,
     *     iname)
         if(iname.gt.0)             goto 180
         if(dictfl.ne.'yes') then
         call hash_find(nametb(1:lnametb),
     *     dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,jdict)
         if(jdict.ne.0) then
CDBG        print *,' found entry ', jdict, dicxtyp(jdict)
           dictype_=dicxtyp(jdict)
           if(dcindex(jdict).ne.0) diccat_=dcname(dcindex(jdict))
           dicname_=nametb(1:lnametb)
           if(aroot(jdict).ne.0) then
             dicname_=dictag(aroot(jdict))
             call hash_find(dicnam(aroot(jdict)),
     *         dname,dchain,NUMBLOCK,nname,dhash,NUMHASH,
     *         iname)
             if(iname.gt.0)      goto 180
           endif
           type_='null'
           ttype_='    '
           depth_=0
           index_=0
           tagname_=' '
           strg_=' '
           long_=1
           go to 1000
         endif
         endif

         continue
         type_='null'
         ttype_='    '
         depth_=0
         index_=0
         dictype_='null'
         diccat_='(none)'
         dicname_=name
         long_=1
         goto 1000
C
C
180      tagname_=dtag(iname)
         if(ddict(iname).ne.0) tagname_=dictag(ddict(iname))
         posnam_=tchar(iname)
         if(tabx_)posnam_=xchar(iname)
         if(nloop(iname).le.0)      goto 500
C
C....... Process loop packet if first item request
C
         if(nloop(iname).ne.loopnl) goto 200
         if(lloop(iname).lt.loopct) goto 300
         if(loop_)                  goto 230
200      loop_=.true.
         depth_ = 0
         loopct=0
         loopnl=nloop(iname)
         nitem=loopni(loopnl)
         npakt=loopnp(loopnl)
         irecd=drecd(iname)-1
         call getlin(flag)
         jchar=max(0,dchar(iname)-1)
CDBG     if(jchar.lt.0) write(6,'(7H dchar ,i5)') jchar
         do 220 i=1,nitem
220      lloop(i+iname-iloop(iname))=0
         goto 240
C
C....... Read a packet of loop items
C
230      nitem=loopni(loopnl)
         npakt=loopnp(loopnl)
         irecd=looprd(nitem+1)-1
         call getlin(flag)
         jchar=loopch(nitem+1)
CDBG     if(jchar.lt.0) write(6,'(7H loopch,i5)') jchar
240      iitem=0
250      iitem=iitem+1
         quote_=' '
         text_=.false.
         if(iitem.le.nitem)     goto 255
         loopch(iitem)=jchar
         looprd(iitem)=irecd
         goto 270
255      call getstr
         loopch(iitem)=jchar-long_
         if(quote_.ne.' ') then
           if (quote_.eq.';') then
             loopch(iitem)=1
           else
             if (quote_.eq.''''''''.or.quote_.eq.'"""') then
               loopch(iitem)=jchar-long_-3
             else
               loopch(iitem)=jchar-long_-1
             end if
           end if
         end if
         loopln(iitem)=long_
         looprd(iitem)=irecd
         
         if (text_ .or.depth_ .gt. 0) then
         if (depth_.gt.0) then
           loopch(iitem)= posbrkstk(1)
           loopln(iitem)= 1
           looprd(iitem)= srecd
         end if
260        call getstr
           if (type_.eq.'fini') call tbxxerr(' Unexpected end of data')
           if (text_.or.depth_ .gt. 0) goto 260
         end if
         goto 250
270      loopct=loopct+1
         if(loopct.lt.npakt)    goto 300
         loop_=.false.
C
C....... Point to the loop data item
C
300      lloop(iname)=lloop(iname)+1
         loopi=iloop(iname)
         irecd=looprd(loopi)-1
         call getlin(flag)
         long_=loopln(loopi)
         kchar=loopch(loopi)
         if ((buffer(kchar:kchar).eq.'(' .and. rdprn_)
     *     .or. (buffer(kchar:kchar).eq.'[' .and. rdbkt_)
     *     .or. (buffer(kchar:kchar).eq.'{' .and. rdbrc_)) then
           if (kchar.gt.1) then
             if (buffer(kchar-1:kchar-1).eq.''''
     *          .or. buffer(kchar-1:kchar-1).eq.'"') goto 550
           end if
           jchar = kchar-1
           call getstr
CDBG         print *,' strg_ ', strg_(1:max(1,long_))
CDBG         print *,' depth_ ', depth_
           itpos=jchar-long_
           posval_=itpos
           posend_=itpos+long_-1
           jchar=kchar+long_
           if(jchar.le.MAXBUF) then
             if(buffer(jchar:jchar).ne.' ' .and.
     *       buffer(jchar:jchar).ne.tab)
     *       jchar=jchar+1
           endif
           if(type_.eq.'numb') then
             call ctonum
             if(posdec_.gt.0) posdec_=posval_+posdec_-1
           endif
           go to 1000
         end if
         goto 550
C
C....... Point to the non-loop data item
C
500      irecd=drecd(iname)-1
         call getlin(flag)
         kchar=dchar(iname)+1
         long_=iloop(iname)
         loop_=.false.
         loopct=0
         loopnl=0
C
C....... Place data item into variable string and make number
C
550      type_=dtype(iname)
         quote_=' '
         text_=.false.
         dictype_=dxtyp(iname)
         diccat_='(none)'
         if(cindex(iname).gt.0) diccat_=dcname(cindex(iname))
         if(cindex(iname).lt.0) diccat_=cname(-cindex(iname))
         if(diccat_.eq.' ') diccat_='(none)'
         dicname_=dtag(iname)
         if(ddict(iname).ne.0) then
           if (aroot(ddict(iname)).ne.0) then
             dicname_=dictag(aroot(ddict(iname)))
           endif
         endif
         strg_=' '
         if(long_.gt.0) then
CDBG       if (kchar.le.0)call tbxxwarn(' kchar le 0')
           strg_(1:long_)=buffer(kchar:kchar+long_-1)
           if ((buffer(kchar:kchar).eq.'('.and.rdprn_)
     *       .or. (buffer(kchar:kchar).eq.'['.and.rdbkt_)
     *       .or. (buffer(kchar:kchar).eq.'{'.and.rdbrc_)) then
             if (kchar.gt.1) then
               if (buffer(kchar-1:kchar-1).eq.''''.or.
     *           buffer(kchar-1:kchar-1).eq.'"') then
                 go to 555
               end if        
             end if
CDBG         print *,' getitm: kchar, irecd ',kchar,irecd
             jchar = kchar-1
CDBG         print *,' strg_ ', strg_(1:max(1,long_))
             call getstr
             if(type_.eq.'numb') then
               call ctonum
               if(posdec_.gt.0) posdec_=posval_+posdec_-1
             endif
CDBG         print *,' strg_ ', strg_(1:max(1,long_))
CDBG         print *,' depth_ ', depth_
             go to 1000
           end if
         endif
555      itpos=kchar
         posval_=itpos
         posend_=itpos+long_-1
         jchar=kchar+long_
         if(jchar.le.MAXBUF) then
           if(buffer(jchar:jchar).ne.' ' .and.
     *       buffer(jchar:jchar).ne.tab) jchar=jchar+1
         endif
         quote_=' '
         if(kchar.gt.1) then
           if(buffer(kchar-1:kchar-1).ne.' ' .and.
     *        buffer(kchar-1:kchar-1).ne.tab) then
             quote_=buffer(kchar-1:kchar-1)
             if (kchar.gt.3.and.rdtq_) then
               if (buffer(kchar-3:kchar-1).eq.
     *           quote_//quote_//quote_) then
                 quote_ = buffer(kchar-3:kchar-1)
               end if
             end if
           endif
         endif
         if(type_.eq.'char' .and. kchar.eq.1 .and.
     *     buffer(1:1).eq.';') then
           type_='text'
           fold_=.false.
           quote_=';'
         endif
         if(type_.eq.'text') then
           if(buffer(1:1).eq.';') then
             quote_=';'
             if (clipt_.or.long_.lt.2) then
             strg_(1:1)=' '
             if (strg_(1:long_).eq.(' '//slash) ) then
               fold_=.true.
               if(unfold_) then
                 strg_(1:long_)=slash
                   long_=1
                 endif
               endif
             else
               do ipp = 2,long_
                 strg_(ipp-1:ipp-1)=strg_(ipp:ipp)
               end do
               long_=long_-1
               if (strg_(1:long_).eq.slash) then
                 fold_=.true.
                 if (unfold_) then
                 long_=1
                 endif
               endif
             endif
           else
             type_='char'
             if (quote_.eq.';') quote_=' '
           endif
         endif
         if(type_.eq.'numb') then
           call ctonum
           if(posdec_.gt.0) posdec_=posval_+posdec_-1
         endif
         if(type_.eq.'char' .and. strg_.eq.' '.and.nblank_)
     *     type_='null'
         if (quote_.ne.' ') goto 1000
         if (long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
         if (long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
         if (tcheck.eq.'yes') then
           call hash_find(nametb(1:lnametb),
     *       dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,jdict)
           if (jdict.gt.0)
     *     call tbxxckv(jdict)
         endif
C
1000     return
         end

C
C
C
C
C
C
C
C >>>>>> Convert string to integer, marking non-digit
C
C
         function tbxxsti(xstr,nondig)
         integer tbxxsti
         character *(*) xstr
         integer nondig, i
         integer sign, digits, kdv

         tbxxsti = 0
         digits = 0
         nondig = 0
         sign = 1
         do i = 1,len(xstr)
           kdv = ichar(xstr(i:i))-ichar('0')
           if (digits.eq.0) then
             if (xstr(i:i).eq.'-') then
               sign = -1
               digits = 1
             else
               if (xstr(i:i).eq.'+') then
                 sign = 1
                 digits = 1
               else
                 if (kdv.ge.0 .and. kdv.le.9) then
                   digits = 1
                   tbxxsti = kdv
                 else
                   if (xstr(i:i).ne.' ') then
                     nondig = i
                     return
                   endif
                 endif
               endif
             endif
           else
             if (kdv.ge.0 .and.kdv.le.9) then
               tbxxsti = tbxxsti*10+kdv
             else
               tbxxsti = sign*tbxxsti
               nondig = i
               return
             endif
           endif
         enddo
         tbxxsti = sign*tbxxsti
         return

         end



C
C
C
C
C
C
C
C >>>>>> Convert string to double, marking non-digit
C
C
         function tbxxstd(xstr,nondig)
         double precision tbxxstd
         integer tbxxsti
         character *(*) xstr
         integer nondig, i
         integer sign, digits, kdv
         integer idp, eval
         tbxxstd = 0.0
         digits = 0
         nondig = 0
         sign = 1
         idp = 0
         do i = 1,len(xstr)
           kdv = ichar(xstr(i:i))-ichar('0')
           if (i.lt.len(xstr)
     *       .and. (xstr(i:i).eq.'e'
     *       .or. xstr(i:i).eq.'E'
     *       .or. xstr(i:i).eq.'d'
     *       .or. xstr(i:i).eq.'D'
     *       .or. xstr(i:i).eq.'q'
     *       .or. xstr(i:i).eq.'Q')) then
             eval = tbxxsti(xstr(i+1:len(xstr)),nondig)
             tbxxstd = sign*tbxxstd*10.**eval
             if (nondig.ne.0) nondig=nondig+i+1
             return
           endif
           if (i.lt.len(xstr) .and. digits .ne.0
     *       .and. (xstr(i:i).eq.'+'
     *       .or. xstr(i:i).eq.'-')) then
             eval = tbxxsti(xstr(i:len(xstr)),nondig)
             tbxxstd = sign*tbxxstd*10.**eval
             if (nondig.ne.0) nondig=nondig+i
             return
           endif
           if (xstr(i:i).eq.'.'.and.idp.eq.0) then
             idp = i
             digits = 1
           endif
           if (digits.eq.0) then
             if (xstr(i:i).eq.'-') then
               sign = -1
               digits = 1
             else
               if (xstr(i:i).eq.'+') then
                 sign = 1
                 digits = 1
               else
                 if (kdv.ge.0 .and. kdv.le.9) then
                   digits = 1
                   tbxxstd = kdv
                 else
                   if (xstr(i:i).ne.' ') then
                     nondig = i
                     return
                   endif
                 endif
               endif
             endif
           else
             if (kdv.ge.0 .and.kdv.le.9) then
               if (idp.eq.0) then
                 tbxxstd = tbxxstd*10.+kdv
               else
                 tbxxstd = tbxxstd+kdv*(10.**(idp-i))
               endif
             else
               if (i.ne.idp) then
                 tbxxstd = sign*tbxxstd
                 nondig = i
                 return
               endif
             endif
           endif
         enddo
         tbxxstd = sign*tbxxstd
         return

         end


C
C
C
C
C
C
C
C >>>>>> Validate the string in strg_(1:long_) of type type_
C        against the dictionary item at jdict
C
C
         subroutine tbxxckv(jdict)
         integer jdict
C
         include   'ciftbx.sys'

         character*(MAXBUF) temp, target, lcvalue
         integer tbxxfstb
         integer tbxxsti
         double precision tbxxstd
         integer tlen
         logical igood, isword, nolo, nohi
         integer fblank,ftab, symop, xlate
         integer yyyy, mm, dd, hr, mi, se, sf, tz
         integer nondig, prevdig, ldt, ldn
         logical enumflg
         integer lastnb
         integer kptr, ltarget, icptr
         integer ilolo, ilohi, ihilo, ihihi
         integer llcvalue

         valid_ = .false.
         igood = .false.
         enumflg = .false.
         kptr = 0
         if (long_ .lt. 1) return

         fblank = index(strg_(1:long_),' ')
         ftab = index(strg_(1:long_),tab)
         ldt = max(1,lastnb(dictype_))
         ldn = max(1,lastnb(dicnam(jdict)))
         isword = .true.
         if (fblank.ne.0 .or. ftab.ne.0) isword =.false.

         if (type_.eq.'null') igood = .true.

         if ((type_.eq.'char' .or. type_.eq. 'numb').and. isword) then
           if (dictype_.eq.'uchar3') then
             if (long_.eq.3.or.
     *            (long_.eq.4.and.strg_(1:1).eq.'+')
     *            ) igood = .true.
             go to 90
           endif

           if (dictype_.eq.'uchar1') then
             if (long_.eq.1.or.
     *            (long_.eq.2.and.strg_(1:1).eq.'+')
     *          ) igood = .true.
             go to 90
           endif

           if (dictype_(1:4).eq.'symo') then
             symop = tbxxsti(strg_(1:long_),nondig)
             xlate = 0
             if (nondig.ne.0.and.nondig.lt.long_) then
               if (strg_(nondig:nondig).eq.'_') then
                 xlate = tbxxsti(strg_(nondig+1:long_),nondig)
               endif
             endif

             if (nondig.eq.0 .and.
     *         symop .ge. 1 .and.
     *         symop .le. 192 .and.
     *         xlate .ge. 0 .and.
     *         xlate .le. 1000) igood =.true.
             go to 90
           endif

           if (dictype_(1:5).eq.'yyyy-') then

             mm=-1
             dd=-1
             hr=0
             mi =0
             se=0
             sf=0
             tz = 0

             yyyy = tbxxsti(strg_(1:long_),nondig)
             if (nondig.ne.0.and.nondig.lt.long_) then
             if (strg_(nondig:nondig).eq.'-') then
               prevdig = nondig
               mm = tbxxsti(strg_(nondig+1:long_),nondig)
               if (nondig.ne.0) nondig=prevdig+nondig
               if (nondig.ne.0.and.nondig.lt.long_) then
               if (strg_(nondig:nondig).eq.'-') then
                 prevdig = nondig
                 dd = tbxxsti(strg_(nondig+1:long_),nondig)
                 if (nondig.ne.0) nondig=prevdig+nondig
                 if (nondig.ne.0.and.nondig.lt.long_) then
                 if (strg_(nondig:nondig).eq.'T'
     *             .or. strg_(nondig:nondig).eq.'t'
     *             .or. strg_(nondig:nondig).eq.':') then
                   prevdig = nondig
                   hr = tbxxsti(strg_(nondig+1:long_),nondig)
                   if (nondig.ne.0) nondig=prevdig+nondig
                   if (nondig.ne.0.and.nondig.lt.long_) then
                   if (strg_(nondig:nondig).eq.':') then
                     prevdig = nondig
                     mi = tbxxsti(strg_(nondig+1:long_),nondig)
                     if (nondig.ne.0) nondig=prevdig+nondig
                     if (nondig.ne.0.and.nondig.lt.long_) then
                     if (strg_(nondig:nondig).eq.':') then
                       prevdig = nondig
                       se = tbxxsti(strg_(nondig+1:long_),nondig)
                       if (nondig.ne.0) nondig=prevdig+nondig
                       if (nondig.ne.0.and.nondig.lt.long_) then
                       if (strg_(nondig:nondig).eq.'.') then
                         prevdig = nondig
                         sf = tbxxsti(strg_(nondig+1:long_),nondig)
                         if (nondig.ne.0) nondig=prevdig+nondig
                       endif
                       endif
                     endif
                     endif
                   endif
                   endif
                 endif
                 endif
               endif
               endif
             endif
             endif
             if (nondig.ne.0) then
             if (strg_(nondig:nondig).eq.'-'
     *         .or. strg_(nondig:nondig).eq.'+') then
               tz = tbxxsti(strg_(nondig+1:long_),nondig)
             endif
             endif
             if (nondig.eq.0
     *         .and. yyyy .ge. 0  .and. yyyy .lt. 10000
     *         .and. mm .gt. 0 .and. mm .lt. 13
     *         .and. dd .gt. 0 .and. dd .lt. 32
     *         .and. hr .ge. 0 .and. hr .lt. 25
     *         .and. mi .ge. 0 .and. mi .lt. 61
     *         .and. se .ge. 0 .and. se .lt. 61
     *         .and. sf .ge. 0
     *         .and. tz .ge. 0 .and. tz .lt. 25 ) igood =.true.
             go to 90
           endif

           if (dictype_(1:4).eq.'char'
     *       .or. dictype_(1:4).eq.'ucha'
     *       .or. dictype_(1:4).eq.'code'
     *       .or. dictype_(1:4).eq.'ucod'
     *       .or. dictype_(1:4).eq.'line'
     *       .or. dictype_(1:4).eq.'ulin'
     *       .or. dictype_(1:3).eq.'any'
     *       .or. dictype_(1:4).eq.'atco'
     *       .or. dictype_(1:4).eq.'phon'
     *       .or. dictype_(1:4).eq.'emai'
     *       .or. dictype_(1:4).eq.'fax'
     *       .or. dictype_(1:4).eq.'text')  then
             igood = .true.
             go to 90
           endif

           if (dictype_(1:4).eq.'numb'
     *       .or. dictype_(1:3).eq.'int'
     *       .or. dictype_(1:4).eq.'floa') then
             tbxxintr = tbxxsti(strg_(1:long_),nondig)
             if (nondig.eq.0) then
               igood = .true.
               go to 90
             endif
             if (strg_(nondig:nondig).eq.'('
     *         .and. nondig .lt. long_) then
               tbxxintr = tbxxsti(strg_(nondig+1:long_),nondig)
               if (nondig.gt.0) then
                 if (strg_(nondig:nondig).eq.')') then
                   igood = .true.
                   go to 90
                 endif
               endif
             endif
             if (dictype_(1:4).eq.'numb'
     *         .or. dictype_(1:4).eq.'floa') then
               if (type_.eq.'numb') igood = .true.
             endif
             go to 90
           endif
           go to 90
         endif

         if (type_.eq.'char') then
           if (dictype_(1:4).eq.'text'
     *       .or. dictype_(1:3).eq.'any'
     *       .or. dictype_(1:4).eq.'line'
     *       .or. dictype_(1:4).eq.'ulin'
     *       .or. dictype_(1:4).eq.'phon'
     *       .or. dictype_(1:4).eq.'atco'
     *       .or. dictype_(1:4).eq.'phon'
     *       .or. dictype_(1:4).eq.'char'
     *       .or. dictype_(1:4).eq.'ucha' ) igood = .true.
           go to 90
         endif

         if (type_.eq.'text') then
           if (dictype_(1:4).eq.'text'
     *       .or. dictype_(1:3).eq.'any'
     *       .or. dictype_(1:4).eq.'char'
     *       .or. dictype_(1:4).eq.'ucha' ) igood = .true.
           go to 90
         endif


 90      continue

         if (.not.igood) then
           call tbxxwarn(' Dictionary type '//dictype_(1:ldt)//
     *     ' for '//dicnam(jdict)(1:ldn)//
     *     ' not matched by '//strg_(1:long_))
           return
         endif
         kptr = deindex(jdict)
         if (kptr.eq.0 .or. type_.eq.'null') then
           valid_ = .true.
           return
         endif



         call tbxxclc(lcvalue,llcvalue,strg_,long_)

100      if (kptr.ne.0) then

           tlen = tbxxfstb(temp,ivtsbp(kptr),.false.)
           if (tlen.gt.0) then
             call tbxxclc(target,ltarget,temp,tlen)
             if (ivtvet(kptr) .eq. 0) then
               enumflg = .true.
               if (target(1:ltarget).eq.lcvalue(1:llcvalue)) then
                 valid_ = .true.
                 return
               endif
               if(type_.eq.'numb'
     *           .and. (dictype_(1:4).eq.'numb'
     *              .or. dictype_(1:3).eq.'int'
     *              .or. dictype_(1:4).eq.'floa')) then
                 if (tbxxstd(target(1:ltarget),nondig)
     *              .eq.numbtb) then
                    valid_= .true.
                    return
                 endif
               endif
             else
               enumflg = .false.
               icptr = index(target(1:ltarget),':')
               ilolo = 1
               ilohi = icptr-1
               ihilo = icptr+1
               ihihi = ltarget
               nolo = .true.
               if (ilohi.ge.ilolo) then
                 nolo = .false.
                 if (target(ilolo:ilohi).eq.'.')
     *             nolo = .true.
               endif
               nohi = .true.
               if (ihihi.ge.ihilo) then
                 nohi = .false.
                 if (target(ihilo:ihihi).eq.'.')
     *             nohi = .true.
               endif
              if (dictype_(1:4).eq.'numb'
     *         .or. dictype_(1:3).eq.'int'
     *         .or. dictype_(1:4).eq.'floa') then
              if (nolo.and.(.not.nohi)) then
                 if ((ivtvet(kptr).gt.0
     *             .and. numbtb .lt.
     *               tbxxstd(target(ihilo:ihihi),nondig)) .or.
     *               (ivtvet(kptr).lt.0
     *             .and. numbtb .le.
     *               tbxxstd(target(ihilo:ihihi),nondig))) then
                   valid_= .true.
                   return
                 endif
               endif
               if (nohi.and.(.not.nolo)) then
                 if ((ivtvet(kptr).gt.0
     *             .and. numbtb .gt.
     *               tbxxstd(target(ilolo:ilohi),nondig)) .or.
     *               (ivtvet(kptr).lt.0
     *             .and. numbtb .ge.
     *               tbxxstd(target(ilolo:ilohi),nondig))) then
                   valid_= .true.
                   return
                 endif
               endif
               if ((.not.nohi).and.(.not.nolo)) then
                 if ((ivtvet(kptr).gt.0
     *             .and. numbtb .lt.
     *               tbxxstd(target(ihilo:ihihi),nondig)
     *             .and. numbtb .gt.
     *               tbxxstd(target(ilolo:ilohi),nondig)) .or.
     *               (ivtvet(kptr).lt.0
     *             .and. numbtb .le.
     *               tbxxstd(target(ihilo:ihihi),nondig)
     *             .and. numbtb .ge.
     *               tbxxstd(target(ilolo:ilohi),nondig))) then
                   valid_= .true.
                   return
                 endif
               endif
               else
               if (nolo.and.(.not.nohi)) then
                 if ((ivtvet(kptr).gt.0
     *             .and. lcvalue(1:llcvalue) .lt.
     *               target(ihilo:ihihi)) .or.
     *               (ivtvet(kptr).lt.0
     *             .and. lcvalue(1:llcvalue) .le.
     *               target(ihilo:ihihi))) then
                   valid_= .true.
                   return
                 endif
               endif
               if (nohi.and.(.not.nolo)) then
                 if ((ivtvet(kptr).gt.0
     *             .and. lcvalue(1:llcvalue) .gt.
     *               target(ilolo:ilohi)) .or.
     *               (ivtvet(kptr).lt.0
     *             .and. lcvalue(1:llcvalue) .ge.
     *               target(ilolo:ilohi))) then
                   valid_= .true.
                   return
                 endif
               endif
               if ((.not.nohi).and.(.not.nolo)) then
                 if ((ivtvet(kptr).gt.0
     *             .and. lcvalue(1:llcvalue) .lt.
     *               target(ihilo:ihihi)
     *             .and. lcvalue(1:llcvalue) .gt.
     *               target(ilolo:ilohi)) .or.
     *               (ivtvet(kptr).lt.0
     *             .and. lcvalue(1:llcvalue) .le.
     *               target(ihilo:ihihi)
     *             .and. lcvalue(1:llcvalue) .ge.
     *               target(ilolo:ilohi))) then
                   valid_= .true.
                   return
                 endif
               endif
              endif
             endif
           endif
           kptr = ivtnxt(kptr)
           go to 100
         endif
         continue
         if (enumflg) then
           call tbxxwarn(' Dictionary type '//dictype_(1:ldt)//
     *      ' for '//dicnam(jdict)(1:ldn)//', '//
     *      strg_(1:long_)//
     *      ' not in dictionary list of values')
         else
         call tbxxwarn(' Dictionary type '//dictype_(1:ldt)//
     *      ' for '//dicnam(jdict)(1:ldn)//
     *      ' range not matched by '//strg_(1:long_))
         endif

         return
         end
C
C
C
C
C
C
C
C
C >>>>>> Test for separator
C
C
       logical function tbxxtsts(c)
       include   'ciftbx.sys'
       character*1 c
       tbxxtsts = .true.
       if (rdrcqt_) return
       if (c.eq.' ') return
       if (c.eq.tab) return
       tbxxtsts = .false.
       if (depth_ .eq. 0) return
       if (c.eq.',' .or.
     *   c.eq.':' .or.
     *   ((c.eq.')'.or.c.eq.'(') .and. rdprn_) .or.
     *   ((c.eq.'}'.or.c.eq.'{') .and. rdbrc_) .or.
     *   ((c.eq.']'.or.c.eq.'[') .and. rdbkt_) ) then
         tbxxtsts=.true.
         return
       end if
       return
       end
C
C
C
C
C
C
C
C
C >>>>>> Test for terminal treble quote
C
C      tests buffer(jchar:lastch) for a terminal treble quote
C      and returns the location in jtloc
C
C
       logical function tbxxtttq(jtloc)
       include   'ciftbx.sys'
       character*1 slash
       logical escaped
       logical tbxxtsts
       integer i,jtloc

       slash = rsolidus(1:1)
       tbxxtttq = .false.
       jtloc = jchar

       if (rdrcqt_) then

C        Process according to CIF2 rules on quotes
         escaped = .false.
         if (jchar .le. lastch-2) then
           do i = jchar,lastch-2
             if (.not.escaped) then
               if (buffer(i:i).eq.slash) then
                 escaped = .true.
               else
                 if (buffer(i:i+2).eq.quote_) then
                   jtloc = i
                   tbxxtttq = .true.
                   return
                 end if
               end if
             else
               escaped=.false.
             end if
           end do
         end if
         return
       else
C
C        Process according to CIF1 rules
         if (jchar .le. lastch-2) then
           do i = jchar,lastch-2
             if (buffer(i:i+2).eq.quote_) then
               if (i.lt.lastch-2) then
                 if (tbxxtsts(buffer(i+3:i+3))) then
                   jtloc = i
                   tbxxtttq = .true.
                   return
                 end if
               else
                 jtloc = i
                 tbxxtttq = .true.
                 return
               end if
             end if
           end do
         end if
         return
       end if
         end


C
C
C
C
C
C
C
C >>>>>> Read the next string from the file
C
C
       subroutine getstr
C
C      On entry, jchar is set to one less than the next character
C      to be read, on the line given by irecd, which is assumed
C      to have been loaded into buffer, with lastch set to the
C      position of the last character
C
C      if depth_ is greater than 0, then statestack(depth_),
C      brackstack(depth_) and indexstack(depth_) give the
C      state of the scan within a list, array, tuple of table
C
C      If the state is 0, we are starting a search for a token
C      if the state is 2, we had a token last pass and need
C      to find a comma, colon or a terminating ) ] or } before looking
C      for the next token.
C
C      On entry, the state of text_ is used to determine if we are continuing
C      a text field or a bracketed construct.  The case of a text field within
C      a bracketed construct is handled by checking quote_.  If it is any of
C      ';', "'''", or '"""', then the next line needs to be read if text_ is
C      true.  If it any any other string, then the next element of the
C      bracketed construct needs to be read
C
C      If the depth is zero, text_ will be cleared one getstr call
C      prior to the last, empty read.  This change is not made for treble
C      quoted strings
C
C      In a bracketed construct text_ will be left set after the read
C      of the last read with text and the next read will return a null
C      type_
C
       include   'ciftbx.sys'
       integer   i,j,jj(11),im,ip
       integer   state,innerdepth
       logical   quoted
       logical   escaped
       character c*1,num*21,flag*4
       character slash*1
       logical   tbxxtsts, tbxxtttq
       integer   jtloc

       data num/'0123456789+-.()EDQedq'/
       slash = rsolidus(1:1)
       im = 0
CDBG   print *,' entering getstr type, text_, quote_, depth_ jchar: '
CDBG   print *,  type_, text_, quote_, depth_, jchar


C      If text_ is true, we may be continuing a text field, a
C      treble-quoted string or a bracketed construct
C
C      We deal with the first 2 cases here
C
       if (text_) then
         if (quote_.eq.';'.or.quote_.eq.'"""'.or.quote_.eq."'''") then
CDBG     print *,' processing next line'
           call getlin(flag)
           if (flag.eq.'fini') then
             type_='fini'
             text_=.false.
             depth_=0
             ttype_='    '
             depth_=0
             quote_=' '
             quoted=.false.
             goto 500
           end if
C
C          Handle the case of a text field
C          This is terminated by \n; which
C          is unconditionally recognized if rdrcqt_ is true
C          or which requires a trailing separator if depth_
C          is greater than 0

CDBG   print *, 'read line:',buffer(1:lastch)
           if (quote_.eq.';') then
CDBG     print *, 'semicolon quote_ detected'
             if (buffer(1:1).eq.';') then
CDBG         print *, 'terminal semicolon detected'
               if (lastch.gt.1) then
                 if (.not.tbxxtsts(buffer(2:2))) goto 10
               end if
               jchar = 2
               long_ = 0
               strg_(1:1) = ' '
               text_=.false.
               quote_ = ' '
               type_='null'
               goto 500
             end if
C            Here is the line we have read is part of the text field
10           continue
CDBG         print *,'processing as text field'
             jchar = lastch+1
             strg_(1:lastch) = buffer(1:lastch)
CDBG         print *, ' buffer before backup ', buffer(1:lastch)
             long_ = lastch
             text_ = .true.
             
             if (depth_.eq.0) then
               call getlin(flag)
               if(flag.eq.'fini') then
                 text_ = .false.
               else
                 if (buffer(1:1).eq.';') then
                   text_ = .false.
                   jchar = 2
                 end if
               end if
               if (text_) then
                 irecd = irecd-1
                 jchar=MAXBUF+1
               else
                 go to 500
               endif
CDBG         print *, ' buffer after backup ', buffer(1:lastch)
             end if
             goto 500
           else
C
C          Handle the case of a treble-quoted string
C          This is terminated by an unescaped quote
C
             if (tbxxtttq(jtloc)) then
               long_ = jtloc-jchar
               if (long_.eq.0) then
                 long_=0
                 strg_(1:1) = ' '
               else
                 strg_(1:long_)=buffer(jchar:jtloc)
               end if
               jchar = jtloc+3
               text_=.false.
               goto 500
             else
               jchar = lastch+1
               strg_(1:lastch) = buffer(1:lastch)
               long_ = lastch
               text_ = .true.
               goto 500
             end if
           end if
         end if
       end if
C
C      Now we are sure we are done with the multiline quoted cases
       quoted=.false.
       quote_=' '
       ttype_=' '
       if (depth_ .gt. 0) then
         ttype_ = typestack(1)
         type_ = typestack(depth_)
         state = statestack(depth_)
         index_= indexstack(depth_)
         go to 3000
       end if

C      We are not in a bracketed construct

       if(irecd.gt.0.and.
     *   jchar.le.1.and.lastch.gt.0) then
         jchar=1
         goto 140
       end if
100    jchar=jchar+1
       if(jchar.le.lastch)     goto 150
C
C....... Read a new line
C
110    call getlin(flag)
       type_='fini'
       dictype_=type_
       diccat_='(none)'
       dicname_=' '
CDBG   write(6,'(/5i5,a)')
CDBG *   irecd,jrecd,lrecd,nrecd,lastch, buffer(1:lastch)
       if(flag.eq.'fini')  goto 500
C
C....... Test if the new line is the start of a text sequence
C
140    if(buffer(1:1).ne.';') goto 150
       type_='text'
       quote_=';'
       jchar=lastch+1
       long_=lastch
       if (clipt_) then
         strg_(1:long_)=buffer(1:long_)
         strg_(1:1)=' '
       else
         if (long_.eq.1) then
           strg_(1:long_) = ' '
           long_ = 0
         else
           long_ = long_-1
           strg_(1:long_)=buffer(2:long_+1)
         endif
       endif
       text_ = .true.
       goto 500
C
C....... Process this character in the line
C
150    c=buffer(jchar:jchar)
       ip = jchar
       if(c.eq.' ')       goto 100
       if(c.eq.tab)       goto 100
       if(c.eq.'#')       goto 110
       if(c.eq.'''')      goto 300
       if(c.eq.'"')       goto 300
       if(c.eq.'('.and.rdprn_)       goto 350
       if(c.eq.'['.and.rdbkt_)       goto 360
       if(c.eq.'{'.and.rdbrc_)       goto 370
       if(c.ne.'_')       goto 200

       type_='name'
       goto 210
C
C....... Span blank delimited token; test if a number or a character
C
200    type_='numb'
       im=0
       quoted=.false.
       quote_=' '
       do 205 i=1,11
205    jj(i)=0
210    ip = jchar
       do 250 i=jchar,lastch
         ip = i
         if(buffer(i:i).eq.' ')       goto 400
         if(buffer(i:i).eq.tab)       goto 400
         if(type_.ne.'numb')          goto 250
         j=index(num,buffer(i:i))
         if(j.eq.0)                 type_='char'
         if(j.le.10) then
           im=im+1
           goto 250
         endif
         if(j.gt.13.and.im.eq.0) type_='char'
         jj(j-10)=jj(j-10)+1
250    continue
       i=lastch+1
       ip = i
       if(type_.ne.'numb') goto 400
       do 270 j=1,5
         if((jj(j).gt.1.and.j.gt.2) .or.
     *     jj(j).gt.2)             type_='char'
270    continue
       goto 400       
C
C....... Span quote delimited token; assume character
C
300    type_='char'
       quoted=.true.
       jchar=jchar+1
       if (rdtq_ .and. jchar+1 .le. lastch
     *   .and. buffer(jchar:jchar+1).eq.c//c) then
         quote_ = c//c//c
         jchar = jchar+2
         if (tbxxtttq(jtloc)) then
           text_=.false.
         else
           jtloc = lastch
           text_=.true.
           type_='text'
         endif
         long_ = jtloc-jchar
         if (long_.eq.0) then
           strg_=' '
         else
           strg_(1:long_)=buffer(jchar:jtloc)
         endif
         goto 500
       end if
       escaped = .false.
       ip = jchar
       do 320 i=jchar,lastch
         ip = i
         if (.not.escaped.and.rdrcqt_) then
           if (c.eq.slash) then
             escaped = .true.
             go to 320
           end if
         end if
         if (escaped) then
           escaped = .false.
           go to 320
         end if
         if(buffer(i:i).ne.c)             goto 320
         if(rdrcqt_)                      goto 400
         if(i+1.ge.lastch)                goto 400
         if (tbxxtsts(buffer(i+1:i+1)))   goto 400
320    continue
CDBG   write(6,'(a,4i5,a)')
CDBG *     '**** ',irecd,lastch,i,jchar,buffer(jchar:i)
       call tbxxwarn(' Quoted string not closed')
       i = lastch+1
       goto 400

C
C...... Here to start a bracketed construct
C
350    type_='tupl'
       go to 390
360    type_='list'
       go to 390
370    type_='tabl'
       if (.not. rdbkt_) type_='list'
390    depth_=1
       srecd=irecd
       ttype_=type_
       typestack(depth_) = type_
       brackstack(depth_) = c
       posbrkstk(depth_) = ip
       if (c.eq.':') brackstack(depth_) = ' '
       delimstack(depth_+1) = c
       posdlmstk(depth_+1) = ip
       recdlmstk(depth_+1) = irecd
       indexstack(depth_) = 1
       statestack(depth_) = 0
       state = 0
       go to 3100

C
C....... Here to process within a bracketed construct
C
3000   continue
CDBG   print *,' Processing in backeted construct '
       if(irecd.gt.0.and.
     *   jchar.le.1.and.lastch.gt.0) then
         jchar=1
         goto 3140
       endif

3100   jchar = jchar+1
       if(jchar.le.lastch)     goto 3150
C
C....... Read a new line
C
3110   call getlin(flag)
       if(flag.eq.'fini')  goto 3500

C
C....... Test if the new line is the start of a text sequence
C
3140   continue
CDBG   print *,buffer(jchar:lastch)
       if (buffer(1:1).ne.';') goto 3150
         type_='text'
         quote_=';'
         jchar=lastch+1
         long_=lastch
         if (long_ .gt. 1) then
         if (clipt_) then
           strg_(2:long_)=buffer(2:long_)
           strg_(1:1) = ' '
         else
           strg_(1:long_-1)=buffer(2:long_)
         long_ = long_-1
         endif
       else
         strg_(1:1)=' '
         long_ = 0
       endif
       state = 2
       statestack(depth_) = 2
       goto 500
C
C..... Process this character in the line
C      within a bracket construct
C
3150   continue
CDBG   print *,buffer(jchar:lastch)
       c=buffer(jchar:jchar)
       ip = jchar
       if(c.eq.' ')       goto 3100
       if(c.eq.tab)       goto 3100
       if(c.eq.'#')       goto 3110
       if(c.eq.'''')      goto 3300
       if(c.eq.'"')       goto 3300
       if(c.eq.'('.and.rdprn_)       goto 3350
       if(c.eq.'['.and.rdbkt_)       goto 3360
       if(c.eq.'{'.and.rdbrc_)       goto 3370
       if(c.eq.':'.and.rdcolon_)     goto 3380
       if(c.eq.',')       goto 3160
       if((c.eq.')' .and. brackstack(depth_).eq.'(') .or.
     *   (c.eq.'}' .and. brackstack(depth_).eq.'{') .or.
     *   (c.eq.']' .and. brackstack(depth_).eq.'[')) goto 3160
       if(c.eq.')' .and.
     *   ((brackstack(depth_).eq.'{').or.
     *   (brackstack(depth_).eq.'[')))
     *     call tbxxwarn(' Unbalanced ) treated as comma')
       if(c.eq.']' .and.
     *   ((brackstack(depth_).eq.'{').or.
     *   (brackstack(depth_).eq.'(')))
     *     call tbxxwarn(' Unbalanced ] treated as comma')
       if(c.eq.'}' .and.
     *   ((brackstack(depth_).eq.'(').or.
     *   (brackstack(depth_).eq.'[')))
     *     call tbxxwarn(' Unbalanced } treated as comma')
       go to 3200
C
C..... Process comma or close bracket found within bracketed construct
C
3160   continue
       if (state .eq. 2) then
         state = 0
         if (c.eq.',') then
           index_ = index_+1
           indexstack(depth_) = index_
           delimstack(depth_+1) = c
           posdlmstk(depth_+1) = ip
           recdlmstk(depth_+1) = irecd
           goto 3100
         endif
       endif
       depth_ = depth_-1
CDBG   print *,' decreasing depth ',depth_, recn_, 
CDBG *  buffer(jchar:lastch)

       type_='null'
       long_ = 1
       strg_(1:long_) = ' '
       if (depth_ .eq.0 ) goto 500
       state = 2
       statestack(depth_) = 2
       delimstack(depth_+1) = c
       posdlmstk(depth_+1) = ip
       recdlmstk(depth_+1) = irecd
       go to 500
C
C..... Process colon found within bracketed construct
C      treat as a comma if already started
C
3380   continue
       if (state .eq. 2) then
         state = 0
         index_ = index_+1
         indexstack(depth_) = index_
         delimstack(depth_) = c
         posdlmstk(depth_+1) = ip
         recdlmstk(depth_+1) = irecd
         goto 3100
       endif
       type_='null'
       long_ = 1
       strg_(1:long_) = ' '
       state = 2
       statestack(depth_) = 2
       delimstack(depth_) = c
       go to 500


C
C..... Span blank delimited token; test if a number or a character
C
3200   type_='numb'
       im=0
       innerdepth = depth_
       do 3205 i=1,11
3205     jj(i)=0
       ip = jchar
       do 3250 i=jchar,lastch
       ip = i
       if((buffer(i:i).eq.'('.and.rdprn_)
     *   .or.(buffer(i:i).eq.'{'.and.rdbrc_)
     *   .or.(buffer(i:i).eq.'['.and.rdbkt_)) then
         if (depth_ .ge. MAXDEPTH) then
           call tbxxerr(' Stack overflow, increase MAXDEPTH')
         end if
         depth_=depth_+1
CDBG   print *,' increasing depth ',depth_, recn_, 
CDBG *   buffer(jchar:lastch)
         typestack(depth_) = type_
         indexstack(depth_) = 1
         brackstack(depth_) = buffer(i:i)
       endif
       if((buffer(i:i).eq.')'.and.brackstack(depth_).eq.'(')
     *   .or.(buffer(i:i).eq.'}'.and.brackstack(depth_).eq.'{')
     *   .or.(buffer(i:i).eq.']'.and.brackstack(depth_).eq.'[')
     *   .or.(buffer(i:i).eq.' '.and.brackstack(depth_).eq.' ')) then
         if (depth_.gt.innerdepth) then
           depth_=innerdepth
           call tbxxwarn(
     * ' Failed to balance brackets in blank-delimited token')
CDBG   print *,' decreasing depth ',depth_, recn_, 
CDBG *  buffer(jchar:lastch)
         else
           go to 3395
         endif
       endif
       if(buffer(i:i).eq.' ')       goto 3400
       if(buffer(i:i).eq.tab)       goto 3400
       if(buffer(i:i).eq.',')       goto 3390
       if(buffer(i:i).eq.':'.and.rdcolon_) go to 3390
       if(buffer(i:i).eq.brackstack(depth_)) goto 3400
       if(type_.ne.'numb')          goto 3250
       j=index(num,buffer(i:i))
       if(j.eq.0)                 type_='char'
       if(j.le.10) then
         im=im+1
         goto 3250
       endif
       if(j.gt.13.and.im.eq.0) type_='char'
       jj(j-10)=jj(j-10)+1
3250   continue
       i=lastch+1
       ip = i
       if(type_.ne.'numb') goto 3400
       do 3270 j=1,5
         if((jj(j).gt.1.and.j.gt.2) .or.
     *     jj(j).gt.2)             type_='char'
3270   continue
       go to 3400

C
C..... Span '\'' or '\"' quote delimited token; assume character
C
3300   type_='char'
       quoted=.true.
       jchar=jchar+1

       if (rdtq_ .and. jchar+1 .le. lastch
     *   .and. buffer(jchar:jchar+1).eq.c//c) then
         quote_ = c//c//c
         jchar = jchar+2
         if (tbxxtttq(jtloc)) then
           text_=.false.
         else
           jtloc = lastch
           text_=.true.
           type_='text'
         endif
         long_ = jtloc-jchar
         if (long_.eq.0) then
           strg_=' '
         else
           strg_(1:long_)=buffer(jchar:jtloc)
         endif
         goto 500
       end if

       escaped = .false.
      
CDBG   print *,'Processing quoted string '
CDBG   print *,buffer(jchar:lastch)

       ip = jchar
       do 3320 i=jchar,lastch
         ip = i
         if (.not.escaped.and.rdrcqt_) then
           if (c.eq.slash) then
             escaped = .true.
             go to 3320
           end if
         end if
         if (escaped) then
           escaped = .false.
           go to 3320
         end if
CDBG     print *,'i,c,buffer(i:i)',i,c,buffer(i:i) 

         if(rdrcqt_)                      goto 3400
         if(buffer(i:i).ne.c)             goto 3320
         if(i+1.gt.lastch)                goto 3400
         if (tbxxtsts(buffer(i+1:i+1)))   goto 3400
3320    continue
       go to 3400
C.....  Span (-delimited tuple
3350    type_ = 'tupl'
       go to 3375
C.....  Span [-delimited list or array
3360    type_ = 'list'
       go to 3375
C.....  Span { delimited table
3370    type_ = 'tabl'
        if (.not. rdbkt_) type_='list'
3375    continue
       if (depth_ .ge. MAXDEPTH) then
         call tbxxerr(' Stack overflow, increase MAXDEPTH')
       end if
       depth_=depth_+1
CDBG   print *,' increasing depth ',depth_, recn_, 
CDBG *  buffer(jchar:lastch)
       typestack(depth_) = type_
       indexstack(depth_) = 1
       brackstack(depth_) = buffer(ip:ip)
       state = 0
       statestack(depth_) = state
       go to 3100


3390    if(depth_ .ne. innerdepth) then
         call tbxxwarn(' failed to close bracketed string')
         depth_ = innerdepth
        endif
        go to 3400


3395   ip = ip-1

3400   state = 2
       statestack(depth_) = state
       go to 400

3500   type_='fini'
       dictype_=type_
       diccat_='(none)'
       dicname_=' '
       if (depth_ .gt.0) then
         call tbxxwarn(
     *     ' File ended in unterminated bracketed construct')
         depth_ = 0
       endif
       go to 500



C
C..... Store the string for the getter
C
400    long_=0
       strg_=' '
         if(ip.gt.jchar) then
         long_=ip-jchar
         strg_(1:long_)=buffer(jchar:ip-1)
       endif
       jchar=ip
       quote_=' '
       if(quoted) then
         quote_=buffer(jchar:jchar)
         if (depth_.eq.0) jchar =jchar+1
       endif
       if(type_.ne.'char'.or.quoted.or.depth_.gt.0) goto 500
       if(strg_(1:5).eq.'data_') then
          type_='data'
          depth_=0
       end if
       if(strg_(1:5).eq.'loop_') then
         type_='loop'
         depth_=0
       end if
CDBG   if (strg_(1:max(1,long_)).eq.'?') print *,long_,strg_(1:1)
       if(long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
       if(long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
       if(strg_(1:5).eq.'save_') then
         type_='save'
         depth_=0
       end if
       if(long_.eq.7.and. strg_(1:7).eq.'global_') then
         type_='glob'
         depth_=0
       end if
C
500    continue
CDBG   print *,' leaving getstr with strg: ', strg_(1:long_)
CDBG   print *,' leaving getstr type, text_, quote_, depth_, jchar: '
CDBG   print *, type_,', ',text_,', ',quote_,', ',depth_,', ',jchar
       return
       end
C
C
C
C
C
C
C >>>>>> Convert a character string into a number and its esd
C
C                                          Q
C                                          D+
C                                          E-
C                                +         +
C           number string        -xxxx.xxxx-xxx(x)
C           component count CCNT 11111222223333444
C           (with at least 1 digit in the mantissa)
C
         subroutine ctonum
C
         integer   lastnb
         include  'ciftbx.sys'
         character test*26,c*1
         integer*4 m,nchar
         integer*4 ccnt,expn,msin,esin,ndec,ids,nmd
         integer*4 nms,ned,nef,nes
         double precision numb,sdev,ntemp,mant
         data test /'0123456789+.-()EDQedq :,]}'/
C
         numbtb=0.D0
         sdevtb=-1.D0
         numb=1.D0
         sdev=0.D0
         ccnt=0
         mant=0.D0
         expn=0.
         msin=+1
         esin=+1
         ndec=0
         ids=0
         nmd=0
         nms=0
         ned=0
         nef=0
         nes=0
         type_='char'
         posdec_=0
         esddig_=0
         if(long_.eq.1.and.
     *     index('0123456789',strg_(1:1)).eq.0) goto 500
         lzero_=.false.
         decp_=.false.
C
C....... Loop over the string and identify components
C
C        The scan works in phases
C          ccnt = 0   processing looking for first digit
C          ccnt = 1   processing before decimal point
C          ccnt = 2   processing after decimal point
C          ccnt = 3   processing exponent
C          ccnt = 4   processing standard deviation
C
         do 400 nchar=1,long_
C
         c=strg_(nchar:nchar)
         m=index(test,c)
         if(m.eq.0)     goto 500
         if(m.gt.10)    goto 300
C
C....... Process the digits
C
         if(ccnt.eq.0)  ccnt=1
         if(ccnt.eq.2)  ndec=ndec+1
         if(ccnt.gt.2)  goto 220
         ntemp=m-1
         if (ndec.eq.0) then
           mant=mant*10.D0+ntemp
         else
           mant=mant+ntemp/10.D0**(ndec)
         endif
         nmd=nmd+1
         if(ccnt.eq.1.and.mant.ne.0.D0) ids=ids+1
         goto 400
220      if(ccnt.gt.3)  goto 240
         expn=expn*10+m-1
         goto 400
240      esddig_=esddig_+1
         ntemp=m-1
         sdev=sdev*10.D0+ntemp
         sdevtb=1.D0
         goto 400
C
C....... Process the characters    . + - ( ) E D Q
C
300      if(c.ne.'.')  goto 320
         decp_=.true.
         if(nchar.gt.1.and.mant.eq.0.d0) then
           if(strg_(nchar-1:nchar-1).eq.'0') lzero_=.true.
         endif
         if(ccnt.gt.1) goto 500
         posdec_=nchar
         ccnt=2
         goto 400
C
320      if(nmd.eq.0.and.m.gt.13) goto 500
         if(c.ne.'(')  goto 340
         if(posdec_.eq.0) posdec_=nchar
         ccnt=4
         goto 400
C
340      if(posdec_.eq.0.and.ccnt.gt.0) posdec_=nchar
         if(c.eq.')' .or. c.eq.' ')  goto 400
         if(ccnt.eq.3 .and. ned.gt.0) goto 500
         if(m.gt.13) then
           if (nef.gt.0) goto 500
           nef = nef+1
           ccnt = 3
           esin = 1
         else
           if(ccnt.gt.0) then
             if (nes.gt.0) goto 500
             nes = nes+1
             ccnt = 3
             esin = 12-m
           else
             if (nms.gt.0) goto 500
             nms = nms+1
             ccnt=1
             msin=12-m
           endif
         endif
C
400      continue
C
         if(posdec_.eq.0) posdec_=lastnb(strg_(1:long_))+1
C
C....... String parsed; construct the numbers
C
         expn=expn*esin
         if(expn+ids.gt.-minexp) then
           call tbxxwarn(' Exponent overflow in numeric input')
           expn=-minexp-ids
         endif
         if(expn.lt.minexp) then
           call tbxxwarn(' Exponent underflow in numeric input')
           expn=minexp
         endif
         if(expn-ndec.lt.0) numb=1./10.D0**abs(expn-ndec)
         if(expn-ndec.gt.0) numb=10.D0**(expn-ndec)
         if(sdevtb.gt.0.0) sdevtb=numb*sdev
         numb=1.D0
         if(expn.lt.0) numb=1./10.D0**abs(expn)
         if(expn.gt.0) numb=10.D0**(expn)
         ntemp=msin
         numbtb=numb*mant*ntemp
         type_='numb'
C
500      return
         end
C
C
C
C
C
C
C >>>>>> Read a new line from the direct access file
C
         subroutine getlin(flag)
C
         include  'ciftbx.sys'
         character flag*4
         integer   kpp,lpp,mpp,npp,ir
         integer   tbxxrld
         integer   lip,mp,kip,ip,mip,mis
         integer   icpos,itpos,ixpos,ixtpos
C
         irecd=irecd+1
         jchar=1
         kpp = 0
         if(irecd.eq.jrecd.and.
     *     irecd.gt.recbeg_.and.
     *     irecd.le.recend_)  goto 200
         if(irecd.le.min(lrecd,recend_))  goto 100
         irecd=min(lrecd,recend_)+1
         buffer(1:1)=' '
         lastch=0
         jchar=MAXBUF+1
         jrecd=-1
         flag='fini'
         goto 200
100      continue
         lpp=-1
         mpp=-1
         npp=kpp
         call tbxxflin(irecd,lip,kpp,mp,kip,ip,mip,mis)
         if (lip.eq.0) then
           buffer(1:1) = ' '
           lastch = 1
           go to 130
         endif
         do ir = 1,NUMPAGE
           if(iabs(mppoint(ir)).eq.kpp) then
             lpp = ir
             goto 120
           endif
           if(mppoint(ir).eq.0) then
             lpp=ir
           else
             if(iabs(iabs(mppoint(ir))-kpp)
     *         .gt.iabs(npp-kpp)) then
               mpp=ir
               npp=iabs(mppoint(ir))
             endif
           endif
         enddo
C
C        failed to find page as resident
C        remove a target page
C
         if(lpp.eq.-1)lpp=mpp
         if(lpp.eq.-1)lpp=1
         if (mppoint(lpp).lt.0) then
           write(dirdev,'(a)',rec=-mppoint(lpp)) pagebuf(lpp)
         endif
         mppoint(lpp)=kpp
         read(dirdev,'(a)',rec=kpp) pagebuf(lpp)
120      lastch = tbxxrld(buffer,pagebuf(lpp)(mp:NUMCPP), .false.)
130      recn_=irecd
         jrecd=irecd
         flag=' '
         if (lastch.gt.0 .and. tabx_) then
           icpos=1
           itpos=1

140        ixpos=index(buffer(icpos:lastch),tab)
           ixtpos=ixpos+itpos-1
           if(ixpos.gt.0.and.ixtpos.le.MAXBUF) then
             ixtpos=((ixtpos+7)/8)*8
             if(ixpos.gt.1) then
               bufntb(itpos:ixtpos)=
     *           buffer(icpos:ixpos+icpos-2)
             else
             bufntb(itpos:ixtpos)=' '
             endif
             itpos=ixtpos+1
             icpos=ixpos+icpos
             goto 140
           else
             bufntb(itpos:min(MAXBUF,itpos+lastch-icpos))=
     *         buffer(icpos:lastch)
           endif
           buffer(1:min(MAXBUF,itpos+lastch-icpos))=
     *        bufntb(1:min(MAXBUF,itpos+lastch-icpos))
           lastch = min(MAXBUF,itpos+lastch-icpos)
         endif
200      return
         end
C
C
C
C
C
C
C >>>>>> Write error message and exit.
C
         subroutine tbxxerr(mess)
         character*(*) mess
         call tbxxcmsg('error',mess)
         stop
         end
C
C
C
C
C
C
C >>>>>> Write warning message and continue.
C
         subroutine tbxxwarn(mess)
         character*(*) mess
         call tbxxcmsg('warning',mess)
         return
         end
C
C
C
C
C
C
C >>>>>> Write a message to the error device
C
         subroutine tbxxcmsg(flag,mess)
C
         integer    lastnb
         include   'ciftbx.sys'
         character*(*) flag
         character*(*) mess
         character*(MAXBUF)  tline
         character*5   btype
         integer       ll,ls,ltry,ii,i
C
         btype = 'data_'
         if(save_) btype = 'save_'
         if(.not.glob_) then
         tline= ' ciftbx '//flag//': '
     *   //file_(1:longf_)//' '//btype
     *   //bloc_(1:max(1,lastnb(bloc_)))//' line:'
         else
         tline= ' ciftbx '//flag//': '
     *   //file_(1:longf_)//' global_'//' line:'
         endif
         ll = max(1,lastnb(tline))
         write(errdev,'(a,i7)')tline(1:ll),irecd
         ll=len(mess)
         ls=1
100      if(ll-ls.le.79) then
           write(errdev,'(1X,a)') mess(ls:ll)
           return
         else
           ltry = min(ll,ls+79)
           do ii = ls+1,ltry
           i = ltry-ii+ls+1
           if(mess(i:i).eq.' ') then
             write(errdev,'(1X,a)') mess(ls:i-1)
             ls=i+1
             if(ls.le.ll) go to 100
             return
           endif
           enddo
           write(errdev,'(1X,a)') mess(ls:ltry)
           ls=ltry+1
           if(ls.le.ll) go to 100
           return
         endif
         end
C
C
C
C
C >>>>>> Create a named file.
C
         function pfile_(fname)
C
         logical   pfile_
         include   'ciftbx.sys'
         logical   test
         integer   lfname
         integer   i
         character fname*(*)
C
C....... Test if a file by this name is already open.
C
         if(pfilef.eq.'yes') call close_
         pfilef='no '
         file_(1:longf_) = ' '
         lfname = len(fname)
         file_(1:lfname)=fname
         do 120 i=1,lfname
         if(file_(i:i).eq.' ') goto 140
120      continue
         i = lfname+1
140      if (i.gt.1) then
           inquire(file=file_(1:i-1),exist=test)
           pfile_=.false.
           longf_ = i-1
           if(test)            goto 200
         else
           file_ = ' '
           pfile_ = .true.
           longf_ = 1
         endif
C
C....... Open up a new CIF
C
         if (file_(1:1) .ne. ' ')  then
         open(unit=outdev,file=file_(1:longf_),status='NEW',
     *                    access='SEQUENTIAL',
     *                    form='FORMATTED')
         precn_=0
         endif
         pfile_=.true.
         pfilef='yes'
         nbloc=0
         pchar=1+lprefx
         pcharl=0
         obuf=prefx
         obuf(pchar:MAXBUF)=' '
200      ploopn = 0
         ploopc = 0
         ploopf = 'no '
         ptextf = 'no '
         pdepth_ = 0
         pdelimstack(1) = ' '
         pposdlmstk(1) = 0
         plcat = ' '
         pdblok = ' '
         plhead(1) = ' '
         if (xmlout_) then
           call tbxxpstr('<?xml version="1.0"?>')
         endif
         return
         end


C
C
C
C
C
C <<<<<< Substitute item in data block XML translation
C
         function tbxxxsub(oblok,xstring)
         include    'ciftbx.sys'
         character   oblok*(*)
         character   xstring*(*)
         character   tbxxxsub*(MAXBUF)
         integer     ii, jj, kk
         integer     lastnb
         jj = 1
         tbxxxsub = ' '
         do ii = 1,lastnb(xstring)
           if(xstring(ii:ii).ne.'%') then
             tbxxxsub(jj:jj) = xstring(ii:ii)
             jj = jj+1
           else
             do kk = 1,lastnb(oblok)
               tbxxxsub(jj:jj) = oblok(kk:kk)
               jj = jj+1
             enddo
           endif
         enddo
         return
         end

C
C
C
C
C
C >>>>>> Store a data block command in the CIF
C        Call with blank name to close current block only
C
         function pdata_(name)
C
         logical   pdata_
         SAVE
         include  'ciftbx.sys'
         character name*(*),temp*(MAXBUF)
         character dbloc(100)*(NUMCHAR)
         character tbxxxsub*(MAXBUF)
         integer   i
         integer   lastnb
C
         pdata_=.true.
         if(ptextf.eq.'yes') call tbxxeot
         if(pdepth_ .gt.0)   call tbxxebkt
         if(ploopn.ne.0)     call tbxxelp

         if(psaveo) then
           pchar=-1
           if(pposval_.ne.0) then
             pchar=lprefx+1
             call tbxxpstr(' ')
             pchar=lprefx+pposval_
             pposval_=0
           endif
           if (xmlout_) then
             call tbxxpxct('save_',' ')
           else
             call tbxxpstr('save_')
           endif
           psaveo=.false.
         endif
         if (pdblok(1:1).ne.' ') then
           if (xmlout_) then
             if (xmdata.eq.0) then
               call tbxxpxct(pdblok,' ')
             else
               call tbxxpxct(tbxxxsub(pdblok,xmlate(xmdata)),' ')
             endif
           endif
           pdblok=' '
         endif
         if(globo_) then
           pchar=-1
           temp='global_'
           pdblok='global_'
           psaveo=.false.
           goto 135
         endif
C
C....... Check for duplicate data name
C
         temp=name
         if(temp.eq.' ')        goto 200
         if(saveo_)             goto 130
         pdata_=.false.
         do 110 i=1,nbloc
         if(temp.eq.dbloc(i))   goto 130
110      continue
         pdata_ = .true.
         goto 125
C
C....... Save block name and put data_ statement
C
125      nbloc=nbloc+1
         if(nbloc.le.100) dbloc(nbloc)=temp(1:min(NUMCHAR,MAXBUF))
         pdblok = temp(1:min(NUMCHAR,MAXBUF))
130      pchar=-1
         temp='data_'//name
         if(saveo_) temp='save_'//name
         if(globo_) temp='global_'
         psaveo=saveo_
135      if(pposnam_.gt.0) then
           pchar=lprefx+1
           call tbxxpstr(' ')
           pchar=lprefx+pposnam_
           pposnam_=0
         endif
         if (xmlout_) then
           if (globo_) then
             call tbxxpxot('global_',' ')
           else
             if (xmdata.eq.0) then
               call tbxxpxot(pdblok,' ')
             else
               call tbxxpxot(tbxxxsub(pdblok,xmlate(xmdata)),' ')
             endif
             if (saveo_) then
               call tbxxpxot('save_',' ')
             endif
           endif
         else
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         pchar=lprefx
         plcat = ' '
         ploopn = 0
C
200      return
         end

C
C
C
C
C
C
C >>>>>> Process a name to extract the category and item
C
         subroutine tbxxgcat(name,type,flag,tflag,mycat,myxcat,
     *     item,xitem,nroot)
C
         character  name*(*),mycat*(*),item*(*),nroot*(*),type*4
         character  myxcat*(*),xitem*(*)
         include   'ciftbx.sys'
         character  xxxtemp*(NUMCHAR)
         logical    flag,tflag
         integer    lastnb,kpl,npl
         character  str1*(NUMCHAR), str2*(NUMCHAR)
         integer    kdc, lmycat


         item  = name
         xitem = ' '
         nroot = name
         mycat = ' '
         myxcat = ' '
         flag = .true.
         tflag = .true.
         if(vcheck.eq.'yes') then
           kdc = 0
           call tbxxdck(name,type,flag,tflag)
           if (xdchk.ne.0) then
             kdc = dcindex(xdchk)
             if (xmindex(xdchk).ne.0) xitem = xmlate(xmindex(xdchk))
           endif
           if (aliaso_.and.xdchk.ne.0) then
             if (aroot(xdchk).ne.0) then
               nroot = dictag(aroot(xdchk))
               kdc = dcindex(aroot(xdchk))
             endif
           endif
           if (kdc.ne.0) then
             mycat = dcname(kdc)
             myxcat = ' '
             if (xmcind(kdc).ne.0) myxcat = xmlate(xmcind(kdc))
           endif
         else
           call tbxxcat(name,mycat,lmycat)
         endif
         kpl = lastnb(mycat)
         npl = lastnb(name)
         call tbxxnlc(str1, mycat)
         call tbxxnlc(str2, name)
         if (mycat(1:1).ne.' ' .and. name(1:1).eq.'_') then
           if(str1(1:kpl).eq.str2(2:kpl+1) .and. npl .gt. kpl+2 .and.
     *       (name(kpl+2:kpl+2).eq.'.' .or.
     *       name(kpl+2:kpl+2).eq.'_') ) then
             item = name(kpl+3:npl)
           else
             item = name(2:npl)
           endif
         else
           if (mycat(1:1).eq.' ' .and. plcat(1:1).ne.' '
     *       .and. name(1:1).eq.'_') then
             call tbxxnlc(str1, plcat)
             kpl = lastnb(plcat)
             if(str1(1:kpl).eq.str2(2:kpl+1) .and. npl .gt. kpl+2 .and.
     *         (name(kpl+2:kpl+2).eq.'.' .or.
     *         name(kpl+2:kpl+2).eq.'_') ) then
               mycat = plcat
               item = name(kpl+3:npl)
             else
               item = name(2:npl)
             endif
           else
             item = name
             if (item(1:1).eq.'_') item = name(2:npl)
           endif
         endif
         if (xmlong_) then
           item = name
           if (item(1:1).eq.'_') item = name(2:npl)
         endif
         call tbxxnupc(xxxtemp,mycat)
         mycat = xxxtemp
         return
         end

C
C
C
C
C
C
C >>>>>> Put a number into the CIF, perhaps with an esd appended
C
         function pnumb_(name,numb,sdev)
C
         logical    pnumb_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
         character  mycat*(NUMCHAR),item*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         real       numb,sdev
         double precision dnumb,dsdev,dprec
         integer    kmn
         integer    lastnb
C
         pnumb_=.true.
         flag  =.true.
         tflag =.true.
         temp=name
         if(ptextf.eq.'yes') call tbxxeot
         if (pdepth_ .gt.0.and.name(1:1).ne.char(0))
     *     call tbxxebkt
C
         if(name(1:1).eq.' '.or.name(1:1).eq.char(0))
     *     goto 110
         call tbxxgcat(name,'numb',flag,tflag,mycat,myxcat,
     *     item,xitem,temp)
         pnumb_=flag
         if(ploopn.ne.0)        call tbxxelp
         if (xmlout_) then
           if (plcat(1:1).ne.' ' .and. plcat.ne.mycat) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             call tbxxpxct(plcat,plxcat)
             plcat = ' '
             plxcat = ' '
           endif
         endif
         pchar=-1
         if(pposnam_.ne.0)pchar=pposnam_+lprefx
         if (xmlout_) then
           if ((plhead(1)(1:1).ne.' '.and.plhead(1).ne.item)
     *       .or. plcat.ne.mycat) then
             call tbxxpxct (plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = mycat
               plxcat = myxcat
               call tbxxpxot(plcat,plxcat)
             endif
             call tbxxpxot (item,xitem)
           else
             if(plhead(1)(1:1).eq.' ') call tbxxpxot (item,xitem)
           endif
           plhead(1) = item
           plxhead(1) = xitem
         else
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         go to 120
C
110      if (xmlout_) then
           if (ploopn.gt.0) then
             kmn = mod(ploopc,ploopn)+2
             if (ploopn.gt.1.or.ploopf.eq.'yes') then
               call tbxxpxot(plhead(kmn),plxhead(kmn))
             endif
           endif
         endif
C
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         dprec=decprc
         dnumb=numb
         dsdev=sdev
         call tbxxpnum(dnumb,dsdev,dprec)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not in dictionary -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not in dictionary')
             call tbxxpstr(char(0))
           endif
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not correct type -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not correct type')
             call tbxxpstr(char(0))
           endif
         endif
         if (xmlout_) then
           if (ploopn.gt.1 .and.ploopc.gt.0) then
             call tbxxpxct(plhead(ploopc+1),plxhead(ploopc+1))
           endif
         endif
C
         pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         pesddig_=0
         return
         end
C
C
C
C
C
C
C >>>>>> Put a double precision number into the CIF, perhaps
C        with an esd appended
C
         function pnumd_(name,numb,sdev)
C
         logical    pnumd_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
         character  mycat*(NUMCHAR),item*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         double precision numb,sdev
         integer    kmn
         integer    lastnb
C
         pnumd_=.true.
         flag  =.true.
         tflag =.true.
         temp=name
         if(ptextf.eq.'yes') call tbxxeot
         if (pdepth_ .gt.0.and.name(1:1).ne.char(0))
     *      call tbxxebkt
C
         if(name(1:1).eq.' '.or.name(1:1).eq.char(0))
     *      goto 110
         call tbxxgcat(name,'numb',flag,tflag,mycat,myxcat,
     *     item,xitem,temp)
         pnumd_=flag
         if(ploopn.ne.0)        call tbxxelp
         if (xmlout_) then
           if (plcat(1:1).ne.' ' .and. plcat.ne.mycat) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             call tbxxpxct(plcat,plxcat)
             plcat = ' '
             plxcat = ' '
           endif
         endif
         pchar=-1
         if(pposnam_.ne.0)pchar=pposnam_+lprefx
         if (xmlout_) then
           if ((plhead(1)(1:1).ne.' '.and.plhead(1).ne.item)
     *       .or. plcat.ne.mycat) then
             call tbxxpxct (plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = mycat
               plxcat = myxcat
               call tbxxpxot(plcat,myxcat)
             endif
             call tbxxpxot (item,xitem)
           else
             if(plhead(1)(1:1).eq.' ') call tbxxpxot (item,xitem)
           endif
           plhead(1) = item
           plxhead(1) = xitem
         else
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         go to 120
C
110      if (xmlout_) then
           if (ploopn.gt.0) then
             kmn = mod(ploopc,ploopn)+2
             if (ploopn.gt.1.or.ploopf.eq.'yes') then
               call tbxxpxot(plhead(kmn),plxhead(kmn))
             endif
           endif
         endif
C
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         call tbxxpnum(numb,sdev,dpprc)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not in dictionary -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not in dictionary')
             call tbxxpstr(char(0))
           endif
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not correct type -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not correct type')
             call tbxxpstr(char(0))
           endif
         endif
         if (xmlout_) then
           if (ploopn.gt.1 .and.ploopc.gt.0) then
             call tbxxpxct(plhead(ploopc+1),plxhead(ploopc+1))
           endif
         endif
C
         pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         pesddig_=0
         return
         end
C
C
C
C
C
C
C >>>>>> Put a character string into the CIF.
C
         function pchar_(name,string)
C
         logical    pchar_
         include   'ciftbx.sys'
         logical    flag,tflag
         logical    pdelim_
         character  name*(*),temp*(NUMCHAR),string*(*)
         character  mycat*(NUMCHAR),item*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         character  line*(MAXBUF),strg*(MAXBUF)
         character*3 tsq,tdq,pqt
         integer    i, j, kfold, pql
         integer    lstring
         integer    lastnb
         integer    kmn, ic
         character*1 slash
C
         slash = rsolidus(1:1)
         pchar_=.true.
         flag  =.true.
         tflag =.true.
         tsq = ''''''''
         tdq = '"""'
         temp  =name
         lstring = lastnb(string)
         if (lstring .gt. MAXBUF) then
           call tbxxwarn(
     *       'Output CIF line longer than MAXBUF, truncated')
           lstring = MAXBUF
         endif
         pqt = pquote_
         pql = lastnb(pqt)
         if(ptextf.eq.'yes') call tbxxeot
         if(pdepth_ .gt.0.and.name(1:1).ne.char(0))
     *     call tbxxebkt
C
         if(name(1:1).eq.' '.or.name(1:1).eq.char(0))
     *     goto 110
         call tbxxgcat(name,'char',flag,tflag,mycat,myxcat,
     *     item,xitem,temp)
         pchar_=flag
         if(ploopn.ne.0)        call tbxxelp
         if (xmlout_) then
           if (plcat(1:1).ne.' ' .and. plcat.ne.mycat) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             call tbxxpxct(plcat,plxcat)
             plcat = ' '
             plxcat = ' '
           endif
         endif
         pchar=-1
         if(pposnam_.gt.0) pchar=posnam_+lprefx
         if (xmlout_) then
           if ((plhead(1)(1:1).ne.' '.and.plhead(1).ne.item)
     *       .or. plcat.ne.mycat) then
             call tbxxpxct (plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = mycat
               plxcat = myxcat
               call tbxxpxot(plcat,plxcat)
             endif
             call tbxxpxot (item,xitem)
           else
             if(plhead(1)(1:1).eq.' ') call tbxxpxot (item,xitem)
           endif
           plhead(1) = item
           plxhead(1) = xitem
         else
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         go to 120
C
110      if (xmlout_) then
           if (ploopn.gt.0) then
             kmn = mod(ploopc,ploopn)+2
             if (ploopn.gt.1.or.ploopf.eq.'yes') then
               call tbxxpxot(plhead(kmn),plxhead(kmn))
             endif
           endif
         endif
C
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         i=1
         if (string(1:1).eq.char(0)) go to 210
         if (xmlout_) then
           do ic = 1,lstring
            if ( string(ic:ic).eq.'&'
     *        .or. string(ic:ic).eq.'<'
     *        .or. string(ic:ic).eq.'>' ) then
              if(i.lt.MAXBUF) then
                line(i:i) = '&'
              endif
              if (i.lt.MAXBUF) then
                if( string(ic:ic).eq.'&' ) then
                  line(i:MAXBUF)='amp;'
                  i = i+4
                endif
                if( string(ic:ic).eq.'<' ) then
                  line(i:MAXBUF)='lt;'
                  i = i+3
                endif
                if( string(ic:ic).eq.'>' ) then
                  line(i:MAXBUF)='gt;'
                  i = i+3
                endif
              endif
              if (i.gt.MAXBUF+1) then
                i = MAXBUF+1
              endif
            else
              if(i.lt.MAXBUF) then
                line(i:i) = string(ic:ic)
                i = i+1
              endif
            endif
           enddo
           if (i.gt.1) i = i-1
           if (i.lt.MAXBUF) line(i+1:MAXBUF) = ' '
         else
           line=string
           i = lstring
         endif
         if(pposval_.ne.0.and.pposend_.ge.pposval_)
     *      i=max(i,pposend_-pposval_+1)
         if(pfold_ .ne. 0 .and. lstring .gt. min(pfold_,line_) )
     *      go to 290
         if (i .gt. MAXBUF) then
           call tbxxwarn(
     *       'Output CIF line longer than MAXBUF, truncated')
           i = MAXBUF
         endif
         if(pquote_.ne.' ')   go to 150
         do 140 j=i,1,-1
         if(line(j:j).eq.' ') go to 150
140      continue
         if((line(1:1).eq.'_'
     *     .or. line(i:i).eq.'_'
     *     .or. line(1:1).eq.''''
     *     .or. line(1:1).eq.'"'
     *     .or. line(1:1).eq.';'
     *     .or. line(1:1).eq.'('
     *     .or. line(1:1).eq.'['
     *     .or. line(1:1).eq.'(')
     *     .and.line(1:i).ne.'''.'''
     *     .and.line(1:i).ne.'''?'''
     *     .and.line(1:i).ne.'"."'
     *     .and.line(1:i).ne.'"?"') go to 150
         strg=line(1:i)
         goto 200
150      if(pqt.eq.';'
     *     .or. pqt.eq. tsq
     *     .or. pqt.eq. tdq
     *     .or. pqt.eq. '('
     *     .or. pqt.eq. '{'
     *     .or. pqt.eq. '[')       go to 190
         if(line(1:i).eq.' '.and.nblanko_) then
           strg = '.'
           i = 1
           if(pposval_.ne.0) then
             pchar=pposval_+lprefx
           endif
           call tbxxpstr(strg(1:i))
           go to 210
         endif
         if(pqt.eq.'"' .or. pqt.eq.'"""')       go to 170
         do 160 j=1,i-1
         if(line(j:j).eq.''''.and.
     *     (line(j+1:j+1).eq.' '.or.line(j+1:j+1).eq.tab))
     *     goto 170
160      continue
         strg=''''//line(1:i)//''''
         i=i+2
         pqt = ''''
         pql = 1
         goto 200
170      do 180 j=1,i-1
         if(line(j:j).eq.'"'.and.
     *     (line(j+1:j+1).eq.' '.or.line(j+1:j+1).eq.tab))
     *     goto 190
180      continue
         strg='"'//line(1:i)//'"'
         i=i+2
         pqt = '"'
         pql = 1
         if(pfold_ .gt. 1 .and. i .gt. min(pfold_,line_) ) go to 290
         goto 200
190      pchar=-1
         if (xmlout_) then
           if (pqt.eq.';') then
           strg = '<![CDATA['//line(1:i)
           i = i+9
         else
             strg = pqt(1:pql)//'<![CDATA['//line(1:i)
             i = i+9+pql
           endif
         else
           if (pclipt_ .and. pqt.eq.';' .and. line(1:1).eq.' ') then
             strg =pqt(1:pql)//line(2:i)
             i = i-1+pql
           else
             strg =pqt(1:pql)//line(1:i)
             i = i+pql
           endif
         endif
         ptextf='yes'
         call tbxxpstr(strg(1:i))
         pchar=-1
         ptextf='no '
         if (pqt.eq.'('
     *     .or.pqt.eq.'{'
     *     .or.pqt.eq.'[') then
           if (pqt.eq.'(') call tbxxpstr(')')
           if (pqt.eq.'{') call tbxxpstr('}')
           if (pqt.eq.'[') call tbxxpstr(']')
         else
           call tbxxpstr(pqt(1:pql))
         endif
         if (xmlout_) then
           call tbxxpstr(']]>')
         endif
         pchar=lprefx
         call tbxxpstr(' ')
         strg =
     *    ' Converted pchar_ output to text for: '//string(1:lstring)
         call tbxxwarn(strg)
         goto 210
C
200      if(pposval_.ne.0) then
           pchar=pposval_+lprefx
           if(pqt.ne.' ') pchar=pchar-pql
         endif
         if (pdepth_.gt.0) then
           if(pstatestack(pdepth_).eq.2)
     *     tbxxrslt = pdelim_(',',.false.,0)
         end if
         call tbxxpstr(strg(1:i))
         if (pdepth_.gt.0) pstatestack(pdepth_) = 2
210      if(.not.flag) then
           pchar = pcharl+4
           if (.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not in dictionary -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not in dictionary')
             call tbxxpstr(char(0))
           endif
         endif
         if((.not.tflag).and.line(1:i).ne.'.'.and.
     *     line(1:i).ne.'?'.and.pqt.eq.' ') then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not correct type -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not correct type')
             call tbxxpstr(char(0))
           endif
         endif
         if (xmlout_) then
           if (ploopn.gt.1 .and. ploopc.gt.0) then
             call tbxxpxct(plhead(ploopc+1),plxhead(ploopc+1))
           endif
         endif


         pposval_=0
         pposdec_=0
         pposnam_=0
         pposend_=0
         pquote_=' '
         return

C
C        fold a string to min(pfold_,line_)
C
290      pchar=-1
         pqt = ';'
         pql = 1
         if (xmlout_) then
           call tbxxpstr('<![CDATA['//slash)
         else
           call tbxxpstr(';'//slash)
         endif
         kfold = min(pfold_,line_)
         call tbxxpfs(string(1:lstring),' ',kfold)
         pchar=-1
         ptextf='no '
         if (xmlout_) then
           call tbxxpstr(']]>')
         else
           call tbxxpstr(';')
         endif
         pchar=lprefx
         call tbxxpstr(' ')
         if (pdepth_.gt.0) pstatestack(pdepth_) = 2
         goto 210

         end
C
C
C
C
C
C >>>>>> Put a comment in the output CIF
C
         function pcmnt_(string)
C
         logical    pcmnt_
         include   'ciftbx.sys'
         character  string*(*), temp*(MAXBUF)
         character*3 pqt
         character*1 slash
         integer lstring, kfold, pql, ic, ik
         integer lastnb
C
         slash = rsolidus(1:1)
         lstring = lastnb(string)
         pqt = pquote_
         pql = lastnb(pquote_)
         kfold = min(pfold_,line_)
         if(ptextf.eq.'yes') call tbxxeot
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         if(string.eq.' '.or.
     *     (string.eq.char(0)) .or.
     *     (string.eq.tab.and.(.not.ptabx_))) then
           if(string.eq.' ') pchar=-1
           if (pqt.eq.'#') then
             temp(1:1+lstring) = pqt(1:pql)//string(1:lstring)
             call tbxxpstr(temp(1:1+lstring))
           else
             call tbxxpstr(string)
           endif
           if(string.eq.' ') call tbxxpstr(char(0))
         else
           if ((kfold .ne. 0) .and.
     *        ((xmlout_ .and. (max(pchar,1)+8+lstring.gt.kfold))
     *        .or.((.not.xmlout_) .and.
     *        ((max(pchar,1)+lstring).gt.kfold)))) then
           if (xmlout_) then
             call tbxxpstr('<!-- '//slash)
             call tbxxpfs(string(1:lstring),' ',kfold)
           else
             call tbxxpstr('#'//slash)
             call tbxxpfs(string(1:lstring),"#",kfold)
           endif
           call tbxxpstr(char(0))
           else
           if (xmlout_) then
             temp(1:5) = '<!-- '
             ik = 6
             do ic = 1,lastnb(string)
               if (string(ic:ic).eq.'-'.and.
     *           temp(ik-1:ik-1).eq.'-'.and.
     *           ik.lt.MAXBUF-4) then
                 temp(ik:ik) = slash
                 ik = ik+1
               endif
               if (ik.lt.MAXBUF-4) then
                 temp(ik:ik) = string(ic:ic)
                 ik = ik+1
               endif
             enddo
             temp(ik:ik+3) = ' -->'
             ik = ik+4
             if (ik.lt.MAXBUF) temp(ik:MAXBUF) = ' '
           else
             temp='#'//string
           endif
           call tbxxpstr(temp(1:lastnb(temp)))
           call tbxxpstr(char(0))
           endif
         endif
         pcmnt_=.true.
         pposnam_=0
         if(string.ne.tab)pchar=lprefx+1
         return
         end
C
C
C
C
C
C
C
C >>>>>> Put a text sequence into the CIF.
C
         function ptext_(name,string)
C
         logical    ptext_
         integer    lastnb
         include   'ciftbx.sys'
         logical    flag,tflag
         integer    ll
         character  name*(*),temp*(NUMCHAR),string*(*),store*(NUMCHAR)
         character  mycat*(NUMCHAR),item*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         character  temp2*(MAXBUF)
         character  slash*1
         character*3 pqt
         integer    pql
         integer    kmn
         integer    kfold
         data       store/'                                '/
C
CDBG     print *,' ptext_, pclipt_ ', pclipt_
         ptext_=.true.
         flag  =.true.
         tflag =.true.
         slash = rsolidus(1:1)
         pqt = pquote_
         if (pqt .eq. ' ') pqt = ';'
         pql = lastnb(pqt)
         ll=lastnb(string)
         temp=name
         if(ptextf.eq.'no ')    goto 100
         if(temp.eq.store)      goto 150
         call tbxxeot
         if (pdepth_ .gt.0)  call tbxxebkt
C
100      if(name(1:1).ne.' ')   goto 110
         if(ptextf.eq.'yes')    goto 150
         goto 120
C
110      if(ploopn.ne.0)        call tbxxelp
         call tbxxgcat(name,'char',flag,tflag,mycat,myxcat,
     *     item,xitem,temp)
         ptext_=flag
         if (xmlout_) then
           if (plcat(1:1).ne.' ' .and. plcat.ne.mycat) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             call tbxxpxct(plcat,plxcat)
             plcat = ' '
             plxcat = ' '
           endif
         endif
         pchar=-1
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         if (xmlout_) then
           if ((plhead(1)(1:1).ne.' '.and.plhead(1).ne.item)
     *       .or. plcat.ne.mycat) then
             call tbxxpxct (plhead(1),plxhead(1))
             plhead(1) = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = mycat
               plxcat = myxcat
               call tbxxpxot(plcat,plxcat)
             endif
             call tbxxpxot (item,xitem)
           else
             if(plhead(1)(1:1).eq.' ') call tbxxpxot (item,xitem)
           endif
           plhead(1) = item
           plxhead(1) = xitem
         else
           call tbxxpstr(temp)
         endif
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not in dictionary -->')
             call tbxxpstr(char(0))
           else
             call tbxxpstr('#< not in dictionary')
             call tbxxpstr(char(0))
           endif
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           if (xmlout_) then
             call tbxxpstr('<!-- not correct type -->')
           else
             call tbxxpstr('#< not correct type')
           endif
         endif
         go to 130
C
120      if (xmlout_) then
           if (ploopn.gt.0) then
             kmn = mod(ploopc,ploopn)+2
             if (ploopn.gt.1.or.ploopf.eq.'yes') then
               call tbxxpxot(plhead(kmn),plxhead(kmn))
             endif
           endif
         endif
C
130      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         ptextf='yes'
         store=temp
         if (pfold_.eq.0) then
           pchar=-1
           if (xmlout_) then
             if (pclipt_ .and. string(1:1).eq.' '.and.ll.gt.1) then
               if (pqt.eq.';') then
             temp2 = '<![CDATA['//string(2:ll)
               else
                 temp2 = '<![CDATA['//pqt(1:pql)//string(2:ll)
               endif
             else
               if (pqt.eq.';') then
                 temp2 = '<![CDATA['//string(1:ll)
           else
                 temp2 = '<![CDATA['//pqt(1:pql)//string(1:ll)
               endif
             endif
           else
             if (pclipt_ .and. string(1:1).eq.' '.and.ll.gt.1) then
               temp2 = pqt(1:pql)//string(2:ll)
             else
               temp2 = pqt(1:pql)//string(1:ll)
             endif
           endif
           call tbxxpstr(temp2(1:lastnb(temp2)))
           pchar=-1
           return
         endif
         pchar=-1
         if (xmlout_) then
           if (pqt.eq.';') then
             call tbxxpstr('<![CDATA['//slash)
           else
             call tbxxpstr('pqt(1:pql)//<![CDATA['//slash)
           endif
         else
           call tbxxpstr(pqt(1:pql)//slash)
         endif
150      pchar=-1
         kfold = 0
         if (pfold_ .gt. 0 ) kfold = pfold_
         if (line_ .gt. 0 .and. pfold_ .gt. line_) then
            kfold = line_
         endif
         if (kfold .gt. 0) then
           call tbxxpfs(string,' ',kfold)
         else
           call tbxxpstr(string(1:max(1,ll)))
         endif
         continue
         pchar=-1
         pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         return
         end
C
C
C
C
C
C
C >>>>>> Put a folded string to the output CIF
C
         subroutine tbxxpfs(string,prefix,kfold)
C
         include   'ciftbx.sys'
         character *(*) string,prefix
         character *(MAXBUF) temp
         character *1 slash
         logical stabl
         integer kfold
         integer sploopn
         integer i, klow, khi, kpref, klen
         integer lastnb

         slash = rsolidus(1:1)
         sploopn = ploopn
         ploopn = -1
         stabl = tabl_
         tabl_ = .false.
         if (kfold .lt. 4) then
           call
     *     tbxxwarn(
     *       'Invalid attempt to fold output line, limit reset to 4')
           pfold_ = 4
           kfold = 4
         endif
         klen = lastnb(string)
         kpref = len(prefix)
         if (prefix.eq.' ') kpref=0
         klow = 1
 100     khi = klen
         if (khi.gt.klow+kfold-1-kpref) then
           khi = klow+kfold-1-kpref-1
           do i = khi,klow+1,-1
             if(string(i:i).eq.' ') then
               khi = i
               go to 120
             endif
           enddo
 120       if (kpref.gt.0) then
             temp(1:kpref+khi-klow+2) = prefix//string(klow:khi)//slash
           else
             temp(1:kpref+khi-klow+2) = string(klow:khi)//slash
           endif
           pchar = -1
           call tbxxpstr(temp(1:kpref+khi-klow+2))
           call tbxxpstr(char(0))
           klow = khi+1
           go to 100
         else
           if (string(khi:khi).eq.slash) then
             if (khi.lt.klow+kfold-1-kpref) then
                if (kpref.gt.0) then
                  temp(1:kpref+khi-klow+2) =
     *              prefix//string(klow:khi)//slash
                  pchar = -1
                  call tbxxpstr(temp(1:kpref+khi-klow+2))
                  call tbxxpstr(char(0))
                  pchar = -1
                  call tbxxpstr(prefix)
                else
                  temp(1:khi-klow+2) = string(klow:khi)//slash
                  pchar = -1
                  call tbxxpstr(temp(1:khi-klow+2))
                  pchar = -1
                  call tbxxpstr(' ')
                endif
                call tbxxpstr(char(0))
             else
                if (kpref.gt.0) then
                  temp(1:kpref+khi-klow+1) = prefix//string(klow:khi)
                  pchar = -1
                  call tbxxpstr(temp(1:kpref+khi-klow+1))
                  call tbxxpstr(char(0))
                  temp(1:kpref+2) = prefix//slash//slash
                  pchar = -1
                  call tbxxpstr(temp(1:kpref+2))
                  call tbxxpstr(char(0))
                  call tbxxpstr(prefix)
                else
                  pchar = -1
                  call tbxxpstr(string(klow:khi))
                  call tbxxpstr(char(0))
                  pchar = -1
                  call tbxxpstr(slash//slash)
                  call tbxxpstr(char(0))
                  call tbxxpstr(' ')
                endif
                call tbxxpstr(char(0))
                pchar = -1
             endif
           else
             pchar = -1
             if (kpref.gt.0) then
               temp(1:kpref+khi-klow+1)=prefix//string(klow:khi)
               call tbxxpstr(temp(1:kpref+khi-klow+1))
             else
               call tbxxpstr(string(klow:khi))
             endif
             call tbxxpstr(char(0))
           endif
         endif
         pchar = -1
         ploopn = sploopn
         tabl_ = stabl
         return
         end
C
C
C
C
C
C
C >>>>>> Put a delimiter symbol into the CIF.
C
         function pdelim_(delim,force,posdlm)
C
         logical    pdelim_
         character*(*) delim
         logical    force
         integer    posdlm
         
         include   'ciftbx.sys'
         pdelim_ = .true.
         if (ptextf.eq.'yes') call tbxxeot
         if (delim.eq.'('
     *     .or. delim.eq.'['
     *     .or. delim.eq.'{') then
           if (pdepth_.gt.0) then
             if (pstatestack(pdepth_).eq.2) then
               pposdlmstk(pdepth_+1) = pchar
               call tbxxpstr(',')
               pdelimstack(pdepth_+1) = ','
               pstatestack(pdepth_) = 1
             end if
           else
             if (ploopf.eq.'yes') ploopc=0
             ploopf='no '
             ploopc = ploopc+1
             if (ploopc.gt.ploopn) ploopc=ploopc-ploopn
           end if
           pdepth_ = pdepth_+1
           pbrackstack(pdepth_) = delim
           pdelimstack(pdepth_+1) = delim
           pstatestack(pdepth_) = 1
           pposbrkstk(pdepth_) = posdlm
           pposdlmstk(pdepth_+1) = posdlm 
           go to 100
         end if
         if (delim.eq.'}') then
           if (pdepth_.eq.0) then
             if (.not.force) pdelim_ = .false.
           else
             if (pbrackstack(pdepth_) .ne.'{'
     *         .and. .not.force) pdelim_ = .false.
           end if
         end if
         if (delim.eq.')') then
           if (pdepth_.eq.0) then
             if (.not.force) pdelim_ = .false.
           else
             if (pbrackstack(pdepth_) .ne.'('
     *         .and. .not.force) pdelim_ = .false.
           end if
         end if
         if (delim.eq.']') then
           if (pdepth_.eq.0) then
             if (.not.force) pdelim_ = .false.
           else
             if (pbrackstack(pdepth_) .ne.'['
     *         .and. .not.force) pdelim_ = .false.
           end if
         end if
         if (delim.eq.':' .and. pdelimstack(pdepth_+1).eq.':') then
            if (.not.force) pdelim_= .false.
         end if
         if (.not.pdelim_) return
         if (delim.eq.'}' .or. delim.eq.']' .or. delim.eq.')') then
           pdepth_ = pdepth_-1
CDBG       print *,' decreasing pdepth ',pdepth_, precn_
           if (pdepth_ .gt.0) pstatestack(pdepth_) = 2
           go to 100
         end if
         if (delim.eq.',' .or. delim.eq.':') then
           pstatestack(pdepth_) = 1
           pdelimstack(pdepth_+1) = delim
           pposdlmstk(pdepth_+1) = posdlm
         end if
100      continue 
         if (posdlm.ne.0) pchar = lprefx+posdlm
         call tbxxpstr(delim)
         return
         end
         
C
C
C
C
C
C
C >>>>>> Put a loop_ data name into the CIF.
C
         function ploop_(name)
C
         logical    ploop_
         integer    kpc,kpl,npl
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*), temp*(NUMCHAR), mycat*(NUMCHAR)
         character  myxcat*(XMLCHAR),xitem*(XMLCHAR)
         character  item*(NUMCHAR), str*(NUMCHAR)
         character  shead*(NUMCHAR),xshead*(XMLCHAR)
         integer    lastnb
C
         ploop_=.true.
         flag  =.true.
         if(ptextf.eq.'yes')    call tbxxeot
         if (pdepth_ .gt.0)     call tbxxebkt
         if(ploopn.ne.0. and. ploopf.ne.'yes'
     *     .and. name(1:1).eq.' ') then
           call tbxxelp
         endif
         temp = ' '
         mycat = ' '
         item = ' '
         shead = plhead(1)
         xshead = plxhead(1)(1:min(XMLCHAR,NUMCHAR))
         str = ' '
         if(name(1:1).eq.' ')   goto 100
C
         call tbxxgcat(name,'    ',flag,tflag,mycat,myxcat,
     *     item,xitem,str)
         ploop_ = flag
         if (ploopn.ne.0. and. ploopf.ne.'yes') then
           if (plcat.eq.mycat) then
             plcat = ' '
             call tbxxelp
             plcat = mycat
             plxcat = myxcat
           else
             call tbxxelp
           endif
         endif
         if (xmlout_) then
           if (plcat(1:1).ne.' '.and.ploopn.eq.0) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
             shead = ' '
             xshead = ' '
             if (plcat.ne.mycat) then
               call tbxxpxct(plcat,plxcat)
               plcat = ' '
               plxcat = ' '
             endif
           endif
         endif
         if(tabl_.and.pposnam_.eq.0) then
           temp='    '//str(1:NUMCHAR-4)
         else
           temp=str
         endif
         plhead(max(ploopn,0)+2) = item
         plxhead(max(ploopn,0)+2) = xitem
100      if(ploopn.ne.0)        goto 120
         ploopf='yes'
         pchar=-1
         if(pposval_.ne.0) then
           pchar=lprefx+1
C           call tbxxpstr(' ')
           pchar=pposval_+lprefx
         else
           if(pposnam_.ne.0) then
             pchar=lprefx+1
             call tbxxpstr(' ')
             pchar=pposnam_+lprefx+1
           endif
         endif
         if (xmlout_) then
           if (shead(1:1).ne.' ') then
             call tbxxpxct (shead,xshead)
           endif
         else
           call tbxxpstr('loop_')
         endif
         pchar=-1
         if(name(1:1).eq.' ') then
           ploopn=-1
           plhead(1) = ' '
           plxhead(1) = ' '
           return
         endif
120      if(ploopn.le.0) then
           if (xmlout_.and.plcat.ne.mycat) then
             call tbxxpxct(plcat,plxcat)
             plcat = mycat
             plxcat = myxcat
             call tbxxpxot(mycat,myxcat)
           endif
         else
           if(xmlout_ .and. plcat.ne.mycat) then
             kpl = lastnb(plcat)
             if(mycat(1:1).eq.' ') then
                mycat = '(none)'
                myxcat = '_NONE_ '
             endif
             npl = lastnb(mycat)
             kpc = pchar
             call tbxxpstr('<!--  mixed categories '//plcat(1:kpl)
     *                      //' and '//mycat(1:npl)//' -->')
             pchar = kpc
           endif
         endif
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         if (.not. xmlout_) then
           call tbxxpstr(temp(1:lastnb(temp)))
         endif
         if(flag)               goto 130
         if(.not.tabl_) pchar=lprefx+57
         if (xmlout_) then
           call tbxxpstr('<!-- '//temp(1:lastnb(temp))//
     *     ' not in dictionary -->')
           call tbxxpstr(char(0))
         else
           call tbxxpstr('#< not in dictionary')
           call tbxxpstr(char(0))
         endif
130      pchar=lprefx+1
         ploopn=max(ploopn,0)+1
         ploopc = 0
C
         return
         end
C
C
C
C
C
C >>>>>> Create or clear a prefix string
C        Any change in the length of the prefix string flushes
C        pending text, if any,  loops and partial output lines
C
         function prefx_(strg,lstrg)
C
         logical    prefx_
         include   'ciftbx.sys'
         character  strg*(*)
         integer    lstrg,mxline
C
         mxline=MAXBUF
         if(line_.gt.0) mxline=min(line_,MAXBUF)
         if(lstrg.ne.lprefx.and.pcharl.gt.0) then
           pchar=-1
           call tbxxpstr(' ')
         endif
         if (lstrg.le.0) then
           prefx=' '
           if(pchar.ge.lprefx+1)pchar=pchar-lprefx
           lprefx=0
         else
           if(lstrg.gt.mxline) then
             call tbxxwarn(' Prefix string truncated')
           endif
           prefx=strg
           if(pchar.ge.lprefx+1)pchar=pchar-lprefx+lstrg
           obuf(1:min(mxline,lstrg))=prefx
           lprefx=lstrg
           if(mxline-lprefx.lt.NUMCHAR) then
             call tbxxwarn(' Output prefix may force line overflow')
           endif
         endif
         prefx_=.true.
         return
         end
C
C
C
C
C
C
C >>>>>> Close the CIF
C
         subroutine close_
C
         include   'ciftbx.sys'
         character  tbxxxsub*(MAXBUF)
C
         if(ptextf.eq.'yes') call tbxxeot
         if (pdepth_ .gt.0)  call tbxxebkt
         if(ploopn.ne.0)     call tbxxelp
         if (xmlout_) then
           if (plhead(1)(1:1).ne.' ')
     *       call tbxxpxct(plhead(1),plxhead(1))
           if (plcat(1:1).ne.' ') call tbxxpxct(plcat,plxcat)
           if (pdblok(1:1).ne.' ')  then
             if (xmdata.eq.0) then
               call tbxxpxct(pdblok,' ')
             else
               call tbxxpxct(tbxxxsub(pdblok,xmlate(xmdata)),' ')
             endif
           endif
         endif
         pdblok = ' '
         plcat = ' '
         plxcat = ' '
         plhead(1) = ' '
         plxhead(1) = ' '
         if(pcharl.ge.lprefx+1) then
           pchar=-1
           call tbxxpstr(' ')
         endif
         if (file_(1:1) .ne. ' ') then
           file_(1:longf_) = ' '
           longf_ = 1
           close(outdev)
           precn_=0
         endif
         return
         end
C
C
C >>>>>  Clean out characters not valid for an XML name
C
C        An XML name may begin with a letter, '_' or ':'
C        and may contain letters, digits, '_', ':', '.' or '-'
C
C        Note that the full Unicode character set would also permit
C        combining characters and extender characters, but these
C        have no representation in a 128 character ASCII set
C
C
         function tbxxxcln(xstring,lstr)
         logical tbxxxcln
         character*(*) xstring
         integer lstr, ii, ix
         character*10 chkstr1
         character*28 chkstr2
         character*28 chkstr3
         character*1 c
         data chkstr1/'0123456789'/
         data chkstr2/'_:ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data chkstr3/'abcdefghijklmnopqrstuvwxyz.-'/
         tbxxxcln = .true.
         do ii = 1,lstr
           c = xstring(ii:ii)
           if (c.eq.' ') return
           ix = index(chkstr2,c)
           if (ix.eq.0) then
             ix = index(chkstr3,c)
             if (ix.eq.0.and.ii.gt.1) then
               ix = index(chkstr1,c)
             endif
             if(ix.eq.0.or.(ix.gt.26.and.ii.eq.1)) then
               xstring(ii:ii) = '_'
               tbxxxcln = .false.
             endif
           endif
         enddo
         return
         end




C
C
C >>>>>> Put out the given string as an xml open tag
C
C        Note that the string may have embedded blanks and
C        parameters.  The second argument is an optional
C        translation to be used if non-blank.
C
         subroutine tbxxpxot(string,xstring)
C
         integer lastnb
         include 'ciftbx.sys'
         character sbuf*(MAXBUF)
         character*(*) string, xstring
         integer ik
         logical tbxxxcln

         if (string(1:1).eq.' ') return
         sbuf(1:1) = '<'
         if (xstring(1:1).eq.' ') then
           ik = lastnb(string)
           sbuf(2:ik+1)=string(1:ik)
         else
           ik = lastnb(xstring)
           sbuf(2:ik+1)=xstring(1:ik)
         endif
         sbuf(ik+2:ik+2) = '>'
         pchar = -1
         if (.not.tbxxxcln(sbuf(2:ik+1),ik)) then
           call tbxxwarn(' XML required remapping for '//sbuf(2:ik+1))
         endif
         call tbxxpstr(sbuf(1:ik+2))
         return
         end


C
C
C >>>>>> Put out the given string as an xml close tag
C
C        Note that the string may have embedded blanks and
C        parameters.  Only the first token will be used for close.
C        The second argument is an optional translation to be
C        used if non-blank
C
         subroutine tbxxpxct(string, xstring)
C
         include 'ciftbx.sys'
         character sbuf*(MAXBUF)
         character*(*) string, xstring
         integer ik
         logical tbxxxcln

         if (string(1:1).eq.' ') return
         sbuf(1:2) = '</'
         if (xstring(1:1).eq.' ') then
           do ik = 1,len(string)
             if (string(ik:ik).eq.' ') go to 100
           enddo
           ik = len(string)+1
 100       ik = ik-1
           sbuf(3:ik+2)=string(1:ik)
         else
           do ik = 1,len(xstring)
             if (xstring(ik:ik).eq.' ') go to 200
           enddo
           ik = len(xstring)+1
 200       ik = ik-1
           sbuf(3:ik+2)=xstring(1:ik)
         endif
         sbuf(ik+3:ik+3) = '>'
         if (.not.tbxxxcln(sbuf(3:ik+2),ik)) then
           call tbxxwarn(' XML required remapping for '//sbuf(3:ik+2))
         endif
         pchar = -1
         call tbxxpstr(sbuf(1:ik+3))
         return
         end

C
C
C
C
C
C >>>>>> Put the string into the output CIF buffer
C
         subroutine tbxxpstr(string)
C
         integer    lastnb
         include   'ciftbx.sys'
         SAVE
         character  string*(*),temp*(MAXBUF),bfill*(MAXBUF)
         character  temp2*(MAXBUF)
         integer    i,ii,mxline,ioffst,ifree,icpos,itpos
         integer    ixpos,ixtpos,it,im,kbin,kpass
         integer    lstring
         logical    pflush,waslop
         data       waslop /.false./

C
CDBG     print *,' entry to tbxxpstr, pchar, pcharl, string'
CDBG     print *, pchar, ', ',pcharl,', ', string
         bfill = ' '
         lstring = min(MAXBUF,lastnb(string))
         mxline=MAXBUF
         if(line_.gt.0) mxline=min(line_,MAXBUF)
         if(pfold_.gt.0) then
           if (pfold_ .lt. lprefx+lstring) then
             call tbxxwarn('Invalid value of pfold_, reset')
             pfold_ = min(line_,lprefx+lstring)
           endif
           mxline=min(mxline,pfold_)
         endif
         temp(1:lstring)=string
         temp2=temp
         pflush=.false.
         if(pchar.lt.0) pflush=.true.
C
         do 100 i=lstring,1,-1
         if(temp(i:i).eq.' ')              goto 100
         if(ptabx_.and.temp(i:i).eq.tab) goto 100
         goto 110
100      continue
         i=0
         it=i
C
C....... Organise the output of loop_ items
C
110      if(i.eq.0)             goto 130
         if(i.eq.1.and.string.eq.tab) goto 130
         if(i.eq.1.and.string.eq.char(0)) then
           pcharl=MAXBUF
           goto 200
         endif
         if((.not.xmlout_).and.temp(1:1).eq.'#')   goto 130
         if(xmlout_.and.temp(1:1).eq.'<') go to 130
         if(ploopf.eq.'yes')    goto 130
         if(ptextf.eq.'yes')    goto 130
         if(ploopn.le.0)        goto 130
         if(pdepth_.gt.0)       goto 130
         ploopc=ploopc+1
         if((align_.or.tabl_).and.pdepth_.eq.0 )then
           if(ploopc.gt.ploopn) then
             if(pcharl.gt.lprefx) pflush=.true.
             ploopc=1
             if(pchar.gt.0) pchar=lprefx+1
           endif
           if(pchar.lt.0)    goto 130
           if(tabl_) then
           kbin=(mxline-lprefx)/8
           if(ploopn.lt.kbin) then
             if(kbin/(ploopn+1).gt.1) then
             pchar=9+lprefx+
     *         (ploopc-1)*8*(kbin/(ploopn+1))
             else
             pchar=1+lprefx+
     *         (ploopc-1)*8*(kbin/ploopn)
             endif
           else
             if(ploopc.le.kbin) then
               pchar=1+lprefx+(ploopc-1)*8
             else
               kpass=(ploopc-kbin-1)/(kbin-1)+1
               pchar=2*kpass+1+lprefx+
     *           mod(ploopc-kbin-1,kbin-1)*8
             endif
           endif
           else
             if(ptabx_) then
             icpos=1
             itpos=1
120          ixpos = 0
             if (icpos.le.i) ixpos=index(temp(icpos:i),tab)
             ixtpos=(pchar+itpos-1+ixpos)
             ixtpos=((ixtpos+7)/8)*8
             if(ixpos.gt.0) then
               if(ixpos.gt.1) then
                 temp2(itpos:ixtpos-pchar+1)=temp(icpos:ixpos-1)
               else
                 temp2(itpos:ixtpos-pchar+1)=' '
               endif
               icpos=ixpos+1
               itpos=ixtpos+2-pchar
               if(icpos.le.i) goto 120
               it=itpos-1
             else
               if(icpos.le.i) then
                 temp2(itpos:itpos+i-icpos)=temp(icpos:i)
                 it=itpos+i-icpos
               endif
             endif
             endif
             if((pchar+i).gt.mxline+1.or.
     *          (ptabx_.and.pchar+it.gt.mxline+1)) then
               if(pcharl.gt.lprefx)pflush=.true.
               pchar=lprefx+1
             endif
           endif
         else
           if(ploopc.le.ploopn)   goto 130
           ploopc=mod(ploopc-1,ploopn)+1
         endif
C
C....... Is the buffer full and needs flushing?
C
130      if(i.eq.1.and.string.eq.tab) then
           if(pcharl.gt.lprefx) then
             if(obuf(pcharl:pcharl).eq.' ') pcharl=pcharl-1
           endif
         endif
         if(pdepth_.gt.0
     *     .and.string(1:1).eq.'#'
     *     .and. pcharl.gt.lprefx) then
           if (obuf(pcharl:pcharl).ne.' ') then
             pcharl = pcharl+1
             obuf(pcharl:pcharl) = ' '
           end if
         end if
         if(pchar.le.pcharl.and.pcharl.gt.lprefx) pflush=.true.
         pchar=max(lprefx+1,pchar)
         if (string.ne.'('
     *     .and.string.ne.'['
     *     .and.string.ne.'{'
     *     .and.string.ne.')'
     *     .and.string.ne.']'
     *     .and.string.ne.'}'
     *     .and.string.ne.','
     *     .and.string.ne.':'
     *     .and.(ploopf.eq.'yes'.or.ploopn.le.0).and.tabl_)
     *     pchar=((pchar-lprefx+6)/8)*8+1+lprefx
         if(ptabx_) then
           icpos=1
           itpos=1
135        ixpos=0
           if(icpos.le.i) ixpos=index(temp(icpos:i),tab)
           ixtpos=(pchar+itpos-1+ixpos)
           ixtpos=((ixtpos+7)/8)*8
           if(ixpos.gt.0) then
             if(ixpos.gt.1) then
               temp2(itpos:ixtpos-pchar+1)=temp(icpos:ixpos-1)
             else
               temp2(itpos:ixtpos-pchar+1)=' '
             endif
             icpos=ixpos+1
             itpos=ixtpos+2-pchar
             if(icpos.le.i) goto 135
             it=itpos-1
           else
             if(icpos.le.i) then
               temp2(itpos:itpos+i-icpos)=temp(icpos:i)
               it=itpos+i-icpos
             endif
           endif
         endif
         if((pchar+i).gt.mxline+1.or.
     *     (ptabx_.and.pchar+it.gt.mxline+1)) then
            pflush=.true.
            pchar=mxline+1-i
            if (xmlout_) pchar = 1
            if (pdepth_.gt.0) then
              pchar = 1
              if (pposbrkstk(pdepth_)+i.lt.mxline)
     *          pchar = pposbrkstk(pdepth_)+1
            end if
            pchar=max(lprefx+1,pchar)
         endif
         if(.not.pflush)  goto 150
         if(pcharl.gt.lprefx) then
           if(waslop.or.(.not.tabl_)) goto 145
           ioffst=0
           pcharl=max(lastnb(obuf(1:pcharl)),lprefx+1)
           ifree=mxline-pcharl
           if(ifree.gt.0) then
           im=numtab+2
           if(numtab.gt.0.and.numtab.le.MAXTAB) then
             if(obuf(itabp(numtab):itabp(numtab)).eq.'#')
     *         im=im-1
           endif
           if(ifree.ge.16.and.im.lt.4.and.
     *       (obuf(1+lprefx:1+lprefx).ne.'#'
     *        .and.((.not.xmlout_).or.(
     *             obuf(1+lprefx:1+lprefx).ne.'<'
     *        .and.obuf(1+lprefx:1+lprefx).ne.']'))
     *        .and.obuf(1+lprefx:1+lprefx).ne.';'
     *        .and.obuf(1+lprefx:1+lprefx).ne.'_'
     *        .and.obuf(1+lprefx:1+lprefx).ne.' '
     *        .and.obuf(1+lprefx:5+lprefx).ne.'data_'
     *        .and.obuf(1+lprefx:5+lprefx).ne.'save_'
     *        .and.obuf(1+lprefx:5).ne.'loop_')) then
             temp(1+lprefx:pcharl)=obuf(1+lprefx:pcharl)
             obuf(1+lprefx:pcharl+8)=
     *         bfill(1:8)//temp(1+lprefx:pcharl)
             ioffst = 8
             ifree=ifree-8
             pcharl=pcharl+8
           endif
           do ii=1,min(MAXTAB,numtab)
             icpos=itabp(ii)+ioffst
             if(icpos.gt.pcharl)   goto 145
             if(im.lt.4) then
             itpos=(max(icpos-lprefx,
     *         ii*(mxline-lprefx)/im)+6)/8
             itpos=itpos*8+1+lprefx
             else
             itpos=(max(icpos-lprefx,
     *         ii*(mxline-lprefx)/im)+4)/6
             itpos=itpos*6+1+lprefx
             endif
             if((obuf(icpos:icpos).eq.''''.or.
     *          obuf(icpos:icpos).eq.'"').and.
     *          itpos.gt.icpos) itpos=itpos-1
             if(itpos-icpos.gt.ifree) itpos=icpos+ifree
             if(itpos.gt.icpos) then
               temp(1:pcharl-icpos+1)=
     *           obuf(icpos:pcharl)
               if(i.lt.numtab) then
                 ixpos=itabp(ii+1)+ioffst
                 if(ixpos.gt.icpos+itpos-icpos+1) then
                   if(obuf(ixpos-(itpos-icpos+1):ixpos-1).eq.
     *               bfill(1:itpos-icpos+1)) then
                     temp(ixpos-itpos+1:pcharl-itpos+1)=
     *               obuf(ixpos:pcharl)
                     pcharl=pcharl-(itpos-icpos)
                   endif
                 endif
               endif
               obuf(icpos:pcharl+itpos-icpos)=
     *           bfill(1:itpos-icpos)//temp(1:pcharl-icpos+1)
               ifree=ifree-(itpos-icpos)
               ioffst=ioffst+itpos-icpos
               pcharl=pcharl+itpos-icpos
             endif
             if(ifree.le.0)      goto 145
           enddo
           endif
145        pcharl=max(1,lastnb(obuf))
           write(outdev,'(a)') obuf(1:pcharl)
         else
           if(precn_.gt.0) then
           if(lprefx.gt.0) then
           write(outdev,'(a)') obuf(1:lprefx)
           else
           write(outdev,'(a)')
           endif
           else
           precn_=precn_-1
           endif
         endif
         waslop=.false.
         precn_=precn_+1
         do ii = 1,MAXTAB
           itabp(ii)=0
         enddo
         numtab=0
         if(lprefx.gt.0) then
           obuf=prefx(1:lprefx)
         else
           obuf=' '
         endif
C
C....... Load the next item into the buffer
C
150      pcharl=pchar+i
         if(ptabx_) pcharl=pchar+it
         waslop= ploopf.eq.'no '.and.ploopn.gt.0.and.align_
         if(i.eq.0) then
           if(pcharl.eq.lprefx+1.and.
     *       obuf(lprefx+1:lprefx+1).eq.' ') pcharl=pcharl-1
             pchar=pcharl+1
             if (pdepth_.gt.0) then
               if(pcharl.gt.lprefx+1
     *           .and.obuf(pcharl:pcharl).eq.' ') then
                 pchar = pcharl
                 pcharl = pcharl-1
               end if
             endif
           goto 200
         endif
         if(ptabx_) then
           obuf(pchar:pcharl)=temp2(1:it)
         else
           if(string.eq.tab) pcharl=pcharl-1
           obuf(pchar:pcharl)=string(1:i)
         endif
         if(pchar.gt.1+lprefx) then
           numtab=numtab+1
           if(numtab.le.MAXTAB) itabp(numtab)=pchar
         endif
         pchar=pcharl+1
         if (pdepth_.gt.0.and.pcharl.gt.lprefx+1
     *     .and.obuf(pcharl:pcharl).eq.' ') then
           pchar = pcharl
           pcharl = pcharl-1
         endif
         if(pchar.gt.mxline+2) then
           if (pfold_.eq.0) then
             call tbxxwarn(' Output CIF line longer than line_')
           else
             call tbxxwarn(
     *         ' Output CIF line longer than line_ or pfold_')
           endif
         endif
C
200      continue
C
CDBG     print *,' exit from tbxxpstr, pchar, pcharl, obuf'
CDBG     print *, pchar, ', ',pcharl,', ', obuf(1:pcharl)

         return
         end
C
C
C
C
C
C >>>>>> Convert the number and esd to string nnnn(m), limited
C        by relative precision prec
C
         subroutine tbxxpnum(numb,sdev,prec)
C
         include   'ciftbx.sys'
         character  string*30,temp*30,c*1,sfmt*8
         double precision numb,sdev,prec,xxnumb,xsdev,slog
         integer    i,iexp,ifp,ii,jj,j,jlnz,jn,kexp,m,ixsdev,islog
         integer    kdecp,ibexp,lexp
C
         kdecp=0
         lexp = 0
         jn = 0
         if (sdev.gt.abs(numb)*prec) then
           if (iabs(esdlim_).ne.esdcac) then
C
C            determine the number of digits set by esdlim_
C
             if (iabs(esdlim_).lt.9 .or.iabs(esdlim_).gt.99999) then
               call tbxxwarn(' Invalid value of esdlim_ reset to 19')
               esdlim_ = 19
             endif
C
C            determine the number of esd digits
C
             esddigx = int(1.+alog10(float(iabs(esdlim_))))
             esdcac = iabs(esdlim_)
           endif
C
C          if esdlim_ < 0, validate pesddig_
C
           if (esdlim_.lt. 0 )then
             if (pesddig_.lt.0 .or. pesddig_.gt.5) then
               call tbxxwarn(' Invalid value of pesddig_ reset to 0')
               pesddig_ = 0
             endif
           endif
C
C          determine kexp, the power of 10 necessary
C          to present sdev as an integer in the range
C          (esdlim_/10,esdlim_] or [1,-esdlim_] if esdlim_ < 0
C
           slog = dlog10(sdev)
           islog = int(slog+1000.)
           islog = islog-1000
           kexp = -islog+esddigx
C
C          Adjust exponent kexp, so that sdev*10**kexp
C          is in the interval (esdlim_/10,esdlim_] or [1,-esdlim_]
C
 20        if (kexp.lt.minexp) then
             call tbxxwarn(' Underflow of esd')
             ixsdev = 0
             go to 30
           endif
           if (kexp.gt.-minexp) then
             call tbxxwarn(' Overflow of esd')
             ixsdev = 99999
             go to 30
           endif
           xsdev = sdev*10.D0**kexp
           ixsdev = int(xsdev+.5)
           if (ixsdev.gt.iabs(esdlim_)) then
             kexp = kexp -1
             go to 20
           endif
           if (ixsdev.lt.(iabs(esdlim_)+5)/10) then
             kexp = kexp+1
             go to 20
           endif
C
C          lexp holds the number of trailing zeros which may be
C          sacrificed in the esd if the number itself has
C          trailing zeros in the fraction which is permitted if
C          esdlim_ is negative
C
C          If esdlim_ is negative and pesddig_ is .gt.0,
C          pesddig_ will be used to force the number of digits
C          in which case lexp has the number of digits that
C          must be sacrificed (lexp > 0) or zeros to add (lexp < 0)
C
           lexp=0
           if(esdlim_.lt.0) then
             if(pesddig_.gt.0) then
25             continue
               if(ixsdev*10**(-lexp).ge.10**(pesddig_))then
                 if(lexp.gt.0)
     *             ixsdev=ixsdev-5*10**(lexp-1)
                 ixsdev=ixsdev+5*10**lexp
                 lexp=lexp+1
                 goto 25
               endif
               if(ixsdev.lt.10**(pesddig_-1+lexp)
     *           .and.lexp.gt.0) then
                 if(ixsdev*10**(-lexp).le.iabs(esdlim_))then
                   lexp =lexp-1
                   if(lexp.ge.0) then
                     ixsdev=ixsdev-5*10**lexp
                   endif
                   if(lexp.gt.0) then
                     ixsdev=ixsdev+5*10**(lexp-1)
                   endif
                   goto 25
                 endif
               endif
               kexp=kexp-lexp
               ixsdev = ixsdev/(10**lexp)
               lexp=0
             else
             do ii = 1,4
               if(mod(ixsdev,10**ii).ne.0) go to 30
               lexp = ii
             enddo
             endif
           endif
C
C          We need to present the number to the same scaling
C          at first, but will adjust to avoid Ennn notation
C          if possible
C
 30        xxnumb = dabs(numb)*10.d0**kexp+.5
           if(xxnumb*prec .gt.1.D0) then
             call tbxxwarn(' ESD less than precision of machine')
             ixsdev=0
           endif
           if(numb.lt.0.d0) xxnumb = -xxnumb
           write(string,ndpfmt)xxnumb
           if(xxnumb.lt.1.d0 .and. xxnumb.ge.0.d0)
     *        string='                         0.0E0'
           if(xxnumb.gt.-1.d0 .and. xxnumb.lt.0.d0)
     *        string='                        -0.0E0'
C
C          Extract the power of 10
C
           iexp = 0
           ibexp = 0
           do ii = 0,4
             i = 30-ii
             c = string(i:i)
             m = index('0123456789',c)
             if (m.gt.0) then
               iexp = iexp+(m-1)*10**(ii-ibexp)
             else
               if (c.eq.' ') then
                 ibexp = ibexp+1
               else
               if (c.eq.'-') iexp=-iexp
               goto 40
               endif
             endif
           enddo
           call tbxxerr(' Internal error in tbxxpnum')
C
C          Scan the rest of the string shifting the
C          decimal point to get an integer
C
40         ifp = 0
           j=1
           do ii = 1,i-1
           c = string(ii:ii)
           if (c.ne.' ')then
             m=index('0123456789+-',c)
             if(m.ne.0) then
               temp(j:j)=c
               if(j.gt.1.or.c.ne.'0')j=j+1
               if(j.eq.3.and.temp(1:2).eq.'-0')j=j-1
               if(ifp.ne.0)then
                 iexp=iexp-1
                 if(iexp.le.0) goto 50
               endif
             else
               if(c.eq.'.') then
                 ifp=1
                 if(iexp.le.0) goto 50
               endif
             endif
           endif
           enddo
C
C          The string from 1 to j-1 has an integer
C          If iexp < 0, we present a 0.  If iexp > 0
C          we pad with zeros
C
50         if(j.eq.2 .and. temp(1:1).eq.'-') then
              temp(1:2)='-0'
              j=3
              iexp=0
           endif
           if(j.eq.1 .or .iexp.lt.0) then
             temp(1:1)='0'
             j=2
             iexp = 0
             if(xxnumb.lt.0.d0) then
               temp(1:2)='-0'
               j=3
             endif
           endif
           if (iexp.gt.0) then
             do ii = 1,iexp
             temp(j:j)='0'
             j=j+1
             enddo
             iexp=0
           endif
           string=temp(1:j-1)
C
C          We have the number for which the presentation
C          would be nnnnnE-kexp.  If kexp is gt 0, we can
C          decrease it and introduce a decimal point
C
           jj=0
           if(index('0123456789',temp(1:1)).eq.0) jj=1
           if(kexp.gt.0.and.kexp.lt.j-jj+8) then
             if(kexp.lt.j-1) then
               if(plzero_ .and.
     *           j-1-kexp.eq.1.and.temp(1:1).eq.'-') then
                 string=temp(1:j-1-kexp)//'0.'//
     *             temp(j-kexp:j-1)
                 j=j+2
               else
                 string=temp(1:j-1-kexp)//'.'//
     *           temp(j-kexp:j-1)
                 j=j+1
               endif
               kexp = 0
             else
               if(jj.ne.0)string(1:1)=temp(1:1)
               if(plzero_) then
                 string(1+jj:2+jj)='0.'
                 do ii=1,kexp-(j-1-jj)
                   string(2+jj+ii:2+jj+ii)='0'
                 enddo
                 string(3+jj+(kexp-(j-1-jj)):30)=
     *             temp(1+jj:j-1)
                 j=j+2+kexp-(j-1-jj)
               else
                 string(1+jj:1+jj)='.'
                 do ii=1,kexp-(j-1-jj)
                   string(1+jj+ii:1+jj+ii)='0'
                 enddo
                 string(2+jj+(kexp-(j-1-jj)):30)=
     *             temp(1+jj:j-1)
                 j=j+1+kexp-(j-1-jj)
               endif
               kexp=0
             endif
           endif
           kdecp=index(string(1:j-1),'.')
           if(kdecp.gt.0.and.kdecp.lt.j-1.and.lexp.gt.0) then
             jj=0
             do ii = 1,min(lexp,j-1-kdecp)
               c = string(j-ii:j-ii)
               if(c.ne.'0') goto 60
               jj=jj+1
             enddo
60           j=j-jj
             ixsdev=ixsdev/10**jj
             if(.not.pdecp_.and.string(j-1:j-1).eq.'.') then
               j=j-1
               kdecp=0
             endif
           endif
           if(kdecp.eq.0) then
             kdecp=j
             if(pdecp_) then
               if(plzero_.and.
     *           (j.eq.1 .or. (j.eq.2.and.string(1:1).eq.'-'))) then
                 string(j:j)='0'
                 j=j+1
               endif
               string(j:j)='.'
               j=j+1
             endif
           endif
           if(kexp.ne.0) then
             write(temp(1:5),'(i5)') -kexp
             string(j:j)='E'
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
           endif
C
C          if there is a standard deviation
C          append it in parentheses
C
           if(ixsdev.ne.0) then
             write(temp(1:5),'(i5)') ixsdev
             string(j:j)='('
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
             string(j:j)=')'
             j=j+1
           endif
         else
C
C          There is no standard deviation, just write numb
C          But limit to the digits implied by prec
C
           slog = dlog10(min(.1D0,max(prec,dpprc)))
           islog = int(slog+1000.5)
           islog = islog-1000
           kexp = -islog
           write(sfmt,'(5h(D30.,i2,1h))') kexp
           write(temp,sfmt)numb
C
C          Now have the number in the form
C          [sign][0].nnnnnnnnDeee
C          which, while sufficient, is not neat
C          we reformat for the case 0<=eee<=kexp
C
C
C          Extract the power of 10
C
           iexp = 0
           ibexp = 0
           do ii = 0,4
             i = 30-ii
             c = temp(i:i)
             m = index('0123456789',c)
             if (m.gt.0) then
               iexp = iexp+(m-1)*10**(ii-ibexp)
             else
               if (c.eq.' ') then
                 ibexp = ibexp+1
               else
               if (c.eq.'-') iexp=-iexp
               goto 140
               endif
             endif
           enddo
           call tbxxerr(' Internal error in tbxxpnum')
C
C          Scan the rest of the string shifting the
C          decimal point to get a number with exponent 0,
C          if possible
C
140        ifp = 0
           j=1
           do ii = 1,i-1
           jn=ii
           c = temp(ii:ii)
           if (c.ne.' ')then
             m=index('0123456789+-',c)
             if(m.ne.0) then
               string(j:j)=c
               if(j.gt.1.or.c.ne.'0')j=j+1
               if(j.eq.3.and.string(1:2).eq.'-0')j=j-1
               if(ifp.ne.0)then
                 iexp=iexp-1
                 if(iexp.le.0) goto 150
               endif
             else
               if(c.eq.'.') then
                 ifp = -1
                 if(iexp.le.0) goto 150
               endif
             endif
           endif
           enddo
150        if(plzero_ .and.
     *       (j.eq.1 .or.(j.eq.2.and.string(1:1).eq.'-'))) then
             string(j:j)='0'
             j=j+1
           endif
           string(j:j)='.'
           ifp = j
           j = j+1
           jlnz = j-1
           do ii = jn+1,i-1
             c = temp(ii:ii)
             if (c.ne.' ')then
               m=index('0123456789',c)
               if(m.ne.0) then
                 string(j:j)=c
                 j=j+1
                 if(m.ne.1)jlnz=j
                 if(m.eq.1.and.ifp.ge.1.and.
     *             pposdec_.ne.0.and.pposend_.ne.0) then
                   if(j-1-ifp-min(iexp,0).le.pposend_-pposdec_)
     *               jlnz=j
                 endif
               else
                 goto 160
               endif
             endif
           enddo
160        j=jlnz
           if(j.eq.1) then
            string(1:1)='0'
            j=2
           endif
           if(iexp.lt.0.and.iexp.gt.-7.and.ifp.lt.j-1.and.
     *       ifp.ne.0.and.j-ifp-iexp.le.kexp) then
             temp(1:ifp)=string(1:ifp)
             do ii = 1,-iexp
               temp(ifp+ii:ifp+ii) = '0'
             enddo
             temp(ifp-iexp+1:j-iexp-1) = string(ifp+1:j-1)
             j = j-iexp
             iexp=0
             string(1:j-1) = temp(1:j-1)
           endif
           kdecp=index(string(1:j-1),'.')
           if(kdecp.eq.0) then
             kdecp=j
             if(pdecp_) then
               string(kdecp:kdecp)='.'
               j=j+1
             endif
           endif
           if(iexp.ne.0) then
             write(temp(1:5),'(i5)')iexp
             string(j:j)='E'
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
           endif
         endif
C
         if(j.lt.1) then
           string(1:1)='0'
           j=2
         endif
         if(kdecp.lt.1)kdecp=j
         if(pposdec_.ne.0) then
           pchar=lprefx+pposdec_-kdecp+1
         else
           if(pposval_.ne.0)pchar=lprefx+pposval_
         endif
         call tbxxpstr(string(1:j-1))
         return
         end
C
C
C
C
C
C >>>>>> Check dictionary for data name validation
C
         subroutine tbxxdck(name,type,flag,tflag)
C
         include   'ciftbx.sys'
         logical    flag,tflag
         integer    nln
         character  name*(*),temp*(NUMCHAR),
     *              type*4
C
         flag=.true.
         tflag=.true.
         nln = min(len(name),len(temp))
         call tbxxnlc(temp(1:nln),name)
         call hash_find(temp(1:nln),
     *     dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,xdchk)
         if(xdchk.eq.0) goto 150
         if(tcheck.eq.'no ')          goto 200
         if(type.eq.dictyp(xdchk))    goto 200
         if(type.eq.'    ')           goto 200
         if(dictyp(xdchk).eq.'text' .and. type.eq.'char') goto 200
         if(dictyp(xdchk).eq.'char' .and. type.eq.'numb') goto 200
         tflag=.false.
         goto 200
150      flag=.false.
200      continue
         return
         end
C
C
C
C
C
C >>>>>> End of text string
C
         subroutine tbxxeot
C
         include   'ciftbx.sys'
C
         character*3 pqt
         integer pql
         integer lastnb

         pqt = pquote_
         pql = lastnb(pqt)
         if (pqt.eq.' ') pqt = ';'
         if (pqt.eq.'(') pqt = ')'
         if (pqt.eq.'{') pqt = '}'
         if (pqt.eq.'[') pqt = ']'

         if(ptextf.ne.'yes') then
           call tbxxwarn(' Out-of-sequence call to end text block')
           return
         endif
         ptextf='no '
         pchar=-1

         if (xmlout_) then
           if (pqt.eq.';') then
           call tbxxpstr(']]>')
           else
             call tbxxpstr(']]>'//pqt(1:pql))
           endif
           if (ploopn.gt.1) then
             call tbxxpxct(plhead(ploopc+1),plxhead(ploopc+1))
           endif
           if (ploopn.le.0) then
             call tbxxpxct(plhead(1),plxhead(1))
             plhead(1) = ' '
             plxhead(1) = ' '
           endif
         else
           call tbxxpstr(pqt(1:pql))
         endif
         if (pqt.eq.';') then
         call tbxxpstr(char(0))
         else
           call tbxxpstr(' ')
         endif
         return
         end
C
C
C
C
C
C >>>>>> End of bracketed structure detected;
C        close all open levels
C
         subroutine tbxxebkt
C
         include   'ciftbx.sys'
         integer   i
         character*1 cd
         if (pdepth_ .eq. 0) return
         do i = 1,pdepth_
           cd = '}'
           if (pbrackstack(1+pdepth_-i).eq.'(' ) cd = ')'
           if (pbrackstack(1+pdepth_-i).eq.'[' ) cd = ']'
           pchar = max(pcharl,lprefx+pposbrkstk(1+pdepth_-i))
           call tbxxpstr(cd)
         end do
         pdepth_ = 0
         return
         end
C
C
C
C
C
C >>>>>> End of loop detected; check integrity and tidy up pointers
C
         subroutine tbxxelp
C
         include   'ciftbx.sys'
         integer   i
C
         if(ploopn.eq.0)          goto 200
         if(ploopn.eq.-1) then
           if (xmlout_) then
             plcat = ' '
             plxcat = ' '
             plhead(1) = 'DUMMY'
             plxhead(1) = ' '
           else
             call tbxxpstr('_DUMMY')
           endif
           ploopn=1
           ploopc=0
           call tbxxwarn(
     *       ' Missing: missing loop_ name set as _DUMMY')
         endif
         if (xmlout_ .and. ploopn.eq.1 .and.
     *     ploopf.ne.'yes') then
           call tbxxpxct(plhead(2),plxhead(2))
         endif
         if(ploopn.eq.ploopc)     goto 200
         do i=ploopc+1,ploopn
         if (xmlout_) then
           call tbxxpxot(plhead(i+1),plxhead(1+1))
           call tbxxpstr('DUMMY')
           call tbxxpxct(plhead(i+1),plxhead(i+1))
         else
           call tbxxpstr('DUMMY')
         endif
         enddo
         call tbxxwarn(
     *         ' Missing: missing loop_ items set as DUMMY')
         plhead(1) = ' '
         plxhead(1) = ' '
C
200      ploopc=0
         ploopn=0
         if (xmlout_) then
           call tbxxpxct(plhead(1),plxhead(1))
           plhead(1) = ' '
           call tbxxpxct(plcat,plxcat)
           plcat = ' '
         endif
         return
         end
C
C
C
C
C
C
C >>>>>> Set common default values
C
         block data
C
         include   'ciftbx.sys'
         data cifdev     /1/
         data outdev     /2/
         data dirdev     /3/
         data errdev     /6/
         data recbeg_    /1/
         data recend_    /0/
         data loopct     /0/
         data nhash      /0/
         data ndict      /0/
         data nname      /0/
         data nbloc      /0/
         data ploopn     /0/
         data ploopc     /0/
         data xmnxlat    /0/
         data xmdata     /0/
         data rsolidus   /'\\'/
         data ploopf     /'no '/
         data ptextf     /'no '/
         data pfilef     /'no '/
         data testfl     /'no '/
         data textfl     /'no '/
         data vcheck     /'no '/
         data tcheck     /'no '/
         data catchk     /'yes'/
         data parchk     /'yes'/
         data align_     /.true./
         data append_    /.false./
         data tabl_      /.true./
         data tabx_      /.true./
         data ptabx_     /.true./
         data text_      /.false./
         data loop_      /.false./
         data ndcname    /0/
         data ncname     /0/
         data rdprn_     /.false./
         data rdbrc_     /.false./
         data rdbkt_     /.false./
         data rdtq_      /.false./
         data rdrcqt_    /.false./
         data rdcolon_   /.false./
         data save_      /.false./
         data saveo_     /.false./
         data psaveo     /.false./
         data glob_      /.false./
         data globo_     /.false./
         data alias_     /.true./
         data aliaso_    /.false./
         data nblank_    /.false./
         data nblanko_   /.false./
         data decp_      /.false./
         data pdecp_     /.false./
         data lzero_     /.false./
         data plzero_    /.false./
         data xmlout_    /.false./
         data catkey     /NUMDICT*.false./
         data xmlong_    /.true./
         data dchash     /NUMHASH*0/
         data dichash    /NUMHASH*0/
         data dhash      /NUMHASH*0/
         data dcchain    /NUMDICT*0/
         data aroot      /NUMDICT*0/
         data keychain   /NUMDICT*0/
         data ccatkey    /NUMDICT*0/
         data cindex     /NUMBLOCK*0/
         data deindex    /NUMDICT*0/
         data dcindex    /NUMDICT*0/
         data line_      /80/
         data lastch     /0/
         data dictype_   /' '/
         data dicname_   /' '/
         data dicver_    /' '/
         data diccat_    /' '/
         data tagname_   /' '/
         data plcat      /' '/
         data plhead     /NUMLP1*' '/
         data prefx      /' '/
         data file_      /' '/
         data longf_     /1/
         data tbxver_    /'CIFtbx version 4.1.0 29 Nov 2009'/
         data lprefx     /0/
         data esdlim_    /19/
         data esddig_    /0/
         data pesddig_   /0/
         data esdcac     /19/
         data esddigx    /2/
         data esdfmt     /'(e12.2)'/
         data edpfmt     /'(d12.2)'/
         data ndpfmt     /'(d30.14)'/
         data decprc     /1.e-6/
         data dpprc      /1.d-14/
         data decmin     /1.e-37/
         data dpmin      /1.d-307/
         data minexp     /-307/
         data itabp      /MAXTAB*0/
         data jrect      /-1/
         data numtab     /0/
         data recn_      /0/
         data precn_     /0/
         data posnam_    /0/
         data posval_    /0/
         data posdec_    /0/
         data posend_    /0/
         data pposnam_   /0/
         data pposval_   /0/
         data pposdec_   /0/
         data pposend_   /0/
         data quote_     /' '/
         data pquote_    /' '/
         data unfold_    /.false./
         data fold_      /.false./
         data clipt_     /.true./
         data pclipt_    /.true./
         data pfold_     /0/
         data ibkmrk     /MAXBOOK*-1,MAXBOOK*-1,
     *                    MAXBOOK*-1,MAXBOOK*-1,
     *                    MAXBOOK*-1,MAXBOOK*-1/
         data lnametb    /1/
         data nametb     /' '/

         end
C
C
C       change the following include to include 'clearfp_sun.f'
C       for use on a SUN
C
        include 'clearfp.f'

