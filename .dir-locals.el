((nil . ((eval . (set (make-local-variable 'my-project-path)
                      (file-name-directory
                       (let ((d (dir-locals-find-file ".")))
                         (if (stringp d) d (car d))))))
 	 (compile-command . (format "make -C %s -j 4" (concat my-project-path "build"))))))
