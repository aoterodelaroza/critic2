((nil . (
         (compile-command . "make -C ../build -j 4")
         (projectile-project-compilation-cmd . "make -C ../build -j 4")
         (eval . (setq-local flycheck-gfortran-args
			     (append (list (concat "-I" (projectile-project-root) "build/src") (concat "-J" (projectile-project-root) "build/src"))
				     (mapcar (lambda (str) (concat "-I" (projectile-project-root) "build/src/" str))
				     (split-string (shell-command-to-string
						    (concat "find " (projectile-project-root) "build/src" " -maxdepth 1 -type d ! -name '*CMakeFiles*' -printf '%P\n'")))))))
)))
