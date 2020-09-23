((nil . (
	 (compile-command . "make -C ../build -j 4")
	 (projectile-project-compilation-cmd . "make -C ../build -j 4")
	 (eval . (setq-local flycheck-gfortran-include-path
			     (list (expand-file-name "build" (projectile-project-root))
				   (expand-file-name "build/src" (projectile-project-root)))))
)))
