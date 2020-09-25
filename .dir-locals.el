((nil . (
         (compile-command . "make -C ../build -j 4")
         (projectile-project-compilation-cmd . "make -C ../build -j 4")
         (eval . (setq-local flycheck-gfortran-args
                             (list (concat "-J" (projectile-project-root) (file-name-as-directory "build/") (file-relative-name (file-name-directory buffer-file-name) (projectile-project-root))))))
)))
