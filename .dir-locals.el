;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

;; ((c++-mode . (flycheck-gcc-language-standard . "c++11"))
((nil
  (projectile-project-run-cmd . "./build/FESR")
  (projectile-project-compilation-cmd . "make -C /Users/knowledge/Developer/PhD/FESR/build/ -j8 default_target")
  (projectile-project-test-cmd . "./build/runUnitTests")))
