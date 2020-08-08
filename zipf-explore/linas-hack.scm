
(use-modules (srfi srfi-27))
(use-modules (srfi srfi-43)) ; for vector-fold
(use-modules (rnrs io ports)) ; for eof-object

 (load "srfi-194-impl.scm")

(use-modules (srfi srfi-64)) ; for test-group
; (use-modules (srfi srfi-158)) ; for generators...
; (load "srfi-194-test.scm")


(load "scaffolding.scm")
; (load "zipf-impl.scm")
(load "zipf-zri.scm")
(load "zipf-test.scm")

