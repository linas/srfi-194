

; for vector-map
(use-modules (srfi srfi-43))

; guile impl of gaussian distribution.
(define (make-normal-generator center width)
	(lambda () (+ center (* width (random:normal)))))

; (load "../srfi/sphere.scm")

  (define dim-sizes '#(10 2))
  (define gaussg-vec
    (vector-map
      (lambda (i size)
        (make-normal-generator 0.0 (sqrt size)))
      dim-sizes))




(define g (make-sphere-generator* '#(10 2)))


(vector-map  (lambda (i size) 42) '#(10 2))
