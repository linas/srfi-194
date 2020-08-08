

(use-modules (srfi srfi-1))
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
        (make-normal-generator 0.0 size))
      dim-sizes))

  (define (l2-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x l)
              (+ sum (/ (* x x)
                        (* l l)
                        )))
            0
            VEC
            dim-sizes)))


(define g (make-sphere-generator* '#(10 2)))

(vector-map  (lambda (i size) 42) '#(10 2))

(define points (map (lambda (x) (g)) (iota 10000)))

(map l2-norm points)

; verify they sit on the surface
(fold (lambda (p sum) (+ sum (abs (- 1 (l2-norm p))))) 0 points)

; Sort into clockwise order
(define (clockwise pts)
	(sort pts (lambda (a b)
		(if (and (< 0 (vector-ref a 1)) (< 0 (vector-ref b 1)))
			(< (vector-ref b 0) (vector-ref a 0))
			(if (and (< (vector-ref a 1) 0) (< (vector-ref b 1) 0))
				(< (vector-ref a 0) (vector-ref b 0))
				(< (vector-ref b 1) (vector-ref a 1)))))))


; test
(define clock (list
	'#(1 1e-3) '#(0.8 0.2) '#(0.2 0.8)
	'#(0 1) '#(-0.2 0.8) '#(-0.8 0.2) '#(-1 1e-3)
	'#(-1 -1e-3) '#(-0.8 -0.2) '#(-0.2 -0.8)
	'#(0 -1) '#(0.2 -0.8) '#(0.8 -0.2) '#(1 -1e-3)))

(equal? (clockwise clock) clock)

(define ordered-points (clockwise points))
