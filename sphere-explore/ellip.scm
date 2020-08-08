

(use-modules (srfi srfi-1))
; for vector-map
(use-modules (srfi srfi-43))

; guile impl of gaussian distribution.
(define (make-normal-generator center width)
	(lambda () (+ center (* width (random:normal)))))

(load "../srfi/sphere.scm")

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

; Place the point sample into order
(define ordered-points (clockwise points))

(define (vector-diff a b)
	(vector-map (lambda (idx ea eb) (- ea eb)) a b))

(vector-diff '#( 2 3) '#(0.5 0.7))

; Compute differences between neighbor points
(define (delta pts rv)
	(if (null? (cdr pts)) (reverse! rv)
		(delta (cdr pts) (cons (vector-diff (car pts) (cadr pts)) rv))))

(define diffs (delta ordered-points '()))

; Compute the distances between neighboring points
(define dists (map l2-norm diffs))

; Debug utility for gnuplot graphing.
; You can use this to dump a list to a tab-delimited file.
(define (list-to-file lst filename)
	(let ((outport (open-file filename "w")))
		(fold
			(lambda (x i) (format outport "~A	~A\n" i x) (+ i 1))
			1
			lst)
		(close outport)))

(list-to-file dists "dists.csv")

