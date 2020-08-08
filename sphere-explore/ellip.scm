

(use-modules (srfi srfi-1))
; for vector-map
(use-modules (srfi srfi-43))

; guile impl of gaussian distribution.
(define (make-normal-generator center width)
	(lambda () (+ center (* width (random:normal)))))

; Avoid heart-ache of porting to guile...
(load "../srfi/sphere.scm")


; Sample ellipse
(define dim-sizes '#(10 2))

  ;;This is what it should be doing....
  (define gaussg-vec
	 (vector-map
		(lambda (i size)
		  (make-normal-generator 0.0 size))
		dim-sizes))

	; Handy-dandy utility
  (define (l2-norm VEC)
	 (sqrt (vector-fold
				(lambda (idx sum x l) (+ sum (/ (* x x) (* l l))))
				0 VEC dim-sizes)))


; Lets do it...
(define g (make-sphere-generator* '#(10 2)))

(define points (map (lambda (x) (g)) (iota 10000)))

(map l2-norm points)

; Verify they sit on the surface
(fold (lambda (p sum) (+ sum (abs (- 1 (l2-norm p))))) 0 points)

; Sort into clockwise order. Only for 2D ellipsoids...
(define (clockwise pts)
	(sort pts (lambda (a b)
		(if (and (< 0 (vector-ref a 1)) (< 0 (vector-ref b 1)))
			(< (vector-ref b 0) (vector-ref a 0))
			(if (and (< (vector-ref a 1) 0) (< (vector-ref b 1) 0))
				(< (vector-ref a 0) (vector-ref b 0))
				(< (vector-ref b 1) (vector-ref a 1)))))))


; Test the clockwise sort...
(define clock (list
	'#(1 1e-3) '#(0.8 0.2) '#(0.2 0.8)
	'#(0 1) '#(-0.2 0.8) '#(-0.8 0.2) '#(-1 1e-3)
	'#(-1 -1e-3) '#(-0.8 -0.2) '#(-0.2 -0.8)
	'#(0 -1) '#(0.2 -0.8) '#(0.8 -0.2) '#(1 -1e-3)))

(equal? (clockwise clock) clock)

; Place the point sample into clockwise order
(define ordered-points (clockwise points))

; Vector differences
; Example usage: (vector-diff '#( 2 3) '#(0.5 0.7))
(define (vector-diff a b)
	(vector-map (lambda (idx ea eb) (- ea eb)) a b))

; Compute differences between neighbor points
(define (delta pts rv)
	(if (null? (cdr pts)) (reverse! rv)
		(delta (cdr pts) (cons (vector-diff (car pts) (cadr pts)) rv))))

(define diffs (delta ordered-points '()))

; Compute the distances between neighboring points
(define dists (map l2-norm diffs))

; Compute sum of a list of numbers
(define (sum lst) (fold (lambda (x sum) (+ sum x)) 0 lst))

; Compute average of a list of numbers
(define (avg lst) (/ (sum lst) (length lst)))

(define pi 3.14159265358979)
; Expect this to be 1.0
(define perimeter (/ (sum dists) (* 2 pi)))

; Debug utility for gnuplot graphing.
; You can use this to dump a list to a tab-delimited file.
(define (list-to-file lst filename)
	(let ((outport (open-file filename "w")))
		(fold
			(lambda (x i) (format outport "~A	~A\n" i x) (+ i 1))
			1
			lst)
		(close outport)))

(list-to-file dists "dists.dat")


; samples
(avg dists) ; 6.283012635109339e-4 == 2pi / num-points
(avg (take dists 300))  ; 6.667184236324033e-4
(avg (take (drop dists 300) 300)) ; 5.841017658792995e-4

; Compute moving average
(define (moving-avg lst window)
	(map (lambda (offset) (avg (take (drop lst offset) window)))
		(iota (- (length lst) window))))

(define moving-300 (moving-avg dists 300))

(list-to-file moving-300 "moving.dat")

(define b (make-sphere-bork* '#(10 2)))
(define borks (map (lambda (x) (b)) (iota 10000)))
(define ordered-borks (clockwise borks))
(define biffs (delta ordered-borks '()))
(define bists (map l2-norm biffs))
(define bork-300 (moving-avg bists 300))
(list-to-file bork-300 "bork.dat")

