

(use-modules (srfi srfi-1))
; for vector-map
(use-modules (srfi srfi-43))

; guile impl of gaussian distribution.
(define (make-normal-generator center width)
	(lambda () (+ center (* width (random:normal)))))

(define (make-uniform-generator)
	(lambda () (random:uniform)))

; Avoid heart-ache of porting to guile...
(load "../srfi/sphere.scm")


; Sample ellipse
(define axes '#(10 2))

  ;;This is what it should be doing....
  (define gaussg-vec
	 (vector-map
		(lambda (i size)
		  (make-normal-generator 0.0 size))
		axes))

	; Handy-dandy utility
  (define (ellipse-norm VEC AXES)
	 (sqrt (vector-fold
				(lambda (idx sum x l) (+ sum (/ (* x x) (* l l))))
				0 VEC AXES)))

  (define (l2-norm VEC)
	 (sqrt (vector-fold
				(lambda (idx sum x) (+ sum (* x x)))
				0 VEC)))


; Lets do it...
(define g (make-sphere-generator* '#(10 2)))

(define points (map (lambda (x) (g)) (iota 10000)))

(map ellipse-norm points axes)

; Verify they sit on the surface
(fold (lambda (p sum) (+ sum (abs (- 1 (ellipse-norm p axes))))) 0 points)

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

; Pseudo-perimeter -- expect this to be 2pi
(define circu (sum (map (lambda (d) (ellipse-norm d axes)) diffs)))

; factorial
(define (fact n rv)
	(if (zero? n) rv (fact (- n 1) (* n rv))))

; Double factorial
; https://en.wikipedia.org/wiki/Double_factorial
(define (double-fact n rv)
	(if (<= n 0) rv (double-fact (- n 2) (* n rv))))

; Complete elliptic integral per wikipedia, see the Ivorty& Bessel
; expansion. Here `a` and `b` are the axes.
; https://en.wikipedia.org/wiki/Ellipse
(define (complete-elliptic a b)
	(define rh (/ (- a b) (+ a b)))
	(define h (* rh rh))

	(define (ivory term n twon hn fact-n dfact-n sum)
		(if (< term 1e-10) (+ sum term)
			; (format #t "yo n= ~A term=~A 2^n=~A h^n=~A n!=~A n!!=~A sum=~A\n"
			; n term twon hn fact-n dfact-n sum)
			(ivory
				(/ (* dfact-n dfact-n hn) (* twon twon fact-n fact-n))
				(+ n 1)
				(* 2 twon)
				(* h hn)
				(* (+ n 1) fact-n)
				(* (- (* 2 n) 1) dfact-n)
				(+ term sum))))

	(* pi (+ a b) (+ 1 (/ h 4)
		(ivory (/ (* h h) 64) 3 8 (* h h h) 6 3 0.0))))

; Sum of the intervals
(define perimeter (sum dists))

; Should equal the complete integral.
(define perim-exact
	(complete-elliptic (vector-ref axes 0) (vector-ref axes 1)))

; Should be less than one
(abs (* (/ (- perimeter perim-exact) perim-exact) (length points)))


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
(define (veclist-to-file lst filename)
	(let ((outport (open-file filename "w")))
		(fold
			(lambda (x i) (format outport "~A	~A	~A\n"
				 i (vector-ref x 0) (vector-ref x 1))
				 (+ i 1))
			1 lst)
		(close outport)))


(define b (make-sphere-drop* '#(2 10)))
(define borks (map (lambda (x) (b)) (iota 1000)))
(define ordered-borks (clockwise borks))
(define biffs (delta ordered-borks '()))
(define bists (map l2-norm biffs))
(define bork-300 (moving-avg bists 30))
(list-to-file bork-300 "bork.dat")
(veclist-to-file ordered-borks "bell.dat")

(define bork-300 (moving-avg bists 3000))
(list-to-file bork-300 "bork.dat")

(define t (make-sphere-drop* '#(2 10 6)))
(define first-two
	(lambda () 
		(define (try)
			(define samp (t))
			(if (and (< (vector-ref samp 2) 0.1) (< -0.1 (vector-ref samp 2)))
				(vector (vector-ref samp 0) (vector-ref samp 1))
				(try)))
		(try)))


; Take a two-2d slice out of ellipsoidhaving `axes`. Thickness of the
; slice is `thickness`. Location of the slice is `where`.  Caution:
; may be extremely slow if a high-dimenional vector is given, or if
; thickness is too thin. Will infinite-loop if `where` is not inside
; the ellipsiod.  Returns an n-dimensional point, of which the first
; two coords are unconstrained.
; 
; Usage:
;
; (define gen (make-slicer '#(2 10 6 4) 0.1  '#(0 0 2 1)))
; (gen)
(define (make-slicer axes thickness where)

	; generator of points on n-dimensional ellipsoid
	; (define elli (make-sphere-drop* axes))
	(define elli (make-ellipsoid-generator axes))

	; return true if point is in the slice.
	(define (accept point)

		;; take the slice off-center
		(define diff
			(vector-map (lambda (idx r s) (- r s)) point where))

		;; return #t if the point is in the slice
		(define (ok vec)
			(vector-fold (lambda (idx pass coord)
				(or (< idx 2)
					(and pass (< (- 0 thickness) coord) (< coord thickness))))
				#t vec))

		; test
		(ok diff))

	(define (try)
		(define sample (elli))
		(if (accept sample) sample (try)))

	; return the looper
	try)

; Like above but return 2D point only.
(define (two-d-slicer axes thickness where)
	(define sli (make-slicer axes thickness where))

	(lambda ()
		(define sample (sli))
		(vector (vector-ref sample 0) (vector-ref sample 1))))

; ---------------------------------------------------------

; iota
(define (generator->list f n)
   (map (lambda (junk) (f)) (make-list n)))


(define (gtake GEN N)
   (generator->list GEN N))

(define (generator-for-each LAM LST) (for-each LAM LST))


(define REPS 10000)
(define egen (make-ellipsoid-generator* '#(2 10)))

(define egen (two-d-slicer '#(4 10 16) 0.04  '#(0 0 2)))
(define egen (two-d-slicer '#(4 10 16) 0.0005  '#(0 0 1)))

(define pts (map (lambda (x) (egen)) (iota REPS)))
(define ordered-pts (clockwise pts))
(define diffs (delta ordered-pts '()))
(define dists (map l2-norm diffs))
(define perimeter (sum dists))
(define exp-dist (/ perimeter REPS))
(define norm-dists (map (lambda (x) (/ x exp-dist)) dists))

(define CENT 100)
(define NBINS (* 8 CENT))

(define counts (do-count norm-dists CENT NBINS))
(vector-to-file counts "counts.dat")


  ; Bin-counter containing accumulated histogram.
(define (do-count FLST SCALE NBINS)
    (let ((bin-counts (make-vector NBINS 0.0))
			 (incr (/ SCALE (length FLST))))
     ; Accumulate samples into the histogram.
     (for-each
       (lambda (SAMP)
         (define off SAMP)
         (define offset (if (< off NBINS) off (- NBINS 1)))
         (vector-set! bin-counts offset (+ incr (vector-ref bin-counts offset))))
       (map (lambda (x) (inexact->exact (floor (* x SCALE)))) FLST))
     bin-counts))

(define (vector-to-file vec filename)
  (define (write-vec)
    (for-each
      (lambda (i)
        (define index (+ i 1))
        (define val (vector-ref vec i))
        (display index)
        (display "  ")
        (display val)
        (newline))
      (iota (vector-length vec))))
  (with-output-to-file filename write-vec))

