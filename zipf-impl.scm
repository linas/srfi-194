;
; Create a Zipf random distribution.
;
; Created by Linas Vepstas 24 June 2020
; Nominated for inclusion in srfi-194
; Not yet unit-tested, may be buggy!
;
; Implementation taken from drobilla's May 24, 2017 answer to
; https://stackoverflow.com/questions/9983239/how-to-generate-zipf-distributed-numbers-efficiently
;
; That code is referenced with this:
; "Rejection-inversion to generate variates from monotone discrete
; distributions", Wolfgang Hörmann and Gerhard Derflinger
; ACM TOMACS 6.3 (1996): 169-184
;
; This code works best for large-N distributions; for small-N one
; can precompute an array and cache it.  For an example, see
; https://github.com/opencog/cogutil/blob/master/opencog/util/zipf.h
; The cached version takes longer to create, but once created, it runs
; 2x or 3x faster ... for small N (say, less than 300).  But for large
; N, the cache is blown (literally, L1 D-cache misses on the CPU core).
;
; Example usage:
;    (define zgen (make-zipf-generator 50 1.01))
;    (generator->list zgen 10)
;
;
; ------------------------------------------------------------

; Defaults should be n==int-max and s==1
(define (make-zipf-generator/new n s q)

	; Epsilon to avoid convergence issues with divide-by-zero.
	; This is to avoid explosions due to scheme not having many
	; common elementary functions (that I know of, at any rate).
	(define epsilon 2.0e-5)

	; (exp(x) - 1) / x 
	; This partly works around the issue of scheme not having a native
	; high-precision exp(x)-1 function that is accurate for small x.
	; (that I know of). Uses the epsilon above.
	(define (expxm1bx x)
		(if (< epsilon (abs x))
			(/ (- (exp x) 1) x)
			(+ 1 (* (/ x 2) (+ 1 (* (/ x 3) (+ 1 (/ x 4)))))))
	)

	; log(1 + x) / x 
	; As before, uses the epsilon to work around the missing high-precision
	; log(1+x) function. Bummer.
	(define (log1pxbx x)
		(if (< epsilon (abs x))
			(/ (log (+ 1 x)) x)
			(- 1 (* x (- (/ 1 2) (* x (- (/ 1 3) (* x (/ 1 4))))))))
	)

	; The hat function h(x) = 1/(x^s) 
	(define (hat x) 
		(expt (+ x q) (- s))
	)

	(define 1ms (- 1 s))
	(define logq (log (+ 1 q)))
	(define big-hq (* logq (expxm1bx (* 1ms logq))))

	; The integral of hat(x)
	; H(x) = log(x) if s == 1, (x^(1-s) - 1)/(1 - s) otherwise.
	; Note the numerator is one less than in the paper order to
	; work with all positive s.
	(define (big-h x)
		(define logx (log (+ x q)))
		(- (* logx (expxm1bx (* 1ms logx))) big-hq)
	)

	; The inverse function of H(x)
	(define (big-h-inv y)
		(define t (+ y big-hq))
		(define u (* t 1ms))
		(- (exp (* t (log1pxbx u))) q)
	)

	; Clamp x to [lo, hi].
	(define (clamp x lo hi)
		(max lo (min x hi))
	)

	; For zero q, have a 4x performance improvement by offsetting.
	(define big-h-x1
		(if (< 0.1 (abs q)) (big-h 0.5) (- (big-h 1.5) 1)))

	(define big-h-n (big-h (+ n 0.5)))

	(define dist (make-random-real-generator big-h-x1 big-h-n))

	; Attempt to hit the dartboard. Return #f if we fail,
	; otherwise return an integer between 1 and n.
	(define (try)
		(define u (dist))
		(define x (big-h-inv u))
		(define flt-k (clamp (round x) 1 n))
		; Convert to integer. This is guile-specific, I don't know
		; generic scheme for this.
		(define k (inexact->exact (floor flt-k)))

		(if (>= u (- (big-h (+ k 0.5)) (hat k))) k #f))

	; Did we hit the dartboard? If not, try again.
	(define (loop-until)
		(define k (try))
		(if k k (loop-until)))

	; Return the generator.
	loop-until
)

*unspecified*
