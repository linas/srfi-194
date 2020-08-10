;
; sphere.scm
; Uniform distributions on a sphere, and a ball.
; Submitted for inclusion in srfi-194
;
; Algorithm based on BoxMeuller as described in
; http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
;


; make-sphere-generator N - return a generator of points uniformly
; distributed on an N-dimensional sphere.
; This implements the BoxMeuller algorithm, that is, of normalizing
; N+1 Gaussian random variables.
(define (make-sphere-generator arg)
  (cond
    ((integer? arg) (make-sphere-generator* (make-vector (+ 1 arg) 1.0)))
    ((vector? arg) (make-sphere-generator* arg))
    (else (error "expected argument to either be a number (dimension), or vector (axis length for the dimensions)"))))

(define (make-sphere-generator* dim-sizes)
  (define gaussg-vec
    (vector-map
      (lambda (idx size)
        (make-normal-generator 0.0 size))
      dim-sizes))
  ; Banach l2-norm aka root-mean-square distance.
  (define (l2-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x l)
              (+ sum (/ (* x x)
                        (* l l)
                        )))
            0
            VEC
            dim-sizes)))

  (lambda ()
    (define vect
      (vector-map
        (lambda (idx gaussg)
          (gaussg))
        gaussg-vec))
    (define norm (/ 1.0 (l2-norm vect)))
    (vector-map (lambda (idx x)
                  (* x norm))
                vect)))

(define (make-sphere-bork* dim-sizes)
  ; Banach l2-norm aka root-mean-square distance.
  (define (l2-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x) (+ sum (* x x)))
            0
            VEC)))
  (define gaussg-vec
    (vector-map
      (lambda (idx size)
        (make-normal-generator 0.0 size ))
      dim-sizes))

  ; Distance to ellipse.
  (define (ellipse-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x l)
              (+ sum (/ (* x x)
                        (* l l)
                        )))
            0
            VEC
            dim-sizes)))

  (define (sgn x) (if (< 0 x) 1 -1))
  (lambda ()
    (define vect
      (vector-map
        (lambda (idx gaussg)
          (define x (gaussg))
          (* 1 x))
        gaussg-vec))
    (define norm (/ 1.0 (l2-norm vect)))
    (vector-map (lambda (idx x)
                  (* x norm (vector-ref dim-sizes idx)))
                vect)))

(define (make-sphere-bark* dim-sizes)
  (define gaussg-vec
    (vector-map
      (lambda (idx size) (make-normal-generator 0 1 ))
      dim-sizes))

  ; Distance to ellipse.
  (define (ellipse-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x l) (+ sum (* x x l l )))
            0
            VEC
            dim-sizes)))

  (lambda ()
    (define vect
      (vector-map
        (lambda (idx gaussg) (gaussg))
        gaussg-vec))

    (define norm (/ 1.0 (ellipse-norm vect)))
    (vector-map
         (lambda (idx x)
              (define a (vector-ref dim-sizes idx))
              (* x norm a a ))
                vect)))

(define (make-sphere-birk* dim-sizes)
  ; Banach l2-norm aka root-mean-square distance.
  (define (l2-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x) (+ sum (* x x)))
            0
            VEC)))
  (define gaussg-vec
    (vector-map
      (lambda (idx size)
        (make-normal-generator 0.0 size ))
      dim-sizes))

  ; Distance to ellipse.
  (define (cube-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x l)
              (+ sum (/ (* x x)
                        (* l l )
                        )))
            0
            VEC
            dim-sizes)))

  (define (ellipse-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x l)
              (+ sum (/ (* x x)
                        (* l l)
                        )))
            0
            VEC
            dim-sizes)))

  (lambda ()
    (define vect
      (vector-map
        (lambda (idx gaussg)
          (define x (gaussg))
          (* 1 x))
        gaussg-vec))
    (define cnorm (/ 1.0 (cube-norm vect)))
  (define cube-vec
    (vector-map (lambda (idx x)
          (define a (vector-ref dim-sizes idx))
                  (* x cnorm (/ 1 (* a)) ))
                vect))
    (define norm (/ 1.0 (ellipse-norm cube-vec)))
    (vector-map (lambda (idx x)
          (define a (vector-ref dim-sizes idx))
                  (* x norm  ))
                cube-vec))
  )

; -----------------------------------------------
(define (make-sphere-drop-exp* dim-sizes)
  (define gaussg-vec
    (vector-map
      (lambda (idx size) (make-normal-generator 0 1))
      dim-sizes))

  ; Banach l2-norm
  (define (l2-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x) (+ sum (* x x)))
            0
            VEC)))

	; Generate point on sphere
	(define (sph)
    (define vect 
      (vector-map
        (lambda (idx gaussg) (gaussg))
        gaussg-vec))
    (define norm (/ 1.0 (l2-norm vect)))
    (vector-map
         (lambda (idx x)
              (* x norm))
                vect))


  ; Distance to ellipse.
  (define (ellipse-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x l)
              (+ sum (/ (* x x)
                        (* l l)
                        )))
            0
            VEC
            dim-sizes)))

	(define vol
    (vector-fold
            (lambda (idx prod l) (* prod l))
            1 dim-sizes))

	(define minor
    (vector-fold
            (lambda (idx mino l) (if (< l mino) l mino))
            1e60 dim-sizes))
	(define uni (make-uniform-generator))

; (format #t "maj=~A vol=~A\n" minor vol)
  ; probability of accept
  (define (keep VEC)
	(define cut (* minor (ellipse-norm VEC)))
	(define ran (uni))
; (format #t "yo cut=~A uni=~A\n" cut ran)
	(< ran cut))

	(define (sample)
		(define vect (sph))
		(if (keep vect) vect (sample)))

  (lambda ()
    (define vect (sample))

    (define norm (/ 1.0 (l2-norm vect)))
    (vector-map
         (lambda (idx x)
              (define a (vector-ref dim-sizes idx))
              (* x norm a ))
                vect))
)

; -----------------------------------------------
(define (make-sphere-drop* axes)

	; A vector of normal gaussian generators
	(define gaussg-vec
		(make-vector (vector-length axes) (make-normal-generator 0 1)))

  ; Banach l2-norm of a vector
  (define (l2-norm VEC)
    (sqrt (vector-fold
            (lambda (idx sum x) (+ sum (* x x)))
            0
            VEC)))

	; Generate one point on a sphere
	(define (sph)
		; Sample a point
    (define point
      (vector-map (lambda (idx gaussg) (gaussg)) gaussg-vec))
		; project it to the unit sphere
    (define norm (/ 1.0 (l2-norm point)))
    (vector-map (lambda (idx x) (* x norm)) point))

  ; Distance from origin to the surface of the
  ; ellipsoid along direction VEC
  (define (ellipsoid-dist VEC)
    (sqrt (vector-fold
            (lambda (idx sum x a) (+ sum (/ (* x x) (* a a))))
            0
            VEC
            axes)))

	; Find the shortest axis
	(define minor
    (vector-fold
            (lambda (idx mino l) (if (< l mino) l mino))
            1e308 axes))

	; uniform generator
	(define uni (make-uniform-generator))

  ; Return #t if the POINT can be kept; else must resample.
  (define (keep POINT)
	(< (uni) (* minor (ellipsoid-dist POINT))))

	; Sample once. The returned sample is a point on a sphere
	; of unit radius (we already normed up above).
	(define (sample)
		(define vect (sph))
		(if (keep vect) vect (sample)))

  (lambda ()
	; Rescale to ellipsoid.
    (vector-map
         (lambda (idx x a) (* x a )) (sample) axes))
)

; make-ball-generator N - return a generator of points uniformly
; distributed inside an N-dimensional ball.
; This implements the Harman-Lacko-Voelker Dropped Coordinate method.
(define (make-ball-generator arg)
  (define dim-sizes
    (cond
      ((integer? arg) (make-vector (+ 2 arg) 1.0))
      ((vector? arg) (vector-append arg (vector 1.0 1.0)))
      (else (error "expected argument to either be a number (dimension), or vector (axis length for the dimensions)"))))
  (define N (- (vector-length dim-sizes) 2))
  (define sphereg (make-sphere-generator (+ N 2)))
  ; Create a vector of N+2 values, and drop the last two.
  ; (The sphere already added one, so we only add one more)
  (lambda ()
    (vector-map
      (lambda (el dim-size _) (* el dim-size))
      (sphereg)
      dim-sizes
      (make-vector N #f))))

