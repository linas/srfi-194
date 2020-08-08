
(define s 1.001)

   (define 1ms (- 1 s))
   (define (trm n u lg) (* lg (+ 1 (/ (* 1ms u) n))))
   (define (exn lg) (trm 2 (trm 3 (trm 4 1 lg) lg) lg))

(define (acn lg) (/ (- (exp (* 1ms lg)) 1) 1ms))

(define (prtex)
(for-each (lambda (p)
	(set! 1ms (expt 0.1 p))
	(define diff (- (exn lg) (acn lg)))
	(format #t "~A ~A ~A ~A\n" 1ms diff (exn lg) (acn lg)))
	(iota 16 1.5)))

(define q 0.1)

   (define (big-h x) (exn (log (+ q x))))

(define (act x) (/ (- (exp (* 1ms (log (+ q x)))) 1)  1ms))

(define (prt)
(for-each (lambda (x)
	(define diff (- (big-h x) (act x)))
	(format #t "~A ~A ~A\n" diff (big-h x) (act x)))
	(iota 10 1.5)))

----------------

(define (alg y) (/ (log (+ 1 (* y 1ms))) 1ms))

   (define (lg y)
      (define yms (* y 1ms))
      (define (trm n u r) (- (/ 1 n) (* u r)))
      (* y (trm 1 yms (trm 2 yms (trm 3 yms (trm 4 yms 0))))))


(define (prta)
(for-each (lambda (x)
	(define diff (- (lg x) (alg x)))
	(format #t "~A ~A ~A\n" diff (lg x) (alg x)))
	(iota 10 1.5)))



-----------------


(define (aiv y) (- (expt (+ 1 (* y 1ms)) (/ 1 1ms)) q ))

   (define (big-h-inv y)
      (- (exp (lg y)) q))


(define (sets es) (set! s es) (set! 1ms (- 1 es)))

(define (prti)
(for-each (lambda (x)
	(define diff (- (big-h-inv x) (aiv x)))
	(format #t "~A ~A ~A\n" diff (big-h-inv x) (aiv x)))
	(iota 10 1.5)))

