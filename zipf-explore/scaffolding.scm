
; iota
(define (generator->list f n)
   (map (lambda (junk) (f)) (make-list n)))


(define (gtake GEN N)
	(generator->list GEN N))

(define (generator-for-each LAM LST) (for-each LAM LST))
