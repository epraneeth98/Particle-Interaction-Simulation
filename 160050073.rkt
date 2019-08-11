#lang racket
(provide buildTree)
(provide calcForces)
(provide moveparticles)
(define (pcm ll)
  (define (helper Nx Ny D LL)
    (if (= (length LL) 0) (list D (cons (/ Nx D) (/ Ny D)))
        (helper (+ Nx (* (vec-x (particle-posn (car LL))) (particle-mass (car LL))))
                (+ Ny (* (vec-y (particle-posn (car LL))) (particle-mass (car LL))))
                (+ D (particle-mass (car LL)))
                (cdr LL))))
  (helper 0 0 0 ll))

(define (assign-num body llx lly rux ruy)
  (if (and (>= (vec-x (particle-posn body)) llx)
           (< (vec-x (particle-posn body)) rux)
           (>= (vec-y (particle-posn body)) lly)
           (< (vec-y (particle-posn body)) ruy)) 1 0))

(define (lessthantwo llx lly rux ruy ll)
  (define (helper llx lly rux ruy num ll)
    (cond [(and (= (length ll) 0) (= num 1)) 1]
          [(and (= (length ll) 0) (= num 0)) 0]
          [(> num 1) #f]
          [#t (helper llx lly rux ruy (+ num (assign-num (car ll) llx lly rux ruy)) (cdr ll))]))
  (helper llx lly rux ruy 0 ll))

(define (filter llx lly rux ruy particles)
  (define (filter-help llx lly rux ruy ans particles)
    (cond [(= 0 (length particles)) ans] 
          [(= 1 (assign-num (car particles) llx lly rux ruy))
        (filter-help llx lly rux ruy (cons (car particles) ans) (cdr particles))]
          [(= 0 (assign-num (car particles) llx lly rux ruy))
           (filter-help llx lly rux ruy ans (cdr particles))]))
  (filter-help llx lly rux ruy '() particles))

(define (make-tree llx lly rux ruy particles)
  (cond [(equal? (lessthantwo llx lly rux ruy particles) #f)
        (gnode (car (pcm particles)) (vec [car (car (cdr (pcm particles)))]
                                          [cdr (car (cdr (pcm particles)))])
               (list (make-tree llx (/ (+ lly ruy) 2) (/ (+ llx rux) 2) ruy (filter llx (/ (+ lly ruy) 2) (/ (+ llx rux) 2) ruy particles))
                     (make-tree (/ (+ llx rux) 2) (/ (+ lly ruy) 2) rux ruy (filter (/ (+ llx rux) 2) (/ (+ lly ruy) 2) rux ruy particles))
                     (make-tree llx lly (/ (+ llx rux) 2) (/ (+ lly ruy) 2) (filter llx lly (/ (+ llx rux) 2) (/ (+ lly ruy) 2) particles))
                     (make-tree (/ (+ llx rux) 2) lly rux (/ (+ lly ruy) 2) (filter (/ (+ llx rux) 2) lly rux (/ (+ lly ruy) 2) particles))))]
        [(equal? (lessthantwo llx lly rux ruy particles) 1) (car (filter llx lly rux ruy particles))]
        [(equal? (lessthantwo llx lly rux ruy particles) 0) '()]))

(define (buildTree initialArea particles)
  (define k (bounding-box particles))
     (make-tree (bbox-llx k) (bbox-lly k) (bbox-rux k) (bbox-ruy k) particles))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (square a) (* a a))
(define (cube a) (* a a a))

(define (sep a b)
  (cond [(gnode? b) (sqrt (+ (square (- (vec-x (particle-posn a)) (vec-x (gnode-posn b))))
                           (square (- (vec-y (particle-posn a)) (vec-y (gnode-posn b))))))]
      [(particle? b) (sqrt (+ (square (- (vec-x (particle-posn a)) (vec-x (particle-posn b))))
                           (square (- (vec-y (particle-posn a)) (vec-y (particle-posn b))))))]))
(define (add a b)
  (vec (+ (vec-x a) (vec-x b)) (+ (vec-y a) (vec-y b))))

(define (add-elements list)
  (define (helper list ans)
    (if (= (length list) 0) ans (helper (cdr list) (add ans (car list)))))
    (helper list (vec 0 0)))

(define (one-vec mult a1 b1 a2 b2)
  (vec (* mult (/ (- a2 a1) (sqrt (+ (square (- a1 a2)) (square (- b1 b2))))))
        (* mult (/ (- b2 b1) (sqrt (+ (square (- a1 a2)) (square (- b1 b2))))))))

;(define (function-list list f)
; (map (lambda (x) (f x)) list))

(define (attraction a b)
  (cond [(particle? b) (one-vec (/ (* g (particle-mass a) (particle-mass b)) (square (sep a b)))
                                (vec-x (particle-posn a)) (vec-y (particle-posn a))
                                (vec-x (particle-posn b)) (vec-y (particle-posn b)))]
        [(gnode? b) (one-vec (/ (* g (particle-mass a) (gnode-mass b)) (square (sep a b)))
                                (vec-x (particle-posn a)) (vec-y (gnode-posn a))
                                (vec-x (particle-posn b)) (vec-y (gnode-posn b)))]))

(define (give-list tree)
  (cond [(gnode? tree) (append* (map (lambda (x) (give-list x)) (gnode-subtrees tree)))]
        [(particle? tree) (list tree)]
        [(equal? '() tree) '()]))

(define (bbl input)
  (- (bbox-rux (bounding-box (give-list input))) (bbox-llx (bounding-box (give-list input)))))

(define (force input object)
  (define qq (bbl input))
  (define (help-me input object)
    (cond [(gnode? input) (if (< (/ (sep object input) qq) theta)
                              (add-elements (map (lambda (x) (help-me x object)) (gnode-subtrees input)))
                              (attraction object input))]
                              
          [(equal? '() input) (vec 0 0)]
          [(particle? input) (if (equal? input object) (vec 0 0)
                                 (attraction object input))]))
  (help-me input object))

(define (calcForces initialArea tree list)
  (define ee (buildTree initialArea list))
  (map (lambda (x) (force ee x)) list))

  (define (accn list)
    (define ee (buildTree list))
    (map (lambda (x) (vec [/ (vec-x (force ee x)) (particle-mass x)] [/ (vec-y (force ee x)) (particle-mass x)])) list))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (moveparticles particles forces)
  (define chettu (buildTree forces particles))
  (map (lambda (x) (particle [particle-mass x]
                             [vec (+ (+ (vec-x (particle-posn x)) {* (vec-x (particle-velocity x)) timeslice})
                                     (* 0.5 (/ (vec-x (force chettu x)) (particle-mass x)) (square timeslice)))
                                  (+ (+ (vec-y (particle-posn x)) {* (vec-y (particle-velocity x)) timeslice})
                                     (* 0.5 (/ (vec-y (force chettu x)) (particle-mass x)) (square timeslice)))]
                             [vec (+ (vec-x (particle-velocity x)) (* timeslice (/ (vec-x (force chettu x)) (particle-mass x))))
                                   (+ (vec-y (particle-velocity x)) (* timeslice (/ (vec-y (force chettu x)) (particle-mass x))))]))
       particles))


; ;; ;;; ;;;; ;;;;; ;;;;;; ;;;;;;; ;;;;;;;;