;; Ordinary Différential Equations Solvers - Serial RLC Example using state2ode converter 
;;
;; Methods based on: Gerald-Wheatley, 2003, Applied Numerical Analysis, Seventh Edition, 
;; Pearson Education, ISBN 0-321-13304-8.
;;
;;
;; Copyright (c) 2009 Stephane Rousseau (stephaner gmail.com). All rights reserved.  
;;
;; The use and distribution terms for this software are covered by the Eclipse Public License 1.0,
;; which can be found in the file epl-v10.html at the root of this distribution. 
;; By using this software in any fashion, you are agreeing to be bound by the terms of this license.
;; You must not remove this notice, or any other, from this software.

(use 'com.konato.ode)

(defn rlcvcstate [R L Cp U T0 TF DT]
"Test the Runge-Kutta order 4 method for vc RLC circuit using State2ode transformation and odesolve.
Example from: csserver.evansville.edu/~blandfor/EE210/MatLab/MLEx5.pdf"
(let [
      A [[(/ (* -1.0 R) L) (/ (* -1.0 R) L)] [(/ 1 (* R Cp)) 0]]
      B [[(/ (* 1.0 R) L)] [0]]
      C [[0 1]]
      D [[0]]
      fps (state2ode A B C D U)
     ]
  (odesolve rk4 (:dotfps fps) {:x1 0.0, :x2 0.0} :t T0 TF DT (:outfps fps))
))

(comment
  (load-file "com_konato_ode/src/examples/rlc.clj")
;; Using data from: csserver.evansville.edu/~blandfor/EE210/MatLab/MLEx5.pdf
  (def resa (rlcvcstate 1000 0.010 0.000001 [#(* 1 (step (:t %)))] 0.0 0.005 0.00001));overdamped
  (def resb (rlcvcstate 200 0.010 0.000001 [#(* 1 (step (:t %)))] 0.0 0.005 0.00001));criticaly damped
  (def resc (rlcvcstate 40 0.010 0.000001 [#(* 1 (step (:t %)))] 0.0 0.005 0.00001));underdamped
)