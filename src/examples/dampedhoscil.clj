;; Ordinary Différential Equations Solvers - Damped Harmonic Oscillator Example 
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

(defn dampedho [m fv k xini xdotini T0 TF DT]
"Test the Runge-Kutta order 4 method for a classic damped harmonic oscilator."
 (let [
       fps  {:x #(:xdot %), :xdot #(/ (+ (* -1.0 k (:x %)) (* -1.0 fv (:xdot %))) m)}
       outfps {:x #(:x %), :t #(:t %)}
      ]
      (odesolve rk4 fps {:xdot 0.0, :x 0.9} :t T0 TF DT outfps)
))

(defn transferdampedho [m fv k xini xdotini T0 TF DT]
"Test the Runge-Kutto order 4 method and the transfer2ode for a classic damped harmonic oscilator."
  (let [
	DEN [1 (/ k m) (/ fv m)]
	NUM [0 1]
	inp [0]
	;inp [#(/ (:t %) m)]
	fps (transfer2ode NUM DEN inp)
       ]
    (odesolve rk4 (:dotfps fps) {:x2 0.0, :x1 0.9} :t T0 TF DT (:outfps fps))
))

(defn statedampedho [m fv k xini xdotini T0 TF DT]
"Test the Runge-Kutto order 4 method and the state2ode for a classic damped harmonic oscilator."
  (let [
	A [[0 1][(* -1 (/ k m)) (* -1 (/ fv m))]]  
	B [[0][1]]
	C [[1 0]]
	D [[0]]
	U [0]
	fps (state2ode A B C D U)
       ]
    (odesolve rk4 (:dotfps fps) {:x2 0.0, :x1 0.9} :t T0 TF DT (:outfps fps))
))


(comment
  (load-file "com_konato_ode/src/examples/dampedhoscil.clj")
  (def res (dampedho 1 2 100 0.9 0.0 0.0 2.0 0.01));G.A. Kornm,1989, Interactive Dynamic System Simulation, example 1.4
  (dampedho 1 5 6 1.0 0.0 0.0 2.0 0.01); Real and unequal roots
  (dampedho 1 25 2 1.0 0.0 0.0 2.0 0.01); Real and unequal roots
  (dampedho 1 1 1 1.0 0.0 0.0 2.0 0.01); Complex roots
  (dampedho 1 4 4 1.0 0.0 0.0 2.0 0.01); Real and equal roots
  (transferdampedho 1 1 1 1.0 0.0 0.0 2.0 0.01); Complex roots from a transfer representation
  (statedampedho 1 1 1 1.0 0.0 0.0 2.0 0.01); Complex roots from a state representation
  (= 0 (reduce + (map #(- %1 %2) (:x (dampedho 1 1 1 1.0 0.0 0.0 2.0 0.01)) (:y1 (transferdampedho 1 1 1 1.0 0.0 0.0 2.0 0.01)))));error study with transfer2ode conversion
  (= 0 (reduce +(map #(- %1 %2) (:x (dampedho 1 1 1 1.0 0.0 0.0 2.0 0.01)) (:y1 (statedampedho 1 1 1 1.0 0.0 0.0 2.0 0.01)))));error study with state2ode conversion
)