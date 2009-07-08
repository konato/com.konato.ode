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

(defn errorstudy []
"Test the Runge-Kutta order 4 method for a classic damped harmonic oscilator, based on G.A. Kornm,1989, Interactive Dynamic System Simulation, example 1.4"
 (let [
       DT 0.02
       T0 0.0
       TF 5.2
       NN 261
       ts (take NN (interval-seq T0 DT))
       fps  {:x1 #(:x2 %), :x2 #(* -1.0 (:x1 %))}
       outfps {:x1 #(:x1 %), :x2 #(:x2 %)}
       rsim1 (odesolve eulerpc fps {:x1 1.0, :x2 0.0} :t T0 TF DT outfps)
       rsim2 (odesolve rk4 fps {:x1 1.0, :x2 0.0} :t T0 TF DT outfps)
       rsim3 (odesolve rkf45 fps {:x1 1.0, :x2 0.0} :t T0 TF DT outfps)
       rexact (map #(Math/cos %) ts) 
       error1 (map #(- %1 %2) rexact  (:x1 rsim1))
       error2 (map #(- %1 %2) rexact  (:x1 rsim2))
       error3 (map #(- %1 %2) rexact  (:x1 rsim3))
      ]
      { :t ts, :sim1 rsim1, :sim2 rsim2, :sim3 rsim3, :exact rexact, :error1 error1, :error2 error2, :error3 error3}
))

(comment
  (load-file "com_konato_ode/src/examples/errorstudy.clj")
  (errorstudy)
)