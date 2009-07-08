;; Ordinary Différential Equations Solvers - Pendulum Example 
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

(defn pendulum []
"Test the Runge-Kutta order 4 method for a classic pendulum dtheta/dt=omega,domega/dt=(-g/L)sin(theta),L=1m,theta(0)=20degrees,omega(0)=0"
 (let [L 1.0
       G 9.807
       PERIOD (* 2 (* Math/PI (Math/sqrt (/ L G))))
       N 200
       DT (/ (* 1.5 PERIOD) (- N 1))
       T0 0.0
       fps  {:theta #(:omega %), :omega #(* -1.0 (* (/ G L) (Math/sin (:theta %))))}
       outfps {:omega #(:omega %), :theta #(:theta %), :t #(:t %)}
      ]
      (odesolve rk4 fps {:omega 0.0, :theta (/ (* 20 Math/PI) 180)} :t T0 (Math/round (* DT N)) DT outfps)
))

(comment
  (load-file "com_konato_ode/src/examples/pendulum.clj")
  (pendulum)
)