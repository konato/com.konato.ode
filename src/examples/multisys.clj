;; Ordinary Différential Equations Solvers - Multi Equations System Example
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

(defn odesolvenonstiff []
"Test the Runge-Kutta order 4 method for a non stiff system dy1/dt=y2*y3;dy2/dt=-y1*y3;dy3/dt=-0.51y1*y2, y1(0)=0, y2(0)=1, y3(0)=1, h=0.1,n=120"
 (let [
       fps  {:y1 #(* (:y2 %) (:y3 %)), :y2 #(* (* -1.0 (:y1 %)) (:y3 %)),  :y3 #(* (* -0.51 (:y1 %)) (:y2 %))}
       outfps {:y1 #(:y1 %), :y2 #(:y2 %), :y3 #(:y3 %), :x #(:x %)}
      ]
   (odesolve rk4 fps {:y1 0.0,:y2 1.0,:y3 1.0} :x 0.0 12 0.01 outfps)
))

(comment
  (load-file "com_konato_ode/src/examples/multisys.clj")
  (odesolvenonstiff)
)
