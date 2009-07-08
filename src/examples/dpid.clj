;; Ordinary Différential Equations Solvers - Digital PID Controller Example 
;;
;; Methods based on: Gerald-Wheatley, 2003, Applied Numerical Analysis, Seventh Edition, 
;; Pearson Education, ISBN 0-321-13304-8.
;; Example based on: G.A. Korn,1989, Interactive Dynamic System Simulation, example 7.9"
;;
;;
;; Copyright (c) 2009 Stephane Rousseau (stephaner gmail.com). All rights reserved.  
;;
;; The use and distribution terms for this software are covered by the Eclipse Public License 1.0,
;; which can be found in the file epl-v10.html at the root of this distribution. 
;; By using this software in any fashion, you are agreeing to be bound by the terms of this license.
;; You must not remove this notice, or any other, from this software.

(use 'com.konato.ode)

(defn sat [val]
"Saturation function"
(cond (> val 1.0) 1
      (< val -1.0) -1
      :else val))

(defn dpid [DT T0 TF TS depvar]
"Example of a digital PID controller using an ODE solver, based on G.A. Korn,1989, Interactive Dynamic System Simulation, example 7.9"
;; (digitalpid 0.01 0 1 0.1 :t)
 (let [TSS (map #(+ (* TS %) T0) (take (inc (Math/round (/ (- TF T0) TS))) (iterate inc 0)));Produce a seq of digital interval
       r  0.75 ;step input
       A  0.99 ;torque limit
       x1 0.0  ;init val. of digital states vars
       x2 0.0  ;...
       torque 0 ; init val of torque
       ivals  {:c 0.0, :cdot 0.0} ; init vals of odes
       KP 2.0  ;Proportional coef
       KI 1.2  ;Integral coef
       KD 0.2  ;Derivative coef
       B0 (/ KD TS); preevals
       B1 (+ (* -1.0 KP) (* 0.5 KI TS) (* -2.0 B0))
       B2 (+ KP (* 0.5 KI TS) B0)
       outfps {:c #(:c %), :cdot #(:cdot %),:t #(:t %)}
      ]
   (loop [rzd  {:torque []}
	  rzc  (conj (into {} (for [[k v] ivals] [k []])) {:t []})
	  ts  TSS
	  fps  {:c #(:cdot %), :cdot #(+ (* 10.0 torque) (* -2.0 (:c %)) (* -3.0 (:cdot %)))}
	  ivs  ivals
	  error (- r (:c ivals)); digital initial eval
	  X1 x2 ;d state vars
	  X2 (+ error x2) ;...
	  y  (+ (* B0 x1) (* B1 X1) (* B2 X2))
	  torque (* A (sat (/ y A)))
	  x1 X1
	  x2 X2
	  
	 ]
     (if (= 1 (count ts))
       {depvar TSS, :torque (:torque rzd), :rzc rzc}
       (let [evs (odesolve rkf45 fps ivs :t (first ts) (second ts)
		      {:szmax 0.02, :szmin 0.001, :szstart DT :maxsteperr 0.00000001 :maxiter 50 :eps 0.000001 } outfps)
	     error (- r (last (:c evs))); digital initial eval
	     X1 x2 ;d state vars
	     X2 (+ error x2) ;...
	     y  (+ (* B0 x1) (* B1 X1) (* B2 X2))
	     torque (* A (sat (/ y A)))
	     x1 X1
	     x2 X2
	     ]
       (recur 
	  {:torque (merge (:torque rzd) torque)}  
	  (into {} (for [[k v] rzc] [k (into [] (concat v (k evs)))]))
	  (rest ts)
	  {:c #(:cdot %), :cdot #(+ (* 10.0 torque) (* -2.0 (:c %)) (* -3.0 (:cdot %)))}
	  (into {} (for [[k v] evs] [k (last v)])); the next initvals are the last 
	  error 
	  X1 
	  X2
	  y 
	  torque
	  x1
	  x2
	)))

)))

(comment
  (load-file "com_konato_ode/src/examples/dpid.clj")
  (dpid 0.01 0 3.5 0.1 :t)
)