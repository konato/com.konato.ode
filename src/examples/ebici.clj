;; Ordinary Différential Equations Solvers - An Electric Bicycle Example 
;;
;; Methods based on: Gerald-Wheatley, 2003, Applied Numerical Analysis, Seventh Edition, 
;; Pearson Education, ISBN 0-321-13304-8.
;; Example based on: Urieli I., Winter 2002, An Approacj to Solving Ordinary Differentail Equations
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

(defn calc_torque [vel vbat imax rw winding kt n]
"Calc acceleration" 
;;(calc-torque vel vbat imax mt mw rw winding kt n) 
(let [
      tw (* (* n (/ kt winding)) (- vbat (* kt n (/ vel rw))))
      twmax (* n kt imax)
      ]
      (if (> tw twmax) twmax tw) ; continous current limiting controler
))

(defn calc_accel [vel vbat imax slope mt mw rw winding kt n cda density cr g]
"Calc acceleration"  
(let [
      jw (* mw (* rw rw)) ; [kg.m^2] moment of inertia of each wheel
      mi (+ mt (/ (* 2 jw) (* rw rw))) ; equivalent inertial mass of bike + rider
      mtw (calc_torque vel vbat imax rw winding kt n)
      drag (/ (* cda density (Math/pow vel 2)) 2)
      resist (* mt g (+ cr slope))
     ]
	  (/ (- (/ mtw rw) (+ drag resist)) mi)
))


(defn ebici [METHOD VBAT IMAX SLOPE DUR TS]
"Example of an electric bicycle parameters are :ode method, battery voltage, continous limiting controller current max, slope, duration of simulation, time sample period. "
;; (ebici rk4 36 32 0 30 0.1)
 (let [TSS (map #(+ (* TS %) 0.0) (take (inc (Math/round (/ DUR TS))) (iterate inc 0)));Produce a seq of digital interval
       ivals  {:dis 0.0, :vel 0.0} ; init vals of odes
       mt 100.0 ; [kg] mass of bike + rider
       mw 1.0 ; [kg] mass of each wheel
       rw 0.2115 ; [m] radius of wheel
       ;; motor/generator specs mg62
       winding 0.395 ; [ohms] resistance
       kt 0.62 ; [Nm/A] torque constant
       ;; wheel/motor gear ratio
       n 30/22 ; wheel torque/ motor torque
       ;; resistive forces
       cda 0.4 ; [m^2] coef of drag * frontal area
       density 1.18 ; [kg/m^3] air density
       cr 0.005 ; coefficient of rolling resistance
       g 9.807 ; [m/s^2] gravity 
       mtw (calc_torque (:vel ivals) VBAT IMAX rw winding kt n); calc_torque initial using current limiter
      ]
   (loop [
	  rzc  {:dis [(:dis ivals)], :vel [(:vel ivals)], :amp [(/ mtw (* n kt))], :pow [(* (:vel ivals) (/ mtw rw))]}
	  ti  TSS
	  ivs ivals
	  fps  {:dis #(:vel %), :vel  #(calc_accel (:vel %) VBAT IMAX SLOPE mt mw rw winding kt n cda density cr g)}
	 ]
     (if (empty? ti)
       (conj rzc {:t TSS})
       (let [
	     evs (euler fps ivs :t (first ti) TS) 
	     mtw (calc_torque (:vel evs) VBAT IMAX rw winding kt n)
	     amps (conj (:amp rzc) (/ mtw (* n kt)))
	     pows (conj (:pow rzc) (* (/ mtw rw) (:vel evs)))
	     ]
       (recur   
	  (merge (into {} (for [[k v] rzc] [k (into [] (conj v (k evs)))])) {:amp amps} {:pow pows })
	  (rest ti)
	  evs; the next initvals are the last
	  {:dis #(:vel %), :vel  #(calc_accel (:vel %) VBAT IMAX SLOPE mt mw rw winding kt n cda density cr g)}
	))))
))

(comment
  (load-file "com_konato_ode/src/examples/ebici.clj")
  (ebici rk4 36 32 0 30 0.1)
)