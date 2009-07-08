;; Ordinary Différential Equations Solvers
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

(ns com.konato.ode)

;; Auxiliary function
(defn interval-seq
"Produce a sequence of number beginning at v0 incremented by stepsz"
[v0 stepsz] 
(map #(+ v0 (* stepsz %))  (iterate inc 0)))
;;(take 10 (interval-seq 0.0 0.1))

;;Auxiliary function
(defn max-vec
  "Returns the max x of a vector."
[vec] (reduce (fn [x y] (if (> x y) x y)) (first  
vec) (rest vec))) 

;;Auxiliary function
(defn keepdecim 
"Keep a number of decimals."
[val digits]  
(/ (Math/round (* val (Math/pow 10.0 digits))) (Math/pow 10.0 digits)))

;;
;; Inputs common functions
;;
(defn step
"Return 1 if val >= 0 else 0"
[val]
(if (>= val 0.0) 1.0 0.0))

(defn dirac
"Return 1 if val = 0 else 0"
[val]
(if (= val 0.0) 1.0 0.0))

(defn ramp
"Return val, could be simpler use the val or macro"
[val]
val)

;;
;; ODE Solvers implementations
;;
(defn euler
"Evaluate an ODE system using Euler method."
[fps initvals depvar depv0 stepsz]
  (into {} (for [[k fp] fps] [k (+ (k initvals) (* stepsz (fp  (merge initvals {depvar depv0}))))])))


(defn eulerpc
"Evaluate an ODE system using an improved Euler method using predictor and corrector."
[fps initvals depvar depv0 stepsz]
  (let [dys  (into {} (for [[k fp] fps] [k (fp  (merge initvals {depvar depv0}))]));  
        y1p  (into {} (for [[k dy] dys] [k (+ (k initvals) (* stepsz dy))]));y(n+1),p   
        dy1p (into {} (for [[k fp] fps] [k (fp (merge y1p {depvar (+ depv0 stepsz)}))]));dy(n+1),p  
       ]
     (into {} (for [[k _] fps] [k (+ (k initvals) (* stepsz (/ (+ (k dys) (k dy1p)) 2)))]))
))

(defn rk4
"Evaluate an ODE system using an Runge-Kutta order 4."
[fps initvals depvar depv0 stepsz]
(let [k1 (into {} (for [[k fp] fps] [k (* stepsz (fp (merge initvals {depvar depv0})))]))
      y1  (into {} (for [[k fp] fps] [k (+ (k initvals) (* 0.5 (k k1)))]))
      k2 (into {} (for [[k fp] fps] [k (* stepsz (fp (merge y1 {depvar (+ depv0 (* 0.5 stepsz))})))]))
      y2  (into {} (for [[k fp] fps] [k (+ (k initvals) (* 0.5 (k k2)))]))
      k3 (into {} (for [[k fp] fps] [k (* stepsz (fp (merge y2 {depvar (+ depv0 (* 0.5 stepsz))})))]))
      y3  (into {} (for [[k fp] fps] [k (+ (k initvals) (k k3))]))
      k4 (into {} (for [[k fp] fps] [k (* stepsz (fp (merge y3 {depvar (+ depv0 stepsz)})))]))  
      ]
  (into {} (for [[k _] fps] [k (+ (k initvals) (* 1/6 (+ (k k1) (* 2 (k k2)) (* 2 (k k3)) (k k4))))]))
))


(defn evalks
"Recursive evaluation of the ks for the Runge-Kutta parametrized method (rkp)."
[alps bets fps initvals depvar depv0 stepsz]
(loop [alps alps
       bets bets 
       ks  (into {} (for [[k fp] fps] [k (vector (* stepsz (fp (merge initvals {depvar depv0}))))]))
      ]
     (if (empty? alps)
       ks ; finished, we return the ks to rpk
       (let [
	     yn (into {} (for [[k fp] fps] [k (+ (k initvals) (reduce + (map * (first bets) (k ks))))])) ;; compute yn+
	     kn (into {} (for [[k fp] fps] [k (* stepsz (fp (merge yn {depvar (+ depv0 (* stepsz (first alps)))})))])) ;; compute kn+
	     ks (into {} (for [[k rz] ks] [k (concat rz (vector (k kn)))])) ; add k to list of ks
	    ]
	     (recur (rest alps) (rest bets) ks )))
))

(defn rkp
"Solve an ODE system using paramatrized Runge-Kutta method, this allows to implements various orders and methods. If :error par a furnished, the function will compute the error for each independant var and will return the maximum error as :err tag into return val."
[par fps initvals depvar depv0 stepsz]
(let [ks (evalks (:alphas par) (:betas par) fps initvals depvar depv0 stepsz);eval the ks
      rz  (into {} (for [[k fp] fps] [k (+ (k initvals) (* (first (:abs par)) (reduce + (map * (second (:abs par)) (k ks)))))]));for each ode eval yn+1 = yn + abfact (ak1+bk2+...)
      err (if (:error par) (into [] (for [[k fp] fps] (Math/abs (reduce + (map * (:error par) (k ks)))))) 0.0);eval the error if :error param is present
     ]
  (if (:error par) (assoc rz :err (max-vec err)) rz) ;return the rz add the maximun error if error params allowed evaluation
))


(defn make-rkp [pars]
"Make a new Runge-Kutta method using parameters." 
  (let [par pars]
    (fn  [fps initvals depvar depv0 stepsz] (rkp par fps initvals depvar depv0 stepsz))))

(def rkp4 (make-rkp {:abs [1/6 [1 2 2 1]], :alphas [1/2 1/2 1], :betas [[1/2][0 1/2][0 0 1]]}))

(def rkf45 (make-rkp {:abs [1 [16/135 0 6656/12825 28561/56430 -9/50 2/55]], :alphas [1/4 3/8 12/13 1 1/2], :betas [[1/4][3/32 9/32][1932/2197 -7200/2197 7296/2197][439/216 -8 3680/513 -845/4104][-8/27 2 -3544/2565 1859/4104 -11/40]], :error [1/360 0 -128/4275 -2197/75240 1/50 2/55]}))

(defmulti odesolve
"Solve an ODE system using either an fixed step or a variable step solver. Map outfps allows selection of outputs, apply formulas."
(fn [method fps initvals depvar depv0 depvf stepsz outfps] [(if (map? stepsz) :vstep :fstep)]))

;; Ordinary differential equations solver fixed step
(defmethod odesolve [:fstep] [method fps initvals depvar depv0 depvf stepsz outfps]
"Solve an ODE system using a fixed method"
   (loop [tvals (take (Math/round (Math/abs (/ (- depvf depv0) stepsz))) (interval-seq depv0 stepsz)) 
          rz (into {} (for [[k fp] outfps] [k (vector (fp (merge initvals {depvar depv0})))]));eval the output fcts with initvals and depvar v0
	  ;rz (conj (into {} (for [[k val] initvals] [k (vector val)])) {depvar (vector depv0)}) ;; initvals are included in the result 
          ivals initvals]
      (if (empty? tvals)
          rz
          (let [vals1 (dissoc (method fps ivals depvar (first tvals) stepsz) :err)];eval and remove :err key not needed
	    (recur 
	     (rest tvals)
	     (into {} (for [[k fp] outfps] [k (merge (k rz) (fp (merge vals1 {depvar (+ (first tvals) stepsz)})))]))
	     ;(conj (into {} (for [[k val] vals1] [k (merge (k rz) val)])) {depvar (merge (depvar rz) (+ (first tvals) stepsz))}) 
	     vals1))))
)

;; required for variable step odesolve
(defn- adjust-stepsz
"Adjust the step size using a popular method described in http://math.fullerton.edu/mathews/n2003/RungeKuttaFehlbergMod.html."
[tol err cursz minsz maxsz]
(let [
      s (if (= err 0) 0 (Math/pow (/ (* tol cursz) (* 2 err)) 1/4));eval the new s
      ss (cond (< s 0.1) 0.1
	       (> s 4.0) 4.0
	       :else s)
      sz (* ss cursz)
      ]
  (cond (< sz minsz) minsz
	(> sz maxsz) maxsz
	:else sz))
)

;; Ordinary differential equations solver using a variable step
(defmethod odesolve [:vstep] [method fps initvals depvar depv0 depvf stepsz outfps]
"Solve an ODE system usign a variable method, method should give a :error ks series to calcule step size. 
Map outfps allows selection of outputs, apply formulas."
(loop [sz (:szstart stepsz)
       depvn depv0
       rz (into {} (for [[k fp] outfps] [k (vector (fp (merge initvals {depvar depv0})))]));eval the output fcts with initvals and depvar v0
       ivals initvals
       iter 1
       ]
  (if (or (>= (+ depvn (:eps stepsz)) depvf) (>= iter (:maxiter stepsz))) 
    rz ;return the results
    (let [ ev (method fps ivals depvar depvn sz)
	   vals1 (dissoc ev :err)
	   err (:err ev)
	 ]
      (if (and (>= err (:maxsteperr stepsz)) (>= sz (* 2 (:szmin stepsz))))
	(let [newsz (adjust-stepsz (:maxsteperr stepsz) err sz (:szmin stepsz) (:szmax stepsz));eval the new step size
	     ]
	  (recur 
	   (if (> (+ depvn newsz) depvf) (- depvf depvn) newsz) 
	   depvn 
	   rz 
	   ivals 
	   iter)); correct if > end of interval
	(recur 
	 (if (> (+ depvn sz) depvf) (- depvf depvn) sz) ; correct if > than end of interval
	 (if (> (+ depvn sz) depvf) depvf (+ depvn sz)) ; correct if > than end of interval
         (into {} (for [[k fp] outfps] [k (merge (k rz) (fp (merge vals1 {depvar (+ depvn sz)})))]))
	 vals1
	 (inc iter)))
    ))
))

;--------- state2ode section ------------

(defmacro xs [pre s]
"Return prefixed map of none zero factors (xs pre [0.0 -1.0 -2.0]) == {:x2 -1.0, :x3 -2.0}"
`(let [ xeq1#  (into {} (map (fn [x# y# z#] [(keyword (str x# y#)) z#]) (repeat (count ~s) ~pre) (take (count ~s) (iterate inc 1)) ~s))
   ]
 (into {} (for [[k# m#] xeq1# :when (not= m# 0.0)] [k# m#]))));remove zero factor
;;(xs "x" [0.0 -1.0 -2.0]) == {:x2 -1.0, :x3 -2.0}

(defmacro toterms [fs ivals]
"Return a partial linear equation using factors in the form {:x2 -1.0, :x3 -2.0}"
`(reduce + (map (fn [fa#] (* ((key fa#) ~ivals) (val fa#))) ~fs)))
;;toterms fonctionne reste avoir avec fct
;(def f1 (fn [ival] (toterms {:x2 -1.0 :x1 -3.0} ival)))
;(f1 {:x2 2.0, :x1 2.0})

(defmacro frominput [ts us ivals]
"Return an equation from the input dependant part."
`(reduce + (map (fn [t# u#] (* t# (cond (number? u#) u# :else (u# ~ivals)))) ~ts ~us)))
;(def f2 (fn [ivals] (frominput  [1 1] [10 #(* 5 (step (:t %)))] ivals)))
;(f2 {:t 1.1}) == 15

(defmacro tofp  [pre s t u ivals]
"Return a linear equation"
`(let [pre# (xs ~pre ~s)]
  (+ (toterms pre# ~ivals) (frominput ~t ~u ~ivals))))
;;(def f3 (fn [ival] (tofp "x" [-1.0 0.0 -2.0] [1 1] [10 #(* 5 (step (:t %)))] ival)))
;;(f3 {:x1 (* (step 0.001)), :x3 1, :t 0.5})

(defn- tofps [prefixres prefix S T U] 
"Convert a system represented by prefixres = Sprefix + TU"
(into {} (map #(hash-map (keyword (str prefixres %1)) (fn [ival] (tofp prefix %2 %3 U ival))) (take (count S) (iterate inc 1)) S T)))

(defn state2ode [A B C D U]
"Convert a State Space representation to an ODE system of equations who can be solved using odesolve"
{:dotfps (tofps "x" "x" A B U), :outfps (tofps "y" "x" C D U)}
)

;-------- transfer2ode section --------
(defmacro singledot [i]
`(keyword (str "x" ~i)))

(defn tf2odedot1 [DEN] 
"Convert part 1 of transfer into a dot eqns "
(into {} (map #(hash-map (keyword (str "x" %1)) (fn [ival] ((singledot (+ 1 %1)) ival))) (take (count (rest (rest DEN))) (iterate inc 1)))))

(defn tf2odedotlast [DEN U]
"Convert to the last row who is the same a DEN minus first item and multiplied by -1"
(let [D (reverse (map #(* -1.0 %) (rest DEN)))] 
  {(keyword (str "x" (int (- (count DEN) 1.0)))) (fn [ival] (tofp "x" D [1] U ival))}
))

(defn tf2odeout [NUM U]
"Convert to the last row who is the same a DEN minus first item and multiplied by -1"
(let [N (reverse NUM)] 
  {:y1 (fn [ival] (tofp "x" N [0] U ival))}
))

(defn transfer2ode [NUM  DEN inp]
"Convert a Transfer function representation to an ODE system of equations who can be solved using odesolve"
; validate length of NUM and DEN
(if (not (= 1 (- (count DEN) (count NUM))))
  (throw (IllegalArgumentException. "Number of coefs in numerator should be 1 less than number of coefs in denominator"))
  ; make first coefs of NUM to unity by dividing all NUM and DEN
  (let [N (map #(/ % (first DEN)) NUM)
	D (map #(/ % (first DEN)) DEN)]
    {:dotfps (into (tf2odedot1 D) (tf2odedotlast D inp)), :outfps (tf2odeout N inp)}
)))

(comment
  (load-file "com_konato_ode/src/com/konato/ode.clj")
  (use 'com.konato.ode)
)

