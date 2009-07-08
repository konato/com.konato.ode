;; Ordinary Différential Equations Solvers Tests library
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

(ns com.konato.ode.tests.test-ode
	(:use clojure.contrib.test-is)
	(:use com.konato.ode))

(deftest testkeepdecim
  (is (= 1.0 (keepdecim 0.99999999 1)))
  (is (= -1.0 (keepdecim -0.99999999 1)))
  (is (= 0.1111(keepdecim 0.11111111 4)))
  (is (= 0.1 (keepdecim 0.1111 1)))
)

(deftest testmaxvec
   (is (= 1.0 (max-vec '(0.034 1.0 -2.0))))
   (is (= 1.0 (max-vec [0.034 1.0 -2.0])))
)

(deftest testeuler
"Test the euler ode method for dy/dx=-2x-y, y(0)=-1, h=0.1"
 (let [fp  {:y #(+ (* (:x %) -2.0) (* (:y %) -1.0)), :y1 #(+ (* (:x %) -2.0) (* (:y %) -1.0))}
      ]
   (is (= {:y1 -0.8, :y -0.9} (euler fp {:y  -1.0, :y1 -0.9 } :x 0.0 0.1)))))

(deftest testeulerpc
"Test the euler predictor corector for ode method for dy/dx=-2x-y, y(0)=-1, h=0.1"
 (let [fp  {:y #(+ (* (:x %) -2.0) (* (:y %) -1.0)), :y1 #(+ (* (:x %) -2.0) (* (:y %) -1.0))} 
      ]
   (is (= {:y1 -0.8300000000000001, :y -0.915} (eulerpc fp {:y  -1.0, :y1 -0.915 } :x 0.0 0.1)))))

(deftest testrk4
"Test the Runge-Kurder 4 method for dy/dx=-2x-y, y(0)=-1, h=0.1"
 (let [fp  {:y #(+ (* (:x %) -2.0) (* (:y %) -1.0))}
      ]
   (is (= {:y -0.9145125} (rk4 fp {:y  -1.0} :x 0.0 0.1)))))

(deftest testevalks
"Test the subfunction used by rkp to evaluate the ks for the classical rk4 and the rk-fehlbergh method"
 (let [fp  {:y #(+ (* (:x %) -2.0) (* (:y %) -1.0))}
       park4 {:abs [1/6 [1 2 2 1]], :alphas [1/2 1/2 1], :betas [[1/2][0 1/2][0 0 1]]}
       parkrfk45  {:abs [1 [16/135 0 6656/12825 28561/56430 -9/50 2/55]], :alphas [1/4 3/8 12/13 1 1/2], :betas [[1/4][3/32 9/32][1932/2197 -7200/2197 7296/2197][439/216 -8 3680/513 -845/4104][-8/27 2 -3544/2565 1859/4104 -11/40]], :error [1/360 0 -128/4275 -2197/75240 1/50 2/55]}
       ]
  (is (= '(0.1 0.085 0.0858 0.0714) (map #(keepdecim % 4) (:y (evalks (:alphas park4) (:betas park4)  fp {:y  -1.0} :x 0.0 0.1)))))
  (is (= '(0.1 0.0925 0.0889609 0.0735157 0.0713736 0.0853872) (map #(keepdecim % 7) (:y (evalks (:alphas parkrfk45) (:betas parkrfk45)  fp {:y  -1.0} :x 0.0 0.1)))))
))

(deftest testrkp
"Test the flexible Runge-Kutta using parameters for dy/dx=-2x-y, y(0)=-1, h=0.1 
with the classical rk4 and the rk-fehlbergh method"
 (let [fp  {:y #(+ (* (:x %) -2.0) (* (:y %) -1.0))}
       park4 {:abs [1/6 [1 2 2 1]], :alphas [1/2 1/2 1], :betas [[1/2][0 1/2][0 0 1]]}
       parkrfk45  {:abs [1 [16/135 0 6656/12825 28561/56430 -9/50 2/55]], :alphas [1/4 3/8 12/13 1 1/2], :betas [[1/4][3/32 9/32][1932/2197 -7200/2197 7296/2197][439/216 -8 3680/513 -845/4104][-8/27 2 -3544/2565 1859/4104 -11/40]], :error [1/360 0 -128/4275 -2197/75240 1/50 2/55]}
       ]
   (is (= {:y -0.9145125} (rkp park4 fp {:y  -1.0} :x 0.0 0.1)))
   (is (= -0.914512251 (keepdecim (:y (rkp parkrfk45 fp {:y  -1.0} :x 0.0 0.1)) 9)))
))


(deftest testodesolve
"Test the odesolve for dy/dx=-2x-y, y(0)=-1, h=0.1 with all the defined methods"
 (let [fp  {:y #(+ (* (:x %) -2.0) (* (:y %) -1.0))}
       outfps {:y #(:y %), :x #(:x %)}]
   (is (= '(-1.0 -0.9 -0.83 -0.787 -0.7683 -0.7715) 
	  (map #(keepdecim % 4) (:y (odesolve euler fp {:y  -1.0} :x 0.0 0.5 0.1 outfps)))))
   (is (= '(-1.0 -0.915 -0.8571 -0.8237 -0.8124 -0.8212) 
	  (map #(keepdecim % 4) (:y (odesolve eulerpc fp {:y  -1.0} :x 0.0 0.5 0.1 outfps)))))
   (is (= '(-1.0 -0.91451 -0.85619 -0.82246 -0.81096 -0.81959 -0.84644) 
	  (map #(keepdecim % 5) (:y (odesolve rk4 fp {:y  -1.0} :x 0.0 0.6 0.1 outfps)))))
   (is (= '(-1.0 -0.91451 -0.85619 -0.82246 -0.81096 -0.81959 -0.84644) 
	  (map #(keepdecim % 5) (:y (odesolve rkp4 fp {:y  -1.0} :x 0.0 0.6 0.1 outfps)))))
   (is (= '(-1.0 -0.914512251 -0.856192254 -0.822454655 -0.81096013 -0.81959197 -0.846434899) 
	  (map #(keepdecim % 9) (:y (odesolve rkf45 fp {:y  -1.0} :x 0.0 0.6 0.1 outfps)))))
   (is (= '(-1.0 -0.914512251 -0.856192254 -0.822454655 -0.81096013 -0.81959197 -0.846434899) 
	  (map #(keepdecim % 9) (:y (odesolve rkf45 fp {:y  -1.0} :x 0.0 0.6 
					      {:szmax 0.250, :szmin 0.001, :szstart 0.1 :maxsteperr 0.0001 :maxiter 10 :eps 0.000001 } outfps)))))
))

(deftest teststep
"Test the step response utilitarian function"
(let [tes #(* 10 (step (:t %)))
     ]
(is (= 0 (tes {:t -1}))) ; test the step response 
(is (= 10 (tes {:t 0})))
(is (= 10 (tes {:t 1})))
))

(deftest teststate2ode
"Test the transformation of a State-Space representation to an ODE system."
(let [
      A [[0 1 0][0 0 1][-24 -26 -9]]  
      B [[0][0][1]]
      C [[2 7 1]]
      D [[0]]
      U [#(* 1 (step (:t %)))]
      res (state2ode A B C D U)
     ]
(is (= 4 ((:x2 (:dotfps res)){:x3 4 :t 1})))
(is (= 1 ((:x1 (:dotfps res)){:x2 1 :t 1})))
(is (= -102 ((:x3 (:dotfps res)) {:x1 1 :x2 2 :x3 3 :t 1}))); -24-52-27+1=-102
(is (= 19   ((:y1 (:outfps res)) {:x1 1 :x2 2 :x3 3 :t 1}))); 2+14+3=19
))

(deftest testtransfer2ode
"Test the transformation of a Transfert function representation to an ODE system."
(let [
      NUM [1 7 2]
      DEN [1 9 26 24]
      inp [#(* 1 (step (:t %)))]
      res (transfer2ode NUM DEN inp)
     ]
(is (= 4 ((:x2 (:dotfps res)){:x3 4 :t 1})))
(is (= 1 ((:x1 (:dotfps res)){:x2 1 :t 1})))
(is (= -102 ((:x3 (:dotfps res)) {:x1 1 :x2 2 :x3 3 :t 1}))); -24-52-27+1=-102
(is (= 19   ((:y1 (:outfps res)) {:x1 1 :x2 2 :x3 3 :t 1}))); 2+14+3=19
))


(comment
(load-file "com_konato_ode/src/com/konato/ode.clj")
(load-file "com_konato_ode/src/com/konato/ode/tests/test_ode.clj")
(in-ns 'com.konato.ode.tests.test-ode)
(run-tests 'com.konato.ode.tests.test-ode)
)

