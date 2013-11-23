{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE MultiParamTypeClasses      #-}

module Hammer.Texture.SH.Rotation where

{--
( rotM
       , rotMReal
       , rotMZ
       , rotMZReal
       , multLStack
       )
--}

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Packed         as M

import           Data.Complex
import           Foreign.Storable      (Storable)
import           Data.Vector           (Vector)
import           Data.Packed.Matrix    ((><), (@@>))
import           Numeric.LinearAlgebra ((<>))
import           Numeric.Container     (Product, add)

-- import           Hammer.Math.Algebra
import           Hammer.Texture.SH.Pyramid
import           Hammer.Texture.SH.HyperSphere
import           Hammer.Texture.SH.SupportFunctions
import           Hammer.Texture.SphericalHarmonics
import           Hammer.Texture.Orientation

import           Debug.Trace
dbg s x = trace (show s L.++ show x) x


-- ==================================== HSH rotation =====================================

-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotHSHMVec :: Int -> SO3 -> SO3 -> Vector (M.Matrix (Complex Double))
rotHSHMVec n r p = let
  lmax = n `quot` 2
  us   = rotSHMVec id        lmax r
  us'  = rotSHMVec conjugate lmax p
  in V.izipWith (\l u u' -> rotHSHM' u u' (2*l)) us us'

-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotHSHM :: Int -> SO3 -> SO3 -> M.Matrix (Complex Double)
rotHSHM n r p = let
  l  = n `quot` 2
  u  = (rotSHMVec id        l r) V.! l
  u' = (rotSHMVec conjugate l p) V.! l
  in rotHSHM' u u' n

-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotHSHM' :: M.Matrix (Complex Double) -> M.Matrix (Complex Double) -> Int -> M.Matrix (Complex Double)
rotHSHM' u u' n = let
  l  = n `quot` 2

  calcAtM (row, col) = let
    p' = ps V.! row
    p  = ps V.! col
    us = [l, l-1 .. (-l)]
    in sum [ func p' p (m', m) (m''', m'')
           | m <- us, m' <- us, m'' <- us, m''' <- us
           ]

  func (lb', mu') (lb, mu) (m', m) (m''', m'') = let
    -- calcClebshGordan j1 m1 j2 m2 j m
    c1 = calcClebshGordan (l, m'' ) (l, m') (lb', mu')
    c2 = calcClebshGordan (l, m''') (l, m ) (lb , mu )
    z  = u  @@> (l - m'  , l - m  )
    z' = u' @@> (l - m''', l - m'')
    in ((c1 * c2) :+ 0) * z * z'

  ps = V.fromList [ (lb, mu) | lb <- [0 .. 2*l], mu <- [lb, lb-1 .. (-lb)] ]
  ms = [ (row, col) | col <- [0 .. s-1], row <- [0 .. s-1] ]
  s  = (2*l+1) * (2*l+1)
  in (s >< s) $ map calcAtM ms


rotHSHRealMVec :: Int -> SO3 -> SO3 -> Vector (M.Matrix Double)
rotHSHRealMVec n r p = let
  l  = n `quot` 2
  vu = rotHSHMVec (2*l) r p
  vt = realizeHSHVec (2*l)
  func i u = let
    t  = M.diagBlock $ V.toList $ V.take (2*i+1) vt
    t' = M.trans $ M.mapMatrix conjugate t
    in M.mapMatrix realPart $ t' <> u <> t
  in V.imap func vu

-- ===================================== SH rotation =====================================

-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotSHMVec :: (Complex Double -> Complex Double) -> Int -> SO3 ->
             Vector (M.Matrix (Complex Double))
rotSHMVec foo llen r = let
  z = genSHFunc (2*llen) r

  func l (m', m) (L lb, MF mu) = let
    ld  = fromIntegral l
    lbd = fromIntegral lb
    k1  = (pi / (2 * ld + 1)) * sqrt (2 * (2 * lbd + 1))
    -- calcClebshGordan j1 m1 j2 m2 j m
    c1 = calcClebshGordan (l, m') (lb, mu) (l, m)
    z1 = foo $ z %! (N (2*l), L lb, MF mu)
    in ((k1 * c1) :+ 0) * z1

  calcAtM l (m', m) = let
    ps = genLinSeq (2*llen)
    in L.foldl' (\acc p -> acc + func l (m', m) p) 0 ps

  buildM l = let
    ms = [l, l-1 .. (-l)]
    ps = [(i, j) | i <- ms, j <- ms]
    in ((2*l+1) >< (2*l+1)) $ map (calcAtM l) ps
  in V.generate (llen+1) buildM

rotSHRealMVec :: (Complex Double -> Complex Double) -> Int -> SO3 ->
                 Vector (M.Matrix Double)
rotSHRealMVec foo l r = let
  vu = rotSHMVec foo l r
  func i u = let
    c  = realizeSH i
    c' = M.trans $ M.mapMatrix conjugate c
    in M.mapMatrix realPart $ c' <> u <> c
  in V.imap func vu

-- ======================================== Map L ========================================

multLStack :: (U.Unbox a, Storable a, Product a)=> (Int -> M.Matrix a) ->
             Pyramid (L, MF) a -> Pyramid (L, MF) a
multLStack getMs pyC = let
  l       = maxStack pyC
  ls      = [0 .. l]
  toVec   = U.fromList . M.toList
  func li = (getMs li) <> (getLSlice pyC (L li))
  pyra    = U.concat $ map (toVec . func) ls
  in pyC { pyramid = pyra }

getLSlice :: (U.Unbox a, Storable a)=> Pyramid (L, MF) a -> L -> M.Vector a
getLSlice pyC l@(L li) = toVec $ slicer $ pyramid pyC
  where
    toVec  = M.fromList . U.toList
    slicer = U.slice (getKeyPos (l, MF li)) (2 * li + 1)

-- ======================================== Map N ========================================

multNStack :: (U.Unbox a, Storable a, Product a)=> (Int -> M.Matrix a) ->
             Pyramid (N, L, MF) a -> Pyramid (N, L, MF) a
multNStack getMs pyC = let
  n       = dbg "max: " $ maxStack pyC
  ls      = [0 .. (n `quot` 2)]
  toVec   = U.fromList . M.toList
  func li = (getMs li) <> (getNSlice pyC (N $ 2*li))
  pyra    = U.concat $ map (toVec . func) ls
  in pyC { pyramid = pyra }

getNSlice :: (U.Unbox a, Storable a)=> Pyramid (N, L, MF) a -> N -> M.Vector a
getNSlice pyC n = toVec $ slicer $ pyramid pyC
  where
    l      = (unN n) `quot` 2
    toVec  = M.fromList . U.toList
    slicer = U.slice (getKeyPos (N (2*l), L 0, MF 0)) ((2*l+1)*(2*l+1))

-- ==================================== Map L over N =====================================

multNLStack :: (U.Unbox a, Storable a, Product a)=> (Int -> M.Matrix a) ->
             Pyramid (N, L, MF) a -> Pyramid (N, L, MF) a
multNLStack getMs pyC = let
  n       = maxStack pyC
  ns      = [0, 2 .. n]
  func ni = multLStack getMs (getNLSlice pyC (N ni))
  pyra    = U.concat $ map (pyramid . func) ns
  in pyC { pyramid = pyra }

getNLSlice :: (U.Unbox a, Storable a)=> Pyramid (N, L, MF) a -> N -> Pyramid (L, MF) a
getNLSlice pyC nt = generatePyramid func ni
  where
    vec = slicer $ pyramid pyC
    n@(N ni) = if even (unN nt) then nt else nt-1
    func lmf = vec U.! getKeyPos (lmf :: (L, MF))
    slicer = U.slice (getKeyPos (n, L 0, MF 0)) ((ni + 1)^(2 :: Int))

-- ============================= Conversion to real matrix ===============================

realizeSH :: Int -> M.Matrix (Complex Double)
realizeSH l = let
  calcAtM (i, j)
    | i == 0 && j == 0 = 1 :+ 0
    | abs i /= abs j   = 0 :+ 0
    | j < 0 && i < 0   = 0 :+ (-s / (sqrt 2))
    | j < 0 && i > 0   = 0 :+ (1 / (sqrt 2))
    | j > 0 && i < 0   = (s / (sqrt 2)) :+ 0
    | otherwise        = (1 / (sqrt 2)) :+ 0     --  j > 0 && i > 0
    where s = if even i then 1 else (-1)
  in ((2*l+1) >< (2*l+1)) $ map calcAtM [(i, j) | i <- [-l..l], j <- [-l..l]]

realizeSHVec :: Int -> Vector (M.Matrix (Complex Double))
realizeSHVec l = V.generate (l+1) realizeSH

realizeHSH :: Int -> M.Matrix (Complex Double)
realizeHSH l = let
  calcAtM (mi, mj)
    | mi == 0 && mj == 0 = (i l)
    | abs mi /= abs mj   = 0 :+ 0
    | mj < 0 && mi < 0   = (i l) * ((s / (sqrt 2)) :+ 0)
    | mj > 0 && mi < 0   = (i l) * ((1 / (sqrt 2)) :+ 0)
    | mj < 0 && mi > 0   = (i $ l-1) * ((-1 / (sqrt 2)) :+ 0)
    | otherwise          = (i $ l-1) * ((s  / (sqrt 2)) :+ 0) --  j > 0 && i > 0
    where
      i = powerComplex
      s = if even mi then 1 else (-1)
  ms = [l, l-1 .. (-l)]
  ps = [(mi, mj) | mi <- ms, mj <- ms]
  in ((2*l+1) >< (2*l+1)) $ map calcAtM ps

realizeHSHVec :: Int -> Vector (M.Matrix (Complex Double))
realizeHSHVec l = V.generate (l+1) realizeHSH

-- ================================ Conversion to real Pyramid ===========================

fromComplexNLMF :: Pyramid (N, L, MF) (Complex Double) -> Pyramid (N, L, MF) Double
fromComplexNLMF pyC = let
  n       = maxStack pyC
  ns      = [0, 2 .. n]
  func ni = fromComplexLMF2 (getNLSlice pyC (N ni))
  pyra    = U.concat $ map (pyramid . func) ns
  in pyC { pyramid = pyra }

fromComplexLMF :: Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) Double
fromComplexLMF pyC = generatePyramid (realPart . func) lmax
  where
    lmax = maxStack pyC
    kr = (1/(sqrt 2)) :+ 0
    ki = 0 :+ (-1/(sqrt 2))
    func (l, mf)
      | mf > 0 = kr * (pyC %! (l, mf)  + s*(pyC %! (l, -mf)))
      | mf < 0 = ki * (pyC %! (l, -mf) - s*(pyC %! (l,  mf)))
      | otherwise = pyC %! (l, mf)
      where s = if even (unMF mf) then 1 else -1

fromComplexLMF2 :: Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) Double
fromComplexLMF2 pyC = generatePyramid (k . realPart . func) lmax
  where
    lmax = maxStack pyC
    k = (/ (sqrt 2))
    i = powerComplex
    func (l, mf)
      | mf > 0 = i (unL l)     * (s*pyC %! (l, mf) + (pyC %! (l, -mf)))
      | mf < 0 = i (unL $ l-1) * (s*pyC %! (l, mf) - (pyC %! (l, -mf)))
      | otherwise = pyC %! (l, mf)
      where s = if even (unMF mf) then 1 else -1

-- ======================================== Test ========================================

testRotSH :: IO ()
testRotSH = let
  foo  = (rotSHRealMVec id 10 (SO3 (pi/4) (pi/4) 0) V.!)
  func = multLStack foo
  in testSH func

testRotSymm :: IO ()
testRotSymm = let
  rot1 = rotHSHMVec 6 (SO3 0 0 0) (SO3 (pi/2) (pi/2) 0)
  rot2 = rotHSHMVec 6 (SO3 0 0 0) (SO3 (-pi/2) (pi/2) 0)
  rot3 = rotHSHMVec 6 (SO3 0 0 0) (SO3 (0)    (0)    0)
  rot  = L.foldl' (V.zipWith add) rot1 [rot2, rot3]
  func = fromComplexNLMF . multNStack (rot V.!)
  in testHSH_C func

testRotHSH2 :: IO ()
testRotHSH2 = let
  foo  = (rotHSHMVec 6 (SO3 0 0 0) (SO3 (-pi/2) (pi/2) 0) V.!)
  func = fromComplexNLMF . multNStack foo
  in testHSH_C func

testRotHSH :: IO ()
testRotHSH = let
  foo  = (rotSHMVec id 6 (SO3 (pi/2) 0 0) V.!)
  func = fromComplexNLMF . multNLStack foo
  in testHSH_C func

testHSH_C :: (Pyramid (N,L,MF) (Complex Double) -> Pyramid (N,L,MF) Double) -> IO ()
testHSH_C func = let
  s1 = SO3 (pi/2) (pi/2) (pi/2)
  s2 = SO3 (pi/4) (pi/4) 0
  xs :: Vector (Complex Double, SO3)
  xs  = V.fromList [(10 :+ 0, s1), (10 :+ 0, s2)]
  c   = findSHCoefWeight 6 xs
  vtk = renderSO3SolidVTK (evalSH $ func c)

  rot = so3ToQuaternion $ SO3 (pi/2) (pi/2) (0)
  q1  = so3ToQuaternion s1
  q2  = so3ToQuaternion s2
  so  = V.fromList [s1, s2]
  sp  = V.fromList [ quaternionToSO3 $ q1 #<= rot
                   , quaternionToSO3 $ q2 #<= rot ]
  sa  = V.fromList [ quaternionToSO3 $ activeRot q1 rot
                   , quaternionToSO3 $ activeRot q2 rot ]
  pso = renderSO3PointsVTK so
  psp = renderSO3PointsVTK sp
  psa = renderSO3PointsVTK sa
  in do
    writeQuater "Rot-SHS-original-points" pso
    writeQuater "Rot-SHS-passive-points"  psp
    writeQuater "Rot-SHS-active-points"   psa
    writeQuater "Rot-SHS-test" vtk

testSHRotC :: Pyramid (L, MF) (Complex Double)
testSHRotC = let
  xs :: Vector (Complex Double, SO2)
  xs = V.fromList [(10 :+ 0, (SO2 0 0)), (10 :+ 0, SO2 (pi/4) (pi/4)), (10 :+ 0, SO2 (pi/2) (pi/2))]
  in findSHCoefWeight 10 xs

testSHRot :: Pyramid (L, MF) (Double)
testSHRot = let
  xs :: Vector (Double, SO2)
  xs = V.fromList [(10, SO2 0 0), (10, SO2 (pi/4) (pi/4)), (10, SO2 (pi/2) (pi/2))]
  in findSHCoefWeight 10 xs

-- ======================================== Trash ========================================

-- TODO  move to Orientation and add test case
-- | Applies an active rotation to an other rotation. It means that the second rotation is
-- applied on the reference frame instead of the rotated frame. It is used for applying
-- sample symmetry.
activeRot :: Quaternion -> Quaternion -> Quaternion
activeRot q1 q2 = let
  (q20, q2v) = splitQuaternion q2
  q2vAct     = activeVecRotation q2v (invert q1)
  in q1 #<= unsafeMergeQuaternion (q20, q2vAct)

testRowCol0 :: N -> SO3 -> SO3 -> Double
testRowCol0 n r p = let
  rp  = quaternionToSO3 ((so3ToQuaternion r) #<= (so3ToQuaternion p))
  pr  = quaternionToSO3 ((so3ToQuaternion p) #<= (so3ToQuaternion r))
  vC1 = rotHSHCol0 n pr
  vC2 = head $ M.toColumns $ rotHSHM (unN n) r p
  vR1 = rotHSHRow0 n rp
  vR2 = head $ M.toRows    $ rotHSHM (unN n) r p
  tc  = V.imap (\i x -> abs $ magnitude $ x - (vC2 M.@> i) ) vC1
  tr  = V.imap (\i x -> abs $ magnitude $ x - (vR2 M.@> i) ) vR1
  in max (V.maximum tc) (V.maximum tr)

rotHSHRow0 :: N -> SO3 -> Vector (Complex Double)
rotHSHRow0 n r = let
  l = (unN n) `quot` 2
  z = genSHFunc (2*l) r

  func (lb, mu) = let
    k1 = sqrt 2 * pi / (fromIntegral $ 2*l+1)
    z1 = conjugate $ z %! (N (2*l), L lb, MF mu)
    in (k1 :+ 0) * z1

  ps = V.fromList [ (lb, mu) | lb <- [0 .. 2*l], mu <- [lb, lb-1 .. (-lb)] ]

  in V.map func ps

rotHSHCol0 :: N -> SO3 -> Vector (Complex Double)
rotHSHCol0 n r = let
  l = (unN n) `quot` 2
  z = genSHFunc (2*l) r

  func (lb, mu) = let
    k1 = (-1)^lb * sqrt 2 * pi / (fromIntegral $ 2*l+1)
    z1 = z %! (N (2*l), L lb, MF mu)
    in (k1 :+ 0) * z1

  ps = V.fromList [ (lb, mu) | lb <- [0 .. 2*l], mu <- [lb, lb-1 .. (-lb)] ]

  in V.map func ps

rotMZ :: Int -> Double -> M.Matrix (Complex Double)
rotMZ l omega = let
  d = M.fromList $ map (\m -> cis (fromIntegral (m) * omega)) [-l .. l]
  in M.diagRect (0 :+ 0) d (2*l+1) (2*l+1)

rotMZReal :: Int -> Double -> M.Matrix Double
rotMZReal l omega = let
  calcAtM (i,j)
    |   i  == 0 && j == 0 = 1
    |   i  == j           = cos w
    | (-i) == j && i < 0  = sin w
    | (-i) == j && i > 0  = -(sin w)
    | otherwise           = 0
    where w = (fromIntegral i)*omega
  in ((2*l+1) >< (2*l+1)) $ map calcAtM [(i, j) | i <- [-l..l], j <- [-l..l]]
