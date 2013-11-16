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
import           Numeric.Container     (Product)

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
  us   = rotSHMVec  lmax r
  us'  = rotSHMVec' lmax p
  in V.izipWith (\l u u' -> rotHSHM' u u' (dbg "n:" $ 2*l)) us us'

-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotHSHM :: Int -> SO3 -> SO3 -> M.Matrix (Complex Double)
rotHSHM n r p = let
  l  = n `quot` 2
  u  = (rotSHMVec  l r) V.! l
  u' = (rotSHMVec' l p) V.! l
  in rotHSHM' u u' n

-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotHSHM' :: M.Matrix (Complex Double) -> M.Matrix (Complex Double) -> Int -> M.Matrix (Complex Double)
rotHSHM' u u' n = let
  l  = n `quot` 2

  calcAtM (row, col) = let
    p' = ps V.! row
    p  = ps V.! col
    us = [(-l) .. l]
    in sum [ func p' p (m', m) (m''', m'')
           | m <- us, m' <- us, m'' <- us, m''' <- us

           --, let m''' = mu  - m
           --, let m''  = mu' - m'
           --, m''' >= (-l) && m''' <= l
           --, m''  >= (-l) && m''  <= l
           ]

  func (lb', mu') (lb, mu) (m', m) (m''', m'') = let
    -- calcClebshGordan j1 m1 j2 m2 j m
    c1 = calcClebshGordan (l, m'' ) (l, m') (lb', mu')
    c2 = calcClebshGordan (l, m''') (l, m ) (lb , mu )
    z  = u  @@> (l + m'  , l + m  )
    z' = u' @@> (l + m''', l + m'')
    in ((c1 * c2) :+ 0) * z * z'

  --ps = V.fromList [ (lb, mu) | lb <- [0 .. 2*l], mu <- [lb, lb-1 .. (-lb)] ]
  ps = V.fromList [ (lb, mu) | lb <- [0 .. 2*l], mu <- [(-lb) .. lb] ]
  ms = [ (row, col) | col <- [0 .. s-1], row <- [0 .. s-1] ]
  s  = (2*l+1) * (2*l+1)
  in (s >< s) $ map calcAtM ms

-- ===================================== SH rotation =====================================

-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotSHMVec' :: Int -> SO3 -> Vector (M.Matrix (Complex Double))
rotSHMVec' llen r = let
  z = genSHFunc (2*llen) r

  func l (m', m) (L lb, MF mu) = let
    ld  = fromIntegral l
    lbd = fromIntegral lb
    k1  = (pi / (2 * ld + 1)) * sqrt (2 * (2 * lbd + 1))
    -- calcClebshGordan j1 m1 j2 m2 j m
    c1 = calcClebshGordan (l, m') (lb, mu) (l, m)
    z1 = conjugate $ z %! (N (2*l), L lb, MF mu)
    in ((k1 * c1) :+ 0) * z1

  calcAtM l (m', m) = let
    ps = genLinSeq (2*llen)
    in L.foldl' (\acc p -> acc + func l (m', m) p) 0 ps
  mL l = ((2*l+1) >< (2*l+1)) $ map (calcAtM l) [(i, j) | i <- [-l .. l], j <- [-l .. l]]
  in V.generate (llen+1) mL

-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotSHMVec :: Int -> SO3 -> Vector (M.Matrix (Complex Double))
rotSHMVec llen r = let
  z = genSHFunc (2*llen) r

  func l (m', m) (L lb, MF mu) = let
    ld  = fromIntegral l
    lbd = fromIntegral lb
    k1  = (pi / (2 * ld + 1)) * sqrt (2 * (2 * lbd + 1))
    -- calcClebshGordan j1 m1 j2 m2 j m
    c1 = calcClebshGordan (l, m') (lb, mu) (l, m)
    z1 = z %! (N (2*l), L lb, MF mu)
    in ((k1 * c1) :+ 0) * z1

  calcAtM l (m', m) = let
    ps = genLinSeq (2*llen)
    in L.foldl' (\acc p -> acc + func l (m', m) p) 0 ps
  mL l = ((2*l+1) >< (2*l+1)) $ map (calcAtM l) [(i, j) | i <- [-l .. l], j <- [-l .. l]]
  in V.generate (llen+1) mL

rotSHRealMVec :: Int -> SO3 -> Vector (M.Matrix Double)
rotSHRealMVec l r = let
  vu = rotSHMVec l r
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
    slicer = U.slice (getKeyPos (l, MF (-li))) (2 * li + 1)

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

-- =================================== From Complex ======================================

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

testRotHSH2 :: IO ()
testRotHSH2 = let
  foo  = (rotHSHMVec 6 (SO3 0 0 0) (SO3 (-pi/2) (pi/4) 0) V.!)
  func = fromComplexNLMF . multNStack foo
  in testHSH_C func

testRotHSH :: IO ()
testRotHSH = let
  foo  = (rotSHMVec 6 (SO3 (pi/4) (pi/4) 0) V.!)
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

  rot = so3ToQuaternion $ SO3 (pi/4) (pi/4) 0
  q1  = so3ToQuaternion s1
  q2  = so3ToQuaternion s2
  ss  = V.fromList [s1, s2, quaternionToSO3 $ q1 #<= rot, quaternionToSO3 $ q2 #<= rot]
  ps  = renderSO3PointsVTK ss
  in do
    writeQuater "Rot-SHS-test-points" ps
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

rotHSHRow0 :: N -> SO3 -> Vector (Complex Double)
rotHSHRow0 n r = let
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
